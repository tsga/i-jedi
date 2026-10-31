#include "ijedi/Io/mom6/WriteMOM6Netcdf.h"

#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/mom6/GeometryMOM6.h"

namespace ijedi {

namespace {
void ncCheck(int status, const std::string & what, const std::string & filepath) {
  if (status != NC_NOERR)
    throw eckit::Exception("writeMOM6Netcdf: " + what + " failed for " + filepath + ": "
                           + nc_strerror(status), Here());
}

// uT, vT and mask (level-major on this task's block, nz levels) at the cells
// bordering the block to the east and north, fetched from their owners.  The
// east neighbour of the last global column is the first (periodic).  The
// north neighbour of the last global row is its tripolar fold image
// (nj-1, ni-1-i), where grid-relative v changes sign; without a fold that
// face is closed (mask 0).
void fetchEastNorth(const GeometryMOM6 & geom, const int nz,
                    const std::vector<double> & uT, const std::vector<double> & vT,
                    const std::vector<double> & mask,
                    std::vector<double> * uEast, std::vector<double> * mEast,
                    std::vector<double> * vNorth, std::vector<double> * mNorth) {
  const eckit::mpi::Comm & comm = geom.comm();
  const int npes = static_cast<int>(comm.size());
  const mom6::Block b = geom.mom6Block();
  const int niG = geom.niEff(), njG = geom.njEff();
  const size_t nLocal = b.size();

  // Slots: east column (j) first, then north row (i)
  struct Request { int rank, localIdx, slot; double vSign; };
  std::vector<Request> requests;
  for (int j = 0; j < b.nj; ++j) {
    const auto owner = geom.mom6Owner((b.i0 + b.ni) % niG, b.j0 + j);
    requests.push_back({owner.first, owner.second, j, 1.0});
  }
  for (int i = 0; i < b.ni; ++i) {
    if (b.j0 + b.nj < njG) {
      const auto owner = geom.mom6Owner(b.i0 + i, b.j0 + b.nj);
      requests.push_back({owner.first, owner.second, b.nj + i, 1.0});
    } else if (geom.hasFold()) {
      const auto owner = geom.mom6Owner(niG - 1 - (b.i0 + i), njG - 1);
      requests.push_back({owner.first, owner.second, b.nj + i, -1.0});
    }
  }
  std::stable_sort(requests.begin(), requests.end(),
                   [](const Request & a, const Request & c) { return a.rank < c.rank; });

  std::vector<int> reqCounts(npes, 0), srvCounts(npes, 0);
  for (const auto & r : requests) reqCounts[r.rank]++;
  comm.allToAll(reqCounts, srvCounts);
  std::vector<int> reqDispl(npes, 0), srvDispl(npes, 0);
  for (int r = 1; r < npes; ++r) {
    reqDispl[r] = reqDispl[r - 1] + reqCounts[r - 1];
    srvDispl[r] = srvDispl[r - 1] + srvCounts[r - 1];
  }
  const int nSrv = srvDispl[npes - 1] + srvCounts[npes - 1];

  std::vector<int> reqIdx(requests.size()), srvIdx(nSrv);
  for (size_t n = 0; n < requests.size(); ++n) reqIdx[n] = requests[n].localIdx;
  comm.allToAllv(reqIdx.data(), reqCounts.data(), reqDispl.data(),
                 srvIdx.data(), srvCounts.data(), srvDispl.data());

  // Serve (uT, vT, mask) x nz per requested point
  const int nv = 3 * nz;
  std::vector<double> srvBuf(static_cast<size_t>(nSrv) * nv);
  for (int n = 0; n < nSrv; ++n)
    for (int k = 0; k < nz; ++k) {
      const size_t src = static_cast<size_t>(k) * nLocal + srvIdx[n];
      double * dst = &srvBuf[static_cast<size_t>(n) * nv + 3 * k];
      dst[0] = uT[src];
      dst[1] = vT[src];
      dst[2] = mask[src];
    }
  auto perValue = [nv](std::vector<int> c) {
    for (auto & x : c) x *= nv;
    return c;
  };
  std::vector<double> reqBuf(requests.size() * nv);
  comm.allToAllv(srvBuf.data(), perValue(srvCounts).data(), perValue(srvDispl).data(),
                 reqBuf.data(), perValue(reqCounts).data(), perValue(reqDispl).data());

  uEast->assign(static_cast<size_t>(nz) * b.nj, 0.0);
  mEast->assign(uEast->size(), 0.0);
  vNorth->assign(static_cast<size_t>(nz) * b.ni, 0.0);
  mNorth->assign(vNorth->size(), 0.0);
  for (size_t n = 0; n < requests.size(); ++n) {
    const Request & r = requests[n];
    for (int k = 0; k < nz; ++k) {
      const double * src = &reqBuf[n * nv + 3 * k];
      if (r.slot < b.nj) {
        (*uEast)[static_cast<size_t>(k) * b.nj + r.slot] = src[0];
        (*mEast)[static_cast<size_t>(k) * b.nj + r.slot] = src[2];
      } else {
        const int i = r.slot - b.nj;
        (*vNorth)[static_cast<size_t>(k) * b.ni + i] = r.vSign * src[1];
        (*mNorth)[static_cast<size_t>(k) * b.ni + i] = src[2];
      }
    }
  }
}
}  // namespace

// ---------------------------------------------------------------------------
void writeMOM6Netcdf(const std::string & filepath,
                     const atlas::FieldSet & x,
                     const std::vector<std::string> & fileVarNames,
                     const int ni, const int nj, const int nz,
                     const eckit::mpi::Comm & comm,
                     const bool nativeVelocity,
                     const mom6::VelocityGridInfo * velGrid) {
  const bool isRoot = (comm.rank() == 0);

  // Velocity pairs to write on the native MOM6 C-grid (grid-relative, u on
  // east faces, v on north faces).  eastIdx -> northIdx.
  std::unordered_map<int, int> nativePairs;
  std::unordered_set<int> nativeNorth;
  if (nativeVelocity) {
    for (const auto & pair : mom6::velocityPairs()) {
      int iE = -1, iN = -1;
      for (int f = 0; f < x.size(); ++f) {
        if (fileVarNames[f].empty()) continue;
        if (x.field(f).name() == pair.first)  iE = f;
        if (x.field(f).name() == pair.second) iN = f;
      }
      if (iE < 0 && iN < 0) continue;
      if (iE < 0 || iN < 0)
        throw eckit::UserError(
            "writeMOM6Netcdf: velocity grid native needs both " + pair.first +
            " and " + pair.second, Here());
      nativePairs[iE] = iN;
      nativeNorth.insert(iN);
    }
    ASSERT_MSG(nativePairs.empty() || !isRoot || velGrid != nullptr,
               "writeMOM6Netcdf: native velocity requires grid rotation info");
  }

  // --- 1. Gather every field to rank 0 using Atlas NodeColumns::gather ---
  const auto fs =
      atlas::functionspace::NodeColumns(x.field(0).functionspace());

  // The JEDI function space is unstructured and only contains active points
  // (ocean + fringe); land points are omitted.  After fs.gather() the global
  // field has nActiveGlobal entries -- NOT ni*nj.  We therefore also gather
  // the global-index field so we know where each gathered node sits on the
  // full structured grid.
  //
  // Convention: global_index = jG * ni + iG + 1  (1-based)

  // 1a. Build a local field holding each node's global index as double,
  //     then gather it to root.
  auto localGidx = atlas::array::make_view<atlas::gidx_t, 1>(
      fs.nodes().global_index());
  atlas::Field gidxField = fs.createField<double>(
      atlas::option::name("__gidx") | atlas::option::levels(1));
  {
    auto gv = atlas::array::make_view<double, 2>(gidxField);
    for (atlas::idx_t n = 0; n < fs.nodes().size(); ++n)
      gv(n, 0) = static_cast<double>(localGidx(n));
  }
  atlas::Field globalGidxField = fs.createField<double>(
      atlas::option::name("__gidx") | atlas::option::levels(1)
      | atlas::option::global());
  fs.gather(gidxField, globalGidxField);

  // 1b. Gather all data fields.
  std::vector<atlas::Field> globalFields;
  globalFields.reserve(x.size());

  int fi = 0;
  for (const auto & field : x) {
    const std::string & varName = fileVarNames[fi++];
    if (varName.empty()) {
      globalFields.emplace_back();          // placeholder
      continue;
    }
    const int nLevels = field.shape(1);
    atlas::Field gf = fs.createField<double>(
        atlas::option::name(varName)
        | atlas::option::levels(nLevels)
        | atlas::option::global());
    fs.gather(field, gf);
    globalFields.push_back(gf);
  }

  // --- 2. Root creates the NetCDF file, defines dims/vars, writes data ---
  if (isRoot) {
    std::remove(filepath.c_str());
    // NetCDF-4: classic formats cap a variable at 2-4 GB, and a 3D double
    // field reaches ~4.9 GB at 1/12 deg (41 x 3297 x 4500).
    int ncid;
    ncCheck(nc_create(filepath.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid), "nc_create", filepath);

    // Define dimensions
    int dim_time, dim_layer, dim_lath, dim_lonh;
    ncCheck(nc_def_dim(ncid, "Time",  NC_UNLIMITED, &dim_time), "nc_def_dim", filepath);
    ncCheck(nc_def_dim(ncid, "Layer", nz,           &dim_layer), "nc_def_dim", filepath);
    ncCheck(nc_def_dim(ncid, "lath",  nj,           &dim_lath), "nc_def_dim", filepath);
    ncCheck(nc_def_dim(ncid, "lonh",  ni,           &dim_lonh), "nc_def_dim", filepath);
    int dim_latq = -1, dim_lonq = -1;
    if (!nativePairs.empty()) {
      // Non-symmetric MOM6 layout: u(I) on the east face, v(J) on the north face
      ncCheck(nc_def_dim(ncid, "latq", nj, &dim_latq), "nc_def_dim", filepath);
      ncCheck(nc_def_dim(ncid, "lonq", ni, &dim_lonq), "nc_def_dim", filepath);
    }

    // Define variables (skip duplicates, e.g. two JEDI names → same file var)
    std::unordered_map<std::string, int> varids;
    fi = 0;
    for (const auto & field : x) {
      const std::string & varName = fileVarNames[fi++];
      if (varName.empty() || varids.count(varName)) continue;
      int varid;
      const int fIdx = fi - 1;
      const int dimY = nativeNorth.count(fIdx) ? dim_latq : dim_lath;
      const int dimX = nativePairs.count(fIdx) ? dim_lonq : dim_lonh;
      if (field.shape(1) > 1) {
        int dimids[4] = {dim_time, dim_layer, dimY, dimX};
        ncCheck(nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 4, dimids, &varid),
                "nc_def_var " + varName, filepath);
      } else {
        int dimids[3] = {dim_time, dimY, dimX};
        ncCheck(nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 3, dimids, &varid),
                "nc_def_var " + varName, filepath);
      }
      varids[varName] = varid;
    }
    ncCheck(nc_enddef(ncid), "nc_enddef", filepath);

    // Retrieve the gathered global indices (only valid on root)
    auto gIdx = atlas::array::make_view<double, 2>(globalGidxField);
    const int nGathered = static_cast<int>(globalGidxField.shape(0));

    // --- Diagnostic: validate gathered global indices ---
    {
      int nBadIdx = 0;
      double minIdx = 1.0e30, maxIdx = -1.0e30;
      for (int n = 0; n < nGathered; ++n) {
        const double g = gIdx(n, 0);
        if (g < minIdx) minIdx = g;
        if (g > maxIdx) maxIdx = g;
        const int flatIdx = static_cast<int>(g) - 1;
        if (flatIdx < 0 || flatIdx >= nj * ni) ++nBadIdx;
      }
      oops::Log::info() << "writeMOM6Netcdf: gathered " << nGathered
                        << " nodes (full grid = " << nj * ni << ")"
                        << ", gIdx range [" << minIdx << ", " << maxIdx << "]"
                        << ", out-of-range count = " << nBadIdx << std::endl;
    }

    // Expand a gathered field onto the full (nLevels, nj, ni) grid.
    auto expand = [&](const atlas::Field & gf) {
      auto gv = atlas::array::make_view<double, 2>(gf);
      const int nL = static_cast<int>(gf.shape(1));
      std::vector<double> full(static_cast<size_t>(nL) * nj * ni, 0.0);
      for (int n = 0; n < static_cast<int>(gf.shape(0)); ++n) {
        const size_t flatIdx = static_cast<size_t>(gIdx(n, 0)) - 1;
        for (int k = 0; k < nL; ++k)
          full[static_cast<size_t>(k) * nj * ni + flatIdx] = gv(n, k);
      }
      return full;
    };
    auto put = [&](const std::string & varName, int nL, const std::vector<double> & buf) {
      if (nL > 1) {
        const size_t start[4] = {0, 0, 0, 0};
        const size_t count[4] = {1, static_cast<size_t>(nL),
                                    static_cast<size_t>(nj),
                                    static_cast<size_t>(ni)};
        ncCheck(nc_put_vara_double(ncid, varids.at(varName), start, count, buf.data()),
                "nc_put_vara_double " + varName, filepath);
      } else {
        const size_t start[3] = {0, 0, 0};
        const size_t count[3] = {1, static_cast<size_t>(nj),
                                    static_cast<size_t>(ni)};
        ncCheck(nc_put_vara_double(ncid, varids.at(varName), start, count, buf.data()),
                "nc_put_vara_double " + varName, filepath);
      }
    };

    // Write each gathered field
    fi = 0;
    for (const auto & field : x) {
      const std::string & varName = fileVarNames[fi];
      const auto & gf = globalFields[fi];
      const int fIdx = fi;
      ++fi;
      if (varName.empty()) continue;

      if (nativeNorth.count(fIdx)) continue;  // written with its east partner
      if (nativePairs.count(fIdx)) {
        // Geographic at h-points -> grid-relative, then h-points -> C-grid faces
        const int iN = nativePairs.at(fIdx);
        const int nL = static_cast<int>(gf.shape(1));
        ASSERT(static_cast<int>(globalFields[iN].shape(1)) == nL);
        std::vector<double> uT = expand(gf);
        std::vector<double> vT = expand(globalFields[iN]);
        mom6::geoToGrid(velGrid->cosRot, velGrid->sinRot, nL, &uT, &vT);
        std::vector<double> u, v;
        mom6::tracerToFaces(uT, vT, velGrid->mask3d, ni, nj, nL, velGrid->hasFold, &u, &v);
        put(varName, nL, u);
        put(fileVarNames[iN], nL, v);
        oops::Log::info() << "writeMOM6Netcdf: " << field.name() << ", "
                          << x.field(iN).name() << " -> C-grid " << varName << ", "
                          << fileVarNames[iN] << " (rotated, interpolated to faces) to "
                          << filepath << std::endl;
        continue;
      }

      auto gView = atlas::array::make_view<double, 2>(gf);
      const int nGlobal = static_cast<int>(gf.shape(0));
      const int nLevels = static_cast<int>(gf.shape(1));

      // The gathered field has nActiveGlobal entries (ocean+fringe only).
      // Use the gathered global_index to place each value on the full grid.
      //   global_index = jG * ni + iG + 1  →  flatIdx = global_index - 1
      if (nLevels > 1) {
        const size_t nHoriz = static_cast<size_t>(nj) * ni;
        std::vector<double> buf(static_cast<size_t>(nz) * nHoriz, 0.0);
        for (int n = 0; n < nGlobal; ++n) {
          const size_t flatIdx = static_cast<size_t>(gIdx(n, 0)) - 1;
          for (int k = 0; k < nLevels; ++k)
            buf[static_cast<size_t>(k) * nHoriz + flatIdx] = gView(n, k);
        }
        const size_t start[4] = {0, 0, 0, 0};
        const size_t count[4] = {1, static_cast<size_t>(nz),
                                    static_cast<size_t>(nj),
                                    static_cast<size_t>(ni)};
        ncCheck(nc_put_vara_double(ncid, varids.at(varName), start, count, buf.data()),
                "nc_put_vara_double " + varName, filepath);
      } else {
        std::vector<double> buf(static_cast<size_t>(nj * ni), 0.0);
        for (int n = 0; n < nGlobal; ++n) {
          const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
          buf[flatIdx] = gView(n, 0);
        }
        const size_t start[3] = {0, 0, 0};
        const size_t count[3] = {1, static_cast<size_t>(nj),
                                    static_cast<size_t>(ni)};
        ncCheck(nc_put_vara_double(ncid, varids.at(varName), start, count, buf.data()),
                "nc_put_vara_double " + varName, filepath);
      }

      oops::Log::info() << "writeMOM6Netcdf: " << field.name()
                        << " -> " << varName << " to " << filepath << std::endl;
    }

    ncCheck(nc_close(ncid), "nc_close", filepath);
  }

  // Ensure the file is fully written before any rank returns (and potentially
  // tries to read it back).
  comm.barrier();
}

// ---------------------------------------------------------------------------
void writeMOM6NetcdfBlocks(const std::string & filepath,
                           const atlas::FieldSet & x,
                           const std::vector<std::string> & fileVarNames,
                           const int nz,
                           const GeometryMOM6 & geom,
                           const bool nativeVelocity) {
  const eckit::mpi::Comm & comm = geom.comm();
  const bool isRoot = (comm.rank() == 0);
  const mom6::Block block = geom.mom6Block();
  const size_t nLocal = block.size();

  // Velocity pairs to write on the native MOM6 C-grid (grid-relative, u on
  // east faces, v on north faces).  eastIdx -> northIdx.
  std::unordered_map<int, int> nativePairs;
  std::unordered_set<int> nativeNorth;
  if (nativeVelocity) {
    for (const auto & pair : mom6::velocityPairs()) {
      int iE = -1, iN = -1;
      for (int f = 0; f < x.size(); ++f) {
        if (fileVarNames[f].empty()) continue;
        if (x.field(f).name() == pair.first)  iE = f;
        if (x.field(f).name() == pair.second) iN = f;
      }
      if (iE < 0 && iN < 0) continue;
      if (iE < 0 || iN < 0)
        throw eckit::UserError(
            "writeMOM6NetcdfBlocks: velocity grid native needs both " + pair.first +
            " and " + pair.second, Here());
      nativePairs[iE] = iN;
      nativeNorth.insert(iN);
    }
  }

  if (isRoot) std::remove(filepath.c_str());
  comm.barrier();

  int ncid;
  ncCheck(nc_create_par(filepath.c_str(), NC_CLOBBER | NC_NETCDF4,
                        MPI_Comm_f2c(comm.communicator()), MPI_INFO_NULL, &ncid),
          "nc_create_par", filepath);

  int dim_time, dim_layer, dim_lath, dim_lonh;
  ncCheck(nc_def_dim(ncid, "Time",  NC_UNLIMITED, &dim_time), "nc_def_dim", filepath);
  ncCheck(nc_def_dim(ncid, "Layer", nz,           &dim_layer), "nc_def_dim", filepath);
  ncCheck(nc_def_dim(ncid, "lath",  geom.njEff(), &dim_lath), "nc_def_dim", filepath);
  ncCheck(nc_def_dim(ncid, "lonh",  geom.niEff(), &dim_lonh), "nc_def_dim", filepath);
  int dim_latq = -1, dim_lonq = -1;
  if (!nativePairs.empty()) {
    // Non-symmetric MOM6 layout: u(I) on the east face, v(J) on the north face
    ncCheck(nc_def_dim(ncid, "latq", geom.njEff(), &dim_latq), "nc_def_dim", filepath);
    ncCheck(nc_def_dim(ncid, "lonq", geom.niEff(), &dim_lonq), "nc_def_dim", filepath);
  }

  // Define variables (skip duplicates, e.g. two JEDI names → same file var)
  std::unordered_map<std::string, int> varids;
  int fi = 0;
  for (const auto & field : x) {
    const std::string & varName = fileVarNames[fi++];
    if (varName.empty() || varids.count(varName)) continue;
    int varid;
    const int fIdx = fi - 1;
    const int dimY = nativeNorth.count(fIdx) ? dim_latq : dim_lath;
    const int dimX = nativePairs.count(fIdx) ? dim_lonq : dim_lonh;
    if (field.shape(1) > 1) {
      int dimids[4] = {dim_time, dim_layer, dimY, dimX};
      ncCheck(nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 4, dimids, &varid),
              "nc_def_var " + varName, filepath);
    } else {
      int dimids[3] = {dim_time, dimY, dimX};
      ncCheck(nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 3, dimids, &varid),
              "nc_def_var " + varName, filepath);
    }
    varids[varName] = varid;
  }
  ncCheck(nc_enddef(ncid), "nc_enddef", filepath);

  // Owned JEDI values → this task's MOM6 block, level-major (land stays 0)
  auto toBlock = [&](const atlas::Field & field) {
    const int nLevels = static_cast<int>(field.shape(1));
    atlas::Field mom6Field = geom.mom6FunctionSpace().createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(nLevels));
    auto vm = atlas::array::make_view<double, 2>(mom6Field);
    vm.assign(0.0);
    geom.gatherFromJedi(field, &mom6Field);
    std::vector<double> buf(static_cast<size_t>(nLevels) * nLocal);
    for (size_t n = 0; n < nLocal; ++n)
      for (int k = 0; k < nLevels; ++k)
        buf[static_cast<size_t>(k) * nLocal + n] = vm(n, k);
    return buf;
  };

  // Collective hyperslab write of this task's block
  auto putBlock = [&](const std::string & varName, const int nLevels,
                      const std::vector<double> & buf) {
    const int varid = varids.at(varName);
    ncCheck(nc_var_par_access(ncid, varid, NC_COLLECTIVE), "nc_var_par_access", filepath);
    if (nLevels > 1) {
      const size_t start[4] = {0, 0, static_cast<size_t>(block.j0),
                               static_cast<size_t>(block.i0)};
      const size_t count[4] = {1, static_cast<size_t>(nLevels),
                               static_cast<size_t>(block.nj), static_cast<size_t>(block.ni)};
      ncCheck(nc_put_vara_double(ncid, varid, start, count, buf.data()),
              "nc_put_vara_double " + varName, filepath);
    } else {
      const size_t start[3] = {0, static_cast<size_t>(block.j0), static_cast<size_t>(block.i0)};
      const size_t count[3] = {1, static_cast<size_t>(block.nj), static_cast<size_t>(block.ni)};
      ncCheck(nc_put_vara_double(ncid, varid, start, count, buf.data()),
              "nc_put_vara_double " + varName, filepath);
    }
  };

  fi = 0;
  for (const auto & field : x) {
    const int fIdx = fi;
    const std::string & varName = fileVarNames[fi++];
    if (varName.empty()) continue;
    const int nLevels = static_cast<int>(field.shape(1));

    if (nativeNorth.count(fIdx)) continue;  // written with its east partner
    if (nativePairs.count(fIdx)) {
      // Geographic at h-points -> grid-relative, then h-points -> C-grid faces
      const int iN = nativePairs.at(fIdx);
      const atlas::Field & fieldN = x.field(iN);
      ASSERT(static_cast<int>(fieldN.shape(1)) == nLevels);
      std::vector<double> uT = toBlock(field);
      std::vector<double> vT = toBlock(fieldN);

      std::vector<double> cosRot(nLocal), sinRot(nLocal);
      auto vCos = atlas::array::make_view<double, 2>(geom.mom6Fields().field("cos_rot"));
      auto vSin = atlas::array::make_view<double, 2>(geom.mom6Fields().field("sin_rot"));
      for (size_t n = 0; n < nLocal; ++n) {
        cosRot[n] = vCos(n, 0);
        sinRot[n] = vSin(n, 0);
      }
      mom6::geoToGrid(cosRot, sinRot, nLevels, &uT, &vT);

      const atlas::Field & fMask = geom.mom6Fields().field("mask3d");
      ASSERT(nLevels <= static_cast<int>(fMask.shape(1)));
      auto vMask = atlas::array::make_view<double, 2>(fMask);
      std::vector<double> mask(static_cast<size_t>(nLevels) * nLocal);
      for (size_t n = 0; n < nLocal; ++n)
        for (int k = 0; k < nLevels; ++k)
          mask[static_cast<size_t>(k) * nLocal + n] = vMask(n, k);

      std::vector<double> uEast, mEast, vNorth, mNorth, u, v;
      fetchEastNorth(geom, nLevels, uT, vT, mask, &uEast, &mEast, &vNorth, &mNorth);
      mom6::tracerToFacesHalo(uT, vT, mask, block.ni, block.nj, nLevels,
                              uEast, mEast, vNorth, mNorth, &u, &v);
      putBlock(varName, nLevels, u);
      putBlock(fileVarNames[iN], nLevels, v);
      if (isRoot)
        oops::Log::info() << "writeMOM6NetcdfBlocks: " << field.name() << ", "
                          << fieldN.name() << " -> C-grid " << varName << ", "
                          << fileVarNames[iN] << " (rotated, interpolated to faces) to "
                          << filepath << std::endl;
      continue;
    }

    putBlock(varName, nLevels, toBlock(field));
    if (isRoot)
      oops::Log::info() << "writeMOM6NetcdfBlocks: " << field.name()
                        << " -> " << varName << " to " << filepath << std::endl;
  }

  ncCheck(nc_close(ncid), "nc_close", filepath);
  comm.barrier();
}

// ---------------------------------------------------------------------------
}  // namespace ijedi
