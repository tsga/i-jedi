#include "ijedi/Io/mom6/ReadMOM6Netcdf.h"

#include <netcdf.h>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/mom6/GeometryMOM6.h"

namespace ijedi {

namespace {

std::string findDimName(int ncid,
                        const std::string & role,
                        const std::vector<std::string> & candidates,
                        size_t * len) {
  int dimid;
  for (const auto & name : candidates) {
    if (nc_inq_dimid(ncid, name.c_str(), &dimid) == NC_NOERR) {
      nc_inq_dimlen(ncid, dimid, len);
      return name;
    }
  }

  std::string tried;
  for (const auto & name : candidates) {
    if (!tried.empty()) tried += ", ";
    tried += name;
  }
  throw eckit::Exception(
      "readMOM6Netcdf: cannot find " + role +
          " dimension. Tried: " + tried,
      Here());
}

// Layout of one horizontal-structured variable, from its own dimensions.
struct VarLayout {
  int ni = 0, nj = 0, nz = 1;
  bool hasTime = false;
  bool hasZ = false;
  bool xStaggered = false;  // on lonq/xq (u points)
  bool yStaggered = false;  // on latq/yq (v points)
};

VarLayout inquireLayout(int ncid, int varid) {
  int ndims = 0;
  nc_inq_varndims(ncid, varid, &ndims);
  std::vector<int> dimids(ndims);
  nc_inq_vardimid(ncid, varid, dimids.data());
  int unlimdim = -1;
  nc_inq_unlimdim(ncid, &unlimdim);

  VarLayout lay;
  const std::set<std::string> xq = {"lonq", "xq"}, yq = {"latq", "yq"};
  std::vector<size_t> lens;
  for (int d = 0; d < ndims; ++d) {
    char name[NC_MAX_NAME + 1];
    size_t len = 0;
    nc_inq_dim(ncid, dimids[d], name, &len);
    const std::string dname(name);
    if (d == 0 && (dimids[d] == unlimdim || dname == "Time" || dname == "time")) {
      lay.hasTime = true;
      continue;
    }
    if (d == ndims - 1) lay.xStaggered = xq.count(dname) > 0;
    if (d == ndims - 2) lay.yStaggered = yq.count(dname) > 0;
    lens.push_back(len);
  }
  ASSERT_MSG(lens.size() == 2 || lens.size() == 3,
             "readMOM6Netcdf: expected ([z,] y, x) dimensions");
  lay.ni = static_cast<int>(lens.back());
  lay.nj = static_cast<int>(lens[lens.size() - 2]);
  if (lens.size() == 3) {
    lay.hasZ = true;
    lay.nz = static_cast<int>(lens[0]);
  }
  return lay;
}

// Read the first time record of rows [j0, j0+nj) and columns [i0, i0+ni) of
// a variable into a [nz][nj][ni] buffer.
std::vector<double> readSlab(int ncid, int varid, const VarLayout & lay,
                             size_t j0, size_t nj, size_t i0, size_t ni) {
  std::vector<double> buf(static_cast<size_t>(lay.nz) * nj * ni);
  std::vector<size_t> start, count;
  if (lay.hasTime) { start.push_back(0); count.push_back(1); }
  if (lay.hasZ) { start.push_back(0); count.push_back(static_cast<size_t>(lay.nz)); }
  start.push_back(j0); count.push_back(nj);
  start.push_back(i0); count.push_back(ni);
  if (!buf.empty())
    ASSERT(nc_get_vara_double(ncid, varid, start.data(), count.data(), buf.data())
           == NC_NOERR);
  return buf;
}

// Read the first time record of a variable into a [nz][nj][ni] buffer.
std::vector<double> readStructured(int ncid, int varid, const VarLayout & lay) {
  return readSlab(ncid, varid, lay, 0, lay.nj, 0, lay.ni);
}

// C-grid faces bounding one block, as a symmetric layout of the block:
// u [nz][nj][ni+1] (west face of the first column first), v [nz][nj+1][ni]
// (south face of the first row first).  A non-symmetric u wraps periodically
// in i; a non-symmetric v takes the face south of row 0 as 0 (closed).
void readBlockFaces(int ncid, int varidU, int varidV,
                    const VarLayout & layU, const VarLayout & layV,
                    bool uSym, bool vSym, int niGlobal, const mom6::Block & b,
                    std::vector<double> * u, std::vector<double> * v) {
  const size_t ni = b.ni, nj = b.nj, i0 = b.i0, j0 = b.j0;
  const int nz = layU.nz;
  if (uSym) {
    *u = readSlab(ncid, varidU, layU, j0, nj, i0, ni + 1);
  } else if (i0 > 0) {
    *u = readSlab(ncid, varidU, layU, j0, nj, i0 - 1, ni + 1);
  } else {
    const std::vector<double> west = readSlab(ncid, varidU, layU, j0, nj, niGlobal - 1, 1);
    const std::vector<double> rest = readSlab(ncid, varidU, layU, j0, nj, 0, ni);
    u->resize(static_cast<size_t>(nz) * nj * (ni + 1));
    for (size_t kj = 0; kj < static_cast<size_t>(nz) * nj; ++kj) {
      (*u)[kj * (ni + 1)] = west[kj];
      std::copy_n(rest.begin() + kj * ni, ni, u->begin() + kj * (ni + 1) + 1);
    }
  }
  if (vSym) {
    *v = readSlab(ncid, varidV, layV, j0, nj + 1, i0, ni);
  } else if (j0 > 0) {
    *v = readSlab(ncid, varidV, layV, j0 - 1, nj + 1, i0, ni);
  } else {
    const std::vector<double> rest = readSlab(ncid, varidV, layV, 0, nj, i0, ni);
    v->assign(static_cast<size_t>(nz) * (nj + 1) * ni, 0.0);
    for (int k = 0; k < nz; ++k)
      std::copy_n(rest.begin() + static_cast<size_t>(k) * nj * ni, nj * ni,
                  v->begin() + static_cast<size_t>(k) * (nj + 1) * ni + ni);
  }
}

}  // namespace

// ---------------------------------------------------------------------------
void readMOM6Netcdf(const std::vector<std::string> & filepaths,
                    atlas::FieldSet & x,
                    const std::vector<std::string> & fileVarNames,
                    const std::vector<double> & scalings,
                    const eckit::mpi::Comm & comm,
                    const mom6::VelocityGridInfo * velGrid) {
  const bool isRoot = (comm.rank() == 0);
  const auto fs =
      atlas::functionspace::NodeColumns(x.field(0).functionspace());

  // --- 1. Gather global indices to root (collective, done once). ---
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
  fs.gather(gidxField, globalGidxField);  // collective

  auto gIdx = atlas::array::make_view<double, 2>(globalGidxField);
  const int nGlobal = static_cast<int>(globalGidxField.shape(0));

  // --- 2. Pre-allocate one global field per output field (all ranks). ---
  std::vector<atlas::Field> globalFields;
  globalFields.reserve(x.size());
  for (const auto & field : x) {
    globalFields.push_back(fs.createField<double>(
        atlas::option::name(field.name())
        | atlas::option::levels(field.shape(1))
        | atlas::option::global()));
  }

  // --- 3. Root reads each file in order; a later file may overwrite an
  //        earlier one, but in practice variables don't overlap. ---
  if (isRoot) {
    for (const auto & filepath : filepaths) {
      int ncid = -1;
      ASSERT_MSG(nc_open(filepath.c_str(), NC_NOWRITE, &ncid) == NC_NOERR,
                 "readMOM6Netcdf: cannot open " + filepath);

      size_t dimNi = 0, dimNj = 0, dimNz = 0;
      const std::string dimX =
        findDimName(ncid, "x (longitude)", {"lonh", "xh", "xaxis_1", "nx", "ni"}, &dimNi);
      const std::string dimY =
        findDimName(ncid, "y (latitude)", {"lath", "yh", "yaxis_1", "ny", "nj"}, &dimNj);

      // z dimension is optional: 2D files (sea ice, fix) may not have one
      bool hasZDim = false;
      for (const auto & zName : std::vector<std::string>{"Layer", "z_l", "zaxis_1", "nz"}) {
        int dimid;
        if (nc_inq_dimid(ncid, zName.c_str(), &dimid) == NC_NOERR) {
          nc_inq_dimlen(ncid, dimid, &dimNz);
          hasZDim = true;
          break;
        }
      }

      const int file_ni = static_cast<int>(dimNi);
      const int file_nj = static_cast<int>(dimNj);
      const int file_nz = static_cast<int>(dimNz);
      const int nStructured = file_ni * file_nj;

      for (int n = 0; n < nGlobal; ++n) {
        const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
        if (flatIdx < 0 || flatIdx >= nStructured) {
          nc_close(ncid);
          throw eckit::BadValue(
              "readMOM6Netcdf: global index " + std::to_string(flatIdx + 1) +
                  " is out of bounds for file horizontal dimensions " +
                  dimX + "=" + std::to_string(file_ni) +
                  ", " + dimY + "=" + std::to_string(file_nj) +
                  " in " + filepath,
              Here());
        }
      }

      // C-grid velocity pairs: average faces to tracer points, then rotate
      // grid-relative -> eastward/northward.  Velocities stored on h-points
      // (e.g. written by i-jedi with "velocity grid: tracer") are already
      // geographic and go through the generic path below.
      std::set<int> handled;
      for (const auto & pair : mom6::velocityPairs()) {
        int iE = -1, iN = -1;
        for (int f = 0; f < x.size(); ++f) {
          if (fileVarNames[f].empty()) continue;
          if (x.field(f).name() == pair.first)  iE = f;
          if (x.field(f).name() == pair.second) iN = f;
        }
        int varidE = -1, varidN = -1;
        const bool foundE = iE >= 0 &&
            nc_inq_varid(ncid, fileVarNames[iE].c_str(), &varidE) == NC_NOERR;
        const bool foundN = iN >= 0 &&
            nc_inq_varid(ncid, fileVarNames[iN].c_str(), &varidN) == NC_NOERR;
        if (!foundE && !foundN) continue;

        VarLayout layE, layN;
        if (foundE) layE = inquireLayout(ncid, varidE);
        if (foundN) layN = inquireLayout(ncid, varidN);
        const bool staggered = (foundE && layE.xStaggered) || (foundN && layN.yStaggered);
        if (!staggered) continue;

        if (!(foundE && foundN)) {
          nc_close(ncid);
          throw eckit::UserError(
              "readMOM6Netcdf: C-grid velocity in " + filepath + " needs both " +
              pair.first + " and " + pair.second + " to be requested and present",
              Here());
        }
        ASSERT_MSG(velGrid != nullptr,
                   "readMOM6Netcdf: C-grid velocity requires grid rotation info");
        const int ni = velGrid->ni, nj = velGrid->nj;
        const bool uSym = (layE.ni == ni + 1), vSym = (layN.nj == nj + 1);
        ASSERT_MSG(layE.xStaggered && layN.yStaggered &&
                   (uSym || layE.ni == ni) && layE.nj == nj &&
                   (vSym || layN.nj == nj) && layN.ni == ni,
                   "readMOM6Netcdf: unexpected C-grid velocity dimensions for " +
                   fileVarNames[iE] + "/" + fileVarNames[iN] + " in " + filepath);

        const int nLevels = x.field(iE).shape(1);
        if (layE.nz != nLevels || layN.nz != nLevels) {
          oops::Log::warning() << "readMOM6Netcdf: skipping C-grid velocity '"
                               << fileVarNames[iE] << "'/'" << fileVarNames[iN]
                               << "' from " << filepath
                               << " (z dimension size mismatch)" << std::endl;
          handled.insert(iE);
          handled.insert(iN);
          continue;
        }

        std::vector<double> u = readStructured(ncid, varidE, layE);
        std::vector<double> v = readStructured(ncid, varidN, layN);
        if (scalings[iE] != 0.0) for (auto & val : u) val *= scalings[iE];
        if (scalings[iN] != 0.0) for (auto & val : v) val *= scalings[iN];

        std::vector<double> uT, vT;
        mom6::facesToTracer(u, v, ni, nj, nLevels, uSym, vSym, &uT, &vT);
        mom6::gridToGeo(velGrid->cosRot, velGrid->sinRot, nLevels, &uT, &vT);

        const size_t nh = static_cast<size_t>(nj) * ni;
        auto gE = atlas::array::make_view<double, 2>(globalFields[iE]);
        auto gN = atlas::array::make_view<double, 2>(globalFields[iN]);
        for (int n = 0; n < nGlobal; ++n) {
          const size_t flatIdx = static_cast<size_t>(gIdx(n, 0)) - 1;
          for (int k = 0; k < nLevels; ++k) {
            gE(n, k) = uT[k * nh + flatIdx];
            gN(n, k) = vT[k * nh + flatIdx];
          }
        }
        handled.insert(iE);
        handled.insert(iN);
        oops::Log::info() << "readMOM6Netcdf: " << pair.first << ", " << pair.second
                          << " <- C-grid " << fileVarNames[iE] << ", " << fileVarNames[iN]
                          << " (averaged to h-points, rotated) from " << filepath
                          << std::endl;
      }

      int fi = 0;
      for (const auto & field : x) {
        if (handled.count(fi)) { ++fi; continue; }
        const std::string & fileVarName = fileVarNames[fi];
        const double scale = scalings[fi];
        const int nLevels = field.shape(1);

        if (!fileVarName.empty()) {
          int varid;
          const bool found =
              (nc_inq_varid(ncid, fileVarName.c_str(), &varid) == NC_NOERR);

          if (found) {
            auto gView = atlas::array::make_view<double, 2>(globalFields[fi]);
            const bool is3D = (nLevels > 1);

            // Detect whether the variable has a leading time dimension.
            // Files like MOM6 restarts use (time, [z,] nj, ni); geometry
            // cache files written by saveStructuredGrid use ([z,] nj, ni).
            int var_ndims = 0;
            nc_inq_varndims(ncid, varid, &var_ndims);
            const bool hasTimeDim = (var_ndims == (is3D ? 4 : 3));

            if (is3D) {
              if (!hasZDim || file_nz != nLevels) {
                oops::Log::warning()
                    << "readMOM6Netcdf: skipping 3D variable '" << fileVarName
                    << "' from " << filepath
                    << " (z dimension absent or size mismatch)" << std::endl;
              } else {
                std::vector<double> buf(
                    static_cast<size_t>(file_nz) * file_nj * file_ni);
                if (hasTimeDim) {
                  const size_t start[4] = {0, 0, 0, 0};
                  const size_t count[4] = {1, static_cast<size_t>(file_nz),
                                              static_cast<size_t>(file_nj),
                                              static_cast<size_t>(file_ni)};
                  nc_get_vara_double(ncid, varid, start, count, buf.data());
                } else {
                  const size_t start[3] = {0, 0, 0};
                  const size_t count[3] = {static_cast<size_t>(file_nz),
                                              static_cast<size_t>(file_nj),
                                              static_cast<size_t>(file_ni)};
                  nc_get_vara_double(ncid, varid, start, count, buf.data());
                }

                if (scale != 0.0)
                  for (auto & v : buf) v *= scale;

                for (int n = 0; n < nGlobal; ++n) {
                  const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
                  const int iG = flatIdx % file_ni;
                  const int jG = flatIdx / file_ni;
                  for (int k = 0; k < nLevels; ++k) {
                    gView(n, k) = buf[static_cast<size_t>(k) * file_nj * file_ni
                                    + static_cast<size_t>(jG) * file_ni
                                    + static_cast<size_t>(iG)];
                  }
                }
                oops::Log::info() << "readMOM6Netcdf: " << field.name()
                                  << " <- " << fileVarName
                                  << " from " << filepath << std::endl;
              }
            } else {
              std::vector<double> buf(
                  static_cast<size_t>(file_nj) * file_ni);
              if (hasTimeDim) {
                const size_t start[3] = {0, 0, 0};
                const size_t count[3] = {1, static_cast<size_t>(file_nj),
                                            static_cast<size_t>(file_ni)};
                nc_get_vara_double(ncid, varid, start, count, buf.data());
              } else {
                const size_t start[2] = {0, 0};
                const size_t count[2] = {static_cast<size_t>(file_nj),
                                            static_cast<size_t>(file_ni)};
                nc_get_vara_double(ncid, varid, start, count, buf.data());
              }

              if (scale != 0.0)
                for (auto & v : buf) v *= scale;

              for (int n = 0; n < nGlobal; ++n) {
                const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
                gView(n, 0) = buf[flatIdx];
              }
              oops::Log::info() << "readMOM6Netcdf: " << field.name()
                                << " <- " << fileVarName
                                << " from " << filepath << std::endl;
            }
          }
        }
        ++fi;
      }

      nc_close(ncid);
    }
  }

  // --- 4. Scatter all global fields to distributed ranks (collective). ---
  // Zero the distributed fields first. atlas::createField leaves memory
  // uninitialised, and scatter only writes nodes whose global index is present
  // in the gathered global field. With the northern-fold / periodic mesh, some
  // owned nodes are not covered by the gather/scatter round-trip; without this
  // they would retain heap garbage, which is then summed by reductions
  // (e.g. normField over exclude_halo) and makes results non-reproducible.
  for (auto & field : x) {
    auto v = atlas::array::make_view<double, 2>(field);
    v.assign(0.0);
  }

  int fi = 0;
  for (auto & field : x) {
    fs.scatter(globalFields[fi], field);
    ++fi;
  }

  // Halo exchange to fill ghost nodes from their owning ranks
  for (auto & field : x) {
    atlas::functionspace::NodeColumns(field.functionspace()).haloExchange(field);
  }
}

// ---------------------------------------------------------------------------
void readMOM6NetcdfBlocks(const std::vector<std::string> & filepaths,
                          atlas::FieldSet & x,
                          const std::vector<std::string> & fileVarNames,
                          const std::vector<double> & scalings,
                          const GeometryMOM6 & geom) {
  const bool isRoot = (geom.comm().rank() == 0);
  const mom6::Block block = geom.mom6Block();
  const size_t nLocal = block.size();

  // Owned JEDI nodes are overwritten by the scatter, ghosts by the halo
  // exchange; zero first so a variable missing from every file reads as 0,
  // as in readMOM6Netcdf.
  for (auto & field : x) {
    auto v = atlas::array::make_view<double, 2>(field);
    v.assign(0.0);
  }

  // Block (level-major) → MOM6-space field → owned JEDI nodes
  auto toJedi = [&](const std::vector<double> & buf, atlas::Field * field) {
    const int nLevels = field->shape(1);
    atlas::Field mom6Field = geom.mom6FunctionSpace().createField<double>(
        atlas::option::name(field->name()) | atlas::option::levels(nLevels));
    auto vm = atlas::array::make_view<double, 2>(mom6Field);
    for (size_t n = 0; n < nLocal; ++n)
      for (int k = 0; k < nLevels; ++k)
        vm(n, k) = buf[static_cast<size_t>(k) * nLocal + n];
    geom.scatterToJedi(mom6Field, field);
  };

  // This block's grid rotation, for C-grid velocities
  std::vector<double> cosRot(nLocal), sinRot(nLocal);
  {
    auto vCos = atlas::array::make_view<double, 2>(geom.mom6Fields().field("cos_rot"));
    auto vSin = atlas::array::make_view<double, 2>(geom.mom6Fields().field("sin_rot"));
    for (size_t n = 0; n < nLocal; ++n) {
      cosRot[n] = vCos(n, 0);
      sinRot[n] = vSin(n, 0);
    }
  }

  for (const auto & filepath : filepaths) {
    int ncid = -1;
    ASSERT_MSG(nc_open(filepath.c_str(), NC_NOWRITE, &ncid) == NC_NOERR,
               "readMOM6NetcdfBlocks: cannot open " + filepath);

    size_t dimNi = 0, dimNj = 0, dimNz = 0;
    findDimName(ncid, "x (longitude)", {"lonh", "xh", "xaxis_1", "nx", "ni"}, &dimNi);
    findDimName(ncid, "y (latitude)", {"lath", "yh", "yaxis_1", "ny", "nj"}, &dimNj);
    if (static_cast<int>(dimNi) != geom.niEff() || static_cast<int>(dimNj) != geom.njEff()) {
      nc_close(ncid);
      throw eckit::BadValue("readMOM6NetcdfBlocks: " + filepath + " is " +
                            std::to_string(dimNi) + " x " + std::to_string(dimNj) +
                            " but the geometry grid is " + std::to_string(geom.niEff()) +
                            " x " + std::to_string(geom.njEff()), Here());
    }
    int zDimId = -1;
    for (const auto & zName : std::vector<std::string>{"Layer", "z_l", "zaxis_1", "nz"}) {
      if (nc_inq_dimid(ncid, zName.c_str(), &zDimId) == NC_NOERR) {
        nc_inq_dimlen(ncid, zDimId, &dimNz);
        break;
      }
      zDimId = -1;
    }

    // C-grid velocity pairs: read the faces bounding this task's block,
    // average them to tracer points, then rotate grid-relative ->
    // eastward/northward.  Velocities stored on h-points (e.g. written by
    // i-jedi with "velocity grid: tracer") are already geographic and go
    // through the generic path below.
    std::set<int> handled;
    for (const auto & pair : mom6::velocityPairs()) {
      int iE = -1, iN = -1;
      for (int f = 0; f < x.size(); ++f) {
        if (fileVarNames[f].empty()) continue;
        if (x.field(f).name() == pair.first)  iE = f;
        if (x.field(f).name() == pair.second) iN = f;
      }
      int varidE = -1, varidN = -1;
      const bool foundE = iE >= 0 &&
          nc_inq_varid(ncid, fileVarNames[iE].c_str(), &varidE) == NC_NOERR;
      const bool foundN = iN >= 0 &&
          nc_inq_varid(ncid, fileVarNames[iN].c_str(), &varidN) == NC_NOERR;
      if (!foundE && !foundN) continue;

      VarLayout layE, layN;
      if (foundE) layE = inquireLayout(ncid, varidE);
      if (foundN) layN = inquireLayout(ncid, varidN);
      const bool staggered = (foundE && layE.xStaggered) || (foundN && layN.yStaggered);
      if (!staggered) continue;

      if (!(foundE && foundN)) {
        nc_close(ncid);
        throw eckit::UserError(
            "readMOM6NetcdfBlocks: C-grid velocity in " + filepath + " needs both " +
            pair.first + " and " + pair.second + " to be requested and present",
            Here());
      }
      const int ni = geom.niEff(), nj = geom.njEff();
      const bool uSym = (layE.ni == ni + 1), vSym = (layN.nj == nj + 1);
      ASSERT_MSG(layE.xStaggered && layN.yStaggered &&
                 (uSym || layE.ni == ni) && layE.nj == nj &&
                 (vSym || layN.nj == nj) && layN.ni == ni,
                 "readMOM6NetcdfBlocks: unexpected C-grid velocity dimensions for " +
                 fileVarNames[iE] + "/" + fileVarNames[iN] + " in " + filepath);

      const int nLevels = x.field(iE).shape(1);
      if (layE.nz != nLevels || layN.nz != nLevels) {
        if (isRoot)
          oops::Log::warning() << "readMOM6NetcdfBlocks: skipping C-grid velocity '"
                               << fileVarNames[iE] << "'/'" << fileVarNames[iN]
                               << "' from " << filepath
                               << " (z dimension size mismatch)" << std::endl;
        handled.insert(iE);
        handled.insert(iN);
        continue;
      }

      std::vector<double> u, v;
      readBlockFaces(ncid, varidE, varidN, layE, layN, uSym, vSym, ni, block, &u, &v);
      if (scalings[iE] != 0.0) for (auto & val : u) val *= scalings[iE];
      if (scalings[iN] != 0.0) for (auto & val : v) val *= scalings[iN];

      std::vector<double> uT, vT;
      mom6::facesToTracer(u, v, block.ni, block.nj, nLevels, true, true, &uT, &vT);
      mom6::gridToGeo(cosRot, sinRot, nLevels, &uT, &vT);
      atlas::Field fieldE = x.field(iE), fieldN = x.field(iN);
      toJedi(uT, &fieldE);
      toJedi(vT, &fieldN);

      handled.insert(iE);
      handled.insert(iN);
      if (isRoot)
        oops::Log::info() << "readMOM6NetcdfBlocks: " << pair.first << ", " << pair.second
                          << " <- C-grid " << fileVarNames[iE] << ", " << fileVarNames[iN]
                          << " (averaged to h-points, rotated) from " << filepath
                          << std::endl;
    }

    int fi = 0;
    for (auto & field : x) {
      if (handled.count(fi)) { ++fi; continue; }
      const std::string & fileVarName = fileVarNames[fi];
      const double scale = scalings[fi];
      ++fi;
      const int nLevels = field.shape(1);
      int varid;
      if (fileVarName.empty() ||
          nc_inq_varid(ncid, fileVarName.c_str(), &varid) != NC_NOERR) continue;

      // Variable layout: ([time,] [z,] j, i).  It has a z dimension when the
      // one before (j, i) is the file's z dimension.
      int ndims = 0;
      nc_inq_varndims(ncid, varid, &ndims);
      std::vector<int> dimids(std::max(ndims, 1));
      nc_inq_vardimid(ncid, varid, dimids.data());
      const bool varHasZ = (ndims >= 3 && zDimId >= 0 && dimids[ndims - 3] == zDimId);

      const bool is3D = (nLevels > 1);
      if (is3D && (!varHasZ || static_cast<int>(dimNz) != nLevels)) {
        if (isRoot)
          oops::Log::warning() << "readMOM6NetcdfBlocks: skipping 3D field " << field.name()
                               << ": variable '" << fileVarName << "' in " << filepath
                               << " has no matching z dimension" << std::endl;
        continue;
      }

      // Hyperslab of this task's block; a leading time dimension is read at
      // index 0, and a 2D field takes level 0 of a variable with levels.
      std::vector<size_t> start(ndims, 0), count(ndims, 1);
      if (is3D) count[ndims - 3] = static_cast<size_t>(nLevels);
      start[ndims - 2] = static_cast<size_t>(block.j0);
      count[ndims - 2] = static_cast<size_t>(block.nj);
      start[ndims - 1] = static_cast<size_t>(block.i0);
      count[ndims - 1] = static_cast<size_t>(block.ni);
      std::vector<double> buf(static_cast<size_t>(nLevels) * nLocal);
      if (nLocal > 0)
        ASSERT_MSG(nc_get_vara_double(ncid, varid, start.data(), count.data(), buf.data())
                   == NC_NOERR, "readMOM6NetcdfBlocks: cannot read " + fileVarName
                   + " from " + filepath);
      if (scale != 0.0)
        for (auto & v : buf) v *= scale;

      toJedi(buf, &field);

      if (isRoot)
        oops::Log::info() << "readMOM6NetcdfBlocks: " << field.name() << " <- "
                          << fileVarName << " from " << filepath << std::endl;
    }
    nc_close(ncid);
  }

  // Halo exchange to fill ghost nodes from their owning tasks
  for (auto & field : x) {
    atlas::functionspace::NodeColumns(field.functionspace()).haloExchange(field);
  }
}

// ---------------------------------------------------------------------------
}  // namespace ijedi
