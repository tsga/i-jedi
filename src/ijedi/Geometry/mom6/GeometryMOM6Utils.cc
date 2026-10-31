#include "ijedi/Geometry/mom6/GeometryMOM6Utils.h"

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"

namespace ijedi {
namespace mom6 {

int computeExtent(int N, int ndivs, int pe) {
  const int base = N / ndivs;
  const int extra = N % ndivs;
  return base + (pe < extra ? 1 : 0);
}

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
      "MOM6 vertical geometry: cannot find " + role +
          " dimension. Tried: " + tried,
      Here());
}

bool parseLeveledNodeDataName(const std::string & nameLine,
                              std::string * baseName,
                              int * level) {
  if (nameLine.size() < 5 || nameLine.front() != '"' || nameLine.back() != '"')
    return false;
  const size_t lb = nameLine.rfind('[');
  const size_t rb = nameLine.rfind(']');
  if (lb == std::string::npos || rb == std::string::npos || lb <= 1 || rb <= lb + 1)
    return false;
  if (rb != nameLine.size() - 2) return false;

  const std::string lvl = nameLine.substr(lb + 1, rb - lb - 1);
  if (!std::all_of(lvl.begin(), lvl.end(), [](char c) { return c >= '0' && c <= '9'; }))
    return false;

  *baseName = nameLine.substr(1, lb - 1);
  *level = std::stoi(lvl);
  return true;
}

// Read the [r0, r0+nr) x [c0, c0+nc) hyperslab of the two trailing (row,
// column) dimensions of 'var'; any leading dimension (time, ntiles) is taken
// at index 0.
void readSlab2D(int ncid, const std::string & path, const std::string & var,
                size_t r0, size_t nr, size_t c0, size_t nc, double * buf) {
  int vid, ndims;
  if (nc_inq_varid(ncid, var.c_str(), &vid) != NC_NOERR ||
      nc_inq_varndims(ncid, vid, &ndims) != NC_NOERR || ndims < 2)
    throw eckit::Exception("MOM6 geometry: cannot find 2D variable '" + var +
                           "' in " + path, Here());
  std::vector<size_t> start(ndims, 0), count(ndims, 1);
  start[ndims - 2] = r0;  count[ndims - 2] = nr;
  start[ndims - 1] = c0;  count[ndims - 1] = nc;
  if (nc_get_vara_double(ncid, vid, start.data(), count.data(), buf) != NC_NOERR)
    throw eckit::Exception("MOM6 geometry: failed reading '" + var + "' from " + path,
                           Here());
}

int openOrThrow(const std::string & path) {
  int ncid;
  if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
    throw eckit::CantOpenFile(path, Here());
  return ncid;
}

}  // namespace

void readHgridBlock(const std::string & path,
                    int coarsenFactor,
                    const Block & block,
                    std::vector<double> * lon,
                    std::vector<double> * lat,
                    std::vector<double> * dxT,
                    std::vector<double> * dyT,
                    std::vector<double> * areaT,
                    std::vector<double> * lonU,
                    std::vector<double> * latU,
                    std::vector<double> * lonV,
                    std::vector<double> * latV,
                    std::vector<double> * cosRot,
                    std::vector<double> * sinRot) {
  const int K = coarsenFactor;
  // Supergrid corner rows/columns covering the block: x, y are (nyp, nxp),
  // dx is (nyp, nx), dy is (ny, nxp) and area is (ny, nx).
  const size_t r0 = static_cast<size_t>(2 * K) * block.j0;
  const size_t c0 = static_cast<size_t>(2 * K) * block.i0;
  const size_t R  = static_cast<size_t>(2 * K) * block.nj + 1;
  const size_t C  = static_cast<size_t>(2 * K) * block.ni + 1;

  std::vector<double> x(R * C), y(R * C), dx(R * (C - 1)), dy((R - 1) * C),
                      a((R - 1) * (C - 1));
  const int ncid = openOrThrow(path);
  readSlab2D(ncid, path, "x",    r0, R,     c0, C,     x.data());
  readSlab2D(ncid, path, "y",    r0, R,     c0, C,     y.data());
  readSlab2D(ncid, path, "dx",   r0, R,     c0, C - 1, dx.data());
  readSlab2D(ncid, path, "dy",   r0, R - 1, c0, C,     dy.data());
  readSlab2D(ncid, path, "area", r0, R - 1, c0, C - 1, a.data());
  nc_close(ncid);

  const size_t nT = block.size();
  for (auto * v : {lon, lat, dxT, dyT, areaT, lonU, latU, lonV, latV, cosRot, sinRot})
    v->resize(nT);

  // Shift x to within half a period of xc (MOM6 modulo_around_point).
  auto moduloAround = [](double x, double xc) {
    return x - 360.0 * std::floor((x - xc) / 360.0 + 0.5);
  };
  const double piOver720 = std::atan(1.0) / 180.0;

  for (int jL = 0; jL < block.nj; ++jL) {
    const size_t jS = static_cast<size_t>(K) * (2 * jL + 1);
    const size_t j2 = static_cast<size_t>(2 * K) * jL;
    for (int iL = 0; iL < block.ni; ++iL) {
      const size_t iS = static_cast<size_t>(K) * (2 * iL + 1);
      const size_t i2 = static_cast<size_t>(2 * K) * iL;
      const size_t n = static_cast<size_t>(jL) * block.ni + iL;

      (*lon)[n]  = x[jS * C + iS];
      (*lat)[n]  = y[jS * C + iS];
      (*lonU)[n] = x[jS * C + (i2 + 2 * K)];
      (*latU)[n] = y[jS * C + (i2 + 2 * K)];
      (*lonV)[n] = x[(jS + K) * C + iS];
      (*latV)[n] = y[(jS + K) * C + iS];

      // Rotation of grid "north" clockwise from true north, from the four
      // cell corners -- MOM6 initialize_grid_rotation_angle (non-bugs branch).
      {
        const size_t jSo = j2, jNo = j2 + 2 * K, iWe = i2, iEa = i2 + 2 * K;
        const double lonC = (*lon)[n];
        const double lonSW = moduloAround(x[jSo * C + iWe], lonC);
        const double lonSE = moduloAround(x[jSo * C + iEa], lonC);
        const double lonNW = moduloAround(x[jNo * C + iWe], lonC);
        const double lonNE = moduloAround(x[jNo * C + iEa], lonC);
        const double latSW = y[jSo * C + iWe];
        const double latSE = y[jSo * C + iEa];
        const double latNW = y[jNo * C + iWe];
        const double latNE = y[jNo * C + iEa];
        const double lonScale = std::cos(piOver720 * ((latSW + latNE) + (latSE + latNW)));
        const double angle = std::atan2(lonScale * ((lonNW - lonSE) + (lonNE - lonSW)),
                                        (latNW - latSE) + (latNE - latSW));
        (*cosRot)[n] = std::cos(angle);
        (*sinRot)[n] = std::sin(angle);
      }

      (*dxT)[n] = 0.0;
      for (int dc = 0; dc < 2 * K; ++dc)
        (*dxT)[n] += dx[jS * (C - 1) + i2 + dc];

      (*dyT)[n] = 0.0;
      for (int dr = 0; dr < 2 * K; ++dr)
        (*dyT)[n] += dy[(j2 + dr) * C + iS];

      (*areaT)[n] = 0.0;
      for (int dr = 0; dr < 2 * K; ++dr)
        for (int dc = 0; dc < 2 * K; ++dc)
          (*areaT)[n] += a[(j2 + dr) * (C - 1) + i2 + dc];
    }
  }
}

void readHgridLonLat(const std::string & path,
                     int niGlobal,
                     int niEff,
                     int njEff,
                     int coarsenFactor,
                     std::vector<double> * lon,
                     std::vector<double> * lat) {
  const int K = coarsenFactor;
  const size_t nxp = static_cast<size_t>(2 * niGlobal) + 1;
  lon->resize(static_cast<size_t>(njEff) * niEff);
  lat->resize(static_cast<size_t>(njEff) * niEff);

  // Read a strip of T rows at a time: bounded memory, few large reads.
  const int strip = 32;
  std::vector<double> xs, ys;
  const int ncid = openOrThrow(path);
  for (int jG0 = 0; jG0 < njEff; jG0 += strip) {
    const int nRows = std::min(strip, njEff - jG0);
    const size_t rowFirst = static_cast<size_t>(K) * (2 * jG0 + 1);
    const size_t nr = static_cast<size_t>(2 * K) * (nRows - 1) + 1;
    xs.resize(nr * nxp);
    ys.resize(nr * nxp);
    readSlab2D(ncid, path, "x", rowFirst, nr, 0, nxp, xs.data());
    readSlab2D(ncid, path, "y", rowFirst, nr, 0, nxp, ys.data());
    for (int jj = 0; jj < nRows; ++jj) {
      const size_t r = static_cast<size_t>(2 * K) * jj;
      for (int iG = 0; iG < niEff; ++iG) {
        const size_t iS = static_cast<size_t>(K) * (2 * iG + 1);
        const size_t n = static_cast<size_t>(jG0 + jj) * niEff + iG;
        (*lon)[n] = xs[r * nxp + iS];
        (*lat)[n] = ys[r * nxp + iS];
      }
    }
  }
  nc_close(ncid);
}

void readTopogBlock(const std::string & path,
                    int coarsenFactor,
                    double minimumDepth,
                    const Block & block,
                    std::vector<double> * depth,
                    std::vector<double> * wet) {
  const int K = coarsenFactor;
  const size_t nr = static_cast<size_t>(K) * block.nj;
  const size_t nc = static_cast<size_t>(K) * block.ni;
  const size_t r0 = static_cast<size_t>(K) * block.j0;
  const size_t c0 = static_cast<size_t>(K) * block.i0;

  std::vector<double> d(nr * nc);
  const int ncid = openOrThrow(path);
  readSlab2D(ncid, path, "depth", r0, nr, c0, nc, d.data());

  if (K == 1) {
    wet->resize(d.size());
    int vid;
    if (nc_inq_varid(ncid, "wet", &vid) == NC_NOERR) {
      readSlab2D(ncid, path, "wet", r0, nr, c0, nc, wet->data());
    } else {
      for (size_t k = 0; k < d.size(); ++k)
        (*wet)[k] = (d[k] > minimumDepth) ? 1.0 : 0.0;
    }
    nc_close(ncid);
    *depth = std::move(d);
    return;
  }
  nc_close(ncid);

  depth->resize(block.size());
  wet->resize(block.size());
  for (int jL = 0; jL < block.nj; ++jL) {
    for (int iL = 0; iL < block.ni; ++iL) {
      double sumD = 0.0;
      for (int jr = 0; jr < K; ++jr)
        for (int ir = 0; ir < K; ++ir)
          sumD += d[static_cast<size_t>(K * jL + jr) * nc + (K * iL + ir)];
      const size_t n = static_cast<size_t>(jL) * block.ni + iL;
      (*depth)[n] = sumD / (K * K);
      (*wet)[n] = ((*depth)[n] > minimumDepth) ? 1.0 : 0.0;
    }
  }
}

void readVerticalGeometryBlock(const std::string & path,
                               int niGlobal,
                               int njGlobal,
                               int numLevels,
                               int coarsenFactor,
                               const Block & block,
                               std::vector<double> * layerThickness,
                               std::vector<double> * layerCenterDepth) {
  const int ncid = openOrThrow(path);

  size_t fileNi = 0, fileNj = 0, fileNz = 0;
  const std::string dimX =
      findDimName(ncid, "x (longitude)", {"lonh", "xh", "xaxis_1", "nx"}, &fileNi);
  const std::string dimY =
      findDimName(ncid, "y (latitude)", {"lath", "yh", "yaxis_1", "ny"}, &fileNj);
  const std::string dimZ =
      findDimName(ncid, "z (layer)", {"Layer", "z_l", "zaxis_1"}, &fileNz);

  if (static_cast<int>(fileNi) != niGlobal ||
      static_cast<int>(fileNj) != njGlobal ||
      static_cast<int>(fileNz) != numLevels) {
    nc_close(ncid);
    throw eckit::BadValue(
        "MOM6 vertical geometry: file dimensions in " + path +
            " are " + dimX + "=" + std::to_string(fileNi) +
            ", " + dimY + "=" + std::to_string(fileNj) +
            ", " + dimZ + "=" + std::to_string(fileNz) +
            " but geometry expects " + std::to_string(niGlobal) + " x " +
            std::to_string(njGlobal) + " x " + std::to_string(numLevels),
        Here());
  }

  int vid, ndims;
  if (nc_inq_varid(ncid, "h", &vid) != NC_NOERR ||
      nc_inq_varndims(ncid, vid, &ndims) != NC_NOERR || ndims < 3) {
    nc_close(ncid);
    throw eckit::Exception(
        "MOM6 vertical geometry: variable 'h' not found in " + path, Here());
  }

  // Fine-grid block (K x the effective block), all levels; any leading
  // (time) dimension is taken at index 0.
  const int K = coarsenFactor;
  const size_t nr = static_cast<size_t>(K) * block.nj;
  const size_t nc = static_cast<size_t>(K) * block.ni;
  std::vector<size_t> start(ndims, 0), count(ndims, 1);
  count[ndims - 3] = static_cast<size_t>(numLevels);
  start[ndims - 2] = static_cast<size_t>(K) * block.j0;  count[ndims - 2] = nr;
  start[ndims - 1] = static_cast<size_t>(K) * block.i0;  count[ndims - 1] = nc;
  std::vector<double> hBlk(static_cast<size_t>(numLevels) * nr * nc);
  if (nc_get_vara_double(ncid, vid, start.data(), count.data(), hBlk.data()) != NC_NOERR) {
    nc_close(ncid);
    throw eckit::Exception(
        "MOM6 vertical geometry: failed reading 'h' from " + path, Here());
  }
  nc_close(ncid);

  const size_t nHoriz = block.size();
  if (K > 1) {
    std::vector<double> hCoarse(static_cast<size_t>(numLevels) * nHoriz, 0.0);
    for (int k = 0; k < numLevels; ++k) {
      const size_t kFine = static_cast<size_t>(k) * nr * nc;
      const size_t kEff = static_cast<size_t>(k) * nHoriz;
      for (int jL = 0; jL < block.nj; ++jL) {
        for (int iL = 0; iL < block.ni; ++iL) {
          double sumH = 0.0;
          for (int jr = 0; jr < K; ++jr)
            for (int ir = 0; ir < K; ++ir)
              sumH += hBlk[kFine + static_cast<size_t>(K * jL + jr) * nc +
                           static_cast<size_t>(K * iL + ir)];
          hCoarse[kEff + static_cast<size_t>(jL) * block.ni + iL] =
              sumH / static_cast<double>(K * K);
        }
      }
    }
    *layerThickness = std::move(hCoarse);
  } else {
    *layerThickness = std::move(hBlk);
  }

  layerCenterDepth->assign(static_cast<size_t>(numLevels) * nHoriz, 0.0);
  for (size_t n = 0; n < nHoriz; ++n) {
    double depthTop = 0.0;
    for (int k = 0; k < numLevels; ++k) {
      const size_t idx = static_cast<size_t>(k) * nHoriz + n;
      const double thickness = (*layerThickness)[idx];
      (*layerCenterDepth)[idx] = depthTop + 0.5 * thickness;
      depthTop += thickness;
    }
  }
}

void rewriteLeveledNodeDataAsTime(const std::string & filename,
                                  const std::unordered_set<std::string> & timeFields) {
  std::ifstream in(filename);
  if (!in)
    throw eckit::CantOpenFile(filename, Here());

  const std::string tmp = filename + ".tmp_levels_as_time";
  std::ofstream out(tmp);
  if (!out)
    throw eckit::CantOpenFile(tmp, Here());

  std::string line;
  while (std::getline(in, line)) {
    if (line != "$NodeData") {
      out << line << '\n';
      continue;
    }

    out << line << '\n';

    std::string nStringTags, nameLine, nRealTags, realTag;
    std::string nIntTags, int1, int2, int3, int4;
    if (!std::getline(in, nStringTags) ||
        !std::getline(in, nameLine) ||
        !std::getline(in, nRealTags) ||
        !std::getline(in, realTag) ||
        !std::getline(in, nIntTags) ||
        !std::getline(in, int1) ||
        !std::getline(in, int2) ||
        !std::getline(in, int3) ||
        !std::getline(in, int4)) {
      throw eckit::Exception("Malformed $NodeData block while rewriting " + filename,
                             Here());
    }

    std::string baseName;
    int level = 0;
    const bool isLeveled = parseLeveledNodeDataName(nameLine, &baseName, &level);
    if (isLeveled && timeFields.count(baseName) > 0) {
      nameLine = "\"" + baseName + "\"";
      if (nRealTags == "1") realTag = std::to_string(level);
      if (nIntTags == "4") int1 = std::to_string(level);
    }

    out << nStringTags << '\n'
        << nameLine << '\n'
        << nRealTags << '\n'
        << realTag << '\n'
        << nIntTags << '\n'
        << int1 << '\n'
        << int2 << '\n'
        << int3 << '\n'
        << int4 << '\n';
  }

  in.close();
  out.close();

  if (std::rename(tmp.c_str(), filename.c_str()) != 0)
    throw eckit::Exception("Failed replacing Gmsh file after level->time rewrite: " +
                               filename,
                           Here());
}

}  // namespace mom6
}  // namespace ijedi
