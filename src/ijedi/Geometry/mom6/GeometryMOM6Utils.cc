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

void readHgrid(const std::string & path,
               int niGlobal,
               int njGlobal,
               int niEff,
               int njEff,
               int coarsenFactor,
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
  int ncid;
  if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
    throw eckit::CantOpenFile(path, Here());

  const int nxp = 2 * niGlobal + 1;
  const int nyp = 2 * njGlobal + 1;
  const int nx = 2 * niGlobal;
  const int nySub = 2 * njGlobal;

  std::vector<double> xFull(static_cast<size_t>(nyp) * nxp);
  std::vector<double> yFull(static_cast<size_t>(nyp) * nxp);
  std::vector<double> dxFull(static_cast<size_t>(nyp) * nx);
  std::vector<double> dyFull(static_cast<size_t>(nySub) * nxp);
  std::vector<double> aFull(static_cast<size_t>(nySub) * nx);

  int vid;
  nc_inq_varid(ncid, "x", &vid);
  nc_get_var_double(ncid, vid, xFull.data());
  nc_inq_varid(ncid, "y", &vid);
  nc_get_var_double(ncid, vid, yFull.data());
  nc_inq_varid(ncid, "dx", &vid);
  nc_get_var_double(ncid, vid, dxFull.data());
  nc_inq_varid(ncid, "dy", &vid);
  nc_get_var_double(ncid, vid, dyFull.data());
  nc_inq_varid(ncid, "area", &vid);
  nc_get_var_double(ncid, vid, aFull.data());
  nc_close(ncid);

  const size_t nT = static_cast<size_t>(njEff) * niEff;
  lon->resize(nT);
  lat->resize(nT);
  dxT->resize(nT);
  dyT->resize(nT);
  areaT->resize(nT);
  lonU->resize(nT);
  latU->resize(nT);
  lonV->resize(nT);
  latV->resize(nT);
  cosRot->resize(nT);
  sinRot->resize(nT);

  // Shift x to within half a period of xc (MOM6 modulo_around_point).
  auto moduloAround = [](double x, double xc) {
    return x - 360.0 * std::floor((x - xc) / 360.0 + 0.5);
  };
  const double piOver720 = std::atan(1.0) / 180.0;

  const int K = coarsenFactor;
  for (int jG = 0; jG < njEff; ++jG) {
    const int jS = K * (2 * jG + 1);
    const int j2 = 2 * K * jG;
    for (int iG = 0; iG < niEff; ++iG) {
      const int iS = K * (2 * iG + 1);
      const int i2 = 2 * K * iG;
      const int n = jG * niEff + iG;

      (*lon)[n] = xFull[jS * nxp + iS];
      (*lat)[n] = yFull[jS * nxp + iS];
      (*lonU)[n] = xFull[jS * nxp + (i2 + 2 * K)];
      (*latU)[n] = yFull[jS * nxp + (i2 + 2 * K)];
      (*lonV)[n] = xFull[(jS + K) * nxp + iS];
      (*latV)[n] = yFull[(jS + K) * nxp + iS];

      // Rotation of grid "north" clockwise from true north, from the four
      // cell corners -- MOM6 initialize_grid_rotation_angle (non-bugs branch).
      {
        const int jSo = j2, jNo = j2 + 2 * K, iWe = i2, iEa = i2 + 2 * K;
        const double lonC = (*lon)[n];
        const double lonSW = moduloAround(xFull[jSo * nxp + iWe], lonC);
        const double lonSE = moduloAround(xFull[jSo * nxp + iEa], lonC);
        const double lonNW = moduloAround(xFull[jNo * nxp + iWe], lonC);
        const double lonNE = moduloAround(xFull[jNo * nxp + iEa], lonC);
        const double latSW = yFull[jSo * nxp + iWe];
        const double latSE = yFull[jSo * nxp + iEa];
        const double latNW = yFull[jNo * nxp + iWe];
        const double latNE = yFull[jNo * nxp + iEa];
        const double lonScale = std::cos(piOver720 * ((latSW + latNE) + (latSE + latNW)));
        const double angle = std::atan2(lonScale * ((lonNW - lonSE) + (lonNE - lonSW)),
                                        (latNW - latSE) + (latNE - latSW));
        (*cosRot)[n] = std::cos(angle);
        (*sinRot)[n] = std::sin(angle);
      }

      (*dxT)[n] = 0.0;
      for (int dc = 0; dc < 2 * K; ++dc)
        (*dxT)[n] += dxFull[jS * nx + i2 + dc];

      (*dyT)[n] = 0.0;
      for (int dr = 0; dr < 2 * K; ++dr)
        (*dyT)[n] += dyFull[(j2 + dr) * nxp + iS];

      (*areaT)[n] = 0.0;
      for (int dr = 0; dr < 2 * K; ++dr)
        for (int dc = 0; dc < 2 * K; ++dc)
          (*areaT)[n] += aFull[(j2 + dr) * nx + i2 + dc];
    }
  }
}

void readTopog(const std::string & path,
               int niGlobal,
               int njGlobal,
               int niEff,
               int njEff,
               int coarsenFactor,
               double minimumDepth,
               std::vector<double> * depth,
               std::vector<double> * wet) {
  int ncid;
  if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
    throw eckit::CantOpenFile(path, Here());

  const size_t n = static_cast<size_t>(njGlobal) * niGlobal;
  depth->resize(n);
  wet->resize(n);
  int vid;
  nc_inq_varid(ncid, "depth", &vid);
  nc_get_var_double(ncid, vid, depth->data());

  if (nc_inq_varid(ncid, "wet", &vid) == NC_NOERR) {
    nc_get_var_double(ncid, vid, wet->data());
  } else {
    for (size_t k = 0; k < n; ++k)
      (*wet)[k] = ((*depth)[k] > minimumDepth) ? 1.0 : 0.0;
  }
  nc_close(ncid);

  if (coarsenFactor > 1) {
    const int K = coarsenFactor;
    const size_t nEff = static_cast<size_t>(njEff) * niEff;
    std::vector<double> depthC(nEff), wetC(nEff);
    for (int jG = 0; jG < njEff; ++jG) {
      for (int iG = 0; iG < niEff; ++iG) {
        double sumD = 0.0;
        for (int jr = 0; jr < K; ++jr)
          for (int ir = 0; ir < K; ++ir)
            sumD += (*depth)[(K * jG + jr) * niGlobal + (K * iG + ir)];
            depthC[jG * niEff + iG] = sumD / (K * K);
            wetC[jG * niEff + iG] =
              (depthC[jG * niEff + iG] > minimumDepth) ? 1.0 : 0.0;
      }
    }
    *depth = std::move(depthC);
    *wet = std::move(wetC);
  }
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

}  // namespace

void readVerticalGeometry(const std::string & path,
                         int niGlobal,
                         int njGlobal,
                         int niEff,
                         int njEff,
                         int numLevels,
                         int coarsenFactor,
                         std::vector<double> * layerThickness,
                         std::vector<double> * layerCenterDepth) {
  int ncid;
  if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
    throw eckit::CantOpenFile(path, Here());

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

  int vid;
  if (nc_inq_varid(ncid, "h", &vid) != NC_NOERR) {
    nc_close(ncid);
    throw eckit::Exception(
        "MOM6 vertical geometry: variable 'h' not found in " + path, Here());
  }

  const size_t nFull = static_cast<size_t>(numLevels) * njGlobal * niGlobal;
  std::vector<double> hFull(nFull);
  const size_t start[4] = {0, 0, 0, 0};
  const size_t count[4] = {1,
                           static_cast<size_t>(numLevels),
                           static_cast<size_t>(njGlobal),
                           static_cast<size_t>(niGlobal)};
  if (nc_get_vara_double(ncid, vid, start, count, hFull.data()) != NC_NOERR) {
    nc_close(ncid);
    throw eckit::Exception(
        "MOM6 vertical geometry: failed reading 'h' from " + path, Here());
  }
  nc_close(ncid);

  if (coarsenFactor > 1) {
    const int K = coarsenFactor;
    const size_t nEff = static_cast<size_t>(numLevels) * njEff * niEff;
    std::vector<double> hCoarse(nEff, 0.0);
    for (int k = 0; k < numLevels; ++k) {
      const size_t kFine = static_cast<size_t>(k) * njGlobal * niGlobal;
      const size_t kEff = static_cast<size_t>(k) * njEff * niEff;
      for (int jG = 0; jG < njEff; ++jG) {
        for (int iG = 0; iG < niEff; ++iG) {
          double sumH = 0.0;
          for (int jr = 0; jr < K; ++jr) {
            for (int ir = 0; ir < K; ++ir) {
              sumH += hFull[kFine + static_cast<size_t>(K * jG + jr) * niGlobal +
                            static_cast<size_t>(K * iG + ir)];
            }
          }
          hCoarse[kEff + static_cast<size_t>(jG) * niEff + iG] =
              sumH / static_cast<double>(K * K);
        }
      }
    }
    *layerThickness = std::move(hCoarse);
  } else {
    *layerThickness = std::move(hFull);
  }

  const size_t nHoriz = static_cast<size_t>(njEff) * niEff;
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
