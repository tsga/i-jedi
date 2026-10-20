#include "ijedi/Io/mom6/VelocityStagger.h"

#include <string>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"

namespace ijedi {
namespace mom6 {

// ---------------------------------------------------------------------------
const std::vector<std::pair<std::string, std::string>> & velocityPairs() {
  static const std::vector<std::pair<std::string, std::string>> pairs = {
    {"eastward_sea_water_velocity", "northward_sea_water_velocity"},
    {"barotropic_eastward_sea_water_velocity", "barotropic_northward_sea_water_velocity"},
  };
  return pairs;
}

// ---------------------------------------------------------------------------
void facesToTracer(const std::vector<double> & u, const std::vector<double> & v,
                   int ni, int nj, int nz, bool uSymmetric, bool vSymmetric,
                   std::vector<double> * uT, std::vector<double> * vT) {
  const int niU = uSymmetric ? ni + 1 : ni;
  const int njV = vSymmetric ? nj + 1 : nj;
  ASSERT(u.size() == static_cast<size_t>(nz) * nj * niU);
  ASSERT(v.size() == static_cast<size_t>(nz) * njV * ni);

  const size_t nh = static_cast<size_t>(nj) * ni;
  uT->assign(static_cast<size_t>(nz) * nh, 0.0);
  vT->assign(static_cast<size_t>(nz) * nh, 0.0);

  for (int k = 0; k < nz; ++k) {
    const size_t offU = static_cast<size_t>(k) * nj * niU;
    const size_t offV = static_cast<size_t>(k) * njV * ni;
    const size_t offT = static_cast<size_t>(k) * nh;
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        // u: west face + east face of cell (i, j)
        double uW, uE;
        if (uSymmetric) {
          uW = u[offU + static_cast<size_t>(j) * niU + i];
          uE = u[offU + static_cast<size_t>(j) * niU + i + 1];
        } else {
          const int iw = (i == 0) ? ni - 1 : i - 1;
          uW = u[offU + static_cast<size_t>(j) * ni + iw];
          uE = u[offU + static_cast<size_t>(j) * ni + i];
        }
        // v: south face + north face of cell (i, j)
        double vS, vN;
        if (vSymmetric) {
          vS = v[offV + static_cast<size_t>(j) * ni + i];
          vN = v[offV + static_cast<size_t>(j + 1) * ni + i];
        } else {
          vS = (j == 0) ? 0.0 : v[offV + static_cast<size_t>(j - 1) * ni + i];
          vN = v[offV + static_cast<size_t>(j) * ni + i];
        }
        const size_t t = offT + static_cast<size_t>(j) * ni + i;
        (*uT)[t] = 0.5 * (uW + uE);
        (*vT)[t] = 0.5 * (vS + vN);
      }
    }
  }
}

// ---------------------------------------------------------------------------
void tracerToFaces(const std::vector<double> & uT, const std::vector<double> & vT,
                   const std::vector<double> & mask, int ni, int nj, int nz,
                   bool hasFold, std::vector<double> * u, std::vector<double> * v) {
  const size_t nh = static_cast<size_t>(nj) * ni;
  ASSERT(uT.size() == static_cast<size_t>(nz) * nh);
  ASSERT(vT.size() == static_cast<size_t>(nz) * nh);
  ASSERT(mask.size() >= static_cast<size_t>(nz) * nh);

  u->assign(static_cast<size_t>(nz) * nh, 0.0);
  v->assign(static_cast<size_t>(nz) * nh, 0.0);

  for (int k = 0; k < nz; ++k) {
    const size_t off = static_cast<size_t>(k) * nh;
    auto at = [&](int j, int i) { return off + static_cast<size_t>(j) * ni + i; };
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        const size_t c = at(j, i);

        // East face: cell (i, j) and its periodic eastern neighbour
        const size_t e = at(j, (i + 1) % ni);
        (*u)[c] = mask[c] * mask[e] * 0.5 * (uT[c] + uT[e]);

        // North face: cell (i, j) and its northern neighbour (fold image on
        // the last row, where grid-relative vectors change sign)
        if (j < nj - 1) {
          const size_t n = at(j + 1, i);
          (*v)[c] = mask[c] * mask[n] * 0.5 * (vT[c] + vT[n]);
        } else if (hasFold) {
          const size_t n = at(j, ni - 1 - i);
          (*v)[c] = mask[c] * mask[n] * 0.5 * (vT[c] - vT[n]);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
void gridToGeo(const std::vector<double> & cosRot, const std::vector<double> & sinRot,
               int nz, std::vector<double> * u, std::vector<double> * v) {
  const size_t nh = cosRot.size();
  ASSERT(sinRot.size() == nh);
  ASSERT(u->size() == static_cast<size_t>(nz) * nh);
  ASSERT(v->size() == static_cast<size_t>(nz) * nh);
  for (int k = 0; k < nz; ++k) {
    for (size_t n = 0; n < nh; ++n) {
      const size_t t = static_cast<size_t>(k) * nh + n;
      const double c = cosRot[n], s = sinRot[n];
      const double ug = (*u)[t], vg = (*v)[t];
      (*u)[t] = c * ug + s * vg;
      (*v)[t] = c * vg - s * ug;
    }
  }
}

// ---------------------------------------------------------------------------
void geoToGrid(const std::vector<double> & cosRot, const std::vector<double> & sinRot,
               int nz, std::vector<double> * u, std::vector<double> * v) {
  const size_t nh = cosRot.size();
  ASSERT(sinRot.size() == nh);
  ASSERT(u->size() == static_cast<size_t>(nz) * nh);
  ASSERT(v->size() == static_cast<size_t>(nz) * nh);
  for (int k = 0; k < nz; ++k) {
    for (size_t n = 0; n < nh; ++n) {
      const size_t t = static_cast<size_t>(k) * nh + n;
      const double c = cosRot[n], s = sinRot[n];
      const double ue = (*u)[t], vn = (*v)[t];
      (*u)[t] = c * ue - s * vn;
      (*v)[t] = c * vn + s * ue;
    }
  }
}

}  // namespace mom6
}  // namespace ijedi
