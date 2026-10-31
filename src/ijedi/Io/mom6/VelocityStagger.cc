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

  // East neighbour of the last column is the first (periodic); north
  // neighbour of the last row is its fold image (nj-1, ni-1-i), where
  // grid-relative vectors change sign.  Without a fold the north face is closed.
  std::vector<double> uEast(static_cast<size_t>(nz) * nj), mEast(uEast.size());
  std::vector<double> vNorth(static_cast<size_t>(nz) * ni, 0.0), mNorth(vNorth.size(), 0.0);
  for (int k = 0; k < nz; ++k) {
    const size_t off = static_cast<size_t>(k) * nh;
    for (int j = 0; j < nj; ++j) {
      uEast[static_cast<size_t>(k) * nj + j] = uT[off + static_cast<size_t>(j) * ni];
      mEast[static_cast<size_t>(k) * nj + j] = mask[off + static_cast<size_t>(j) * ni];
    }
    if (hasFold) {
      for (int i = 0; i < ni; ++i) {
        const size_t img = off + static_cast<size_t>(nj - 1) * ni + (ni - 1 - i);
        vNorth[static_cast<size_t>(k) * ni + i] = -vT[img];
        mNorth[static_cast<size_t>(k) * ni + i] = mask[img];
      }
    }
  }
  tracerToFacesHalo(uT, vT, mask, ni, nj, nz, uEast, mEast, vNorth, mNorth, u, v);
}

// ---------------------------------------------------------------------------
void tracerToFacesHalo(const std::vector<double> & uT, const std::vector<double> & vT,
                       const std::vector<double> & mask, int ni, int nj, int nz,
                       const std::vector<double> & uEast, const std::vector<double> & mEast,
                       const std::vector<double> & vNorth, const std::vector<double> & mNorth,
                       std::vector<double> * u, std::vector<double> * v) {
  const size_t nh = static_cast<size_t>(nj) * ni;
  ASSERT(uT.size() == static_cast<size_t>(nz) * nh);
  ASSERT(vT.size() == static_cast<size_t>(nz) * nh);
  ASSERT(mask.size() >= static_cast<size_t>(nz) * nh);
  ASSERT(uEast.size() == static_cast<size_t>(nz) * nj && mEast.size() == uEast.size());
  ASSERT(vNorth.size() == static_cast<size_t>(nz) * ni && mNorth.size() == vNorth.size());

  u->assign(static_cast<size_t>(nz) * nh, 0.0);
  v->assign(static_cast<size_t>(nz) * nh, 0.0);

  for (int k = 0; k < nz; ++k) {
    const size_t off = static_cast<size_t>(k) * nh;
    auto at = [&](int j, int i) { return off + static_cast<size_t>(j) * ni + i; };
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        const size_t c = at(j, i);

        // East face: cell (i, j) and its eastern neighbour
        double uE, mE;
        if (i < ni - 1) {
          uE = uT[at(j, i + 1)];
          mE = mask[at(j, i + 1)];
        } else {
          uE = uEast[static_cast<size_t>(k) * nj + j];
          mE = mEast[static_cast<size_t>(k) * nj + j];
        }
        (*u)[c] = mask[c] * mE * 0.5 * (uT[c] + uE);

        // North face: cell (i, j) and its northern neighbour
        double vN, mN;
        if (j < nj - 1) {
          vN = vT[at(j + 1, i)];
          mN = mask[at(j + 1, i)];
        } else {
          vN = vNorth[static_cast<size_t>(k) * ni + i];
          mN = mNorth[static_cast<size_t>(k) * ni + i];
        }
        (*v)[c] = mask[c] * mN * 0.5 * (vT[c] + vN);
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
