#pragma once

#include <string>
#include <utility>
#include <vector>

namespace ijedi {
namespace mom6 {

/// Velocity component pairs (eastward, northward) that live on the MOM6
/// C-grid in model files but at tracer points (geographic) in i-jedi.
const std::vector<std::pair<std::string, std::string>> & velocityPairs();

/// Structured-grid context needed to convert velocities between the MOM6
/// C-grid and i-jedi tracer points.  Arrays are full global, row-major
/// (jG * ni + iG); mask3d is [nz][nj][ni].  Only meaningful on the rank that
/// does the netCDF I/O.
struct VelocityGridInfo {
  int ni = 0;
  int nj = 0;
  int nz = 0;
  bool hasFold = false;
  std::vector<double> cosRot;
  std::vector<double> sinRot;
  std::vector<double> mask3d;
};

// All functions below act on full structured arrays laid out [nz][nj][ni]
// (faces: [nz][nj][ni(+1)] for u, [nz][nj(+1)][ni] for v).

/// Average C-grid faces to tracer points, as MOM6's ssu_east diagnostic:
///   uT(i) = 1/2 (u(I-1) + u(I)),  vT(j) = 1/2 (v(J-1) + v(J))
/// Non-symmetric u wraps periodically in i; non-symmetric v takes v(-1) = 0
/// (closed southern boundary).  A symmetric array carries the extra
/// western/southern face explicitly.
void facesToTracer(const std::vector<double> & u, const std::vector<double> & v,
                   int ni, int nj, int nz, bool uSymmetric, bool vSymmetric,
                   std::vector<double> * uT, std::vector<double> * vT);

/// Interpolate tracer-point components to non-symmetric C-grid faces:
///   u(I) = m(i) m(i+1) 1/2 (uT(i) + uT(i+1))   (periodic in i)
///   v(J) = m(j) m(j+1) 1/2 (vT(j) + vT(j+1))
/// On the last row the northern neighbour is the tripolar fold image
/// (nj-1, ni-1-i) with vT sign-flipped; without a fold that face is 0.
/// mask is [nz][nj][ni] (only the first nz levels are used).
void tracerToFaces(const std::vector<double> & uT, const std::vector<double> & vT,
                   const std::vector<double> & mask, int ni, int nj, int nz,
                   bool hasFold, std::vector<double> * u, std::vector<double> * v);

/// Grid-relative -> geographic at tracer points (in place):
///   uE = c u + s v,  vN = c v - s u
void gridToGeo(const std::vector<double> & cosRot, const std::vector<double> & sinRot,
               int nz, std::vector<double> * u, std::vector<double> * v);

/// Geographic -> grid-relative at tracer points (in place):
///   u = c uE - s vN,  v = c vN + s uE
void geoToGrid(const std::vector<double> & cosRot, const std::vector<double> & sinRot,
               int nz, std::vector<double> * u, std::vector<double> * v);

}  // namespace mom6
}  // namespace ijedi
