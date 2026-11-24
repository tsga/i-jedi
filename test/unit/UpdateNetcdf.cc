/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

// updateNetcdf on an MPAS-like restart: cells numbered 0..nCells-1 (as the atlas global index - 1),
// variables (Time, nCells, nVertLevels) with the levels innermost, a 2-D (Time, nCells) variable
// and one the update must not touch. The cells are distributed over the MPI tasks (PointCloud,
// equal_regions). The values are a function of the cell position, so a wrong mapping of the
// nodes to the cells fails the test.

#include <netcdf.h>

#include <cmath>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid.h"
#include "atlas/library/Library.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "ijedi/Io/utils/UpdateNetcdf.h"

namespace {

const int nCells = 2562;      // as the MPAS 480 km mesh
const int nVertLevels = 5;
const std::string fileName = "mpas_like_restart.nc";

// Cell centers: quasi-uniform points on the sphere (Fibonacci lattice)
std::vector<atlas::PointXY> cellCenters() {
  std::vector<atlas::PointXY> points;
  const double golden = M_PI * (3.0 - std::sqrt(5.0));
  for (int c = 0; c < nCells; ++c) {
    const double z = 1.0 - 2.0 * (c + 0.5) / nCells;
    const double lon = std::fmod(c * golden * 180.0 / M_PI, 360.0);
    points.emplace_back(lon, std::asin(z) * 180.0 / M_PI);
  }
  return points;
}

// Values of the file before the update, functions of the cell position
double theta(const double lon, const double lat, const int k) {
  return 300.0 + 10.0 * std::cos(lat * M_PI / 180.0) + 0.01 * lon - 5.0 * k;
}
double surfacePressure(const double lon, const double lat) {
  return 1.0e5 + 100.0 * std::sin(lat * M_PI / 180.0) + lon;
}
double untouched(const int c) { return 1.0 + c; }

void check(const int status) { ASSERT(status == NC_NOERR); }

// Rank 0 writes the file before the update
void createFile(const std::vector<atlas::PointXY> & centers) {
  int ncid, dimTime, dimCells, dimLev, vTheta, vPs, vKeep;
  check(nc_create(fileName.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &ncid));
  check(nc_def_dim(ncid, "Time", NC_UNLIMITED, &dimTime));
  check(nc_def_dim(ncid, "nCells", nCells, &dimCells));
  check(nc_def_dim(ncid, "nVertLevels", nVertLevels, &dimLev));
  const int dims3[3] = {dimTime, dimCells, dimLev};
  const int dims2[2] = {dimTime, dimCells};
  check(nc_def_var(ncid, "theta", NC_DOUBLE, 3, dims3, &vTheta));
  check(nc_def_var(ncid, "surface_pressure", NC_DOUBLE, 2, dims2, &vPs));
  check(nc_def_var(ncid, "untouched", NC_DOUBLE, 2, dims2, &vKeep));
  check(nc_enddef(ncid));
  std::vector<double> t(nCells * nVertLevels), ps(nCells), keep(nCells);
  for (int c = 0; c < nCells; ++c) {
    for (int k = 0; k < nVertLevels; ++k)
      t[c * nVertLevels + k] = theta(centers[c].x(), centers[c].y(), k);
    ps[c] = surfacePressure(centers[c].x(), centers[c].y());
    keep[c] = untouched(c);
  }
  const size_t start[3] = {0, 0, 0};
  const size_t count3[3] = {1, nCells, nVertLevels};
  const size_t count2[2] = {1, nCells};
  check(nc_put_vara_double(ncid, vTheta, start, count3, t.data()));
  check(nc_put_vara_double(ncid, vPs, start, count2, ps.data()));
  check(nc_put_vara_double(ncid, vKeep, start, count2, keep.data()));
  check(nc_close(ncid));
}

std::vector<double> readVariable(const std::string & name) {
  int ncid, varid, ndims, dimids[NC_MAX_VAR_DIMS];
  check(nc_open(fileName.c_str(), NC_NOWRITE, &ncid));
  check(nc_inq_varid(ncid, name.c_str(), &varid));
  check(nc_inq_varndims(ncid, varid, &ndims));
  check(nc_inq_vardimid(ncid, varid, dimids));
  size_t size = 1;
  for (int d = 0; d < ndims; ++d) {
    size_t len;
    check(nc_inq_dimlen(ncid, dimids[d], &len));
    size *= len;
  }
  std::vector<double> values(size);
  check(nc_get_var_double(ncid, varid, values.data()));
  check(nc_close(ncid));
  return values;
}

// The distributed MPAS-like cells and the FieldSet to write: 3 x the file values
struct Setup {
  std::vector<atlas::PointXY> centers = cellCenters();
  atlas::functionspace::PointCloud fs;
  atlas::FieldSet x;
  Setup()
    : fs(atlas::UnstructuredGrid(centers), atlas::grid::Partitioner("equal_regions")) {
    atlas::Field t = fs.createField<double>(atlas::option::name("air_potential_temperature")
                                            | atlas::option::levels(nVertLevels));
    atlas::Field ps = fs.createField<double>(atlas::option::name("air_pressure_at_surface")
                                             | atlas::option::levels(1));
    const auto lonlat = atlas::array::make_view<double, 2>(fs.lonlat());
    auto tv = atlas::array::make_view<double, 2>(t);
    auto pv = atlas::array::make_view<double, 2>(ps);
    for (atlas::idx_t n = 0; n < fs.size(); ++n) {
      for (int k = 0; k < nVertLevels; ++k) tv(n, k) = 3.0 * theta(lonlat(n, 0), lonlat(n, 1), k);
      pv(n, 0) = 3.0 * surfacePressure(lonlat(n, 0), lonlat(n, 1));
    }
    x.add(t);
    x.add(ps);
  }
};

CASE("updateNetcdf/mpas_like_restart") {
  const eckit::mpi::Comm & comm = eckit::mpi::comm();
  Setup s;
  if (comm.rank() == 0) createFile(s.centers);
  comm.barrier();

  // Levels innermost, a 2-D variable without nVertLevels in the same call
  ijedi::updateNetcdf(fileName, s.x, {"theta", "surface_pressure"},
                      ijedi::NetcdfLayout{{"nCells"}, {"nVertLevels"}}, comm);

  if (comm.rank() == 0) {
    const auto t = readVariable("theta");
    const auto ps = readVariable("surface_pressure");
    const auto keep = readVariable("untouched");
    double errT = 0.0, errPs = 0.0, errKeep = 0.0;
    for (int c = 0; c < nCells; ++c) {
      const double lon = s.centers[c].x(), lat = s.centers[c].y();
      for (int k = 0; k < nVertLevels; ++k)
        errT = std::max(errT, std::abs(t[c * nVertLevels + k] - 3.0 * theta(lon, lat, k)));
      errPs = std::max(errPs, std::abs(ps[c] - 3.0 * surfacePressure(lon, lat)));
      errKeep = std::max(errKeep, std::abs(keep[c] - untouched(c)));
    }
    eckit::Log::info() << "theta error " << errT << ", surface_pressure error " << errPs
                       << ", untouched error " << errKeep << std::endl;
    EXPECT(errT < 1.0e-9);
    EXPECT(errPs < 1.0e-6);
    EXPECT(errKeep == 0.0);
  }
}

CASE("updateNetcdf/layout_not_in_file_throws_on_all_tasks") {
  const eckit::mpi::Comm & comm = eckit::mpi::comm();
  Setup s;
  if (comm.rank() == 0) createFile(s.centers);
  comm.barrier();
  // Wrong level dimension name: rank 0 finds it, every rank throws (no hang in the gathers)
  EXPECT_THROWS_AS(ijedi::updateNetcdf(fileName, s.x, {"theta", "surface_pressure"},
                                       ijedi::NetcdfLayout{{"nCells"}, {"nLevels"}}, comm),
                   eckit::UserError);
}

}  // namespace

int main(int argc, char ** argv) {
  atlas::initialize(argc, argv);
  const int result = eckit::testing::run_tests(argc, argv, false);
  atlas::finalize();
  return result;
}
