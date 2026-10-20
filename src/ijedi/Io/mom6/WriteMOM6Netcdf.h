#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"

#include "eckit/mpi/Comm.h"

#include "ijedi/Io/mom6/VelocityStagger.h"

namespace ijedi {

/// Write an Atlas FieldSet to a MOM6-format NetCDF file.
///
/// Each field in @p x is gathered to rank 0 using Atlas NodeColumns::gather(),
/// unpacked into a structured (Time, Layer/lath, lonh) buffer, and written to
/// the file at @p filepath.
///
/// @param filepath      Output NetCDF path (created on rank 0).
/// @param x             Distributed FieldSet (all ranks participate in gather).
/// @param fileVarNames  Per-field file variable name (empty ⇒ skip).
/// @param ni            Global i-dimension (lonh).
/// @param nj            Global j-dimension (lath).
/// @param nz            Number of vertical levels (Layer).
/// @param comm          MPI communicator.
/// @param nativeVelocity If true, velocity pairs are rotated to grid-relative
///                       and interpolated to the MOM6 C-grid (u on lonq, v on
///                       latq); otherwise written eastward/northward on h-points.
/// @param velGrid        Grid rotation/mask info (root only), required when
///                       nativeVelocity is set and a velocity pair is written.
void writeMOM6Netcdf(const std::string & filepath,
                     const atlas::FieldSet & x,
                     const std::vector<std::string> & fileVarNames,
                     int ni, int nj, int nz,
                     const eckit::mpi::Comm & comm,
                     bool nativeVelocity = false,
                     const mom6::VelocityGridInfo * velGrid = nullptr);

}  // namespace ijedi
