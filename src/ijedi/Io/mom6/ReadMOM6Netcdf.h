#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"

#include "eckit/mpi/Comm.h"

#include "ijedi/Io/mom6/VelocityStagger.h"

namespace ijedi {

class GeometryMOM6;

/// Read one or more MOM6-format NetCDF files into an Atlas FieldSet.
///
/// Rank 0 reads each file in order (ocean, sea-ice, fix, …).  For every file
/// the root maps the structured grid into Atlas global fields via gathered
/// global indices; a later file may overwrite fields from an earlier one, but
/// in practice each variable lives in exactly one file.  After all files are
/// processed a single collective scatter distributes data to all ranks,
/// followed by a halo exchange on ghost nodes.
///
/// @param filepaths    Ordered list of input NetCDF paths (rank 0 reads).
/// @param x            Distributed FieldSet to populate.
/// @param fileVarNames Per-field file variable name (empty => skip field).
/// @param scalings     Per-field scaling factor (0 => no scaling).
/// @param comm         MPI communicator.
/// @param velGrid      Grid rotation/mask info (root only). Required when a
///                     velocity pair is stored on the C-grid (lonq/latq dims):
///                     faces are averaged to h-points and rotated to
///                     eastward/northward.
void readMOM6Netcdf(const std::vector<std::string> & filepaths,
                    atlas::FieldSet & x,
                    const std::vector<std::string> & fileVarNames,
                    const std::vector<double> & scalings,
                    const eckit::mpi::Comm & comm,
                    const mom6::VelocityGridInfo * velGrid = nullptr);

/// Distributed variant: each task reads only its MOM6 block of each variable
/// (hyperslab) and the values are scattered to the JEDI points it owns, then
/// halo-exchanged.  No task holds a global field.  Every file must be on the
/// geometry's (effective) horizontal grid.
void readMOM6NetcdfBlocks(const std::vector<std::string> & filepaths,
                          atlas::FieldSet & x,
                          const std::vector<std::string> & fileVarNames,
                          const std::vector<double> & scalings,
                          const GeometryMOM6 & geom);

}  // namespace ijedi
