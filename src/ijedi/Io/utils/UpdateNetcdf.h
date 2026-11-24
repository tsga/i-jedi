#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"

#include "eckit/mpi/Comm.h"

namespace ijedi {

/// How the nodes and levels of a FieldSet map to the dimensions of a NetCDF variable.
struct NetcdfLayout {
  /// File dimensions of the horizontal position, outer to inner: the atlas global index - 1 of a
  /// node is its C-order flat index over them, e.g. {"nj", "ni"} (MOM6, CICE),
  /// {"tile", "yaxis_1", "xaxis_1"} (FV3 cubed sphere) or {"nCells"} (MPAS).
  std::vector<std::string> horizontalDims;
  /// File dimensions holding the levels of a field, outer to inner (their sizes multiply to the
  /// number of levels); empty for single-level fields, e.g. {"ncat"}, {"zaxis_1"},
  /// {"nVertLevels"}.
  std::vector<std::string> levelDims;
};

/// Update, in place, the variables of an existing NetCDF file from a distributed FieldSet.
///
/// The dimensions of a variable may come in any order (e.g. (Time, nCells, nVertLevels) or
/// (time, tile, z, y, x)); its dimensions that are not in the layout (e.g. a time dimension) must
/// have length 1. Each field in turn is gathered to rank 0, which reads the variable, overwrites
/// it at the gathered nodes and writes it back: the grid points without a node, the variables
/// not in the FieldSet and the other attributes of the file keep their values. Works with any
/// atlas function space (generic gather and global index), and only rank 0 does i/o, so classic
/// and 64-bit offset files can be updated as well as NetCDF-4 ones.
///
/// The file, the variables and their shapes are checked on rank 0 before any collective call,
/// and all ranks throw together on a failure.
///
/// @param filepath      Existing NetCDF file.
/// @param x             Distributed FieldSet (all ranks participate in the gathers).
/// @param fileVarNames  File variable name of each field of x (empty => the field is skipped).
/// @param layout        Horizontal and level dimensions of the variables.
/// @param comm          MPI communicator.
void updateNetcdf(const std::string & filepath,
                  const atlas::FieldSet & x,
                  const std::vector<std::string> & fileVarNames,
                  const NetcdfLayout & layout,
                  const eckit::mpi::Comm & comm);

}  // namespace ijedi
