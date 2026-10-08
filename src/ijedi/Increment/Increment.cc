#include "ijedi/Increment/Increment.h"

#include <optional>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>

#include "atlas/field.h"
#include "atlas/util/Earth.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Utilities/PrintHelper.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/generic/ProximitySearch.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/for_each.h"

namespace ijedi {

  // -----------------------------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom, const oops::Variables & vars,
                       const util::DateTime & time)
      : mist::Increment(geom, vars, time), geom_(geom) {}

  // -----------------------------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom, const Increment & other, const bool ad)
      : mist::Increment(geom, other, ad), geom_(geom) {}

  // -----------------------------------------------------------------------------------------------

  Increment::Increment(const Increment & other, const bool copy)
      : mist::Increment(other, copy), geom_(other.geom_) {}

  // -----------------------------------------------------------------------------------------------

  Increment::~Increment() = default;

  // -----------------------------------------------------------------------------------------------

  Increment & Increment::operator=(const Increment & rhs) {
    mist::Increment::operator=(rhs);
    return *this;
  }

  // -----------------------------------------------------------------------------------------------

  void Increment::read(const eckit::Configuration & config) {
    oops::Log::trace() << "ijedi::Increment::read starting" << std::endl;

    // Create a Parameters object
    IncrementParameters params;
    params.deserialize(swapIoMember(config));

    // Check that there are IO parameters
    if (params.io.value() == boost::none ||
        params.io.value()->ioParameters.value() == nullptr)
    {
      throw eckit::BadParameter("ijedi::Increment::read: No IO parameters provided", Here());
    }

    // Get the polymorphic IO parameters
    const IoParametersBase &ioParams = *params.io.value()->ioParameters.value();

    // Create the IO object to use
    // ---------------------------
    std::unique_ptr<IoBase> io(IoFactory::create(geom_, ioParams));

    // Call read method of child
    // -------------------------
    io->readBase(this->fieldSet());

    oops::Log::trace() << "ijedi::Increment::read done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  void Increment::write(const eckit::Configuration & config) const {
    oops::Log::trace() << "ijedi::Increment::write starting" << std::endl;

    // Create a Parameters object
    IncrementWriteParameters params;
    params.deserialize(swapIoMember(config));

    // Check that there are IO parameters
    if (params.io.value() == boost::none ||
        params.io.value()->ioParameters.value() == nullptr)
    {
      throw eckit::BadParameter("ijedi::Increment::write: No IO parameters provided", Here());
    }

    // Get the polymorphic IO parameters
    const IoParametersBase &ioParams = *params.io.value()->ioParameters.value();

    // Create the IO object to use
    // ---------------------------
    std::unique_ptr<IoBase> io(IoFactory::create(geom_, ioParams));

    // Call write method of child
    // --------------------------
    io->writeBase(this->fieldSet());

    oops::Log::trace() << "ijedi::Increment::write done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  void Increment::dirac(const eckit::Configuration & config) {
    oops::Log::trace() << "ijedi::Increment::dirac starting" << std::endl;

    // Create a Parameters object
    DiracParameters params;
    params.deserialize(config);

    // Grid-agnostic lon/lat dirac: the globally-nearest owned grid node is found
    // via the geometry's shared KD-tree (the same tree used by interpolation).
    // Validate that all vectors are the same length.
    const std::vector<double> & lon = params.lon.value();
    const std::vector<double> & lat = params.lat.value();
    const std::vector<int> & level = params.level.value();
    const std::vector<std::string> & vars = params.variable.value();
    const size_t nDiracs = lon.size();
    ASSERT_MSG(lat.size() == nDiracs, "Dirac: 'lat' inconsistent length");
    ASSERT_MSG(level.size() == nDiracs, "Dirac: 'level' inconsistent length");
    ASSERT_MSG(vars.size() == nDiracs, "Dirac: 'variable' inconsistent length");

    const auto & comm = this->geom_.comm();
    const auto & geomData = this->geom_.geometryData();
    const oops::ProximitySearch & proximitySearch =
        oops::getProximitySearch(geomData.functionSpace(), geomData.comm());
    const auto lonlatView = atlas::array::make_view<double, 2>(geomData.functionSpace().lonlat());

    // Search radius for the nearest owned node. The global tree returns the
    // closest point within this chord distance (meters); a value larger than
    // any grid spacing guarantees a hit. Use a quarter of the Earth's
    // circumference so even the coarsest grids resolve.
    const double searchRadius = 0.25 * 2.0 * M_PI * atlas::util::Earth::radius();

    // Start from zero so only the requested points are nonzero.
    this->zero();

    for (size_t jdir = 0; jdir < nDiracs; ++jdir) {
      atlas::Field field = this->fieldSet().field(vars[jdir]);

      // The MPI task owning the globally-nearest node, then (on that task) the
      // task-local functionspace index of that node.
      const int localTask = proximitySearch.taskOwningClosestPoint(lat[jdir], lon[jdir]);

      // level input is 1-based -> 0-based array index.
      const int lev = level[jdir] - 1;
      double lonDir = 0.0;
      double latDir = 0.0;
      if (static_cast<size_t>(localTask) == comm.rank()) {
        const std::optional<int> index =
            proximitySearch.closestPointWithinRadius(lat[jdir], lon[jdir], searchRadius);
        ASSERT_MSG(index.has_value(), "Dirac: no owned grid node found near requested point");
        auto view = atlas::array::make_view<double, 2>(field);
        view(*index, lev) = 1.0;
        lonDir = lonlatView(*index, 0);
        latDir = lonlatView(*index, 1);
      }

      // Log the resolved location (sum reduction picks up the owning task's value).
      comm.allReduceInPlace(lonDir, eckit::mpi::sum());
      comm.allReduceInPlace(latDir, eckit::mpi::sum());
      oops::Log::info() << "ijedi::Increment::dirac point #" << jdir << " (" << vars[jdir]
                        << "): requested " << lon[jdir] << "/" << lat[jdir]
                        << ", placed at " << lonDir << "/" << latDir
                        << ", level " << level[jdir] << std::endl;
    }

    oops::Log::trace() << "ijedi::Increment::dirac done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  void Increment::print(std::ostream & os) const {
    os << std::endl
       << "  Valid time: " << this->validTime()
       << ", nFields = " << this->variables().size();

    const auto & comm = geom_.comm();
    const auto & fs   = this->fieldSet();
    size_t maxNameLen = 0;
    for (const auto & var : this->variables()) {
      maxNameLen = std::max(maxNameLen, var.name().size());
    }
    for (const auto & var : this->variables()) {
      const atlas::Field & field            = fs.field(var.name());
      const auto[globalMin, globalMax, rms] = fieldMinMaxRMS(comm, field);
      os << std::endl
         << std::left << std::setw(maxNameLen) << var.name()
         << " : " << std::scientific << std::setprecision(10)
         << "Min=" << globalMin << ", Max=" << globalMax << ", RMS=" << rms;
    }
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi
