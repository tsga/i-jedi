#include "ijedi/State/State.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>

#include "atlas/field.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Utilities/PrintHelper.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/for_each.h"

#include "ijedi/Io/IoBase.h"

namespace ijedi {

  // -----------------------------------------------------------------------------------------------

  State::State(const Geometry &geom, const eckit::Configuration &config)
      : mist::State(geom, oops::Variables(config, "state variables"),
                    util::DateTime(config.getString("date")), false),
        geom_(geom)
  {
    // If config has 'analytic init' then call analytic_init, else if config has 'io' then call read
    if (config.has("analytic init")) {
      analytic_init(config);
    } else if (config.has("io")) {
      read(config);
    } else {
      throw eckit::BadParameter("ijedi::State: config must have 'io' or 'analytic init'",
                                Here());
    }
    setAtlasFieldMetadata();
  }

  // -----------------------------------------------------------------------------------------------

  State::State(const Geometry &geom, const oops::Variables &vars, const util::DateTime &time,
               bool initToZero)
      : mist::State(geom, vars, time, initToZero), geom_(geom) {
    setAtlasFieldMetadata();
  }

  // -----------------------------------------------------------------------------------------------

  State::State(const Geometry &geom, const State &other)
      : mist::State(geom, other), geom_(geom) {
    setAtlasFieldMetadata();
  }

  // -----------------------------------------------------------------------------------------------

  State::State(const oops::Variables &vars, const State &other)
      : mist::State(vars, other), geom_(other.geom_) {
    setAtlasFieldMetadata();
  }

  State::State(const State &other) : mist::State(other), geom_(other.geom_) {
    setAtlasFieldMetadata();
  }

  // -----------------------------------------------------------------------------------------------

  State::~State() = default;

  // -----------------------------------------------------------------------------------------------

  State &State::operator=(const State &rhs)
  {
    mist::State::operator=(rhs);
    return *this;
  }

  // -----------------------------------------------------------------------------------------------

  // Required so LocalEnsembleDA links for the current non-inline LETKF path.
  // A real implementation is only needed for inline LETKF runs (`Run Inline: true`),
  // where forecast states must be redistributed onto the DA-local patch layout.
  void State::transpose(const State &, const eckit::mpi::Comm &, int ensNum, int transNum)
  {
    throw eckit::NotImplemented("ijedi::State::transpose is not implemented. "
                                "LETKF inline forecast transposition is unsupported in I-JEDI. "
                                "The current LETKF hookup supports the non-inline path only "
                                "(ensNum=" + std::to_string(ensNum)
                                + ", transNum=" + std::to_string(transNum) + ").", Here());
  }

  // -----------------------------------------------------------------------------------------------

  void State::read(const eckit::Configuration &config)
  {
    oops::Log::trace() << "ijedi::State::read starting" << std::endl;

    // Create a Parameters object
    StateParameters params;
    params.deserialize(swapIoMember(config));

    // Check that there are IO parameters
    if (params.io.value() == boost::none ||
        params.io.value()->ioParameters.value() == nullptr)
    {
      throw eckit::BadParameter("ijedi::State::read: No IO parameters provided", Here());
    }

    // Get the polymorphic IO parameters
    const IoParametersBase &ioParams = *params.io.value()->ioParameters.value();

    // Create the IO object to use
    // ---------------------------
    std::unique_ptr<IoBase> io(IoFactory::create(geom_, ioParams));

    // Call read method of child
    // -------------------------
    io->readBase(this->fieldSet());

    oops::Log::trace() << "ijedi::State::read done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  void State::write(const eckit::Configuration &config) const
  {
    oops::Log::trace() << "ijedi::State::write starting" << std::endl;

    // Create a Parameters object
    StateWriteParameters params;
    params.deserialize(swapIoMember(config));

    // Check that there are IO parameters
    if (params.io.value() == boost::none ||
        params.io.value()->ioParameters.value() == nullptr)
    {
      throw eckit::BadParameter("ijedi::State::write: No IO parameters provided", Here());
    }

    // Get the polymorphic IO parameters
    const IoParametersBase &ioParams = *params.io.value()->ioParameters.value();

    // Create the IO object to use
    // ---------------------------
    std::unique_ptr<IoBase> io(IoFactory::create(geom_, ioParams));

    // Call write method of child
    // --------------------------
    io->writeBase(this->fieldSet());

    oops::Log::trace() << "ijedi::State::write done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  void State::analytic_init(const eckit::Configuration &config)
  {
    oops::Log::trace() << "ijedi::State::analytic_init starting" << std::endl;
    // TODO(someone): implement analytic init
    oops::Log::trace() << "ijedi::State::analytic_init done" << std::endl;
  }

  void State::setAtlasFieldMetadata() {
    for (auto & field : this->fieldSet()) {
      field.metadata().set("interp_type", "default");
      // A temporary hack for interpolation masks for MOM6.
      // This should be replaced by using a proper mask field for different fields
      // (probably coming from FieldsMetaData)
      if (geom_.fields().has("mask2d")) {
        field.metadata().set("mask", "mask2d");
      }
    }
  }

  void State::print(std::ostream &os) const
  {
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
        const atlas::Field & field = fs.field(var.name());
        if (geom_.fields().has("owned")) {
          const atlas::Field owned = geom_.fields().field("owned");
          const auto[globalMin, globalMax, rms] = fieldMinMaxRMS(comm, field, &owned);
          os << std::endl
             << std::left << std::setw(maxNameLen) << var.name()
             << " : " << std::scientific << std::setprecision(10)
             << "Min=" << globalMin << ", Max=" << globalMax << ", RMS=" << rms;
        } else {
          const auto[globalMin, globalMax, rms] = fieldMinMaxRMS(comm, field);
          os << std::endl
             << std::left << std::setw(maxNameLen) << var.name()
             << " : " << std::scientific << std::setprecision(10)
             << "Min=" << globalMin << ", Max=" << globalMax << ", RMS=" << rms;
        }
    }
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi
