#include <ostream>
#include <string>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/fv3-restart/IoFV3Restart.h"
#include "ijedi/Io/fv3-restart/IoFV3Restart.interface.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoFV3Restart> makerIoFV3Restart_("fv3 restart");
  // -------------------------------------------------------------------------------------------------
  namespace {
  eckit::LocalConfiguration makeRuntimeConfig(const IoFV3Restart::Parameters_ &params,
                                              const Geometry &geom,
                                              const atlas::FieldSet &x)
  {
    eckit::LocalConfiguration runtime(params.toConfiguration());

    std::vector<std::string> activeFields;
    std::vector<std::string> tracerFields;
    activeFields.reserve(x.size());

    for (const auto &field : x)
    {
      const std::string fieldName = field.name();
      activeFields.push_back(fieldName);
      if (geom.getFieldMetadata().getFieldMetadata(fieldName).getIsTracer())
      {
        tracerFields.push_back(fieldName);
      }
    }

    runtime.set("active fields", activeFields);
    runtime.set("tracer fields", tracerFields);
    return runtime;
  }
  }  // namespace
  // -------------------------------------------------------------------------------------------------
  IoFV3Restart::IoFV3Restart(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()),
        geom_(geom),
        parameters_(params)
  {
    util::Timer timer(classname(), "IoFV3Restart");
    oops::Log::trace() << classname() << " constructor starting" << std::endl;
    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  IoFV3Restart::~IoFV3Restart()
  {
    util::Timer timer(classname(), "~IoFV3Restart");
    oops::Log::trace() << classname() << " destructor starting" << std::endl;
    oops::Log::trace() << classname() << " destructor done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoFV3Restart::read(atlas::FieldSet &x,
                          const eckit::LocalConfiguration &fileionames,
                          const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read state");
    oops::Log::trace() << classname() << " read state starting" << std::endl;

    const eckit::LocalConfiguration runtimeConfig = makeRuntimeConfig(parameters_, geom_, x);
    ijedi_io_fv3_restart_read_f90(runtimeConfig, geom_.modelData(), x.get(),
                                  fileionames, fileioscaling);

    // Halo-exchange the fields the Fortran reader just filled, as the MOM6 reader already does.
    // Needed because the reader populates owned points only, so the halo keeps its allocated
    // value and consumers with cross-subdomain stencils (the SABER interpolation blocks) read it.
    for (auto &field : x) {
      if (field.functionspace()) {
        field.set_dirty();
        field.haloExchange();
      }
    }

    oops::Log::trace() << classname() << " read state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoFV3Restart::write(const atlas::FieldSet &x,
                           const eckit::LocalConfiguration &fileionames,
                           const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;

    const eckit::LocalConfiguration runtimeConfig = makeRuntimeConfig(parameters_, geom_, x);
    ijedi_io_fv3_restart_write_f90(runtimeConfig, geom_.modelData(), x.get(),
                                   fileionames, fileioscaling);

    oops::Log::trace() << classname() << " write state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoFV3Restart::print(std::ostream &os) const
  {
    os << classname() << " Io for FMS Restarts";
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi
