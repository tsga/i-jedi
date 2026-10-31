#pragma once

#include <ostream>
#include <string>

#include "mist/base/State.h"
#include "oops/util/ObjectCounter.h"

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Io/IoBase.h"

namespace eckit
{
  class Configuration;
}  // namespace eckit

namespace oops
{
  class Variables;
}  // namespace oops

namespace util
{
  class DateTime;
}  // namespace util

namespace ijedi
{

  class Geometry;

  // -----------------------------------------------------------------------------------------------

  class AnalyticICParameters : public oops::Parameters
  {
    OOPS_CONCRETE_PARAMETERS(AnalyticICParameters, Parameters)
   public:
    // Analytic initial condition parameters
    oops::RequiredParameter<std::string> method{"method", this};
  };

  // -----------------------------------------------------------------------------------------------

  class StateParameters : public oops::Parameters
  {
    OOPS_CONCRETE_PARAMETERS(StateParameters, Parameters)
   public:
    oops::OptionalParameter<util::DateTime> datetime{"datetime", this};
    oops::OptionalParameter<oops::Variables> stateVariables{"state variables", this};
    // Analytic initial condition parameters
    oops::OptionalParameter<AnalyticICParameters> analytic{"analytic init", this};
    // Io parameters wrapper for polymorphic IO parameters (nested under "io" key)
    oops::OptionalParameter<IoParametersWrapper> io{"io", this};
  };

  // -----------------------------------------------------------------------------------------------

  class StateWriteParameters : public oops::Parameters
  {
    OOPS_CONCRETE_PARAMETERS(StateWriteParameters, Parameters)
   public:
    // Io parameters for writing (nested under "io" key)
    oops::OptionalParameter<IoParametersWrapper> io{"io", this};
  };

  // -----------------------------------------------------------------------------------------------

  class State : public mist::State, private util::ObjectCounter<State>
  {
   public:
    static std::string classname() { return "ijedi::State"; }

    State(const Geometry &, const eckit::Configuration &);
    State(const Geometry &, const oops::Variables &, const util::DateTime &,
          bool initToZero = true);
    State(const Geometry &, const State &);
    State(const oops::Variables &, const State &);
    State(const State &);

    ~State();

    State &operator=(const State &);

    void transpose(const State &, const eckit::mpi::Comm &, int ensNum, int transNum);

    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &) const;

    void setAtlasFieldMetadata();

   private:
    void analytic_init(const eckit::Configuration &);

    void print(std::ostream &os) const override;

    const Geometry &geom_;
  };

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi
