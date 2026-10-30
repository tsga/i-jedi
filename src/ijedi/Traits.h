#pragma once

#include <string>

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/LinearModel/LinearModel.h"
#include "ijedi/LinearVariableChange/LinearVariableChange.h"
#include "ijedi/Model/Model.h"
#include "ijedi/State/State.h"
#include "ijedi/VariableChange/VariableChange.h"

#include "mist/base/ErrorCovariance.h"
#include "mist/base/ModelAuxControl.h"
#include "mist/base/ModelAuxCovariance.h"
#include "mist/base/ModelAuxIncrement.h"
#include "mist/base/ModelData.h"

#include "oops/generic/UnstructuredInterpolator.h"
#include "oops/generic/AtlasInterpolator.h"
#include "ijedi/Geometry/landVector/LocalInterpolatorLandVector.h"

namespace ufo {
template<typename ITERATOR> class ObsLocalization;
}

namespace ijedi
{

  struct Traits
  {
    static std::string name() { return "ijedi"; }
    static std::string nameCovar() { return "ijediError"; }

    typedef ijedi::Geometry                  Geometry;
    typedef mist::GeometryIterator           GeometryIterator;
    typedef ijedi::State                     State;
    typedef ijedi::Increment                 Increment;
    typedef mist::ModelData                  ModelData;
    typedef ijedi::VariableChange            VariableChange;
    typedef ijedi::LinearVariableChange      LinearVariableChange;
    typedef ijedi::Model                     Model;
    typedef ijedi::LinearModel               LinearModel;
    typedef mist::ErrorCovariance            Covariance;
    typedef mist::ModelAuxControl            ModelAuxControl;
    typedef mist::ModelAuxIncrement          ModelAuxIncrement;
    typedef mist::ModelAuxCovariance         ModelAuxCovariance;
    //typedef oops::UnstructuredInterpolator   LocalInterpolator;
    // Bind custom KNN local interpolator here!
    //typedef ijedi::LocalInterpolatorLandVector LocalInterpolator;
    typedef oops::AtlasInterpolator          LocalInterpolator;
    typedef ufo::ObsLocalization<mist::GeometryIterator> ObsLocalization;
  };

  struct TraitsAtm
  {
    // The name here is not meaningful because ijedi supports both atmosphere
    // and ocean models. This name is used in yamls for coupled applications, and
    // "atmosphere" and "ocean" are chosen solely for user-friendliness
    static std::string name() { return "atmosphere"; }
    static std::string nameCovar() { return "ijediError"; }

    typedef ijedi::Geometry                  Geometry;
    typedef mist::GeometryIterator           GeometryIterator;
    typedef ijedi::State                     State;
    typedef ijedi::Increment                 Increment;
    typedef mist::ModelData                  ModelData;
    typedef ijedi::VariableChange            VariableChange;
    typedef ijedi::LinearVariableChange      LinearVariableChange;
    typedef ijedi::Model                     Model;
    typedef ijedi::LinearModel               LinearModel;
    typedef mist::ErrorCovariance            Covariance;
    typedef mist::ModelAuxControl            ModelAuxControl;
    typedef mist::ModelAuxIncrement          ModelAuxIncrement;
    typedef mist::ModelAuxCovariance         ModelAuxCovariance;
    //typedef oops::UnstructuredInterpolator   LocalInterpolator;
    //typedef ijedi::LocalInterpolatorLandVector LocalInterpolator;
    typedef oops::AtlasInterpolator          LocalInterpolator;
    typedef ufo::ObsLocalization<mist::GeometryIterator> ObsLocalization;
  };

  struct TraitsOcn
  {
    static std::string name() { return "ocean"; }
    static std::string nameCovar() { return "ijediError"; }

    typedef ijedi::Geometry                  Geometry;
    typedef mist::GeometryIterator           GeometryIterator;
    typedef ijedi::State                     State;
    typedef ijedi::Increment                 Increment;
    typedef mist::ModelData                  ModelData;
    typedef ijedi::VariableChange            VariableChange;
    typedef ijedi::LinearVariableChange      LinearVariableChange;
    typedef ijedi::Model                     Model;
    typedef ijedi::LinearModel               LinearModel;
    typedef mist::ErrorCovariance            Covariance;
    typedef mist::ModelAuxControl            ModelAuxControl;
    typedef mist::ModelAuxIncrement          ModelAuxIncrement;
    typedef mist::ModelAuxCovariance         ModelAuxCovariance;
    //typedef oops::UnstructuredInterpolator   LocalInterpolator;
    typedef oops::AtlasInterpolator          LocalInterpolator;
    //typedef ijedi::LocalInterpolatorLandVector LocalInterpolator;
    typedef ufo::ObsLocalization<mist::GeometryIterator> ObsLocalization;
  };

}  // namespace ijedi
