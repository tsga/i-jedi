
/*To use landVector seamlessly with JEDI algorithms (e.g., oops::Variational, oops::HofX), 
register it inside your model traits file (e.g., src/ijedi/TraitsLandVector.h):*/

#pragma once

#include "ijedi/Geometry/landVector/GeometryLandVector.h"
#include "ijedi/State/landVector/StateLandVector.h"
#include "ijedi/Increment/landVector/IncrementLandVector.h"

namespace ijedi {

struct LandVectorTraits {
  static std::string name() { return "LandVectorTraits"; }

  using Geometry  = landVector::GeometryLandVector;
  using State     = landVector::StateLandVector;
  using Increment = landVector::IncrementLandVector;
  // Add LinearVariableChange, Model, etc. as needed
};

}  // namespace ijedi