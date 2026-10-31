#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "eckit/mpi/Comm.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "mist/base/Geometry.h"

#include "ijedi/Geometry/base/GeometryBase.h"

#include "ijedi/FieldMetadata/FieldsMetadata.h"

// Forward declarations
namespace eckit
{
  class Configuration;
}

namespace oops
{
  class Variables;
}

namespace ijedi
{

  // -----------------------------------------------------------------------------
  // Geometry handles geometry.

  class Geometry : public mist::Geometry, private util::ObjectCounter<Geometry>
  {
   public:
    static const std::string classname() { return "ijedi::Geometry"; }

    Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
    ~Geometry();

    // This might need to change with mist.
    const int &numLevels() const { return numberLevels_; }

    // Function to access field metadata
    const FieldsMetadata &getFieldMetadata() const { return *fieldsMeta_; }

    // Add ingredient fields required by some Vader recipes that are not state
    // variables. All sourcing is model-specific, via the per-model hook
    // GeometryBase::addModelVaderIngredients (default no-op).
    void addVaderIngredients(atlas::FieldSet &fset) const
      { geometryImpl_->addModelVaderIngredients(fields(), fset, numberLevels_); }

    // Model-specific implementation, for model IO that needs more than the
    // generic fieldset (e.g. the MOM6 structured decomposition).
    const GeometryBase &geometryImpl() const { return *geometryImpl_; }

   private:
    Geometry &operator=(const Geometry &);
    void print(std::ostream &) const;
    std::shared_ptr<FieldsMetadata> fieldsMeta_;
    std::shared_ptr<GeometryBase> geometryImpl_;
    int numberLevels_;
    std::string type_;
  };
  // -----------------------------------------------------------------------------

}  // namespace ijedi
