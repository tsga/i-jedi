#include <iostream>
#include <map>
#include <utility>

#include "eckit/exception/Exceptions.h"
#include "ijedi/FieldMetadata/FieldsMetadata.h"
#include "ijedi/FieldMetadata/FieldsMetadataDefault.h"
#include "ijedi/FieldMetadata/FieldsMetadataMarine.h"

// -------------------------------------------------------------------------------------------------

namespace ijedi
{

  // -----------------------------------------------------------------------------------------------

  FieldsMetadata::FieldsMetadata(const int nlev, const eckit::Configuration &modelData)
    : longNames_()
  {
    // Set the default metadata
    // ------------------------
    setMetadata(fieldsMetadata_, nlev);

    // Marine fields that depend on the geometry (CICE per-category sea ice)
    // ---------------------------------------------------------------------
    setMetadataMarine(fieldsMetadata_, nlev, modelData);

    // Create vector of the field long names
    // -------------------------------------
    for (const auto &field : fieldsMetadata_)
    {
      longNames_.push_back(field.first);
    }
  }

  // -----------------------------------------------------------------------------------------------

  FieldMetadata FieldsMetadata::getFieldMetadata(const std::string &longName) const
  {
    // Check that fieldsMetadata_ has longName in the keys and abort if not
    ASSERT_MSG(fieldsMetadata_.find(longName) != fieldsMetadata_.end(),
               "FieldMetadata error. Field \"" + longName + "\" not found in map. Ensure that " +
                   "the field is listed in FieldMetadataDefault.h");
    // Return Field Metadata
    return fieldsMetadata_.find(longName)->second;
  }

  // -----------------------------------------------------------------------------------------------

  std::unordered_map<std::string, size_t> FieldsMetadata::levelsPerVariable() const {
    std::unordered_map<std::string, size_t> levelsMap;
    for (const auto &field : fieldsMetadata_) {
      levelsMap[field.first] = field.second.getNumLevls();
    }
    return levelsMap;
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi

// -------------------------------------------------------------------------------------------------
