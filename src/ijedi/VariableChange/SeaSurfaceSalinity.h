/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

class SeaSurfaceSalinity_AParameters : public vader::RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(SeaSurfaceSalinity_AParameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{"recipe name", this};
};

/*! \brief SeaSurfaceSalinity_A class defines a recipe for sea surface salinity
 *
 *  \details For now, just use the first level of sea water salinity
 */
class SeaSurfaceSalinity_A : public vader::RecipeBase {
 public:
    static const char Name[];
    static const oops::Variables Ingredients;

    typedef SeaSurfaceSalinity_AParameters Parameters_;

    SeaSurfaceSalinity_A(const Parameters_ &, const vader::VaderConfigVars &);

    // Recipe base class overrides
    bool hasTLAD() const override { return true; }
    std::string name() const override;
    oops::Variable product() const override;
    oops::Variables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    void executeNL(atlas::FieldSet &) override;
    void executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    void executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
    const vader::VaderConfigVars & configVariables_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi
