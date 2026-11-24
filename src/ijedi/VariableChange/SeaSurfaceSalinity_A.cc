/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/VariableChange/SeaSurfaceSalinity.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "oops/util/for_each.h"
#include "oops/util/Logger.h"

namespace ijedi
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char SeaSurfaceSalinity_A::Name[] = "SeaSurfaceSalinity_A";
const oops::Variables SeaSurfaceSalinity_A::Ingredients{std::vector<std::string>{
                                                           "sea_water_salinity"}};

// Register the maker
static vader::RecipeMaker<SeaSurfaceSalinity_A> makerSSSaltToSSSalt_(
                    SeaSurfaceSalinity_A::Name);

SeaSurfaceSalinity_A::SeaSurfaceSalinity_A(const Parameters_ & params,
                             const vader::VaderConfigVars & configVariables):
                                            configVariables_{configVariables}
{
    oops::Log::trace() << "SeaSurfaceSalinity_A::SeaSurfaceSalinity_A(params)"
        << std::endl;
}

std::string SeaSurfaceSalinity_A::name() const
{
    return SeaSurfaceSalinity_A::Name;
}

oops::Variable SeaSurfaceSalinity_A::product() const
{
    return oops::Variable{"sea_surface_salinity"};
}

oops::Variables SeaSurfaceSalinity_A::ingredients() const
{
    return SeaSurfaceSalinity_A::Ingredients;
}

size_t SeaSurfaceSalinity_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return 1;
}

atlas::FunctionSpace SeaSurfaceSalinity_A::productFunctionSpace
                                                (const atlas::FieldSet & afieldset) const
{
    return afieldset.field("sea_water_salinity").functionspace();
}

// -------------------------------------------------------------------------------------------------

void SeaSurfaceSalinity_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering SeaSurfaceSalinity_A::executeNL function"
        << std::endl;

    // Get fields
    atlas::Field salinity = afieldset.field("sea_water_salinity");
    atlas::Field surface_salinity = afieldset.field("sea_surface_salinity");

    util::for_each_column(
        [&](
            const auto s,
            auto sss) {
                // For now, just use the first level of sea water salinity
                // as sea surface salinity
                sss(0) = s(0);
            },
        salinity,
        surface_salinity);

    oops::Log::trace() << "leaving SeaSurfaceSalinity_A::executeNL function" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void SeaSurfaceSalinity_A::executeTL(atlas::FieldSet & afieldsetTL,
                                        const atlas::FieldSet & afieldsetTraj) {
    oops::Log::trace() << "entering SeaSurfaceSalinity_A::executeTL function"
        << std::endl;

    // Get fields
    atlas::Field salinity_tl = afieldsetTL.field("sea_water_salinity");
    atlas::Field sss_tl = afieldsetTL.field("sea_surface_salinity");

    util::for_each_column(
        [&](
            const auto ds,
            auto dsss) {
                // For now, just use the first level of sea water salinity
                // as sea surface salinity
                dsss(0) = ds(0);
            },
        salinity_tl,
        sss_tl);

    oops::Log::trace() << "leaving SeaSurfaceSalinity_A::executeTL function" << std::endl;
}

// -------------------------------------------------------------------------------------------------


void SeaSurfaceSalinity_A::executeAD(atlas::FieldSet & afieldsetAD,
                                        const atlas::FieldSet & afieldsetTraj) {
    oops::Log::trace() << "entering SeaSurfaceSalinity_A::executeAD function"
        << std::endl;

    // Get fields
    atlas::Field salinity_ad = afieldsetAD.field("sea_water_salinity");
    atlas::Field sss_ad = afieldsetAD.field("sea_surface_salinity");

    util::for_each_column(
        [&](
            auto ds,
            auto dsss) {
                // For now, just use the first level of sea water salinity
                // as sea surface salinity
                ds(0) += dsss(0);
                dsss(0) = 0.;
            },
        salinity_ad,
        sss_ad);

    oops::Log::trace() << "leaving SeaSurfaceSalinity_A::executeAD function" << std::endl;
}

}  // namespace ijedi
