/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/VariableChange/VaderCookbook.h"

#include <map>
#include <string>
#include <vector>

namespace detail {
static std::map<std::string, std::vector<std::string>> cookbook() {
  return {
      // ps: from delp
      {"air_pressure_at_surface",    {"SurfaceAirPressure_A"}},
      // P: from delp, from ps (and ak/bk)
      {"air_pressure_levels",        {"AirPressureAtInterface_B",
                                      "AirPressureAtInterface_A"}},
      // p: from pe
      {"air_pressure",               {"AirPressure_A"}},
      // delp: from p
      {"air_pressure_thickness",     {"AirPressureThickness_A"}},
      // ln(p): from p
      {"ln_air_pressure",            {"LnAirPressure_A"}},
      // ln(pe): from pe
      {"ln_air_pressure_at_interface", {"LnAirPressureAtInterface_A"}},
      // p^kappa: from pe and ln(p)
      {"air_pressure_to_kappa",        {"AirPressureToKappa_A"}},
      // rh: from t q and qsat
      {"relative_humidity",            {"RelativeHumidity_A"}},
      // r: from q
      {"water_vapor_mixing_ratio_wrt_dry_air", {"WaterVaporMixingRatioWrtDryAir_C"}},
      // es: from t
      {"svp",                                  {"SaturationVaporPressure_A"}},
      // qsat: from t, and es
      {"saturation_water_vapor_mixing_ratio_wrt_moist_air", {"SaturationSpecificHumidity_A"}},
      // phi_e from phi tv ln(p) and ln(pe)
      {"geopotential_levels",                  {"GeopotentialLevels_A"}},
      // z: from phi
      {"geopotential_height",                  {"GeopotentialHeight_A"}},
      // z_surf: from phi_surf
      {"geopotential_height_at_surface",       {"GeopotentialHeightAtSurface_A"}},
      // ze: from phi_e
      {"geopotential_height_levels",           {"GeopotentialHeightLevels_A"}},
      // h_amsl: from z_surf
      {"height_above_mean_sea_level_at_surface",  {"HeightAboveMeanSeaLevelAtSurface_A"}},
      // pt: from t and pkz, from t and ps
      {"air_potential_temperature",               {"AirPotentialTemperature_B",
                                                   "AirPotentialTemperature_A"}},
      // tv: from t and q
      {"virtual_temperature",                     {"AirVirtualTemperature_A"}},
      // t: from tv and q; needed when tv is the control variable
      {"air_temperature",                         {"AirTemperature_B"}},
      // o3mr: from o3 mole fraction; needed when the control variable is ppmv
      // TODO(samueldegelia): re-enable once vader #421 is merged
      // {"ozone_mass_mixing_ratio",                 {"OzoneMassMixingRatio_A"}},
      // o3 mole fraction: from o3mr; lets the trajectory hold a variable the state does not
      {"mole_fraction_of_ozone_in_air",           {"MoleFractionOfOzoneInAir_A"}},
      // sst: from t
      {"sea_surface_temperature",                 {"SeaSurfaceTemperature_A"}},
      // t_insitu: NL from thetao (sea_water_potential_temperature), salinity, depth, lat/lon
      //           (SeaWaterTemperature_A); TL/AD via linearized TEOS-10 Jacobian
      //           (SeaWaterTemperature_B). Both are needed so the variable change
      //           works in the variational inner loop as well as the NL observer.
      // TODO(guillaumevernieres): re-enable once the vader insitu (GSW ocean) recipe PR is merged.
      // {"sea_water_temperature",                  {"SeaWaterTemperature_A",
      //                                             "SeaWaterTemperature_B"}},
  };
}
}  // namespace detail

namespace ijedi {

std::map<std::string, std::vector<std::string>> vaderDefaultCookbook() {
  return detail::cookbook();
}

}  // namespace ijedi
