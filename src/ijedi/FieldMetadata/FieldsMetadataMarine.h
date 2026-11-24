#pragma once

#include <map>
#include <string>

#include "eckit/config/Configuration.h"

#include "ijedi/FieldMetadata/FieldsMetadata.h"
#include "ijedi/FieldMetadata/FieldsMetadataDefault.h"

namespace ijedi
{

  // -----------------------------------------------------------------------------------------------

  // CICE per-category sea-ice fields (CICE restarts, io "cice"). The level dimension is a physical
  // vertical level, so the category is part of the name and the levels are the thermodynamic
  // layers of ice or snow. With k = 1..number of categories:
  //
  //   sea_ice_category<k>_area_fraction           1 level            (CICE aicen)
  //   sea_ice_category<k>_volume                  1 level            (vicen)
  //   sea_ice_snow_category<k>_volume             1 level            (vsnon)
  //   sea_ice_category<k>_surface_temperature     1 level            (Tsfcn)
  //   sea_ice_category<k>_pond_area_fraction      1 level            (apnd)
  //   sea_ice_category<k>_pond_depth              1 level            (hpnd)
  //   sea_ice_category<k>_pond_lid_thickness      1 level            (ipnd)
  //   sea_ice_category<k>_enthalpy                sea ice layers     (qice001, qice002, ...)
  //   sea_ice_category<k>_ice_salinity            sea ice layers     (sice001, sice002, ...)
  //   sea_ice_snow_category<k>_enthalpy           snow layers        (qsno001, ...)
  //
  // The layers are top-down, as the MOM6 levels (levels_are_top_down): level 1 is CICE layer 1
  // (qice001, ...), the top layer of the ice (snow) and the last level is the bottom layer, as in
  // Icepack (icepack_therm_mushy/bl99: k = 1 top ice layer, k = nilyr bottom ice layer, snow
  // layers numbered from the top).
  //
  // The numbers of categories and layers are set by the geometry ("sea ice categories",
  // "sea ice layers", "snow layers", see the MOM6 geometry); the fields are only defined for a
  // geometry that provides them. Masked cells keep the values of the file (bctype none).
  void setMetadataMarine(std::map<std::string, FieldMetadata> &fieldsmetadata, const int nlev,
                         const eckit::Configuration &modelData)
  {
    if (!modelData.has("sea ice categories")) return;
    const int nCategories = modelData.getInt("sea ice categories");
    const std::string iceLayers = std::to_string(modelData.getInt("sea ice layers"));
    const std::string snowLayers = std::to_string(modelData.getInt("snow layers"));

    struct metadataStruct md;
    setMetadataStruct(md);

    for (int k = 1; k <= nCategories; ++k) {
      const std::string ice = "sea_ice_category" + std::to_string(k);
      const std::string snow = "sea_ice_snow_category" + std::to_string(k);

      md.longName = ice + "_area_fraction";
      md.units = "none";
      md.kind = "double";
      md.tracer = "false";
      md.levels = "1";
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);

      md.longName = ice + "_volume";
      md.units = "m";
      md.kind = "double";
      md.tracer = "false";
      md.levels = "1";
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);

      md.longName = snow + "_volume";
      md.units = "m";
      md.kind = "double";
      md.tracer = "false";
      md.levels = "1";
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);

      md.longName = ice + "_surface_temperature";
      md.units = "degC";
      md.kind = "double";
      md.tracer = "false";
      md.levels = "1";
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);

      md.longName = ice + "_pond_area_fraction";
      md.units = "none";
      md.kind = "double";
      md.tracer = "false";
      md.levels = "1";
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);

      md.longName = ice + "_pond_depth";
      md.units = "m";
      md.kind = "double";
      md.tracer = "false";
      md.levels = "1";
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);

      md.longName = ice + "_pond_lid_thickness";
      md.units = "m";
      md.kind = "double";
      md.tracer = "false";
      md.levels = "1";
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);

      md.longName = ice + "_enthalpy";
      md.units = "J m-3";
      md.kind = "double";
      md.tracer = "false";
      md.levels = iceLayers;
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);

      md.longName = ice + "_ice_salinity";
      md.units = "ppt";
      md.kind = "double";
      md.tracer = "false";
      md.levels = iceLayers;
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);

      md.longName = snow + "_enthalpy";
      md.units = "J m-3";
      md.kind = "double";
      md.tracer = "false";
      md.levels = snowLayers;
      md.vtype = "magnitude";
      addFieldMetadata(fieldsmetadata, nlev, md);
    }
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi
