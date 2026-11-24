#include "ijedi/Io/cice/IoCICE.h"

#include <netcdf.h>

#include <algorithm>
#include <cstdio>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Geometry/mom6/GeometryMOM6.h"
#include "ijedi/Io/utils/UpdateNetcdf.h"
#include "ijedi/Io/mom6/FillMaskedCells.h"
#include "ijedi/Io/mom6/ReadMOM6Netcdf.h"

namespace ijedi
{
  namespace
  {
  // -------------------------------------------------------------------------------------------------
  // CICE restart variables of the per-category fields (FieldsMetadataMarine.h): field
  // <prefix><k><suffix> is category k of <cice>, or, for a layered quantity, has one level per
  // layer l, category k of <cice><lll>. Level l-1 is CICE layer l, so the levels are top-down
  // (CICE layer 1 is the top layer of the ice or snow), as the MOM6 geometry's levels
  enum class Layers { none, ice, snow };
  struct CiceQuantity {
    std::string prefix;
    std::string suffix;
    std::string cice;
    Layers layers;
  };
  const std::vector<CiceQuantity> & ciceQuantities() {
    static const std::vector<CiceQuantity> quantities = {
      {"sea_ice_category", "_area_fraction", "aicen", Layers::none},
      {"sea_ice_category", "_volume", "vicen", Layers::none},
      {"sea_ice_snow_category", "_volume", "vsnon", Layers::none},
      {"sea_ice_category", "_surface_temperature", "Tsfcn", Layers::none},
      {"sea_ice_category", "_pond_area_fraction", "apnd", Layers::none},
      {"sea_ice_category", "_pond_depth", "hpnd", Layers::none},
      {"sea_ice_category", "_pond_lid_thickness", "ipnd", Layers::none},
      {"sea_ice_category", "_enthalpy", "qice", Layers::ice},
      {"sea_ice_category", "_ice_salinity", "sice", Layers::ice},
      {"sea_ice_snow_category", "_enthalpy", "qsno", Layers::snow}};
    return quantities;
  }

  // A per-category field: its category (1-based) and the CICE variable of each of its levels
  struct CategoryField {
    std::string name;
    int category;
    std::vector<std::string> ciceNames;
  };

  // How the fields of a FieldSet map to the restart variables
  struct CiceLayout {
    std::vector<CategoryField> categoryFields;
    std::vector<std::string> ciceVariables;   // CICE variables of the categoryFields (unique)
    atlas::FieldSet otherFields;              // mapped by "field io names"
    std::vector<std::string> otherNames;
    std::vector<double> otherScalings;
  };

  // Parses sea_ice_[snow_]category<k>_<quantity>; false if the name is not a per-category field
  bool parseCategoryField(const std::string & name, const int nIceLayers, const int nSnowLayers,
                          CategoryField & cf) {
    for (const auto & q : ciceQuantities()) {
      if (name.size() <= q.prefix.size() + q.suffix.size()) continue;
      if (name.compare(0, q.prefix.size(), q.prefix) != 0) continue;
      if (name.compare(name.size() - q.suffix.size(), q.suffix.size(), q.suffix) != 0) continue;
      const std::string k = name.substr(q.prefix.size(),
                                        name.size() - q.prefix.size() - q.suffix.size());
      if (k.find_first_not_of("0123456789") != std::string::npos) continue;
      cf.name = name;
      cf.category = std::stoi(k);
      cf.ciceNames.clear();
      const int nLayers = q.layers == Layers::ice ? nIceLayers
                        : q.layers == Layers::snow ? nSnowLayers : 0;
      if (nLayers == 0) {
        cf.ciceNames.push_back(q.cice);
      } else {
        for (int l = 1; l <= nLayers; ++l) {
          char lll[8];
          std::snprintf(lll, sizeof(lll), "%03d", l);
          cf.ciceNames.push_back(q.cice + lll);
        }
      }
      return true;
    }
    return false;
  }

  CiceLayout ciceLayout(const atlas::FieldSet & x, const eckit::LocalConfiguration & fileionames,
                        const eckit::LocalConfiguration & fileioscaling,
                        const GeometryMOM6 & geom) {
    CiceLayout layout;
    for (const auto & field : x) {
      CategoryField cf;
      if (!parseCategoryField(field.name(), geom.nSeaIceLayers(), geom.nSnowLayers(), cf)) {
        layout.otherFields.add(field);
        layout.otherNames.push_back(fileionames.has(field.name())
                                    ? fileionames.getString(field.name()) : "");
        layout.otherScalings.push_back(fileioscaling.has(field.name())
                                       ? fileioscaling.getDouble(field.name()) : 0.0);
        continue;
      }
      if (cf.category < 1 || cf.category > geom.nSeaIceCategories()) {
        throw eckit::UserError("IoCICE: " + cf.name + ": category out of range (the geometry "
                               "has " + std::to_string(geom.nSeaIceCategories())
                               + " sea ice categories)", Here());
      }
      if (field.shape(1) != static_cast<int>(cf.ciceNames.size())) {
        throw eckit::UserError("IoCICE: " + cf.name + " has " + std::to_string(field.shape(1))
                               + " levels, expected " + std::to_string(cf.ciceNames.size()),
                               Here());
      }
      for (const auto & v : cf.ciceNames) {
        if (std::find(layout.ciceVariables.begin(), layout.ciceVariables.end(), v)
            == layout.ciceVariables.end()) layout.ciceVariables.push_back(v);
      }
      layout.categoryFields.push_back(cf);
    }
    return layout;
  }

  // Rank 0 checks that the restart has the geometry's number of categories and every variable;
  // all ranks throw together on a failure
  void checkCiceRestart(const std::string & filepath, const std::vector<std::string> & variables,
                        const int nCategories, const eckit::mpi::Comm & comm) {
    std::string error;
    if (comm.rank() == 0) {
      int ncid = -1, dimid = -1, varid = -1;
      size_t ncat = 0;
      if (nc_open(filepath.c_str(), NC_NOWRITE, &ncid) != NC_NOERR) {
        error = "cannot open " + filepath;
      } else {
        if (nc_inq_dimid(ncid, "ncat", &dimid) != NC_NOERR
            || nc_inq_dimlen(ncid, dimid, &ncat) != NC_NOERR) {
          error = "no ncat dimension in " + filepath;
        } else if (static_cast<int>(ncat) != nCategories) {
          error = filepath + " has " + std::to_string(ncat) + " categories, the geometry "
                  + std::to_string(nCategories);
        }
        for (const auto & v : variables) {
          if (!error.empty()) break;
          if (nc_inq_varid(ncid, v.c_str(), &varid) != NC_NOERR)
            error = "no variable " + v + " in " + filepath;
        }
        nc_close(ncid);
      }
    }
    int failed = error.empty() ? 0 : 1;
    comm.broadcast(failed, 0);
    if (failed) {
      throw eckit::UserError("IoCICE: " + (comm.rank() == 0 ? error
                             : std::string("restart check failed on rank 0")), Here());
    }
  }

  // The CICE variables, read from the restart: one field per variable, one level per category
  atlas::FieldSet readCiceVariables(const std::string & filepath,
                                    const std::vector<std::string> & variables,
                                    const atlas::FunctionSpace & functionSpace,
                                    const GeometryMOM6 & geom) {
    checkCiceRestart(filepath, variables, geom.nSeaIceCategories(), geom.comm());
    const atlas::functionspace::NodeColumns fs(functionSpace);
    atlas::FieldSet ciceFs;
    for (const auto & v : variables) {
      ciceFs.add(fs.createField<double>(atlas::option::name(v)
                                        | atlas::option::levels(geom.nSeaIceCategories())));
    }
    readMOM6NetcdfBlocks({filepath}, ciceFs, variables,
                         std::vector<double>(variables.size(), 0.0), geom, {"ncat"});
    return ciceFs;
  }

  const GeometryMOM6 & mom6Geometry(const Geometry & geom, const std::string & classname) {
    const auto * mom6Geom = dynamic_cast<const GeometryMOM6 *>(&geom.geometryImpl());
    if (mom6Geom == nullptr)
      throw eckit::UserError(classname + ": the cice io needs a MOM6 geometry", Here());
    return *mom6Geom;
  }
  // -------------------------------------------------------------------------------------------------
  }  // namespace

  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoCICE> makerIoCICE_("cice");
  // -------------------------------------------------------------------------------------------------
  IoCICE::IoCICE(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()),
        geom_(geom)
  {
    oops::Log::trace() << classname() << " constructor starting" << std::endl;
    filepath_ = params.datapath.value() + "/" + params.restart_file.value();
    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoCICE::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                    const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read state");
    oops::Log::trace() << classname() << " read state starting" << std::endl;

    const GeometryMOM6 & mom6Geom = mom6Geometry(geom_, classname());
    CiceLayout layout = ciceLayout(x, fileionames, fileioscaling, mom6Geom);

    // Per-category fields: each task reads its MOM6 block of the CICE variables (one level per
    // category), then level l of a field is category k of the variable of layer l
    if (!layout.categoryFields.empty()) {
      const atlas::FieldSet ciceFs = readCiceVariables(filepath_, layout.ciceVariables,
                                                       x.field(0).functionspace(), mom6Geom);
      for (const auto & cf : layout.categoryFields) {
        auto view = atlas::array::make_view<double, 2>(x.field(cf.name));
        const double scale = fileioscaling.has(cf.name) ? fileioscaling.getDouble(cf.name) : 0.0;
        for (size_t l = 0; l < cf.ciceNames.size(); ++l) {
          const auto ciceView = atlas::array::make_view<const double, 2>(
            ciceFs.field(cf.ciceNames[l]));
          for (atlas::idx_t n = 0; n < view.shape(0); ++n) {
            view(n, l) = ciceView(n, cf.category - 1);
            if (scale != 0.0) view(n, l) *= scale;
          }
        }
      }
    }

    // Other fields: one variable each
    if (layout.otherFields.size() > 0) {
      readMOM6NetcdfBlocks({filepath_}, layout.otherFields, layout.otherNames,
                           layout.otherScalings, mom6Geom, {"ncat"});
    }

    // Boundary conditions of masked cells, as for the MOM6 io (the per-category CICE restart
    // fields have none: their masked cells keep the values of the file)
    applyBoundaryConditions(x, geom_.fields().field("mask3d"), geom_.getFieldMetadata());

    oops::Log::trace() << classname() << " read state done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoCICE::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                     const eckit::LocalConfiguration &) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;

    // The restart is updated in place: it must already exist (staged by the workflow)
    const GeometryMOM6 & mom6Geom = mom6Geometry(geom_, classname());
    CiceLayout layout = ciceLayout(x, fileionames, eckit::LocalConfiguration(), mom6Geom);

    atlas::FieldSet out;
    std::vector<std::string> outNames;

    // Per-category fields: the CICE variables are read first, so that the categories JEDI does
    // not hold keep their values, then category k of the variable of layer l is overwritten by
    // level l of the field
    if (!layout.categoryFields.empty()) {
      atlas::FieldSet ciceFs = readCiceVariables(filepath_, layout.ciceVariables,
                                                 x.field(0).functionspace(), mom6Geom);
      for (const auto & cf : layout.categoryFields) {
        const auto view = atlas::array::make_view<const double, 2>(x.field(cf.name));
        for (size_t l = 0; l < cf.ciceNames.size(); ++l) {
          auto ciceView = atlas::array::make_view<double, 2>(ciceFs.field(cf.ciceNames[l]));
          for (atlas::idx_t n = 0; n < view.shape(0); ++n) {
            ciceView(n, cf.category - 1) = view(n, l);
          }
        }
      }
      for (const auto & v : layout.ciceVariables) {
        out.add(ciceFs.field(v));
        outNames.push_back(v);
      }
    }

    // Other fields: one variable each
    for (int i = 0; i < layout.otherFields.size(); ++i) {
      out.add(layout.otherFields.field(i));
      outNames.push_back(layout.otherNames[i]);
    }

    // In-place update of the restart: horizontal position (nj, ni), levels on the category
    // dimension for the per-category variables, none for the 2-D ones
    if (out.size() > 0) {
      updateNetcdf(filepath_, out, outNames, NetcdfLayout{{"nj", "ni"}, {"ncat"}},
                   geom_.getComm());
    }

    oops::Log::trace() << classname() << " write state done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoCICE::print(std::ostream &os) const
  {
    os << classname() << " Io for CICE restarts on the MOM6 grid: " << filepath_;
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi
