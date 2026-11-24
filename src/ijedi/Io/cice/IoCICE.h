#pragma once

#include <ostream>
#include <string>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/IoBase.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  class IoCICEParameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoCICEParameters, IoParametersBase)

   public:
    // Path prepended to all files
    oops::Parameter<std::string> datapath{"datapath",
                                          "path to location of files to be read",
                                          "./", this};

    // CICE restart file: read, and updated in place on write (the file must already exist:
    // staging it, e.g. copying the background restart, is left to the workflow)
    oops::RequiredParameter<std::string> restart_file{"restart_file",
                                                      "CICE restart file name",
                                                      this};
  };

  // -------------------------------------------------------------------------------------------------

  /// Io for CICE restarts on the MOM6 (tripolar) grid of the MOM6 geometry.
  ///
  /// The per-category CICE fields (see FieldsMetadataMarine.h) are mapped to the restart
  /// variables without configuration: sea_ice_category<k>_<quantity> is category k of the CICE
  /// variable(s) of that quantity, and its levels are the ice or snow layers, each a separate
  /// CICE variable (e.g. level l of sea_ice_category<k>_enthalpy is qice<lll>(k, :, :)); the
  /// levels are top-down, as in CICE (layer 1 is the top layer) and the MOM6 geometry. The
  /// numbers of categories and layers come from the MOM6 geometry. Other fields are mapped by
  /// "field io names" to 2-D (nj, ni) variables.
  ///
  /// A write updates the existing variables of an existing restart (at the geometry's mesh nodes
  /// only), so the variables, categories and grid points that JEDI does not hold keep their values.
  class IoCICE : public IoBase, private util::ObjectCounter<IoCICE>
  {
   public:
    static const std::string classname() { return "ijedi::IoCICE"; }

    typedef IoCICEParameters Parameters_;

    IoCICE(const Geometry &, const Parameters_ &);
    ~IoCICE() = default;
    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry & geom_;
    std::string filepath_;
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi
