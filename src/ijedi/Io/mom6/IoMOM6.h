#pragma once

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/IoBase.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  class IoMOM6Parameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoMOM6Parameters, IoParametersBase)

   public:
    // Path prepended to all files
    oops::Parameter<std::string> datapath{"datapath",
                                          "path to location of files to be read",
                                          "./", this};

    // Ocean background file
    oops::RequiredParameter<std::string> ocn_file{"ocn_file",
                                                  "ocean background file name",
                                                  this};

    // Sea ice file
    oops::OptionalParameter<std::string> ice_file{"ice_file",
                                                  "sea ice file name",
                                                  this};

    // Fix file (decorrelation length scales, distance from coast, etc.)
    oops::OptionalParameter<std::string> fix_file{"fix_file",
                                                  "fix fields file name",
                                                  this};

    // Layout of written velocities: "tracer" (eastward/northward at h-points)
    // or "native" (grid-relative on the MOM6 C-grid, e.g. for an IAU increment).
    // Reading auto-detects the layout from the file dimensions.
    oops::Parameter<std::string> velocity_grid{"velocity grid",
                                               "tracer or native",
                                               "tracer", this};

    // Fields read from a file variable that holds a per-cell volume (thickness times
    // sea_ice_area_fraction), e.g. the CICE history hi_h and hs_h: field -> file variable.
    // On read, these fields are divided by sea_ice_area_fraction (0 where there is no ice),
    // which must be a variable of the same read. Writing is not affected.
    oops::OptionalParameter<std::map<std::string, std::string>> thickness_from_volume{
        "thickness from volume",
        "fields read as per-cell volumes and divided by sea_ice_area_fraction: "
        "field -> file variable", this};
  };

  // -------------------------------------------------------------------------------------------------

  class IoMOM6 : public IoBase, private util::ObjectCounter<IoMOM6>
  {
   public:
    static const std::string classname() { return "ijedi::IoMOM6"; }

    typedef IoMOM6Parameters Parameters_;

    IoMOM6(const Geometry &, const Parameters_ &);
    ~IoMOM6() = default;
    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;
    // Divide the "thickness from volume" fields of x by sea_ice_area_fraction
    void divideVolumesByIceFraction(atlas::FieldSet &) const;

    const Geometry & geom_;
    std::vector<std::string> readFilepaths_;  // ordered: ocn, ice, fix
    std::string writeFilepath_;               // ocean file for write
    bool nativeVelocity_ = false;             // write velocities on the C-grid
    std::map<std::string, std::string> thicknessFromVolume_;  // field -> file volume variable
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi
