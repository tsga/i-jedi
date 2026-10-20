#pragma once

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

    const Geometry & geom_;
    std::vector<std::string> readFilepaths_;  // ordered: ocn, ice, fix
    std::string writeFilepath_;               // ocean file for write
    bool nativeVelocity_ = false;             // write velocities on the C-grid
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi
