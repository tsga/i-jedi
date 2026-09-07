/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include <unordered_map>
#include <vector>

#include "atlas/field.h"

#include "oops/generic/GlobalInterpolator.h"

#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Io/IoBase.h"

#include "ijedi/FieldMetadata/FieldsMetadata.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/State/State.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

class IOLandVectorParameters : public IoParametersBase {
  OOPS_CONCRETE_PARAMETERS(IOLandVectorParameters, IoParametersBase)

 public:

  oops::Parameter<std::string> mode{"mode", "read/write", "write", this};

  // Filenames provided as a list for input
  oops::OptionalParameter<std::vector<std::string>> filenames{"filenames",
                                                              "names of the files to be read",
                                                              this};
  // Filename of output
  oops::Parameter<std::string> filename{"filename", "filename", 
                                        "landVector_%Y%m%dT%H%M%S.nc4", this};

  // Filename of geom input (for reading lat/lon/elevation from a NetCDF file)
  oops::OptionalParameter<std::string> geomfilename{"geom filename",
                                                     "geometry NetCDF filename for read()",
                                                     this};

    // Path prepended to all files
  oops::Parameter<std::string> datapath{"datapath", "path to location of files to be read",
                                        "./", this};

  oops::OptionalParameter<std::string> dateTime{"date time", 
                                                 "date time for read/write", this};

  //***TODO: handle this
  // Interpolator type
  oops::Parameter<std::string> interpolator{"local interpolator type", "local interpolator type",
                                            "oops unstructured grid interpolator",
                                            this};

  // Optionally config may contain member
  oops::OptionalParameter<int> member{"member", "ensemble member number", this};

  // Floating point precision in bytes for NetCDF write
  oops::Parameter<int> floatPrecision{"float precision in bytes", "float precision in bytes", 8,
                                      this};

  // Geometry var names
  oops::Parameter<std::string> locationName{"location dime name", "location dim name", "location", this};
  oops::Parameter<std::string> layerName{"layer dim name", "layer dim name", "soil_levels", this};
  oops::Parameter<std::string> timeName{"time dim name", "time dim name", "time", this};

};

// -------------------------------------------------------------------------------------------------
class IOLandVector : public IoBase, private util::ObjectCounter<IOLandVector> {
 public:
  static const std::string classname() {return "ijedi::IOLandVector";}

  typedef IOLandVectorParameters Parameters_;

  IOLandVector(const Geometry &, const Parameters_ &);
  ~IOLandVector();

  /*void read(State &, const eckit::LocalConfiguration &,
            const eckit::LocalConfiguration &) const override;
  void read(Increment &, const eckit::LocalConfiguration &,
            const eckit::LocalConfiguration &) const override;
  void write(const State &, const eckit::LocalConfiguration &,
             const eckit::LocalConfiguration &) const override;
  void write(const Increment &, const eckit::LocalConfiguration &,
             const eckit::LocalConfiguration &) const override;*/
  
  void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
                  const eckit::LocalConfiguration &) const override;
  void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
                   const eckit::LocalConfiguration &) const override;

 private:
  // Methods
  void print(std::ostream &) const override;
  template <typename T>
  void writeVector(const T & obj, const std::string & label,
                      const eckit::LocalConfiguration & fileionames,
                      const eckit::LocalConfiguration & fileioscaling) const;
  void writeVectorFields(const atlas::FieldSet &, //const util::DateTime &,
                             size_t num_points,
                             const eckit::LocalConfiguration &,
                             const eckit::LocalConfiguration &) const;
  void readVectorFields(std::string pathFile,
                            atlas::FieldSet &, const util::DateTime &,
                            size_t, size_t,
                            const eckit::LocalConfiguration &,
                            const eckit::LocalConfiguration &) const;


  const Geometry & geom_;
  std::string vectorStr_;
  Parameters_ params_;
 
  atlas::Field readLonLat_;

   //std::unique_ptr<oops::GlobalInterpolator> interpolator_;
  //mutable std::unique_ptr<oops::GlobalInterpolator> interpolator_;
   //std::unique_ptr<atlas::functionspace::StructuredColumns> functionSpace_;
  //mutable std::unique_ptr<atlas::functionspace::StructuredColumns> functionSpace_;  

};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi
