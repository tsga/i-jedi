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

class IOStructuredGridParameters : public IoParametersBase {
  OOPS_CONCRETE_PARAMETERS(IOStructuredGridParameters, IoParametersBase)

 public:
  // Type of structured grid to write
  oops::Parameter<std::string> outputGridType{"gridtype", "gridtype", "F12", this};
  oops::Parameter<std::string> mode{"mode", "read/write", "write", this};

  // Filenames provided as a list for input
  oops::OptionalParameter<std::vector<std::string>> filenames{"filenames",
                                                              "names of the files to be read",
                                                              this};
  // Filename of output
  oops::Parameter<std::string> filename{"filename", "filename",  
                                        "cube_to_geometric_%Y%m%dT%H%M%S.nc4", this};

  // Filename of geom input (for reading external structured-grid files)
  oops::OptionalParameter<std::string> geomfilename{"geom filename",
                                                     "geometry NetCDF filename for read()",
                                                     this};

  // Filename of input (for reading external structured-grid files)
  oops::OptionalParameter<std::string> inputFilename{"input filename",
                                                     "input NetCDF filename for read()",
                                                     this};

    // Path prepended to all files
  oops::Parameter<std::string> datapath{"datapath", "path to location of files to be read",
                                        "./", this};

  // Flag to indicate whether to remap vertical coordinates based on orography
  oops::Parameter<bool> doVerticalRemapping{"do vertical remapping",
                                            "do vertical remapping", false, this};

  // Orography filename
  oops::OptionalParameter<std::string> orographyFilename{"orography filename",
                                                         "orography filename", this};

  oops::OptionalParameter<std::string> dateTime{"date time", 
                                                 "date time for read/write", this};

  // Interpolator type
  oops::Parameter<std::string> interpolator{"local interpolator type", "local interpolator type",
                                            "oops unstructured grid interpolator",
                                            this};

  // Optionally config domain region boundaries
  oops::OptionalParameter<float> lon_min{"lon_min", "minimum longitude for read in",this};
  oops::OptionalParameter<float> lon_max{"lon_max", "maximum longitude for read in",this};
  oops::OptionalParameter<float> lat_min{"lat_min", "minimum latitude for read in",this};
  oops::OptionalParameter<float> lat_max{"lat_max", "maximum latitude for read in",this};

  // Optionally config may contain member
  oops::OptionalParameter<int> member{"member", "ensemble member number", this};

  // Floating point precision in bytes for NetCDF write
  oops::Parameter<int> floatPrecision{"float precision in bytes", "float precision in bytes", 8,
                                      this};

  // Dimension names
  oops::Parameter<std::string> latName{"latitude dim name", "latitude dim name", "grid_yt", this};
  oops::Parameter<std::string> lonName{"longitude dim name", "longitude dim name", "grid_xt", this};
  oops::Parameter<std::string> latvarName{"latitude var name", "latitude var name", "lat", this};
  oops::Parameter<std::string> lonvarName{"longitude var name", "longitude var name", "lon", this};
  oops::Parameter<std::string> levName{"level dim name", "level dim name", "pfull", this};
  oops::Parameter<std::string> edgName{"edge dim name", "edge dim name", "phalf", this};  //TODO: this may be wrong
  oops::Parameter<std::string> forName{"four level dim name", "four level dim name", "four", this};
  oops::Parameter<std::string> timName{"time dim name", "time dim name", "time", this};
};

// -------------------------------------------------------------------------------------------------
class IOStructuredGrid : public IoBase, private util::ObjectCounter<IOStructuredGrid> {
 public:
  static const std::string classname() {return "ijedi::IOStructuredGrid";}

  typedef IOStructuredGridParameters Parameters_;

  IOStructuredGrid(const Geometry &, const Parameters_ &);
  ~IOStructuredGrid();

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
  void writeGauss(const T & obj, const std::string & label,
                      const eckit::LocalConfiguration & fileionames,
                      const eckit::LocalConfiguration & fileioscaling) const;
  void writeStructuredFields(const atlas::FieldSet &, const util::DateTime &,
                             const eckit::LocalConfiguration &,
                             const eckit::LocalConfiguration &) const;
  void readStructuredFields(std::string pathFile,
                            atlas::FieldSet &, const util::DateTime &,
                            const eckit::LocalConfiguration &,
                            const eckit::LocalConfiguration &) const;

  // Data
  //std::unique_ptr<oops::GlobalInterpolator> interpolator_;
  //mutable std::unique_ptr<oops::GlobalInterpolator> readInterpolator_;
  //mutable std::unique_ptr<oops::GlobalInterpolator> interpolatorBack_;  // mutable: created in read()
  const Geometry & geom_;
  std::string gridStr_;
  Parameters_ params_;
  //std::unique_ptr<atlas::functionspace::StructuredColumns> writeFunctionSpace_;
  //mutable std::unique_ptr<atlas::functionspace::StructuredColumns> readFunctionSpace_;  // mutable: used in read()
  atlas::Field readLonLat_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
