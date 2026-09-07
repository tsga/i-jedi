/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

 // A geometry class that only needs lat/lon/elevation information--for land (snow/soil) DA over independent grids.
 // Written with the intent of being a minimal geometry class that can be used for land DA over independent grids, without requiring a full model geometry.
 // Written with the help of Gemeneni, fully reviewed and tested by land DA team.


/*Sample YAML file to configure the landVector geometry:
geometry:
  geometry type: landVector
  lats: [34.05, 36.16, 39.73]
  lons: [-118.24, -115.13, -104.99]
  elevations: [89.0, 610.0, 1609.0
*/

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/grid/Grid.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ijedi/Geometry/base/GeometryBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi {

class GeometryLandVector : public GeometryBase {
 public:
  GeometryLandVector(const eckit::Configuration &, const eckit::mpi::Comm &,
                  eckit::LocalConfiguration &, atlas::FunctionSpace &,
                  atlas::FieldSet &, bool &, int &);  
  GeometryLandVector(const GeometryLandVector &);
  ~GeometryLandVector();

  static const std::string classname() { return "ijedi::GeometryLandVector"; }

  std::vector<double> verticalCoord(std::string &) const override;

  void print(std::ostream &) const override;

  const eckit::mpi::Comm & getComm() const { return comm_; }
  const atlas::FunctionSpace & functionSpace() const { return functionSpace_; }
  const atlas::FieldSet & fields() const { return fields_; }

  // Spatial metadata accessors
  size_t globalNodeCount() const;
  size_t localNodeCount() const;

  
 private:

  void readLatLonElevFromFile(const std::string, const std::string, const std::string, const std::string,
    std::vector<double> &, std::vector<double> &, std::vector<double> &);
    
  const eckit::mpi::Comm & comm_;
  atlas::FunctionSpace functionSpace_;
  atlas::FieldSet fields_;
 
  atlas::Grid grid_;

  int numLevels_ = 1;

  // Atlas field containers for coordinates
  //atlas::Field lonlatField_;
  //atlas::Field globalIndexField_;
  atlas::Field elevationField_;
  atlas::Field lat_;
  atlas::Field lon_;

  atlas::Field lonlat_;
  atlas::Field partition_; 
  mutable atlas::Field ghost_;
  atlas::Field global_index_;
  mutable atlas::Field remote_index_;
  atlas::Field owned_;
  atlas::Field area_;
  atlas::Field vertical_;

};

}  // namespace ijedi
