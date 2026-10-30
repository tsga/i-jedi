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
#include "atlas/parallel/GatherScatter.h"

#include "atlas/util/KDTree.h"

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

  // Enable default C++ shallow copy semantics for OOPS cloning/ensemble creation
  GeometryLandVector(const GeometryLandVector & other) = default;
  GeometryLandVector & operator=(const GeometryLandVector & other) = default;

  ~GeometryLandVector() = default;

  // Accessors delegating directly to shared data
  const atlas::FunctionSpace & functionSpace() const {
    return sharedData_->functionSpace_;
  }
  //const atlas::FunctionSpace & functionSpace() const { return functionSpace_; }

  const atlas::Grid & grid() const {
    return sharedData_->grid_;
  }

  const eckit::mpi::Comm & getComm() const {
    return sharedData_->comm_;
  }
  //const eckit::mpi::Comm & getComm() const { return comm_; }

  const atlas::FieldSet & fields() const { return sharedData_->fields_; }

  // Static cache manager functions
  static void clearCache();
  static size_t cacheSize();

  static const std::string classname() { return "ijedi::GeometryLandVector"; }

  std::vector<double> verticalCoord(std::string &) const override;

  void print(std::ostream &) const override;

  // Spatial metadata accessors
  //size_t globalNodeCount() const;
  //size_t localNodeCount() const;

  
 private:

  // Heavy, member-invariant Atlas resources wrapped inside a shared pointer
  struct SharedData {
    SharedData(const eckit::Configuration &, const eckit::mpi::Comm &,
                  eckit::LocalConfiguration &, atlas::FunctionSpace &,
                  atlas::FieldSet &, bool &, int &);
   
    void readLatLonElevFromFile(const std::string, const std::string, const std::string, const std::string,
      std::vector<double> &, std::vector<double> &, std::vector<double> &);
  
    // Spatial metadata accessors
    size_t globalNodeCount() const;
    size_t localNodeCount() const;

    const eckit::mpi::Comm & comm_;
    atlas::FunctionSpace functionSpace_;
    atlas::FieldSet fields_;
   
    atlas::Grid grid_;
    atlas::Mesh mesh_;
    
    //atlas::util::IndexKDTree localCellCenterTree_; 
    std::unique_ptr<atlas::util::IndexKDTree> localCellCenterTree_;

    int numLevels_ = 1;
    eckit::LocalConfiguration geomVariables_;
    bool levelsAreTopDown_;

    // Atlas field containers for coordinates
    atlas::Field elevationField_;
    atlas::Field lat_;
    atlas::Field lon_;

    atlas::Field lonlat_;
    atlas::Field partition_; 
    //mutable 
    atlas::Field ghost_;
    atlas::Field global_index_;
    //mutable 
    atlas::Field remote_index_;
    atlas::Field owned_;
    atlas::Field area_;
    atlas::Field vertical_;

    int ownedCount_ = 0;   // nG/npe (+1 for some)
    int ownedMin_ = 0;
    int ownedMax_ = 0;
  };

  // Helper to construct a unique lookup string from geometry settings
  static std::string makeCacheKey(const eckit::Configuration & config,
                                 const eckit::mpi::Comm & comm);

  // Static registry & mutex across all GeometryLandVector instances
  static std::map<std::string, std::shared_ptr<SharedData>> cache_;
  static std::mutex mutex_;

  // Shared pointer handle across shallow-copied & cached instances
  std::shared_ptr<SharedData> sharedData_;

  //atlas::parallel::GatherScatter gather_scatter_;

};

}  // namespace ijedi
