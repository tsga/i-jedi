/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

 // A geometry class that only needs lat/lon/elevation information--for land (snow/soil) DA over independent grids.
 // Written with the intent of being a minimal geometry class that can be used for land DA over independent grids, without requiring a full model geometry.
 // Written with the help of Gemeni, fully reviewed and tested by land DA team.
  
#include <netcdf.h>

#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/UnstructuredGrid.h"

#include "atlas/util/Config.h"
#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/util/Metadata.h"
#include "oops/util/FieldSetHelpers.h"

#include "ijedi/Geometry/landVector/GeometryLandVector.h"

namespace ijedi {

static inline void nc_rc(const int return_code, const std::string & operation) {
  if (return_code != NC_NOERR) {
    std::string errMsg = "GeometryLandVector netCDF operation '" + operation +
                         "' failed with error: " + nc_strerror(return_code);
    oops::Log::error() << errMsg << std::endl;
    throw eckit::Exception(errMsg, Here());
  }
}

// -----------------------------------------------------------------------------
GeometryLandVector::GeometryLandVector(const eckit::Configuration &config,
                                       const eckit::mpi::Comm &comm,
                                       eckit::LocalConfiguration & /*geomVariables*/,
                                       atlas::FunctionSpace &functionSpace,
                                       atlas::FieldSet &geomFields,
                                       bool &levelsAreTopDown, int &numLevels)
    : comm_(comm) {
  oops::Log::trace() << "GeometryLandVector::GeometryLandVector starting..." << std::endl;
 
  numLevels_ = config.getInt("nlevels", 1);
  numLevels = numLevels_;
  levelsAreTopDown = false;

  // 1. Read coordinates (lat, lon, elevation) from configuration or input file
  std::string latlonFile = config.getString("latlon_file", "");
  std::string latName = config.getString("lat_name", "latitude");
  std::string lonName = config.getString("lon_name", "longitude");
  std::string elevName = config.getString("elev_name", "elevation");

   std::vector<double> lats;
   std::vector<double> lons; 
   std::vector<double> elevations;

  if(!latlonFile.empty() && config.has("lat_name") && config.has("lon_name") && config.has("elev_name")) {
    oops::Log::info() << "Reading lat/lon/elevation from file: " << latlonFile << std::endl;
    readLatLonElevFromFile(latlonFile, latName, lonName, elevName, lats, lons, elevations);
  } else if (config.has("lats") && config.has("lons") && config.has("elevations")) {
    oops::Log::info() << "Reading lat/lon/elevation from configuration." << std::endl;
    std::vector<double> lats = config.getDoubleVector("lats");
    std::vector<double> lons = config.getDoubleVector("lons");
    std::vector<double> elevations = config.getDoubleVector("elevations");
  }
  else {
    throw eckit::BadParameter("GeometryLandVector: Must specify either latlon_file or lats/lons/elevations in configuration.", Here());
  }

  ASSERT(lats.size() == lons.size());
  ASSERT(lats.size() == elevations.size());

  int myRank = comm_.rank();

  const size_t numPoints = lats.size();
  std::vector<atlas::PointXY> global_pts(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    global_pts[i] = atlas::PointXY(lons[i], lats[i]);
  }

  // 2. Wrap global points into an UnstructuredGrid and partition them across tasks
  //atlas::UnstructuredGrid grid(global_pts);
  grid_ = atlas::UnstructuredGrid(global_pts);
  atlas::grid::Partitioner partitioner("equal_regions", comm_.size());
  atlas::grid::Distribution distribution = partitioner.partition(grid_);

  // 3. Construct PointCloud functionspace WITH grid and distribution!
  // THIS automatically configures the internal gather/scatter engine in Atlas
  eckit::LocalConfiguration fs_config;
  fs_config.set("mpi_comm", comm_.name());
  functionSpace = atlas::functionspace::PointCloud(grid_, partitioner, fs_config);

  // 4. Create metadata fields on functionSpace
  const size_t localSize = functionSpace.size();

  lonlat_ = functionSpace.createField<double>(
      atlas::option::name("lonlat") | atlas::option::variables(2));
  
  elevationField_ = functionSpace.createField<double>(
      atlas::option::name("height"));

  auto lonlatView = atlas::array::make_view<double, 2>(lonlat_);
  auto elevView   = atlas::array::make_view<double, 1>(elevationField_);

  // Extract owned local points from functionSpace
  auto pc = atlas::functionspace::PointCloud(functionSpace);
  size_t localIdx = 0;
  for (size_t i = 0; i < numPoints; ++i) {
    if (distribution.partition(i) == comm_.rank()) {
      lonlatView(localIdx, 0) = global_pts[i].x();
      lonlatView(localIdx, 1) = global_pts[i].y();
      elevView(localIdx)      = elevations[i];
      localIdx++;
    }
  }
  
  // Explicit 1D latitude and longitude fields for OOPS KD-tree search
  auto lat_ = functionSpace.createField<double>(atlas::option::name("latitude"));
  auto lon_ = functionSpace.createField<double>(atlas::option::name("longitude"));
  
  auto latView = atlas::array::make_view<double, 1>(lat_);
  auto lonView = atlas::array::make_view<double, 1>(lon_);
  for (size_t i = 0; i < localSize; ++i) {
    lonView(i) = lonlatView(i, 0);
    latView(i) = lonlatView(i, 1);
  }

    // Create and populate global_index field on the FunctionSpace
  global_index_ = functionSpace.createField<atlas::gidx_t>(
      atlas::option::name("global_index"));

  auto gidxView = atlas::array::make_view<atlas::gidx_t, 1>(global_index_);

  localIdx = 0;
  for (size_t i = 0; i < numPoints; ++i) {
    if (distribution.partition(i) == myRank) {
      // Note: Atlas uses 1-based global indexing by convention
      gidxView(localIdx++) = static_cast<atlas::gidx_t>(i + 1); 
    }
  }

  // Create remote_index field (0-based local index on the remote owning rank)
  remote_index_ = functionSpace.createField<int>(atlas::option::name("remote_index"));
  auto ridxView  = atlas::array::make_view<int, 1>(remote_index_);

  // Create ghost field (0 = owned node, 1 = ghost node)
  ghost_ = functionSpace.createField<int>(atlas::option::name("ghost"));
  auto ghostView  = atlas::array::make_view<int, 1>(ghost_);

  // Create partition field (myRank for all local points)
  partition_ = functionSpace.createField<int>(atlas::option::name("partition"));
  auto partView  = atlas::array::make_view<int, 1>(partition_);

  localIdx = 0;
  for (size_t i = 0; i < numPoints; ++i) {
    if (distribution.partition(i) == myRank) {
      ridxView(localIdx)  = static_cast<int>(localIdx);        // Index on owning rank
      partView(localIdx) = comm_.rank();
      ghostView(localIdx) = 0; // All points in local_pts are owned locally
      localIdx++;
    }
  }

  // 1. Create 'owned' field (1 for owned point, 0 for ghost point)
  owned_ = functionSpace.createField<int>(atlas::option::name("owned"));
  auto ownedView  = atlas::array::make_view<int, 1>(owned_);
  area_ = functionSpace.createField<double>(atlas::option::name("area"));
  auto areaView  = atlas::array::make_view<double, 1>(area_);

  for (size_t i = 0; i < localSize; ++i) {
    ownedView(i) = 1;          // EVERY local point in this rank is owned
    areaView(i) = 1.0;
  }

  // Add ALL of them to geomFields
  geomFields.add(lonlat_);
  geomFields.add(elevationField_);
  geomFields.add(lat_);
  geomFields.add(lon_); 
  geomFields.add(global_index_); 
  geomFields.add(remote_index_);
  geomFields.add(ghost_);
  geomFields.add(partition_);
  geomFields.add(area_);
  geomFields.add(owned_);

  functionSpace_ = functionSpace;
  fields_ = geomFields;

  oops::Log::info() << "Registered fields in GeometryLandVector:" << std::endl;

  oops::Log::info() << "=== GeometryData Compatibility Check ===" << std::endl;
  oops::Log::info() << "Fields present in geom.fields():" << std::endl;
  for (const auto & f : fields_) {
    oops::Log::info() << "  - " << f.name() 
                      << " | Rank: " << f.rank() 
                      << " | Shape: [" << f.shape(0) 
                      << (f.rank() > 1 ? ", " + std::to_string(f.shape(1)) : "") << "]" 
                      << std::endl;
  }

  oops::Log::info() << "=== GeometryLandVector Diagnostic ===" << std::endl;
  oops::Log::info() << "  Local node count:  " << localNodeCount() << std::endl;
  oops::Log::info() << "  Global node count: " << globalNodeCount() << std::endl;

  if (fields_.has("lonlat")) {
    auto llView = atlas::array::make_view<double, 2>(fields_["lonlat"]);
    if (localNodeCount() > 0) {
      oops::Log::info() << "  Sample Point 0: Lon = " << llView(0, 0) 
                        << ", Lat = " << llView(0, 1) << std::endl;
      // Check for NaN or 0,0 coordinates
      ASSERT(!std::isnan(llView(0,0)) && !std::isnan(llView(0,1)));   //,"Geometry coordinates contain NaN!");
    }
  } else {
    oops::Log::error() << "CRITICAL: 'lonlat' field missing in fields_!" << std::endl;
  }

  // Test if Atlas can gather points globally across MPI tasks
  atlas::Field globalLonLat = functionSpace_.createField<double>(
      atlas::option::name("global_lonlat") | atlas::option::variables(2) | atlas::option::global());

  try {
    functionSpace_.gather(lonlat_, globalLonLat);
    oops::Log::info() << "Atlas gather successful. Global shape: " 
                      << globalLonLat.shape(0) << " x " << globalLonLat.shape(1) << std::endl;
  } catch (const eckit::Exception & e) {
    oops::Log::error() << "Atlas gather failed: " << e.what() << std::endl;
  }

  // Diagnostic: Force test GeometryData's exact gather operation
  atlas::Field local_latlon = geomFields["lonlat"];
  atlas::Field global_latlon = functionSpace.createField<double>(
      atlas::option::name("global_lonlat") | 
      atlas::option::variables(2) | 
      atlas::option::global());

  functionSpace.gather(local_latlon, global_latlon);

  oops::Log::info() << ">>> GEOMETRY GATHER TEST RESULT <<<" << std::endl;
  oops::Log::info() << "Global gathered shape: " << global_latlon.shape(0) 
                    << " x " << global_latlon.shape(1) << std::endl;

  if (comm_.rank() == 0) {
    ASSERT(global_latlon.shape(0) > 0); //, "Atlas functionSpace.gather() returned 0 points! GeometryData cannot build globalNodeTree_.");
  }

  oops::Log::trace() << "GeometryLandVector::GeometryLandVector completed with "
                     << numPoints << " global and " << localSize << " local points." << std::endl;

}

// -----------------------------------------------------------------------------
/* GeometryLandVector::GeometryLandVector(const eckit::Configuration &config,
                                       const eckit::mpi::Comm &comm,
                                       eckit::LocalConfiguration &,    //geomVariables,
                                       atlas::FunctionSpace &functionSpace,
                                       atlas::FieldSet &geomFields,
                                       bool &levelsAreTopDown, int &numLevels)
    : comm_(comm) {
  oops::Log::trace() << "GeometryLandVector::GeometryLandVector starting..." << std::endl;
 
  numLevels_ = config.getInt("nlevels", 1);
  numLevels = numLevels_;

  // 1. Read coordinates (lat, lon, elevation) from configuration or input file
  std::string latlonFile = config.getString("latlon_file", "");
  std::string latName = config.getString("lat_name", "latitude");
  std::string lonName = config.getString("lon_name", "longitude");
  std::string elevName = config.getString("elev_name", "elevation");

   std::vector<double> lats;
   std::vector<double> lons; 
   std::vector<double> elevations;

  if(!latlonFile.empty() && config.has("lat_name") && config.has("lon_name") && config.has("elev_name")) {
    oops::Log::info() << "Reading lat/lon/elevation from file: " << latlonFile << std::endl;
    readLatLonElevFromFile(latlonFile, latName, lonName, elevName, lats, lons, elevations);
  } else if (config.has("lats") && config.has("lons") && config.has("elevations")) {
    oops::Log::info() << "Reading lat/lon/elevation from configuration." << std::endl;
    std::vector<double> lats = config.getDoubleVector("lats");
    std::vector<double> lons = config.getDoubleVector("lons");
    std::vector<double> elevations = config.getDoubleVector("elevations");
  }
  else {
    throw eckit::BadParameter("GeometryLandVector: Must specify either latlon_file or lats/lons/elevations in configuration.", Here());
  }

  ASSERT(lats.size() == lons.size());
  ASSERT(lats.size() == elevations.size());

  size_t numPoints = lats.size();

  // Global point cloud geometry
  std::vector<atlas::PointXY> global_pts(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    global_pts[i] = atlas::PointXY(lons[i], lats[i]);
  }

  atlas::UnstructuredGrid grid(global_pts);

  // Create partitioner (EqualRegions splits points evenly across MPI tasks)
  atlas::grid::Partitioner partitioner("equal_regions", comm_.size());
  atlas::grid::Distribution distribution = partitioner.partition(grid);

  // Extract points owned by rank
  int myRank = comm_.rank();
  std::vector<atlas::PointXY> local_pts;
  std::vector<double> local_elevations;

  for (size_t i = 0; i < numPoints; ++i) {
    if (distribution.partition(i) == myRank) {
      local_pts.push_back(global_pts[i]);
      local_elevations.push_back(elevations[i]);
    }
  }

  // Initialize PointCloud functionspace across local MPI tasks
  functionSpace_ = atlas::functionspace::PointCloud(local_pts);

  // Create and populate global_index field on the FunctionSpace
  global_index_ = functionSpace_.createField<atlas::gidx_t>(
      atlas::option::name("global_index"));

  auto gidxView = atlas::array::make_view<atlas::gidx_t, 1>(global_index_);

  size_t localIdx = 0;
  for (size_t i = 0; i < numPoints; ++i) {
    if (distribution.partition(i) == myRank) {
      // Note: Atlas uses 1-based global indexing by convention
      gidxView(localIdx++) = static_cast<atlas::gidx_t>(i + 1); 
    }
  }

  // Register standard JEDI coordinate fields (lonlat and vertical height/elevation)
  lonlat_ = functionSpace_.createField<double>(
      atlas::option::name("lonlat") | atlas::option::variables(2));
  
  elevationField_ = functionSpace_.createField<double>(
      atlas::option::name("height") | atlas::option::variables(1));

  // Populate field values
  auto lonlatView = atlas::array::make_view<double, 2>(lonlat_);
  auto elevView   = atlas::array::make_view<double, 2>(elevationField_);

  for (size_t i = 0; i < local_pts.size(); ++i) {
    lonlatView(i, 0) = local_pts[i].x();   //lons;
    lonlatView(i, 1) = local_pts[i].y();   //lats;
    elevView(i, 0)      = local_elevations[i];
  }

  //Create explicit scalar fields for latitude and longitude
  lat_ = functionSpace_.createField<double>(atlas::option::name("latitude"));
  lon_ = functionSpace_.createField<double>(atlas::option::name("longitude"));

  auto latView = atlas::array::make_view<double, 1>(lat_);
  auto lonView = atlas::array::make_view<double, 1>(lon_);

  for (size_t i = 0; i < local_pts.size(); ++i) {
    lonView(i) = local_pts[i].x();
    latView(i) = local_pts[i].y();
  }

  // Create ghost field (0 = owned node, 1 = ghost node)
  auto ghost_ = functionSpace_.createField<int>(atlas::option::name("ghost"));
  auto ghostView  = atlas::array::make_view<int, 1>(ghost_);
  for (size_t i = 0; i < local_pts.size(); ++i) {
    ghostView(i) = 0; // All points in local_pts are owned locally
  }

  // Create partition field (myRank for all local points)
  auto partition_ = functionSpace_.createField<int>(atlas::option::name("partition"));
  auto partView  = atlas::array::make_view<int, 1>(partition_);
  for (size_t i = 0; i < local_pts.size(); ++i) {
    partView(i) = comm_.rank();
  }

  // Create remote_index field (0-based local index on the remote owning rank)
  remote_index_ = functionSpace_.createField<int>(atlas::option::name("remote_index"));
  auto ridxView  = atlas::array::make_view<int, 1>(remote_index_);

  localIdx = 0;
  for (size_t i = 0; i < numPoints; ++i) {
    if (distribution.partition(i) == myRank) {
      ridxView(localIdx)  = static_cast<int>(localIdx);        // Index on owning rank
      localIdx++;
    }
  }

  // Add fields to geometry fieldset so JEDI ObsOperators / Increment objects can query them
  fields_.add(lonlat_);
  fields_.add(elevationField_);
  fields_.add(global_index_);
  fields_.add(lat_);
  fields_.add(lon_);
  // 3. Add to geometry fields_
  fields_.add(ghost_);
  fields_.add(partition_);
  fields_.add(remote_index_);

  functionSpace = functionSpace_;
  geomFields = fields_;

  // Explicitly construct Atlas GatherScatter structures for PointCloud
  auto pc = atlas::functionspace::PointCloud(functionSpace_);
  pc.setupGatherScatter();  

  oops::Log::info() << "Registered fields in GeometryLandVector:" << std::endl;
  for (const auto & f : fields_) {
      oops::Log::info() << "  - Field name: " << f.name() << std::endl;
  }

  oops::Log::info() << "=== GeometryLandVector Diagnostic ===" << std::endl;
  oops::Log::info() << "  Local node count:  " << localNodeCount() << std::endl;
  oops::Log::info() << "  Global node count: " << globalNodeCount() << std::endl;

  if (fields_.has("lonlat")) {
    auto view = atlas::array::make_view<double, 2>(fields_["lonlat"]);
    if (localNodeCount() > 0) {
      oops::Log::info() << "  Sample Point 0: Lon = " << view(0, 0) 
                        << ", Lat = " << view(0, 1) << std::endl;
    }
  } else {
    oops::Log::error() << "CRITICAL: 'lonlat' field missing in fields_!" << std::endl;
  }

  // Test if Atlas can gather points globally across MPI tasks
  atlas::Field globalLonLat = functionSpace_.createField<double>(
      atlas::option::name("global_lonlat") | atlas::option::variables(2) | atlas::option::global());

  try {
    functionSpace_.gather(lonlat_, globalLonLat);
    oops::Log::info() << "Atlas gather successful. Global shape: " 
                      << globalLonLat.shape(0) << " x " << globalLonLat.shape(1) << std::endl;
  } catch (const eckit::Exception & e) {
    oops::Log::error() << "Atlas gather failed: " << e.what() << std::endl;
  }

  oops::Log::trace() << "GeometryLandVector::GeometryLandVector completed with "
                     << numPoints << " global and " << local_pts.size() << " local points." << std::endl;

}*/

// -----------------------------------------------------------------------------
GeometryLandVector::GeometryLandVector(const GeometryLandVector & other)
    : comm_(other.comm_),
      functionSpace_(other.functionSpace_),
      fields_(other.fields_),
      lonlat_(other.lonlat_),
      elevationField_(other.elevationField_),
      global_index_(other.global_index_),
      lat_(other.lat_),
      lon_(other.lon_),
      ghost_(other.ghost_),
      partition_(other.partition_),
      remote_index_(other.remote_index_),
      owned_(other.owned_),
      area_(other.area_)  {}

// -----------------------------------------------------------------------------
GeometryLandVector::~GeometryLandVector() {}

// -----------------------------------------------------------------------------
size_t GeometryLandVector::globalNodeCount() const {
  size_t localSize = functionSpace_.size();
  size_t globalSize = 0;
  // Sum local sizes across all MPI ranks in comm_
  comm_.allReduce(localSize, globalSize, eckit::mpi::Operation::SUM);

  return globalSize;
}

// -----------------------------------------------------------------------------
size_t GeometryLandVector::localNodeCount() const {
  return functionSpace_.size();
}

// -----------------------------------------------------------------------------
void GeometryLandVector::print(std::ostream & os) const {
  os << "GeometryLandVector: [ Count = " << globalNodeCount() << " points ]";
}

// -----------------------------------------------------------------------------------------------

std::vector<double> GeometryLandVector::verticalCoord(std::string & /*vcUnits*/) const {
  std::stringstream errorMsg;
  errorMsg << "GeometryLandVector::verticalCoord is not implemented" << std::endl;
  throw eckit::NotImplemented(errorMsg.str(), Here());
}

void GeometryLandVector::readLatLonElevFromFile(const std::string pathFile, const std::string latName, const std::string lonName, const std::string elevName,
    std::vector<double> & lats, std::vector<double> & lons, std::vector<double> & elevations)
{
    // NetCDF IDs
    int fileId;

    oops::Log::info() << "Reading file " << pathFile << std::endl;
    nc_rc(nc_open(pathFile.c_str(), NC_NOWRITE, &fileId), "nc_open " + pathFile);

    // Get the variable ID for this field
    int varId;
    nc_rc(nc_inq_varid(fileId, latName.c_str(), &varId), "nc_inq_varid " + latName);

    // Get number of dimensions + their IDs
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    nc_rc(nc_inq_var(fileId, varId,
                     nullptr,   // var name (unused)
                     nullptr,   // type (unused)
                     &ndims,
                     dimids,
                     nullptr),  // attributes (unused)
          "nc_inq_var");

    // Ensure that we are working with a vector (1D variables)
    ASSERT(ndims == 1);
    // Read the variable data
    size_t dimlen;
    nc_rc(nc_inq_dimlen(fileId, dimids[0], &dimlen), "nc_get_dimlen " + latName);
    lats.resize(dimlen);
    lons.resize(dimlen);
    elevations.resize(dimlen);
    
    nc_rc(nc_get_var_double(fileId, varId, lats.data()), "nc_get_var_double " + latName);

    nc_rc(nc_inq_varid(fileId, lonName.c_str(), &varId), "nc_inq_varid " + lonName);
    nc_rc(nc_get_var_double(fileId, varId, lons.data()), "nc_get_var_double " + lonName);

    nc_rc(nc_inq_varid(fileId, elevName.c_str(), &varId), "nc_inq_varid " + elevName);
    nc_rc(nc_get_var_double(fileId, varId, elevations.data()), "nc_get_var_double " + elevName);
    
    // Close file
    nc_rc(nc_close(fileId), "nc_close "+ pathFile);

    oops::Log::info() << "Done reading file " << pathFile << std::endl;

}

// -----------------------------------------------------------------------------------------------

}  // namespace ijedi