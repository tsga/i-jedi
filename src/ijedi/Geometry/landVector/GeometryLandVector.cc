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
#include "atlas/grid/Grid.h"
#include "atlas/grid/UnstructuredGrid.h"

#include "atlas/util/Config.h"
#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/util/Metadata.h"
#include "oops/util/FieldSetHelpers.h"

#include "ijedi/Geometry/landVector/GeometryLandVector.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include "atlas/util/Point.h"


/*class CustomPartitionerGenerator : public atlas::grid::PartitionerGenerator {
public:
    CustomPartitionerGenerator(const eckit::Configuration& config = atlas::util::NoConfig())
        : PartitionerGenerator(config) {}

    void populate(const atlas::Grid& grid, atlas::grid::DistributionImpl::partition_t part[]) const override {
        size_t N = grid.size();

        // Populate the `part` array where part[i] is the target MPI rank for point i
        for (size_t i = 0; i < N; ++i) {
            part[i] = get_custom_rank_for_point(i); // Your custom rank assignment logic
        }
    }

    std::string type() const override { return "custom"; }

private:
    int get_custom_rank_for_point(size_t index) const {
        // Return target rank for global point index
        return index % 4; // Example assignment
    }
};
*/

namespace {

// Standard 2D Hilbert point-to-index conversion (2D box of size N x N, where N = 2^bits)
uint64_t hilbertEncode(uint32_t x, uint32_t y, int bits) {
    uint64_t d = 0;
    for (int s = 1 << (bits - 1); s > 0; s >>= 1) {
        uint32_t rx = (x & s) > 0;
        uint32_t ry = (y & s) > 0;
        d += static_cast<uint64_t>(s) * s * ((3 * rx) ^ ry);

        // Rotate/flip quadrant
        if (ry == 0) {
            if (rx == 1) {
                x = (1U << bits) - 1 - x;
                y = (1U << bits) - 1 - y;
            }
            // Swap x and y
            std::swap(x, y);
        }
    }
    return d;
}

// Convert Lon (x) [-180, 180] or [0, 360] and Lat (y) [-90, 90] to integer grid coordinates
uint64_t pointToHilbertKey(const atlas::PointXY& pt, int bits = 16) {
    double lon = pt.x(); // degrees
    double lat = pt.y(); // degrees

    // Normalize Longitude to [0, 1)
    while (lon < 0.0) lon += 360.0;
    while (lon >= 360.0) lon -= 360.0;
    double norm_lon = lon / 360.0;

    // Normalize Latitude to [0, 1]
    double norm_lat = (lat + 90.0) / 180.0;
    norm_lat = std::max(0.0, std::min(1.0, norm_lat));

    // Scale to integer grid [0, 2^bits - 1]
    uint32_t max_val = (1U << bits) - 1;
    uint32_t x = static_cast<uint32_t>(norm_lon * max_val);
    uint32_t y = static_cast<uint32_t>(norm_lat * max_val);

    return hilbertEncode(x, y, bits);
}

} // anonymous namespace

// Call this function before passing landPoints to atlas::UnstructuredGrid
void sortLandPointsByHilbert(std::vector<atlas::PointXY>& landPoints) {
    // 1. Precompute Hilbert keys for fast comparison
    struct ScoredPoint {
        atlas::PointXY point;
        uint64_t key;
    };

    std::vector<ScoredPoint> scoredPoints;
    scoredPoints.reserve(landPoints.size());

    for (const auto& pt : landPoints) {
        scoredPoints.push_back({pt, pointToHilbertKey(pt, 16)}); // 16-bit resolution per axis
    }

    // 2. Sort by Hilbert curve key
    std::sort(scoredPoints.begin(), scoredPoints.end(),
              [](const ScoredPoint& a, const ScoredPoint& b) {
                  return a.key < b.key;
              });

    // 3. Copy back to landPoints array
    for (size_t i = 0; i < landPoints.size(); ++i) {
        landPoints[i] = scoredPoints[i].point;
    }
}

 bool is_point_in_vector(const atlas::PointXY& target, const std::vector<atlas::PointXY>& points, double epsilon = 1e-9) {
    return std::any_of(points.begin(), points.end(), [&target, epsilon](const atlas::PointXY& p) {
        return std::abs(p.x() - target.x()) < epsilon &&
               std::abs(p.y() - target.y()) < epsilon;
    });
}


namespace ijedi {

static inline void nc_rc(const int return_code, const std::string & operation) {
  if (return_code != NC_NOERR) {
    std::string errMsg = "GeometryLandVector netCDF operation '" + operation +
                         "' failed with error: " + nc_strerror(return_code);
    oops::Log::error() << errMsg << std::endl;
    throw eckit::Exception(errMsg, Here());
  }
}

// Helper to convert degrees to radians
constexpr double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

// Computes the Haversine distance in meters between two lon/lat points
double haversine_distance(const atlas::PointXY& p1, const atlas::PointXY& p2) {
    constexpr double EARTH_RADIUS_METERS = 6371000.0; // Mean Earth radius

    double lon1 = deg2rad(p1.x()); // x represents longitude
    double lat1 = deg2rad(p1.y()); // y represents latitude
    double lon2 = deg2rad(p2.x());
    double lat2 = deg2rad(p2.y());

    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    double a = std::sin(dlat / 2.0) * std::sin(dlat / 2.0) +
               std::cos(lat1) * std::cos(lat2) *
               std::sin(dlon / 2.0) * std::sin(dlon / 2.0);

    double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));

    return EARTH_RADIUS_METERS * c;
}

// Computes distances from all points in the vector to the target point
std::vector<double> compute_distances_to_target(const atlas::PointXY& target, const std::vector<atlas::PointXY>& points) {
    std::vector<double> distances;
    distances.reserve(points.size());

    for (const auto& point : points) {
        distances.push_back(haversine_distance(target, point));
    }

    return distances;
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
  double R_halo = config.getDouble("halo_radius", 125.0); //km
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

  const size_t myRank = static_cast<size_t>(comm_.rank());
  const size_t npes = static_cast<size_t>(comm.size());
  oops::Log::info() << "Proc  " << myRank << " running with " << npes << " procs" << std::endl;

  const size_t numPoints = lats.size();
  std::vector<atlas::PointXY> global_pts(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    global_pts[i] = atlas::PointXY(lons[i], lats[i]);
  }
  
  // Sort points using Hilbert curve
  sortLandPointsByHilbert(global_pts);
  //TODO: Sort corresponding elevation 
  
/*  // assign contiguous chunks to ranks
  std::vector<int> cust_partition(numPoints);
  size_t points_per_rank = numPoints / npes;
  size_t remainder = numPoints % npes;

  size_t current_idx = 0;
  for (int rank = 0; rank < npes; ++rank) {
      size_t count = points_per_rank + (rank < remainder ? 1 : 0);
      for (size_t i = 0; i < count; ++i) {
          cust_partition[current_idx++] = rank;
      }
  }

  auto isOwnedByThisRank = [&](atlas::PointXY target, std::vector<atlas::PointXY> vPoints) -> bool {
    return is_point_in_vector(target, vPoints);
  };

  auto distanceToRankDomain = [&] (atlas::PointXY target, std::vector<atlas::PointXY> vPoints) -> double {
    auto dist_vec = compute_distances_to_target(target, vPoints);
    if(!dist_vec.empty()){
      auto min_dist = *std::min_element(dist_vec.begin(), dist_vec.end());  
      return 0.001*min_dist;  // km 
    } else {
      throw eckit::BadParameter("GeometryLandVector: Error in distance computation.", Here());
    }
  };
  
  size_t start_idx = myRank * points_per_rank + std::min<size_t>(myRank, remainder);
  size_t local_count = points_per_rank + (myRank < remainder ? 1 : 0);

  std::vector<atlas::PointXY> rankPoints;
  //(global_pts.begin() + start_idx, global_pts.begin() + start_idx + local_count);
  std::vector<int> ghostMask;  //(ownedCount_,0); // 0 for owned, 1 for ghost
  std::vector<atlas::gidx_t> globalIndices; // Original 1-based global index
  for (size_t g_idx = 0; g_idx < global_pts.size(); ++g_idx) {
    if (cust_partition[g_idx] == myRank) {
        rankPoints.push_back(global_pts[g_idx]);
        ghostMask.push_back(0); // Owned
        globalIndices.push_back(g_idx + 1); // 1-based global ID
    }
  }
  
  std::vector<atlas::PointXY> locPoints(rankPoints);
  //ownedCount_ = rankPoints.size();  //local_count
   size_t numGhost = 0;
  // append distance-based ghost points within R_halo of owned points
  for (size_t g_idx = 0; g_idx < global_pts.size(); ++g_idx) {
      auto pt = global_pts[g_idx];  //for (const auto& pt : global_pts) {
      if (!isOwnedByThisRank(pt, locPoints) && distanceToRankDomain(pt, locPoints) <= R_halo) {
          rankPoints.push_back(pt); // Add as ghost point
	  ghostMask.push_back(1); // Ghost
	  numGhost++;
	  globalIndices.push_back(g_idx + 1);
      }
  }
  
  // use this spatial 2D Hilbert distribution
  //atlas::grid::Distribution distribution(npes, numPoints, cust_partition.data())
  // Construct PointCloud purely from local owned + halo points
  //functionSpace = atlas::functionspace::PointCloud(rankPoints);
*/

  eckit::LocalConfiguration fs_config;
  fs_config.set("mpi_comm", comm_.name());
  fs_config.set("halo_radius", R_halo);

  grid_ = atlas::UnstructuredGrid(global_pts);
  atlas::grid::Partitioner partitioner("equal_regions", fs_config); 
  atlas::grid::Distribution distribution(grid_, partitioner);

  functionSpace = atlas::functionspace::PointCloud(grid_, partitioner, fs_config);

  const size_t localSize = functionSpace.size();

/*  lonlat_ = functionSpace.createField<double>(
      atlas::option::name("lonlat") | atlas::option::variables(2));
  auto lonlatView = atlas::array::make_view<double, 2>(lonlat_);
  global_index_ = functionSpace.createField<atlas::gidx_t>(
      atlas::option::name("global_index") | atlas::option::levels(0));
  auto gidxView = atlas::array::make_view<atlas::gidx_t, 1>(global_index_);
  ghost_ = functionSpace.createField<int>(atlas::option::name("ghost")
                  | atlas::option::levels(0));
  auto ghostView  = atlas::array::make_view<int, 1>(ghost_);
  // Create remote_index field (0-based local index on the remote owning rank)
  remote_index_ = functionSpace.createField<atlas::idx_t>(atlas::option::name("remote_index")
                  | atlas::option::levels(0));
  auto ridxView  = atlas::array::make_view<atlas::idx_t, 1>(remote_index_);
  // Create partition field (myRank for all local points)
  partition_ = functionSpace.createField<int>(atlas::option::name("partition")
                  | atlas::option::levels(0));
  auto partView  = atlas::array::make_view<int, 1>(partition_);
*/
/*
  for (size_t l_idx = 0; l_idx < rankPoints.size(); ++l_idx) {
    lonlatView(l_idx, 0) = rankPoints[l_idx].x();  //rankPoints[localIdx].x();
    lonlatView(l_idx, 1) = rankPoints[l_idx].y();
  }
   for (size_t i = 0; i < rankPoints.size(); ++i) {
    gidxView(i) = static_cast<atlas::gidx_t>(globalIndices[i]);
    ghostView(i) = ghostMask[i]; // 0 for i < ownedCount_, 1 for ghost
    size_t g_idx = globalIndices[i] - 1;
    int owner_rank = cust_partition[g_idx];
    partView(i) = owner_rank;
    if (owner_rank == myRank) {
      // Owned point: local index on this rank is simply i
      ridxView(i) = static_cast<int>(i); //index on owning rank
    } else {
      // 0-based local index on owner_rank for ghost
      // (Assuming rank offset calculation based on points_per_rank and remainder)
      size_t start_idx = owner_rank * points_per_rank + std::min<size_t>(owner_rank, remainder);
      size_t local_idx_on_owner = g_idx - start_idx;
      ridxView(i) = static_cast<int>(local_idx_on_owner);
    }
  }

oops::Log::info() << "here 6" << std::endl;

  elevationField_ = functionSpace.createField<double>(
      atlas::option::name("height")
      | atlas::option::levels(1));
  auto elevView   = atlas::array::make_view<double, 2>(elevationField_);
  for (size_t l_idx = 0; l_idx < rankPoints.size(); ++l_idx) {
    size_t g_idx = globalIndices[l_idx] - 1; // 0-based lookup into global elevation array
    elevView(l_idx, 0) = elevations[g_idx];
  }
comm_.barrier();
oops::Log::info() << "here 7" << std::endl;

  // Explicit 1D latitude and longitude fields for OOPS KD-tree search
  lat_ = functionSpace.createField<double>(
		  atlas::option::name("lat")
		  | atlas::option::levels(1));
  lon_ = functionSpace.createField<double>(
		  atlas::option::name("lon")
		  | atlas::option::levels(1));
comm_.barrier();  
oops::Log::info() << "here 8" << std::endl;
  auto latView = atlas::array::make_view<double, 2>(lat_);
  auto lonView = atlas::array::make_view<double, 2>(lon_);
comm_.barrier();
oops::Log::info() << "here 8.5" << std::endl;
  for (size_t i = 0; i < localSize; ++i) {
    lonView(i, 0) = rankPoints[i].x();
    latView(i, 0) = rankPoints[i].y();
  }
comm_.barrier();
oops::Log::info() << "here 9" << std::endl;

  // Create 'owned' field (1 for owned point, 0 for ghost point)
  owned_ = functionSpace.createField<int>(atlas::option::name("owned")
		  | atlas::option::levels(0));
  auto owned_view  = atlas::array::make_view<int, 1>(owned_);
  area_ = functionSpace.createField<double>(atlas::option::name("area")
		  | atlas::option::levels(0));
  auto area_view  = atlas::array::make_view<double, 1>(area_);

  for (size_t i = 0; i < rankPoints.size(); ++i) {
    size_t g_idx = globalIndices[i] - 1;
    int owner_rank = cust_partition[g_idx];
    owned_view(i) = (owner_rank == myRank) ? 1 : 0;  //ghostView(i) == 0) ? 1 : 0;
    area_view(i) = (owner_rank == myRank) ? 1.0 : -1;
  }
*/

 lonlat_ = functionSpace.lonlat();
 ghost_ = functionSpace.ghost();
 global_index_ = functionSpace.global_index();
 remote_index_ = functionSpace.remote_index();
 partition_ = functionSpace.partition();
 //if (functionSpace.has_mask()) mask_ = functionSpace.mask();

  comm_.barrier();
oops::Log::info() << "here 10" << std::endl;  

  // Add fields to geomFields
  geomFields.add(lonlat_);
/*  geomFields.add(lat_);
  geomFields.add(lon_); 
  geomFields.add(area_);
  geomFields.add(owned_);*/
  geomFields.add(global_index_); 
  geomFields.add(remote_index_);
  geomFields.add(ghost_);
  geomFields.add(partition_);
//  geomFields.add(elevationField_);
comm_.barrier();
oops::Log::info() << "here 11" << std::endl;
/*  geomFields.add(functionSpace.lonlat());
  geomFields.add(functionSpace.ghost());
  geomFields.add(functionSpace.global_index());
  geomFields.add(functionSpace.remote_index());
  geomFields.add(functionSpace.partition());
  //geomFields.add(functionSpace.mask());
*/

  functionSpace_ = functionSpace;
  fields_ = geomFields; 
comm_.barrier();
oops::Log::info() << "here 11" << std::endl;

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
comm_.barrier();
oops::Log::info() << "here 15" << std::endl;

  oops::Log::info() << "=== GeometryLandVector Diagnostic ===" << std::endl;
  oops::Log::info() << "  Local node count:  " << localNodeCount() << std::endl;
  oops::Log::info() << "  Global node count: " << globalNodeCount() << std::endl;
  auto loccnt=localNodeCount();
  auto glbcnt=globalNodeCount();
  if (comm_.rank() == (comm_.size() - 1)) {
  std::cout << "  Local node count:  " << loccnt << std::endl;
  std::cout << "  Global node count: " << glbcnt << std::endl;
}
  comm_.barrier();
oops::Log::info() << "here 16" << std::endl;

//  if (fields_.has("lonlat")) {
    auto llView = atlas::array::make_view<double, 2>(fields_["lonlat"]);
    if (localNodeCount() > 0) {
      oops::Log::info() << "  Sample Point 0: Lon = " << llView(0, 0) 
                        << ", Lat = " << llView(0, 1) << std::endl;
      // Check for NaN or 0,0 coordinates
      ASSERT(!std::isnan(llView(0,0)) && !std::isnan(llView(0,1)));   //,"Geometry coordinates contain NaN!");
    }
//  } else {
//    oops::Log::error() << "CRITICAL: 'lonlat' field missing in fields_!" << std::endl;
//  }
comm_.barrier();
oops::Log::info() << "here 17" << std::endl;
  // Diagnostic: Force test GeometryData's exact gather operation
  atlas::Field f = geomFields["lonlat"];
  oops::Log::info() << "  - " << f.name()
                      << " | Rank: " << f.rank()
                      << " | Shape: [" << f.shape(0)
                      << (f.rank() > 1 ? ", " + std::to_string(f.shape(1)) : "") << "]"
                      << std::endl;
comm_.barrier();
oops::Log::info() << "here 18" << std::endl;
  auto nGlobal = static_cast<atlas::idx_t>(lats.size());
  atlas::Field global_latlon = functionSpace_.createField<double>(
      atlas::option::name("global_latlon")
      | atlas::option::variables(2) 
      | atlas::option::global());
     // | atlas::option::shape(std::initializer_list<atlas::idx_t>{nGlobal, 2}));
oops::Log::info() << "  - " << global_latlon.name()
                      << " | Rank: " << global_latlon.rank()
                      << " | Shape: [" << global_latlon.shape(0)
                      << (global_latlon.rank() > 1 ? ", " + std::to_string(global_latlon.shape(1)) : "") << "]"
                      << std::endl;
  comm_.barrier();
oops::Log::info() << "here 18.2" << std::endl;
  atlas::Field loc_latlon = functionSpace_.createField<double>(
      atlas::option::name("loc_latlon")
      | atlas::option::variables(2));
auto locllView = atlas::array::make_view<double, 2>(loc_latlon);
  for (size_t i = 0; i < localSize; ++i) {
    locllView(i, 0) = llView(i, 0);
    locllView(i, 1) = llView(i, 0);;
  }
oops::Log::info() << "  - " << loc_latlon.name()
                      << " | Rank: " << loc_latlon.rank()
                      << " | Shape: [" << loc_latlon.shape(0)
                      << (loc_latlon.rank() > 1 ? ", " + std::to_string(loc_latlon.shape(1)) : "") << "]"
                      << std::endl;
  atlas::functionspace::PointCloud fsHalo(functionSpace_);
  fsHalo.haloExchange(loc_latlon);
  fsHalo.haloExchange(global_latlon);
  comm_.barrier();
oops::Log::info() << "here 19" << std::endl;
  functionSpace_.gather(loc_latlon, global_latlon);
comm_.barrier();
oops::Log::info() << "here 23" << std::endl;
  oops::Log::info() << ">>> GEOMETRY GATHER TEST RESULT <<<" << std::endl;
  oops::Log::info() << "Global gathered shape: " << global_latlon.shape(0) 
                    << " x " << global_latlon.shape(1) << std::endl;

  if (comm_.rank() == 0) {
    ASSERT(global_latlon.shape(0) > 0); //, "Atlas functionSpace.gather() returned 0 points! GeometryData cannot build globalNodeTree_.");
  }
comm_.barrier();
oops::Log::info() << "here 25" << std::endl;
}

// -----------------------------------------------------------------------------
GeometryLandVector::GeometryLandVector(const GeometryLandVector & other)
    : comm_(other.comm_),
      functionSpace_(other.functionSpace_),
      fields_(other.fields_){}
/*      lonlat_(other.lonlat_),
      elevationField_(other.elevationField_),
      global_index_(other.global_index_),
      lat_(other.lat_),
      lon_(other.lon_),
      ghost_(other.ghost_),
      partition_(other.partition_),
      remote_index_(other.remote_index_),
      owned_(other.owned_),
      area_(other.area_)  {}
*/
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
