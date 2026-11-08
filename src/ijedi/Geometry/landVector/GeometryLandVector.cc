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
#include "atlas/mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/mesh/actions/BuildHalo.h"

#include "atlas/meshgenerator.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
//#include "atlas/mesh/actions/BuildDelaunay.h"
//#include "atlas/mesh/actions/Dispatch.h"

#include "atlas/output/Gmsh.h"
#include "atlas/util/Earth.h"
#include "atlas/util/KDTree.h"

#include "atlas/util/Config.h"
#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/util/Metadata.h"
#include "oops/util/FieldSetHelpers.h"
#include "atlas/util/Earth.h"
#include "atlas/util/KDTree.h"

#include "ijedi/Geometry/landVector/GeometryLandVector.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include "atlas/util/Point.h"

#include <unordered_map>
#include <vector>
#include <cmath>
#include <array>
#include "atlas/array.h"
#include "atlas/util/Constants.h"

void compute_quad_nodes(size_t ownedCount, std::vector<atlas::PointXY> global_pts, std::vector<atlas::PointXY> rank_points,
		std::vector<double> rank_lons, std::vector<double> rank_lats,
		std::vector<atlas::gidx_t> globalIdx,
                std::vector<std::array<atlas::gidx_t, 4>> &quadNodes,
		std::vector<atlas::gidx_t> &quadGidx,
		 std::vector<std::array<atlas::gidx_t, 3>> &triNodes,
                std::vector<atlas::gidx_t> &triGidx,
		const eckit::mpi::Comm &comm)
{

  // 1. Calculate bounding box and resolution step sizes (dx, dy)
  double minLon = 360.0, maxLon = -360.0;
  double minLat = 90.0,  maxLat = -90.0;
  
  for (size_t i = 0; i < rank_lons.size(); ++i) {
      minLon = std::min(minLon, rank_lons[i]);
      maxLon = std::max(maxLon, rank_lons[i]);
      minLat = std::min(minLat, rank_lats[i]);
      maxLat = std::max(maxLat, rank_lats[i]);
  }
  
  // Global reduce to get consistent min/max across all MPI ranks
  comm.allReduceInPlace(minLon, eckit::mpi::Operation::MIN);
  comm.allReduceInPlace(maxLon, eckit::mpi::Operation::MAX);
  comm.allReduceInPlace(minLat, eckit::mpi::Operation::MIN);
  comm.allReduceInPlace(maxLat, eckit::mpi::Operation::MAX);
  
  // 2. Define effective resolution based on minimum point distance (dLat / dLon)
  const double dLon = 0.25; // Target grid spacing in degrees
  const double dLat = 0.25;
  
  const int niEff = static_cast<int>(std::round((maxLon - minLon) / dLon)) + 1;
  const int njEff = static_cast<int>(std::round((maxLat - minLat) / dLat)) + 1;

  // 1. Calculate total grid size
  const size_t totalGridSize = static_cast<size_t>(niEff) * static_cast<size_t>(njEff);
  
  // 3. Fast O(1) flat lookup table mapping (iG, jG) -> 1-based globalIdx
  // Populated using ALL local nodes on the rank (owned AND ghost nodes)
  std::vector<atlas::gidx_t> pointGridTable(totalGridSize, 0);
  
  for (size_t k = 0; k < globalIdx.size(); ++k) {
    int iG = static_cast<int>(std::round((rank_lons[k] - minLon) / dLon));
    int jG = static_cast<int>(std::round((rank_lats[k] - minLat) / dLat));
  
    iG = (iG % niEff + niEff) % niEff; // Longitude wrapping
  
    if (jG >= 0 && jG < njEff) {
      size_t flatIdx = static_cast<size_t>(jG) * niEff + iG;
      pointGridTable[flatIdx] = globalIdx[k]; // 1-based node ID matching globalIdx[k]
    }
  }
  
  auto is_local_node = [&](atlas::gidx_t id) {
    return std::find(globalIdx.begin(), globalIdx.end(), id) != globalIdx.end();
  };
  
  for (size_t n = 0; n < ownedCount; ++n) {
    int iG = static_cast<int>(std::round((rank_lons[n] - minLon) / dLon));
    int jG = static_cast<int>(std::round((rank_lats[n] - minLat) / dLat));
  
    if (jG + 1 >= njEff) continue;
  
    int iG1 = (iG + 1) % niEff;
    int jG1 = jG + 1;
  
    size_t k_SW = static_cast<size_t>(jG)  * niEff + iG;
    size_t k_SE = static_cast<size_t>(jG)  * niEff + iG1;
    size_t k_NE = static_cast<size_t>(jG1) * niEff + iG1;
    size_t k_NW = static_cast<size_t>(jG1) * niEff + iG;
  
    atlas::gidx_t id_SW = pointGridTable[k_SW];
    atlas::gidx_t id_SE = pointGridTable[k_SE];
    atlas::gidx_t id_NE = pointGridTable[k_NE];
    atlas::gidx_t id_NW = pointGridTable[k_NW];

    bool has_SW = (id_SW > 0) && is_local_node(id_SW);
    bool has_SE = (id_SE > 0) && is_local_node(id_SE);
    bool has_NE = (id_NE > 0) && is_local_node(id_NE);
    bool has_NW = (id_NW > 0) && is_local_node(id_NW);
  
    int valid_count = has_SW + has_SE + has_NE + has_NW;
  
    // Case A: Full Quad
    if (valid_count == 4) {
      quadNodes.push_back({id_SW, id_SE, id_NE, id_NW});
    }
    // Case B: Fallback Triangle for Boundaries
    else if (valid_count == 3) {
      if (!has_NW) triNodes.push_back({id_SW, id_SE, id_NE});
      else if (!has_NE) triNodes.push_back({id_SW, id_SE, id_NW});
      else if (!has_SE) triNodes.push_back({id_SW, id_NE, id_NW});
      else if (!has_SW) triNodes.push_back({id_SE, id_NE, id_NW});
    }
   // // Confirm all 4 corner IDs exist locally in globalIdx before registration
   // if (id_SW > 0 && id_SE > 0 && id_NE > 0 && id_NW > 0) {
   //   if (is_local_node(id_SW) && is_local_node(id_SE) &&
   //       is_local_node(id_NE) && is_local_node(id_NW)) {
   //     quadNodes.push_back({id_SW, id_SE, id_NE, id_NW});
   //   }
   // }
  }

  // ============================================================================
  // 5. Fallback Dummy Cell Injection for Cell-Empty Ranks
  // ============================================================================
  // If an isolated coastal/island rank produces 0 quads and 0 triangles, inject 
  // a single dummy triangular cell using the rank's first 3 local nodes.
  if (quadNodes.empty() && triNodes.empty() && globalIdx.size() >= 3) {
    triNodes.push_back({globalIdx[0], globalIdx[1], globalIdx[2]});
    oops::Log::warning() << "Rank " << comm.rank() 
                       << " injected a fallback dummy triangle to satisfy oops::MeshTriangulation." 
                       << std::endl;
  } 

  // 5. Compute globally unique 1-based quad element IDs via allGather
  {
    /*const int nQ = static_cast<int>(quadNodes.size());
    std::vector<int> nPerRank(comm.size(), 0);
  
    comm.allGather(nQ, nPerRank.begin(), nPerRank.end());
  
    int offset = 0;
    for (int r = 0; r < comm.rank(); ++r) {
      offset += nPerRank[r];
    }
  
    quadGidx.resize(nQ);
    for (int q = 0; q < nQ; ++q) {
      quadGidx[q] = static_cast<atlas::gidx_t>(offset + q + 1);
    }*/
    const int nQ = static_cast<int>(quadNodes.size());
  const int nT = static_cast<int>(triNodes.size());

  std::vector<int> nQ_perRank(comm.size(), 0), nT_perRank(comm.size(), 0);
  comm.allGather(nQ, nQ_perRank.begin(), nQ_perRank.end());
  comm.allGather(nT, nT_perRank.begin(), nT_perRank.end());

  int quadOffset = 0;
  for (int r = 0; r < comm.rank(); ++r) quadOffset += nQ_perRank[r];

  int totalQuads = 0;
  for (int r = 0; r < comm.size(); ++r) totalQuads += nQ_perRank[r];

  int triOffset = totalQuads;
  for (int r = 0; r < comm.rank(); ++r) triOffset += nT_perRank[r];

  quadGidx.resize(nQ);
  for (int q = 0; q < nQ; ++q) quadGidx[q] = static_cast<atlas::gidx_t>(quadOffset + q + 1);

  triGidx.resize(nT);
  for (int t = 0; t < nT; ++t) triGidx[t]  = static_cast<atlas::gidx_t>(triOffset + t + 1);
  }
  
}

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

static inline void nc_rc(const int return_code, const std::string & operation) {
  if (return_code != NC_NOERR) {
    std::string errMsg = "GeometryLandVector netCDF operation '" + operation +
                         "' failed with error: " + nc_strerror(return_code);
    oops::Log::error() << errMsg << std::endl;
    throw eckit::Exception(errMsg, Here());
  }
}


namespace ijedi {


// Initialize static member variables
std::map<std::string, std::shared_ptr<GeometryLandVector::SharedData>>
    GeometryLandVector::cache_;
std::mutex GeometryLandVector::mutex_;

//---------------cache clear usage e.g. in tests--------------
/*#include "gtest/gtest.h"
#include "ijedi/Geometry/landVector/GeometryLandVector.h"

class GeometryLandVectorTest : public ::testing::Test {
 protected:
  void TearDown() override {
    // Purge cached geometries between tests
    ijedi::GeometryLandVector::clearCache();
    EXPECT_EQ(ijedi::GeometryLandVector::cacheSize(), 0u);
  }
};*/

// -----------------------------------------------------------------------------
// Clear Cache Function: Call this in C++ Test Fixtures (e.g. TearDown)
// -----------------------------------------------------------------------------
void GeometryLandVector::clearCache() {
  std::lock_guard<std::mutex> lock(mutex_);
  // Clearing the map drops reference counts.
  // If no GeometryLandVector objects exist in the test harness,
  // underlying Atlas grid/PointCloud memory is immediately freed.
  cache_.clear();
}

// -----------------------------------------------------------------------------
// Cache Size Inspector (Useful for assertion checks in unit tests)
// -----------------------------------------------------------------------------
size_t GeometryLandVector::cacheSize() {
  std::lock_guard<std::mutex> lock(mutex_);
  return cache_.size();
}

// -----------------------------------------------------------------------------
// Cache Key Generator
// -----------------------------------------------------------------------------
std::string GeometryLandVector::makeCacheKey(const eckit::Configuration & config,
                                             const eckit::mpi::Comm & comm) {
  // Use unique key parameters driving the PointCloud (e.g. filename, resolution)
  std::string grid_resol = config.getString("grid_resol", "C384.mx025");
  return "landVectorGeom_" + grid_resol + "_" + comm.name() + "_" + std::to_string(comm.size());
}

// -----------------------------------------------------------------------------

// Primary Constructor
// -----------------------------------------------------------------------------
GeometryLandVector::GeometryLandVector(const eckit::Configuration &config,
                                       const eckit::mpi::Comm &comm,
                                       eckit::LocalConfiguration & geomVariables,
                                       atlas::FunctionSpace &functionSpace,
                                       atlas::FieldSet &geomFields,
                                       bool &levelsAreTopDown, int &numLevels)
{
  // Read toggle flag from YAML (defaults to true if key is omitted)
  const bool useCache = config.getBool("geometry_cache", true);

  if (useCache) {
    std::string key = makeCacheKey(config, comm);

    std::lock_guard<std::mutex> lock(mutex_);

    auto it = cache_.find(key);
    if (it != cache_.end()) {
      // 1. INSTANCE CREATED FROM SCRATCH MATCHES EXISTING CACHE:
      // Reuse existing shared data handle across distinct OOPS factory calls
      sharedData_ = it->second;
    } else {
      // 2. FIRST TIME THIS GEOMETRY SPECIFICATION IS CREATED:
      // Build expensive spatial structures once and store handle in cache map
      sharedData_ = std::make_shared<SharedData>(config, comm, geomVariables, functionSpace,
                                         geomFields, levelsAreTopDown, numLevels);
      cache_[key] = sharedData_;
    }
  } else {
    // 3. CACHE BYPASSED: Construct dedicated instance directly
    sharedData_ = std::make_shared<SharedData>(config, comm, geomVariables, functionSpace,
                                               geomFields, levelsAreTopDown, numLevels);
  }
  // 4: populate output parameters
  functionSpace = sharedData_->functionSpace_;
  geomFields = sharedData_->fields_;
  //geomVariables = sharedData_->geomVariables;
  //levelsAreTopDown = sharedData_->levelsAreTopDown;
  //numLevels = sharedData_->numLevels_;

}

// -----------------------------------------------------------------------------
GeometryLandVector::SharedData::SharedData(const eckit::Configuration &config,
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

  const size_t my_rank = static_cast<size_t>(comm_.rank());
  const size_t nproc = static_cast<size_t>(comm.size());
  oops::Log::info() << "landGeom NC: Proc  " << my_rank << " running with " << nproc << " procs" << std::endl;

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
  
  const size_t num_points = lats.size();
  std::vector<atlas::PointXY> global_pts(num_points);
  for (size_t i = 0; i < num_points; ++i) {
    global_pts[i] = atlas::PointXY(lons[i], lats[i]);
  }
  
  // Sort points using Hilbert curve
//***** doesn't seem to make a difference---sortLandPointsByHilbert(global_pts);
//TODO: Sort corresponding elevation 
 

  size_t points_per_rank = num_points / nproc;
  size_t remainder = num_points % nproc;
  
  size_t start_idx = my_rank * points_per_rank + std::min<size_t>(my_rank, remainder);
  size_t local_count = points_per_rank + (my_rank < remainder ? 1 : 0);

  eckit::LocalConfiguration fs_config;
  fs_config.set("mpi_comm", comm_.name());
  fs_config.set("halo_radius", R_halo);

  //grid_ = atlas::UnstructuredGrid(global_pts);
  //atlas::grid::Partitioner partitioner("equal_regions", fs_config);
  //atlas::grid::Distribution distribution(grid_, partitioner);
  const atlas::UnstructuredGrid ugrid(global_pts);
  const atlas::grid::Distribution distribution =
      atlas::grid::Partitioner("equal_regions", nproc).partition(ugrid);

  std::vector<int> rank_owning(num_points);
  for (int i = 0; i < num_points; ++i){
    rank_owning[i] = distribution.partition(i);
  }

  std::vector<atlas::PointXY> rank_points;
  std::vector<double> rank_lons;
  std::vector<double> rank_lats;
  std::vector<int> ghosts;  //(ownedCount_,0); // 0 for owned, 1 for ghost
  std::vector<atlas::gidx_t> globalIdx; // Original 1-based global index
  std::vector<atlas::idx_t> remoteIdx;      // use 1-based
  //atlas::idx_t rem_indx = 1;
  std::vector<int> ownerLocalIdx(num_points);  // each active point's local index on its owner rank.
  std::vector<int> ownerCount(nproc, 0);
  for (int i = 0; i < num_points; ++i){
    int owner_rank = rank_owning[i];  // cust_partition[g_idx];
    ownerLocalIdx[i] = ownerCount[owner_rank]++;
  }
  //atlas::idx_t rem_indx = 1;
  for (size_t g_idx = 0; g_idx < global_pts.size(); ++g_idx) {
    if (rank_owning[g_idx] == my_rank){  //cust_partition[g_idx] == my_rank) {
        rank_points.push_back(global_pts[g_idx]);
        ghosts.push_back(0); // Owned
        globalIdx.push_back(g_idx + 1); // 1-based global ID
	rank_lons.push_back(global_pts[g_idx].x());
        rank_lats.push_back(global_pts[g_idx].y());
	remoteIdx.push_back(static_cast<atlas::idx_t>(ownerLocalIdx[g_idx] + 1));  //rem_indx++); 
    }
  }
  
  ownedCount_ = rank_points.size(); 
  int ownedCountd = static_cast<int>(distribution.nb_pts()[my_rank]);
  if (ownedCount_ != ownedCountd) {
    std::string errMsg = "GeometryLandVector: ownedCount_ " + std::to_string(ownedCount_) 
	    + " != ownedCount from dist " + std::to_string(ownedCountd);
    oops::Log::error() << errMsg << std::endl;
    throw eckit::Exception(errMsg, Here());
  }
  int ownedMin_   = static_cast<int>(distribution.min_pts());
  int ownedMax_   = static_cast<int>(distribution.max_pts());
  oops::Log::info() << "min/max owned count " << ownedMin_ << " / " << ownedMax_ << std::endl;

  // append distance-based ghost points within R_halo of owned points
// TODO: the number of added cells along a direction might be 0, 1, or more 
// do we need to limit to specific num e.g., 1?
//             __/
//       /1 \_/_\/
//           \|__/_
//           |  |
  std::vector<int> partitions(rank_points.size(), my_rank);
  std::vector<atlas::PointXY> locPoints(rank_points);
  for (size_t g_idx = 0; g_idx < global_pts.size(); ++g_idx) {
      auto pt = global_pts[g_idx];  //for (const auto& pt : global_pts) {
      if (!isOwnedByThisRank(pt, locPoints) && distanceToRankDomain(pt, locPoints) <= R_halo) {
          rank_points.push_back(pt); // Add as ghost point
	  ghosts.push_back(1); // Ghost
	  globalIdx.push_back(g_idx + 1);
	  rank_lons.push_back(pt.x());
	  rank_lats.push_back(pt.y());
          remoteIdx.push_back(static_cast<atlas::idx_t>(ownerLocalIdx[g_idx] + 1));  //local_idx_on_owner);
	  partitions.push_back(rank_owning[g_idx]);
      }
  }
  size_t nghost = std::accumulate(ghosts.begin(), ghosts.end(), 0); 
  std::vector<size_t> num_ghost_vec(nproc);
  comm_.allGather(nghost, num_ghost_vec.begin(), num_ghost_vec.end());
  auto min_ng = *std::min_element(num_ghost_vec.begin(), num_ghost_vec.end());
  auto max_ng = *std::max_element(num_ghost_vec.begin(), num_ghost_vec.end());
  oops::Log::info() << "proc 0 nghost " << nghost << " min/max numGhost " 
	  << min_ng << " / " << max_ng << std::endl;  
  
  std::vector<std::array<atlas::gidx_t, 4>> quadNodes;
  std::vector<atlas::gidx_t> quadGidx;
  std::vector<std::array<atlas::gidx_t, 3>> triNodes;
  std::vector<atlas::gidx_t> triGidx;
  compute_quad_nodes(ownedCount_, global_pts, rank_points,
                rank_lons, rank_lats,
		globalIdx,
                quadNodes, quadGidx, triNodes, triGidx,
                comm);

  eckit::LocalConfiguration meshConf;
  meshConf.set("mpi_comm", comm_.name());

  const atlas::mesh::MeshBuilder meshBuilder;
  atlas::Mesh mesh = meshBuilder(rank_lons, rank_lats, ghosts,
      globalIdx, remoteIdx, /*remote_index_base=*/1, partitions,
      /*triag_nodes_global=*{},*/ /*triag_global_index={},*/ triNodes, triGidx,
      /* quad_nodes_global=*/quadNodes,/*quad_global_index=*/quadGidx,
      meshConf);

  // 2. Compute parallel node connectivity and owner structures
  //atlas::mesh::actions::build_parallel_fields(mesh);
  // 3. Generate Delaunay 2D triangular cells on local node sets
  //atlas::mesh::actions::build_delaunay(mesh);

//  atlas::MeshGenerator meshGenerator("delaunay");
//  atlas::Mesh mesh = meshGenerator.generate(ugrid, distribution);

  // 4. Rebuild the mesh halo so adjacent cells across ranks are linked
  atlas::mesh::actions::build_halo(mesh, 1);
  
  // Diagnostic log to confirm cells exist across all ranks
  if (mesh.cells().size() < 1){
  std::cout << "Rank " << comm_.rank()
                  << " | Nodes: " << mesh.nodes().size()
                  << " | Quad Cells: " << mesh.cells().size() << std::endl;
  }

  // Guarantee cells are present so MeshTriangulation populates its KD-tree
  ASSERT(mesh.cells().size() > 0);

  functionSpace_ = atlas::functionspace::NodeColumns(mesh, atlas::util::Config("mpi_comm", comm.name()));  //fs_config);

  lonlat_ = functionSpace_.lonlat();
  ghost_ = functionSpace_.ghost();
  auto glb_idx = functionSpace_.global_index();
  auto r_idx = functionSpace_.remote_index();
  partition_ = functionSpace_.partition();

  auto lonlatView = atlas::array::make_view<double, 2>(lonlat_);
  auto glbidxView = atlas::array::make_view<atlas::gidx_t, 1>(glb_idx);
  auto partView = atlas::array::make_view<int,1>(partition_);
  auto ridxView = atlas::array::make_view<int, 1>(r_idx);

  // add global_index field: current name is 'glb_idx
  global_index_ = functionSpace_.createField<atlas::gidx_t>(atlas::option::name("global_index"));
  // Create remote_index field (1-based local index on the remote owning rank)
  remote_index_ = functionSpace_.createField<int>(atlas::option::name("remote_index"));  
  auto global_index_view = atlas::array::make_view<atlas::gidx_t, 1>(global_index_);
  auto remote_index_view  = atlas::array::make_view<int, 1>(remote_index_);
 
  // Explicit 1D latitude and longitude fields, eg for KD-tree search
  lat_ = functionSpace_.createField<double>(atlas::option::name("lat") | atlas::option::levels(1));
  lon_ = functionSpace_.createField<double>(atlas::option::name("lon") | atlas::option::levels(1));
  auto latView = atlas::array::make_view<double, 2>(lat_);
  auto lonView = atlas::array::make_view<double, 2>(lon_);

  elevationField_ = functionSpace_.createField<double>(atlas::option::name("height") | atlas::option::levels(1));
  auto elevView   = atlas::array::make_view<double, 2>(elevationField_);
  
  for (size_t i = 0; i < functionSpace_.size(); ++i){
    if (partView(i) == my_rank){
      atlas::gidx_t g_idx = glbidxView(i) - 1;; // 0-based lookup into global elevation array
      elevView(i, 0) = elevations[g_idx];  //TODO: use sorted elevation
      latView(i, 0) = lonlatView(i, 1);
      lonView(i,0) = lonlatView(i, 0);
      global_index_view(i) = glbidxView(i);
      remote_index_view(i) = ridxView(i);
    }
  }

  // Create 'owned' field (1 for owned point, 0 for ghost point)
  owned_ = functionSpace_.createField<int>(atlas::option::name("owned") | atlas::option::levels(1));
  auto owned_view = atlas::array::make_view<int, 2>(owned_);
  area_ = functionSpace_.createField<double>(atlas::option::name("area") | atlas::option::levels(1));
  auto area_view = atlas::array::make_view<double, 2>(area_); 
  for (size_t i = 0; i < functionSpace_.size(); ++i) {
      if (partView(i) == my_rank) {
          owned_view(i, 0) = 1;
          area_view(i, 0)  = 1.0;
      } else {
          owned_view(i, 0) = 0;
          area_view(i, 0)  = -1.0;
      }
  }
  
  // do halo exchange
  functionSpace_.haloExchange(elevationField_);
  functionSpace_.haloExchange(lat_);
  functionSpace_.haloExchange(lon_);
  functionSpace_.haloExchange(owned_);
  functionSpace_.haloExchange(area_);
  functionSpace_.haloExchange(global_index_);
  functionSpace_.haloExchange(remote_index_);

  std::vector<atlas::Point2> treePoints(functionSpace_.size());
  std::vector<size_t> payloads(functionSpace_.size()); 
  
  for (size_t i = 0; i < functionSpace_.size(); ++i) {
    treePoints.emplace_back(atlas::Point2(lonlatView(i, 0), lonlatView(i, 1)));
    payloads.push_back(i); // local functionspace index as payload
  }
  localCellCenterTree_ = std::make_unique<atlas::util::IndexKDTree>();  // atlas::util::IndexKDTree();
  localCellCenterTree_->build(treePoints, payloads);

  oops::Log::info() << "localCellCenterTree_.size() " << localCellCenterTree_->size() 
	<< " functionSpace_.size() " << functionSpace_.size() << std::endl;

  fields_.add(lonlat_);
  fields_.add(global_index_);
  fields_.add(remote_index_);
  fields_.add(ghost_);
  fields_.add(partition_);
  fields_.add(elevationField_);  
  fields_.add(lat_);
  fields_.add(lon_); 
  fields_.add(area_);
  fields_.add(owned_);

  functionSpace = functionSpace_;
  geomFields = fields_; 

  oops::Log::info() << "Fields present in geom.fields():" << std::endl;
  for (const auto & f : fields_) {
    oops::Log::info() << "  - " << f.name() 
                      << " | Rank: " << f.rank() 
                      << " | Shape: [" << f.shape(0) 
                      << (f.rank() > 1 ? ", " + std::to_string(f.shape(1)) : "") << "]" 
                      << std::endl;
  }

  oops::Log::info() << "=== GeometryLandVector Diagnostic ===" << std::endl;
  auto loccnt = localNodeCount();
  auto glbcnt = globalNodeCount();
  if (comm_.rank() == 0 || comm_.rank() == comm_.size() - 1) {
    std::cout << "Rank " << comm_.rank() << "  Local node count:  " << loccnt << std::endl;
    std::cout << "Rank " << comm_.rank() << "  Global node count: " << glbcnt << std::endl;
  }

  if (!fields_.has("lonlat")) {
    oops::Log::error() << "CRITICAL: 'lonlat' field missing in fields_!" << std::endl;
  }
  auto llView = atlas::array::make_view<double, 2>(fields_["lonlat"]);
  if (localNodeCount() > 0) {
      oops::Log::info() << "  Sample Point 0: Lon = " << llView(0, 0) 
                        << ", Lat = " << llView(0, 1) << std::endl;
      // Check for NaN or 0,0 coordinates
      ASSERT(!std::isnan(llView(0,0)) && !std::isnan(llView(0,1)));   //,"Geometry coordinates contain NaN!");
  }

  atlas::Field loc_latlon = functionSpace_.createField<double>(
      atlas::option::name("loc_latlon") | atlas::option::variables(2));

  atlas::Field global_latlon = functionSpace_.createField<double>(atlas::option::name("global_latlon") 
		  | atlas::option::variables(2) | atlas::option::global(0));  //

  auto locllView = atlas::array::make_view<double, 2>(loc_latlon);
  for (size_t i = 0; i < functionSpace_.size(); ++i) {
    locllView(i, 0) = llView(i, 0);
    locllView(i, 1) = llView(i, 1);
  }

  atlas::functionspace::NodeColumns fsHalo(functionSpace_);
  fsHalo.haloExchange(loc_latlon);
  functionSpace_.gather(loc_latlon, global_latlon);

  oops::Log::info() << ">>> GEOMETRY GATHER TEST RESULT <<<" << std::endl;
  oops::Log::info() << "Global gathered shape: " << global_latlon.shape(0) 
                    << " x " << global_latlon.shape(1) << std::endl;
  if (comm_.rank() == 0) {
    ASSERT(global_latlon.shape(0) > 0); //, "Atlas functionSpace.gather() returned 0 points! GeometryData cannot build globalNodeTree_.");
  }
  oops::Log::info() << "GeometryLandVector::GeometryLandVector created" << std::endl;

}

// -----------------------------------------------------------------------------
/*GeometryLandVector::GeometryLandVector(const GeometryLandVector & other)
      :sharedData_(other.sharedData_) { // Copy the shared_ptr (shallow copy
      :comm_(other.comm_),
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
      area_(other.area_)  {}*/

// -----------------------------------------------------------------------------
//GeometryLandVector::~GeometryLandVector() {}

// -----------------------------------------------------------------------------
size_t GeometryLandVector::SharedData::globalNodeCount() const {
  size_t localSize = functionSpace_.size();
  size_t globalSize = 0;
  // Sum local sizes across all MPI ranks in comm_
  comm_.allReduce(localSize, globalSize, eckit::mpi::Operation::SUM);

  return globalSize;
}

// -----------------------------------------------------------------------------
size_t GeometryLandVector::SharedData::localNodeCount() const {
  return functionSpace_.size();
}

// -----------------------------------------------------------------------------
void GeometryLandVector::print(std::ostream & os) const {
  os << "GeometryLandVector: [ Count = " << sharedData_->globalNodeCount() << " points ]";
}

// -----------------------------------------------------------------------------------------------

std::vector<double> GeometryLandVector::verticalCoord(std::string & /*vcUnits*/) const {
  std::stringstream errorMsg;
  errorMsg << "GeometryLandVector::verticalCoord is not implemented" << std::endl;
  throw eckit::NotImplemented(errorMsg.str(), Here());
}

void GeometryLandVector::SharedData::readLatLonElevFromFile(const std::string pathFile, const std::string latName, const std::string lonName, const std::string elevName,
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
