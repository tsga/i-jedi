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
#include <iomanip>  // Required for std::setprecision and std::fixed
#include <iostream>

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

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <array>
#include "atlas/array.h"
#include "atlas/util/Constants.h"


#include <tuple>
#include "eckit/mpi/Comm.h"


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

// -----------------------------------------------------------------------------
// POD Struct for transferring halo node metadata between MPI ranks
// -----------------------------------------------------------------------------
struct GhostNodeData {
  atlas::gidx_t globalIdx;
  int tile;
  int i;
  int j;
  double lon;
  double lat;
  atlas::idx_t remoteIdx;
  int ownerRank;
};

// -----------------------------------------------------------------------------
// Helper: Pack (tile, i, j) coordinates into a 64-bit lookup hash key
// -----------------------------------------------------------------------------
inline uint64_t pack_tile_ij(int tile, int i, int j) {
  return (static_cast<uint64_t>(tile) << 32) |
         (static_cast<uint64_t>(i & 0xFFFF) << 16) |
         static_cast<uint64_t>(j & 0xFFFF);
}

// -----------------------------------------------------------------------------
// Helper: Resolve (tile, i, j) cross-panel seam boundary transformations
// Uses 1-based indexing for tile (1..6) and coordinate indices (1..N).
// -----------------------------------------------------------------------------
inline std::tuple<int, int, int, bool> resolve_fv3_seam_neighbor(
    int tile, int i, int j, int N)
{
  // 1. Within current panel boundaries
  if (i >= 1 && i <= N && j >= 1 && j <= N) {
    return {tile, i, j, true};
  }

  // 2. Corner singularities (touching 3 panels diagonally; omitted for quad edges)
  if ((i < 1 || i > N) && (j < 1 || j > N)) {
    return {-1, -1, -1, false};
  }

  // 3. Tile 1
  if (tile == 1) {
    if (i < 1) return {4, N, j, true};           // West -> Tile 4 East
    if (i > N) return {2, 1, j, true};           // East -> Tile 2 West
    if (j < 1) return {6, i, N, true};           // South -> Tile 6 North
    if (j > N) return {5, i, 1, true};           // North -> Tile 5 South
  }

  // 4. Tile 2
  if (tile == 2) {
    if (i < 1) return {1, N, j, true};           // West -> Tile 1 East
    if (i > N) return {3, 1, j, true};           // East -> Tile 3 West
    if (j < 1) return {6, N, N - i + 1, true};   // South -> Tile 6 East (+90 deg)
    if (j > N) return {5, N, i, true};           // North -> Tile 5 East (-90 deg)
  }

  // 5. Tile 3
  if (tile == 3) {
    if (i < 1) return {2, N, j, true};           // West -> Tile 2 East
    if (i > N) return {4, 1, j, true};           // East -> Tile 4 West
    if (j < 1) return {6, N - i + 1, 1, true};   // South -> Tile 6 South (180 deg)
    if (j > N) return {5, N - i + 1, N, true};   // North -> Tile 5 North (180 deg)
  }

  // 6. Tile 4
  if (tile == 4) {
    if (i < 1) return {3, N, j, true};           // West -> Tile 3 East
    if (i > N) return {1, 1, j, true};           // East -> Tile 1 West
    if (j < 1) return {6, 1, i, true};           // South -> Tile 6 West (-90 deg)
    if (j > N) return {5, 1, N - i + 1, true};   // North -> Tile 5 West (+90 deg)
  }

  // 7. Tile 5 (North Pole)
  if (tile == 5) {
    if (i < 1) return {4, N - j + 1, N, true};   // West -> Tile 4 North
    if (i > N) return {2, j, N, true};           // East -> Tile 2 North
    if (j < 1) return {1, i, N, true};           // South -> Tile 1 North
    if (j > N) return {3, N - i + 1, N, true};   // North -> Tile 3 North
  }

  // 8. Tile 6 (South Pole)
  if (tile == 6) {
    if (i < 1) return {4, j, 1, true};           // West -> Tile 4 South
    if (i > N) return {2, N - j + 1, 1, true};   // East -> Tile 2 South
    if (j < 1) return {3, N - i + 1, 1, true};   // South -> Tile 3 South
    if (j > N) return {1, i, 1, true};           // North -> Tile 1 South
  }

  return {-1, -1, -1, false};
}

// -----------------------------------------------------------------------------
// Helper: Construct pure quads across panels using FV3 seam neighbor resolution
// -----------------------------------------------------------------------------
void compute_quad_nodes(
    size_t ownedCount_,
    int N,
    const std::vector<int>& rank_tiles,
    const std::vector<int>& rank_i,
    const std::vector<int>& rank_j,
    const std::vector<atlas::gidx_t>& globalIdx,
    std::vector<std::array<atlas::gidx_t, 4>>& quadNodes,
    std::vector<atlas::gidx_t>& quadGidx,
    const eckit::mpi::Comm& comm)
{
  quadNodes.clear();
  quadGidx.clear();

  // 1. Build local node table: (tile, i, j) -> 1-based globalIdx
  std::unordered_map<uint64_t, atlas::gidx_t> nodeLookup;
  nodeLookup.reserve(globalIdx.size());

  for (size_t k = 0; k < globalIdx.size(); ++k) {
    nodeLookup[pack_tile_ij(rank_tiles[k], rank_i[k], rank_j[k])] = globalIdx[k];
  }

  // 2. Identify full quads where this rank owns the SW corner node.
  // Standard counter-clockwise orientation:
  // SW: (t, i, j), SE: (t, i+1, j), NE: (t, i+1, j+1), NW: (t, i, j+1)
  for (size_t n = 0; n < ownedCount_; ++n) {
    int t = rank_tiles[n];
    int i = rank_i[n];
    int j = rank_j[n];

    atlas::gidx_t id_SW = globalIdx[n];

    auto [t_SE, i_SE, j_SE, v_SE] = resolve_fv3_seam_neighbor(t, i + 1, j, N);
    auto [t_NE, i_NE, j_NE, v_NE] = resolve_fv3_seam_neighbor(t, i + 1, j + 1, N);
    auto [t_NW, i_NW, j_NW, v_NW] = resolve_fv3_seam_neighbor(t, i, j + 1, N);

    if (!v_SE || !v_NE || !v_NW) continue;

    auto it_SE = nodeLookup.find(pack_tile_ij(t_SE, i_SE, j_SE));
    if (it_SE == nodeLookup.end()) continue;

    auto it_NE = nodeLookup.find(pack_tile_ij(t_NE, i_NE, j_NE));
    if (it_NE == nodeLookup.end()) continue;

    auto it_NW = nodeLookup.find(pack_tile_ij(t_NW, i_NW, j_NW));
    if (it_NW == nodeLookup.end()) continue;

    quadNodes.push_back({id_SW, it_SE->second, it_NE->second, it_NW->second});
  }

  // 3. Compute non-overlapping global cell indices across all MPI partitions
  const int nQ = static_cast<int>(quadNodes.size());
  std::vector<int> nQ_perRank(comm.size(), 0);
  comm.allGather(nQ, nQ_perRank.begin(), nQ_perRank.end());

  int quadOffset = 0;
  for (size_t r = 0; r < comm.rank(); ++r) {
    quadOffset += nQ_perRank[r];
  }

  quadGidx.resize(nQ);
  for (int q = 0; q < nQ; ++q) {
    quadGidx[q] = static_cast<atlas::gidx_t>(quadOffset + q + 1);
  }
}


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
  geomVariables = sharedData_->geomVariables_;
  levelsAreTopDown = sharedData_->levelsAreTopDown_;
  numLevels = sharedData_->numLevels_;
  //numberLevels_ = sharedData_->numLevels_;

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
  levelsAreTopDown = config.getBool("levels_are_topdown", false);
  levelsAreTopDown_ = levelsAreTopDown;

  // 1. Read coordinates (lat, lon, elevation) from configuration or input file
  std::string latlonFile = config.getString("latlon_file", "");
  std::string latName = config.getString("lat_name", "latitude");
  std::string lonName = config.getString("lon_name", "longitude");
  std::string elevName = config.getString("elev_name", "elevation");

  std::vector<double> lats;
  std::vector<double> lons; 
  std::vector<double> elevations;
  std::vector<int> global_tile;
  std::vector<int> global_i;
  std::vector<int> global_j;

  if(!latlonFile.empty() && config.has("lat_name") && config.has("lon_name") && config.has("elev_name")) {
    oops::Log::info() << "Reading lat/lon/elevation from file: " << latlonFile << std::endl;
    readLatLonElevFromFile(latlonFile, latName, lonName, elevName, lats, lons, elevations, global_tile, global_i, global_j);
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
  const size_t nproc = static_cast<size_t>(comm_.size());
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
  // Warning: only if also sorting other coord vars--tile, i, j
//***** doesn't seem to make a difference---sortLandPointsByHilbert(global_pts);
//TODO: Sort corresponding elevation 
 

  size_t points_per_rank = num_points / nproc;
  size_t remainder = num_points % nproc;
  
  size_t start_idx = my_rank * points_per_rank + std::min<size_t>(my_rank, remainder);
  size_t local_count = points_per_rank + (my_rank < remainder ? 1 : 0);

  eckit::LocalConfiguration fs_config;
  fs_config.set("mpi_comm", comm_.name());
  //fs_config.set("halo_radius", R_halo);

  /*grid_ = atlas::UnstructuredGrid(global_pts);
  atlas::grid::Partitioner partitioner("checkerboard", fs_config);
  const atlas::grid::Distribution distribution(grid_, partitioner);
  //const atlas::UnstructuredGrid ugrid(global_pts);
  //const atlas::grid::Distribution distribution = atlas::grid::Partitioner("checkerboard", nproc).partition(ugrid);

  std::vector<int> rank_owning(num_points, 0);
  for (int i = 0; i < num_points; ++i){
    rank_owning[i] = distribution.partition(i);
  }*/
  
  
  std::vector<int> rank_owning(num_points, 0);
  // 2. Assign balanced chunks across MPI ranks
  for (size_t i = 0; i < num_points; ++i) {
    rank_owning[i] = static_cast<int>((i * nproc) / num_points);
  }
  
  // 3. Optional: Wrap in an explicit Atlas distribution if downstream APIs require one
  std::vector<int> part_list(rank_owning.begin(), rank_owning.end());
  atlas::grid::Distribution distribution(comm_.size(), num_points, part_list.data());

  std::vector<atlas::PointXY> rank_points;
  std::vector<int> rank_tiles;
  std::vector<int> rank_i;
  std::vector<int> rank_j;

  std::vector<double> rank_lons;
  std::vector<double> rank_lats;
  std::vector<int> ghosts;  //(ownedCount_,0); // 0 for owned, 1 for ghost
  std::vector<atlas::gidx_t> globalIdx; // Original 1-based global index
  std::vector<atlas::idx_t> remoteIdx;      // use 1-based
  std::vector<int> partitions;

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

	rank_tiles.push_back(global_tile[g_idx]);
        rank_i.push_back(global_i[g_idx]);
        rank_j.push_back(global_j[g_idx]);
        partitions.push_back(my_rank);
    }
  }
  
  ownedCount_ = rank_points.size(); 
  std::vector<size_t> num_owned_vec(nproc);
  comm_.allGather(ownedCount_, num_owned_vec.begin(), num_owned_vec.end());
  auto min_no = *std::min_element(num_owned_vec.begin(), num_owned_vec.end());
  auto max_no = *std::max_element(num_owned_vec.begin(), num_owned_vec.end());
  oops::Log::info() << "proc " << my_rank << " ownedCount_ " << ownedCount_ << " min/max ownedCount_ " 
	  << min_no << " / " << max_no << std::endl;  

  size_t min_owned = ownedCount_;
  size_t max_owned = ownedCount_;
  comm_.allReduceInPlace(min_owned, eckit::mpi::Operation::MIN);
  comm_.allReduceInPlace(max_owned, eckit::mpi::Operation::MAX);

  oops::Log::info() << "proc " << my_rank << " min/max owned " << min_owned << " / " << max_owned << std::endl;

  /*// 3. Determine panel resolution N across all ranks
  int local_N = 0;
  for (size_t n = 0; n < ownedCount_; ++n) {
    local_N = std::max(local_N, std::max(rank_i[n], rank_j[n]));
  }
  int N = local_N;
  comm.allReduceInPlace(N, eckit::mpi::Operation::MAX);*/

  const int N = config.getInt("resolution", 384);

  // 4. Build global spatial index to discover owners of missing boundary vertices
  std::unordered_map<uint64_t, size_t> globalCoordToIdx;
  globalCoordToIdx.reserve(num_points);
  for (size_t g = 0; g < num_points; ++g) {
    globalCoordToIdx[pack_tile_ij(global_tile[g], global_i[g], global_j[g])] = g;
  }

  // 5. Determine neighbor nodes required to close boundary quads across FV3 seams
  std::vector<std::vector<atlas::gidx_t>> requests(nproc);
  const std::array<std::pair<int, int>, 3> offsets = {{{1, 0}, {1, 1}, {0, 1}}};

  for (size_t n = 0; n < ownedCount_; ++n) {
    int t = rank_tiles[n];
    int i = rank_i[n];
    int j = rank_j[n];

    for (const auto& off : offsets) {
      auto [target_t, target_i, target_j, valid] =
          resolve_fv3_seam_neighbor(t, i + off.first, j + off.second, N);

      if (!valid) continue;

      uint64_t key = pack_tile_ij(target_t, target_i, target_j);
      auto it = globalCoordToIdx.find(key);
      if (it != globalCoordToIdx.end()) {
        size_t neighbor_g_idx = it->second;
        int owner = rank_owning[neighbor_g_idx];

        if (owner != my_rank) {
          atlas::gidx_t req_gid = static_cast<atlas::gidx_t>(neighbor_g_idx + 1);
          requests[owner].push_back(req_gid);
        }
      }
    }
  }

  // Deduplicate request lists per destination rank
  for (int r = 0; r < nproc; ++r) {
    std::sort(requests[r].begin(), requests[r].end());
    requests[r].erase(std::unique(requests[r].begin(), requests[r].end()), requests[r].end());
  }

  // 6. Non-blocking MPI halo exchange (Phase 1: Exchange node ID requests)
  std::vector<int> send_req_counts(nproc, 0);
  for (int r = 0; r < nproc; ++r) {
    send_req_counts[r] = static_cast<int>(requests[r].size());
  }

  std::vector<int> recv_req_counts(nproc, 0);
  comm_.allToAll(send_req_counts, recv_req_counts);

  std::vector<std::vector<atlas::gidx_t>> incoming_requests(nproc);
  std::vector<eckit::mpi::Request> req_stage1;

  for (int r = 0; r < nproc; ++r) {
    if (r != my_rank && recv_req_counts[r] > 0) {
      incoming_requests[r].resize(recv_req_counts[r]);
      req_stage1.push_back(comm_.iReceive(incoming_requests[r].data(),
                                         recv_req_counts[r],
                                         r, /*tag=*/100));
    }
  }

  for (int r = 0; r < nproc; ++r) {
    if (r != my_rank && send_req_counts[r] > 0) {
      req_stage1.push_back(comm_.iSend(requests[r].data(),
                                      send_req_counts[r],
                                      r, /*tag=*/100));
    }
  }

  comm_.waitAll(req_stage1);

  // 7. Non-blocking MPI halo exchange (Phase 2: Pack and send ghost node payloads)
  std::vector<std::vector<GhostNodeData>> outgoing_ghosts(nproc);
  for (int r = 0; r < nproc; ++r) {
    if (r == my_rank) continue;

    outgoing_ghosts[r].reserve(recv_req_counts[r]);
    for (atlas::gidx_t gid : incoming_requests[r]) {
      size_t g_idx = static_cast<size_t>(gid - 1);
      GhostNodeData node{};  //{} zero-initialize
      node.globalIdx = gid;
      node.tile      = global_tile[g_idx];
      node.i         = global_i[g_idx];
      node.j         = global_j[g_idx];
      node.lon       = global_pts[g_idx].x();
      node.lat       = global_pts[g_idx].y();
      node.remoteIdx = static_cast<atlas::idx_t>(ownerLocalIdx[g_idx] + 1);
      node.ownerRank = my_rank;
      outgoing_ghosts[r].push_back(node);
    }
  }

  std::vector<std::vector<GhostNodeData>> incoming_ghosts(nproc);
  std::vector<eckit::mpi::Request> req_stage2;

  for (int r = 0; r < nproc; ++r) {
    if (r != my_rank && send_req_counts[r] > 0) {
      incoming_ghosts[r].resize(send_req_counts[r]);
      req_stage2.push_back(comm_.iReceive(
          reinterpret_cast<char*>(incoming_ghosts[r].data()),
          static_cast<int>(send_req_counts[r] * sizeof(GhostNodeData)),
          r, /*tag=*/200));
    }
  }

  for (int r = 0; r < nproc; ++r) {
    if (r != my_rank && recv_req_counts[r] > 0) {
      req_stage2.push_back(comm_.iSend(
          reinterpret_cast<const char*>(outgoing_ghosts[r].data()),
          static_cast<int>(outgoing_ghosts[r].size() * sizeof(GhostNodeData)),
          r, /*tag=*/200));
    }
  }

  comm_.waitAll(req_stage2);

  // 8. Append received ghost/halo nodes to local structures
  // Keep track of all global indices already registered locally (both owned and ghosts)
  std::unordered_set<atlas::gidx_t> registeredGidx(globalIdx.begin(), globalIdx.end());
  for (int r = 0; r < nproc; ++r) {
    for (const auto& g : incoming_ghosts[r]) {
      // If this ghost node was already added by another quad or incoming message, SKIP IT
      if (!registeredGidx.insert(g.globalIdx).second) {
        continue;
      }
      rank_points.push_back(atlas::PointXY(g.lon, g.lat));
      rank_lons.push_back(g.lon);
      rank_lats.push_back(g.lat);
      rank_tiles.push_back(g.tile);
      rank_i.push_back(g.i);
      rank_j.push_back(g.j);

      ghosts.push_back(1); // 1 = Ghost node
      globalIdx.push_back(g.globalIdx);
      remoteIdx.push_back(g.remoteIdx);
      partitions.push_back(g.ownerRank);
    }
  }

  // 9. Build pure quads across local and halo nodes using seam transformations
  std::vector<std::array<atlas::gidx_t, 4>> quadNodes;
  std::vector<atlas::gidx_t> quadGidx;

  compute_quad_nodes(
      ownedCount_,
      N,
      rank_tiles, rank_i, rank_j,
      globalIdx,
      quadNodes, quadGidx,
      comm_);

  // 1. Verify owned node count matches the expected input exactly across all ranks
  size_t local_owned = 0;
  for (int g : ghosts) {
    if (g == 0) local_owned++;
  }
  
  size_t total_owned = local_owned;
  comm_.allReduceInPlace(total_owned, eckit::mpi::Operation::SUM);
  
  oops::Log::info() << "Total globally owned nodes (sum of ghosts == 0): "
                    << total_owned
                    << " (Expected: " << num_points << ")" << std::endl;
  
  // 2. Check maximum globalIdx assigned across all ranks
  atlas::gidx_t max_gid = 0;
  for (atlas::gidx_t gid : globalIdx) {
    max_gid = std::max(max_gid, gid);
  }
  comm.allReduceInPlace(max_gid, eckit::mpi::Operation::MAX);
  
  oops::Log::info() << "Maximum globalIdx across all ranks: "
                  << max_gid << std::endl;

  size_t nghost = std::accumulate(ghosts.begin(), ghosts.end(), 0);
  std::vector<size_t> num_ghost_vec(nproc);
  comm_.allGather(nghost, num_ghost_vec.begin(), num_ghost_vec.end());
  auto min_ng = *std::min_element(num_ghost_vec.begin(), num_ghost_vec.end());
  auto max_ng = *std::max_element(num_ghost_vec.begin(), num_ghost_vec.end());
  oops::Log::info() << "proc 0 nghost " << nghost << " min/max numGhost "
          << min_ng << " / " << max_ng << std::endl;
  
    // 1. Check for duplicate globalIdx in local rank tables
  std::unordered_set<atlas::gidx_t> seen_gidx;
  size_t local_duplicates = 0;
  for (atlas::gidx_t gid : globalIdx) {
    if (!seen_gidx.insert(gid).second) {
      local_duplicates++;
    }
  }
  
  size_t total_duplicates = local_duplicates;
  comm_.allReduceInPlace(total_duplicates, eckit::mpi::Operation::SUM);
  oops::Log::info() << "Total duplicate node registrations across local vectors: " 
                    << total_duplicates << std::endl;
  
  // 2. Validate ghost remoteIdx against owner capacity
  // Share each rank's owned node vector size across all ranks
  std::vector<int> localOwnedCounts(nproc, 0);
  int myOwned = static_cast<int>(ownedCount_);
  comm_.allGather(myOwned, localOwnedCounts.begin(), localOwnedCounts.end());
  
  size_t out_of_bounds_remoteIdx = 0;
  for (size_t k = 0; k < ghosts.size(); ++k) {
    if (ghosts[k] == 1) {
      int owner = partitions[k];
      atlas::idx_t rem = remoteIdx[k];
      // remoteIdx is 1-based, so it must be between 1 and localOwnedCounts[owner]
      if (rem < 1 || rem > localOwnedCounts[owner]) {
        out_of_bounds_remoteIdx++;
      }
    }
  }
  
  size_t total_bad_remotes = out_of_bounds_remoteIdx;
  comm_.allReduceInPlace(total_bad_remotes, eckit::mpi::Operation::SUM);
  oops::Log::info() << "Total out-of-bounds ghost remoteIdx values: " 
                  << total_bad_remotes << std::endl;

  const size_t total_nodes = rank_points.size();

  ASSERT(rank_lons.size()   == total_nodes);
  ASSERT(rank_lats.size()   == total_nodes);
  ASSERT(ghosts.size()      == total_nodes);
  ASSERT(globalIdx.size()   == total_nodes);
  ASSERT(remoteIdx.size()   == total_nodes);
  ASSERT(partitions.size()  == total_nodes);
  
  // Check for uninitialized / invalid values
  for (size_t k = 0; k < total_nodes; ++k) {
    ASSERT(globalIdx[k] > 0);
    ASSERT(partitions[k] >= 0 && partitions[k] < comm_.size());
    ASSERT(ghosts[k] == 0 || ghosts[k] == 1);
    if (ghosts[k] == 0) {
      ASSERT(partitions[k] == comm_.rank());
    } else {
      ASSERT(partitions[k] != comm_.rank());
    }
  }

  // 10. Build Atlas Mesh with zero triangle structures
  eckit::LocalConfiguration meshConf;
  meshConf.set("mpi_comm", comm.name());
  // Inform Atlas that parallel fields are already validly constructed
  meshConf.set("halo", 0);

  const atlas::mesh::MeshBuilder meshBuilder;
  mesh_ = meshBuilder(
      rank_lons, rank_lats, ghosts,
      globalIdx, remoteIdx, /*remote_index_base=*/1, partitions,
      /*triag_nodes_global=*/{}, /*triag_global_index=*/{},
      /*quad_nodes_global=*/quadNodes, /*quad_global_index=*/quadGidx,
      meshConf);

  // 11. Global verification that at least one rank holds cells
  size_t local_cells = mesh_.cells().size();
  size_t global_cells = local_cells;
  comm_.allReduceInPlace(global_cells, eckit::mpi::Operation::SUM);
  ASSERT(global_cells > 0);

  // Guarantee cells are present so MeshTriangulation populates its KD-tree
  ASSERT(mesh_.cells().size() > 0);

  // 12. Instantiate Atlas NodeColumns function space
  functionSpace_ = atlas::functionspace::NodeColumns(mesh_, atlas::util::Config("mpi_comm", comm.name()));  //fs_config);

  lonlat_ = functionSpace_.lonlat();
  ghost_ = functionSpace_.ghost();
  auto glb_idx = functionSpace_.global_index();
  auto r_idx = functionSpace_.remote_index();
  partition_ = functionSpace_.partition();

  auto lonlatView = atlas::array::make_view<double, 2>(lonlat_);
  auto glbidxView = atlas::array::make_view<atlas::gidx_t, 1>(glb_idx);
  auto partView = atlas::array::make_view<int,1>(partition_);
  auto ridxView = atlas::array::make_view<int, 1>(r_idx);
  
  auto ghostView = atlas::array::make_view<int, 1>(ghost_);

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
    if (partView(i) == my_rank){  //Populate only owned nodes from the local input arrays
      atlas::gidx_t g_idx = glbidxView(i) - 1;; // 0-based lookup into global elevation array
      elevView(i, 0) = elevations[g_idx];  //TODO: use sorted elevation
      latView(i, 0) = lonlatView(i, 1);
      lonView(i,0) = lonlatView(i, 0);
      global_index_view(i) = glbidxView(i);
      remote_index_view(i) = ridxView(i);
      
      /*const double lon = lons[g_idx];
      const double lat = lats[g_idx];
      lonView(i, 0)    = lon;
      latView(i, 0)    = lat;
      lonlatView(i, 0) = lon;
      lonlatView(i, 1) = lat;

    } else {  //Clear / zero ghost slots before the exchange to ensure no stale data remains
      lonView(i, 0)    = 0.0;
      latView(i, 0)    = 0.0;
      lonlatView(i, 0) = 0.0;
      lonlatView(i, 1) = 0.0;*/
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

  // 1. Ensure ghost field matches your manual ghost array
  for (atlas::idx_t j = 0; j < functionSpace_.size(); ++j) {
    int g_expected = (j < ownedCount_) ? 0 : 1;
    ASSERT(ghostView(j) == g_expected);
  }
  
  // 2. Add gmask explicitly for BUMP universe routines
  atlas::Field gmaskField = functionSpace_.createField<int>(
      atlas::option::name("gmask") | atlas::option::levels(1));
  auto gmaskView = atlas::array::make_view<int, 2>(gmaskField);
  for (atlas::idx_t j = 0; j < functionSpace_.size(); ++j) {
    gmaskView(j, 0) = (j < ownedCount_) ? 1 : 0;
  }
 
  // 2. Add mask explicitly for BUMP universe routines
  atlas::Field maskField = functionSpace_.createField<int>(
      atlas::option::name("mask") | atlas::option::levels(1));
  auto maskView = atlas::array::make_view<int, 2>(maskField);
  for (atlas::idx_t j = 0; j < functionSpace_.size(); ++j) {
    maskView(j, 0) = (j < ownedCount_) ? 1 : 0;
  }
  

  // do halo exchange
  atlas::functionspace::NodeColumns fsHalo(functionSpace_);

  fsHalo.haloExchange(elevationField_);
  fsHalo.haloExchange(lat_);
  fsHalo.haloExchange(lon_);
  fsHalo.haloExchange(lonlat_);
  //fsHalo.haloExchange(owned_); do not fieldexchange owned
  fsHalo.haloExchange(area_);
  fsHalo.haloExchange(global_index_);
  fsHalo.haloExchange(remote_index_);
  fsHalo.haloExchange(gmaskField);
  fsHalo.haloExchange(maskField);

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
  fields_.add(gmaskField);
  fields_.add(maskField);

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

  //atlas::functionspace::NodeColumns fsHalo(functionSpace_);
  fsHalo.haloExchange(loc_latlon);
  functionSpace_.gather(loc_latlon, global_latlon);

  oops::Log::info() << ">>> GEOMETRY GATHER TEST RESULT <<<" << std::endl;
  oops::Log::info() << "Global gathered shape: " << global_latlon.shape(0) 
                    << " x " << global_latlon.shape(1) << std::endl;
  if (comm_.rank() == 0) {
    ASSERT(global_latlon.shape(0) > 0); //, "Atlas functionSpace.gather() returned 0 points! GeometryData cannot build globalNodeTree_.");
  }

  // Gather all owned coordinates to rank 0 to test for global duplicates
  std::vector<double> local_owned_lons;
  for (size_t n = 0; n < ownedCount_; ++n) 
    local_owned_lons.push_back(rank_lons[n]);

  // Count owned points across all ranks
  std::vector<size_t> recvcounts(comm_.size());
  comm_.allGather(local_owned_lons.size(), recvcounts.begin(), recvcounts.end());
  
  // Inspect whether the sum of owned points equals the global active count exactly
  size_t total_rec = 0;
  for (size_t c : recvcounts) total_rec += c;
  oops::Log::info() << "Total global owned points gathered: " 
                    << total_rec << " (Target: " << glbcnt << ")" << std::endl;
  ASSERT(total_rec == glbcnt);


// -----------------------------------------------------------------------------
// Coordinate Audit Using eckit::mpi gatherv
// -----------------------------------------------------------------------------
{
  const size_t nLocal = functionSpace_.size();
  auto lonlatView = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
  auto ghostView  = atlas::array::make_view<int, 1>(functionSpace_.ghost());

  const size_t fieldsPerNode = 4;
  std::vector<double> localBuffer(nLocal * fieldsPerNode);

  for (atlas::idx_t j = 0; j < static_cast<atlas::idx_t>(nLocal); ++j) {
    localBuffer[j * fieldsPerNode + 0] = lonlatView(j, 0);
    localBuffer[j * fieldsPerNode + 1] = lonlatView(j, 1);
    localBuffer[j * fieldsPerNode + 2] = static_cast<double>(ghostView(j));
    localBuffer[j * fieldsPerNode + 3] = static_cast<double>(j);
  }

  const int nRanks = comm_.size();
  const int myRank = comm_.rank();
  const size_t rootRank = 0;

  // eckit MPI wrapper expects int for counts and displs
  std::vector<int> counts(nRanks, 0);
  int localCount = static_cast<int>(localBuffer.size());
  
  // Gather counts from all ranks to root rank
  comm_.gather(localCount, counts, rootRank);

  std::vector<int> displs(nRanks, 0);
  int totalElements = 0;
  if (myRank == 0) {
    for (int r = 0; r < nRanks; ++r) {
      displs[r] = totalElements;
      totalElements += counts[r];
    }
  }

  std::vector<double> globalBuffer(totalElements);

  // eckit::mpi::Comm::gatherv signature:
  // gatherv(const double* sendbuf, size_t sendcount, double* recvbuf,
  //         const int* recvcounts, const int* displs, size_t root)
  comm_.gatherv(localBuffer.data(), localBuffer.size(),
                globalBuffer.data(), counts.data(), displs.data(), rootRank);

  if (myRank == 0) {
    struct NodeMeta {
      int rank;
      int localIdx;
      int isGhost;
      double lon;
      double lat;
    };

    auto quantize = [](double val) -> int64_t {
      return static_cast<int64_t>(std::round(val * 1e6));
    };

    std::map<std::pair<int64_t, int64_t>, std::vector<NodeMeta>> coordMap;
    size_t totalOwnedNodes = 0;
    size_t totalGhostNodes = 0;
    int offset = 0;

    for (int r = 0; r < nRanks; ++r) {
      const int nodeCount = counts[r] / fieldsPerNode;
      for (int n = 0; n < nodeCount; ++n) {
        const double lon   = globalBuffer[offset + n * fieldsPerNode + 0];
        const double lat   = globalBuffer[offset + n * fieldsPerNode + 1];
        const int isGhost  = static_cast<int>(globalBuffer[offset + n * fieldsPerNode + 2]);
        const int localIdx = static_cast<int>(globalBuffer[offset + n * fieldsPerNode + 3]);

        if (isGhost == 0) totalOwnedNodes++;
        else totalGhostNodes++;

        auto key = std::make_pair(quantize(lon), quantize(lat));
        coordMap[key].push_back({r, localIdx, isGhost, lon, lat});
      }
      offset += counts[r];
    }

    oops::Log::info() << "==================================================" << std::endl;
    oops::Log::info() << " GLOBAL COORDINATE AUDIT REPORT" << std::endl;
    oops::Log::info() << " Total Nodes Gathered: " << (totalElements / fieldsPerNode) << std::endl;
    oops::Log::info() << " Total Owned Nodes:    " << totalOwnedNodes << " (Target: 272699)" << std::endl;
    oops::Log::info() << " Total Ghost Nodes:    " << totalGhostNodes << std::endl;
    oops::Log::info() << " Unique Coordinate Bins: " << coordMap.size() << std::endl;

    size_t ownedOwnedCollisions = 0;
    size_t ownedGhostMatches = 0;

    for (const auto &entry : coordMap) {
      const auto &nodes = entry.second;
      if (nodes.size() > 1) {
        int ownedCountAtCoord = 0;
        for (const auto &node : nodes) {
          if (node.isGhost == 0) ownedCountAtCoord++;
        }

        if (ownedCountAtCoord > 1) {
          ownedOwnedCollisions++;
          oops::Log::error() << ">>> [OWNED-OWNED COLLISION] at lon: " 
                             << std::setprecision(10) << std::fixed
                             << nodes[0].lon << ", lat: " << nodes[0].lat << std::endl;
          for (const auto &node : nodes) {
            oops::Log::error() << "    -> Rank " << node.rank
                               << " | LocalIdx: " << node.localIdx 
                               << " | Ghost: " << node.isGhost << std::endl;
          }
        } else if (ownedCountAtCoord == 1 && nodes.size() > 1) {
          ownedGhostMatches++;
        }
      }
    }

    oops::Log::info() << " Normal Owned-to-Ghost Seam Matches: " << ownedGhostMatches << std::endl;
    oops::Log::info() << " Fatal Owned-Owned Collisions:      " << ownedOwnedCollisions << std::endl;
  }
  

  size_t northPolePoints = 0;
  size_t southPolePoints = 0;
  for (size_t n = 0; n < ownedCount_; ++n) {
    if (std::abs(lats[n] - 90.0) < 1e-4) northPolePoints++;
    if (std::abs(lats[n] + 90.0) < 1e-4) southPolePoints++;
  }
  size_t globalNorth = northPolePoints, globalSouth = southPolePoints;
  comm_.allReduceInPlace(globalNorth, eckit::mpi::Operation::SUM);
  comm_.allReduceInPlace(globalSouth, eckit::mpi::Operation::SUM);
  
  oops::Log::info() << "Points at North Pole (90 deg): " << globalNorth
                    << " | Points at South Pole (-90 deg): " << globalSouth << std::endl;
}

// -----------------------------------------------------------------------------
// Intercept BUMP Universe Coordinates on Rank 0
// -----------------------------------------------------------------------------
{
  const int myRank = comm_.rank();
  const int nRanks = comm_.size();
  const size_t nLocal = functionSpace_.size();

  auto lonlatView = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
  auto ghostView  = atlas::array::make_view<int, 1>(functionSpace_.ghost());

  // Package only what BUMP would inspect (active/owned points, or all points)
  // Let's inspect all points that have ghost == 0 vs all points in functionSpace
  std::vector<double> localCoords;
  for (atlas::idx_t j = 0; j < static_cast<atlas::idx_t>(nLocal); ++j) {
    // If BUMP ignores ghost, this will catch cross-rank duplicates:
    localCoords.push_back(lonlatView(j, 0)); // lon
    localCoords.push_back(lonlatView(j, 1)); // lat
    localCoords.push_back(static_cast<double>(ghostView(j))); // 0=owned, 1=ghost
    localCoords.push_back(static_cast<double>(j));            // local index
    localCoords.push_back(static_cast<double>(myRank));       // owning rank
  }

  const int fieldsPerPoint = 5;
  int localCount = static_cast<int>(localCoords.size());
  std::vector<int> counts(nRanks, 0);
  comm_.gather(localCount, counts, 0);

  std::vector<int> displs(nRanks, 0);
  int totalElements = 0;
  if (myRank == 0) {
    for (int r = 0; r < nRanks; ++r) {
      displs[r] = totalElements;
      totalElements += counts[r];
    }
  }

  std::vector<double> globalCoords(totalElements);
  comm_.gatherv(localCoords.data(), localCoords.size(),
                globalCoords.data(), counts.data(), displs.data(), 0);

  if (myRank == 0) {
    const int totalPoints = totalElements / fieldsPerPoint;
    oops::Log::info() << ">>> AUDITING " << totalPoints << " COORDINATES FOR BUMP UNIVERSE DUPLICATES <<<" << std::endl;

    // Convert degrees to radians on the unit sphere for great-circle chord distance
    const double deg2rad = M_PI / 180.0;
    struct Point3D {
      double x, y, z;
      double lon, lat;
      int isGhost;
      int localIdx;
      int rank;
    };

    std::vector<Point3D> pts(totalPoints);
    int poleCount = 0;

    for (int i = 0; i < totalPoints; ++i) {
      double lon = globalCoords[i * fieldsPerPoint + 0];
      double lat = globalCoords[i * fieldsPerPoint + 1];
      int ghost  = static_cast<int>(globalCoords[i * fieldsPerPoint + 2]);
      int lIdx   = static_cast<int>(globalCoords[i * fieldsPerPoint + 3]);
      int rk     = static_cast<int>(globalCoords[i * fieldsPerPoint + 4]);

      // Check for singular poles
      if (std::abs(std::abs(lat) - 90.0) < 1e-5) {
        poleCount++;
      }

      double rlon = lon * deg2rad;
      double rlat = lat * deg2rad;
      pts[i] = {
        std::cos(rlat) * std::cos(rlon),
        std::cos(rlat) * std::sin(rlon),
        std::sin(rlat),
        lon, lat, ghost, lIdx, rk
      };
    }

    oops::Log::info() << ">>> Polar points (|lat| = 90 deg): " << poleCount << std::endl;

    // Check pairwise distance using BUMP's zero-distance threshold (chord distance squared < 1e-12)
    const double eps2 = 1.0e-12;
    int duplicateCount = 0;

    for (int i = 0; i < totalPoints; ++i) {
      for (int j = i + 1; j < totalPoints; ++j) {
        // Only evaluate if both are considered or if testing owned vs owned
        double dx = pts[i].x - pts[j].x;
        double dy = pts[i].y - pts[j].y;
        double dz = pts[i].z - pts[j].z;
        double d2 = dx*dx + dy*dy + dz*dz;

        if (d2 < eps2) {
          duplicateCount++;
          if (duplicateCount <= 10) { // Print first 10 collisions
            oops::Log::error() << ">>> BUMP COLLISION #" << duplicateCount << " (dist^2 = " << d2 << ") <<<" << std::endl;
            oops::Log::error() << "  Point A: Rank=" << pts[i].rank << " LocalIdx=" << pts[i].localIdx
                               << " Ghost=" << pts[i].isGhost << " Lon=" << pts[i].lon << " Lat=" << pts[i].lat << std::endl;
            oops::Log::error() << "  Point B: Rank=" << pts[j].rank << " LocalIdx=" << pts[j].localIdx
                               << " Ghost=" << pts[j].isGhost << " Lon=" << pts[j].lon << " Lat=" << pts[j].lat << std::endl;
          }
        }
      }
    }
    oops::Log::info() << ">>> Total Identical Coordinate Collisions on Sphere: " << duplicateCount << std::endl;
  }
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
  /*size_t localSize = functionSpace_.size();
  size_t globalSize = 0;
  // Sum local sizes across all MPI ranks in comm_
  comm_.allReduce(localSize, globalSize, eckit::mpi::Operation::SUM);

  return globalSize;*/
  atlas::functionspace::NodeColumns fsHalo(functionSpace_);
  return fsHalo.nb_nodes_global();
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
    std::vector<double> & lats, std::vector<double> & lons, std::vector<double> & elevations, 
    std::vector<int>& global_tile, std::vector<int>& global_i, std::vector<int>& global_j)
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
    global_tile.resize(dimlen);
    global_i.resize(dimlen);
    global_j.resize(dimlen);

    nc_rc(nc_get_var_double(fileId, varId, lats.data()), "nc_get_var_double " + latName);

    nc_rc(nc_inq_varid(fileId, lonName.c_str(), &varId), "nc_inq_varid " + lonName);
    nc_rc(nc_get_var_double(fileId, varId, lons.data()), "nc_get_var_double " + lonName);

    nc_rc(nc_inq_varid(fileId, elevName.c_str(), &varId), "nc_inq_varid " + elevName);
    nc_rc(nc_get_var_double(fileId, varId, elevations.data()), "nc_get_var_double " + elevName);
    
    nc_rc(nc_inq_varid(fileId, "cube_tile", &varId), "nc_inq_varid cube_tile ");  // + elevName);
    nc_rc(nc_get_var_int(fileId, varId, global_tile.data()), "nc_get_var_int cube_tile "); // + elevName);

    nc_rc(nc_inq_varid(fileId, "cube_i", &varId), "nc_inq_varid cube_i" );  
    nc_rc(nc_get_var_int(fileId, varId, global_i.data()), "nc_get_var_int cube_j ");

    nc_rc(nc_inq_varid(fileId, "cube_j", &varId), "nc_inq_varid cube_j ");
    nc_rc(nc_get_var_int(fileId, varId, global_j.data()), "nc_get_var_int cube_j");

    // Close file
    nc_rc(nc_close(fileId), "nc_close "+ pathFile);

    oops::Log::info() << "Done reading file " << pathFile << std::endl;

}

// -----------------------------------------------------------------------------------------------

}  // namespace ijedi
