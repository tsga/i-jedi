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

  const size_t my_rank = static_cast<size_t>(comm_.rank());
  const size_t nproc = static_cast<size_t>(comm.size());
  oops::Log::info() << "Proc  " << my_rank << " running with " << nproc << " procs" << std::endl;

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
  sortLandPointsByHilbert(global_pts);
//TODO: Sort corresponding elevation 
 

  size_t points_per_rank = num_points / nproc;
  size_t remainder = num_points % nproc;
  
  size_t start_idx = my_rank * points_per_rank + std::min<size_t>(my_rank, remainder);
  size_t local_count = points_per_rank + (my_rank < remainder ? 1 : 0);

  eckit::LocalConfiguration fs_config;
  fs_config.set("mpi_comm", comm_.name());
  fs_config.set("halo_radius", R_halo);

  grid_ = atlas::UnstructuredGrid(global_pts);
  atlas::grid::Partitioner partitioner("equal_regions", fs_config);
  atlas::grid::Distribution distribution(grid_, partitioner);

  // assign contiguous chunks to ranks
  // using custom partitioner
  std::vector<int> cust_partition(num_points);
  size_t current_idx = 0;
  for (int rank = 0; rank < nproc; ++rank) {
      size_t count = points_per_rank + (rank < remainder ? 1 : 0);
      for (size_t i = 0; i < count; ++i) {
          cust_partition[current_idx++] = rank;
      }
  }
  // using 'equal regions' partitioner
  std::vector<int> rank_owning(num_points);
  for (int i = 0; i < num_points; ++i){
    rank_owning[i] = distribution.partition(i);
  }

  if (cust_partition != rank_owning) {
    
    oops::Log::warning() << "Warning! cust_partition not= rank_owning from partition " << std::endl;
/*    auto it1 = cust_partition.begin();
    auto it2 = rank_owning.begin();

    oops::Log::info() << "Mismatched indices: " << std::endl;
    while (true) {
        // Find the next mismatch
        auto p = std::mismatch(it1, cust_partition.end(), it2);
        
        // If it reached the end of v1, break
        if (p.first == cust_partition.end()) break;

        // Calculate and print the index using std::distance
        oops::Log::info() << std::distance(cust_partition.begin(), p.first) << " ";

        // Move past the current mismatch to keep searching
        it1 = p.first + 1;
        it2 = p.second + 1;
    }
    oops::Log::info() << "\n";
    std::string errMsg = "GeometryLandVector: cust_partition != rank_owning";
    oops::Log::error() << errMsg << std::endl;
    throw eckit::Exception(errMsg, Here());
*/
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
          /*// 1-based local index on owner_rank for ghost
          // (Assuming rank offset calculation based on points_per_rank and remainder)
          int owner_rank = rank_owning[g_idx];  // cust_partition[g_idx];
	  size_t start_idx = owner_rank * points_per_rank + std::min<size_t>(owner_rank, remainder);
          size_t local_idx_on_owner = g_idx - start_idx + 1;
	  remoteIdx.push_back(local_idx_on_owner);*/
          
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
  
  eckit::LocalConfiguration meshConf;
  meshConf.set("mpi_comm", comm_.name());

  const atlas::mesh::MeshBuilder meshBuilder;
  mesh_ = meshBuilder(rank_lons, rank_lats, ghosts,
      globalIdx, remoteIdx, /*remote_index_base=*/1, partitions,
      /*tri_boundary_nodes=*/{}, /*tri_global_indices=*/{}, /*quadNodes*/{}, /*quadGidx*/{},
      meshConf);
  atlas::mesh::actions::build_halo(mesh_, 1);

  functionSpace_ = atlas::functionspace::NodeColumns(mesh_, fs_config);
  //functionSpace_ = atlas::functionspace::PointCloud(grid_, partitioner, fs_config);

  lonlat_ = functionSpace_.lonlat();
  ghost_ = functionSpace_.ghost();
  global_index_ = functionSpace_.global_index();
  remote_index_ = functionSpace_.remote_index();
  partition_ = functionSpace_.partition();

  auto lonlatView = atlas::array::make_view<double, 2>(lonlat_);
  auto gidxView = atlas::array::make_view<atlas::gidx_t, 1>(global_index_);
  auto partView = atlas::array::make_view<int,1>(partition_);
 
  // Explicit 1D latitude and longitude fields, eg for KD-tree search
  lat_ = functionSpace_.createField<double>(atlas::option::name("lat") | atlas::option::levels(1));
  lon_ = functionSpace_.createField<double>(atlas::option::name("lon") | atlas::option::levels(1));
  auto latView = atlas::array::make_view<double, 2>(lat_);
  auto lonView = atlas::array::make_view<double, 2>(lon_);

  elevationField_ = functionSpace_.createField<double>(atlas::option::name("height") | atlas::option::levels(1));
  auto elevView   = atlas::array::make_view<double, 2>(elevationField_);
  
  for (size_t i = 0; i < functionSpace_.size(); ++i){
    if (partView(i) == my_rank){
      atlas::gidx_t g_idx = gidxView(i) - 1;; // 0-based lookup into global elevation array
      elevView(i, 0) = elevations[g_idx];  //TODO: use sorted elevation
      latView(i, 0) = lonlatView(i, 1);
      lonView(i,0) = lonlatView(i, 0);
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
  auto loccnt=localNodeCount();
  auto glbcnt=globalNodeCount();
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

  atlas::Field global_latlon = functionSpace_.createField<double>(atlas::option::name("global_latlon") 
		  | atlas::option::variables(2) | atlas::option::global());
  atlas::Field loc_latlon = functionSpace_.createField<double>(
      atlas::option::name("loc_latlon") | atlas::option::variables(2));
  auto locllView = atlas::array::make_view<double, 2>(loc_latlon);
  for (size_t i = 0; i < functionSpace_.size(); ++i) {
    locllView(i, 0) = llView(i, 0);
    locllView(i, 1) = llView(i, 0);;
  }

  atlas::functionspace::PointCloud fsHalo(functionSpace_);
  fsHalo.haloExchange(loc_latlon);
  fsHalo.haloExchange(global_latlon);
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
