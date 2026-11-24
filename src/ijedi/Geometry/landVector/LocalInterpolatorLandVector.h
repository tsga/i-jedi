#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/interpolation/Interpolation.h"
#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ijedi/Geometry/landVector/GeometryLandVector.h"
// Include i-JEDI state and increment headers
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/State/State.h"
#include "oops/generic/SourceProximityPartitioner.h"
#include "oops/base/GeometryData.h"

// Helper to normalize any angle into [0, 360) degrees
inline double normalizeLon(double lon) {
  double norm = std::fmod(lon, 360.0);
  if (norm < 0.0) norm += 360.0;
  return norm;
}

namespace ijedi {

// -----------------------------------------------------------------------------
// Target Partitioner Wrapper exposing interpolatingTask() expected by OOPS
// -----------------------------------------------------------------------------
/*class LandVectorTargetPartitioner {
 public:
  explicit LandVectorTargetPartitioner(const Geometry & geom)
      : partitioner_("equal_regions") {}

  // Method called by OOPS line 217 to assign target obs locations to MPI ranks
  int interpolatingTask(double lat, double lon) const {
    atlas::Point2 pt(lon, lat);
    return partitioner_.partition(geo.functionspace().grid()); //pt);
  }

 private:
  atlas::grid::Partitioner partitioner_;
};*/	


class LocalInterpolatorLandVector {
 public:
  // ---------------------------------------------------------------------------
  // 1. Static Partitioner Factory required by oops::LocalInterpolator
  // ---------------------------------------------------------------------------
  /*static auto makeTargetPartitioner(const Geometry & geom) { 
    // Default to Atlas EqualRegions or matching partitioner for point cloud distribution
    //return atlas::grid::Partitioner("equal_regions");
    return LandVectorTargetPartitioner(geom);
  }*/

  // Partition target points based on proximity to a source grid node.
  static oops::SourceProximityPartitioner makeTargetPartitioner(const Geometry & geom) {
    return oops::makeSourceProximityPartitioner(geom.functionSpace(), geom.comm());
  }

  // ---------------------------------------------------------------------------
  // 2. Static Preprocess Hooks (No-ops required by OOPS interface branches)
  // ---------------------------------------------------------------------------
  static void preprocess(const State &) {}
  static void preprocess(State &) {}
  static void preprocess(const Increment &) {}
  static void preprocess(Increment &) {}
  static void preprocess(atlas::FieldSet &) {}
  static void preprocessAD(Increment &) {}
  static void preprocessAD(atlas::FieldSet &) {}

  // ---------------------------------------------------------------------------
  // 3. Constructor: Pre-computes KNN Atlas Interpolation Matrix
  // ---------------------------------------------------------------------------
  LocalInterpolatorLandVector(const eckit::Configuration & conf,
                              const Geometry & geom,
                              const std::vector<double> & target_lats,
                              const std::vector<double> & target_lons)
      : numTargets_(target_lats.size()) {
    // Force KNN interpolation parameters directly
    atlas::util::Config interpConfig;

    std::string interpType; // = interpConf.getString("type", "k-nearest-neighbours");
    int kNN; // = interpConf.getInt("k-nearest-neighbours", 1);
    double maxDist; // = interpConf.getDouble("max_distance", 0.0);
    double bufferDeg; //nterpConf.getDouble("bounding box buffer deg", 0.5);


    eckit::LocalConfiguration interpConf;
    if (conf.has("interpolation method")) {   //local interpolator")) {
      // Read settings from interpConf (now correctly unwrapped!)
      interpConf = conf.getSubConfiguration("interpolation method");     //local interpolator");
      interpType = interpConf.getString("type", "k-nearest-neighbours");
      kNN = interpConf.getInt("k-nearest-neighbours", 1);
      maxDist = interpConf.getDouble("max_distance", 0.0);
      bufferDeg = interpConf.getDouble("bounding box buffer deg", 0.5);
    }
    else {
      interpType = conf.getString("type", "k-nearest-neighbours");
      kNN = conf.getInt("k-nearest-neighbours", 1);
      maxDist = conf.getDouble("max_distance", 0.0);
      bufferDeg = conf.getDouble("bounding box buffer deg", 0.5);
    }

    interpConfig.set("type", interpType); 
    interpConfig.set("k-nearest-neighbours", kNN);
    interpConfig.set("max_distance", maxDist);
    
    //if (geom.getComm().rank() == 0){
    //oops::Log::info() << "LocalInterpolator configuration: " <<  conf << std::endl;
    //oops::Log::info() << "interp config used: " <<  interpConfig << std::endl;
    //}

    if (numTargets_ == 0) return;

    // 1. Compute Local Bounding Box with a safety buffer (~0.5 deg halo)
    double minLat = 90.; //std::numeric_limits<double>::max();
    double maxLat = -90.; //std::numeric_limits<double>::lowest();
    //double minLon = 360.; //std::numeric_limits<double>::max();
    //double maxLon = -360.; //std::numeric_limits<double>::lowest();
    bool hasPoints = false;

    // Extract lon/lat coordinates from local Geometry functionSpace
    auto ghostView = atlas::array::make_view<int, 1>(geom.functionSpace().ghost());
    auto lonlatView = atlas::array::make_view<double, 2>(geom.functionSpace().lonlat());

    for (size_t i = 0; i < geom.functionSpace().size(); ++i) {
      // Include both owned and halo/ghost points for accurate local extent
      //double lon = lonlatView(i, 0);
      double lat = lonlatView(i, 1);

      //minLon = std::min(minLon, lon);
      //maxLon = std::max(maxLon, lon);
      minLat = std::min(minLat, lat);
      maxLat = std::max(maxLat, lat);
      hasPoints = true;
    }
    
    // Add search halo margin (~50km / 0.5 degrees)
    minLat = std::min(-90., minLat - bufferDeg);
    maxLat = std::max(90., maxLat + bufferDeg);

    // Wrap target coordinates into a temporary PointCloud FunctionSpace
    std::vector<atlas::PointXY> pts;
    pts.reserve(numTargets_);
    localTargetIndices_.reserve(numTargets_);
    
    for (size_t i = 0; i < numTargets_; ++i) {
      double targetLat = target_lats[i];
    
      // 1. Check latitude first (cheapest check)
      if (targetLat < minLat || targetLat > maxLat) continue;
    
      double targetLon = normalizeLon(target_lons[i]);
      bool isNearDomain = false;
    
      // 2. Longitude Check: Check distance to ANY local geometry node within radius
      // For LandVector grids, checking angular distance to local nodes prevents boundary errors
      for (size_t j = 0; j < geom.functionSpace().size(); ++j) {
        double nodeLon = normalizeLon(lonlatView(j, 0));
        double nodeLat = lonlatView(j, 1);
    
        // Compute direct lat/lon deltas considering wrapping
        double dLat = std::abs(targetLat - nodeLat);
        double dLon = std::abs(targetLon - nodeLon);
        if (dLon > 180.0) dLon = 360.0 - dLon;
    
        // Fast Manhattan distance check against buffer threshold
        if (dLat <= bufferDeg && dLon <= bufferDeg) {
          isNearDomain = true;
          break;
        }
      }
    
      if (isNearDomain) {
        pts.emplace_back(target_lons[i], target_lats[i]);
        localTargetIndices_.push_back(i);
      }
    }

    pts.shrink_to_fit();
    localTargetIndices_.shrink_to_fit();

    numFilteredTargets_ = pts.size();
    if (numFilteredTargets_ == 0) return;

    targetFs_ = atlas::functionspace::PointCloud(pts);

    // Build KD-tree matrix W once upon initialization
    interpolator_ = atlas::Interpolation(interpConfig,
                                                         geom.functionSpace(),
                                                         targetFs_);
  }

  ~LocalInterpolatorLandVector() = default;

 void print(std::ostream & os) const {
    os << "LocalInterpolatorLandVector: numTargets = " << numTargets_
       << ", type = PointCloud KNN";
  }
  
  
  // ---------------------------------------------------------------------------
  // 4. State & Increment Overloads for apply()
  // ---------------------------------------------------------------------------
  void apply(const oops::Variables & vars,
             const State & state,
             const std::vector<bool> & mask,
             std::vector<double> & buffer) const {
    applyFieldSet(vars, state.fieldSet(), mask, buffer);
  }

  void apply(const oops::Variables & vars,
             const Increment & dx,
             const std::vector<bool> & mask,
             std::vector<double> & buffer) const {
    applyFieldSet(vars, dx.fieldSet(), mask, buffer);
  }
 
 // ---------------------------------------------------------------------------
  // Adjoint Apply: Map flat std::vector<double> buffer back to Increment FieldSet
  // ---------------------------------------------------------------------------
  void applyAD(const oops::Variables & vars,
               const State & state,
               const std::vector<bool> & mask,
               const std::vector<double> & buffer) const {
    applyFieldSetAD(vars, state.fieldSet(), mask, buffer);
  }

  void applyAD(const oops::Variables & vars,
               const Increment & dx,
               const std::vector<bool> & mask,
               const std::vector<double> & buffer) const {
    applyFieldSetAD(vars, dx.fieldSet(), mask, buffer);
  }
  

 private:
  void applyFieldSetAD(const oops::Variables & vars,
                       const atlas::FieldSet & adjFieldSet,
                       const std::vector<bool> & mask,
                       const std::vector<double> & buffer) const 
  {
    
    eckit::Log::warning() << "LocalInterpolatorLandVector::applyAD not implemented" << std::endl;	  
    /*if (numTargets_ == 0) return;

    // 1. Allocate an adjoint target FieldSet matching requested variables
    atlas::FieldSet targetFieldSet;
    for (const auto & varName : vars.variables()) {
      if (adjFieldSet.has(varName)) {
        const auto & srcField = adjFieldSet[varName];
        size_t levels = srcField.levels() > 0 ? srcField.levels() : 1;

        auto targetField = targetFs_.createField<double>(
            atlas::option::name(varName) | atlas::option::levels(levels));

        // Zero-initialize the target field buffer before accumulation
        auto view1D = atlas::array::make_view<double, 1>(targetField);
        view1D.assign(0.0);

        targetFieldSet.add(targetField);
      }
    }

    // 2. Unpack flat std::vector<double> buffer into targetFieldSet
    size_t bufferIdx = 0;
    for (const auto & varName : vars.variables()) {
      if (!targetFieldSet.has(varName)) continue;

      auto & field = targetFieldSet[varName];
      size_t levels = field.levels() > 0 ? field.levels() : 1;

      if (field.rank() == 1) {
        auto view = atlas::array::make_view<double, 1>(field);
        for (size_t i = 0; i < numTargets_; ++i) {
          if (mask.empty() || mask[i]) {
            view(i) = buffer[bufferIdx++];
          }
        }
      } else if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        for (size_t lev = 0; lev < levels; ++lev) {
          for (size_t i = 0; i < numTargets_; ++i) {
            if (mask.empty() || mask[i]) {
              view(i, lev) = buffer[bufferIdx++];
            }
          }
        }
      }
    }

    // 3. Execute Adjoint Interpolation: adjFieldSet += W^T * targetFieldSet
    //interpolator_.execute_adjoint(adjFieldSet, targetFieldSet);
    // 3. Execute Adjoint via Matrix Transpose
    // Access internal matrix W and apply W^T to targetFieldSet -> adjFieldSet
    const auto & matrix = interpolator_->matrix();

    for (const auto & varName : vars.variables()) 
    {
      if (!adjFieldSet.has(varName) || !targetFieldSet.has(varName)) continue;

      auto & srcField = adjFieldSet[varName];
      const auto & tgtField = targetFieldSet[varName];
      size_t levels = srcField.levels() > 0 ? srcField.levels() : 1;

      if (srcField.rank() == 1) {
        auto srcView = atlas::array::make_view<double, 1>(srcField);
        auto tgtView = atlas::array::make_view<double, 1>(tgtField);

        // Accumulate transpose matrix product: srcView += W^T * tgtView
        matrix.execute_transpose(tgtView, srcView);
      } else if (srcField.rank() == 2) {
        for (size_t lev = 0; lev < levels; ++lev) {
          auto srcView = atlas::array::make_view<double, 2>(srcField);
          auto tgtView = atlas::array::make_view<double, 2>(tgtField);

          // Slice per level for 2D fields
          matrix.execute_transpose(tgtView.slice(atlas::array::Range::all(), lev),
                                  srcView.slice(atlas::array::Range::all(), lev));
        }
      }
    }*/
  
  }
  
  void applyFieldSet(const oops::Variables & vars,
                     const atlas::FieldSet & srcFieldSet,
                     const std::vector<bool> & mask,
                     std::vector<double> & buffer) const {
    // Safety Guard 1: Early return if no targets exist on this MPI task
    if (numTargets_ == 0 || numFilteredTargets_ == 0) return;

    atlas::FieldSet cleanSrcFieldSet;
    atlas::FieldSet targetFieldSet;
    for (const auto & varName : vars.variables()) {
      if (srcFieldSet.has(varName)) {
        const auto & srcField = srcFieldSet[varName];
	cleanSrcFieldSet.add(srcField);

        size_t levels = srcField.levels() > 0 ? srcField.levels() : 1;

        auto targetField = targetFs_.createField<double>(
            atlas::option::name(varName) | atlas::option::levels(levels));
        targetFieldSet.add(targetField);
      }
    }

    // Execute matrix product on filtered subset
    //interpolator_.execute(srcFieldSet, targetFieldSet);
    if (targetFieldSet.size() > 0) {
      interpolator_.execute(cleanSrcFieldSet, targetFieldSet);
    }
    // Safety Guard 2: Sequentially unpack using flat buffer index tracking
    size_t bufIdx = 0;
    for (const auto & varName : vars.variables()) {
      if (!targetFieldSet.has(varName)) continue;

      const auto & field = targetFieldSet[varName];
      size_t levels = field.levels() > 0 ? field.levels() : 1;

      if (field.rank() == 1) {
        auto view = atlas::array::make_view<double, 1>(field);

        // Zero out entire buffer block for this variable first
        for (size_t i = 0; i < numTargets_; ++i) {
          if (bufIdx + i < buffer.size()) buffer[bufIdx + i] = 0.0;
        }

        // Fill filtered predictions at their exact original indices
        for (size_t k = 0; k < numFilteredTargets_; ++k) {
          size_t origIdx = localTargetIndices_[k];
          if ((mask.empty() || mask[origIdx]) && (bufIdx + origIdx < buffer.size())) {
            buffer[bufIdx + origIdx] = view(k);
          }
        }
        bufIdx += numTargets_;

      } else if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        for (size_t lev = 0; lev < levels; ++lev) {

          // Zero out buffer block for this level
          for (size_t i = 0; i < numTargets_; ++i) {
            if (bufIdx + i < buffer.size()) buffer[bufIdx + i] = 0.0;
          }

          // Fill filtered predictions
          for (size_t k = 0; k < numFilteredTargets_; ++k) {
            size_t origIdx = localTargetIndices_[k];
            if ((mask.empty() || mask[origIdx]) && (bufIdx + origIdx < buffer.size())) {
              buffer[bufIdx + origIdx] = view(k, lev);
            }
          }
          bufIdx += numTargets_;
        }
      }
    }
  }

  size_t numTargets_{0};
  size_t numFilteredTargets_{0};
  std::vector<size_t> localTargetIndices_;
  atlas::FunctionSpace targetFs_;
  atlas::Interpolation interpolator_;
};
 // Stream insertion operator expected by OOPS interface
inline std::ostream & operator<<(std::ostream & os, const LocalInterpolatorLandVector & interpolator) {
  interpolator.print(os);
  return os;
}

}  // namespace ijedi
