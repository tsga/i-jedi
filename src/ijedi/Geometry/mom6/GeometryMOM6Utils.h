#pragma once

#include <cstddef>
#include <string>
#include <unordered_set>
#include <vector>

namespace ijedi {
namespace mom6 {

int computeExtent(int N, int ndivs, int pe);

/// Rectangular block of the effective (possibly coarsened) T grid:
/// columns [i0, i0+ni) and rows [j0, j0+nj), 0-based.  Block readers return
/// row-major arrays of size ni*nj (times numLevels for 3D, level-major).
struct Block {
  int i0 = 0, j0 = 0, ni = 0, nj = 0;
  size_t size() const { return static_cast<size_t>(ni) * nj; }
};

/// T/U/V coordinates and T-cell metrics of one block, read from the
/// supergrid sub-rectangle covering it (no global arrays).
void readHgridBlock(const std::string & path,
                    int coarsenFactor,
                    const Block & block,
                    std::vector<double> * lon,
                    std::vector<double> * lat,
                    std::vector<double> * dxT,
                    std::vector<double> * dyT,
                    std::vector<double> * areaT,
                    std::vector<double> * lonU,
                    std::vector<double> * latU,
                    std::vector<double> * lonV,
                    std::vector<double> * latV,
                    std::vector<double> * cosRot,
                    std::vector<double> * sinRot);

/// Global T-point lon/lat (niEff*njEff), read a strip of supergrid rows at a
/// time.  Same values as readHgrid's lon/lat.
void readHgridLonLat(const std::string & path,
                     int niGlobal,
                     int niEff,
                     int njEff,
                     int coarsenFactor,
                     std::vector<double> * lon,
                     std::vector<double> * lat);

/// Depth and wet mask of one block (coarsened like readTopog).
void readTopogBlock(const std::string & path,
                    int coarsenFactor,
                    double minimumDepth,
                    const Block & block,
                    std::vector<double> * depth,
                    std::vector<double> * wet);

/// Layer thickness 'h' and layer-centre depth of one block (coarsened like
/// readVerticalGeometry).
void readVerticalGeometryBlock(const std::string & path,
                               int niGlobal,
                               int njGlobal,
                               int numLevels,
                               int coarsenFactor,
                               const Block & block,
                               std::vector<double> * layerThickness,
                               std::vector<double> * layerCenterDepth);

void rewriteLeveledNodeDataAsTime(const std::string & filename,
                                  const std::unordered_set<std::string> & timeFields);

}  // namespace mom6
}  // namespace ijedi
