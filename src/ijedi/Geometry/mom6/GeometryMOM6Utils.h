#pragma once

#include <string>
#include <unordered_set>
#include <vector>

namespace ijedi {
namespace mom6 {

int computeExtent(int N, int ndivs, int pe);

void readHgrid(const std::string & path,
               int niGlobal,
               int njGlobal,
               int niEff,
               int njEff,
               int coarsenFactor,
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

void readTopog(const std::string & path,
               int niGlobal,
               int njGlobal,
               int niEff,
               int njEff,
               int coarsenFactor,
               double minimumDepth,
               std::vector<double> * depth,
               std::vector<double> * wet);

void readVerticalGeometry(const std::string & path,
                         int niGlobal,
                         int njGlobal,
                         int niEff,
                         int njEff,
                         int numLevels,
                         int coarsenFactor,
                         std::vector<double> * layerThickness,
                         std::vector<double> * layerCenterDepth);

void rewriteLeveledNodeDataAsTime(const std::string & filename,
                                  const std::unordered_set<std::string> & timeFields);

}  // namespace mom6
}  // namespace ijedi
