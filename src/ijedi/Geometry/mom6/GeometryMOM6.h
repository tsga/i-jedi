#pragma once

#include <cstdint>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "ijedi/Geometry/base/GeometryBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class GeometryMOM6 : public GeometryBase
  {
   public:
    GeometryMOM6(const eckit::Configuration &, const eckit::mpi::Comm &,
                 eckit::LocalConfiguration &, atlas::FunctionSpace &,
                 atlas::FieldSet &, bool &, int &);
    void print(std::ostream &) const override;
    std::vector<double> verticalCoord(std::string &) const override;
    void addModelVaderIngredients(const atlas::FieldSet & geomFields,
                                  atlas::FieldSet & fset, int nlevels) const override;

    // Accessors for MOM6-specific (structured) decomposition.
    // Used by ModelMOM6 for state injection — not exposed via GeometryBase.
    const atlas::FunctionSpace & mom6FunctionSpace() const
      { return mom6FunctionSpace_; }
    const atlas::FieldSet      & mom6Fields() const
      { return mom6Fields_; }

    // Scatter/gather map: for each JEDI-owned ocean point, the owning MOM6
    // rank and the row-major local index on that rank.  Built once at init;
    // used by ModelMOM6 to exchange state before/after MOM6::run().
    struct ScatterMap {
      std::vector<int> mom6Rank;     // [0..ownedCount_) → MOM6 rank
      std::vector<int> mom6LocalIdx;  // [0..ownedCount_) → row-major idx
    };
    const ScatterMap & scatterMap() const { return scatterMap_; }

    // Move a field between the MOM6 structured decomposition and the JEDI
    // unstructured one.  ExchangePlan is precomputed once in the constructor;
    // scatter/gather are cheap thereafter.
    void scatterToJedi(const atlas::Field & mom6Field,
                       atlas::Field * jediField) const;
    void gatherFromJedi(const atlas::Field & jediField,
                        atlas::Field * mom6Field) const;

   private:
    // --- Function space builders ---
    // MOM6 structured space: replicates FMS rectangular tile decomposition
    // → mom6FunctionSpace_
    void buildMom6FunctionSpace(const eckit::mpi::Comm & comm,
                                const std::vector<double> & lonGlobal,
                                const std::vector<double> & latGlobal);

    // JEDI unstructured space: ocean + configurable land-fringe points,
    // load-balanced
    // → functionSpace_
    void buildJediFunctionSpace(const eckit::mpi::Comm & comm,
                                const std::vector<double> & wetGlobal,
                                const std::vector<double> & lonGlobal,
                                const std::vector<double> & latGlobal);

    // --- Field builders ---
    // Populate mom6Fields_ on mom6FunctionSpace_ from global arrays
    // (extracts local slice)
    void buildMom6Fields(const std::vector<double> & lonGlobal,
                         const std::vector<double> & latGlobal,
                         const std::vector<double> & depthGlobal,
                         const std::vector<double> & wetGlobal,
                         const std::vector<double> & layerThicknessGlobal,
                         const std::vector<double> & layerCenterDepthGlobal,
                         const std::vector<double> & mask3dGlobal,
                         const std::vector<double> & dxTGlobal,
                         const std::vector<double> & dyTGlobal,
                         const std::vector<double> & areaTGlobal,
                         const std::vector<double> & lonUGlobal,
                         const std::vector<double> & latUGlobal,
                         const std::vector<double> & lonVGlobal,
                         const std::vector<double> & latVGlobal,
                         const std::vector<double> & cosRotGlobal,
                         const std::vector<double> & sinRotGlobal);

    // Populate fields_ on functionSpace_ (JEDI unstructured) from global
    // arrays via jediPoints_
    void buildFields(const std::vector<double> & lonGlobal,
                     const std::vector<double> & latGlobal,
                     const std::vector<double> & depthGlobal,
                     const std::vector<double> & wetGlobal,
                     const std::vector<double> & layerThicknessGlobal,
                     const std::vector<double> & layerCenterDepthGlobal,
                     const std::vector<double> & mask3dGlobal,
                     const std::vector<double> & dxTGlobal,
                     const std::vector<double> & dyTGlobal,
                     const std::vector<double> & areaTGlobal,
                     const std::vector<double> & lonUGlobal,
                     const std::vector<double> & latUGlobal,
                     const std::vector<double> & lonVGlobal,
                     const std::vector<double> & latVGlobal,
                     const std::vector<double> & cosRotGlobal,
                     const std::vector<double> & sinRotGlobal,
                     bool buildVerticalGeometry);

    // Distance from each JEDI point to the nearest land cell (metres)
    void buildDistFromCoast(const std::vector<double> & lonGlobal,
                            const std::vector<double> & latGlobal,
                            const std::vector<double> & wetGlobal,
                            const std::vector<double> & mask3dGlobal);

    // --- Scatter map & exchange plan ---
    void buildScatterMap();
    void buildExchangePlan();
    void checkScatterMap();

    // --- Output helpers ---
    void saveGrid(const std::string & filename,
                  const eckit::mpi::Comm & comm) const;
    void saveStructuredGrid(const std::string & filename,
                            const eckit::mpi::Comm & comm) const;
    void saveDebugMesh(const std::string & prefix,
                       const eckit::mpi::Comm & comm) const;
    void saveGmsh(const std::string & filename,
                  const eckit::mpi::Comm & comm) const;

    // Grid parameters (from YAML config)
    int niGlobal_  = 0;
    int njGlobal_  = 0;
    int numLevels_ = 0;
    int layoutX_   = 1;
    int layoutY_   = 1;
    double minimumDepth_     = 0.0;
    double minimumThickness_ = 1e-6;  // layer thickness below which a cell is masked
    int coarsenFactor_ = 1;
    int niEff_         = 0;   // niGlobal_ / coarsenFactor_
    int njEff_         = 0;   // njGlobal_ / coarsenFactor_
    int fringeWidth_   = 2;   // retained land-fringe layers around ocean
    bool hasFold_      = false;  // true for tripolar grids with a northern fold

    // Local MOM6 compute domain (1-based, inclusive)
    int iStart_ = 1, iCount_ = 0;
    int jStart_ = 1, jCount_ = 0;

    // MOM6 structured function space and fields (local compute domain)
    atlas::FunctionSpace mom6FunctionSpace_;
    atlas::FieldSet      mom6Fields_;

    // JEDI function space and fields
    atlas::FunctionSpace functionSpace_;
    atlas::FieldSet      fields_;

    // JEDI unstructured partition
    int nActiveGlobal_ = 0;   // global active point count (ocean + fringe)
    int ownedCount_    = 0;   // number of active points owned by this rank
    int ownedMin_      = 0;   // min owned count across ranks (load balance)
    int ownedMax_      = 0;   // max owned count across ranks (load balance)

    // (iGlobal, jGlobal) for each JEDI node in functionSpace_ order
    // (owned first, then ghost)
    std::vector<std::pair<int, int>> jediPoints_;

    // Scatter/gather map: JEDI unstructured ↔ MOM6 structured
    ScatterMap scatterMap_;

    // Communication plan precomputed from ScatterMap — drives allToAllv
    // in scatter/gather.
    struct ExchangePlan {
      std::vector<int> recvCounts;    // [npes] vals this rank recvs per MOM6
      std::vector<int> recvDispl;     // [npes] prefix sum of recvCounts
      std::vector<int> recvOrder;     // [ownedCount_] recv-buf → JEDI index
      std::vector<int> sendCounts;    // [npes] vals this rank sends per JEDI
      std::vector<int> sendDispl;     // [npes] prefix sum of sendCounts
      std::vector<int> sendLocalIdx;  // [sum(sendCounts)] mom6FS_ local idx
    };
    ExchangePlan exchangePlan_;

    const eckit::mpi::Comm & comm_;
  };

}  // namespace ijedi
