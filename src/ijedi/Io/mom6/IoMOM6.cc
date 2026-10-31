#include <ostream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Geometry/mom6/GeometryMOM6.h"
#include "ijedi/Io/mom6/FillMaskedCells.h"
#include "ijedi/Io/mom6/IoMOM6.h"
#include "ijedi/Io/mom6/ReadMOM6Netcdf.h"
#include "ijedi/Io/mom6/VelocityStagger.h"
#include "ijedi/Io/mom6/WriteMOM6Netcdf.h"

namespace ijedi
{
  namespace
  {
    // True if x holds any eastward/northward velocity component.
    bool hasVelocity(const atlas::FieldSet & x)
    {
      for (const auto & pair : mom6::velocityPairs())
        if (x.has(pair.first) || x.has(pair.second)) return true;
      return false;
    }

    // Gather the geometry rotation and mask onto the full structured grid on
    // root (collective).  Points outside the JEDI mesh are land: mask 0, no
    // rotation.
    mom6::VelocityGridInfo gatherVelocityGrid(const Geometry & geom)
    {
      const auto modelData = geom.modelData();
      mom6::VelocityGridInfo info;
      info.ni = modelData.getInt("ni");
      info.nj = modelData.getInt("nj");
      info.nz = modelData.getInt("nz");
      info.hasFold = modelData.getBool("has northern fold");

      const atlas::FieldSet & gf = geom.fields();
      const auto fs = atlas::functionspace::NodeColumns(gf.field("cos_rot").functionspace());
      auto localGidx = atlas::array::make_view<atlas::gidx_t, 1>(fs.nodes().global_index());
      atlas::Field gidx = fs.createField<double>(atlas::option::name("__gidx")
                                                 | atlas::option::levels(1));
      {
        auto v = atlas::array::make_view<double, 2>(gidx);
        for (atlas::idx_t n = 0; n < fs.nodes().size(); ++n)
          v(n, 0) = static_cast<double>(localGidx(n));
      }
      auto gather = [&](const atlas::Field & f) {
        atlas::Field g = fs.createField<double>(atlas::option::name(f.name())
                                                | atlas::option::levels(f.shape(1))
                                                | atlas::option::global());
        fs.gather(f, g);
        return g;
      };
      const atlas::Field gIdx = gather(gidx);
      const atlas::Field gCos = gather(gf.field("cos_rot"));
      const atlas::Field gSin = gather(gf.field("sin_rot"));
      const atlas::Field gMask = gather(gf.field("mask3d"));

      if (geom.getComm().rank() == 0) {
        const size_t nh = static_cast<size_t>(info.nj) * info.ni;
        info.cosRot.assign(nh, 1.0);
        info.sinRot.assign(nh, 0.0);
        info.mask3d.assign(static_cast<size_t>(info.nz) * nh, 0.0);
        auto vIdx = atlas::array::make_view<double, 2>(gIdx);
        auto vCos = atlas::array::make_view<double, 2>(gCos);
        auto vSin = atlas::array::make_view<double, 2>(gSin);
        auto vMask = atlas::array::make_view<double, 2>(gMask);
        for (atlas::idx_t n = 0; n < gIdx.shape(0); ++n) {
          const size_t flat = static_cast<size_t>(vIdx(n, 0)) - 1;
          info.cosRot[flat] = vCos(n, 0);
          info.sinRot[flat] = vSin(n, 0);
          for (int k = 0; k < info.nz; ++k)
            info.mask3d[static_cast<size_t>(k) * nh + flat] = vMask(n, k);
        }
      }
      return info;
    }
  }  // namespace

  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoMOM6> makerIoMOM6_("mom6");
  // -------------------------------------------------------------------------------------------------
  IoMOM6::IoMOM6(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()),
        geom_(geom)
  {
    util::Timer timer(classname(), "IoMOM6");
    oops::Log::trace() << classname() << " constructor starting" << std::endl;

    const std::string datapath = params.datapath.value();

    writeFilepath_ = datapath + "/" + params.ocn_file.value();
    readFilepaths_.push_back(writeFilepath_);

    const auto addFile = [&](const boost::optional<std::string> & opt) {
      if (opt) readFilepaths_.push_back(datapath + "/" + *opt);
    };
    addFile(params.ice_file.value());
    addFile(params.fix_file.value());

    const std::string velocityGrid = params.velocity_grid.value();
    if (velocityGrid != "tracer" && velocityGrid != "native")
      throw eckit::UserError(classname() + ": velocity grid must be tracer or native, got "
                             + velocityGrid, Here());
    nativeVelocity_ = (velocityGrid == "native");

    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoMOM6::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                    const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read state");
    oops::Log::trace() << classname() << " read state starting" << std::endl;

    // Resolve variable names and scalings for all fields
    std::vector<std::string> fileVarNames;
    std::vector<double> scalings;
    for (const auto & field : x) {
      const std::string jediName = field.name();
      fileVarNames.push_back(fileionames.has(jediName)
                             ? fileionames.getString(jediName) : "");
      scalings.push_back(fileioscaling.has(jediName)
                         ? fileioscaling.getDouble(jediName) : 0.0);
    }

    // Delegate to the NetCDF reader.  On a MOM6 geometry every task reads only
    // its MOM6 block and scatters it to its JEDI points; otherwise root reads
    // whole variables and scatters them.
    const auto * mom6Geom = dynamic_cast<const GeometryMOM6 *>(&geom_.geometryImpl());
    if (mom6Geom != nullptr) {
      readMOM6NetcdfBlocks(readFilepaths_, x, fileVarNames, scalings, *mom6Geom);
    } else if (hasVelocity(x)) {
      const mom6::VelocityGridInfo velGrid = gatherVelocityGrid(geom_);
      readMOM6Netcdf(readFilepaths_, x, fileVarNames, scalings, geom_.getComm(), &velGrid);
    } else {
      readMOM6Netcdf(readFilepaths_, x, fileVarNames, scalings, geom_.getComm());
    }

    // Apply boundary conditions to masked (land) cells: nearest-neighbour
    // extrapolation for tracers, zero for non-tracers.
    applyBoundaryConditions(x, geom_.fields().field("mask3d"),
                          geom_.getFieldMetadata());

    oops::Log::trace() << classname() << " read state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoMOM6::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                     const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;

    // Get geometry info
    const auto modelData = geom_.modelData();
    const int ni = modelData.getInt("ni");
    const int nj = modelData.getInt("nj");
    const int nz = modelData.getInt("nz");

    // Resolve variable names for all fields
    std::vector<std::string> fileVarNames;
    for (const auto & field : x) {
      const std::string jediName = field.name();
      fileVarNames.push_back(fileionames.has(jediName)
                             ? fileionames.getString(jediName) : "");
    }

    // Delegate to the NetCDF writer.  On a MOM6 geometry every task writes its
    // own MOM6 block (parallel NetCDF-4); otherwise root gathers and writes.
    const auto * mom6Geom = dynamic_cast<const GeometryMOM6 *>(&geom_.geometryImpl());
    if (mom6Geom != nullptr) {
      writeMOM6NetcdfBlocks(writeFilepath_, x, fileVarNames, nz, *mom6Geom, nativeVelocity_);
    } else if (nativeVelocity_ && hasVelocity(x)) {
      const mom6::VelocityGridInfo velGrid = gatherVelocityGrid(geom_);
      writeMOM6Netcdf(writeFilepath_, x, fileVarNames, ni, nj, nz, geom_.getComm(),
                      true, &velGrid);
    } else {
      writeMOM6Netcdf(writeFilepath_, x, fileVarNames, ni, nj, nz, geom_.getComm());
    }

    oops::Log::trace() << classname() << " write state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoMOM6::print(std::ostream &os) const
  {
    os << classname() << " Io for MOM6 ocean, sea ice, and fix files";
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi
