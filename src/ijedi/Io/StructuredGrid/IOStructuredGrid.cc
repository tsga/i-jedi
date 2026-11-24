/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <netcdf.h>

#include <map>
#include <algorithm>
#include <cmath>
#include <ostream>
#include <string>
#include <vector>


#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/stringFunctions.h"

//#include "oops/base/GeometryData.h"

#include "ijedi/FieldMetadata/FieldsMetadata.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/State/State.h"

//#include "fv3jedi/Utilities/fv3jedi_vertical_remap.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/StructuredGrid/IOStructuredGrid.h"

namespace ijedi {
// -------------------------------------------------------------------------------------------------
static IoMaker<IOStructuredGrid> makerIOStructuredGrid_("structured grid");
static IoMaker<IOStructuredGrid> makerIOAuxGrid_("auxgrid");
// -------------------------------------------------------------------------------------------------
static inline void nc_rc(const int return_code, const std::string & operation) {
  if (return_code != NC_NOERR) {
    ABORT("IOStructuredGrid netCDF operation \'" + operation + "\' failed with error: "
          + nc_strerror(return_code));
  }
}

IOStructuredGrid::IOStructuredGrid(const Geometry & geom, const Parameters_ & params) 
  : IoBase(geom, params.toConfiguration()), geom_(geom), params_(params) {
    
    // 1. All ranks participate natively in the parallel layout built during initialization
    // 2. No `writeFunctionSpace_` or `interpolator_` members needed!
    oops::Log::trace() << classname() << " Constructor configured for pure geographic space." << std::endl;
}
// -------------------------------------------------------------------------------------------------
IOStructuredGrid::~IOStructuredGrid() {
  util::Timer timer(classname(), "~IOStructuredGrid");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}

void IOStructuredGrid::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                            const eckit::LocalConfiguration &fileioscaling) const 
{
  util::Timer timer(classname(), "read");
  oops::Log::trace() << classname() << " read started" << std::endl;

  // 1. Get your parallel FunctionSpace handle
  auto readFunctionSpace_ = atlas::functionspace::StructuredColumns(geom_.functionSpace());
  const atlas::Grid & grid = readFunctionSpace_.grid();

  // 2. Create a temporary serial distribution where ALL points live on Rank 0
  std::vector<int> zeros(grid.size(), 0);
  atlas::grid::Distribution serialDist(geom_.getComm().size(), grid.size(), zeros.data());
  
  eckit::LocalConfiguration atlas_conf;
  atlas_conf.set("mpi_comm", geom_.getComm().name());
  atlas::functionspace::StructuredColumns serialFunctionSpace(grid, serialDist, atlas_conf);

  // 3. Allocate fieldsSerial using the serial function space factory
  // This guarantees all ranks have perfectly matching metadata, ranks, and levels
  atlas::FieldSet fieldsSerial;
  for (const auto & field : x) {
    // Pass the matching metadata/options from the target field
    atlas::Field fSerial = serialFunctionSpace.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(field.levels()));
    fieldsSerial.add(fSerial);
  }

  // 4. Loop through and read filenames on Rank 0
  const auto & filenamesOpt = params_.filenames.value();
  if (filenamesOpt != boost::none) {
    for (const auto & filename : filenamesOpt.value()) {
      if (geom_.getComm().rank() == 0) {
        util::DateTime dummyTime; 
        this->readStructuredFields(params_.datapath.value() + "/" + filename, 
                                   fieldsSerial, dummyTime, fileionames, fileioscaling);
      }
    }
  } else {
    ABORT("IOStructuredGrid::read: 'filenames' parameter is missing.");
  }

  // 5. Scatter the data globally
  // Since both function spaces share the exact same grid, this will now succeed seamlessly
  readFunctionSpace_.scatter(fieldsSerial, x);
  
  oops::Log::trace() << classname() << " read done" << std::endl;
}

void IOStructuredGrid::write(const atlas::FieldSet & fieldsGeographic,
                             const eckit::LocalConfiguration & fileionames,
                             const eckit::LocalConfiguration & fileioscaling) const {
  util::Timer timer(classname(), "write");
  oops::Log::trace() << classname() << " write started" << std::endl;

  if (!params_.doVerticalRemapping.value()) {
    auto readFunctionSpace_ = atlas::functionspace::StructuredColumns(geom_.functionSpace());
    const atlas::Grid & grid = readFunctionSpace_.grid();
        
    // Create the temporary serial helper space for gathering data back to Rank 0
    std::vector<int> zeros(grid.size(), 0);
    atlas::grid::Distribution serialDist(geom_.getComm().size(), grid.size(), zeros.data());
    
    eckit::LocalConfiguration atlas_conf;
    atlas_conf.set("mpi_comm", geom_.getComm().name());
    atlas::functionspace::StructuredColumns serialFunctionSpace(grid, serialDist, atlas_conf);

    // Allocate the serial target containers using the factory
    atlas::FieldSet fieldsSerial;
    for (const auto & field : fieldsGeographic) {
      atlas::Field fSerial = serialFunctionSpace.createField<double>(
          atlas::option::name(field.name()) | atlas::option::levels(field.levels()));
      fieldsSerial.add(fSerial);
    }

    // Gather distributed data smoothly from all parallel ranks back into Rank 0
    readFunctionSpace_.gather(fieldsGeographic, fieldsSerial);

    // Resolve the valid time dynamically
    util::DateTime validTime;
    if (fieldsGeographic.metadata().has("time")) {
      validTime = util::DateTime(fieldsGeographic.metadata().get<std::string>("time"));
    } else if (fieldsGeographic.metadata().has("timestamp")) {
      validTime = util::DateTime(fieldsGeographic.metadata().get<std::string>("timestamp"));
    } else if (params_.dateTime.value() != boost::none) {
      validTime = util::DateTime(params_.dateTime.value().value());
    } else {
      oops::Log::warning() << "Using default DateTime placeholder for file writing." << std::endl;
      validTime = util::DateTime();
    }

    // Write to disk exclusively on rank 0
    if (geom_.getComm().rank() == 0) {
      //const util::DateTime dateTime(datTimeString);
      this->writeStructuredFields(fieldsSerial, validTime, fileionames, fileioscaling);
    }
  } else {
    
    ABORT("IOStructuredGrid::write doVerticalRemapping not implemented yet");
    
    /*ASSERT(params_.orographyFilename.value() != boost::none);
    
    // Define orography variables
    atlas::FieldSet fieldsOrog;
    atlas::Field zsOrogNew = fieldsGeographic["geopotential_height_at_surface"].clone();
    fieldsOrog.add(zsOrogNew);

    // Write to disk if rank 0
    if (geom_.getComm().rank() == 0) {
      // Read structured-grid orography from file
      const std::string orogFilename = params_.orographyFilename.value().value();
      this -> readStructuredFields(orogFilename, fieldsOrog, obj.validTime(),
                                   fileionames, fileioscaling);

      // Remap the vertical coordinates to account for orography changes
      fv3jedi::VertRemap vert_remap(geom_, fieldsOrog);
      atlas::FieldSet fieldsGeographicRemap = vert_remap.remap(fieldsGeographic);

      // Write to disk
      this->writeStructuredFields(fieldsGeographicRemap, obj.validTime(),
                                  fileionames, fileioscaling);*/
  }
  
  oops::Log::trace() << classname() << " write done" << std::endl;
}

/*void IOStructuredGrid::write(const atlas::FieldSet & fieldsGeographic,
                             const eckit::LocalConfiguration & fileionames,
                             const eckit::LocalConfiguration & fileioscaling) const {
  util::Timer timer(classname(), "write");
  oops::Log::trace() << classname() << " write started"<< std::endl;

  if (!params_.doVerticalRemapping.value()) {
    const auto & baseSpace = geom_.functionSpace();
    const std::string spaceType = baseSpace.type();
        
    atlas::FieldSet fieldsSerial;

    // 1. Check if the function space is natively structured
    if (spaceType == "StructuredColumns") {
      auto readFunctionSpace = atlas::functionspace::StructuredColumns(baseSpace);
      const atlas::Grid & grid = readFunctionSpace.grid();
          
      std::vector<int> zeros(grid.size(), 0);
      atlas::grid::Distribution serialDist(geom_.getComm().size(), grid.size(), zeros.data());
      
      eckit::LocalConfiguration atlas_conf;
      atlas_conf.set("mpi_comm", geom_.getComm().name());
      atlas::functionspace::StructuredColumns serialFunctionSpace(grid, serialDist, atlas_conf);

      // Allocate the structured serial target containers
      for (const auto & field : fieldsGeographic) {
        atlas::Field fSerial = serialFunctionSpace.createField<double>(
            atlas::option::name(field.name()) | atlas::option::levels(field.levels()));
        fieldsSerial.add(fSerial);
      }

      // Gather distributed data smoothly across all parallel ranks back into Rank 0
      readFunctionSpace.gather(fieldsGeographic, fieldsSerial);
    } 
    else {
      // 2. FALLBACK FOR GAUSSIAN/NODE COLUMNS: Avoid bad_cast completely
      // Allocate global flat placeholder fields on Rank 0 only
      if (geom_.getComm().rank() == 0) {
        for (const auto & field : fieldsGeographic) {
          atlas::Field fSerial = atlas::Field(field.name(), atlas::array::DataType::real64(), 
                                              atlas::array::make_shape(baseSpace.size(), field.levels()));
          fieldsSerial.add(fSerial);
        }
      }

      // Note: If running on a single task (serial), you can map local directly to global here:
      if (geom_.getComm().size() == 1) {
        fieldsSerial = fieldsGeographic; 
      } else {
        // For multi-rank MPI setups using unstructured spaces, the parallel data must be 
        // manually gathered via MPI collectives into Rank 0 before writing.
        oops::Log::warning() << "Gathering on unstructured " << spaceType 
                             << " space requires manual MPI reduction." << std::endl;
      }
    }
    
    // Resolve the valid time dynamically (using your multi-tier fallback logic)
    util::DateTime validTime;
    if (fieldsGeographic.metadata().has("time")) {
      validTime = util::DateTime(fieldsGeographic.metadata().get<std::string>("time"));
    } else if (fieldsGeographic.metadata().has("timestamp")) {
      validTime = util::DateTime(fieldsGeographic.metadata().get<std::string>("timestamp"));
    } else if (params_.dateTime.value() != boost::none) {
      validTime = util::DateTime(params_.dateTime.value().value());
    } else {
      validTime = util::DateTime();
    }

    // 3. Write to disk exclusively on Rank 0
    if (geom_.getComm().rank() == 0) {
      this->writeStructuredFields(fieldsSerial, validTime, fileionames, fileioscaling);
    }
  } else {
    ABORT("IOStructuredGrid::write doVerticalRemapping not implemented yet");
  }
  oops::Log::trace() << classname() << " write done" << std::endl;
}*/

void IOStructuredGrid::writeStructuredFields(const atlas::FieldSet & fields,
                                             const util::DateTime & time,
                                             const eckit::LocalConfiguration & ioNames,
                                             const eckit::LocalConfiguration & ioScaling) const {
  // NetCDF IDs
  // ----------
  int fileId;

  // Dimension indices
  int latId;
  int lonId;
  int levId;
  int edgId;
  int forId;
  int timId;

  // Variable indices
  int fIv;
  std::map<std::string, int> fieldIvs;

/*   // Get ak/bk for writing
  // ---------------------
  std::vector<double> ak = geom_.ak();
  std::vector<double> bk = geom_.bk(); */
  
  // Get the name of the file and adjust with datetime
  // -------------------------------------------------
  std::string pathFile = params_.filename.value();

  // For backward compatibility add some things to the filename if not already present
  if (pathFile.find("%Y") == std::string::npos) {
    pathFile += "%Y%m%d_%H%M%Sz";
  }
  if (pathFile.find(".nc") == std::string::npos) {
    pathFile += ".nc4";
  }

  // Format the datetime string
  pathFile = time.formatString(pathFile);

  // Replace member number (ensemble applciaitons)
  util::stringfunctions::swapNameMember(params_.toConfiguration(), pathFile);

  // Create a file to write fields into
  // ----------------------------------
  nc_rc(nc_create(pathFile.c_str(), NC_CLOBBER | NC_NETCDF4, &fileId), "nc_create" + pathFile);
  oops::Log::warning() << "nc created" << std::endl;
  // Create regular grid for determining lat/lon values
  // --------------------------------------------------
  auto writeFunctionSpace_ = atlas::functionspace::StructuredColumns(geom_.functionSpace());
  const atlas::RegularGrid regGrid(writeFunctionSpace_.grid());
  // Get grid dimensions
  // -------------------
  const int nLat = regGrid.ny();
  const int nLon = regGrid.nx();
  const int nLev = geom_.numLevels();    //writeFunctionSpace_.levels(); 
  const int nEdg = nLev + 1;             //regGrid.npz() + 1;  //geom_.npz();
  const int nFor = 4;
  const int nTim = 1;

  nc_rc(nc_def_dim(fileId, params_.latName.value().c_str(), nLat, &latId), "nc_def_dim (lat)");
  nc_rc(nc_def_dim(fileId, params_.lonName.value().c_str(), nLon, &lonId), "nc_def_dim (lon)");
  nc_rc(nc_def_dim(fileId, params_.levName.value().c_str(), nLev, &levId), "nc_def_dim (lev)");
  nc_rc(nc_def_dim(fileId, params_.edgName.value().c_str(), nEdg, &edgId), "nc_def_dim (edg)");
  nc_rc(nc_def_dim(fileId, params_.forName.value().c_str(), nFor, &forId), "nc_def_dim (for)");
  nc_rc(nc_def_dim(fileId, params_.timName.value().c_str(), nTim, &timId), "nc_def_dim (tim)");

  // Define the dimensions variables in the file
  // -------------------------------------------
  std::vector<double> latArr(nLat);
  std::vector<double> lonArr(nLon);
  std::vector<int> levArr(nLev);
  std::vector<int> edgArr(nEdg);
  std::vector<int> forArr(nFor);
  std::vector<int> timArr(nTim);

  for (int i = 0; i < nLat; ++i) {
    latArr[i] = regGrid.y(nLat - 1 - i);
  }
  for (int i = 0; i < nLon; ++i) {
    lonArr[i] = regGrid.x(i);
  }
  for (int i = 0; i < nLev; ++i) {
    levArr[i] = i + 1;
  }
  for (int i = 0; i < nEdg; ++i) {
    edgArr[i] = i + 1;
  }
  for (int i = 0; i < nFor; ++i) {
    forArr[i] = i + 1;
  }
  for (int i = 0; i < nTim; ++i) {
    timArr[i] = i + 1;
  }

  // Write the dimension variables (and attributes) to the file
  // ----------------------------------------------------------
  nc_rc(nc_def_var(fileId, params_.latName.value().c_str(), NC_DOUBLE, 1, &latId, &fIv),
        "nc_def_var (lat)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("degrees_north"), "degrees_north"),
        "nc_put_att_text (lat)");
  fieldIvs[params_.latName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.lonName.value().c_str(), NC_DOUBLE, 1, &lonId, &fIv),
        "nc_def_var (lon)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("degrees_east"), "degrees_east"),
        "nc_put_att_text (lon)");
  fieldIvs[params_.lonName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.levName.value().c_str(), NC_INT, 1, &levId, &fIv),
        "nc_def_var (lev)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("1"), "1"),
        "nc_put_att_text (lev)");
  fieldIvs[params_.levName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.edgName.value().c_str(), NC_INT, 1, &edgId, &fIv),
        "nc_def_var (edg)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("1"), "1"),
        "nc_put_att_text (edg)");
  fieldIvs[params_.edgName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.forName.value().c_str(), NC_INT, 1, &forId, &fIv),
        "nc_def_var (for)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("1"), "1"),
        "nc_put_att_text (for)");
  fieldIvs[params_.forName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.timName.value().c_str(), NC_INT, 1, &timId, &fIv),
        "nc_def_var (tim)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("1"), "1"),
        "nc_put_att_text (tim)");
  fieldIvs[params_.timName.value()] = fIv;

  // Define some categories of dimension IDs for fields
  // --------------------------------------------------
  std::map<int, std::vector<int>> fieldDims;
  fieldDims[nLev] = {timId, levId, latId, lonId};  // Fields at levels
  fieldDims[nEdg] = {timId, edgId, latId, lonId};  // Fields at edges
  fieldDims[4] = {timId, forId, latId, lonId};     // Fields at four levels
  fieldDims[1] = {timId, latId, lonId};            // Fields at surface
  fieldDims[0] = {timId, latId, lonId};            // Fields at surface

  // Set float precision for fields
  // ------------------------------
  const int floatPrecision = params_.floatPrecision.value();
  const int ncPrec = (floatPrecision == 4) ? NC_FLOAT : NC_DOUBLE;

  // Define all the fields that will be written
  // ------------------------------------------
  for (auto& field : fields) {
    // Get number of levels for this field
    const int nLevField = field.shape(1);

    // Get dimensions for this field from map
    auto it = fieldDims.find(nLevField);
    if (it == fieldDims.end()) {
      std::ostringstream oss;
      oss << "IOStructuredGrid::writeStructuredFields: "
          << "No entry in fieldDims for field '" << field.name()
          << "' with " << nLevField << " levels.";
      ABORT(oss.str());
    }
    const auto &dims = it->second;

    // Look for fieldname in the iofile configuration and use the value if key found
    const std::string fieldLong = field.name();
    const char * fieldLongC = fieldLong.c_str();
    
    std::string fieldName = fieldLong;
    if (ioNames.has(fieldName)) {
      fieldName = ioNames.getString(fieldLong);
    }

    // Define the field in the file
    nc_rc(nc_def_var(fileId, fieldName.c_str(), ncPrec, dims.size(), dims.data(), &fIv),
          "nc_def_var " + fieldName);

    // Fallback defaults if metadata keys are missing
    std::string unitsStr = "unknown";
    std::string longNameStr = fieldLong;

    // Extract values dynamically if Atlas has them populated
    if (field.metadata().has("units")) {
      unitsStr = field.metadata().get<std::string>("units");
    }
    if (field.metadata().has("long_name")) {
      longNameStr = field.metadata().get<std::string>("long_name");
    }

    const char * units = unitsStr.c_str();
    const char * longNameC = longNameStr.c_str();

    // Write to NetCDF
    nc_rc(nc_put_att_text(fileId, fIv, "units", strlen(units), units),
          "nc_put_att_text " + fieldName + " units");
    nc_rc(nc_put_att_text(fileId, fIv, "long_name", strlen(longNameC), longNameC),
          "nc_put_att_text " + fieldName + " long_name");

    // Insert field into the fieldIvs map
    fieldIvs[field.name()] = fIv;
  }

 /*  // Write ak/bk to the file as global attributes
  // --------------------------------------------
  nc_rc(nc_put_att_double(fileId, NC_GLOBAL, "ak", NC_DOUBLE, ak.size(), ak.data()),
          "nc_put_att_double (ak)");
  nc_rc(nc_put_att_double(fileId, NC_GLOBAL, "bk", NC_DOUBLE, bk.size(), bk.data()),
          "nc_put_att_double (bk)"); */
          
  // Write grid type and dimensions to the file as global attributes
  nc_rc(nc_put_att_text(fileId, NC_GLOBAL, "grid", strlen(gridStr_.c_str()), gridStr_.c_str()),
          "nc_put_att_text (grid)");
  nc_rc(nc_put_att_int(fileId, NC_GLOBAL, "im", NC_INT, 1, &nLon),
          "nc_put_att_int (im)");
  nc_rc(nc_put_att_int(fileId, NC_GLOBAL, "jm", NC_INT, 1, &nLat),
          "nc_put_att_int (im)");

  // End definition mode
  // -------------------
  nc_rc(nc_enddef(fileId), "nc_enddef");

  // Write coordinate data into the file
  // -----------------------------------
  nc_rc(nc_put_var_double(fileId, fieldIvs[params_.latName.value()], latArr.data()),
        "nc_put_var_double (lat)");
  nc_rc(nc_put_var_double(fileId, fieldIvs[params_.lonName.value()], lonArr.data()),
        "nc_put_var_double (lon)");
  nc_rc(nc_put_var_int(fileId, fieldIvs[params_.levName.value()], levArr.data()),
        "nc_put_var_int (lev)");
  nc_rc(nc_put_var_int(fileId, fieldIvs[params_.edgName.value()], edgArr.data()),
        "nc_put_var_int (edg)");
  nc_rc(nc_put_var_int(fileId, fieldIvs[params_.timName.value()], timArr.data()),
        "nc_put_var_int (tim)");

  // Write the fields into the file
  // ------------------------------
  for (auto& field : fields) {
    // Get number of levels for this field
    const int nLevField = field.shape(1);

    // Create a rank 2 view of the field
    const auto fieldView = atlas::array::make_view<double, 2>(field);

    // Vector to hold the packed field
    std::vector<double> values(nLat*nLon*nLevField);

    // Loop over dimensions and pack the field
    for (size_t k = 0; k < nLevField; ++k) {
      for (size_t j = 0; j < nLat; ++j) {
        for (size_t i = 0; i < nLon; ++i) {
          values[k*nLat*nLon + j*nLon + i] = fieldView((nLat - 1 - j) * nLon + i, k);
        }
      }
    }

    // Write the field to the file
    nc_rc(nc_put_var_double(fileId, fieldIvs[field.name()], values.data()),
          "nc_put_var_double " + field.name());
  }

  // Close netCDF file
  // -----------------
  nc_rc(nc_close(fileId), "nc_close");
}

// -------------------------------------------------------------------------------------------------

void IOStructuredGrid::readStructuredFields(const std::string pathFile,
                                            atlas::FieldSet & fields,
                                            const util::DateTime & time,
                                            const eckit::LocalConfiguration & ioNames,
                                            const eckit::LocalConfiguration & ioScaling) const {
  // NetCDF IDs
  int fileId;

  // Open a file to read fields from
  // -------------------------------
  oops::Log::info() << "Reading file " << pathFile << std::endl;
  nc_rc(nc_open(pathFile.c_str(), NC_NOWRITE, &fileId), "nc_open " + pathFile);

  // Get file number of dimensions + their IDs
  // -----------------------------------------
  int ndims;
  nc_rc(nc_inq_ndims(fileId, &ndims), "nc_inq_ndims");

  std::vector<int> dimids(ndims);
  nc_rc(nc_inq_dimids(fileId, &ndims, dimids.data(), 0), "nc_inq_dimids");

  // Create regular grid for determining lat/lon values
  // --------------------------------------------------
  auto writeFunctionSpace_ = atlas::functionspace::StructuredColumns(geom_.functionSpace());
  const atlas::RegularGrid regGrid(writeFunctionSpace_.grid());
  // Get grid dimensions
  // -------------------
  const int nLat = regGrid.ny();
  const int nLon = regGrid.nx();
  const int nLev = geom_.numLevels();    //writeFunctionSpace_.levels(); 
  const int nEdg = nLev + 1;
  const int nFor = 4;
  const int nTim = 1;

  // Ensure that the lat and lon dimensions are found and have the correct lengths
  size_t dimSize;
  bool hasLat = false;
  bool hasLon = false;
  bool hasLev = false;
  bool hasEdg = false;
  bool hasTim = false;
  int latId;
  int lonId;
  int levId;
  int edgId;
  int timId;
  for (int i = 0; i < ndims; ++i) {
    // Get the name and size of the dimension
    char dimName[NC_MAX_NAME + 1];
    nc_rc(nc_inq_dim(fileId, dimids[i], dimName, &dimSize), "nc_inq_dim");

    if (std::string(dimName) == params_.latName.value().c_str()) {
      oops::Log::info() << params_.latName.value().c_str() << " " << nLat << std::endl;
      hasLat = true;
      latId = dimids[i];
      ASSERT(dimSize == nLat);
    } else if (std::string(dimName) == params_.lonName.value().c_str()) {
      oops::Log::info() << params_.lonName.value().c_str() << " " << nLon << std::endl;
      hasLon = true;
      lonId = dimids[i];
      ASSERT(dimSize == nLon);
    } else if (std::string(dimName) == params_.levName.value().c_str()) {
      oops::Log::info() << params_.levName.value().c_str() << " " << nLev << std::endl;
      hasLev = true;
      levId = dimids[i];
      ASSERT(dimSize == nLev);
      // FIX: reference params_ instead of geom_
      //ASSERT(dimSize == static_cast<size_t>(params_.nlevels.value()));
    } else if (std::string(dimName) == params_.edgName.value().c_str()) {
      oops::Log::info() << params_.edgName.value().c_str() << " " << nEdg << std::endl;
      hasEdg = true;
      edgId = dimids[i];
      ASSERT(dimSize == nEdg);
    } else if (std::string(dimName) == params_.timName.value().c_str()) {
      oops::Log::info() << params_.timName.value().c_str() << " " << nTim << std::endl;
      hasTim = true;
      timId = dimids[i];
      ASSERT(dimSize == nTim);
    }
  }

  // Ensure required dimensions were found
  ASSERT(hasLat);
  ASSERT(hasLon);
  ASSERT(hasLev);
  ASSERT(hasEdg);
  ASSERT(hasTim);

  // Read the fields from the file
  // -------------------------------
  for (auto & field : fields) {
    // Get IO name for this field
    std::string fieldName = field.name();
    if (ioNames.has(fieldName)) {
      fieldName = ioNames.getString(field.name());
    }
    oops::Log::info() << "Field " << fieldName << std::endl;
    // Get the variable ID for this field
    int varId;
    int status = nc_inq_varid(fileId, fieldName.c_str(), &varId);
    if (status == NC_ENOTVAR) {
      // Variable is not in this file; skip it silently (or log an info message)
      oops::Log::info() << "Field " << fieldName << " not found in this file, skipping..." << std::endl;
      continue; 
    } else {
      // Check for any other unexpected NetCDF errors
      nc_rc(status, "nc_inq_varid " + fieldName);
    }
    // nc_rc(nc_inq_varid(fileId, fieldName.c_str(), &varId), "nc_inq_varid " + fieldName);


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

    // Ensure that the field has either 3 or 4 dimensions
    ASSERT(ndims == 3 || ndims == 4);

    // Ensure that the dimensions are in the expected order
    size_t nLevField = 0;
    if ( ndims == 3 ) {
      ASSERT(dimids[0] == timId &&
             dimids[1] == latId &&
             dimids[2] == lonId);
      nLevField = 1;
    } else if ( ndims == 4 ) {
      ASSERT((dimids[0] == timId &&
              dimids[1] == levId &&
              dimids[2] == latId &&
              dimids[3] == lonId) ||
             (dimids[0] == timId &&
              dimids[1] == edgId &&
              dimids[2] == latId &&
              dimids[3] == lonId));

      nLevField = field.shape(1);
      if ( dimids[1] == edgId ) {
        ASSERT(nLevField == nEdg);
      } else {
        ASSERT(nLevField == nLev);
      }
    }

    // Read the variable data
    std::vector<double> values(field.size());

    nc_rc(nc_get_var_double(fileId, varId, values.data()), "nc_get_var_double " + fieldName);
/*
    // Create field and unpack data into it
    auto fieldView = atlas::array::make_view<double, 2>(field);

    // Corrected unpacking loop in readStructuredFields:
    for (size_t k = 0; k < nLevField; ++k) {
      for (size_t j = 0; j < nLat; ++j) {
        for (size_t i = 0; i < nLon; ++i) {
          // REMOVED "(nLat - 1 - j)" to match the constructor layout
          //fieldView((nLat - 1 - j) * nLon + i, k) = values[ k*nLat*nLon + j*nLon + i ];
          fieldView(j * nLon + i, k) = values[ k*nLat*nLon + j*nLon + i ];
        }
      }
    }
*/
    // Create field and unpack data into it
    if (field.rank() == 2) {
      // Standard multi-level or Rank-2 surface field [Points, Levels]
      auto fieldView = atlas::array::make_view<double, 2>(field);

      for (size_t k = 0; k < nLevField; ++k) {
        for (size_t j = 0; j < nLat; ++j) {
          for (size_t i = 0; i < nLon; ++i) {
            fieldView(j * nLon + i, k) = values[ k*nLat*nLon + j*nLon + i ];
          }
        }
      }
    } else if (field.rank() == 1) {
      // Pure Rank-1 surface field [Points]
      auto fieldView = atlas::array::make_view<double, 1>(field);
      ASSERT(nLevField == 1); // Double check it's a surface field

      for (size_t j = 0; j < nLat; ++j) {
        for (size_t i = 0; i < nLon; ++i) {
          fieldView(j * nLon + i) = values[ j*nLon + i ];
        }
      }
    } else {
      ABORT("IOStructuredGrid::readStructuredFields - Unsupported field rank: " + std::to_string(field.rank()));
    }

    oops::Log::info() << "Done reading Field " << fieldName << std::endl;
  }
  // Close file
  nc_rc(nc_close(fileId), "nc_close");
}

// -------------------------------------------------------------------------------------------------

void IOStructuredGrid::print(std::ostream & os) const {
  os << classname() << " IO using Atlas Structured Grid";
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jed
