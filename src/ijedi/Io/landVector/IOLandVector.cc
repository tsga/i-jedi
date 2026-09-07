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

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
//#include "oops/base/GeometryData.h"

#include "ijedi/FieldMetadata/FieldsMetadata.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/State/State.h"

//#include "fv3jedi/Utilities/fv3jedi_vertical_remap.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/landVector/IOLandVector.h"



namespace ijedi {
// -------------------------------------------------------------------------------------------------
static IoMaker<IOLandVector> makerIOLandVector_("landVector");
// -------------------------------------------------------------------------------------------------
static inline void nc_rc(const int return_code, const std::string & operation) {
  if (return_code != NC_NOERR) {
    ABORT("IOLandVector netCDF operation \'" + operation + "\' failed with error: "
          + nc_strerror(return_code));
  }
}

IOLandVector::IOLandVector(const Geometry & geom, const Parameters_ & params) 
  : IoBase(geom, params.toConfiguration()), geom_(geom), params_(params) {
    
    // 1. All ranks participate natively in the parallel layout built during initialization
    // 2. No `writeFunctionSpace_` or `interpolator_` members needed!
    oops::Log::trace() << classname() << " Constructor configured for land vector." << std::endl;
}
// -------------------------------------------------------------------------------------------------
IOLandVector::~IOLandVector() {
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}

/*void IOLandVector::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                        const eckit::LocalConfiguration &fileioscaling) const 
{
  util::Timer timer(classname(), "read");
  oops::Log::trace() << classname() << " read started" << std::endl;

  auto readFunctionSpace = atlas::functionspace::PointCloud(geom_.functionSpace());
  const size_t localSize = readFunctionSpace.size();
  // Calculate global count of grid points
  size_t globalSize = 0;
  geom_.getComm().allReduce(localSize, globalSize, eckit::mpi::Operation::SUM);

  // Create serial FieldSet on Rank 0
  atlas::FieldSet fieldsSerial;
  // Construct point cloud containing all global points on Rank 0
  std::vector<atlas::PointXY> globalPoints(globalSize); 
  atlas::functionspace::PointCloud serialFunctionSpace(globalPoints);

  for (const auto & field : x) {
    atlas::Field fSerial = serialFunctionSpace.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(field.levels()));
    fieldsSerial.add(fSerial);
  }

  // Read NetCDF data into fieldsSerial on Rank 0
  size_t layer_index = 0;
  const auto & filenamesOpt = params_.filenames.value();
  if (filenamesOpt != boost::none) {
      if (geom_.getComm().rank() == 0) {
        for (const auto & filename : filenamesOpt.value()) {
          util::DateTime dummyTime; 
          this->readVectorFields(params_.datapath.value() + "/" + filename, 
                                  fieldsSerial, dummyTime, globalSize, layer_index, 
                                  fileionames, fileioscaling);
        }
    }
  } else {
    ABORT("IOLandVector::read: 'filenames' parameter is missing.");
  }
  
  // 2. Extract global indices stored in Geometry to scatter values to local fields
  //auto pc = atlas::functionspace::PointCloud(geom_.functionSpace());  
  //const auto & gidxField = pc.field("global_index");  
  const auto & gidxField = geom_.fields()["global_index"];
  auto gidxView = atlas::array::make_view<atlas::gidx_t, 1>(gidxField);
  // 3. Scatter fields from Rank 0 to all ranks for each field in x
  for (auto & field : x) {
    const std::string name = field.name();
    const size_t levels = field.levels();

    if (field.rank() == 2) {
      // 2D Field: (nodes, levels)
      std::vector<double> globalBuffer;
      if (geom_.getComm().rank() == 0) {
        auto serialView = atlas::array::make_view<double, 2>(fieldsSerial[name]);
        globalBuffer.resize(globalSize * levels);
        for (size_t g = 0; g < globalSize; ++g) {
          for (size_t lev = 0; lev < levels; ++lev) {
            globalBuffer[g * levels + lev] = serialView(g, lev);
          }
        }
      }

      // Broadcast global data from Rank 0 to all processes
      size_t totalDoubles = globalSize * levels;
      geom_.getComm().broadcast(globalBuffer.data(), totalDoubles, 0);

      // Extract local points using 1-based global index mapping
      auto localView = atlas::array::make_view<double, 2>(field);
      for (size_t i = 0; i < localSize; ++i) {
        atlas::gidx_t gidx = gidxView(i) - 1; // Convert 1-based to 0-based
        for (size_t lev = 0; lev < levels; ++lev) {
          localView(i, lev) = globalBuffer[gidx * levels + lev];
        }
      }
    } 
    else if (field.rank() == 1) {
      // 1D Field: (nodes)
      std::vector<double> globalBuffer(globalSize);
      if (geom_.getComm().rank() == 0) {
        auto serialView = atlas::array::make_view<double, 1>(fieldsSerial[name]);
        for (size_t g = 0; g < globalSize; ++g) {
          globalBuffer[g] = serialView(g);
        }
      }

      geom_.getComm().broadcast(globalBuffer.data(), globalSize, 0);

      auto localView = atlas::array::make_view<double, 1>(field);
      for (size_t i = 0; i < localSize; ++i) {
        atlas::gidx_t gidx = gidxView(i) - 1;
        localView(i) = globalBuffer[gidx];
      }
    }
  }

  oops::Log::info() << classname() << " read done" << std::endl;
}*/

void IOLandVector::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                            const eckit::LocalConfiguration &fileioscaling) const 
{
  util::Timer timer(classname(), "read");
  oops::Log::trace() << classname() << " read started" << std::endl;

  // 1. Get your parallel FunctionSpace handle
  auto readFunctionSpace_ = atlas::functionspace::PointCloud(geom_.functionSpace());
  const size_t localSize = readFunctionSpace_.size();
  std::vector<atlas::PointXY> localPoints(localSize);
  const auto & lonlatField = readFunctionSpace_.lonlat();
  auto lonlatView = atlas::array::make_view<double, 2>(lonlatField);
  for (size_t i = 0; i < localSize; ++i) {
    localPoints[i] = atlas::PointXY(lonlatView(i, 0), lonlatView(i, 1));
  }
  /*for (atlas::idx_t i = 0; i < localSize; ++i) {
    atlas::PointLonLat p = readFunctionSpace_.lonlat(i);
    localPoints[i] = atlas::PointXY(p.lon(), p.lat());
  }*/
  
  // 3. Gather all local point counts to determine rank offsets
  std::vector<size_t> localSizes(geom_.getComm().size(), 0);
  geom_.getComm().allGather(localSize, localSizes.begin(), localSizes.end());

  // Convert points to contiguous doubles (2 values per point: lon, lat) for MPI Gatherv
  std::vector<double> localCoords;
  localCoords.reserve(localSize * 2);
  for (const auto &p : localPoints) {
    localCoords.push_back(p.x());
    localCoords.push_back(p.y());
  }

  std::vector<int> counts(geom_.getComm().size());
  std::vector<int> displs(geom_.getComm().size(), 0);
  int totalDoubles = 0;

  for (size_t r = 0; r < geom_.getComm().size(); ++r) {
    counts[r] = static_cast<int>(localSizes[r] * 2);
    if (r > 0) {
      displs[r] = displs[r - 1] + counts[r - 1];
    }
    totalDoubles += counts[r];
  }

  std::vector<double> globalCoords(totalDoubles);
  geom_.getComm().gatherv(localCoords.data(), localCoords.size(),
                          globalCoords.data(), counts.data(), displs.data(), 0);

  // 4. Construct globalPoints vector (populated on rank 0, empty on others)
  std::vector<atlas::PointXY> globalPoints;
  if (geom_.getComm().rank() == 0) {
    globalPoints.reserve(totalDoubles / 2);
    for (size_t i = 0; i < globalCoords.size(); i += 2) {
      globalPoints.emplace_back(globalCoords[i], globalCoords[i + 1]);
    }
  }

  // 5. Create serial FunctionSpace on rank 0
  // On rank 0, globalPoints contains all points; on other ranks, it is empty.
  atlas::functionspace::PointCloud serialFunctionSpace(globalPoints);

  // 3. Allocate fieldsSerial using the serial function space factory
  // This guarantees all ranks have perfectly matching metadata, ranks, and levels
  atlas::FieldSet fieldsSerial;
  for (const auto & field : x) {
    // Pass the matching metadata/options from the target field
    atlas::Field fSerial = serialFunctionSpace.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(field.levels()));
    fieldsSerial.add(fSerial);
  }

  size_t layer_index = 0;

  // 4. Loop through and read filenames on Rank 0
  const auto & filenamesOpt = params_.filenames.value();
  if (filenamesOpt != boost::none) {
    for (const auto & filename : filenamesOpt.value()) {
      if (geom_.getComm().rank() == 0) {
        util::DateTime dummyTime; 
        this->readVectorFields(params_.datapath.value() + "/" + filename, 
                                   fieldsSerial, dummyTime, globalPoints.size(), layer_index, fileionames, fileioscaling);
      }
    }
  } else {
    ABORT("IOLandVector::read: 'filenames' parameter is missing.");
  }

  // 5. Scatter the data globally
  // Since both function spaces share the exact same grid, this will now succeed seamlessly
  readFunctionSpace_.scatter(fieldsSerial, x);
  
  oops::Log::info() << classname() << " read done" << std::endl;
}


// -------------------------------------------------------------------------------------------------

void IOLandVector::readVectorFields(const std::string pathFile,
                                            atlas::FieldSet & fields,
                                            const util::DateTime & time,
                                            size_t num_points,
                                            size_t layer_index,
                                            const eckit::LocalConfiguration & ioNames,
                                            const eckit::LocalConfiguration & ioScaling) const {
  // NetCDF IDs
  int fileId;

  // Open a file to read fields from
  // -------------------------------
  oops::Log::info() << "Reading file " << pathFile << std::endl;
  nc_rc(nc_open(pathFile.c_str(), NC_NOWRITE, &fileId), "nc_open " + pathFile);

  //int num_points = geom_.globalNodeCount();

  // Get file number of dimensions + their IDs
  // -----------------------------------------
  int ndims;
  // soil fields in form: soil_liquid_vol(time, soil_levels, location) ;
  nc_rc(nc_inq_ndims(fileId, &ndims), "nc_inq_ndims");

  std::vector<int> dimids(ndims);
  nc_rc(nc_inq_dimids(fileId, &ndims, dimids.data(), 0), "nc_inq_dimids");

  
  // Ensure that the lat and lon dimensions are found and have the correct lengths
  size_t dimSize;
  size_t numLayers = 1;
  bool hasLoc = false;
  bool hasLayer = false;
  bool hasTim = false;
  int locId;
  int layerId;
  int timId;
  for (int i = 0; i < ndims; ++i) {
    // Get the name and size of the dimension
    char dimName[NC_MAX_NAME + 1];
    nc_rc(nc_inq_dim(fileId, dimids[i], dimName, &dimSize), "nc_inq_dim");

    if (std::string(dimName) == params_.locationName.value().c_str()) {
      hasLoc = true;
      locId = dimids[i];
      ASSERT(dimSize == num_points);
    } else if (std::string(dimName) == params_.layerName.value().c_str()) {
      hasLayer = true;
      layerId = dimids[i];
      ASSERT(dimSize > layer_index);  //TODO: Check if this is the correct assertion for layer_index
    } else if (std::string(dimName) == params_.timeName.value().c_str()) {
      hasTim = true;
      timId = dimids[i];
      ASSERT(dimSize == 1); // Only one time step is expected for this read
    }
  }

  // Ensure required dimensions were found
  ASSERT(hasLoc);
  ASSERT(hasTim);
  //ASSERT(hasLayer || layer_index == -1); // no layer dimension
  
  // Read the fields from the file
  // -------------------------------
  for (auto & field : fields) {
    // Get IO name for this field
    std::string fieldName = field.name();
    if (ioNames.has(fieldName)) {
      fieldName = ioNames.getString(field.name());
    }
    oops::Log::info() << "Field " << fieldName << std::endl;
    oops::Log::info() << "Field Name: " << field.name() << "\n"
                  << "  Rank:  " << field.rank() << "\n"
                  << "  Size:  " << field.size() << "\n"
                  << "  Shape: [";
    for (atlas::idx_t i = 0; i < field.rank(); ++i) {
        oops::Log::info() << field.shape(i) << (i + 1 < field.rank() ? ", " : "");
    }
    oops::Log::info() << "]" << std::endl;
    /*
    Field Name: sheleg
      0:   Rank:  2
      0:   Size:  18320
      0:   Shape: [18320, 1]
      0: Field Name: air_temperature_at_2m
      0:   Rank:  2
      0:   Size:  18320
      0:   Shape: [18320, 1]
      0: Field Name: specfic_humidity_at_2m
      0:   Rank:  2
      0:   Size:  18320
      0:   Shape: [18320, 1]
      0: Field Name: stc
      0:   Rank:  2
      0:   Size:  73280
      0:   Shape: [18320, 4]
    */

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
    int vardimids[NC_MAX_VAR_DIMS];
    nc_rc(nc_inq_var(fileId, varId,
                     nullptr,   // var name (unused)
                     nullptr,   // type (unused)
                     &ndims,
                     vardimids,
                     nullptr),  // attributes (unused)
          "nc_inq_var");

    ASSERT(ndims == 2 || ndims == 3);

     std::vector<double> values(field.size());
    // Ensure that the dimensions are in the expected order
    if ( ndims == 2 ) {
      ASSERT(vardimids[0] == timId && vardimids[1] == locId);
      //size_t start[2] = {0, 0};
      //size_t count[2] = {1, num_points};
      //nc_rc(nc_get_vara_double(fileId, varId, start, count, values.data()), "nc_get_var_double " + fieldName);
    } else if ( ndims == 3 ) {
      ASSERT(vardimids[0] == timId && vardimids[1] == layerId && vardimids[2] == locId && hasLayer && layer_index >= 0);
      //size_t start[3] = {0, 0, 0};  //layer_index, 0};
      //size_t count[3] = {1, 1, num_points};
      //nc_rc(nc_get_vara_double(fileId, varId, start, count, values.data()), "nc_get_var_double " + fieldName);
    }

    nc_rc(nc_get_var_double(fileId, varId, values.data()), "nc_get_var_double " + fieldName);
    
    numLayers = field.shape(1);
    // Create field and unpack data into it
    if (field.rank() == 2) {
      // Standard multi-level or Rank-2 surface field [Points, layers]
      auto fieldView = atlas::array::make_view<double, 2>(field);

      for (size_t k = 0; k < numLayers; ++k) {
          for (size_t i = 0; i < num_points; ++i) {
            fieldView(i, k) = values[ k*num_points + i ];
          }
        }
    } else if (field.rank() == 1) {
      // Pure Rank-1 surface field [Points]
      auto fieldView = atlas::array::make_view<double, 1>(field);
      //ASSERT(numLayers == 1);
      for (size_t j = 0; j < num_points; ++j) {
          fieldView(j) = values[j];
      }
    } else {
      ABORT("IOLandVector::readVectorFields - Unsupported field rank: " + std::to_string(field.rank()));
    }
    oops::Log::info() << "Done reading Field " << fieldName << std::endl;
  }
  // Close file
  nc_rc(nc_close(fileId), "nc_close");
}

void IOLandVector::write(const atlas::FieldSet & fieldsVector,
                             const eckit::LocalConfiguration & fileionames,
                             const eckit::LocalConfiguration & fileioscaling) const 
{

  util::Timer timer(classname(), "write");
  oops::Log::trace() << classname() << " write started" << std::endl;
  oops::Log::info() << classname() << " write started" << std::endl;
  // 1. Get your parallel FunctionSpace handle
  auto readFunctionSpace_ = atlas::functionspace::PointCloud(geom_.functionSpace());
  const size_t localSize = readFunctionSpace_.size();
  std::vector<atlas::PointXY> localPoints(localSize);
  const auto & lonlatField = readFunctionSpace_.lonlat();
  auto lonlatView = atlas::array::make_view<double, 2>(lonlatField);
  for (size_t i = 0; i < localSize; ++i) {
    localPoints[i] = atlas::PointXY(lonlatView(i, 0), lonlatView(i, 1));
  }
  
   // 3. Gather all local point counts to determine rank offsets
  std::vector<size_t> localSizes(geom_.getComm().size(), 0);
  geom_.getComm().allGather(localSize, localSizes.begin(), localSizes.end());
  oops::Log::info() << classname() << " 1" << std::endl;
  // Convert points to contiguous doubles (2 values per point: lon, lat) for MPI Gatherv
  std::vector<double> localCoords;
  localCoords.reserve(localSize * 2);
  for (const auto &p : localPoints) {
    localCoords.push_back(p.x());
    localCoords.push_back(p.y());
  }

  std::vector<int> counts(geom_.getComm().size());
  std::vector<int> displs(geom_.getComm().size(), 0);
  int totalDoubles = 0;

  for (size_t r = 0; r < geom_.getComm().size(); ++r) {
    counts[r] = static_cast<int>(localSizes[r] * 2);
    if (r > 0) {
      displs[r] = displs[r - 1] + counts[r - 1];
    }
    totalDoubles += counts[r];
  }

  std::vector<double> globalCoords(totalDoubles);
  geom_.getComm().gatherv(localCoords.data(), localCoords.size(),
                          globalCoords.data(), counts.data(), displs.data(), 0);
  oops::Log::info() << classname() << " 2" << std::endl;
  // 4. Construct globalPoints vector (populated on rank 0, empty on others)
  std::vector<atlas::PointXY> globalPoints;
  if (geom_.getComm().rank() == 0) {
    globalPoints.reserve(totalDoubles / 2);
    for (size_t i = 0; i < globalCoords.size(); i += 2) {
      globalPoints.emplace_back(globalCoords[i], globalCoords[i + 1]);
    }
  }

  // 5. Create serial FunctionSpace on rank 0
  // On rank 0, globalPoints contains all points; on other ranks, it is empty.
  atlas::functionspace::PointCloud serialFunctionSpace(globalPoints);
  oops::Log::info() << classname() << " 3" << std::endl;
  // 3. Allocate fieldsSerial using the serial function space factory
  // This guarantees all ranks have perfectly matching metadata, ranks, and levels
  atlas::FieldSet fieldsSerial;
  for (const auto & field : fieldsVector) {
    // Pass the matching metadata/options from the target field
    atlas::Field fSerial = serialFunctionSpace.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(field.levels()));
    fieldsSerial.add(fSerial);
  }

  /* Allocate fieldsSerial using the serial function space factory
  atlas::FieldSet fieldsSerial;
  if (geom_.getComm().rank() == 0) {
    for (const auto & field : fieldsVector) {
      // Pass the matching metadata/options from the target field
      atlas::Field fSerial = serialFunctionSpace.createField<double>(
          atlas::option::name(field.name()) | atlas::option::levels(field.levels()));
      fieldsSerial.add(fSerial);
    }
  }  */

  // Gather distributed data smoothly from all parallel ranks back into Rank 0
  readFunctionSpace_.gather(fieldsVector, fieldsSerial);
  oops::Log::info() << classname() << " 4" << std::endl;
  // Resolve the valid time dynamically
  util::DateTime validTime;
  if (fieldsVector.metadata().has("time")) {
    validTime = util::DateTime(fieldsVector.metadata().get<std::string>("time"));
  } else if (fieldsVector.metadata().has("timestamp")) {
    validTime = util::DateTime(fieldsVector.metadata().get<std::string>("timestamp"));
  } else if (params_.dateTime.value() != boost::none) {
    validTime = util::DateTime(params_.dateTime.value().value());
  } else {
    oops::Log::warning() << "Using default DateTime placeholder for file writing." << std::endl;
    validTime = util::DateTime("2026-07-01T12:00:00Z");
  }

  // Write to disk exclusively on rank 0
  if (geom_.getComm().rank() == 0) {
    //const util::DateTime dateTime(datTimeString);
    this->writeVectorFields(fieldsSerial, globalPoints.size(), fileionames, fileioscaling);
  }
  oops::Log::trace() << classname() << " write done" << std::endl;

}

void IOLandVector::writeVectorFields(const atlas::FieldSet & fields,
                                             //const util::DateTime & time,
                                             size_t num_locations,
                                             const eckit::LocalConfiguration & ioNames,
                                             const eckit::LocalConfiguration & ioScaling) const {
  
  // NetCDF IDs
  // ----------
  int fileId, fIv, locId, layerId, timId;
  std::map<std::string, int> fieldIvs;
  int nTim = 1;
  int num_layers = 4;

  //int num_locations = geom_.globalNodeCount();

  // Get the name of the file and adjust with datetime
  // -------------------------------------------------
  /*std::string pathFile = params_.filename.value();
   // TODO: remove this later--for now to deal with default name
  if (pathFile.find("%Y") == std::string::npos) {
    pathFile += "%Y%m%d_%H%M%Sz";
  }
  if (pathFile.find(".nc") == std::string::npos) {
    pathFile += ".nc4";
  }

  // Format the datetime string
  pathFile = time.formatString(pathFile);*/

  std::string pathFile = params_.datapath.value() + "/" + params_.filename.value();

  // Replace member number (ensemble applciaitons)
  util::stringfunctions::swapNameMember(params_.toConfiguration(), pathFile);

  // Create a file to write fields into
  // ----------------------------------
  nc_rc(nc_create(pathFile.c_str(), NC_CLOBBER | NC_NETCDF4, &fileId), "nc_create" + pathFile);
  oops::Log::warning() << "nc created" << std::endl;

  // Set float precision for fields
  // ------------------------------
  const int floatPrecision = params_.floatPrecision.value();
  const int ncPrec = (floatPrecision == 4) ? NC_FLOAT : NC_DOUBLE;

  const auto & dateTimeOpt = params_.dateTime.value();
  util::DateTime time;
  if (dateTimeOpt != boost::none) {
    time = util::DateTime(dateTimeOpt.value());
  } else {
    ABORT("invalid datetime for write");
  }

  //Get time in seconds since epoch
  const util::DateTime epoch("1970-01-01T00:00:00Z");
  const util::Duration duration = time - epoch;
  int seconds_since_epoch = duration.toSeconds();

  nc_rc(nc_def_dim(fileId, params_.locationName.value().c_str(), num_locations, &locId), "nc_def_dim (location)");
  nc_rc(nc_def_dim(fileId, params_.layerName.value().c_str(), num_layers, &layerId), "nc_def_dim (layer)");
  nc_rc(nc_def_dim(fileId, params_.timeName.value().c_str(), nTim, &timId), "nc_def_dim (time)");

  // Write the dimension variables: only time for now
  nc_rc(nc_def_var(fileId, params_.timeName.value().c_str(), NC_INT, 1, &timId, &fIv),
        "nc_def_var (tim)");
  nc_rc(nc_put_att_text(fileId, fIv, "long name", strlen("time"), "time"),
        "nc_put_att_text (long name time)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("seconds since 1970-01-01 00:00:00"), "seconds since 1970-01-01 00:00:00"),
        "nc_put_att_text (units time)");
  fieldIvs[params_.timeName.value()] = fIv;

  // Define some categories of dimension IDs for fields
  // --------------------------------------------------
  
  // Define all the fields that will be written
  // ------------------------------------------
  oops::Log::info() << "In num locations " << num_locations << std::endl;
  for (auto& field : fields) {
    
    // Get IO name for this field
    std::string fieldNameI = field.name();
    if (ioNames.has(fieldNameI)) {
      fieldNameI = ioNames.getString(field.name());
    }
    oops::Log::info() << "Field " << fieldNameI << std::endl;
    oops::Log::info() << "Field Name: " << field.name() << "\n"
                  << "  Rank:  " << field.rank() << "\n"
                  << "  Size:  " << field.size() << "\n"
                  << "  Shape: [";
    for (atlas::idx_t i = 0; i < field.rank(); ++i) {
        oops::Log::info() << field.shape(i) << (i + 1 < field.rank() ? ", " : "");
    }
    oops::Log::info() << "]" << std::endl;

    ASSERT(field.shape(0) == num_locations);  // Ensure the field has the expected number of locations

    // Get dimensions for this field 
    //const auto &dims = field.shape();

    std::vector<int> fieldDims; // = {timId, locId};  // only one layer written out
    // Create field and unpack data into it
    if (field.rank() == 2) {
      fieldDims = {timId, layerId, locId};
    } else if (field.rank() == 1) {
      fieldDims = {timId, locId};
    } else {
      ABORT("IOLandVector::readVectorFields - Unsupported field rank: " + std::to_string(field.rank()));
    }

    // Look for fieldname in the iofile configuration and use the value if key found
    const std::string fieldLong = field.name();
    const char * fieldLongC = fieldLong.c_str();
    
    std::string fieldName = fieldLong;
    if (ioNames.has(fieldName)) {
      fieldName = ioNames.getString(fieldLong);
    }

    // Define the field in the file
    nc_rc(nc_def_var(fileId, fieldName.c_str(), ncPrec, fieldDims.size(), fieldDims.data(), &fIv), "nc_def_var " + fieldName);

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
    nc_rc(nc_put_att_text(fileId, fIv, "units", strlen(units), units), "nc_put_att_text " + fieldName + " units");
    nc_rc(nc_put_att_text(fileId, fIv, "long_name", strlen(longNameC), longNameC), 
          "nc_put_att_text " + fieldName + " long_name");

    // Insert field into the fieldIvs map
    fieldIvs[field.name()] = fIv;
  }
  // End definition mode
  // -------------------
  nc_rc(nc_enddef(fileId), "nc_enddef");

  // Write coordinate data 
  nc_rc(nc_put_var_int(fileId, fieldIvs[params_.timeName.value()], &seconds_since_epoch), "nc_put_var_int (time)");

  // Write the fields into the file
  // ------------------------------
  for (auto& field : fields) {

    oops::Log::info() << "Writing Field " << field.name() << std::endl;
 
    // Create field and unpack data into it
    if (field.rank() == 2) {
      // int numLayers = field.shape(1);
      // Standard multi-level or Rank-2 surface field [Points, layers]
      auto fieldView = atlas::array::make_view<double, 2>(field);
      // Vector to hold the packed field
      std::vector<double> values(field.size());
      for (size_t k = 0; k < field.shape(1); ++k) {
          for (size_t i = 0; i < field.shape(0); ++i) {
            values[ k*num_locations + i ] = fieldView(i, k);
          }
      }
      nc_rc(nc_put_var_double(fileId, fieldIvs[field.name()], values.data()), "nc_put_var_double " + field.name());
    } else if (field.rank() == 1) {
      // Pure Rank-1 surface field [Points]
      auto fieldView = atlas::array::make_view<double, 1>(field);
      // Vector to hold the packed field
      std::vector<double> values(field.size());
      for (size_t j = 0; j < field.shape(0); ++j) {
          values[j] = fieldView(j);
      }
      nc_rc(nc_put_var_double(fileId, fieldIvs[field.name()], values.data()), "nc_put_var_double " + field.name());
    } else {
      ABORT("IOLandVector::readVectorFields - Unsupported field rank: " + std::to_string(field.rank()));
    }
    // Write the field to the file
    //nc_rc(nc_put_var_double(fileId, fieldIvs[field.name()], values.data()), "nc_put_var_double " + field.name());
  }

  // Close netCDF file
  // -----------------
  nc_rc(nc_close(fileId), "nc_close");
  oops::Log::info() << "IOLandVector::writeVectorFields done " << std::endl;
}

// -------------------------------------------------------------------------------------------------

void IOLandVector::print(std::ostream & os) const {
  os << classname() << " IO for land vector using Atlas PointCloud FunctionSpace";
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jed
