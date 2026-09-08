#pragma once
#include <string>
#include <vector>
#include "eckit/config/LocalConfiguration.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

#include "soil_increments_cpp_interface.h"

//namespace gdasapp {
  class runAddLandIncrement {
   public:
    runAddLandIncrement(const eckit::Configuration & config, const eckit::mpi::Comm & comm)
      : config_(config), comm_(comm) {}

   private:
    const eckit::Configuration & config_;
    const eckit::mpi::Comm & comm_;

    static constexpr std::array<float, 4> zsoil = ({ -0.1, -0.4, -1.0, -2.0 });
    static constexpr int veg_type_landice = 15;
    static constexpr int lsoil = 4;     // zsoil is hard-coded for 4 layers
    static constexpr int ivegsrc = 1;   // The NOAHMP LSM expects that the ivegsrc physics parameter is 1
    static constexpr int isot = 1;      // Noahmp expects 1
    // hard coded defaults--unlikely to change
    bool frac_grid = true;
    float fice_threshold = 0.0;
    float lfrac_threshold = 0.0001;

   public:
    static const std::string classname() {return "runAddLandIncrement";}

    void run(){

      int myrank = comm_.rank();

      // We assume both state and increment are in the same geometry
      const ijedi::Geometry<ijedi::Traits> geom_(eckit::LocalConfiguration(config_, "geometry"),
                                                      comm_);

      // Read state
      ijedi::State<ijedi::Traits> xx(geom_, eckit::LocalConfiguration(config_, "background state"));
      oops::Log::test() << "State: " << xx << std::endl;

      // Read increment
      const eckit::LocalConfiguration incParams(config_, "increment");
      oops::Variables incVars(incParams, "variables");
      ijedi::Increment<ijedi::Traits> dx(geom_, incVars, xx.validTime());
      dx.read(incParams);
      oops::Log::test() << "Increment: " << dx << std::endl;

      // Scale increment
      if (incParams.has("scaling factor")) {
        dx *= incParams.getDouble("scaling factor");
        oops::Log::test() << "Scaled Increment: " << dx << std::endl;
      }

      atlas::FieldSet bkg_fs;
      xx.toFieldSet(bkg_fs);

      // check required fields exist in the fieldset
      if (!bkg_fs.has("sheleg") ) {
          oops::Log::error() << "Missing required fields SWE in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkg_swe = atlas::array::make_view<double, 2>(bkg_fs["sheleg"]);
      if (!bkg_fs.has("vtype") ) {
          oops::Log::error() << "Missing required fields vtype in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkg_vtype = atlas::array::make_view<float, 2>(bkg_fs["vtype"]);
      if (!bkg_fs.has("stype") ) {
          oops::Log::error() << "Missing required fields stype in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkg_stype = atlas::array::make_view<float, 2>(bkg_fs["stype"]);

      if (!bkg_fs.has("stc") ) {
          oops::Log::error() << "Missing required fields stc in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkg_stc = atlas::array::make_view<double, 2>(bkg_fs["stc"]);
      if (!bkg_fs.has("slc") ) {
          oops::Log::error() << "Missing required fields slc in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkg_slc = atlas::array::make_view<double, 2>(bkg_fs["slc"]);
      if (!bkg_fs.has("smc") ) {
          oops::Log::error() << "Missing required fields smc in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkg_smc = atlas::array::make_view<double, 2>(bkg_fs["smc"]);
         
      bool upd_stc = false, upd_slc = false;
      atlas::FieldSet inc_fs;
      dx.toFieldSet(inc_fs);
      if (bkg_fs.has("stc_inc")) {
        auto stc_inc = atlas::array::make_view<double, 2>(bkg_fs["stc_inc"]);
        upd_stc = true;
        oops::Log::trace << "Updating stc" << std::endl;
      }
      if (inc_fs.has("slc_inc")) {
        auto slc_inc = atlas::array::make_view<double, 2>(inc_fs["slc_inc"]);
        upd_slc = true;
        oops::Log::trace << "Updating slc" << std::endl;
      }
      
      // read/construct mask for landice and snow tiles 
      int lsoil_incr = 2;
      config_.get("lsoil_incr", lsoil_incr);
      int len_land_vec = bkg_fs["sheleg"].shape(0);
      // std::vector<int> mask_landice(geom_.nlevsfc(), 0);
      std::vector<int> soil_mask(len_land_vec, 0);
      std::vector<int> ivtype(len_land_vec, -1);
      std::vector<int> istype(len_land_vec, -1);
      std::vector<std::vector<float>> bk_bkg_stc(len_land_vec, std::vector<float>(lsoil, 0.0));
      for (int i = 0; i < len_land_vec; ++i) {
        ivtype[i] = static_cast<int>(bkg_vtype(i, 0));
        istype[i] = static_cast<int>(bkg_stype(i, 0));
        for (int j = 0; j < lsoil; ++j) {
            bk_bkg_stc[i][j] = bkg_stc(i, j);
        }  
      }
      
      // TODO: check if landfrac and icefrac are relevant for mask
      SoilIncrements::calculateLandIncrementMask(bkg_swe, ivtype, istype, 
                              len_land_vec, veg_type_landice, soil_mask);
      // zero out increments for mask not equal to 1
      for (int i = 0; i < len_land_vec; ++i) {
          if (soil_mask[i] != 1) {
              for (int j = 0; j < lsoil_incr; ++j) {
                  stc_inc(i, j) = 0.0;
                  slc_inc(i, j) = 0.0;
              }
          }
      }

      // Add increment to state
      // TODO: see if we can sync xx and dx configs, zero-out dx, then xx += dx;
      bool print_summary = false, print_debug = false;
      config_.get("print_summary", print_summary);
      config_.get("print_debug", print_debug);
      std::vector<int> stc_updated(len_land_vec, 0);
      std::vector<int> slc_updated(len_land_vec, 0);
      SoilIncrements::addIncrementSoil(
          myrank, lsoil, lsoil_incr, len_land_vec, soil_mask, 
          upd_stc, upd_slc, print_summary, print_debug,
          bkg_stc, bkg_slc, bkg_smc, stc_inc, slc_inc,       
          stc_updated, slc_updated  
      );
       
      // post-increment adjustments to ensure consistency b/n soil T and soil M
      SoilIncrements::applyLandDAadjustmentsSoil(
          lsoil_incr, isot, ivegsrc, len_land_vec, lsoil,
          istype, soil_mask, 
          bk_bkg_stc, bkg_stc, bkg_smc, bkg_slc,
          stc_updated, slc_updated, zsoil,
          upd_stc, upd_slc, myrank, print_summary, print_debug,
      );

      // updated state
      xx.fromFieldSet(bkg_fs);      
      oops::Log::test() << "Updated State: " << xx << std::endl;

      // Write updated state to file
      xx.write(eckit::LocalConfiguration(config_, "output state"));

    }
    

  };
//}  // namespace gdasapp