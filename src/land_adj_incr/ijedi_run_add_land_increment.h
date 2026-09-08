#pragma once
#include <string>
#include <vector>
#include "eckit/config/LocalConfiguration.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

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

      // We assume both state and increment are in the same geometry
      const ijedi::Geometry<ijedi::Traits> geom_(eckit::LocalConfiguration(config_, "geometry"),
                                                      comm_);

      // Read state
      ijedi::State<ijedi::Traits> xx(geom_, eckit::LocalConfiguration(config_, "state"));
      oops::Log::test() << "State: " << xx << std::endl;

      // Read increment
      const eckit::LocalConfiguration incParams(config_, "increment");
      oops::Variables addedVars(incParams, "added variables");
      ijedi::Increment<ijedi::Traits> dx(geom_, addedVars, xx.validTime());
      dx.read(incParams);
      oops::Log::test() << "Increment: " << dx << std::endl;

      // Scale increment
      if (incParams.has("scaling factor")) {
        dx *= incParams.getDouble("scaling factor");
        oops::Log::test() << "Scaled Increment: " << dx << std::endl;
      }

      // Add increment to state
      call add_increment_soil(lsoil_incr,noahmp_state%stc_inc,noahmp_state%slc_inc, &
               noahmp_state%stc,noahmp_state%smc,noahmp_state%slc,&
               noahmp_state%stc_updated,noahmp_state%slc_updated,noahmp_state%soilsnow_tile,noahmp_state%soilsnow_tile,&
               len_land_vec,lsoil,myrank, upd_stc, upd_slc, print_summary, print_debug)

            call apply_land_da_adjustments_soil(lsoil_incr, isot, ivegsrc, len_land_vec, &
                 lsoil, noahmp_state%stype, noahmp_state%soilsnow_tile,noahmp_state%stc_bkg, &
                 noahmp_state%stc,noahmp_state%smc,noahmp_state%slc, &
                 noahmp_state%stc_updated,noahmp_state%slc_updated, zsoil, upd_stc, upd_slc, myrank, print_summary, print_debug)
            

      xx += dx;
      oops::Log::test() << "Updated State: " << xx << std::endl;

      // Write updated state to file
      xx.write(eckit::LocalConfiguration(config_, "output state"));

    }
    

  };
//}  // namespace gdasapp