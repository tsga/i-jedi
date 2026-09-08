#pragma once

#include <cstring>
#include <vector>
#include <array>

// Forward declarations of C-interoperable Fortran routines
extern "C" {
    void c_add_increment_soil(
        int* myrank,
        int* lsoil,
        int* lsoil_incr,
        int* lensfc,
        int* soilsnow_tile,
        int* upd_stc,
        int* upd_slc,
        int* print_summary,
        int* print_debug
        float* stc_state,
        float* slc_state,
        float* smc_state,
        float* stcinc,
        float* slcinc,
        int* stc_updated,
        int* slc_updated,
    );

    void c_apply_land_da_adjustments_soil(
        int* lsoil_incr,
        int* isot,
        int* ivegsrc,
        int* lensfc,
        int* lsoil,
        int* isoiltype,
        int* mask,
        float* stc_bck,
        float* stc_adj,
        float* smc_adj,
        float* slc_adj,
        int* stc_updated,
        int* slc_updated,
        float* zsoil,
        int* upd_stc,
        int* upd_slc,
        int* myrank,
        int* print_summary,
        int* print_debug
    );

    void c_calculate_landinc_mask(
        float* swe,
        int* vtype,
        int* stype,
        int* lensfc,
        int* veg_type_landice,
        int* mask
    );
}

// C++ wrapper class for convenient calling from C++
class SoilIncrements {
public:
    static void calculateLandIncrementMask(
        const std::vector<float>& swe_state,       // (lensfc, lsoil)
        const std::vector<int>& vtype,  // (lensfc)
        const std::vector<int>& stype,  // (lensfc)
        int lensfc,
        int veg_type_landice,
        std::vector<int>& mask           // (lensfc)
    ) {
        int i_lensfc = lensfc;
        int i_veg_type_landice = veg_type_landice;
        c_calculate_landinc_mask(
            const_cast<float*>(swe_state.data()),
            const_cast<int*>(vtype.data()),
            const_cast<int*>(stype.data()),
            &i_lensfc,
            &i_veg_type_landice,
            mask.data()
        );
    }   

    static void addIncrementSoil(
        int myrank,
        int lsoil,
        int lsoil_incr,
        int lensfc,
        const std::vector<int>& soilsnow_tile,  // (lensfc)
        bool upd_stc,
        bool upd_slc,
        bool print_summary,
        bool print_debug,
        std::vector<float>& stc_state,          // (lensfc, lsoil)
        std::vector<float>& slc_state,          // (lensfc, lsoil)
        std::vector<float>& smc_state,          // (lensfc, lsoil)
        const std::vector<float>& stcinc,       // (lensfc, lsoil)
        const std::vector<float>& slcinc,       // (lensfc, lsoil)
        std::vector<int>& stc_updated,          // (lensfc)
        std::vector<int>& slc_updated,          // (lensfc)
    ) {
        int i_myrank = myrank;
        int i_lsoil = lsoil;
        int i_lsoil_incr = lsoil_incr;
        int i_lensfc = lensfc;
        int i_upd_stc = upd_stc ? 1 : 0;
        int i_upd_slc = upd_slc ? 1 : 0;
        int i_print_summary = print_summary ? 1 : 0;
        int i_print_debug = print_debug ? 1 : 0;

        c_add_increment_soil(
            &i_myrank,
            &i_lsoil,
            &i_lsoil_incr,
            &i_lensfc,
            const_cast<int*>(soilsnow_tile.data()),
            &i_upd_stc,
            &i_upd_slc,
            &i_print_summary,
            &i_print_debug
            stc_state.data(),
            slc_state.data(),
            smc_state.data(),
            const_cast<float*>(stcinc.data()),
            const_cast<float*>(slcinc.data()),
            stc_updated.data(),
            slc_updated.data()
        );
    }

    static void applyLandDAadjustmentsSoil(
        int lsoil_incr,
        int isot,
        int ivegsrc,
        int lensfc,
        int lsoil,
        const std::vector<int>& isoiltype,      // (lensfc)
        const std::vector<int>& mask,           // (lensfc)
        const std::vector<float>& stc_bck,      // (lensfc, lsoil)
        std::vector<float>& stc_adj,            // (lensfc, lsoil)
        std::vector<float>& smc_adj,            // (lensfc, lsoil)
        std::vector<float>& slc_adj,            // (lensfc, lsoil)
        const std::vector<int>& stc_updated,    // (lensfc)
        const std::vector<int>& slc_updated,    // (lensfc)
        const std::array<float, 4>& zsoil,      // (lsoil)
        bool upd_stc,
        bool upd_slc,
        int myrank,
        bool print_summary,
        bool print_debug
    ) {
        int i_lsoil_incr = lsoil_incr;
        int i_isot = isot;
        int i_ivegsrc = ivegsrc;
        int i_lensfc = lensfc;
        int i_lsoil = lsoil;
        int i_myrank = myrank;
        int i_upd_stc = upd_stc ? 1 : 0;
        int i_upd_slc = upd_slc ? 1 : 0;
        int i_print_summary = print_summary ? 1 : 0;
        int i_print_debug = print_debug ? 1 : 0;

        // Create mutable copy of zsoil for Fortran
        float zsoil_copy[4];
        std::copy(zsoil.begin(), zsoil.end(), zsoil_copy);

        c_apply_land_da_adjustments_soil(
            &i_lsoil_incr,
            &i_isot,
            &i_ivegsrc,
            &i_lensfc,
            &i_lsoil,
            const_cast<int*>(isoiltype.data()),
            const_cast<int*>(mask.data()),
            const_cast<float*>(stc_bck.data()),
            stc_adj.data(),
            smc_adj.data(),
            slc_adj.data(),
            const_cast<int*>(stc_updated.data()),
            const_cast<int*>(slc_updated.data()),
            zsoil_copy,
            &i_upd_stc,
            &i_upd_slc,
            &i_myrank,
            &i_print_summary,
            &i_print_debug
        );
    }
};