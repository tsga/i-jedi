!> @file
!! @brief Routines for applyng soil DA increments
!! copied from land_increments
!! @author Clara Draper ESRL/PSL
!! Tseganeh ZG April 2026 bring soil specific parts close to snow DA incrments code in GDASApp

module soil_increments

    private

    public add_increment_soil
    public calculate_landinc_mask
    public apply_land_da_adjustments_soil

    integer, parameter       :: lsm_noahmp=2      !< flag for NOAHMP land surface model
    real, parameter          :: tfreez=273.16 !< con_t0c  in physcons
    
contains


 !> Read in soil state increments (on the cubed-sphere 
 !! grid),and add to the soil states. Adapted from original add_gsi_increment_soil routine.
 !!
 !! @param[in] SLCINC Liquid soil moisture increments on the cubed-sphere tile
 !! @param[in] STCINC Soil temperature increments on the cubed-sphere tile
 !! @param[inout] STC_STATE Soil temperature state vector
 !! @param[inout] SMC_STATE Soil moisture (liquid plus solid) state vector
 !! @param[inout] SLC_STATE Liquid soil moisture state vector
 !! @param[out] stc_updated Integer to record whether STC in each grid cell was updated 
 !! @param[out] slc_updated Integer to record whether SMC in each grid cell was updated 
 !! @param[in] SOILSNOW_TILE Land mask for increments on the cubed-sphere tile
 !! @param[in] SOILSNOW_FG_TILE First guess land mask for increments on the
 !!             cubed-sphere tile
 !! @param[in] LENSFC Number of land points on a tile
 !! @param[in] LSOIL Number of soil layers
 !! @param[in] LSOIL_INCR Number of soil layers (from top) to apply soil increments to
 !! @param[in] MYRANK MPI rank number
 !!
 !! @author Yuan Xue. 11/2023
 !! updates April 2024, Tseganeh ZG: (calling from jedi-apply_lnd_inc)
 !! stc_upd, slc_upd user inputs
 !! mpi init in calling rourine
 !! applying only for NoahMP

subroutine add_increment_soil(lsoil_incr,stcinc,slcinc,stc_state,smc_state,slc_state,stc_updated,&
              slc_updated,soilsnow_tile,soilsnow_fg_tile,lensfc,lsoil,myrank, upd_stc, upd_slc, print_summary, print_debug)

    use mpi

    implicit none

    integer, intent(in)      :: lsoil_incr, lensfc, lsoil, myrank

    integer, intent(in)      :: soilsnow_tile(lensfc), soilsnow_fg_tile(lensfc)
    real, intent(inout)      :: stc_state(lensfc, lsoil)
    real, intent(inout)      :: slc_state(lensfc, lsoil)
    real, intent(inout)      :: smc_state(lensfc, lsoil)
    integer, intent(out)     :: stc_updated(lensfc), slc_updated(lensfc)
    logical, intent(in)      :: upd_stc, upd_slc, print_summary, print_debug
    
    integer                  :: ij
    integer                  :: mask_tile, mask_fg_tile

    real                     :: stcinc(lensfc,lsoil)
    real                     :: slcinc(lensfc,lsoil)

    integer                  :: k, nother, nsnowupd
    integer                  :: nstcupd, nslcupd,  nfrozen, nfrozen_upd
    logical                  :: soil_freeze, soil_ice

    stc_updated=0
    slc_updated=0

    if (print_summary .and. myrank==0) then
      print*
      print*,'add soil increments on cubed-sphere tiles'
      print*,'updating soil temps', upd_stc
      print*,'updating soil moisture', upd_slc
      print*,'adding to first ', lsoil_incr, ' surface layers only'
    endif
    ! initialize variables for counts statitics to be zeros
    nother = 0 ! grid cells not land
    nsnowupd = 0  ! grid cells with snow (temperature not yet updated)
    nslcupd = 0 ! grid cells that are updated
    nstcupd = 0 ! grid cells that are updated
    nfrozen = 0 ! not update as frozen soil
    nfrozen_upd = 0 ! not update as frozen soil

    ij_loop : do ij = 1, lensfc

        mask_tile    = soilsnow_tile(ij)
        mask_fg_tile = soilsnow_fg_tile(ij)

        !----------------------------------------------------------------------
        ! mask: 1  - soil, 2 - snow, 0 - land-ice, -1 - not land
        !----------------------------------------------------------------------

        if (mask_tile <= 0) then ! skip if neither soil nor snow
         nother = nother + 1
         cycle ij_loop
        endif

        !----------------------------------------------------------------------
        ! if snow is present before or after snow update, skip soil analysis
        !----------------------------------------------------------------------

        if (mask_fg_tile == 2 .or. mask_tile == 2) then
         nsnowupd = nsnowupd + 1
         cycle ij_loop
        endif

        !----------------------------------------------------------------------
        !  do update to soil temperature grid cells
        !----------------------------------------------------------------------

        if (mask_tile == 1) then

           !----------------------------------------------------------------------
           !  add the interpolated increment to the background
           !----------------------------------------------------------------------

           soil_freeze=.false.
           soil_ice=.false.
           do k = 1, lsoil_incr

             if ( stc_state(ij,k) < tfreez)  soil_freeze=.true.
             if ( smc_state(ij,k) - slc_state(ij,k) > 0.001 )  soil_ice=.true.

             if (upd_stc) then
                stc_state(ij,k) = stc_state(ij,k) + stcinc(ij,k)     !TODO: do not add if < min_inc
                if (k==1) then
                    stc_updated(ij) = 1
                    nstcupd = nstcupd + 1
                endif
             endif

             if ( (stc_state(ij,k) < tfreez) .and. (.not. soil_freeze) .and. (k==1) )&
                   nfrozen_upd = nfrozen_upd + 1

             ! do not do updates if this layer or any above is frozen
             if ( (.not. soil_freeze ) .and. (.not. soil_ice ) ) then
                if (upd_slc) then
                if (k==1) then
                    nslcupd = nslcupd + 1
                    slc_updated(ij) = 1
                endif
                   ! apply zero limit here (higher, model-specific limits are
                   ! later)
                   slc_state(ij,k) = max(slc_state(ij,k) + slcinc(ij,k), 0.0)
                   smc_state(ij,k) = max(smc_state(ij,k) + slcinc(ij,k), 0.0)
                endif
             else
                if (k==1) nfrozen = nfrozen+1
             endif

           enddo

        endif ! if soil/snow point

   enddo ij_loop
   
   if (print_summary .and. myrank==0) then 
     write(*,'(a,i2)') ' statistics of grids number processed for rank : ', myrank
     write(*,'(a,i8)') ' soil grid total', lensfc
     write(*,'(a,i8)') ' soil grid cells slc updated = ',nslcupd
     write(*,'(a,i8)') ' soil grid cells stc updated = ',nstcupd
     write(*,'(a,i8)') ' soil grid cells not updated, frozen = ',nfrozen
     write(*,'(a,i8)') ' soil grid cells update, became frozen = ',nfrozen_upd
     write(*,'(a,i8)') ' (not updated yet) snow grid cells = ', nsnowupd
     write(*,'(a,i8)') ' grid cells, without soil or snow = ', nother
   endif

   if (print_debug) then
     write(*,'(a,i2)') ' statistics of grids number processed for rank : ', myrank
     write(*,'(a,i8)') ' soil grid total', lensfc
     write(*,'(a,i8)') ' soil grid cells slc updated = ',nslcupd
     write(*,'(a,i8)') ' soil grid cells stc updated = ',nstcupd
     write(*,'(a,i8)') ' soil grid cells not updated, frozen = ',nfrozen
     write(*,'(a,i8)') ' soil grid cells update, became frozen = ',nfrozen_upd
     write(*,'(a,i8)') ' (not updated yet) snow grid cells = ', nsnowupd
     write(*,'(a,i8)') ' grid cells, without soil or snow = ', nother
   endif


end subroutine add_increment_soil


!> Calculate soil mask for land on model grid.
!! Output is 1  - soil, 2 - snow-covered, 0 - land ice, -1  not land.
!!
!! @param[in] lensfc  Number of land points for this tile 
!! @param[in] veg_type_landice Value of vegetion class that indicates land-ice
!! @param[in] stype Soil type
!! @param[in] swe Model snow water equivalent
!! @param[in] vtype Model vegetation type
!! @param[out] mask Land mask for increments
!! @author Clara Draper @date March 2021
!! @author Yuan Xue: introduce stype to make the mask calculation more generic
subroutine calculate_landinc_mask(swe,vtype,stype,lensfc,veg_type_landice,mask)
 
    implicit none

    integer, intent(in)           :: lensfc, veg_type_landice
    real, intent(in)              :: swe(lensfc)
    integer, intent(in)           :: vtype(lensfc),stype(lensfc)
    integer, intent(out)          :: mask(lensfc)

    integer :: i

    mask = -1 ! not land

    ! land (but not land-ice)
    do i=1,lensfc
        if (stype(i) .GT. 0) then
          if (swe(i) .GT. 0.001) then ! snow covered land
                mask(i) = 2
          else                        ! non-snow covered land
                mask(i) = 1
          endif
        end if ! else should work here too
        if ( vtype(i) ==  veg_type_landice  ) then ! land-ice
                mask(i) = 0
        endif
    end do

end subroutine calculate_landinc_mask

!> Make adjustments to dependent variables after applying land increments.
!! These adjustments are model-dependent, and are currently only coded
!! if full for Noah LSM. 
!! For Noah LSM, copy relevent code blocks from model code (same as has
!! been done in sfc_sub).
!! For Noah-MP, the adjustment scheme shown below as of 11/09/2023:
!! Case 1: frozen ==> frozen, recalculate slc following opt_frz=1, smc remains
!! Case 2: unfrozen ==> frozen, recalculate slc following opt_frz=1, smc remains
!! Case 3: frozen ==> unfrozen, melt all soil ice (if any)
!! Case 4: unfrozen ==> unfrozen along with other cases, (e.g., soil temp=tfrz),do nothing
!! Note: For Case 3, Yuan Xue thoroughly evaluated a total of four options and
!! current option is found to be the best as of 11/09/2023
!! @param[in] isot Integer code for the soil type data set
!! @param[in] ivegsrc Integer code for the vegetation type data set
!! @param[in] lensfc Number of land points for this tile
!! @param[in] lsoil Number of soil layers
!! @param[in] lsoil_incr Number of soil layers (from top) to apply soil increments to
!! @param[in] mask Mask indicating surface type
!! @param[in] stc_bck Background soil temperature states
!! @param[in] stc_adj Analysis soil temperature states
!! @param[inout] smc_adj Analysis soil moisture states
!! @param[inout] slc_adj Analysis liquid soil moisture states
!! @param[in] stc_updated Integer to record whether STC in each grid cell was updated
!! @param[in] slc_updated Integer to record whether SLC in each grid cell was updated
!! @param[in] zsoil Depth of bottom of each soil layer
!! @author Clara Draper @date April 2021
!! updates April 2024, Tseganeh ZG: (calling from jedi-apply_lnd_inc) 
!! stc_upd, slc_upd user inputs
!! mpi init in calling rourine
!! applying only for NoahMP
!! @param[in] isoiltype (soil types) dtype=integer (directly using in jedi-apply_lnd_inc)

subroutine apply_land_da_adjustments_soil(lsoil_incr, isot, ivegsrc,lensfc, &
                 lsoil, isoiltype, mask, stc_bck, stc_adj, smc_adj, slc_adj, &
                 stc_updated, slc_updated, zsoil, upd_stc, upd_slc, myrank, print_summary, print_debug)

    use mpi
    use set_soilveg_snippet_mod, only: set_soilveg_noah,set_soilveg_noahmp
    use sflx_snippet,    only: frh2o

    implicit none
 
    integer, intent(in)           :: lsoil_incr, lensfc, lsoil, isot, ivegsrc, myrank
    integer, intent(in)           :: isoiltype(lensfc) ! soil types, integer 4.17.26
    integer, intent(in)           :: mask(lensfc)
    real, intent(in)              :: stc_bck(lensfc, lsoil)
    integer, intent(in)           :: stc_updated(lensfc), slc_updated(lensfc)
    real, intent(inout)           :: smc_adj(lensfc,lsoil), slc_adj(lensfc,lsoil) 
    real, intent(inout)           :: stc_adj(lensfc, lsoil)
    real(kind=4), intent(in)      :: zsoil(lsoil)
    logical, intent(in)           :: upd_slc, upd_stc, print_summary, print_debug

    logical                       :: frzn_bck, frzn_anl
    logical                       :: soil_freeze, soil_ice

    integer                       :: i, l, n_freeze, n_thaw, ierr
    integer                       :: soiltype, iret, n_stc, n_slc
    
    real                          :: slc_new

    real, parameter               :: tfreez=273.16 !< con_t0c  in physcons
    real, dimension(30)           :: maxsmc, bb, satpsi
    real, dimension(4)            :: dz ! layer thickness

    real, parameter          :: hfus=0.3336e06 !< latent heat of fusion(j/kg)
    real, parameter          :: grav=9.80616   !< gravity accel.(m/s2)
    real                     :: smp !< for computing supercooled water 


    if (upd_stc) then

      call set_soilveg_noahmp(isot, ivegsrc, maxsmc, bb, satpsi, iret)
      if (iret < 0) then
           print *, 'FATAL ERROR: problem in set_soilveg_noahmp'
           call mpi_abort(mpi_comm_world, 10, ierr)
      endif

      n_stc = 0
      n_slc = 0

      do i=1,lensfc
        if (stc_updated(i) == 1 ) then 
            n_stc = n_stc+1
            soiltype = isoiltype(i)    
            do l = 1, lsoil_incr
               !case 1: frz ==> frz, recalculate slc, smc remains
               !case 2: unfrz ==> frz, recalculate slc, smc remains
               !both cases are considered in the following if case
               if (stc_adj(i,l) .LT. tfreez )then
                  !recompute supercool liquid water,smc_anl remain unchanged
                  smp = hfus*(tfreez-stc_adj(i,l))/(grav*stc_adj(i,l)) !(m)
                  slc_new=maxsmc(soiltype)*(smp/satpsi(soiltype))**(-1./bb(soiltype))
                  slc_adj(i,l) = max( min( slc_new, smc_adj(i,l)), 0.0 )
               endif
               !case 3: frz ==> unfrz, melt all soil ice (if any)
               if (stc_adj(i,l) .GT. tfreez )then !do not rely on stc_bck
                  slc_adj(i,l)=smc_adj(i,l)
               endif
            enddo
        endif
      enddo

    endif  

    if (upd_slc) then

      dz(1) = -zsoil(1)
      do l = 2,lsoil 
          dz(l) = -zsoil(l) + zsoil(l-1) 
      enddo 
      print *, 'Applying soil moisture mins ' 

      do i=1,lensfc
      if (slc_updated(i) == 1 ) then 
          n_slc = n_slc+1
          ! apply SM bounds (later: add upper SMC limit)
          do l = 1, lsoil_incr
            ! noah-mp minimum is 1 mm per layer (in SMC)
            ! no need to maintain frozen amount, would be v. small.
            slc_adj(i,l) = max( 0.001/dz(l), slc_adj(i,l) )
            smc_adj(i,l) = max( 0.001/dz(l), smc_adj(i,l) )
          enddo
       endif
      enddo
    endif

    if (print_summary .and. myrank==0) then
      write(*,'(a,i2)') 'statistics of grids number processed for rank : ', myrank
      write(*,'(a,i8)') ' soil grid total', lensfc
      write(*,'(a,i8)') ' soil grid cells with slc update', n_slc
      write(*,'(a,i8)') ' soil grid cells with stc update', n_stc
    endif

    if (print_debug) then 
      write(*,'(a,i2)') 'statistics of grids number processed for rank : ', myrank 
      write(*,'(a,i8)') ' soil grid total', lensfc
      write(*,'(a,i8)') ' soil grid cells with slc update', n_slc
      write(*,'(a,i8)') ' soil grid cells with stc update', n_stc
    endif

end subroutine apply_land_da_adjustments_soil

end module soil_increments
