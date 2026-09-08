program apply_incr_noahmp_soil

 use netcdf

 use soil_increments, only : add_increment_soil, calculate_landinc_mask, apply_land_da_adjustments_soil

 implicit none

 include 'mpif.h'


 type noahmp_type
   double precision, allocatable :: stc_bkg     (:,:)
   double precision, allocatable :: stc         (:,:)
   double precision, allocatable :: smc         (:,:)
   double precision, allocatable :: slc         (:,:)
   double precision, allocatable :: stc_inc     (:,:)
   double precision, allocatable :: slc_inc     (:,:)
   integer, allocatable          :: stc_updated        (:)
   integer, allocatable          :: slc_updated        (:)
   integer, allocatable          :: soilsnow_tile      (:)
   double precision, allocatable :: swe                (:)
   integer, allocatable          :: vtype              (:)
   integer, allocatable          :: stype              (:)
 end type noahmp_type 

 type(noahmp_type)               :: noahmp_state

 integer          :: res, len_land_vec 
 character(len=8) :: date_str 
 character(len=2) :: hour_str

 ! index to map between tile and vector space 
 integer, allocatable          :: tile2vector(:,:) 

 integer :: ierr, irank, nprocs, myrank, lunit, ncid, n, nn
 integer :: ntiles, ens_size, ens_mem, tile_num
 character(len=3) :: ens_str
 logical :: file_exists

 character(len=512) :: orog_path, rst_path_full, inc_path_full, static_filename
 character(len=256) :: rst_path, inc_path
 character(len=20)  :: otype ! orography filename stub. For atm only, oro_C${RES}, for atm/ocean oro_C${RES}.mx100

 character(len=512) :: restart_file
 character(len=1)   :: tilech
 character(len=512) :: ioerrmsg

 logical            :: print_summary, print_debug, frac_grid

 double precision   :: fice_threshold, lfrac_threshold

 integer            :: lsoil, lsoil_incr, isot, ivegsrc
 integer            :: veg_type_landice   !, lsm
 character(len=256) :: inc_prefix, stype_prefix
 double precision, allocatable :: styper(:,:)
 logical                       :: upd_stc
 logical                       :: upd_slc
 logical                       :: csg_increment   ! if true, read increments from cube sphere file instead of fv3 increment file. 
 logical                       :: vector_grids    ! if true, increments and bkg in vector grids
 real(kind=4)       :: zsoil(4) = (/ -0.1, -0.4, -1.0, -2.0 /)   

 namelist /noahmp_soil/ date_str, hour_str, res, rst_path, inc_path, orog_path, otype, ntiles, ens_size, &
                        print_summary, print_debug, lsoil_incr, inc_prefix, stype_prefix, &
                        upd_stc, upd_slc, csg_increment, vector_grids, static_filename          

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    if (myrank==0) print*, "starting apply_incr_noahmp_soil program on ", nprocs, " procs"

    ! SET NAMELIST DEFAULTS
    rst_path = './'
    inc_path = './'
    ntiles = 6
    ens_size = 1
    print_summary = .true.
    print_debug = .false.
    lsoil_incr = 2
    upd_stc = .false.
    upd_slc = .false.
    inc_prefix = ""
    csg_increment = .false.
    vector_grids = .false.
    static_filename = ""
    !stype_prefix = "C96.mx100.soil_type"
    
    ! hard coded defaults--unlikely to change
    frac_grid = .true.
    fice_threshold=0.0
    lfrac_threshold=0.0001
    veg_type_landice = 15
    lsoil = 4     ! zsoil is hard-coded for 4 layers
    ivegsrc = 1   ! The NOAHMP LSM expects that the ivegsrc physics parameter is 1
    isot = 1      ! Noahmp expects 1

    ! READ NAMELIST 
     inquire (file='apply_incr_nml', exist=file_exists) 

    if (.not. file_exists) then
        write (6, *) 'ERROR: apply_incr_nml does not exist'
        call mpi_abort(mpi_comm_world, 10)  
    end if

    open (action='read', file='apply_incr_nml', iostat=ierr, newunit=lunit, iomsg=ioerrmsg)
    read (nml=noahmp_soil, iostat=ierr, unit=lunit)
    close (lunit)
    if (ierr /= 0) then
        print*, "Error code from namelist read", ierr
        write(6,*) trim(ioerrmsg)         
        call mpi_abort(mpi_comm_world, 10)  
    end if

    if (vector_grids) then

        if (myrank==0) print*, "vector_grids: both increments and bkg"

        do irank=myrank, ens_size - 1, nprocs

            write(ens_str, '(I3.3)') irank
            !TBCL: keep the default for ens_size=1
            if(ens_size > 1) then 
                rst_path_full = trim(rst_path)//"/mem"//ens_str//"/"
                inc_path_full = trim(inc_path)//"/mem"//ens_str//"/"
            else
                rst_path_full = trim(rst_path)      
                inc_path_full = trim(inc_path)      
            endif

            call read_vector_restart(rst_path_full, static_filename, date_str//hour_str, noahmp_state, &
                               len_land_vec, lsoil, lsoil_incr)
            noahmp_state%stc_bkg = noahmp_state%stc

            call read_vector_increment(inc_path_full, inc_prefix, len_land_vec, &
                 noahmp_state, lsoil_incr, upd_stc, upd_slc)

            call calculate_landinc_mask(noahmp_state%swe,noahmp_state%vtype,noahmp_state%stype,&
                len_land_vec, veg_type_landice, noahmp_state%soilsnow_tile) 

            call add_increment_soil(lsoil_incr,noahmp_state%stc_inc,noahmp_state%slc_inc, &
               noahmp_state%stc,noahmp_state%smc,noahmp_state%slc,&
               noahmp_state%stc_updated,noahmp_state%slc_updated,noahmp_state%soilsnow_tile,noahmp_state%soilsnow_tile,&
               len_land_vec,lsoil,myrank, upd_stc, upd_slc, print_summary, print_debug)

            call apply_land_da_adjustments_soil(lsoil_incr, isot, ivegsrc, len_land_vec, &
                 lsoil, noahmp_state%stype, noahmp_state%soilsnow_tile,noahmp_state%stc_bkg, &
                 noahmp_state%stc,noahmp_state%smc,noahmp_state%slc, &
                 noahmp_state%stc_updated,noahmp_state%slc_updated, zsoil, upd_stc, upd_slc, myrank, print_summary, print_debug)
            
            ! WRITE OUT ADJUSTED RESTART
            call write_vector_restart(rst_path_full, date_str//hour_str, noahmp_state, &
                               len_land_vec, lsoil)

            if (allocated(noahmp_state%stc_bkg )) deallocate(noahmp_state%stc_bkg           )
            if (allocated(noahmp_state%stc     )) deallocate(noahmp_state%stc               )
            if (allocated(noahmp_state%smc     )) deallocate(noahmp_state%smc               )
            if (allocated(noahmp_state%slc     )) deallocate(noahmp_state%slc               )
            if (allocated(noahmp_state%stc_inc )) deallocate(noahmp_state%stc_inc           )
            if (allocated(noahmp_state%slc_inc )) deallocate(noahmp_state%slc_inc           )
            if (allocated(noahmp_state%stc_updated )) deallocate(noahmp_state%stc_updated       )
            if (allocated(noahmp_state%slc_updated )) deallocate(noahmp_state%slc_updated       )
            if (allocated(noahmp_state%soilsnow_tile )) deallocate(noahmp_state%soilsnow_tile     )
            if (allocated(noahmp_state%swe )) deallocate(noahmp_state%swe               )
            if (allocated(noahmp_state%vtype )) deallocate(noahmp_state%vtype             )
            if (allocated(noahmp_state%stype )) deallocate(noahmp_state%stype             )

        enddo 
         
        if (myrank==0) print*, "apply_incr_noahmp_soil finishing"
        call mpi_finalize(ierr)

        return 

    endif

    allocate(styper(res, res))

    do irank=myrank, ntiles*ens_size - 1, nprocs
        ens_mem = irank/ntiles + 1            !ensemble member
        tile_num = MOD(irank, ntiles) + 1      !tile number

        write(ens_str, '(I3.3)') ens_mem

!TBCL: keep the default for ens_size=1
        if(ens_size > 1) then 
            rst_path_full = trim(rst_path)//"/mem"//ens_str//"/"
            inc_path_full = trim(inc_path)//"/mem"//ens_str//"/"
        else
            rst_path_full = trim(rst_path)      
            inc_path_full = trim(inc_path)      
        endif
        
        ! Calculate MAPPING INDEX based on land fraction
        call get_fv3_mapping_lfrac(tile_num, rst_path_full, date_str, hour_str, res, &
             orog_path,otype,frac_grid,lfrac_threshold,fice_threshold,len_land_vec,tile2vector,stype_prefix,styper)

        ! SET-UP THE NOAH-MP STATE  AND INCREMENT        
        ! The allocations are inside the loop because different ensemble members could have different len_land_vec
        allocate(noahmp_state%stc_bkg         (len_land_vec, lsoil)) 
        allocate(noahmp_state%stc             (len_land_vec, lsoil))
        allocate(noahmp_state%smc             (len_land_vec, lsoil))  
        allocate(noahmp_state%slc             (len_land_vec, lsoil)) 
        allocate(noahmp_state%stc_inc         (len_land_vec, lsoil_incr))
        allocate(noahmp_state%slc_inc         (len_land_vec, lsoil_incr))
        !allocate(noahmp_state%smc_inc         (len_land_vec, lsoil_incr))
        allocate(noahmp_state%stc_updated     (len_land_vec))
        allocate(noahmp_state%slc_updated     (len_land_vec))
        allocate(noahmp_state%soilsnow_tile   (len_land_vec))
        !allocate(noahmp_state%soilsnow_fg_tile(len_land_vec))
        !allocate(noahmp_state%land_frac       (len_land_vec)) 
        allocate(noahmp_state%swe             (len_land_vec))
        allocate(noahmp_state%vtype           (len_land_vec))
        allocate(noahmp_state%stype           (len_land_vec))
       
        ! map soil types
        do nn=1,len_land_vec
            noahmp_state%stype = nint(styper(tile2vector(nn,1), tile2vector(nn,2)))
        enddo

        ! READ RESTART FILE 
        write(tilech, '(i1.1)') (tile_num)
        restart_file = trim(rst_path_full)//"/"//date_str//"."//hour_str//"0000.sfc_data.tile"//tilech//".nc"

        call read_fv3_restart(trim(restart_file), res, ncid, len_land_vec, tile2vector, noahmp_state, lsoil)
        noahmp_state%stc_bkg = noahmp_state%stc

        ! READ soil DA increments
        if (csg_increment) then
             call read_csg_increment(tile_num, inc_path_full, inc_prefix, res, &
                len_land_vec, tile2vector, noahmp_state, lsoil_incr, upd_stc, upd_slc)
        else
             call read_fv3_increment(tile_num, inc_path_full, date_str, hour_str, res, &
                    len_land_vec, tile2vector, inc_prefix, noahmp_state, lsoil_incr, upd_stc, upd_slc)
        endif

        call calculate_landinc_mask(noahmp_state%swe,noahmp_state%vtype,noahmp_state%stype,&
                len_land_vec, veg_type_landice, noahmp_state%soilsnow_tile)                ! soilsnow_fg_tile

        call add_increment_soil(lsoil_incr,noahmp_state%stc_inc,noahmp_state%slc_inc, &
               noahmp_state%stc,noahmp_state%smc,noahmp_state%slc,&
               noahmp_state%stc_updated,noahmp_state%slc_updated,noahmp_state%soilsnow_tile,noahmp_state%soilsnow_tile,&
               len_land_vec,lsoil,myrank, upd_stc, upd_slc, print_summary, print_debug)
        
       !call calculate_landinc_mask(noahmp_state%swe,noahmp_state%vtype,noahmp_state%stype,&
       !         len_land_vec, veg_type_landice, noahmp_state%soilsnow_tile)

        call apply_land_da_adjustments_soil(lsoil_incr, isot, ivegsrc, len_land_vec, &
                 lsoil, noahmp_state%stype, noahmp_state%soilsnow_tile,noahmp_state%stc_bkg, &
                 noahmp_state%stc,noahmp_state%smc,noahmp_state%slc, &
                 noahmp_state%stc_updated,noahmp_state%slc_updated, zsoil, upd_stc, upd_slc, myrank, print_summary, print_debug)
            
        ! WRITE OUT ADJUSTED RESTART
        call write_fv3_restart(trim(restart_file),noahmp_state,res,ncid,len_land_vec,tile2vector,lsoil) 
        
        ! CLOSE RESTART FILE 
        ierr = nf90_close(ncid)
        call netcdf_err( ierr, "closing restart file "//trim(restart_file) )
        
        ! Deallocate. These are required incase a single process loops through multiple tiles with different mapping     
        if (allocated(tile2vector)) deallocate(tile2vector)   
                
        deallocate(noahmp_state%stc_bkg           )
        deallocate(noahmp_state%stc               )
        deallocate(noahmp_state%smc               )
        deallocate(noahmp_state%slc               )
        deallocate(noahmp_state%stc_inc           )
        deallocate(noahmp_state%slc_inc           )
        !deallocate(noahmp_state%smc_inc          )
        deallocate(noahmp_state%stc_updated       )
        deallocate(noahmp_state%slc_updated       )
        deallocate(noahmp_state%soilsnow_tile     )
        !deallocate(noahmp_state%soilsnow_fg_tile )
        !deallocate(noahmp_state%land_frac        )
        deallocate(noahmp_state%swe               )
        deallocate(noahmp_state%vtype             )
        deallocate(noahmp_state%stype             )

    enddo

    deallocate(styper)

    if (myrank==0) print*, "apply_incr_noahmp_soil finishing"
    call mpi_finalize(ierr)

 contains 

!--------------------------------------------------------------
! if at netcdf call returns an error, print out a message and stop processing.
!--------------------------------------------------------------
 subroutine netcdf_err( err, string )

        implicit none

        include 'mpif.h'

        integer, intent(in) :: err
        character(len=*), intent(in) :: string
        character(len=80) :: errmsg
        integer           :: ierr

        if( err == nf90_noerr )return
        errmsg = nf90_strerror(err)
        print*,''
        print*,'fatal error: ', trim(string), ': ', trim(errmsg)
        print*,'stop.'
        call mpi_abort(mpi_comm_world, 999)

        return
 end subroutine netcdf_err

  subroutine write_vector_restart(vector_restart_path, restart_date, noahmp_state, &
                               vector_length, lsoil)
  
    use netcdf

    type(noahmp_type)   :: noahmp_state
    character(*)        :: vector_restart_path
    character*10        :: restart_date
    integer             :: vector_length, lsoil

    character*19        :: date
    character*256       :: vector_filename, filename
    integer             :: ncid, dimid, varid, status
    logical             :: file_exists

    character*19        :: date
    integer             :: yyyy,mm,dd,hh,nn,ss
    integer             :: itile, ix, iy, iloc

    read(restart_date( 1: 4),'(i4.4)') yyyy
    read(restart_date( 5: 6),'(i2.2)') mm
    read(restart_date( 7: 8),'(i2.2)') dd
    read(restart_date( 9:10),'(i2.2)') hh

    write(date,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
        yyyy, "-", mm, "-", dd, "_", hh, "-00-00"  !, nn, "-", ss
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Create vector file name
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(vector_filename,'(a17,a19,a3)') "ufs_land_restart.", date, ".nc"

    filename = trim(vector_restart_path)//trim(vector_filename)
    
    inquire(file=filename, exist=file_exists)
    
    if(.not.file_exists) then 
        print*, trim(filename), " does not exist2"
        print*, "Check paths and file name"
        stop 10 
    end if
        
    print*, "Writing vector file: ", trim(filename)

    status = nf90_open(filename, nf90_write, ncid)
    if (status /= nf90_noerr) then 
        call handle_err(status) 
    endif

    status = nf90_inq_varid(ncid, "temperature_soil", varid)
    if (status /= nf90_noerr) then 
            print *, 'temperature_soil variable missing from vector file' 
            call handle_err(status) 
    endif
    status = nf90_put_var(ncid, varid , noahmp_state%stc , &
        start = (/1            , 1, 1/)                , &
        count = (/vector_length, lsoil, 1/))
    if (status /= nf90_noerr) then  
        call handle_err(status) 
    endif

    status = nf90_inq_varid(ncid, "soil_liquid_vol", varid)
    if (status /= nf90_noerr) then 
        print *, 'soil_liquid_vol variable missing from vector file' 
        call handle_err(status) 
    endif
    status = nf90_put_var(ncid, varid , noahmp_state%slc , &
        start = (/1            , 1, 1/)                , &
        count = (/vector_length, lsoil/))
    if (status /= nf90_noerr) then  
        call handle_err(status) 
    endif
   
    status = nf90_inq_varid(ncid, "soil_moisture_vol", varid)
    if (status /= nf90_noerr) then 
        print *, 'soil_moisture_vol variable missing from vector file' 
        call handle_err(status) 
    endif
    status = nf90_put_var(ncid, varid , noahmp_state%smc, &
        start = (/1            , 1, 1/)                , &
        count = (/vector_length, lsoil, 1/))
    if (status /= nf90_noerr) then  
        call handle_err(status) 
    endif
    
    status = nf90_close(ncid)
    if (status /= nf90_noerr) then  
        call handle_err(status) 
    endif
    
  end subroutine write_vector_restart 

  subroutine read_vector_restart(vector_restart_path, static_filename, restart_date, noahmp_state, &
                               vector_length, lsoil, lsoil_incr)
  
    use netcdf

    type(noahmp_type)   :: noahmp_state
    character(*)        :: vector_restart_path, static_filename
    character*10        :: restart_date
    integer             :: vector_length

    character*19        :: date
    character*256       :: vector_filename, filename
    integer             :: ncid, dimid, varid, status
    logical             :: file_exists

    character*19        :: date
    integer             :: yyyy,mm,dd,hh,nn,ss
    integer             :: itile, ix, iy, iloc
    integer             :: ncid, dimid, varid, status
    logical             :: file_exists
    read(restart_date( 1: 4),'(i4.4)') yyyy
    read(restart_date( 5: 6),'(i2.2)') mm
    read(restart_date( 7: 8),'(i2.2)') dd
    read(restart_date( 9:10),'(i2.2)') hh

    write(date,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
        yyyy, "-", mm, "-", dd, "_", hh, "-00-00"  !, nn, "-", ss
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Create vector file name
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(vector_filename,'(a17,a19,a3)') "ufs_land_restart.", date, ".nc"

    filename = trim(vector_restart_path)//trim(vector_filename)
    
    inquire(file=filename, exist=file_exists)
    
    if(.not.file_exists) then 
        print*, trim(filename), " does not exist2"
        print*, "Check paths and file name"
        stop 10 
    end if
    
    print*, "Reading vector file: ", trim(filename)

    status = nf90_open(filename, NF90_NOWRITE, ncid)

    status = nf90_inq_dimid(ncid, "location", dimid)
    status = nf90_inquire_dimension(ncid, dimid, len = vector_length)

    allocate(noahmp_state%stc_bkg         (vector_length, lsoil)) 
    allocate(noahmp_state%stc             (vector_length, lsoil))
    allocate(noahmp_state%smc             (vector_length, lsoil))  
    allocate(noahmp_state%slc             (vector_length, lsoil)) 
    allocate(noahmp_state%stc_inc         (vector_length, lsoil_incr))
    allocate(noahmp_state%slc_inc         (vector_length, lsoil_incr))
    allocate(noahmp_state%stc_updated     (vector_length))
    allocate(noahmp_state%slc_updated     (vector_length))
    allocate(noahmp_state%soilsnow_tile   (vector_length))
    allocate(noahmp_state%swe             (vector_length))
    allocate(noahmp_state%vtype           (vector_length))
    allocate(noahmp_state%stype           (vector_length))
  
    ! Read the vector fields

    status = nf90_inq_varid(ncid, "snow_water_equiv", varid)
    if (status /= nf90_noerr) then 
            print *, 'snow_water_equiv variable missing from vector file' 
            call handle_err(status) 
    endif
    status = nf90_get_var(ncid, varid , noahmp_state%swe   , &
        start = (/1,1/), count = (/vector_length, 1/))

    
    status = nf90_inq_varid(ncid, "temperature_soil", varid)
    if (status /= nf90_noerr) then 
            print *, 'temperature_soil variable missing from vector file' 
            call handle_err(status) 
    endif
    status = nf90_get_var(ncid, varid , noahmp_state%stc , &
        start = (/1            , 1, 1/)                , &
        count = (/vector_length, lsoil, 1/))

    status = nf90_inq_varid(ncid, "soil_moisture_vol", varid)
    if (status /= nf90_noerr) then 
            print *, 'soil_moisture_vol variable missing from vector file' 
            call handle_err(status) 
    endif
    status = nf90_get_var(ncid, varid , noahmp_state%smc, &
        start = (/1            , 1, 1/)                , &
        count = (/vector_length, lsoil, 1/))

    status = nf90_inq_varid(ncid, "soil_liquid_vol", varid)
    if (status /= nf90_noerr) then
            print *, 'soil_liquid_vol variable missing from vector file'
            call handle_err(status)
    endif
    status = nf90_get_var(ncid, varid , noahmp_state%slc , &
        start = (/1            , 1, 1/)                , &
        count = (/vector_length, lsoil, 1/))

    status = nf90_close(ncid)
    if (status /= nf90_noerr) then
            call handle_err(status)
    endif

    ! read vegetation and soil type from static file 

    filename = trim(vector_restart_path)//trim(static_filename)
    
    inquire(file=filename, exist=file_exists)
    
    if(.not.file_exists) then 
        print*, trim(filename), " does not exist3"
        print*, "Check paths and file name"
        stop 10 
    end if
    
    status = nf90_open(filename, NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid, "vegetation_category", varid)
    if (status /= nf90_noerr) then 
            print *, 'vegetation_category missing from static file' 
            call handle_err(status) 
    endif
    status = nf90_get_var(ncid, varid , noahmp_state%vtype, &
        start = (/1,1/), count = (/vector_length, 1/))

    status = nf90_inq_varid(ncid, "soil_category", varid)
    if (status /= nf90_noerr) then 
            print *, 'soil_category missing from static file' 
            call handle_err(status) 
    endif
    status = nf90_get_var(ncid, varid , noahmp_state%stype, &
        start = (/1,1/), count = (/vector_length, 1/))
    if (status /= nf90_noerr) then
            call handle_err(status)
    endif
        
    status = nf90_close(ncid)
    if (status /= nf90_noerr) then
            call handle_err(status)
    endif

  end subroutine read_vector_restart

!--------------------------------------------------------------
!  read in soil increments from cube sphere file
!--------------------------------------------------------------
 subroutine read_vector_increment(inc_path, inc_prefix, len_land_vec, &
                 noahmp_state, lsoil_incr, upd_stc, upd_slc)

    implicit none 

    include 'mpif.h'

    integer, intent(in) :: len_land_vec, lsoil_incr
    character(len=*), intent(in) :: inc_path, inc_prefix 
    type(noahmp_type), intent(inout)  :: noahmp_state
    logical, intent(in)               :: upd_stc, upd_slc

    character(len=512) :: incr_file
    logical :: file_exists
    integer :: ierr 
    integer :: id_dim, id_var, fres, ncid
    integer :: nn, nl
    character(len=20) :: var_name
    character(len=1)  :: layerch

    ! OPEN FILE
    incr_file = trim(inc_path)//"/"//trim(inc_prefix)//".nc"

    inquire(file=trim(incr_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'incr_file does not exist, ', &
                    trim(incr_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ierr=nf90_open(trim(incr_file),nf90_nowrite,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(incr_file) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'location', id_dim)
    call netcdf_err(ierr, 'reading location dim from '//trim(incr_file) )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading location length from '//trim(incr_file) )

    if ( fres /= len_land_vec) then
       print*,'fatal error: dimension ',fres, ' in '//trim(incr_file), ' not equal to ',len_land_vec, ' in vector restart'
       call mpi_abort(mpi_comm_world, ierr)
    endif
 
    noahmp_state%stc_inc = 0.0  !0 if no inc exists. TODO: need to do liau type "no update on 0"?
    noahmp_state%slc_inc = 0.0
    
    if (upd_stc) then
        status = nf90_inq_varid(ncid, "stc", varid)
        if (status /= nf90_noerr) then 
                print *, 'temperature_soil variable missing from increment file' 
                call handle_err(status) 
        endif
        status = nf90_get_var(ncid, varid , noahmp_state%stc_inc , &
            start = (/1            , 1, 1/)                , &
            count = (/len_land_vec, lsoil_incr, 1/))
    endif

    if (upd_slc) then
        status = nf90_inq_varid(ncid, "slc", varid)
        if (status /= nf90_noerr) then 
                print *, 'soil_liquid_vol variable missing from increment file' 
                call handle_err(status) 
        endif
        status = nf90_get_var(ncid, varid , noahmp_state%slc_inc , &
            start = (/1            , 1, 1/)                , &
            count = (/len_land_vec, lsoil_incr, 1/))

    endif

    ierr=nf90_close(ncid)
    call netcdf_err(ierr, 'closing file: '//trim(incr_file) )

end subroutine read_vector_increment  

!--------------------------------------------------------------
! create index for mapping from tiles (FV3 UFS restart) to vector
! of land locations (offline Noah-MP restart) based on fraction of land 
! (land_frac) field from the oro_grid files.
! !> mask = 1 (land) if: land frac >= lfrac_threshold = 0.01%
!                        && fice (sea ice) not > fice_threshold = 0.0
!                        && veg type not 15 (land ice)
!
! Note: These masks do NOT have exclusion of islands. 
!--------------------------------------------------------------

 subroutine get_fv3_mapping_lfrac(tile_num, rst_path, date_str, hour_str, res, & 
            orog_path, otype, fice_grid, lfrac_thold, fice_fhold, len_land_vec, tile2vector,stype_prefix,styper)

 implicit none 

 include 'mpif.h'

 integer, intent(in)          :: tile_num, res
 character(len=*), intent(in) :: rst_path
 character(len=8), intent(in) :: date_str 
 character(len=2), intent(in) :: hour_str 
 character(len=*), intent(in)   :: orog_path, stype_prefix
 character(len=20), intent(in)  :: otype
 logical, intent(in)            :: fice_grid
 double precision, intent(in)      :: lfrac_thold, fice_fhold
 integer, intent(out)              :: len_land_vec
 integer, allocatable, intent(out) :: tile2vector(:,:)
 double precision, intent(out)     :: styper(res,res) 

 character(len=512) :: restart_file, filename
 character(len=1)   :: rankch
 logical :: file_exists
 integer :: ierr, ncid
 integer :: id_dim, id_var, fres

 integer            :: slmsk_rest(res,res), slmsk_lfrac(res,res)
 double precision   :: fice(res,res)
 double precision   :: vtype(res,res)     ! stored as double in restart files
 double precision   :: land_frac(res,res)
 integer, parameter :: vtype_landice=15   !, vtype_water=17
 integer            :: i, j, nn, len_land_vec_rest, diff_count
 integer, allocatable  :: tile2vector_rest(:,:), tile2vector_diff(:,:)

    ! OPEN FILE
    write(rankch, '(i1.1)') (tile_num)
    filename =trim(orog_path)//"/"//trim(otype)//".tile"//rankch//".nc"

    inquire(file=trim(filename), exist=file_exists)

    if (.not. file_exists) then
            print *, 'filename does not exist, ', &
                    trim(filename) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ierr=nf90_open(trim(filename),nf90_nowrite,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(filename) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'lon', id_dim)
    call netcdf_err(ierr, 'reading lon id from '//trim(filename) )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading lon from '//trim(filename) )

    if ( fres /= res) then
       print*,'fatal error: dimensions wrong in file '//trim(filename)
       call mpi_abort(mpi_comm_world, ierr)
    endif

    ! READ land frac 
    ierr=nf90_inq_varid(ncid, "land_frac", id_var)
    call netcdf_err(ierr, 'reading land_frac id' )
    ierr=nf90_get_var(ncid, id_var, land_frac)
    call netcdf_err(ierr, 'reading land_frac' )

    ! close file 
    ierr=nf90_close(ncid)
    call netcdf_err(ierr, 'closing file: '//trim(filename) )

    ! read soil type 
    filename =trim(orog_path)//"/sfc/"//trim(stype_prefix)//".tile"//rankch//".nc"

    inquire(file=trim(filename), exist=file_exists)

    if (.not. file_exists) then
            print *, 'filename does not exist, ', &
                    trim(filename) , ' exiting'
            call mpi_abort(mpi_comm_world, 10)
    endif

    ierr=nf90_open(trim(filename),nf90_nowrite,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(filename) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'nx', id_dim)
    call netcdf_err(ierr, 'reading nx id from '//trim(filename) )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading nx from '//trim(filename) )

    if ( fres /= res) then
       print*,'fatal error: dimensions wrong in file '//trim(filename)
       call mpi_abort(mpi_comm_world, ierr)
    endif

    ierr=nf90_inq_varid(ncid, "soil_type", id_var)
    call netcdf_err(ierr, 'reading soil_type id' )
    ierr=nf90_get_var(ncid, id_var, styper)
    call netcdf_err(ierr, 'reading soil_type' )

    ! close file
    ierr=nf90_close(ncid)
    call netcdf_err(ierr, 'closing file: '//trim(filename) )

    ! Use vtype to exclude glaciers
    ! OPEN FILE
    write(rankch, '(i1.1)') (tile_num)
    restart_file = trim(rst_path)//"/"//date_str//"."//hour_str//"0000.sfc_data.tile"//rankch//".nc"

    inquire(file=trim(restart_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'restart_file does not exist, ', &
                    trim(restart_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ierr=nf90_open(trim(restart_file),nf90_write,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(restart_file) )
 
    ! READ MASK from restart
    ierr=nf90_inq_varid(ncid, "slmsk", id_var)
    call netcdf_err(ierr, 'reading slmsk id' )
    ierr=nf90_get_var(ncid, id_var, slmsk_rest)
    call netcdf_err(ierr, 'reading slmsk from restart' )

    ! REMOVE GLACIER GRID POINTS
    ierr=nf90_inq_varid(ncid, "vtype", id_var)
    call netcdf_err(ierr, 'reading vtype id' )
    ierr=nf90_get_var(ncid, id_var, vtype)
    call netcdf_err(ierr, 'reading vtype' )

    if (fice_grid) then    
      ierr=nf90_inq_varid(ncid, "fice", id_var)
      call netcdf_err(ierr, 'reading fice id' )
      ierr=nf90_get_var(ncid, id_var, fice)
      call netcdf_err(ierr, 'reading fice' )   
    endif

    ! close file
    ierr=nf90_close(ncid)
    call netcdf_err(ierr, 'closing file: '//trim(restart_file) )



    slmsk_lfrac = 0
    do i = 1, res
        do j = 1, res
            if ( land_frac(i,j) >= lfrac_thold ) slmsk_lfrac(i,j) = 1
        enddo
    enddo

    ! remove land grid cells if ice is present
    if (fice_grid) then
        write (6, *) 'ammending mask to exclude sea ice from', trim(restart_file)
        do i = 1, res
            do j = 1, res
                if (fice(i,j) > fice_fhold ) then 
                        slmsk_lfrac(i,j) = 0
                        slmsk_rest(i,j) = 0
                endif
            enddo
        enddo
    endif

    ! remove land grid cells if glacier land type
    do i = 1, res
        do j = 1, res
            if ( nint(vtype(i,j)) ==  vtype_landice) then  ! vtype is integer, but stored as double
                    slmsk_lfrac(i,j) = 0 
                    slmsk_rest(i,j) = 0
            endif
        enddo
    enddo

    ! get number of land points
    len_land_vec = 0
    do i = 1, res 
        do j = 1, res 
             if ( slmsk_lfrac(i,j) == 1)  len_land_vec = len_land_vec + 1  
        enddo 
    enddo
    
    len_land_vec_rest = 0
    do i = 1, res
        do j = 1, res
             if ( slmsk_rest(i,j) > 0)  len_land_vec_rest = len_land_vec_rest + 1
        enddo
    enddo

    if (len_land_vec .ne. len_land_vec_rest) then 
        print*, "number of land points from "//trim(filename)//" not consitent with those from "//trim(restart_file)
        print*, "orog land points = ",len_land_vec," restart land points = ", len_land_vec_rest 
        call mpi_abort(mpi_comm_world, 10)
    endif

    allocate(tile2vector(len_land_vec,2)) 
    allocate(tile2vector_rest(len_land_vec,2))
    allocate(tile2vector_diff(len_land_vec,2))

    nn=0
    do i = 1, res 
        do j = 1, res 
             if ( slmsk_lfrac(i,j) == 1)   then 
                nn=nn+1
                tile2vector(nn,1) = i 
                tile2vector(nn,2) = j 
             endif
        enddo 
    enddo

    nn=0
    do i = 1, res
        do j = 1, res
             if ( slmsk_rest(i,j) > 0) then  ! some land points are marked 2 (snow) during ufs calls
                nn=nn+1
                tile2vector_rest(nn,1) = i
                tile2vector_rest(nn,2) = j
             endif
        enddo
    enddo

    !check mask consistency from restart and orog
    tile2vector_diff = tile2vector_rest - tile2vector
    diff_count = count(abs(tile2vector_diff) > 0)
    if (diff_count > 0) then 
        print*, diff_count, " differences between land mask from "//trim(filename)//" and from "//trim(restart_file)
        call mpi_abort(mpi_comm_world, 10)
    endif

    deallocate(tile2vector_rest, tile2vector_diff)

end subroutine get_fv3_mapping_lfrac

!--------------------------------------------------------------
! open fv3 restart, and read in required variables
! file is opened as read/write and remains open
!--------------------------------------------------------------
 subroutine read_fv3_restart(restart_file, res, ncid, & 
                len_land_vec,tile2vector, noahmp_state, lsoil)

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: res, len_land_vec   !tile_num, 
 character(len=*), intent(in) :: restart_file
 integer, intent(in)  :: tile2vector(len_land_vec,2)
 integer, intent(in)  :: lsoil
 integer, intent(out) :: ncid
 type(noahmp_type), intent(inout)  :: noahmp_state

 logical :: file_exists
 integer :: ierr, id_dim, fres
 integer :: nn
 double precision :: vstyper(len_land_vec)

    ! OPEN FILE
    inquire(file=trim(restart_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'restart_file does not exist, ', &
                    trim(restart_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ierr=nf90_open(trim(restart_file),nf90_write,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(restart_file) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'xaxis_1', id_dim)
    call netcdf_err(ierr, 'reading xaxis_1 in '//trim(restart_file) )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading xaxis_1 in '//trim(restart_file) )

    if ( fres /= res) then
       print*,'fatal error: dimensions wrong in '//trim(restart_file)
       call mpi_abort(mpi_comm_world, ierr)
    endif

  
    call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 0, & 
                        'weasdl', noahmp_state%swe)
    
    call read_nc_var3D(ncid, trim(restart_file), len_land_vec, res, lsoil,  tile2vector, &
                        'stc', noahmp_state%stc)

    call read_nc_var3D(ncid, trim(restart_file), len_land_vec, res, lsoil, tile2vector, & 
                        'smc', noahmp_state%smc)
 
    call read_nc_var3D(ncid, trim(restart_file), len_land_vec, res, lsoil, tile2vector, & 
                        'slc', noahmp_state%slc)

    call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 0, &
                        'vtype', vstyper)
    noahmp_state%vtype = nint(vstyper)
    
    !call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 0, &
    !                    'stype', vstyper)
    !noahmp_state%stype = nint(vstyper)

end subroutine read_fv3_restart


!--------------------------------------------------------------
!  read in soil increments from jedi increment file
!  file format is same as restart file
!--------------------------------------------------------------
 subroutine read_fv3_increment(tile_num, inc_path, date_str, hour_str, res, & 
                len_land_vec,tile2vector, inc_prefix, noahmp_state, lsoil_incr, upd_stc, upd_slc)

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: tile_num, res, len_land_vec, lsoil_incr
 character(len=*), intent(in) :: inc_path
 character(len=8), intent(in) :: date_str 
 character(len=2), intent(in) :: hour_str 
 integer, intent(in) :: tile2vector(len_land_vec,2)
 type(noahmp_type), intent(inout)  :: noahmp_state
 character(len=*), intent(in)      :: inc_prefix
 logical, intent(in)               :: upd_stc, upd_slc

 character(len=512) :: incr_file
 character(len=1) :: rankch
 logical :: file_exists
 integer :: ierr 
 integer :: id_dim, id_var, fres, ncid
 integer :: nn, nl

    ! OPEN FILE
    write(rankch, '(i1.1)') (tile_num)
    incr_file = trim(inc_path)//"/"//trim(inc_prefix)//date_str//"."//hour_str//"0000.sfc_data.tile"//rankch//".nc"

    inquire(file=trim(incr_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'incr_file does not exist, ', &
                    trim(incr_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ierr=nf90_open(trim(incr_file),nf90_nowrite,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(incr_file) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'xaxis_1', id_dim)
    call netcdf_err(ierr, 'reading xaxis_1 '//trim(incr_file) )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading xaxis_1 '//trim(incr_file) )

    if ( fres /= res) then
       print*,'fatal error: dimension fres ',fres, ' in '//trim(incr_file), ' not equal to res ',res
       call mpi_abort(mpi_comm_world, ierr)
    endif
 
    noahmp_state%stc_inc = 0.0  !0 if no inc exists. TODO: need to do liau type "no update on 0"?
    if (upd_stc) call read_nc_var3D(ncid, trim(incr_file), len_land_vec, res, lsoil_incr,  tile2vector, &
                        'stc   ', noahmp_state%stc_inc)

    noahmp_state%slc_inc = 0.0
    if (upd_slc) call read_nc_var3D(ncid, trim(incr_file), len_land_vec, res, lsoil_incr,  tile2vector, &
                        'slc   ', noahmp_state%slc_inc)

    ierr=nf90_close(ncid)
    call netcdf_err(ierr, 'closing file: '//trim(incr_file) )

end subroutine  read_fv3_increment

!--------------------------------------------------------------
!  read in soil increments from cube sphere file
!--------------------------------------------------------------
 subroutine read_csg_increment(tile_num, inc_path, inc_prefix, res, & 
                len_land_vec, tile2vector, noahmp_state, lsoil_incr, upd_stc, upd_slc)

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: tile_num, res, len_land_vec, lsoil_incr
 character(len=*), intent(in) :: inc_path, inc_prefix 
 integer, intent(in) :: tile2vector(len_land_vec,2)
 type(noahmp_type), intent(inout)  :: noahmp_state
 logical, intent(in)               :: upd_stc, upd_slc

 character(len=512) :: incr_file
 logical :: file_exists
 integer :: ierr 
 integer :: id_dim, id_var, fres, ncid
 integer :: nn, nl
 character(len=20) :: var_name
 character(len=1)  :: layerch

    ! OPEN FILE
    incr_file = trim(inc_path)//"/"//trim(inc_prefix)//".nc"

    inquire(file=trim(incr_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'incr_file does not exist, ', &
                    trim(incr_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ierr=nf90_open(trim(incr_file),nf90_nowrite,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(incr_file) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'grid_xt', id_dim)
    call netcdf_err(ierr, 'reading grid_xt '//trim(incr_file) )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading grid_xt '//trim(incr_file) )

    if ( fres /= res) then
       print*,'fatal error: dimension fres ',fres, ' in '//trim(incr_file), ' not equal to res ',res
       call mpi_abort(mpi_comm_world, ierr)
    endif
 
    noahmp_state%stc_inc = 0.0  !0 if no inc exists. TODO: need to do liau type "no update on 0"?
    noahmp_state%slc_inc = 0.0
    do nl=1, lsoil_incr
        write(layerch, '(i1.1)') nl
        if (upd_stc) then
            var_name = 'soilt'//layerch
            call read_nc_var2D(ncid, trim(incr_file), len_land_vec, res, tile2vector, &
                            6, var_name, noahmp_state%stc_inc(:, nl), tile_num)
        endif
        if (upd_slc) then
            var_name = 'soill'//layerch
            call read_nc_var2D(ncid, trim(incr_file), len_land_vec, res, tile2vector, &
                            6, var_name, noahmp_state%slc_inc(:, nl), tile_num)
        endif
    enddo
    ierr=nf90_close(ncid)
    call netcdf_err(ierr, 'closing file: '//trim(incr_file) )

end subroutine read_csg_increment

!--------------------------------------------------------
! Subroutine to read in a 2D variable from netcdf file, 
! and save to noahmp vector 
!--------------------------------------------------------
subroutine read_nc_var2D(ncid, file_name, len_land_vec, res, tile2vector, in3D_vdim,  & 
                         var_name, data_vec, l_dim)

    integer, intent(in)             :: ncid, len_land_vec, res 
    ! file name added for inspection. In case of failure, filename has all info about proc rank, ens member, and tile number
    character(len=*), intent(in)    :: file_name  
    character(len=*), intent(in)    :: var_name
    integer, intent(in)             :: tile2vector(len_land_vec,2)
    integer, intent(in)             :: in3D_vdim ! 0 - input is 2D, 
                                                 ! >0, gives dim of 3rd dimension
    double precision, intent(out)   :: data_vec(len_land_vec) 
    integer, intent(in), optional   :: l_dim  ! if variable has a level dimension, provide the level to read in (starting from 1). Only used when in3D_vdim > 0

    double precision :: dummy2D(res, res) 
    double precision :: dummy3D(res, res, in3D_vdim)  
    integer          :: nn, ierr, id_var, l_dim_local

    if(present(l_dim)) then 
        l_dim_local = l_dim
    else
        l_dim_local = 1
    endif

    ierr=nf90_inq_varid(ncid, trim(var_name), id_var)
    call netcdf_err(ierr, 'reading '//trim(var_name)//' id in '//trim(file_name) )
    if (in3D_vdim==0) then
        ierr=nf90_get_var(ncid, id_var, dummy2D)
        call netcdf_err(ierr, 'reading '//trim(var_name)//' data in '//trim(file_name) )
    else  ! special case for reading in 3D variable, and retaining only 
          ! level 1
        ierr=nf90_get_var(ncid, id_var, dummy3D)
        call netcdf_err(ierr, 'reading '//trim(var_name)//' data in '//trim(file_name) )
        dummy2D=dummy3D(:,:,l_dim_local) 
    endif

    do nn=1,len_land_vec 
        data_vec(nn) = dummy2D(tile2vector(nn,1), tile2vector(nn,2))
    enddo

end subroutine read_nc_var2D

!--------------------------------------------------------
! Subroutine to read in a 3D variable from netcdf file, 
! and save to noahmp vector
!--------------------------------------------------------
subroutine read_nc_var3D(ncid, file_name, len_land_vec, res, vdim,  & 
                tile2vector, var_name, data_vec)

    integer, intent(in)             :: ncid, len_land_vec, res, vdim
    character(len=*), intent(in)    :: file_name
    character(len=*), intent(in)    :: var_name
    integer, intent(in)             :: tile2vector(len_land_vec,2)
    double precision, intent(out)   :: data_vec(len_land_vec, vdim)

    double precision :: dummy3D(res, res, vdim) 
    integer          :: nn, ierr, id_var

    ierr=nf90_inq_varid(ncid, trim(var_name), id_var)
    call netcdf_err(ierr, 'reading '//trim(var_name)//' id in '//trim(file_name) )
    ierr=nf90_get_var(ncid, id_var, dummy3D)
    call netcdf_err(ierr, 'reading '//trim(var_name)//' data in '//trim(file_name) )

    do nn=1,len_land_vec 
        data_vec(nn,:) = dummy3D(tile2vector(nn,1), tile2vector(nn,2), :) 
    enddo

end subroutine read_nc_var3D

!--------------------------------------------------------------
! write updated fields tofv3_restarts  open on ncid
!--------------------------------------------------------------
 subroutine write_fv3_restart(file_name, noahmp_state, res, ncid, len_land_vec, tile2vector, lsoil) 

 implicit none 

 integer, intent(in) :: ncid, res, len_land_vec, lsoil
 character(len=*), intent(in)    :: file_name
 type(noahmp_type), intent(in) :: noahmp_state
 integer, intent(in) :: tile2vector(len_land_vec,2)
 

    call write_nc_var3D(ncid, trim(file_name), len_land_vec, res, lsoil, tile2vector, & 
                        'stc', noahmp_state%stc)

    call write_nc_var3D(ncid, trim(file_name), len_land_vec, res, lsoil, tile2vector, & 
                        'smc', noahmp_state%smc)

    call write_nc_var3D(ncid, trim(file_name), len_land_vec, res, lsoil, tile2vector, &
                        'slc', noahmp_state%slc)

 end subroutine write_fv3_restart


!--------------------------------------------------------
! Subroutine to write a 2D variable to the netcdf file 
!--------------------------------------------------------
subroutine write_nc_var2D(ncid, file_name, len_land_vec, res, tile2vector,   & 
                in3D_vdim, var_name, data_vec)

    integer, intent(in)             :: ncid, len_land_vec, res
    character(len=*), intent(in)    :: file_name
    character(len=*), intent(in)    :: var_name
    integer, intent(in)             :: tile2vector(len_land_vec,2)
    integer, intent(in)             :: in3D_vdim ! 0 - input is 2D, 
                                                 ! >0, gives dim of 3rd dimension
    double precision, intent(in)    :: data_vec(len_land_vec)

    double precision :: dummy2D(res, res) 
    double precision :: dummy3D(res, res, in3D_vdim)
    integer          :: nn, ierr, id_var

    ierr=nf90_inq_varid(ncid, trim(var_name), id_var)
    call netcdf_err(ierr, 'reading '//trim(var_name)//' id in '//trim(file_name) )
    if (in3D_vdim==0) then 
        ierr=nf90_get_var(ncid, id_var, dummy2D)
        call netcdf_err(ierr, 'reading '//trim(var_name)//' data in '//trim(file_name) )
    else  ! special case for reading in multi-level variable, and 
          ! retaining only first level.
        ierr=nf90_get_var(ncid, id_var, dummy3D)
        call netcdf_err(ierr, 'reading '//trim(var_name)//' data in '//trim(file_name) )
        dummy2D = dummy3D(:,:,1)
    endif
    
    ! sub in updated locations (retain previous fields for non-land)  
    do nn=1,len_land_vec 
        dummy2D(tile2vector(nn,1), tile2vector(nn,2)) = data_vec(nn) 
    enddo

    ! overwrite
    if (in3D_vdim==0) then 
        ierr = nf90_put_var( ncid, id_var, dummy2D)
        call netcdf_err(ierr, 'writing '//trim(var_name)//' in '//trim(file_name) )
    else 
        dummy3D(:,:,1) = dummy2D 
        ierr = nf90_put_var( ncid, id_var, dummy3D)
        call netcdf_err(ierr, 'writing '//trim(var_name)//' in '//trim(file_name) )
    endif
    call remove_checksum(ncid, id_var)
 
end subroutine write_nc_var2D

!--------------------------------------------------------
! Subroutine to write a 3D variable to the netcdf file 
!--------------------------------------------------------
subroutine write_nc_var3D(ncid, file_name, len_land_vec, res, vdim, & 
                tile2vector, var_name, data_vec)

    integer, intent(in)             :: ncid, len_land_vec, res, vdim
    character(len=*), intent(in)    :: file_name
    character(len=*), intent(in)    :: var_name
    integer, intent(in)             :: tile2vector(len_land_vec,2)
    double precision, intent(in)    :: data_vec(len_land_vec, vdim)

    double precision :: dummy3D(res, res, vdim)
    integer          :: nn, ierr, id_var

    ierr=nf90_inq_varid(ncid, trim(var_name), id_var)
    call netcdf_err(ierr, 'reading '//trim(var_name)//' id in '//trim(file_name) )
    ierr=nf90_get_var(ncid, id_var, dummy3D)
    call netcdf_err(ierr, 'reading '//trim(var_name)//' data in '//trim(file_name) )
    
    ! sub in updated locations (retain previous fields for non-land)  
    do nn=1,len_land_vec 
        dummy3D(tile2vector(nn,1), tile2vector(nn,2),:) = data_vec(nn,:)
    enddo

    ! overwrite
    ierr = nf90_put_var( ncid, id_var, dummy3D)
    call netcdf_err(ierr, 'writing '//trim(var_name)//' in '//trim(file_name) )
    call remove_checksum(ncid, id_var)
 
end subroutine write_nc_var3D

!> Remove the checksum attribute from a netcdf record.
!!
!! @param[in] ncid netcdf file id
!! @param[in] id_var netcdf variable id.
!!
!! @author George Gayno NCEP/EMC
 subroutine remove_checksum(ncid, id_var)

 implicit none

 integer, intent(in)       :: ncid, id_var

 integer                   :: error

 error=nf90_inquire_attribute(ncid, id_var, 'checksum')

 if (error == 0) then ! attribute was found

   error = nf90_redef(ncid)
   call netcdf_err(error, 'entering define mode' )

   error=nf90_del_att(ncid, id_var, 'checksum')
   call netcdf_err(error, 'deleting checksum' )

   error= nf90_enddef(ncid)
   call netcdf_err(error, 'ending define mode' )

 endif

 end subroutine remove_checksum

 end program apply_incr_noahmp_soil
