module ijedi_fv3_restart_mod

use atlas_module,               only: atlas_field, atlas_fieldset
use fckit_configuration_module, only: fckit_configuration
use fms2_io_mod,                only: FmsNetcdfDomainFile_t, close_file, get_dimension_names, &
                                      get_dimension_size, get_num_dimensions, &
                                      get_variable_dimension_names, get_variable_num_dimensions, &
                                      is_dimension_registered, open_file, read_restart, &
                                      register_axis, register_field, register_restart_field, &
                                      register_variable_attribute, unlimited, write_restart
use ijedi_fv3_geom_mod,         only: fv3_geom_nodes_to_atlas_nodes, fv3_geom_setup_domain
use ijedi_kinds_mod,            only: kind_real
use mpp_domains_mod,            only: center, domain2D, mpp_deallocate_domain, mpp_update_domains
use mpp_mod,                    only: mpp_pe, mpp_root_pe, mpp_sync

implicit none
private
external abor1_ftn

integer, parameter :: numfiles = 9
integer, parameter :: index_core = 1
integer, parameter :: index_trcr = 2
integer, parameter :: index_sfcd = 3
integer, parameter :: index_sfcw = 4
integer, parameter :: index_cplr = 5
integer, parameter :: index_spec = 6
integer, parameter :: index_phys = 7
integer, parameter :: index_orog = 8
integer, parameter :: index_cold = 9
integer, parameter :: strlen = 1024

public :: ijedi_fv3_restart_read
public :: ijedi_fv3_restart_write

type field_buffer
  character(len=256) :: field_name = ''
  character(len=256) :: io_name = ''
  integer :: file_index = 0
  logical :: copy_to_fieldset = .true.
  real(kind=kind_real), allocatable :: array(:,:,:)
end type field_buffer

contains

subroutine ijedi_fv3_restart_read(conf, geom_conf, afieldset, field_io_names, field_io_scaling)

type(fckit_configuration), intent(in)    :: conf
type(fckit_configuration), intent(in)    :: geom_conf
type(atlas_fieldset),      intent(inout) :: afieldset
type(fckit_configuration), intent(in)    :: field_io_names
type(fckit_configuration), intent(in)    :: field_io_scaling

type(domain2D) :: domain
type(FmsNetcdfDomainFile_t) :: fileobj(numfiles)
type(atlas_field) :: afield
character(len=:), allocatable :: active_fields(:), tracer_fields(:)
character(len=128) :: datapath, filenames(numfiles), filenames_conf(numfiles)
character(len=128) :: prefix
character(len=256) :: io_name
logical :: rstflag(numfiles)
logical :: ps_in_file, skip_coupler, prepend_date, has_prefix, ignore_checksum
logical :: have_delp, want_ps
integer :: n, npx, npy, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, ntiles, ngrid
integer :: layout_x, layout_y, indexrst, ifield, nbuf, ibuf
real(kind=kind_real) :: ptop
real(kind=kind_real), allocatable :: delp(:,:,:)
real(kind=kind_real), allocatable :: pswork(:,:,:)
real(kind=kind_real), pointer :: atlas_ptr(:,:)
type(field_buffer), allocatable :: buffers(:)

call load_runtime(conf, .false., datapath, filenames, filenames_conf, prefix, has_prefix, &
                  prepend_date, ps_in_file, skip_coupler, ignore_checksum)
call load_geom(geom_conf, npx, npy, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
               ntiles, ngrid, layout_x, layout_y, ptop)
call conf%get_or_die('active fields', active_fields)
call conf%get_or_die('tracer fields', tracer_fields)

call fv3_geom_setup_domain(domain, npx-1, npy-1, ntiles, (/layout_x, layout_y/), (/1, 1/), 3)

rstflag = .false.
have_delp = afieldset%has_field('air_pressure_thickness')
want_ps = afieldset%has_field('air_pressure_at_surface') .and. (.not. ps_in_file)
if (want_ps) allocate(delp(isd:ied, jsd:jed, npz))

nbuf = size(active_fields)
if (want_ps .and. (.not. have_delp)) nbuf = nbuf + 1
allocate(buffers(nbuf))

ibuf = 0
do ifield = 1, size(active_fields)
  if (trim(active_fields(ifield)) == 'air_pressure_at_surface' .and. (.not. ps_in_file)) cycle
  ibuf = ibuf + 1
  buffers(ibuf)%field_name = trim(active_fields(ifield))
  buffers(ibuf)%io_name = field_name_for_io(trim(active_fields(ifield)), field_io_names)
  afield = afieldset%field(trim(active_fields(ifield)))
  call get_io_file(trim(active_fields(ifield)), afield%levels(), tracer_fields, buffers(ibuf)%file_index)
  allocate(buffers(ibuf)%array(isd:ied, jsd:jed, max(1, afield%levels())))
  buffers(ibuf)%array = 0.0_kind_real
  call afield%final()
end do

if (want_ps .and. (.not. have_delp)) then
  ibuf = ibuf + 1
  buffers(ibuf)%field_name = 'air_pressure_thickness'
  buffers(ibuf)%io_name = field_name_for_io('air_pressure_thickness', field_io_names)
  buffers(ibuf)%copy_to_fieldset = .false.
  call get_io_file('air_pressure_thickness', npz, tracer_fields, buffers(ibuf)%file_index)
  allocate(buffers(ibuf)%array(isd:ied, jsd:jed, npz))
  buffers(ibuf)%array = 0.0_kind_real
end if

nbuf = ibuf

do n = 1, numfiles
  do ibuf = 1, nbuf
    if (buffers(ibuf)%file_index /= n) cycle
    if (.not. rstflag(n)) then
      if (.not. open_file(fileobj(n), trim(datapath)//'/'//trim(filenames(n)), 'read', domain, &
                          is_restart=.true., dont_add_res_to_filename=.true.)) then
        call abor1_ftn('ijedi_fv3_restart_read: failed to open '//trim(datapath)//'/'//trim(filenames(n)))
      end if
      rstflag(n) = .true.
    end if
    call register_restart_array(fileobj(n), trim(buffers(ibuf)%io_name), buffers(ibuf)%array)
  end do
  if (rstflag(n)) then
    call read_restart(fileobj(n), ignore_checksum=ignore_checksum)
    call close_file(fileobj(n))
  end if
end do

! Populate halo cells from neighbouring PE compute domains so that atlas
! ghost nodes receive valid data when copy_fv3_to_atlas is called below.
do ibuf = 1, nbuf
  call mpp_update_domains(buffers(ibuf)%array, domain)
end do

do ibuf = 1, nbuf
  call scale_array(buffers(ibuf)%array, trim(buffers(ibuf)%field_name), field_io_scaling)
  if (.not. buffers(ibuf)%copy_to_fieldset) then
    delp = buffers(ibuf)%array
    cycle
  end if
  afield = afieldset%field(trim(buffers(ibuf)%field_name))
  call afield%data(atlas_ptr)
  call copy_fv3_to_atlas(buffers(ibuf)%array, atlas_ptr, npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                         ntile, ntiles, ngrid)
  call afield%set_dirty(.true.)
  call afield%final()
  if (trim(buffers(ibuf)%field_name) == 'air_pressure_thickness') delp = buffers(ibuf)%array
end do

if (want_ps) then
  if (have_delp) then
    afield = afieldset%field('air_pressure_thickness')
    call afield%data(atlas_ptr)
    delp = 0.0_kind_real
    call copy_atlas_to_fv3(atlas_ptr, delp, npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                 ntile, ntiles, ngrid)
    call afield%final()
  end if
  afield = afieldset%field('air_pressure_at_surface')
  call afield%data(atlas_ptr)
  allocate(pswork(isd:ied, jsd:jed, 1))
  pswork = 0.0_kind_real
  pswork(:,:,1) = ptop + sum(delp, dim=3)
  call copy_fv3_to_atlas(pswork, atlas_ptr, npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                         ntile, ntiles, ngrid)
  deallocate(pswork)
  call afield%set_dirty(.true.)
  call afield%final()
  deallocate(delp)
end if

do ibuf = 1, nbuf
  if (allocated(buffers(ibuf)%array)) deallocate(buffers(ibuf)%array)
end do
deallocate(buffers)
call mpp_deallocate_domain(domain)

end subroutine ijedi_fv3_restart_read

subroutine ijedi_fv3_restart_write(conf, geom_conf, afieldset, field_io_names, field_io_scaling)

type(fckit_configuration), intent(in) :: conf
type(fckit_configuration), intent(in) :: geom_conf
type(atlas_fieldset),      intent(in) :: afieldset
type(fckit_configuration), intent(in) :: field_io_names
type(fckit_configuration), intent(in) :: field_io_scaling

type(domain2D) :: domain
type(FmsNetcdfDomainFile_t) :: fileobj(numfiles)
type(atlas_field) :: afield
character(len=:), allocatable :: active_fields(:), tracer_fields(:)
character(len=128) :: datapath, filenames(numfiles), filenames_conf(numfiles)
character(len=128) :: prefix
character(len=256) :: io_name
logical :: rstflag(numfiles)
logical :: ps_in_file, skip_coupler, prepend_date, has_prefix, ignore_checksum
integer :: n, npx, npy, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, ntiles, ngrid
integer :: layout_x, layout_y, indexrst, ifield, ibuf, date(6), calendar_type
real(kind=kind_real) :: ptop
real(kind=kind_real), pointer :: atlas_ptr(:,:)
type(field_buffer), allocatable :: buffers(:)

call load_runtime(conf, .true., datapath, filenames, filenames_conf, prefix, has_prefix, &
                  prepend_date, ps_in_file, skip_coupler, ignore_checksum)
call load_geom(geom_conf, npx, npy, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
               ntiles, ngrid, layout_x, layout_y, ptop)
call conf%get_or_die('active fields', active_fields)
call conf%get_or_die('tracer fields', tracer_fields)
call load_calendar(conf, calendar_type)

call fv3_geom_setup_domain(domain, npx-1, npy-1, ntiles, (/layout_x, layout_y/), (/1, 1/), 3)

rstflag = .false.

! Remove any stale output files so FMS open_file('overwrite') starts clean.
! FMS does not unconditionally clobber pre-existing domain-decomposed files.
block
  integer :: del_stat
  if (mpp_pe() == mpp_root_pe()) then
    do n = 1, numfiles
      if (trim(filenames(n)) == 'null') cycle
      open(unit=91, file=trim(datapath)//'/'//trim(filenames(n)), &
           status='old', iostat=del_stat)
      if (del_stat == 0) close(unit=91, status='delete')
    end do
  end if
  call mpp_sync()
end block

allocate(buffers(size(active_fields)))
do ifield = 1, size(active_fields)
  buffers(ifield)%field_name = trim(active_fields(ifield))
  buffers(ifield)%io_name = field_name_for_io(trim(active_fields(ifield)), field_io_names)
  afield = afieldset%field(trim(active_fields(ifield)))
  call afield%data(atlas_ptr)

  allocate(buffers(ifield)%array(isd:ied, jsd:jed, max(1, afield%levels())))
  buffers(ifield)%array = 0.0_kind_real
  call copy_atlas_to_fv3(atlas_ptr, buffers(ifield)%array, npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                         ntile, ntiles, ngrid)
  call unscale_array(buffers(ifield)%array, trim(active_fields(ifield)), field_io_scaling)
  call get_io_file(trim(active_fields(ifield)), afield%levels(), tracer_fields, buffers(ifield)%file_index)
  call afield%final()
end do

do n = 1, numfiles
  do ibuf = 1, size(buffers)
    if (buffers(ibuf)%file_index /= n) cycle
    if (.not. rstflag(n)) then
      if (.not. open_file(fileobj(n), trim(datapath)//'/'//trim(filenames(n)), 'overwrite', domain, &
                          is_restart=.true., dont_add_res_to_filename=.true.)) then
        call abor1_ftn('ijedi_fv3_restart_write: failed to open '//trim(datapath)//'/'//trim(filenames(n)))
      end if
      if (.not. is_dimension_registered(fileobj(n), 'Time')) then
        call register_axis(fileobj(n), 'Time', unlimited)
        call register_field(fileobj(n), 'Time', 'double', (/'Time'/))
        call register_variable_attribute(fileobj(n), 'Time', 'long_name', 'Time', str_len=4)
        call register_variable_attribute(fileobj(n), 'Time', 'units', 'time level', str_len=10)
        call register_variable_attribute(fileobj(n), 'Time', 'cartesian_axis', 'T', str_len=1)
      end if
      rstflag(n) = .true.
    end if
    call register_restart_array(fileobj(n), trim(buffers(ibuf)%io_name), buffers(ibuf)%array)
  end do
  if (rstflag(n)) then
    call write_restart(fileobj(n))
    call close_file(fileobj(n))
  end if
end do

do ibuf = 1, size(buffers)
  if (allocated(buffers(ibuf)%array)) deallocate(buffers(ibuf)%array)
end do
deallocate(buffers)

date = 0

if (mpp_pe() == mpp_root_pe() .and. (.not. skip_coupler)) then
  open(101, file=trim(adjustl(datapath))//'/'//trim(adjustl(filenames(index_cplr))), form='formatted')
  write(101, '(i6,8x,a)') calendar_type, &
    '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
  write(101, '(6i6,8x,a)') date, 'Model start time:   year, month, day, hour, minute, second'
  write(101, '(6i6,8x,a)') date, 'Current model time: year, month, day, hour, minute, second'
  close(101)
end if

call mpp_deallocate_domain(domain)

end subroutine ijedi_fv3_restart_write

subroutine load_runtime(conf, is_write, datapath, filenames, filenames_conf, prefix, has_prefix, &
                        prepend_date, ps_in_file, skip_coupler, ignore_checksum)

type(fckit_configuration), intent(in)  :: conf
logical,                   intent(in)  :: is_write
character(len=128),        intent(out) :: datapath
character(len=128),        intent(out) :: filenames(numfiles)
character(len=128),        intent(out) :: filenames_conf(numfiles)
character(len=128),        intent(out) :: prefix
logical,                   intent(out) :: has_prefix, prepend_date, ps_in_file, skip_coupler, ignore_checksum

character(len=:), allocatable :: str
logical :: filename_is_templated
integer :: n

call conf%get_or_die('datapath', str)
datapath = trim(str)
deallocate(str)

filenames_conf(index_core) = 'fv_core.res.nc'
filenames_conf(index_trcr) = 'fv_tracer.res.nc'
filenames_conf(index_sfcd) = 'sfc_data.nc'
filenames_conf(index_sfcw) = 'fv_srf_wnd.res.nc'
filenames_conf(index_cplr) = 'coupler.res'
filenames_conf(index_spec) = 'null'
filenames_conf(index_phys) = 'phy_data.nc'
filenames_conf(index_orog) = 'oro_data.nc'
filenames_conf(index_cold) = 'gfs_data.nc'

call load_filename(conf, 'filename_core', filenames_conf(index_core))
call load_filename(conf, 'filename_trcr', filenames_conf(index_trcr))
call load_filename(conf, 'filename_sfcd', filenames_conf(index_sfcd))
call load_filename(conf, 'filename_sfcw', filenames_conf(index_sfcw))
call load_filename(conf, 'filename_cplr', filenames_conf(index_cplr))
call load_filename(conf, 'filename_spec', filenames_conf(index_spec))
call load_filename(conf, 'filename_phys', filenames_conf(index_phys))
call load_filename(conf, 'filename_orog', filenames_conf(index_orog))
call load_filename(conf, 'filename_cold', filenames_conf(index_cold))

filename_is_templated = .false.
if (conf%has('filename is datetime templated')) call conf%get_or_die('filename is datetime templated', filename_is_templated)
prepend_date = .true.
if (conf%has('prepend files with date')) call conf%get_or_die('prepend files with date', prepend_date)
ps_in_file = .false.
if (conf%has('psinfile')) call conf%get_or_die('psinfile', ps_in_file)
skip_coupler = .false.
if (conf%has('skip coupler file')) call conf%get_or_die('skip coupler file', skip_coupler)
ignore_checksum = .true.
if (conf%has('ignore checksum')) call conf%get_or_die('ignore checksum', ignore_checksum)
has_prefix = conf%has('prefix')
if (has_prefix) then
  call conf%get_or_die('prefix', str)
  prefix = trim(str)
  deallocate(str)
else
  prefix = ''
end if

do n = 1, numfiles
  filenames(n) = trim(filenames_conf(n))
  if (filename_is_templated) call abor1_ftn('ijedi_fv3_restart_mod: datetime templated filenames require DateTime input')
  if (is_write .and. prepend_date) filenames(n) = trim(filenames(n))
  if (has_prefix) filenames(n) = trim(prefix)//'.'//trim(filenames_conf(n))
end do

end subroutine load_runtime

subroutine load_geom(conf, npx, npy, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, ntiles, &
                     ngrid, layout_x, layout_y, ptop)

type(fckit_configuration), intent(in) :: conf
integer, intent(out) :: npx, npy, npz, isc, iec, jsc, jec, isd, ied, jsd, jed
integer, intent(out) :: ntile, ntiles, ngrid, layout_x, layout_y
real(kind=kind_real), intent(out) :: ptop

call conf%get_or_die('npx', npx)
call conf%get_or_die('npy', npy)
call conf%get_or_die('nLevels', npz)
call conf%get_or_die('isc', isc)
call conf%get_or_die('iec', iec)
call conf%get_or_die('jsc', jsc)
call conf%get_or_die('jec', jec)
call conf%get_or_die('isd', isd)
call conf%get_or_die('ied', ied)
call conf%get_or_die('jsd', jsd)
call conf%get_or_die('jed', jed)
call conf%get_or_die('ntile', ntile)
call conf%get_or_die('ntiles', ntiles)
call conf%get_or_die('ngrid', ngrid)
call conf%get_or_die('layout_x', layout_x)
call conf%get_or_die('layout_y', layout_y)
call conf%get_or_die('air_pressure_at_top_of_atmosphere_model', ptop)

end subroutine load_geom

subroutine load_calendar(conf, calendar_type)

type(fckit_configuration), intent(in) :: conf
integer, intent(out) :: calendar_type

calendar_type = 2
if (conf%has('calendar type')) call conf%get_or_die('calendar type', calendar_type)

end subroutine load_calendar

subroutine load_filename(conf, key, value)

type(fckit_configuration), intent(in) :: conf
character(len=*),          intent(in) :: key
character(len=*),          intent(inout) :: value

character(len=:), allocatable :: str

if (conf%has(trim(key))) then
  call conf%get_or_die(trim(key), str)
  value = trim(str)
  deallocate(str)
end if

end subroutine load_filename

function field_name_for_io(field_name, field_io_names) result(io_name)

character(len=*),          intent(in) :: field_name
type(fckit_configuration), intent(in) :: field_io_names
character(len=256)                    :: io_name
character(len=:), allocatable         :: str

io_name = trim(field_name)
if (field_io_names%has(trim(field_name))) then
  call field_io_names%get_or_die(trim(field_name), str)
  io_name = trim(str)
  if (allocated(str)) deallocate(str)
end if

end function field_name_for_io

subroutine scale_array(array, field_name, field_io_scaling)

real(kind=kind_real),      intent(inout) :: array(:,:,:)
character(len=*),          intent(in)    :: field_name
type(fckit_configuration), intent(in)    :: field_io_scaling
real(kind=kind_real) :: scale

scale = 1.0_kind_real
if (field_io_scaling%has(trim(field_name))) call field_io_scaling%get_or_die(trim(field_name), scale)
array = scale * array

end subroutine scale_array

subroutine unscale_array(array, field_name, field_io_scaling)

real(kind=kind_real),      intent(inout) :: array(:,:,:)
character(len=*),          intent(in)    :: field_name
type(fckit_configuration), intent(in)    :: field_io_scaling
real(kind=kind_real) :: scale

scale = 1.0_kind_real
if (field_io_scaling%has(trim(field_name))) call field_io_scaling%get_or_die(trim(field_name), scale)
if (scale /= 0.0_kind_real) array = array / scale

end subroutine unscale_array

subroutine get_io_file(field_name, nlevels, tracer_fields, indexrst)

character(len=*),          intent(in)  :: field_name
integer,                   intent(in)  :: nlevels
character(len=:), allocatable, intent(in) :: tracer_fields(:)
integer,                   intent(out) :: indexrst

character(len=32) :: io_file

io_file = 'core'
if (is_tracer_field(field_name, tracer_fields)) io_file = 'tracer'
if (nlevels == 1) io_file = 'surface'
if (trim(field_name) == 'air_pressure_at_surface') io_file = 'surface'
if (trim(field_name) == 'geopotential_at_surface') io_file = 'core'
if (trim(field_name) == 'eastward_wind_at_surface') io_file = 'surface_wind'
if (trim(field_name) == 'northward_wind_at_surface') io_file = 'surface_wind'
if (index(trim(field_name), 'orog') /= 0) io_file = 'orography'
if (index(trim(field_name), 'fraction_of_land') /= 0) io_file = 'orography'
if (index(trim(field_name), 'cold') /= 0) io_file = 'cold'
if (trim(field_name) == 'stc') io_file = 'surface'
if (trim(field_name) == 'slc') io_file = 'surface'
if (trim(field_name) == 'soilMoistureVolumetric') io_file = 'surface'
if (trim(field_name) == 'tslb') io_file = 'surface'
if (trim(field_name) == 'smois') io_file = 'surface'
if (trim(field_name) == 'equivalent_reflectivity_factor') io_file = 'physics'

select case (trim(io_file))
case ('core')
  indexrst = index_core
case ('tracer')
  indexrst = index_trcr
case ('surface')
  indexrst = index_sfcd
case ('surface_wind')
  indexrst = index_sfcw
case ('physics')
  indexrst = index_phys
case ('orography')
  indexrst = index_orog
case ('cold')
  indexrst = index_cold
case default
  indexrst = index_core
end select

end subroutine get_io_file

logical function is_tracer_field(field_name, tracer_fields)

character(len=*), intent(in) :: field_name
character(len=:), allocatable, intent(in) :: tracer_fields(:)
integer :: n

is_tracer_field = .false.
do n = 1, size(tracer_fields)
  if (trim(tracer_fields(n)) == trim(field_name)) then
    is_tracer_field = .true.
    return
  end if
end do

end function is_tracer_field

subroutine copy_atlas_to_fv3(atlas_ptr, fv3_array, npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                             ntile, ntiles, ngrid)

real(kind=kind_real),         intent(in)    :: atlas_ptr(:,:)
real(kind=kind_real),         intent(inout) :: fv3_array(isd:, jsd:, :)
integer,                       intent(in)    :: npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed
integer,                       intent(in)    :: ntile, ntiles, ngrid
integer :: jl, a, b, ncopy
logical :: at_lower_left_corner, at_upper_left_corner, at_lower_right_corner
logical :: at_right_edge, at_upper_edge
logical :: halo_w, halo_e, halo_s, halo_n, halo_sw, halo_nw, halo_ne, halo_se, halo_nw3, halo_se6

do jl = 1, size(fv3_array, 3)
  halo_w = .true.
  halo_s = .true.
  halo_sw = .true.
  halo_e = .false.
  halo_n = .false.
  halo_nw = .false.
  halo_ne = .false.
  halo_se = .false.
  halo_nw3 = .false.
  halo_se6 = .false.

  if (ntiles == 6) then
    at_lower_left_corner = (isc == 1 .and. jsc == 1)
    at_upper_left_corner = (isc == 1 .and. jec == npy-1)
    at_lower_right_corner = (iec == npx-1 .and. jsc == 1)

    if (at_lower_left_corner) halo_sw = .false.
    if (at_upper_left_corner .and. (ntile == 3)) halo_nw3 = .true.
    if (at_lower_right_corner .and. (ntile == 6)) halo_se6 = .true.
  else if (ntiles == 1) then
    at_right_edge = (iec == npx-1)
    at_upper_edge = (jec == npy-1)

    if (at_upper_edge) then
      halo_n = .true.
      halo_nw = .true.
    end if
    if (at_right_edge) then
      halo_e = .true.
      halo_se = .true.
      if (at_upper_edge) halo_ne = .true.
    end if
  else
    call abor1_ftn('ijedi_fv3_restart_mod: copy_atlas_to_fv3 requires ntiles == 1 or 6')
  end if

  a = 1
  b = ngrid
  fv3_array(isc:iec, jsc:jec, jl) = reshape(atlas_ptr(jl, a:b), (/iec-isc+1, jec-jsc+1/))

  ncopy = jec - jsc + 1
  if (halo_w) then
    a = b + 1
    b = b + ncopy
    fv3_array(isc-1, jsc:jec, jl) = atlas_ptr(jl, a:b)
  end if
  if (halo_e) then
    a = b + 1
    b = b + ncopy
    fv3_array(iec+1, jsc:jec, jl) = atlas_ptr(jl, a:b)
  end if

  ncopy = iec - isc + 1
  if (halo_s) then
    a = b + 1
    b = b + ncopy
    fv3_array(isc:iec, jsc-1, jl) = atlas_ptr(jl, a:b)
  end if
  if (halo_n) then
    a = b + 1
    b = b + ncopy
    fv3_array(isc:iec, jec+1, jl) = atlas_ptr(jl, a:b)
  end if

  if (halo_sw) then
    a = b + 1
    b = b + 1
    fv3_array(isc-1, jsc-1, jl) = atlas_ptr(jl, a)
  end if
  if (halo_nw) then
    a = b + 1
    b = b + 1
    fv3_array(isc-1, jec+1, jl) = atlas_ptr(jl, a)
  end if
  if (halo_ne) then
    a = b + 1
    b = b + 1
    fv3_array(iec+1, jec+1, jl) = atlas_ptr(jl, a)
  end if
  if (halo_se) then
    a = b + 1
    b = b + 1
    fv3_array(iec+1, jsc-1, jl) = atlas_ptr(jl, a)
  end if
  if (halo_nw3) then
    a = b + 1
    b = b + 1
    fv3_array(isc, jec+1, jl) = atlas_ptr(jl, a)
  end if
  if (halo_se6) then
    a = b + 1
    b = b + 1
    fv3_array(iec+1, jsc, jl) = atlas_ptr(jl, a)
  end if
end do
end subroutine copy_atlas_to_fv3

subroutine copy_fv3_to_atlas(fv3_array, atlas_ptr, npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                             ntile, ntiles, ngrid)

real(kind=kind_real),          intent(in)    :: fv3_array(isd:, jsd:, :)
real(kind=kind_real),          intent(inout) :: atlas_ptr(:,:)
integer,                       intent(in)    :: npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed
integer,                       intent(in)    :: ntile, ntiles, ngrid
integer :: jl

atlas_ptr = 0.0_kind_real
do jl = 1, size(fv3_array, 3)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                     ntiles, ngrid, fv3_array(:,:,jl), atlas_ptr(jl,:))
end do

end subroutine copy_fv3_to_atlas

subroutine register_restart_array(fileobj, io_name, array)

type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
character(len=*),            intent(in)    :: io_name
real(kind=kind_real),        intent(in)    :: array(:,:,:)

logical :: is_registered
integer :: ndims, idim, num_zaxes, nz_dim, nz_field
character(len=8) :: xdim_name, ydim_name, zdim_name
character(len=8), allocatable :: dim_names(:)

if (fileobj%is_readonly) then
  ndims = get_variable_num_dimensions(fileobj, trim(io_name))
  allocate(dim_names(ndims))
  call get_variable_dimension_names(fileobj, trim(io_name), dim_names)

  xdim_name = ''
  ydim_name = ''
  do idim = 1, ndims
    if (dim_names(idim)(1:6) == 'xaxis_') xdim_name = trim(dim_names(idim))
    if (dim_names(idim)(1:6) == 'yaxis_') ydim_name = trim(dim_names(idim))
  end do
  if (len_trim(xdim_name) == 0) xdim_name = trim(dim_names(1))
  if (len_trim(ydim_name) == 0 .and. ndims >= 2) ydim_name = trim(dim_names(2))

  if (.not. is_dimension_registered(fileobj, trim(xdim_name))) then
    call register_axis(fileobj, trim(xdim_name), 'x', domain_position=center)
  end if
  if (len_trim(ydim_name) > 0 .and. .not. is_dimension_registered(fileobj, trim(ydim_name))) then
    call register_axis(fileobj, trim(ydim_name), 'y', domain_position=center)
  end if
  call register_restart_field(fileobj, trim(io_name), array)
  deallocate(dim_names)
else
  is_registered = .false.
  do idim = 1, fileobj%nx
    if (fileobj%xdims(idim)%pos == center) then
      is_registered = .true.
      xdim_name = trim(fileobj%xdims(idim)%varname)
      exit
    end if
  end do
  if (.not. is_registered) then
    write(xdim_name, '(A,I0)') 'xaxis_', fileobj%nx + 1
    call register_axis(fileobj, trim(xdim_name), 'x', domain_position=center)
    call register_field(fileobj, trim(xdim_name), 'double', (/trim(xdim_name)/))
    call register_variable_attribute(fileobj, trim(xdim_name), 'long_name', trim(xdim_name), str_len=len_trim(xdim_name))
    call register_variable_attribute(fileobj, trim(xdim_name), 'units', 'none', str_len=4)
    call register_variable_attribute(fileobj, trim(xdim_name), 'cartesian_axis', 'X', str_len=1)
  end if

  is_registered = .false.
  do idim = 1, fileobj%ny
    if (fileobj%ydims(idim)%pos == center) then
      is_registered = .true.
      ydim_name = trim(fileobj%ydims(idim)%varname)
      exit
    end if
  end do
  if (.not. is_registered) then
    write(ydim_name, '(A,I0)') 'yaxis_', fileobj%ny + 1
    call register_axis(fileobj, trim(ydim_name), 'y', domain_position=center)
    call register_field(fileobj, trim(ydim_name), 'double', (/trim(ydim_name)/))
    call register_variable_attribute(fileobj, trim(ydim_name), 'long_name', trim(ydim_name), str_len=len_trim(ydim_name))
    call register_variable_attribute(fileobj, trim(ydim_name), 'units', 'none', str_len=4)
    call register_variable_attribute(fileobj, trim(ydim_name), 'cartesian_axis', 'Y', str_len=1)
  end if

  nz_field = size(array, 3)
  if (nz_field > 1) then
    ndims = get_num_dimensions(fileobj)
    allocate(dim_names(ndims))
    call get_dimension_names(fileobj, dim_names)
    num_zaxes = 0
    is_registered = .false.
    do idim = 1, ndims
      if (dim_names(idim)(1:6) == 'zaxis_') then
        call get_dimension_size(fileobj, trim(dim_names(idim)), nz_dim)
        if (nz_dim == nz_field) then
          is_registered = .true.
          zdim_name = trim(dim_names(idim))
          exit
        end if
        num_zaxes = num_zaxes + 1
      end if
    end do
    if (.not. is_registered) then
      write(zdim_name, '(A,I0)') 'zaxis_', num_zaxes + 1
      call register_axis(fileobj, trim(zdim_name), nz_field)
      call register_field(fileobj, trim(zdim_name), 'double', (/trim(zdim_name)/))
      call register_variable_attribute(fileobj, trim(zdim_name), 'long_name', trim(zdim_name), str_len=len_trim(zdim_name))
      call register_variable_attribute(fileobj, trim(zdim_name), 'units', 'none', str_len=4)
      call register_variable_attribute(fileobj, trim(zdim_name), 'cartesian_axis', 'Z', str_len=1)
    end if
    deallocate(dim_names)
  end if

  if (nz_field > 1) then
    call register_restart_field(fileobj, trim(io_name), array, (/xdim_name, ydim_name, zdim_name, 'Time    '/))
  else
    call register_restart_field(fileobj, trim(io_name), array, (/xdim_name, ydim_name, 'Time    '/))
  end if
end if

end subroutine register_restart_array

end module ijedi_fv3_restart_mod
