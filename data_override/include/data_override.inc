!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

! This file contains the body of the data_override_r4 and data_override_r8
! modules. These modules are not intended to be used directly - they should be
! used through the data_override_mod API. See data_override.F90 for details.

use platform_mod, only: r4_kind, r8_kind
use yaml_parser_mod
use constants_mod, only: DEG_TO_RAD
use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE, stdout, stdlog, mpp_max
use mpp_mod, only : input_nml_file
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_new, horiz_interp_type, &
                             assignment(=)
use time_interp_external2_mod, only: time_interp_external_init, &
                                     time_interp_external, &
                                     time_interp_external_bridge, get_time_axis, &
                                     init_external_field, &
                                     get_external_field_size, &
                                     set_override_region, &
                                     reset_src_data_region, &
                                     NO_REGION, INSIDE_REGION, OUTSIDE_REGION, &
                                     get_external_fileobj
use fms_mod, only: write_version_number, lowercase, check_nml_error
use axis_utils2_mod,  only : nearest_index, axis_edges
use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, NULL_DOMAIN2D,operator(.NE.),operator(.EQ.)
use mpp_domains_mod, only : mpp_get_global_domain, mpp_get_data_domain
use mpp_domains_mod, only : domainUG, mpp_pass_SG_to_UG, mpp_get_UG_SG_domain, NULL_DOMAINUG
use time_manager_mod, only: time_type, OPERATOR(>), OPERATOR(<)
use fms2_io_mod,     only : FmsNetcdfFile_t, open_file, close_file, &
                            read_data, fms2_io_init, variable_exists, &
                            get_mosaic_tile_file, file_exists
use get_grid_version_mod, only: get_grid_version_1, get_grid_version_2
use fms_string_utils_mod, only: string

implicit none
private

! Include variable "version" to be written to log file.
#include<file_version.h>

!> Private type for holding field and grid information from a data table
!> @ingroup data_override_mod
type data_type
   character(len=3)   :: gridname
   character(len=128) :: fieldname_code  !< fieldname used in user's code (model)
   character(len=128) :: fieldname_file  !< fieldname used in the netcdf data file
   character(len=512) :: file_name       !< name of netCDF data file
   character(len=128) :: interpol_method !< interpolation method (default "bilinear")
   real(FMS_DATA_OVERRIDE_KIND_)        :: factor          !< For unit conversion, default=1, see OVERVIEW above
   real(FMS_DATA_OVERRIDE_KIND_)        :: lon_start, lon_end, lat_start, lat_end
   integer            :: region_type
   logical            :: multifile = .false.
   character(len=512) :: prev_file_name   !< name of netCDF data file for previous segment
   character(len=512) :: next_file_name   !< name of netCDF data file for next segment
   type(time_type), dimension(:), allocatable :: time_records
   type(time_type), dimension(:), allocatable :: time_prev_records
   type(time_type), dimension(:), allocatable :: time_next_records
end type data_type

!> Private type for holding various data fields for performing data overrides
!> @ingroup data_override_mod
type override_type
   character(len=3)                 :: gridname
   character(len=128)               :: fieldname
   integer                          :: t_index            !< index for time interp
   integer                          :: pt_index           !< previous index for time interp
   integer                          :: nt_index           !< next index for time interp
   type(horiz_interp_type), allocatable :: horz_interp(:) !< index for horizontal spatial interp
   integer                          :: dims(4)            !< dimensions(x,y,z,t) of the field in filename
   integer                          :: comp_domain(4)     !< istart,iend,jstart,jend for compute domain
   integer                          :: numthreads
   real(FMS_DATA_OVERRIDE_KIND_), allocatable :: lon_in(:)
   real(FMS_DATA_OVERRIDE_KIND_), allocatable :: lat_in(:)
   logical, allocatable             :: need_compute(:)
   integer                          :: numwindows
   integer                          :: window_size(2)
   integer                          :: is_src, ie_src, js_src, je_src
end type override_type

integer, parameter :: lkind = FMS_DATA_OVERRIDE_KIND_
integer, parameter :: max_table=100, max_array=100

integer            :: table_size !< actual size of data table
logical            :: module_is_initialized = .FALSE.

type(domain2D) :: ocn_domain,atm_domain,lnd_domain, ice_domain
type(domainUG) :: lnd_domainUG

real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:), target, allocatable :: lon_local_ocn, lat_local_ocn
real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:), target, allocatable :: lon_local_atm, lat_local_atm
real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:), target, allocatable :: lon_local_ice, lat_local_ice
real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:), target, allocatable :: lon_local_lnd, lat_local_lnd
real(FMS_DATA_OVERRIDE_KIND_) :: min_glo_lon_ocn, max_glo_lon_ocn
real(FMS_DATA_OVERRIDE_KIND_) :: min_glo_lon_atm, max_glo_lon_atm
real(FMS_DATA_OVERRIDE_KIND_) :: min_glo_lon_lnd, max_glo_lon_lnd
real(FMS_DATA_OVERRIDE_KIND_) :: min_glo_lon_ice, max_glo_lon_ice
integer :: num_fields = 0 !< number of fields in override_array already processed

type(data_type), dimension(:), allocatable :: data_table !< user-provided data table

type(data_type)                            :: default_table
type(override_type), dimension(max_array)  :: override_array !< to store processed fields
type(override_type)                        :: default_array
logical                                    :: debug_data_override
logical                                    :: grid_center_bug = .false.
logical                                    :: reproduce_null_char_bug = .false. !> Flag indicating
                                              !! to reproduce the mpp_io bug where lat/lon_bnd were
                                              !! not read correctly if null characters are present in
                                              !! the netcdf file
logical                                    :: use_data_table_yaml = .false.

namelist /data_override_nml/ debug_data_override, grid_center_bug, reproduce_null_char_bug, use_data_table_yaml

public :: DATA_OVERRIDE_INIT_IMPL_, DATA_OVERRIDE_UNSET_ATM_, DATA_OVERRIDE_UNSET_OCN_, &
        & DATA_OVERRIDE_UNSET_LND_, DATA_OVERRIDE_UNSET_ICE_, DATA_OVERRIDE_0D_, &
        & DATA_OVERRIDE_2D_, DATA_OVERRIDE_3D_, DATA_OVERRIDE_UG_1D_, &
        & DATA_OVERRIDE_UG_2D_

contains

!> @brief Assign default values for default_table, get domain of component models,
!! get global grids of component models.
!! Users should call data_override_init before calling data_override
!!
!! This subroutine should be called by data_override_init.
!!
!! Data_table is initialized here with default values. Users should provide "real" values
!! that will override the default values. Real values can be specified in either data_table
!! or data_table.yaml. Each line of data_table contains one data_entry. Items of data_entry
!! are comma-separated.
subroutine DATA_OVERRIDE_INIT_IMPL_(Atm_domain_in, Ocean_domain_in, Ice_domain_in, Land_domain_in, Land_domainUG_in)
  type (domain2d), intent(in), optional :: Atm_domain_in    !> Atmosphere domain
  type (domain2d), intent(in), optional :: Ocean_domain_in  !> Ocean domain
  type (domain2d), intent(in), optional :: Ice_domain_in    !> Ice domain
  type (domain2d), intent(in), optional :: Land_domain_in   !> Land domain
  type(domainUG) , intent(in), optional :: Land_domainUG_in !> Land domain, unstructured grid

  character(len=128)    :: grid_file = 'INPUT/grid_spec.nc'
  integer               :: is,ie,js,je,use_get_grid_version
  integer               :: i, iunit, io_status, ierr
  logical               :: atm_on, ocn_on, lnd_on, ice_on, lndUG_on
  logical               :: file_open
  type(FmsNetcdfFile_t) :: fileobj

  debug_data_override = .false.

  read (input_nml_file, data_override_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'data_override_nml')
  iunit = stdlog()
  write(iunit, data_override_nml)

! grid_center_bug is no longer supported.
if (grid_center_bug) then
  call mpp_error(FATAL, "data_override_init: You have overridden the default value of " // &
     "grid_center_bug and set it to .true. in data_override_nml.  This was a temporary workaround " // &
     "that is no longer supported. Please remove this namelist variable.")
endif

if (use_data_table_yaml) then
   call mpp_error(NOTE, "You are using YAML.")
else
   call mpp_error(NOTE, "You are using the legacy table.")
end if

  atm_on = PRESENT(Atm_domain_in)
  ocn_on = PRESENT(Ocean_domain_in)
  lnd_on = PRESENT(Land_domain_in)
  ice_on = PRESENT(Ice_domain_in)
  lndUG_on = PRESENT(Land_domainUG_in)
  if(.not. module_is_initialized) then
    atm_domain = NULL_DOMAIN2D
    ocn_domain = NULL_DOMAIN2D
    lnd_domain = NULL_DOMAIN2D
    ice_domain = NULL_DOMAIN2D
    lnd_domainUG = NULL_DOMAINUG
  end if
  if (atm_on) atm_domain = Atm_domain_in
  if (ocn_on) ocn_domain = Ocean_domain_in
  if (lnd_on) lnd_domain = Land_domain_in
  if (ice_on) ice_domain = Ice_domain_in
  if (lndUG_on) lnd_domainUG = Land_domainUG_in

  if(.not. module_is_initialized) then
    call horiz_interp_init
    call write_version_number("DATA_OVERRIDE_MOD", version)

!  Initialize user-provided data table
    default_table%gridname = 'non'
    default_table%fieldname_code = 'none'
    default_table%fieldname_file = 'none'
    default_table%file_name = 'none'
    default_table%factor = 1._lkind
    default_table%interpol_method = 'bilinear'
    default_table%multifile = .false.
    default_table%prev_file_name = ''
    default_table%next_file_name = ''

#ifdef use_yaml
    if (use_data_table_yaml) then
       if (file_exists("data_table")) &
         call mpp_error(FATAL, "You cannot have the legacy data_table if use_data_table_yaml=.true.")
       call read_table_yaml(data_table)
    else
       if (file_exists("data_table.yaml"))&
         call mpp_error(FATAL, "You cannot have the yaml data_table if use_data_table_yaml=.false.")
       allocate(data_table(max_table))
       do i = 1, max_table
          data_table(i) = default_table
       enddo
       call read_table(data_table)
    end if
#else
    if (file_exists("data_table.yaml"))&
         call mpp_error(FATAL, "You cannot have the yaml data_table if use_data_table_yaml=.false.")

    if (use_data_table_yaml) then
       call mpp_error(FATAL, "You cannot have use_data_table_yaml=.true. without compiling with -Duse_yaml")
    else

       allocate(data_table(max_table))
       do i = 1, max_table
          data_table(i) = default_table
       enddo
       call read_table(data_table)
    end if
#endif

!  Initialize override array
    default_array%gridname = 'NONE'
    default_array%fieldname = 'NONE'
    default_array%t_index = -1
    default_array%dims = -1
    default_array%comp_domain = -1
    do i = 1, max_array
       override_array(i) = default_array
    enddo
    call time_interp_external_init
 end if

 module_is_initialized = .TRUE.

 if ( .NOT. (atm_on .or. ocn_on .or. lnd_on .or. ice_on .or. lndUG_on)) return
 if (table_size .eq. 0) then
    call mpp_error(NOTE, "data_table is empty, not doing any data_overrides")
    return
 endif
 call fms2_io_init

! Test if grid_file is already opened
 inquire (file=trim(grid_file), opened=file_open)
 if(file_open) call mpp_error(FATAL, trim(grid_file)//' already opened')

 if(.not. open_file(fileobj, grid_file, 'read' )) then
   call mpp_error(FATAL, 'data_override_mod: Error in opening file '//trim(grid_file))
 endif

 if(variable_exists(fileobj, "x_T" ) .OR. variable_exists(fileobj, "geolon_t" ) ) then
   use_get_grid_version = 1
   call close_file(fileobj)
 else if(variable_exists(fileobj, "ocn_mosaic_file" ) .OR. variable_exists(fileobj, "gridfiles" ) ) then
   use_get_grid_version = 2
   if(variable_exists(fileobj, "gridfiles" ) ) then
     if(count_ne_1((ocn_on .OR. ice_on), lnd_on, atm_on)) call mpp_error(FATAL, 'data_override_mod: the grid file ' //&
          'is a solo mosaic, one and only one of atm_on, lnd_on or ice_on/ocn_on should be true')
   end if
 else
   call mpp_error(FATAL, 'data_override_mod: none of x_T, geolon_t, ocn_mosaic_file or gridfiles exist in '// &
                  & trim(grid_file))
 endif

 if(use_get_grid_version .EQ. 1) then
    if (atm_on .and. .not. allocated(lon_local_atm) ) then
       call mpp_get_compute_domain( atm_domain,is,ie,js,je)
       allocate(lon_local_atm(is:ie,js:je), lat_local_atm(is:ie,js:je))
       call get_grid_version_1(grid_file, 'atm', atm_domain, is, ie, js, je, lon_local_atm, lat_local_atm, &
          min_glo_lon_atm, max_glo_lon_atm )
    endif
    if (ocn_on .and. .not. allocated(lon_local_ocn) ) then
       call mpp_get_compute_domain( ocn_domain,is,ie,js,je)
       allocate(lon_local_ocn(is:ie,js:je), lat_local_ocn(is:ie,js:je))
       call get_grid_version_1(grid_file, 'ocn', ocn_domain, is, ie, js, je, lon_local_ocn, lat_local_ocn, &
          min_glo_lon_ocn, max_glo_lon_ocn )
    endif

    if (lnd_on .and. .not. allocated(lon_local_lnd) ) then
       call mpp_get_compute_domain( lnd_domain,is,ie,js,je)
       allocate(lon_local_lnd(is:ie,js:je), lat_local_lnd(is:ie,js:je))
       call get_grid_version_1(grid_file, 'lnd', lnd_domain, is, ie, js, je, lon_local_lnd, lat_local_lnd, &
          min_glo_lon_lnd, max_glo_lon_lnd )
    endif

    if (ice_on .and. .not. allocated(lon_local_ice) ) then
       call mpp_get_compute_domain( ice_domain,is,ie,js,je)
       allocate(lon_local_ice(is:ie,js:je), lat_local_ice(is:ie,js:je))
       call get_grid_version_1(grid_file, 'ice', ice_domain, is, ie, js, je, lon_local_ice, lat_local_ice, &
          min_glo_lon_ice, max_glo_lon_ice )
    endif
 else
   if (atm_on .and. .not. allocated(lon_local_atm) ) then
       call mpp_get_compute_domain(atm_domain,is,ie,js,je)
       allocate(lon_local_atm(is:ie,js:je), lat_local_atm(is:ie,js:je))
       call get_grid_version_2(fileobj, 'atm', atm_domain, is, ie, js, je, lon_local_atm, lat_local_atm, &
                               min_glo_lon_atm, max_glo_lon_atm )
   endif

   if (ocn_on .and. .not. allocated(lon_local_ocn) ) then
       call mpp_get_compute_domain( ocn_domain,is,ie,js,je)
       allocate(lon_local_ocn(is:ie,js:je), lat_local_ocn(is:ie,js:je))
       call get_grid_version_2(fileobj, 'ocn', ocn_domain, is, ie, js, je, lon_local_ocn, lat_local_ocn, &
                               min_glo_lon_ocn, max_glo_lon_ocn )
   endif

   if (lnd_on .and. .not. allocated(lon_local_lnd) ) then
       call mpp_get_compute_domain( lnd_domain,is,ie,js,je)
       allocate(lon_local_lnd(is:ie,js:je), lat_local_lnd(is:ie,js:je))
       call get_grid_version_2(fileobj, 'lnd', lnd_domain, is, ie, js, je, lon_local_lnd, lat_local_lnd, &
                               min_glo_lon_lnd, max_glo_lon_lnd )
   endif

   if (ice_on .and. .not. allocated(lon_local_ice) ) then
       call mpp_get_compute_domain( ice_domain,is,ie,js,je)
       allocate(lon_local_ice(is:ie,js:je), lat_local_ice(is:ie,js:je))
       call get_grid_version_2(fileobj, 'ocn', ice_domain, is, ie, js, je, lon_local_ice, lat_local_ice, &
                               min_glo_lon_ice, max_glo_lon_ice )
   endif
 end if
 if(use_get_grid_version .EQ. 2) then
   call close_file(fileobj)
 end if
end subroutine DATA_OVERRIDE_INIT_IMPL_

!> @brief Implementation of the following truth table:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Arg 1    Arg 2    Arg 3   | Result  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! .true.   .true.   .true.  | .true.  !!
!! .true.   .true.   .false. | .true.  !!
!! .true.   .false.  .true.  | .true.  !!
!! .true.   .false.  .false. | .false. !!
!! .false.  .true.   .true.  | .true.  !!
!! .false.  .true.   .false. | .false. !!
!! .false.  .false.  .true.  | .false. !!
!! .false.  .false.  .false. | .true.  !!
function count_ne_1(in_1, in_2, in_3)
  logical, intent(in)  :: in_1 !< Argument 1
  logical, intent(in)  :: in_2 !< Argument 2
  logical, intent(in)  :: in_3 !< Argument 3
  logical :: count_ne_1

  count_ne_1 = .not.(in_1.NEQV.in_2.NEQV.in_3) .OR. (in_1.AND.in_2.AND.in_3)
end function count_ne_1

subroutine read_table(data_table)
    type(data_type), dimension(max_table), intent(inout) :: data_table

    integer :: ntable
    integer :: ntable_lima
    integer :: ntable_new

    integer :: iunit
    integer :: io_status
    integer :: index_1col, index_2col
    character(len=256)    :: record
    type(data_type)  :: data_entry

    logical               :: ongrid
    logical               :: table_exists !< Flag indicating existence of data_table
    character(len=128)    :: region, region_type

    integer :: sunit

!  Read coupler_table
    inquire(file='data_table', EXIST=table_exists)
    if (.not. table_exists) then
      call mpp_error(NOTE, 'data_override_mod: File data_table does not exist.')
      table_size = 0
      return
    end if

    open(newunit=iunit, file='data_table', action='READ', iostat=io_status)
    if(io_status/=0) call mpp_error(FATAL, 'data_override_mod: Error in opening file data_table.')

    ntable = 0
    ntable_lima = 0
    ntable_new = 0

    do while (ntable <= max_table)
       read(iunit,'(a)',end=100) record
       if (record(1:1) == '#') cycle
       if (record(1:10) == '          ') cycle
       ntable=ntable+1
       if(index(lowercase(record), "inside_region") .ne. 0 .or. index(lowercase(record), "outside_region") .ne. 0) then
          if(index(lowercase(record), ".false.") .ne. 0 .or. index(lowercase(record), ".true.") .ne. 0 ) then
             ntable_lima = ntable_lima + 1
             read(record,*,err=99) data_entry%gridname, data_entry%fieldname_code, data_entry%fieldname_file, &
                                   data_entry%file_name, ongrid, data_entry%factor, region, region_type
             if(ongrid) then
                data_entry%interpol_method = 'none'
             else
                data_entry%interpol_method = 'bilinear'
             endif
          else
             ntable_new=ntable_new+1
             read(record,*,err=99) data_entry%gridname, data_entry%fieldname_code, data_entry%fieldname_file, &
                                   data_entry%file_name, data_entry%interpol_method, data_entry%factor, region, &
                                 & region_type

             if (index(data_entry%file_name, ":") .ne. 0) then
               data_entry%multifile = .true.
               index_1col = index(data_entry%file_name, ":")
               index_2col = index(data_entry%file_name(index_1col+1:), ":")
               if (index_2col .eq. 0) call mpp_error(FATAL, "data_override: when bridging over forcing files, " &
                 // "central forcing files must be preceded AND followed by the column (:) separator")
               data_entry%prev_file_name = data_entry%file_name(1:index_1col-1)
               data_entry%next_file_name = data_entry%file_name(index_1col+index_2col+1:)
               ! once previous/next files are filled in, overwrite current
               data_entry%file_name = data_entry%file_name(index_1col+1:index_1col+index_2col-1)
             else
               data_entry%multifile = .false.
               data_entry%prev_file_name = ""
               data_entry%next_file_name = ""
             endif
             if (data_entry%interpol_method == 'default') then
                data_entry%interpol_method = default_table%interpol_method
             endif
             if (.not.(data_entry%interpol_method == 'default'  .or. &
                  data_entry%interpol_method == 'bicubic'  .or. &
                  data_entry%interpol_method == 'bilinear' .or. &
                  data_entry%interpol_method == 'none')) then
                sunit = stdout()
                write(sunit,*)" gridname is ", trim(data_entry%gridname)
                write(sunit,*)" fieldname_code is ", trim(data_entry%fieldname_code)
                write(sunit,*)" fieldname_file is ", trim(data_entry%fieldname_file)
                write(sunit,*)" file_name is ", trim(data_entry%file_name)
                write(sunit,*)" factor is ", data_entry%factor
                write(sunit,*)" interpol_method is ", trim(data_entry%interpol_method)
                call mpp_error(FATAL, 'data_override_mod: invalid last entry in data_override_table, ' &
                     //'its value should be "default", "bicubic", "bilinear" or "none" ')
             endif
          endif
          if( trim(region_type) == "inside_region" ) then
             data_entry%region_type = INSIDE_REGION
          else if( trim(region_type) == "outside_region" ) then
             data_entry%region_type = OUTSIDE_REGION
          else
             call mpp_error(FATAL, 'data_override_mod: region type should be inside_region or outside_region')
          endif
          if (data_entry%file_name == "") call mpp_error(FATAL, &
              "data_override: filename not given in data_table when region_type is not NO_REGION")
          if(data_entry%fieldname_file == "") call mpp_error(FATAL, &
             "data_override: fieldname_file must be specified in data_table when region_type is not NO_REGION")
          if( trim(data_entry%interpol_method) == 'none') call mpp_error(FATAL, &
             "data_override(data_override_init): ongrid must be false when region_type is not NO_REGION")
          read(region,*) data_entry%lon_start, data_entry%lon_end, data_entry%lat_start, data_entry%lat_end
          !--- make sure data_entry%lon_end > data_entry%lon_start and data_entry%lat_end > data_entry%lat_start
          if(data_entry%lon_end .LE. data_entry%lon_start) call mpp_error(FATAL, &
             "data_override: lon_end should be greater than lon_start")
          if(data_entry%lat_end .LE. data_entry%lat_start) call mpp_error(FATAL, &
             "data_override: lat_end should be greater than lat_start")
       ! old format
       else if (index(lowercase(record), ".false.") .ne. 0 .or. index(lowercase(record), ".true.") .ne. 0 ) then
          ntable_lima = ntable_lima + 1
          read(record,*,err=99) data_entry%gridname, data_entry%fieldname_code, data_entry%fieldname_file, &
                                   data_entry%file_name, ongrid, data_entry%factor
          if (index(data_entry%file_name, ":") .ne. 0) then
            data_entry%multifile = .true.
            index_1col = index(data_entry%file_name, ":")
            index_2col = index(data_entry%file_name(index_1col+1:), ":")
            if (index_2col .eq. 0) call mpp_error(FATAL, "data_override: when bridging over forcing files, " &
              // "central forcing files must be preceded AND followed by the column (:) separator")
            data_entry%prev_file_name = data_entry%file_name(1:index_1col-1)
            data_entry%next_file_name = data_entry%file_name(index_1col+index_2col+1:)
            ! once previous/next files are filled in, overwrite current
            data_entry%file_name = data_entry%file_name(index_1col+1:index_1col+index_2col-1)
          else
            data_entry%multifile = .false.
            data_entry%prev_file_name = ""
            data_entry%next_file_name = ""
          endif
          if(ongrid) then
             data_entry%interpol_method = 'none'
          else
             data_entry%interpol_method = 'bilinear'
          endif
          data_entry%lon_start = 0.0_lkind
          data_entry%lon_end   = -1.0_lkind
          data_entry%lat_start = 0.0_lkind
          data_entry%lat_end   = -1.0_lkind
          data_entry%region_type = NO_REGION
       else                                      ! new format
          ntable_new=ntable_new+1
          read(record,*,err=99) data_entry%gridname, data_entry%fieldname_code, data_entry%fieldname_file, &
                                data_entry%file_name, data_entry%interpol_method, data_entry%factor
          if (index(data_entry%file_name, ":") .ne. 0) then
            index_1col = index(data_entry%file_name, ":")
            index_2col = index(data_entry%file_name(index_1col+1:), ":")
            data_entry%multifile = .true.
            if (index_2col .eq. 0) call mpp_error(FATAL, "data_override: when bridging over forcing files, " &
              // "central forcing files must be preceded AND followed by the column (:) separator")
            data_entry%prev_file_name = data_entry%file_name(1:index_1col-1)
            data_entry%next_file_name = data_entry%file_name(index_1col+index_2col+1:)
            ! once previous/next files are filled in, overwrite current
            data_entry%file_name = data_entry%file_name(index_1col+1:index_1col+index_2col-1)
          else
            data_entry%multifile = .false.
            data_entry%prev_file_name = ""
            data_entry%next_file_name = ""
          endif
          if (data_entry%interpol_method == 'default') then
            data_entry%interpol_method = default_table%interpol_method
          endif
          if (.not.(data_entry%interpol_method == 'default'  .or. &
                    data_entry%interpol_method == 'bicubic'  .or. &
                    data_entry%interpol_method == 'bilinear' .or. &
                    data_entry%interpol_method == 'none')) then
             sunit = stdout()
             write(sunit,*)" gridname is ", trim(data_entry%gridname)
             write(sunit,*)" fieldname_code is ", trim(data_entry%fieldname_code)
             write(sunit,*)" fieldname_file is ", trim(data_entry%fieldname_file)
             write(sunit,*)" file_name is ", trim(data_entry%file_name)
             write(sunit,*)" factor is ", data_entry%factor
             write(sunit,*)" interpol_method is ", trim(data_entry%interpol_method)
             call mpp_error(FATAL, 'data_override_mod: invalid last entry in data_override_table, ' &
                               //'its value should be "default", "bicubic", "bilinear" or "none" ')
          endif
          data_entry%lon_start = 0.0_lkind
          data_entry%lon_end   = -1.0_lkind
          data_entry%lat_start = 0.0_lkind
          data_entry%lat_end   = -1.0_lkind
          data_entry%region_type = NO_REGION
       endif
       data_table(ntable) = data_entry
    enddo
    call mpp_error(FATAL,'too many enries in data_table')
99  call mpp_error(FATAL,'error in data_table format')
100 continue
    table_size = ntable
    if(ntable_new*ntable_lima /= 0) call mpp_error(FATAL, &
       'data_override_mod: New and old formats together in same data_table not supported')
    close(iunit, iostat=io_status)
    if(io_status/=0) call mpp_error(FATAL, 'data_override_mod: Error in closing file data_table')
end subroutine read_table

#ifdef use_yaml
!> @brief Read and parse the data_table.yaml
subroutine read_table_yaml(data_table)
    type(data_type), dimension(:), allocatable, intent(out) :: data_table !< Contents of the data_table.yaml

    integer, allocatable :: entry_id(:)
    integer :: nentries
    integer :: i
    character(len=50) :: buffer
    integer :: file_id

    file_id = open_and_parse_file("data_table.yaml")
    if (file_id==999) then
      nentries = 0
    else
      nentries = get_num_blocks(file_id, "data_table")
      allocate(data_table(nentries))
      allocate(entry_id(nentries))
      call get_block_ids(file_id, "data_table", entry_id)

      do i = 1, nentries
         call get_value_from_key(file_id, entry_id(i), "gridname", data_table(i)%gridname)
         call check_for_valid_gridname(data_table(i)%gridname)
         call get_value_from_key(file_id, entry_id(i), "fieldname_code", data_table(i)%fieldname_code)

         data_table(i)%fieldname_file = ""
         call get_value_from_key(file_id, entry_id(i), "fieldname_file", data_table(i)%fieldname_file, &
           & is_optional=.true.)

         data_table(i)%multifile = .false.
         call get_value_from_key(file_id, entry_id(i), "is_multi_file", data_table(i)%multifile, &
           & is_optional=.true.)

         if (data_table(i)%multifile) then
          data_table(i)%prev_file_name = ""
          data_table(i)%next_file_name = ""
          call get_value_from_key(file_id, entry_id(i), "prev_file_name", data_table(i)%prev_file_name, &
            & is_optional=.true.)
          call get_value_from_key(file_id, entry_id(i), "next_file_name", data_table(i)%next_file_name, &
            & is_optional=.true.)
          if (trim(data_table(i)%prev_file_name) .eq. "" .and. trim(data_table(i)%next_file_name) .eq. "") &
            call mpp_error(FATAL, "The prev_file_name and next_file_name must be present if is_multi_file. "//&
              "Check your data_table.yaml entry for field:"//trim(data_table(i)%gridname)//":"//&
              trim(data_table(i)%fieldname_code))
         endif

         data_table(i)%file_name = ""
         call get_value_from_key(file_id, entry_id(i), "file_name", data_table(i)%file_name, &
           & is_optional=.true.)

         data_table(i)%interpol_method = "none"
         call get_value_from_key(file_id, entry_id(i), "interpol_method", data_table(i)%interpol_method, &
           & is_optional=.true.)
         call check_interpol_method(data_table(i)%interpol_method, data_table(i)%file_name, &
            data_table(i)%fieldname_file)

         call get_value_from_key(file_id, entry_id(i), "factor", data_table(i)%factor)
         buffer = ""
         call get_value_from_key(file_id, entry_id(i), "region_type", buffer, is_optional=.true.)
         call check_and_set_region_type(buffer, data_table(i)%region_type)

        if (data_table(i)%region_type .ne. NO_REGION) then
          call get_value_from_key(file_id, entry_id(i), "lon_start", data_table(i)%lon_start, is_optional=.true.)
          call get_value_from_key(file_id, entry_id(i), "lon_end", data_table(i)%lon_end, is_optional=.true.)
          call get_value_from_key(file_id, entry_id(i), "lat_start", data_table(i)%lat_start, is_optional=.true.)
          call get_value_from_key(file_id, entry_id(i), "lat_end", data_table(i)%lat_end, is_optional=.true.)
          call check_valid_lat_lon(data_table(i)%lon_start, data_table(i)%lon_end, &
            data_table(i)%lat_start, data_table(i)%lat_end)
        endif
      end do

    end if
    table_size = nentries !< Because one variable is not enough
end subroutine read_table_yaml

!> @brief Check if a grid name is valid, crashes if it is not
subroutine check_for_valid_gridname(gridname)
   character(len=*), intent(in) :: gridname !< Gridname

   select case(trim(gridname))
   case ("OCN", "ATM", "LND", "ICE")
   case default
      call mpp_error(FATAL, trim(gridname)//" is not a valid gridname. "//&
         "The acceptable values are OCN ATM LND and ICE. Check your data_table.yaml")
   end select
end subroutine check_for_valid_gridname

!> @brief Check if the interpol method is correct, crashes if it is not
subroutine check_interpol_method(interp_method, filename, fieldname)
   character(len=*), intent(in) :: interp_method !< The interpo_method
   character(len=*), intent(in) :: filename      !< The filename
   character(len=*), intent(in) :: fieldname     !< The fieldname in the file

   select case(trim(interp_method))
   case ("bicubic", "bilinear")
      if (trim(filename) .eq. "" .or. trim(fieldname) .eq. "") call mpp_error(FATAL, &
         "The file_name and the fieldname_file must be set if using the bicubic or bilinear interpolation method."//&
         " Check your data_table.yaml")
   case ("none")
      if (trim(filename) .ne. "" ) then
         if (trim(fieldname) .eq. "") call mpp_error(FATAL, &
            "If the interpol_method is none and file_name is specified (ongrid case), "//&
            "you must also specify the fieldname_file")
      endif
   case default
      call mpp_error(FATAL, trim(interp_method)//" is not a valid interp method. "//&
         "The acceptable values are bilinear and bicubic")
   end select
end subroutine check_interpol_method

!> @brief Check if a region_type is valid, crashes if it is not. Otherwise it sets the
!! correct integer parameter.
subroutine check_and_set_region_type(region_type_str, region_type_int)
   character(len=*), intent(in)  :: region_type_str !< The region type as defined in the data.yaml
   integer,          intent(out) :: region_type_int !< The region type as an integer parameter

   select case(trim(region_type_str))
   case ("inside_region")
      region_type_int = INSIDE_REGION
   case ("outside_region")
      region_type_int = OUTSIDE_REGION
   case ("")
      region_type_int = NO_REGION
   case default
      call mpp_error(FATAL, trim(region_type_str)//" is not a valid region type. "//&
         "The acceptable values are inside_region and outside_regioon. Check your data_table.yaml")
   end select
end subroutine check_and_set_region_type

!> @brief Check if a region lon_start, lon_end, lat_start and lat_end is valid.
!! Crashes if it is not.
subroutine check_valid_lat_lon(lon_start, lon_end, lat_start, lat_end)
   real(FMS_DATA_OVERRIDE_KIND_), intent(in) :: lon_start !< Starting longitude of the data_override region
   real(FMS_DATA_OVERRIDE_KIND_), intent(in) :: lon_end   !< Ending longitude of the data_override region
   real(FMS_DATA_OVERRIDE_KIND_), intent(in) :: lat_start !< Starting lattiude of the data_override region
   real(FMS_DATA_OVERRIDE_KIND_), intent(in) :: lat_end   !< Ending lattiude of the data_override region

   if (lon_start > lon_end) call mpp_error(FATAL, &
      "lon_start:"//string(lon_start)//" is greater than lon_end"//string(lon_end)//&
      ". Check your data_table.yaml.")

   if (lat_start > lat_end) call mpp_error(FATAL, &
      "lat_start:"//string(lat_start)//" is greater than lat_end:"//string(lat_end)//&
      ". Check your data_table.yaml.")
end subroutine check_valid_lat_lon
#endif

subroutine DATA_OVERRIDE_UNSET_ATM_
  atm_domain = NULL_DOMAIN2D
  if (allocated(lon_local_atm)) deallocate(lon_local_atm)
  if (allocated(lat_local_atm)) deallocate(lat_local_atm)
end subroutine

subroutine DATA_OVERRIDE_UNSET_OCN_
  ocn_domain = NULL_DOMAIN2D
  if (allocated(lon_local_ocn)) deallocate(lon_local_ocn)
  if (allocated(lat_local_ocn)) deallocate(lat_local_ocn)
end subroutine

subroutine DATA_OVERRIDE_UNSET_LND_
  lnd_domain = NULL_DOMAIN2D
  if (allocated(lon_local_lnd)) deallocate(lon_local_lnd)
  if (allocated(lat_local_lnd)) deallocate(lat_local_lnd)
end subroutine

subroutine DATA_OVERRIDE_UNSET_ICE_
  ice_domain = NULL_DOMAIN2D
  if (allocated(lon_local_ice)) deallocate(lon_local_ice)
  if (allocated(lat_local_ice)) deallocate(lat_local_ice)
end subroutine

!> @brief Given a gridname, this routine returns the working domain associated with this gridname
subroutine get_domain(gridname, domain, comp_domain)
  character(len=3), intent(in) :: gridname
  type(domain2D), intent(inout) :: domain
  integer, intent(out), optional :: comp_domain(4) !< istart,iend,jstart,jend for compute domain

  domain = NULL_DOMAIN2D
  select case (gridname)
     case('OCN')
        domain = ocn_domain
     case('ATM')
        domain = atm_domain
     case('LND')
        domain = lnd_domain
     case('ICE')
        domain = ice_domain
     case default
        call mpp_error(FATAL,'error in data_override get_domain')
  end select
  if(domain .EQ. NULL_DOMAIN2D) call mpp_error(FATAL,'data_override: failure in get_domain')
  if(present(comp_domain)) &
     call mpp_get_compute_domain(domain,comp_domain(1),comp_domain(2),comp_domain(3),comp_domain(4))
end subroutine get_domain

!> @brief Given a gridname, this routine returns the working domain associated with this gridname
subroutine get_domainUG(gridname, UGdomain, comp_domain)
  character(len=3), intent(in) :: gridname
  type(domainUG), intent(inout) :: UGdomain
  integer, intent(out), optional :: comp_domain(4) !< istart,iend,jstart,jend for compute domain
  type(domain2D), pointer :: SGdomain => NULL()

  UGdomain = NULL_DOMAINUG
  select case (gridname)
     case('LND')
        UGdomain = lnd_domainUG
     case default
        call mpp_error(FATAL,'error in data_override get_domain')
  end select
!  if(UGdomain .EQ. NULL_DOMAINUG) call mpp_error(FATAL,'data_override: failure in get_domain')
  if(present(comp_domain)) &
     call mpp_get_UG_SG_domain(UGdomain,SGdomain)
     call mpp_get_compute_domain(SGdomain,comp_domain(1),comp_domain(2),comp_domain(3),comp_domain(4))
end subroutine get_domainUG
!===============================================================================================

!> @brief Routine to perform data override for scalar fields
subroutine DATA_OVERRIDE_0D_(gridname,fieldname_code,data_out,time,override,data_index)
  character(len=3), intent(in) :: gridname !< model grid ID (ocn,ice,atm,lnd)
  character(len=*), intent(in) :: fieldname_code !< field name as used in the model (may be
                                                 !! different from the name in NetCDF data file)
  logical, intent(out), optional :: override !< true if the field has been overriden succesfully
  type(time_type), intent(in) :: time !< (target) model time
  real(FMS_DATA_OVERRIDE_KIND_), intent(out) :: data_out !< output data array returned by this call
  integer, intent(in), optional :: data_index

  type(time_type)    :: first_record !< first record of "current" file
  type(time_type)    :: last_record !< last record of "current" file
  character(len=512) :: filename !< file containing source data
  character(len=512) :: prevfilename !< file containing previous source data, when using multiple files
  character(len=512) :: nextfilename !< file containing next source data, when using multiple files
  character(len=128) :: fieldname !< fieldname used in the data file
  integer :: index1 !< field index in data_table
  integer            :: dims(4)
  integer            :: prev_dims(4) !< dimensions of previous source data, when using multiple files
  integer            :: next_dims(4) !< dimensions of next source data, when using multiple files
  integer :: id_time !< index for time interp in override array
  integer            :: id_time_prev=-1 !< time index for previous file, when using multiple files
  integer            :: id_time_next=-1 !< time index for next file, when using multiple files
  integer :: curr_position !< position of the field currently processed in override_array
  integer :: i
  real(FMS_DATA_OVERRIDE_KIND_) :: factor
  logical :: multifile !< use multiple consecutive files for override

  if(.not.module_is_initialized) &
       call mpp_error(FATAL,'Error: need to call data_override_init first')

!1  Look  for the data file in data_table
  if(PRESENT(override)) override = .false.
  if (present(data_index)) then
    index1 = data_index
  else
    index1 = -1
    do i = 1, table_size
       if( trim(gridname) /= trim(data_table(i)%gridname)) cycle
       if( trim(fieldname_code) /= trim(data_table(i)%fieldname_code)) cycle
       index1 = i                               ! field found
       exit
    enddo
    if(index1 .eq. -1) then
       if(debug_data_override) &
         call mpp_error(WARNING,'this field is NOT found in data_table: '//trim(fieldname_code))
       return  ! NO override was performed
    endif
  endif

  fieldname = data_table(index1)%fieldname_file ! fieldname in netCDF data file
  factor = data_table(index1)%factor
  multifile = data_table(index1)%multifile

  if(fieldname == "") then
     data_out = factor
     if(PRESENT(override)) override = .true.
     return
  else
     filename = data_table(index1)%file_name
     if (filename == "") call mpp_error(FATAL,'data_override: filename not given in data_table')
     if (multifile) prevfilename = data_table(index1)%prev_file_name
     if (multifile) nextfilename = data_table(index1)%next_file_name
  endif

!3 Check if fieldname has been previously processed
!$OMP SINGLE
  curr_position = -1
  if(num_fields > 0 ) then
     do i = 1, num_fields
        if(trim(override_array(i)%gridname) /= trim(gridname))   cycle
        if(trim(override_array(i)%fieldname) /= trim(fieldname_code)) cycle
        curr_position = i
        exit
     enddo
  endif

  if(curr_position < 0) then ! the field has not been processed previously
     num_fields = num_fields + 1
     curr_position = num_fields
     ! record fieldname, gridname in override_array
     override_array(curr_position)%fieldname = fieldname_code
     override_array(curr_position)%gridname = gridname
     id_time = init_external_field(filename,fieldname,verbose=.false.)
     if(id_time<0) call mpp_error(FATAL,'data_override:field not found in init_external_field 1')
     override_array(curr_position)%t_index = id_time
  else !curr_position >0
     !9 Get id_time  previously stored in override_array
     id_time = override_array(curr_position)%t_index
  endif !if curr_position < 0


  ! if using consecutive files for data_override, get time axis for previous and next files
  ! and check spatial dims for consistency
  if_multi1: if (multifile) then
    id_time_prev = -1
    if_prev1: if (trim(prevfilename) /= '') then
      id_time_prev = init_external_field(prevfilename,fieldname,verbose=.false.)
      dims = get_external_field_size(id_time)
      prev_dims = get_external_field_size(id_time_prev)
      ! check consistency of spatial dims
      if ((prev_dims(1) .ne. dims(1)) .or. (prev_dims(2) .ne. dims(2)) .or. &
          (prev_dims(3) .ne. dims(3))) then
        call mpp_error(FATAL, 'data_override: dimensions mismatch between consecutive forcing files')
      endif
      allocate(data_table(index1)%time_prev_records(prev_dims(4)))
      call get_time_axis(id_time_prev,data_table(index1)%time_prev_records)
    endif if_prev1
    id_time_next = -1
    if_next1: if (trim(nextfilename) /= '') then
      id_time_next = init_external_field(nextfilename,fieldname,verbose=.false.)
      dims = get_external_field_size(id_time)
      next_dims = get_external_field_size(id_time_next)
      ! check consistency of spatial dims
      if ((next_dims(1) .ne. dims(1)) .or. (next_dims(2) .ne. dims(2)) .or. &
          (next_dims(3) .ne. dims(3))) then
        call mpp_error(FATAL, 'data_override: dimensions mismatch between consecutive forcing files')
      endif
      allocate(data_table(index1)%time_next_records(next_dims(4)))
      call get_time_axis(id_time_next,data_table(index1)%time_next_records)
    endif if_next1
  endif if_multi1


  !10 do time interp to get data in compute_domain

  ! if using consecutive files, allow to perform time interpolation between the last record of previous
  ! file and first record of current file OR between the last record of current file and first record of
  ! next file hence "bridging" over files.
  if_multi2: if (multifile) then
    dims = get_external_field_size(id_time)
    if (.not. allocated(data_table(index1)%time_records)) allocate(data_table(index1)%time_records(dims(4)))
    call get_time_axis(id_time,data_table(index1)%time_records)

    first_record = data_table(index1)%time_records(1)
    last_record = data_table(index1)%time_records(dims(4))

    if_time2: if (time<first_record) then
      if (id_time_prev<0) call mpp_error(FATAL,'data_override:previous file needed with multifile')
      prev_dims = get_external_field_size(id_time_prev)
      if (time<data_table(index1)%time_prev_records(prev_dims(4))) call mpp_error(FATAL, &
          'data_override: time_interp_external_bridge should only be called to bridge with previous file')
      call time_interp_external_bridge(id_time_prev, id_time,time,data_out,verbose=.false.)
    elseif (time>last_record) then
      if (id_time_next<0) call mpp_error(FATAL,'data_override:next file needed with multifile')
      if (time>data_table(index1)%time_next_records(1)) call mpp_error(FATAL, &
          'data_override: time_interp_external_bridge should only be called to bridge with next file')
      call time_interp_external_bridge(id_time, id_time_next,time,data_out,verbose=.false.)
    else ! first_record < time < last_record, do not use bridge
      call time_interp_external(id_time,time,data_out,verbose=.false.)
    endif if_time2
  else ! standard behavior
     call time_interp_external(id_time,time,data_out,verbose=.false.)
  endif if_multi2


  data_out = data_out*factor
!$OMP END SINGLE

  if(PRESENT(override)) override = .true.

end subroutine DATA_OVERRIDE_0D_

!> @brief This routine performs data override for 2D fields.
subroutine DATA_OVERRIDE_2D_(gridname,fieldname,data_2D,time,override, is_in, ie_in, js_in, je_in)
  character(len=3), intent(in) :: gridname !< model grid ID
  character(len=*), intent(in) :: fieldname !< field to override
  logical, intent(out), optional :: override !< true if the field has been overriden succesfully
  type(time_type), intent(in) :: time !<  model time
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:), intent(inout) :: data_2D !< data returned by this call
  integer,           optional,  intent(in) :: is_in, ie_in, js_in, je_in
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:,:), allocatable ::  data_3D
  integer       :: index1
  integer       :: i

!1  Look  for the data file in data_table
  if(PRESENT(override)) override = .false.
  index1 = -1
  do i = 1, table_size
     if( trim(gridname) /= trim(data_table(i)%gridname)) cycle
     if( trim(fieldname) /= trim(data_table(i)%fieldname_code)) cycle
     index1 = i                               ! field found
     exit
  enddo
  if(index1 .eq. -1) return  ! NO override was performed

  allocate(data_3D(size(data_2D,1),size(data_2D,2),1))
  data_3D(:,:,1) = data_2D
  call DATA_OVERRIDE_3D_(gridname,fieldname,data_3D,time,override,data_index=index1,&
                       is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in)

  data_2D(:,:) = data_3D(:,:,1)
  deallocate(data_3D)
end subroutine DATA_OVERRIDE_2D_

!> @brief This routine performs data override for 3D fields
subroutine DATA_OVERRIDE_3D_(gridname,fieldname_code,return_data,time,override,data_index, is_in, ie_in, js_in, je_in)
  character(len=3),             intent(in) :: gridname !< model grid ID
  character(len=*),             intent(in) :: fieldname_code !< field name as used in the model
  logical,           optional, intent(out) :: override !< true if the field has been overriden succesfully
  type(time_type),              intent(in) :: time !< (target) model time
  integer,           optional,  intent(in) :: data_index
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:,:), intent(inout) :: return_data !< data returned by this call
  integer,           optional,  intent(in) :: is_in, ie_in, js_in, je_in
  logical, dimension(:,:,:),   allocatable :: mask_out

  character(len=512) :: filename !< file containing source data
  character(len=512) :: filename2 !< file containing source data
  character(len=512) :: prevfilename !< file containing source data for previous file
  character(len=512) :: prevfilename2 !< file containing source data for previous file
  character(len=512) :: nextfilename !< file containing source data for next file
  character(len=512) :: nextfilename2 !< file containing source data for next file
  character(len=128) :: fieldname !< fieldname used in the data file
  integer            :: i,j
  integer            :: dims(4)
  integer            :: prev_dims(4) !< dimensions of previous source data, when using multiple files
  integer            :: next_dims(4) !< dimensions of next source data, when using multiple files
  integer            :: index1 !< field index in data_table
  integer            :: id_time !< index for time interp in override array
  integer            :: id_time_prev=-1 !< time index for previous file, when using multiple files
  integer            :: id_time_next=-1 !< time index for next file, when using multiple files
  integer            :: axis_sizes(4)
  character(len=32)  :: axis_names(4)
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:), pointer :: lon_local =>NULL() !< of output (target) grid cells
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:), pointer :: lat_local =>NULL() !< of output (target) grid cells
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:), allocatable :: lon_tmp, lat_tmp

  logical :: data_file_is_2D = .false.  !< data in netCDF file is 2D
  logical :: ongrid, use_comp_domain
  type(domain2D) :: domain
  type(time_type)    :: first_record !< first record of "current" file
  type(time_type)    :: last_record !< last record of "current" file
  integer :: curr_position !< position of the field currently processed in override_array
  real(FMS_DATA_OVERRIDE_KIND_) :: factor
  integer, dimension(4) :: comp_domain = 0  !< istart,iend,jstart,jend for compute domain
  integer :: nxd, nyd, nxc, nyc, nwindows
  integer :: nwindows_x, ipos, jpos, window_size(2)
  integer :: istart, iend, jstart, jend
  integer :: isw, iew, jsw, jew
  integer :: omp_get_num_threads, window_id
  logical :: need_compute
  real(FMS_DATA_OVERRIDE_KIND_) :: lat_min, lat_max
  integer :: is_src, ie_src, js_src, je_src
  logical :: exists
  logical :: multifile !< use multiple consecutive files for override
  type(FmsNetcdfFile_t) :: fileobj
  integer :: startingi !< Starting x index for the compute domain relative to the input buffer
  integer :: endingi !< Ending x index for the compute domain relative to the input buffer
  integer :: startingj !< Starting y index for the compute domain relative to the input buffer
  integer :: endingj !< Ending y index for the compute domain relative to the input buffer
  integer :: nhalox !< Number of halos in the x direction
  integer :: nhaloy !< Number of halos in the y direction

  use_comp_domain = .false.
  if(.not.module_is_initialized) &
       call mpp_error(FATAL,'Error: need to call data_override_init first')

!1  Look  for the data file in data_table
  if(PRESENT(override)) override = .false.
  if (present(data_index)) then
    index1 = data_index
  else
    index1 = -1
    do i = 1, table_size
       if( trim(gridname) /= trim(data_table(i)%gridname)) cycle
       if( trim(fieldname_code) /= trim(data_table(i)%fieldname_code)) cycle
       index1 = i                               ! field found
       exit
    enddo
    if(index1 .eq. -1) then
       if(debug_data_override) &
         call mpp_error(WARNING,'this field is NOT found in data_table: '//trim(fieldname_code))
       return  ! NO override was performed
    endif
  endif

  fieldname = data_table(index1)%fieldname_file ! fieldname in netCDF data file
  factor = data_table(index1)%factor
  multifile = data_table(index1)%multifile

  if(fieldname == "") then
     return_data = factor
     if(PRESENT(override)) override = .true.
     return
  else
     filename = data_table(index1)%file_name
     if (filename == "") call mpp_error(FATAL,'data_override: filename not given in data_table')
     if (multifile) prevfilename = data_table(index1)%prev_file_name
     if (multifile) nextfilename = data_table(index1)%next_file_name
  endif

  ongrid = (data_table(index1)%interpol_method == 'none')

!3 Check if fieldname has been previously processed
!$OMP CRITICAL
  curr_position = -1
  if(num_fields > 0 ) then
     do i = 1, num_fields
        if(trim(override_array(i)%gridname) /= trim(gridname))   cycle
        if(trim(override_array(i)%fieldname) /= trim(fieldname_code)) cycle
        curr_position = i
        exit
     enddo
  endif

  if(curr_position < 0) then ! the field has not been processed previously
     num_fields = num_fields + 1
     curr_position = num_fields

! Get working domain from model's gridname
     call get_domain(gridname,domain,comp_domain)
     call mpp_get_data_domain(domain, xsize=nxd, ysize=nyd)
     nxc = comp_domain(2)-comp_domain(1) + 1
     nyc = comp_domain(4)-comp_domain(3) + 1

! record fieldname, gridname in override_array
     override_array(curr_position)%fieldname = fieldname_code
     override_array(curr_position)%gridname = gridname
     override_array(curr_position)%comp_domain = comp_domain
! get number of threads
     override_array(curr_position)%numthreads = 1
#if defined(_OPENMP)
     override_array(curr_position)%numthreads = omp_get_num_threads()
#endif
!--- data_override may be called from physics windows. The following are possible situations
!--- 1. size(return_data,1) == nxd and size(return_data,2) == nyd
!---    (on return_data domain and there is only one window).
!--- 2. nxc is divisible by size(return_data,1), nyc is divisible by size(return_data,2),
!---    nwindow = (nxc/size(return_data(1))*(nyc/size(return_data,2)),
!---    also we require nwindows is divisible by nthreads.
!---    The another restrition is that size(return_data,1) == ie_in - is_in + 1,
!---                                   size(return_data,2) == je_in - js_in + 1
     nwindows = 1
     if( nxd == size(return_data,1) .AND. nyd == size(return_data,2) ) then  !
        use_comp_domain = .false.
     else if ( mod(nxc, size(return_data,1)) ==0 .AND. mod(nyc, size(return_data,2)) ==0 ) then
        use_comp_domain = .true.
        nwindows = (nxc/size(return_data,1))*(nyc/size(return_data,2))
     else
        call mpp_error(FATAL, &
                     & "data_override: data is not on data domain and compute domain is not divisible by size(data)")
     endif
     override_array(curr_position)%window_size(1) = size(return_data,1)
     override_array(curr_position)%window_size(2) = size(return_data,2)

     window_size = override_array(curr_position)%window_size
     override_array(curr_position)%numwindows = nwindows
     if( mod(nwindows, override_array(curr_position)%numthreads) .NE. 0 ) then
        call mpp_error(FATAL, "data_override: nwindow is not divisible by nthreads")
     endif
     allocate(override_array(curr_position)%need_compute(nwindows))
     override_array(curr_position)%need_compute = .true.

!4 get index for time interp
     if(ongrid) then
        if( data_table(index1)%region_type .NE. NO_REGION ) then
           call mpp_error(FATAL,'data_override: ongrid must be false when region_type .NE. NO_REGION')
        endif

!  Allow on-grid data_overrides on cubed sphere grid
        inquire(file=trim(filename),EXIST=exists)
        if (.not. exists) then
           call get_mosaic_tile_file(filename,filename2,.false.,domain)
           filename = filename2
        endif

        ! if using consecutive files for data_override, get file names
        if_multi3: if (multifile) then
          if_prev3: if (trim(prevfilename) /= '') then
            inquire(file=trim(prevfilename),EXIST=exists)
            if (.not. exists) then
               call get_mosaic_tile_file(prevfilename,prevfilename2,.false.,domain)
               prevfilename = prevfilename2
            endif
          endif if_prev3
          if_next3: if (trim(nextfilename) /= '') then
            inquire(file=trim(nextfilename),EXIST=exists)
            if (.not. exists) then
               call get_mosaic_tile_file(nextfilename,nextfilename2,.false.,domain)
               nextfilename = nextfilename2
            endif
          endif if_next3
        endif if_multi3

        !--- we always only pass data on compute domain
        id_time = init_external_field(filename,fieldname,domain=domain,verbose=.false., &
                                    use_comp_domain=use_comp_domain, nwindows=nwindows, ongrid=ongrid)

        ! if using consecutive files for data_override, get time axis for previous and next files
        ! and check spatial dims for consistency
        if_multi4: if (multifile) then
          id_time_prev = -1
          if_prev4:if (trim(prevfilename) /= '') then
            id_time_prev = init_external_field(prevfilename,fieldname,domain=domain, &
               verbose=.false.,use_comp_domain=use_comp_domain, &
               nwindows = nwindows, ongrid=ongrid)
            dims = get_external_field_size(id_time)
            prev_dims = get_external_field_size(id_time_prev)
            ! check consistency of spatial dims
            if ((prev_dims(1) .ne. dims(1)) .or. (prev_dims(2) .ne. dims(2)) .or. &
                (prev_dims(3) .ne. dims(3))) then
              call mpp_error(FATAL, 'data_override: dimensions mismatch between consecutive forcing files')
            endif
            allocate(data_table(index1)%time_prev_records(prev_dims(4)))
            call get_time_axis(id_time_prev,data_table(index1)%time_prev_records)
          endif if_prev4
          id_time_next = -1
          if_next4: if (trim(nextfilename) /= '') then
            id_time_next = init_external_field(nextfilename,fieldname,domain=domain, &
               verbose=.false.,use_comp_domain=use_comp_domain, &
               nwindows = nwindows, ongrid=ongrid)
            dims = get_external_field_size(id_time)
            next_dims = get_external_field_size(id_time_next)
            ! check consistency of spatial dims
            if ((next_dims(1) .ne. dims(1)) .or. (next_dims(2) .ne. dims(2)) .or. &
                (next_dims(3) .ne. dims(3))) then
              call mpp_error(FATAL, 'data_override: dimensions mismatch between consecutive forcing files')
            endif
            allocate(data_table(index1)%time_next_records(next_dims(4)))
            call get_time_axis(id_time_next,data_table(index1)%time_next_records)
          endif if_next4
        endif if_multi4

        dims = get_external_field_size(id_time)
        override_array(curr_position)%dims = dims
        if(id_time<0) call mpp_error(FATAL,'data_override:field not found in init_external_field 1')
        override_array(curr_position)%t_index = id_time
        override_array(curr_position)%pt_index = id_time_prev
        override_array(curr_position)%nt_index = id_time_next
     else !ongrid=false
        id_time = init_external_field(filename,fieldname,domain=domain, axis_names=axis_names,&
            axis_sizes=axis_sizes, verbose=.false.,override=.true.,use_comp_domain=use_comp_domain, &
            nwindows = nwindows)

        ! if using consecutive files for data_override, get time axis for previous and next files
        ! and check spatial dims for consistency
        if_multi5: if (multifile) then
          id_time_prev = -1
          if_prev5: if (trim(prevfilename) /= '') then
            id_time_prev = init_external_field(prevfilename,fieldname,domain=domain, axis_names=axis_names,&
               axis_sizes=axis_sizes, verbose=.false.,override=.true.,use_comp_domain=use_comp_domain, &
               nwindows = nwindows)
            prev_dims = get_external_field_size(id_time_prev)
            allocate(data_table(index1)%time_prev_records(prev_dims(4)))
            call get_time_axis(id_time_prev,data_table(index1)%time_prev_records)
          endif if_prev5
          id_time_next = -1
          if_next5: if (trim(nextfilename) /= '') then
            id_time_next = init_external_field(nextfilename,fieldname,domain=domain, axis_names=axis_names,&
               axis_sizes=axis_sizes, verbose=.false.,override=.true.,use_comp_domain=use_comp_domain, &
               nwindows = nwindows)
            next_dims = get_external_field_size(id_time_next)
            allocate(data_table(index1)%time_next_records(next_dims(4)))
            call get_time_axis(id_time_next,data_table(index1)%time_next_records)
          endif if_next5
        endif if_multi5

        dims = get_external_field_size(id_time)
        override_array(curr_position)%dims = dims
        if(id_time<0) call mpp_error(FATAL,'data_override:field not found in init_external_field 2')
        override_array(curr_position)%t_index = id_time
        override_array(curr_position)%pt_index = id_time_prev
        override_array(curr_position)%nt_index = id_time_next

        !  get lon and lat of the input (source) grid, assuming that axis%data contains
        !  lat and lon of the input grid (in degrees)

        allocate(override_array(curr_position)%horz_interp(nwindows))
        allocate(override_array(curr_position)%lon_in(axis_sizes(1)+1))
        allocate(override_array(curr_position)%lat_in(axis_sizes(2)+1))
        if(get_external_fileobj(filename, fileobj)) then
           call axis_edges(fileobj, axis_names(1), override_array(curr_position)%lon_in, &
              reproduce_null_char_bug_flag=reproduce_null_char_bug)
           call axis_edges(fileobj, axis_names(2), override_array(curr_position)%lat_in, &
              reproduce_null_char_bug_flag=reproduce_null_char_bug)
        else
           call mpp_error(FATAL,'data_override: file '//trim(filename)//' is not opened in time_interp_external')
        end if
! convert lon_in and lat_in from deg to radian
        override_array(curr_position)%lon_in = override_array(curr_position)%lon_in * real(DEG_TO_RAD, lkind)
        override_array(curr_position)%lat_in = override_array(curr_position)%lat_in * real(DEG_TO_RAD, lkind)

        !--- find the region of the source grid that cover the local model grid.
        !--- currently we only find the index range for j-direction because
        !--- of the cyclic condition in i-direction. The purpose of this is to
        !--- decrease the memory usage and increase the IO performance.
        select case(gridname)
        case('OCN')
           lon_local => lon_local_ocn; lat_local => lat_local_ocn
        case('ICE')
           lon_local => lon_local_ice; lat_local => lat_local_ice
        case('ATM')
           lon_local => lon_local_atm; lat_local => lat_local_atm
        case('LND')
           lon_local => lon_local_lnd; lat_local => lat_local_lnd
        case default
           call mpp_error(FATAL,'error: gridname not recognized in data_override')
        end select

        lat_min = minval(lat_local)
        lat_max = maxval(lat_local)
        is_src = 1
        ie_src = axis_sizes(1)
        js_src = 1
        je_src = axis_sizes(2)
        ! do j = 1, axis_sizes(2)+1
        !    if( override_array(curr_position)%lat_in(j) > lat_min ) exit
        !    js_src = j
        ! enddo
        ! do j = 1, axis_sizes(2)+1
        !    je_src = j
        !    if( override_array(curr_position)%lat_in(j) > lat_max ) exit
        ! enddo

        !--- bicubic interpolation need one extra point in each direction. Also add
        !--- one more point for because lat_in is in the corner but the interpolation
        !--- use center points.
        select case (data_table(index1)%interpol_method)
        case ('bilinear')
           js_src = max(1, js_src-1)
           je_src = min(axis_sizes(2), je_src+1)
        case ('bicubic')
           js_src = max(1, js_src-2)
           je_src = min(axis_sizes(2), je_src+2)
        end select
        override_array(curr_position)%is_src = is_src
        override_array(curr_position)%ie_src = ie_src
        override_array(curr_position)%js_src = js_src
        override_array(curr_position)%je_src = je_src
        call reset_src_data_region(id_time, is_src, ie_src, js_src, je_src)
        if (multifile) then
          if (trim(prevfilename) /= '') then
            call reset_src_data_region(id_time_prev, is_src, ie_src, js_src, je_src)
          endif
          if (trim(nextfilename) /= '') then
            call reset_src_data_region(id_time_next, is_src, ie_src, js_src, je_src)
          endif
        endif

!       Find the index of lon_start, lon_end, lat_start and lat_end in the input grid (nearest points)
        if( data_table(index1)%region_type .NE. NO_REGION ) then
           allocate( lon_tmp(axis_sizes(1)), lat_tmp(axis_sizes(2)) )
           call read_data(fileobj, axis_names(1), lon_tmp)
           call read_data(fileobj, axis_names(2), lat_tmp)
           ! limit lon_start, lon_end are inside lon_in
           !       lat_start, lat_end are inside lat_in
           if(data_table(index1)%lon_start < lon_tmp(1) .OR. data_table(index1)%lon_start .GT. lon_tmp(axis_sizes(1)))&
              call mpp_error(FATAL, "data_override: lon_start is outside lon_T")
           if( data_table(index1)%lon_end < lon_tmp(1) .OR. data_table(index1)%lon_end .GT. lon_tmp(axis_sizes(1))) &
              call mpp_error(FATAL, "data_override: lon_end is outside lon_T")
           if(data_table(index1)%lat_start < lat_tmp(1) .OR. data_table(index1)%lat_start .GT. lat_tmp(axis_sizes(2)))&
              call mpp_error(FATAL, "data_override: lat_start is outside lat_T")
           if( data_table(index1)%lat_end < lat_tmp(1) .OR. data_table(index1)%lat_end .GT. lat_tmp(axis_sizes(2))) &
              call mpp_error(FATAL, "data_override: lat_end is outside lat_T")
           istart = nearest_index(data_table(index1)%lon_start, lon_tmp)
           iend   = nearest_index(data_table(index1)%lon_end,   lon_tmp)
           jstart = nearest_index(data_table(index1)%lat_start, lat_tmp)
           jend   = nearest_index(data_table(index1)%lat_end,   lat_tmp)
           ! adjust the index according to is_src and js_src
           istart = istart - is_src + 1
           iend   = iend   - is_src + 1
           jstart = jstart - js_src + 1
           jend   = jend   - js_src + 1
           call set_override_region(id_time, data_table(index1)%region_type, istart, iend, jstart, jend)
           if (multifile) then
             if (trim(prevfilename) /= '') then
               call set_override_region(id_time_prev, data_table(index1)%region_type, istart, iend, jstart, jend)
             endif
             if (trim(nextfilename) /= '') then
               call set_override_region(id_time_next, data_table(index1)%region_type, istart, iend, jstart, jend)
             endif
           endif
           deallocate(lon_tmp, lat_tmp)
        endif

     endif !ongrid
  else !curr_position >0
     dims = override_array(curr_position)%dims
     comp_domain = override_array(curr_position)%comp_domain
     nxc = comp_domain(2)-comp_domain(1) + 1
     nyc = comp_domain(4)-comp_domain(3) + 1
     is_src      = override_array(curr_position)%is_src
     ie_src      = override_array(curr_position)%ie_src
     js_src      = override_array(curr_position)%js_src
     je_src      = override_array(curr_position)%je_src
     window_size = override_array(curr_position)%window_size
     !---make sure data size match window_size
     if( window_size(1) .NE. size(return_data,1) .OR. window_size(2) .NE. size(return_data,2) ) then
        call mpp_error(FATAL, "data_override: window_size does not match size(data)")
     endif
!9 Get id_time  previously stored in override_array
     id_time = override_array(curr_position)%t_index
     id_time_prev = override_array(curr_position)%pt_index
     id_time_next = override_array(curr_position)%nt_index
  endif
!$OMP END CRITICAL

  if( override_array(curr_position)%numwindows > 1 ) then
      if( .NOT. PRESENT(is_in) .OR. .NOT. PRESENT(is_in) .OR. .NOT. PRESENT(is_in) .OR. .NOT. PRESENT(is_in) ) then
          call mpp_error(FATAL, "data_override: is_in, ie_in, js_in, je_in must be present when nwindows > 1")
      endif
  endif

  isw = comp_domain(1)
  iew = comp_domain(2)
  jsw = comp_domain(3)
  jew = comp_domain(4)
  window_id = 1
  if( override_array(curr_position)%numwindows > 1 ) then
     nxc = comp_domain(2) - comp_domain(1) + 1
     nwindows_x = nxc/window_size(1)
     ipos = (is_in-1)/window_size(1) + 1
     jpos = (js_in-1)/window_size(2)

     window_id = jpos*nwindows_x + ipos
     isw = isw + is_in - 1
     iew = isw + ie_in - is_in
     jsw = jsw + js_in - 1
     jew = jsw + je_in - js_in
  endif

  if( ongrid ) then
     need_compute = .false.
  else
     !--- find the index for windows.
     need_compute=override_array(curr_position)%need_compute(window_id)
  endif

  !--- call horiz_interp_new is not initialized

  if( need_compute ) then
     select case(gridname)
     case('OCN')
        lon_local => lon_local_ocn; lat_local => lat_local_ocn
     case('ICE')
        lon_local => lon_local_ice; lat_local => lat_local_ice
     case('ATM')
        lon_local => lon_local_atm; lat_local => lat_local_atm
     case('LND')
        lon_local => lon_local_lnd; lat_local => lat_local_lnd
     case default
        call mpp_error(FATAL,'error: gridname not recognized in data_override')
     end select

     select case (data_table(index1)%interpol_method)
     case ('bilinear')
        call horiz_interp_new (override_array(curr_position)%horz_interp(window_id), &
             override_array(curr_position)%lon_in(is_src:ie_src+1),                  &
             override_array(curr_position)%lat_in(js_src:je_src+1),                  &
             lon_local(isw:iew,jsw:jew), lat_local(isw:iew,jsw:jew), interp_method="bilinear")
     case ('bicubic')
        call horiz_interp_new (override_array(curr_position)%horz_interp(window_id), &
             override_array(curr_position)%lon_in(is_src:ie_src+1),                  &
             override_array(curr_position)%lat_in(js_src:je_src+1),                  &
             lon_local(isw:iew,jsw:jew), lat_local(isw:iew,jsw:jew), interp_method="bicubic")
     end select
     override_array(curr_position)%need_compute(window_id) = .false.
  endif

  ! Determine if  data in netCDF file is 2D or not
  data_file_is_2D = .false.
  if((dims(3) == 1) .and. (size(return_data,3)>1)) data_file_is_2D = .true.

  if(dims(3) .NE. 1 .and. (size(return_data,3) .NE. dims(3))) &
      call mpp_error(FATAL, "data_override: dims(3) .NE. 1 and size(return_data,3) .NE. dims(3)")


  dims = get_external_field_size(id_time)
  if (.not. allocated(data_table(index1)%time_records)) allocate(data_table(index1)%time_records(dims(4)))
  call get_time_axis(id_time,data_table(index1)%time_records)

  first_record = data_table(index1)%time_records(1)
  last_record = data_table(index1)%time_records(dims(4))

  if(ongrid) then
    if (.not. use_comp_domain) then
        !< Determine the size of the halox and the part of `data` that is in the compute domain
        nhalox = (size(return_data,1) - nxc)/2
        nhaloy = (size(return_data,2) - nyc)/2
        startingi = lbound(return_data,1) + nhalox
        startingj = lbound(return_data,2) + nhaloy
        endingi   = ubound(return_data,1) - nhalox
        endingj   = ubound(return_data,2) - nhaloy
    end if

!10 do time interp to get data in compute_domain
    if(data_file_is_2D) then
      if (use_comp_domain) then

        ! if using consecutive files, allow to perform time interpolation between the last record of previous
        ! file and first record of current file OR between the last record of current file and first record of
        ! next file hence "bridging" over files.
        if_multi6: if (multifile) then
          if_time6: if (time<first_record) then
            ! previous file must be init and time must be between last record of previous file and
            ! first record of current file
            if (id_time_prev<0) call mpp_error(FATAL,'data_override:previous file needed with multifile')
            prev_dims = get_external_field_size(id_time_prev)
            if (time<data_table(index1)%time_prev_records(prev_dims(4))) call mpp_error(FATAL, &
                'data_override: time_interp_external_bridge should only be called to bridge with previous file')
            ! bridge with previous file
            call time_interp_external_bridge(id_time_prev,id_time,time,return_data(:,:,1),verbose=.false., &
                                      is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          elseif (time>last_record) then
            ! next file must be init and time must be between last record of current file and
            ! first record of next file
            if (id_time_next<0) call mpp_error(FATAL,'data_override:next file needed with multifile')
            if (time>data_table(index1)%time_next_records(1)) call mpp_error(FATAL, &
                'data_override: time_interp_external_bridge should only be called to bridge with next file')
            ! bridge with next file
            call time_interp_external_bridge(id_time,id_time_next,time,return_data(:,:,1),verbose=.false., &
                                      is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          else  ! first_record <= time <= last_record, do not use bridge
            call time_interp_external(id_time,time,return_data(:,:,1),verbose=.false., &
                                      is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          endif if_time6
        else  ! standard behavior
          call time_interp_external(id_time,time,return_data(:,:,1),verbose=.false., &
                                    is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
        endif if_multi6

      else
        !> If this in an ongrid case and you are not in the compute domain, send in `data` to be the correct
        !! size

        ! if using consecutive files, allow to perform time interpolation between the last record of previous
        ! file and first record of current file OR between the last record of current file and first record of
        ! next file hence "bridging" over files.
        if_multi7: if (multifile) then
          if_time7: if (time<first_record) then
            ! previous file must be init and time must be between last record of previous file and
            ! first record of current file
            if (id_time_prev<0) call mpp_error(FATAL,'data_override:previous file needed with multifile')
            prev_dims = get_external_field_size(id_time_prev)
            if (time<data_table(index1)%time_prev_records(prev_dims(4))) call mpp_error(FATAL, &
                'data_override: time_interp_external_bridge should only be called to bridge with previous file')
            ! bridge with previous file
            call time_interp_external_bridge(id_time_prev,id_time,time,&
                                             return_data(startingi:endingi,startingj:endingj,1),verbose=.false., &
                                             is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          elseif (time>last_record) then
            ! next file must be init and time must be between last record of current file and
            ! first record of next file
            if (id_time_next<0) call mpp_error(FATAL,'data_override:next file needed with multifile')
            if (time>data_table(index1)%time_next_records(1)) call mpp_error(FATAL, &
                'data_override: time_interp_external_bridge should only be called to bridge with next file')
            ! bridge with next file
            call time_interp_external_bridge(id_time,id_time_next,time,&
                                             return_data(startingi:endingi,startingj:endingj,1),verbose=.false., &
                                             is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          else  ! first_record <= time <= last_record, do not use bridge
            call time_interp_external(id_time,time,return_data(startingi:endingi,startingj:endingj,1),verbose=.false., &
                                      is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          endif if_time7
        else  ! standard behavior
          call time_interp_external(id_time,time,return_data(startingi:endingi,startingj:endingj,1),verbose=.false., &
                                    is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
        endif if_multi7

      end if
      return_data(:,:,1) = return_data(:,:,1)*factor
      do i = 2, size(return_data,3)
        return_data(:,:,i) = return_data(:,:,1)
      end do
    else
      if (use_comp_domain) then

        ! if using consecutive files, allow to perform time interpolation between the last record of previous
        ! file and first record of current file OR between the last record of current file and first record of
        ! next file hence "bridging" over files.
        if_multi8: if (multifile) then
          if_time8: if (time<first_record) then
            if (id_time_prev<0) call mpp_error(FATAL,'data_override:previous file needed with multifile')
            prev_dims = get_external_field_size(id_time_prev)
            if (time<data_table(index1)%time_prev_records(prev_dims(4))) call mpp_error(FATAL, &
                'data_override: time_interp_external_bridge should only be called to bridge with previous file')
            call time_interp_external_bridge(id_time_prev,id_time,time,return_data,verbose=.false., &
                                             is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          elseif (time>last_record) then
            if (id_time_next<0) call mpp_error(FATAL,'data_override:next file needed with multifile')
            if (time>data_table(index1)%time_next_records(1)) call mpp_error(FATAL, &
                'data_override: time_interp_external_bridge should only be called to bridge with next file')
            call time_interp_external_bridge(id_time,id_time_next,time,return_data,verbose=.false., &
                                             is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          else ! first_record <= time <= last_record, do not use bridge
            call time_interp_external(id_time,time,return_data,verbose=.false., &
                                      is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          endif if_time8
        else ! standard behavior
          call time_interp_external(id_time,time,return_data,verbose=.false., &
                                    is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
       endif if_multi8

      else
        !> If this in an ongrid case and you are not in the compute domain, send in `data` to be the correct
        !! size

        ! if using consecutive files, allow to perform time interpolation between the last record of previous
        ! file and first record of current file OR between the last record of current file and first record of
        ! next file hence "bridging" over files.
        if_multi9: if (multifile) then
          if_time9: if (time<first_record) then
            if (id_time_prev<0) call mpp_error(FATAL,'data_override:previous file needed with multifile')
            prev_dims = get_external_field_size(id_time_prev)
            if (time<data_table(index1)%time_prev_records(prev_dims(4))) call mpp_error(FATAL, &
                'data_override: time_interp_external_bridge should only be called to bridge with previous file')
            call time_interp_external_bridge(id_time_prev,id_time,time,&
                                             return_data(startingi:endingi,startingj:endingj,:),verbose=.false., &
                                             is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          elseif (time>last_record) then
            if (id_time_next<0) call mpp_error(FATAL,'data_override:next file needed with multifile')
            if (time>data_table(index1)%time_next_records(1)) call mpp_error(FATAL, &
                'data_override: time_interp_external_bridge should only be called to bridge with next file')
            call time_interp_external_bridge(id_time,id_time_next,time,&
                                             return_data(startingi:endingi,startingj:endingj,:),verbose=.false., &
                                             is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          else ! first_record <= time <= last_record, do not use bridge
            call time_interp_external(id_time,time,return_data(startingi:endingi,startingj:endingj,:),verbose=.false., &
                                     is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
          endif if_time9
        else ! standard behavior
          call time_interp_external(id_time,time,return_data(startingi:endingi,startingj:endingj,:),verbose=.false., &
                                   is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
        endif if_multi9

      end if
      return_data = return_data*factor
    endif
  else  ! off grid case
! do time interp to get global data
     if(data_file_is_2D) then
        if( data_table(index1)%region_type == NO_REGION ) then

           ! if using consecutive files, allow to perform time interpolation between the last record of previous
           ! file and first record of current file OR between the last record of current file and first record of
           ! next file hence "bridging" over files.
           if_multi10: if (multifile) then
             if_time10: if (time<first_record) then
               if (id_time_prev<0) call mpp_error(FATAL,'data_override:previous file needed with multifile')
               prev_dims = get_external_field_size(id_time_prev)
               if (time<data_table(index1)%time_prev_records(prev_dims(4))) call mpp_error(FATAL, &
                   'data_override: time_interp_external_bridge should only be called to bridge with previous file')
               call time_interp_external_bridge(id_time_prev,id_time,time,return_data(:,:,1),verbose=.false., &
                                                horz_interp=override_array(curr_position)%horz_interp(window_id), &
                                                is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             elseif (time>last_record) then
               if (id_time_next<0) call mpp_error(FATAL,'data_override:next file needed with multifile')
               if (time>data_table(index1)%time_next_records(1)) call mpp_error(FATAL, &
                   'data_override: time_interp_external_bridge should only be called to bridge with next file')
               call time_interp_external_bridge(id_time,id_time_next,time,return_data(:,:,1),verbose=.false., &
                                                horz_interp=override_array(curr_position)%horz_interp(window_id), &
                                                is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             else ! first_record <= time <= last_record, do not use bridge
               call time_interp_external(id_time,time,return_data(:,:,1),verbose=.false., &
                                         horz_interp=override_array(curr_position)%horz_interp(window_id), &
                                         is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             endif if_time10
           else ! standard behavior
             call time_interp_external(id_time,time,return_data(:,:,1),verbose=.false., &
                                       horz_interp=override_array(curr_position)%horz_interp(window_id), &
                                       is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
           endif if_multi10

           return_data(:,:,1) = return_data(:,:,1)*factor
           do i = 2, size(return_data,3)
             return_data(:,:,i) = return_data(:,:,1)
           enddo
        else
           allocate(mask_out(size(return_data,1), size(return_data,2),1))
           mask_out = .false.

           ! if using consecutive files, allow to perform time interpolation between the last record of previous
           ! file and first record of current file OR between the last record of current file and first record of
           ! next file hence "bridging" over files.
           if_multi11: if (multifile) then
             if_time11: if (time<first_record) then
               if (id_time_prev<0) call mpp_error(FATAL,'data_override:previous file needed with multifile')
               prev_dims = get_external_field_size(id_time_prev)
               if (time<data_table(index1)%time_prev_records(prev_dims(4))) call mpp_error(FATAL, &
                   'data_override: time_interp_external_bridge should only be called to bridge with previous file')
               call time_interp_external_bridge(id_time_prev,id_time,time,return_data(:,:,1),verbose=.false., &
                     horz_interp=override_array(curr_position)%horz_interp(window_id),      &
                     mask_out   =mask_out(:,:,1), &
                     is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             elseif (time>last_record) then
               if (id_time_next<0) call mpp_error(FATAL,'data_override:next file needed with multifile')
               if (time>data_table(index1)%time_next_records(1)) call mpp_error(FATAL, &
                   'data_override: time_interp_external_bridge should only be called to bridge with next file')
               call time_interp_external_bridge(id_time,id_time_next,time,return_data(:,:,1),verbose=.false., &
                     horz_interp=override_array(curr_position)%horz_interp(window_id),      &
                     mask_out   =mask_out(:,:,1), &
                     is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             else ! first_record <= time <= last_record, do not use bridge
               call time_interp_external(id_time,time,return_data(:,:,1),verbose=.false., &
                    horz_interp=override_array(curr_position)%horz_interp(window_id),      &
                    mask_out   =mask_out(:,:,1), &
                    is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             endif if_time11
           else ! standard behavior
             call time_interp_external(id_time,time,return_data(:,:,1),verbose=.false., &
                   horz_interp=override_array(curr_position)%horz_interp(window_id),      &
                   mask_out   =mask_out(:,:,1), &
                   is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
           endif if_multi11

           where(mask_out(:,:,1))
              return_data(:,:,1) = return_data(:,:,1)*factor
           end where
           do i = 2, size(return_data,3)
              where(mask_out(:,:,1))
              return_data(:,:,i) = return_data(:,:,1)
              end where
           enddo
           deallocate(mask_out)
        endif
     else
        if( data_table(index1)%region_type == NO_REGION ) then

           ! if using consecutive files, allow to perform time interpolation between the last record of previous
           ! file and first record of current file OR between the last record of current file and first record of
           ! next file hence "bridging" over files.
           if_multi12: if (multifile) then
             if_time12: if (time<first_record) then
               if (id_time_prev<0) call mpp_error(FATAL,'data_override:previous file needed with multifile')
               prev_dims = get_external_field_size(id_time_prev)
               if (time<data_table(index1)%time_prev_records(prev_dims(4))) call mpp_error(FATAL, &
                   'data_override: time_interp_external_bridge should only be called to bridge with previous file')
               call time_interp_external_bridge(id_time_prev,id_time,time,return_data,verbose=.false.,      &
                                                horz_interp=override_array(curr_position)%horz_interp(window_id), &
                                                is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             elseif (time>last_record) then
               if (id_time_next<0) call mpp_error(FATAL,'data_override:next file needed with multifile')
               if (time>data_table(index1)%time_next_records(1)) call mpp_error(FATAL, &
                   'data_override: time_interp_external_bridge should only be called to bridge with next file')
               call time_interp_external_bridge(id_time,id_time_next,time,return_data,verbose=.false.,      &
                                                horz_interp=override_array(curr_position)%horz_interp(window_id), &
                                                is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             else ! first_record <= time <= last_record, do not use bridge
               call time_interp_external(id_time,time,return_data,verbose=.false.,      &
                    horz_interp=override_array(curr_position)%horz_interp(window_id), &
                    is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             endif if_time12
           else ! standard behavior
             call time_interp_external(id_time,time,return_data,verbose=.false.,      &
                  horz_interp=override_array(curr_position)%horz_interp(window_id), &
                  is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
           endif if_multi12

           return_data = return_data*factor
        else
           allocate(mask_out(size(return_data,1), size(return_data,2), size(return_data,3)) )
           mask_out = .false.

           ! if using consecutive files, allow to perform time interpolation between the last record of previous
           ! file and first record of current file OR between the last record of current file and first record of
           ! next file hence "bridging" over files.
           if_multi13: if (multifile) then
             if_time13: if (time<first_record) then
               if (id_time_prev<0) call mpp_error(FATAL,'data_override:previous file needed with multifile')
               prev_dims = get_external_field_size(id_time_prev)
               if (time<data_table(index1)%time_prev_records(prev_dims(4))) call mpp_error(FATAL, &
                   'data_override: time_interp_external_bridge should only be called to bridge with previous file')
               call time_interp_external_bridge(id_time_prev,id_time,time,return_data,verbose=.false.,      &
                                                horz_interp=override_array(curr_position)%horz_interp(window_id), &
                                                mask_out   =mask_out, &
                                                is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             elseif (time>last_record) then
               if (id_time_next<0) call mpp_error(FATAL,'data_override:next file needed with multifile')
               if (time>data_table(index1)%time_next_records(1)) call mpp_error(FATAL, &
                   'data_override: time_interp_external_bridge should only be called to bridge with next file')
               call time_interp_external_bridge(id_time,id_time_next,time,return_data,verbose=.false.,      &
                                                horz_interp=override_array(curr_position)%horz_interp(window_id), &
                                                mask_out   =mask_out, &
                                                is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             else ! first_record <= time <= last_record, do not use bridge
               call time_interp_external(id_time,time,return_data,verbose=.false.,      &
                    horz_interp=override_array(curr_position)%horz_interp(window_id),    &
                    mask_out   =mask_out, &
                    is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
             endif if_time13
           else ! standard behavior
             call time_interp_external(id_time,time,return_data,verbose=.false.,      &
                  horz_interp=override_array(curr_position)%horz_interp(window_id),    &
                  mask_out   =mask_out, &
                  is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
           endif if_multi13

           where(mask_out)
           return_data = return_data*factor
           end where
           deallocate(mask_out)
        endif
     endif

  endif

  if(PRESENT(override)) override = .true.
end subroutine DATA_OVERRIDE_3D_

!> @brief Data override for 1D unstructured grids
subroutine DATA_OVERRIDE_UG_1D_(gridname,fieldname,return_data,time,override)
  character(len=3),   intent(in) :: gridname !< model grid ID
  character(len=*),   intent(in) :: fieldname !< field to override
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:), intent(inout) :: return_data !< data returned by this call
  type(time_type),    intent(in) :: time !<  model time
  logical, intent(out), optional :: override !< true if the field has been overriden succesfully
  !local vars
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:), allocatable ::  data_SG
  type(domainUG) :: UG_domain
  integer       :: index1
  integer       :: i
  integer, dimension(4) :: comp_domain = 0  !< istart,iend,jstart,jend for compute domain

  !1  Look  for the data file in data_table
  if(PRESENT(override)) override = .false.
  index1 = -1
  do i = 1, table_size
     if( trim(gridname) /= trim(data_table(i)%gridname)) cycle
     if( trim(fieldname) /= trim(data_table(i)%fieldname_code)) cycle
     index1 = i                               ! field found
     exit
  enddo
  if(index1 .eq. -1) return  ! NO override was performed

  call get_domainUG(gridname,UG_domain,comp_domain)
  allocate(data_SG(comp_domain(1):comp_domain(2),comp_domain(3):comp_domain(4)))

  call DATA_OVERRIDE_2D_(gridname,fieldname,data_SG,time,override)

  call mpp_pass_SG_to_UG(UG_domain, data_SG(:,:), return_data(:))

  deallocate(data_SG)
end subroutine DATA_OVERRIDE_UG_1D_

!> @brief Data override for 2D unstructured grids
subroutine DATA_OVERRIDE_UG_2D_(gridname,fieldname,return_data,time,override)
  character(len=3),     intent(in) :: gridname !< model grid ID
  character(len=*),     intent(in) :: fieldname !< field to override
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:), intent(inout) :: return_data !< data returned by this call
  type(time_type),      intent(in) :: time !<  model time
  logical, intent(out), optional :: override !< true if the field has been overriden succesfully
  !local vars
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:,:), allocatable ::  data_SG
  real(FMS_DATA_OVERRIDE_KIND_), dimension(:,:),   allocatable ::  data_UG
  type(domainUG) :: UG_domain
  integer       :: index1
  integer       :: i, nlevel, nlevel_max
  integer, dimension(4) :: comp_domain = 0  !< istart,iend,jstart,jend for compute domain

!1  Look  for the data file in data_table
  if(PRESENT(override)) override = .false.
  index1 = -1
  do i = 1, table_size
     if( trim(gridname) /= trim(data_table(i)%gridname)) cycle
     if( trim(fieldname) /= trim(data_table(i)%fieldname_code)) cycle
     index1 = i                               ! field found
     exit
  enddo
  if(index1 .eq. -1) return  ! NO override was performed

  nlevel = size(return_data,2)
  nlevel_max = nlevel
  call mpp_max(nlevel_max)

  call get_domainUG(gridname,UG_domain,comp_domain)
  allocate(data_SG(comp_domain(1):comp_domain(2),comp_domain(3):comp_domain(4),nlevel_max))
  allocate(data_UG(size(return_data,1), nlevel_max))
  data_SG = 0._lkind
  call DATA_OVERRIDE_3D_(gridname,fieldname,data_SG,time,override)

  call mpp_pass_SG_to_UG(UG_domain, data_SG(:,:,:), data_UG(:,:))
  return_data(:,1:nlevel) = data_UG(:,1:nlevel)

  deallocate(data_SG, data_UG)
end subroutine DATA_OVERRIDE_UG_2D_
