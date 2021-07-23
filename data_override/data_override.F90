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
!> @defgroup data_override_mod data_override_mod
!> @ingroup data_override
!! @brief Routines to get data in a file whose path is described in a user-provided data_table
!! and do spatial and temporal interpolation if necessary to convert data to model's grid and time.
!! @author Z. Liang, M.J. Harrison, M. Winton
!!
!! Before using @ref data_override a data_table must be created with the following entries:
!! gridname, fieldname_code, fieldname_file, file_name, ongrid, factor.
!!
!! More explainations about data_table entries can be found in the source code (defining data_type)
!!
!! If user wants to override fieldname_code with a const, set fieldname_file in data_table = ""
!! and factor = const
!!
!! If user wants to override fieldname_code with data from a file, set fieldname_file = name in
!! the netCDF data file, factor then will be for unit conversion (=1 if no conversion required)
!!
!! A field can be overriden globally (by default) or users can specify one or two regions in which
!! data_override will take place, field values outside the region will not be affected.

!> @file
!> @brief File for @ref data_override_mod

module data_override_mod
use constants_mod, only: PI
use mpp_mod, only : mpp_error, FATAL, WARNING, stdout, stdlog, mpp_max
use mpp_mod, only : input_nml_file, get_unit
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_new, horiz_interp_type, &
                             assignment(=)
use time_interp_external2_mod, only:time_interp_external_init, &
                                   time_interp_external, &
                                   init_external_field, &
                                   get_external_field_size, &
                                   set_override_region, &
                                   reset_src_data_region, &
                                   NO_REGION, INSIDE_REGION, OUTSIDE_REGION,     &
                                   get_external_fileobj
use fms_mod, only: write_version_number, field_exist, lowercase, check_nml_error
use axis_utils_mod, only: get_axis_bounds
use axis_utils2_mod,  only : nearest_index, axis_edges
use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, NULL_DOMAIN2D,operator(.NE.),operator(.EQ.)
use mpp_domains_mod, only : mpp_get_global_domain, mpp_get_data_domain
use mpp_domains_mod, only : domainUG, mpp_pass_SG_to_UG, mpp_get_UG_SG_domain, NULL_DOMAINUG
use time_manager_mod, only: time_type
use fms2_io_mod,     only : FmsNetcdfFile_t, open_file, close_file, &
                            read_data, fms2_io_init, variable_exists, &
                            get_mosaic_tile_file
use get_grid_version_mod, only: get_grid_version_1, get_grid_version_2

implicit none
private

! Include variable "version" to be written to log file.
#include<file_version.h>

!> Private type for holding field and grid information from a data table
!> @ingroup data_override_mod
type data_type
   character(len=3)   :: gridname
   character(len=128) :: fieldname_code !< fieldname used in user's code (model)
   character(len=128) :: fieldname_file !< fieldname used in the netcdf data file
   character(len=512) :: file_name   !< name of netCDF data file
   character(len=128) :: interpol_method   !< interpolation method (default "bilinear")
   real               :: factor !< For unit conversion, default=1, see OVERVIEW above
   real               :: lon_start, lon_end, lat_start, lat_end
   integer            :: region_type
end type data_type

!> Private type for holding various data fields for performing data overrides
!> @ingroup data_override_mod
type override_type
   character(len=3)                 :: gridname
   character(len=128)               :: fieldname
   integer                          :: t_index                 !< index for time interp
   type(horiz_interp_type), allocatable :: horz_interp(:) !< index for horizontal spatial interp
   integer                          :: dims(4)                 !< dimensions(x,y,z,t) of the field in filename
   integer                          :: comp_domain(4)          !< istart,iend,jstart,jend for compute domain
   integer                          :: numthreads
   real, allocatable                :: lon_in(:)
   real, allocatable                :: lat_in(:)
   logical, allocatable             :: need_compute(:)
   integer                          :: numwindows
   integer                          :: window_size(2)
   integer                          :: is_src, ie_src, js_src, je_src
end type override_type

!> Interface for inserting and interpolating data into a file
!! for a model's grid and time. Data path must be described in
!! a user-provided data_table, see @ref data_override_mod "module description"
!! for more information.
!> @ingroup data_override_mod
interface data_override
     module procedure data_override_0d
     module procedure data_override_2d
     module procedure data_override_3d
end interface

!> Version of @ref data_override for unstructured grids
!> @ingroup data_override_mod
interface data_override_UG
     module procedure data_override_UG_1d
     module procedure data_override_UG_2d
end interface

!> @addtogroup data_override_mod
!> @{
 integer, parameter :: max_table=100, max_array=100
 real, parameter    :: deg_to_radian=PI/180.
 integer            :: table_size !< actual size of data table
 logical            :: module_is_initialized = .FALSE.

type(domain2D),save :: ocn_domain,atm_domain,lnd_domain, ice_domain
type(domainUG),save :: lnd_domainUG

real, dimension(:,:), target, allocatable :: lon_local_ocn, lat_local_ocn
real, dimension(:,:), target, allocatable :: lon_local_atm, lat_local_atm
real, dimension(:,:), target, allocatable :: lon_local_ice, lat_local_ice
real, dimension(:,:), target, allocatable :: lon_local_lnd, lat_local_lnd
real                                      :: min_glo_lon_ocn, max_glo_lon_ocn
real                                      :: min_glo_lon_atm, max_glo_lon_atm
real                                      :: min_glo_lon_lnd, max_glo_lon_lnd
real                                      :: min_glo_lon_ice, max_glo_lon_ice
integer:: num_fields = 0 !< number of fields in override_array already processed
type(data_type), dimension(max_table)           :: data_table !< user-provided data table
type(data_type)                                 :: default_table
type(override_type), dimension(max_array), save :: override_array !< to store processed fields
type(override_type), save                       :: default_array
logical                                         :: atm_on, ocn_on, lnd_on, ice_on
logical                                         :: lndUG_on
logical                                         :: debug_data_override
logical                                         :: grid_center_bug = .false.
logical                                         :: reproduce_null_char_bug = .false. !> Flag indicating
                                                   !! to reproduce the mpp_io bug where lat/lon_bnd were
                                                   !! not read correctly if null characters are present in
                                                   !! the netcdf file

namelist /data_override_nml/ debug_data_override, grid_center_bug, reproduce_null_char_bug


public :: data_override_init, data_override, data_override_unset_domains
public :: data_override_UG

contains
function count_ne_1(in_1, in_2, in_3)
  logical, intent(in)  :: in_1, in_2, in_3
  logical :: count_ne_1

  count_ne_1 = .not.(in_1.NEQV.in_2.NEQV.in_3) .OR. (in_1.AND.in_2.AND.in_3)
end function count_ne_1

!> @brief Assign default values for default_table, get domain of component models,
!! get global grids of component models.
!! Users should call data_override_init before calling data_override
!!
!! This subroutine should be called in coupler_init after
!! (ocean/atmos/land/ice)_model_init have been called.
!!
!! data_override_init can be called more than once, in one call the user can pass
!! up to 4 domains of component models, at least one domain must be present in
!! any call
!!
!! Data_table is initialized here with default values. Users should provide "real" values
!! that will override the default values. Real values can be given using data_table, each
!! line of data_table contains one data_entry. Items of data_entry are comma separated.
subroutine data_override_init(Atm_domain_in, Ocean_domain_in, Ice_domain_in, Land_domain_in, Land_domainUG_in)
  type (domain2d), intent(in), optional :: Atm_domain_in
  type (domain2d), intent(in), optional :: Ocean_domain_in, Ice_domain_in
  type (domain2d), intent(in), optional :: Land_domain_in
  type(domainUG) , intent(in), optional :: Land_domainUG_in

  character(len=128)    :: grid_file = 'INPUT/grid_spec.nc'
  integer               :: is,ie,js,je,use_get_grid_version
  integer               :: i, iunit, ntable, ntable_lima, ntable_new, unit,io_status, ierr
  character(len=256)    :: record
  logical               :: file_open
  logical               :: ongrid
  character(len=128)    :: region, region_type
  type(FmsNetcdfFile_t) :: fileobj

  type(data_type)  :: data_entry

  debug_data_override = .false.

  read (input_nml_file, data_override_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'data_override_nml')
  unit = stdlog()
  write(unit, data_override_nml)

!  if(module_is_initialized) return

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
    default_table%gridname = 'none'
    default_table%fieldname_code = 'none'
    default_table%fieldname_file = 'none'
    default_table%file_name = 'none'
    default_table%factor = 1.
    default_table%interpol_method = 'bilinear'
    do i = 1,max_table
       data_table(i) = default_table
    enddo

!  Read coupler_table
    iunit = get_unit()
    open(iunit, file='data_table', action='READ', iostat=io_status)
    if(io_status/=0) call mpp_error(FATAL, 'data_override_mod: Error in opening file data_table')

    ntable = 0
    ntable_lima = 0
    ntable_new = 0

    do while (ntable <= max_table)
       read(iunit,'(a)',end=100) record
       if (record(1:1) == '#') cycle
       if (record(1:10) == '          ') cycle
       ntable=ntable+1
       if (index(lowercase(record), "inside_region") .ne. 0 .or. index(lowercase(record), "outside_region") .ne. 0) then
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
                                   data_entry%file_name, data_entry%interpol_method, data_entry%factor, region, region_type
             if (data_entry%interpol_method == 'default') then
                data_entry%interpol_method = default_table%interpol_method
             endif
             if (.not.(data_entry%interpol_method == 'default'  .or. &
                  data_entry%interpol_method == 'bicubic'  .or. &
                  data_entry%interpol_method == 'bilinear' .or. &
                  data_entry%interpol_method == 'none')) then
                unit = stdout()
                write(unit,*)" gridname is ", trim(data_entry%gridname)
                write(unit,*)" fieldname_code is ", trim(data_entry%fieldname_code)
                write(unit,*)" fieldname_file is ", trim(data_entry%fieldname_file)
                write(unit,*)" file_name is ", trim(data_entry%file_name)
                write(unit,*)" factor is ", data_entry%factor
                write(unit,*)" interpol_method is ", trim(data_entry%interpol_method)
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
       else if (index(lowercase(record), ".false.") .ne. 0 .or. index(lowercase(record), ".true.") .ne. 0 ) then ! old format
          ntable_lima = ntable_lima + 1
          read(record,*,err=99) data_entry%gridname, data_entry%fieldname_code, data_entry%fieldname_file, &
                                   data_entry%file_name, ongrid, data_entry%factor
          if(ongrid) then
             data_entry%interpol_method = 'none'
          else
             data_entry%interpol_method = 'bilinear'
          endif
          data_entry%lon_start = 0.0
          data_entry%lon_end   = -1.0
          data_entry%lat_start = 0.0
          data_entry%lat_end   = -1.0
          data_entry%region_type = NO_REGION
       else                                      ! new format
          ntable_new=ntable_new+1
          read(record,*,err=99) data_entry%gridname, data_entry%fieldname_code, data_entry%fieldname_file, &
                                data_entry%file_name, data_entry%interpol_method, data_entry%factor
          if (data_entry%interpol_method == 'default') then
            data_entry%interpol_method = default_table%interpol_method
          endif
          if (.not.(data_entry%interpol_method == 'default'  .or. &
                    data_entry%interpol_method == 'bicubic'  .or. &
                    data_entry%interpol_method == 'bilinear' .or. &
                    data_entry%interpol_method == 'none')) then
             unit = stdout()
             write(unit,*)" gridname is ", trim(data_entry%gridname)
             write(unit,*)" fieldname_code is ", trim(data_entry%fieldname_code)
             write(unit,*)" fieldname_file is ", trim(data_entry%fieldname_file)
             write(unit,*)" file_name is ", trim(data_entry%file_name)
             write(unit,*)" factor is ", data_entry%factor
             write(unit,*)" interpol_method is ", trim(data_entry%interpol_method)
             call mpp_error(FATAL, 'data_override_mod: invalid last entry in data_override_table, ' &
                               //'its value should be "default", "bicubic", "bilinear" or "none" ')
          endif
          data_entry%lon_start = 0.0
          data_entry%lon_end   = -1.0
          data_entry%lat_start = 0.0
          data_entry%lat_end   = -1.0
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
     if(count_ne_1((ocn_on .OR. ice_on), lnd_on, atm_on)) call mpp_error(FATAL, 'data_override_mod: the grid file ' // &
          'is a solo mosaic, one and only one of atm_on, lnd_on or ice_on/ocn_on should be true')
   end if
 else
   call mpp_error(FATAL, 'data_override_mod: none of x_T, geolon_t, ocn_mosaic_file or gridfiles exist in '//trim(grid_file))
 endif

 if(use_get_grid_version .EQ. 1) then
    if (atm_on .and. .not. allocated(lon_local_atm) ) then
       call mpp_get_compute_domain( atm_domain,is,ie,js,je)
       allocate(lon_local_atm(is:ie,js:je), lat_local_atm(is:ie,js:je))
       call get_grid_version_1(grid_file, 'atm', atm_domain, is, ie, js, je, lon_local_atm, lat_local_atm, &
          min_glo_lon_atm, max_glo_lon_atm, grid_center_bug )
    endif
    if (ocn_on .and. .not. allocated(lon_local_ocn) ) then
       call mpp_get_compute_domain( ocn_domain,is,ie,js,je)
       allocate(lon_local_ocn(is:ie,js:je), lat_local_ocn(is:ie,js:je))
       call get_grid_version_1(grid_file, 'ocn', ocn_domain, is, ie, js, je, lon_local_ocn, lat_local_ocn, &
          min_glo_lon_ocn, max_glo_lon_ocn, grid_center_bug )
    endif

    if (lnd_on .and. .not. allocated(lon_local_lnd) ) then
       call mpp_get_compute_domain( lnd_domain,is,ie,js,je)
       allocate(lon_local_lnd(is:ie,js:je), lat_local_lnd(is:ie,js:je))
       call get_grid_version_1(grid_file, 'lnd', lnd_domain, is, ie, js, je, lon_local_lnd, lat_local_lnd, &
          min_glo_lon_lnd, max_glo_lon_lnd, grid_center_bug )
    endif

    if (ice_on .and. .not. allocated(lon_local_ice) ) then
       call mpp_get_compute_domain( ice_domain,is,ie,js,je)
       allocate(lon_local_ice(is:ie,js:je), lat_local_ice(is:ie,js:je))
       call get_grid_version_1(grid_file, 'ice', ice_domain, is, ie, js, je, lon_local_ice, lat_local_ice, &
          min_glo_lon_ice, max_glo_lon_ice, grid_center_bug )
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

end subroutine data_override_init

!> @brief Unset domains that had previously been set for use by data_override.
!!
!! This subroutine deallocates any data override domains that have been set.
subroutine data_override_unset_domains(unset_Atm, unset_Ocean, &
                                      unset_Ice, unset_Land, must_be_set)
  logical, intent(in), optional :: unset_Atm, unset_Ocean, unset_Ice, unset_Land
  logical, intent(in), optional :: must_be_set

  logical :: fail_if_not_set

  fail_if_not_set = .true. ; if (present(must_be_set)) fail_if_not_set = must_be_set

  if (.not.module_is_initialized) call mpp_error(FATAL, &
     "data_override_unset_domains called with an unititialized data_override module.")

  if (PRESENT(unset_Atm)) then ; if (unset_Atm) then
    if (fail_if_not_set .and. .not.atm_on) call mpp_error(FATAL, &
      "data_override_unset_domains attempted to work on an Atm_domain that had not been set.")
    atm_domain = NULL_DOMAIN2D
    atm_on = .false.
    if (allocated(lon_local_atm)) deallocate(lon_local_atm)
    if (allocated(lat_local_atm)) deallocate(lat_local_atm)
  endif ; endif
  if (PRESENT(unset_Ocean)) then ; if (unset_Ocean) then
    if (fail_if_not_set .and. .not.ocn_on) call mpp_error(FATAL, &
      "data_override_unset_domains attempted to work on an Ocn_domain that had not been set.")
    ocn_domain = NULL_DOMAIN2D
    ocn_on = .false.
    if (allocated(lon_local_ocn)) deallocate(lon_local_ocn)
    if (allocated(lat_local_ocn)) deallocate(lat_local_ocn)
  endif ; endif
  if (PRESENT(unset_Land)) then ; if (unset_Land) then
    if (fail_if_not_set .and. .not.lnd_on) call mpp_error(FATAL, &
      "data_override_unset_domains attempted to work on a Land_domain that had not been set.")
    lnd_domain = NULL_DOMAIN2D
    lnd_on = .false.
    if (allocated(lon_local_lnd)) deallocate(lon_local_lnd)
    if (allocated(lat_local_lnd)) deallocate(lat_local_lnd)
  endif ; endif
  if (PRESENT(unset_Ice)) then ; if (unset_Ice) then
    if (fail_if_not_set .and. .not.ice_on) call mpp_error(FATAL, &
      "data_override_unset_domains attempted to work on an Ice_domain that had not been set.")
    ice_domain = NULL_DOMAIN2D
    ice_on = .false.
    if (allocated(lon_local_ice)) deallocate(lon_local_ice)
    if (allocated(lat_local_ice)) deallocate(lat_local_ice)
  endif ; endif

end subroutine data_override_unset_domains

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

!> @brief This routine performs data override for 2D fields; for usage, see data_override_3d.
subroutine data_override_2d(gridname,fieldname,data_2D,time,override, is_in, ie_in, js_in, je_in)
  character(len=3), intent(in) :: gridname !< model grid ID
  character(len=*), intent(in) :: fieldname !< field to override
  logical, intent(out), optional :: override !< true if the field has been overriden succesfully
  type(time_type), intent(in) :: time !<  model time
  real, dimension(:,:), intent(inout) :: data_2D !< data returned by this call
  integer,           optional,  intent(in) :: is_in, ie_in, js_in, je_in
!  real, dimension(size(data_2D,1),size(data_2D,2),1) :: data_3D
  real, dimension(:,:,:), allocatable ::  data_3D
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
  call data_override_3d(gridname,fieldname,data_3D,time,override,data_index=index1,&
                       is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in)

  data_2D(:,:) = data_3D(:,:,1)
  deallocate(data_3D)
end subroutine data_override_2d

!> @brief This routine performs data override for 3D fields
subroutine data_override_3d(gridname,fieldname_code,data,time,override,data_index, is_in, ie_in, js_in, je_in)
  character(len=3),             intent(in) :: gridname !< model grid ID
  character(len=*),             intent(in) :: fieldname_code !< field name as used in the model
  logical,           optional, intent(out) :: override !< true if the field has been overriden succesfully
  type(time_type),              intent(in) :: time !< (target) model time
  integer,           optional,  intent(in) :: data_index
  real, dimension(:,:,:),    intent(inout) :: data !< data returned by this call
  integer,           optional,  intent(in) :: is_in, ie_in, js_in, je_in
  logical, dimension(:,:,:),   allocatable :: mask_out

  character(len=512) :: filename !< file containing source data
  character(len=512) :: filename2 !< file containing source data
  character(len=128) :: fieldname !< fieldname used in the data file
  integer            :: i,j
  integer            :: dims(4)
  integer            :: index1 !< field index in data_table
  integer            :: id_time !< index for time interp in override array
  integer            :: axis_sizes(4)
  character(len=32)  :: axis_names(4)
  real, dimension(:,:), pointer :: lon_local =>NULL() !< of output (target) grid cells
  real, dimension(:,:), pointer :: lat_local =>NULL() !< of output (target) grid cells
  real, dimension(:), allocatable :: lon_tmp, lat_tmp

  logical :: data_file_is_2D = .false.  !< data in netCDF file is 2D
  logical :: ongrid, use_comp_domain
  type(domain2D) :: domain
  integer :: curr_position !< position of the field currently processed in override_array
  real :: factor
  integer, dimension(4) :: comp_domain = 0  !< istart,iend,jstart,jend for compute domain
  integer :: nxd, nyd, nxc, nyc, nwindows
  integer :: nwindows_x, ipos, jpos, window_size(2)
  integer :: istart, iend, jstart, jend
  integer :: isw, iew, jsw, jew, n
  integer :: omp_get_num_threads, omp_get_thread_num, thread_id, window_id
  logical :: need_compute
  real    :: lat_min, lat_max
  integer :: is_src, ie_src, js_src, je_src
  logical :: exists
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

  if(fieldname == "") then
     data = factor
     if(PRESENT(override)) override = .true.
     return
  else
     filename = data_table(index1)%file_name
     if (filename == "") call mpp_error(FATAL,'data_override: filename not given in data_table')
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
!--- 1. size(data,1) == nxd and size(data,2) == nyd   ( on data domain and there is only one window).
!--- 2. nxc is divisible by size(data,1), nyc is divisible by size(data,2),
!---    nwindow = (nxc/size(data(1))*(nyc/size(data,2)), also we require nwindows is divisible by nthreads.
!---    The another restrition is that size(data,1) == ie_in - is_in + 1,
!---                                   size(data,2) == je_in - js_in + 1
     nwindows = 1
     if( nxd == size(data,1) .AND. nyd == size(data,2) ) then  !
        use_comp_domain = .false.
     else if ( mod(nxc, size(data,1)) ==0 .AND. mod(nyc, size(data,2)) ==0 ) then
        use_comp_domain = .true.
        nwindows = (nxc/size(data,1))*(nyc/size(data,2))
     else
        call mpp_error(FATAL, "data_override: data is not on data domain and compute domain is not divisible by size(data)")
     endif
     override_array(curr_position)%window_size(1) = size(data,1)
     override_array(curr_position)%window_size(2) = size(data,2)

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

        !--- we always only pass data on compute domain
        id_time = init_external_field(filename,fieldname,domain=domain,verbose=.false., &
                                    use_comp_domain=use_comp_domain, nwindows=nwindows, ongrid=ongrid)
        dims = get_external_field_size(id_time)
        override_array(curr_position)%dims = dims
        if(id_time<0) call mpp_error(FATAL,'data_override:field not found in init_external_field 1')
        override_array(curr_position)%t_index = id_time
     else !ongrid=false
        id_time = init_external_field(filename,fieldname,domain=domain, axis_names=axis_names,&
            axis_sizes=axis_sizes, verbose=.false.,override=.true.,use_comp_domain=use_comp_domain, &
            nwindows = nwindows)
        dims = get_external_field_size(id_time)
        override_array(curr_position)%dims = dims
        if(id_time<0) call mpp_error(FATAL,'data_override:field not found in init_external_field 2')
        override_array(curr_position)%t_index = id_time

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
        override_array(curr_position)%lon_in = override_array(curr_position)%lon_in * deg_to_radian
        override_array(curr_position)%lat_in = override_array(curr_position)%lat_in * deg_to_radian

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
        do j = 1, axis_sizes(2)+1
           if( override_array(curr_position)%lat_in(j) > lat_min ) exit
           js_src = j
        enddo
        do j = 1, axis_sizes(2)+1
           je_src = j
           if( override_array(curr_position)%lat_in(j) > lat_max ) exit
        enddo

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

!       Find the index of lon_start, lon_end, lat_start and lat_end in the input grid (nearest points)
        if( data_table(index1)%region_type .NE. NO_REGION ) then
           allocate( lon_tmp(axis_sizes(1)), lat_tmp(axis_sizes(2)) )
           call read_data(fileobj, axis_names(1), lon_tmp)
           call read_data(fileobj, axis_names(2), lat_tmp)
           ! limit lon_start, lon_end are inside lon_in
           !       lat_start, lat_end are inside lat_in
           if( data_table(index1)%lon_start < lon_tmp(1) .OR. data_table(index1)%lon_start .GT. lon_tmp(axis_sizes(1))) &
              call mpp_error(FATAL, "data_override: lon_start is outside lon_T")
           if( data_table(index1)%lon_end < lon_tmp(1) .OR. data_table(index1)%lon_end .GT. lon_tmp(axis_sizes(1))) &
              call mpp_error(FATAL, "data_override: lon_end is outside lon_T")
           if( data_table(index1)%lat_start < lat_tmp(1) .OR. data_table(index1)%lat_start .GT. lat_tmp(axis_sizes(2))) &
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
           deallocate(lon_tmp, lat_tmp)
        endif

     endif
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
     if( window_size(1) .NE. size(data,1) .OR. window_size(2) .NE. size(data,2) ) then
        call mpp_error(FATAL, "data_override: window_size does not match size(data)")
     endif
!9 Get id_time  previously stored in override_array
     id_time = override_array(curr_position)%t_index
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
  if((dims(3) == 1) .and. (size(data,3)>1)) data_file_is_2D = .true.

  if(dims(3) .NE. 1 .and. (size(data,3) .NE. dims(3))) &
      call mpp_error(FATAL, "data_override: dims(3) .NE. 1 and size(data,3) .NE. dims(3)")

  if(ongrid) then
    if (.not. use_comp_domain) then
        !< Determine the size of the halox and the part of `data` that is in the compute domain
        nhalox = (size(data,1) - nxc)/2
        nhaloy = (size(data,2) - nyc)/2
        startingi = lbound(data,1) + nhalox
        startingj = lbound(data,2) + nhaloy
        endingi   = ubound(data,1) - nhalox
        endingj   = ubound(data,2) - nhaloy
    end if

!10 do time interp to get data in compute_domain
    if(data_file_is_2D) then
      if (use_comp_domain) then
        call time_interp_external(id_time,time,data(:,:,1),verbose=.false., &
                                  is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
      else
        !> If this in an ongrid case and you are not in the compute domain, send in `data` to be the correct
        !! size
        call time_interp_external(id_time,time,data(startingi:endingi,startingj:endingj,1),verbose=.false., &
                               is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
      end if
      data(:,:,1) = data(:,:,1)*factor
      do i = 2, size(data,3)
        data(:,:,i) = data(:,:,1)
      end do
    else
      if (use_comp_domain) then
        call time_interp_external(id_time,time,data,verbose=.false., &
                                  is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
      else
        !> If this in an ongrid case and you are not in the compute domain, send in `data` to be the correct
        !! size
        call time_interp_external(id_time,time,data(startingi:endingi,startingj:endingj,:),verbose=.false., &
                               is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
      end if
      data = data*factor
    endif
  else  ! off grid case
! do time interp to get global data
     if(data_file_is_2D) then
        if( data_table(index1)%region_type == NO_REGION ) then
           call time_interp_external(id_time,time,data(:,:,1),verbose=.false., &
                 horz_interp=override_array(curr_position)%horz_interp(window_id), &
                 is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
           data(:,:,1) = data(:,:,1)*factor
           do i = 2, size(data,3)
             data(:,:,i) = data(:,:,1)
           enddo
        else
           allocate(mask_out(size(data,1), size(data,2),1))
           mask_out = .false.
           call time_interp_external(id_time,time,data(:,:,1),verbose=.false., &
                 horz_interp=override_array(curr_position)%horz_interp(window_id),      &
                 mask_out   =mask_out(:,:,1), &
                 is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
           where(mask_out(:,:,1))
              data(:,:,1) = data(:,:,1)*factor
           end where
           do i = 2, size(data,3)
              where(mask_out(:,:,1))
                data(:,:,i) = data(:,:,1)
              end where
           enddo
           deallocate(mask_out)
        endif
     else
        if( data_table(index1)%region_type == NO_REGION ) then
           call time_interp_external(id_time,time,data,verbose=.false.,      &
              horz_interp=override_array(curr_position)%horz_interp(window_id), &
              is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
           data = data*factor
        else
           allocate(mask_out(size(data,1), size(data,2), size(data,3)) )
           mask_out = .false.
           call time_interp_external(id_time,time,data,verbose=.false.,      &
              horz_interp=override_array(curr_position)%horz_interp(window_id),    &
              mask_out   =mask_out, &
              is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
           where(mask_out)
              data = data*factor
           end where
           deallocate(mask_out)
        endif
     endif

  endif

  if(PRESENT(override)) override = .true.

end subroutine data_override_3d

!> @brief Routine to perform data override for scalar fields
subroutine data_override_0d(gridname,fieldname_code,data,time,override,data_index)
  character(len=3), intent(in) :: gridname !< model grid ID (ocn,ice,atm,lnd)
  character(len=*), intent(in) :: fieldname_code !< field name as used in the model (may be
                                                 !! different from the name in NetCDF data file)
  logical, intent(out), optional :: override !< true if the field has been overriden succesfully
  type(time_type), intent(in) :: time !< (target) model time
  real,             intent(out) :: data !< output data array returned by this call
  integer, intent(in), optional :: data_index

  character(len=512) :: filename !< file containing source data
  character(len=128) :: fieldname !< fieldname used in the data file
  integer :: index1 !< field index in data_table
  integer :: id_time !< index for time interp in override array
  integer :: curr_position !< position of the field currently processed in override_array
  integer :: i
  real :: factor

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

  if(fieldname == "") then
     data = factor
     if(PRESENT(override)) override = .true.
     return
  else
     filename = data_table(index1)%file_name
     if (filename == "") call mpp_error(FATAL,'data_override: filename not given in data_table')
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

  !10 do time interp to get data in compute_domain
  call time_interp_external(id_time, time, data, verbose=.false.)
  data = data*factor
!$OMP END SINGLE

  if(PRESENT(override)) override = .true.

end subroutine data_override_0d

!> @brief Data override for 2D unstructured grids
subroutine data_override_UG_1d(gridname,fieldname,data,time,override)
  character(len=3),   intent(in) :: gridname !< model grid ID
  character(len=*),   intent(in) :: fieldname !< field to override
  real, dimension(:), intent(inout) :: data !< data returned by this call
  type(time_type),    intent(in) :: time !<  model time
  logical, intent(out), optional :: override !< true if the field has been overriden succesfully
  !local vars
  real, dimension(:,:), allocatable ::  data_SG
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

  call data_override_2d(gridname,fieldname,data_SG,time,override)

  call mpp_pass_SG_to_UG(UG_domain, data_SG(:,:), data(:))

  deallocate(data_SG)

end subroutine data_override_UG_1d

!> @brief Data override for 2D unstructured grids
subroutine data_override_UG_2d(gridname,fieldname,data,time,override)
  character(len=3),     intent(in) :: gridname !< model grid ID
  character(len=*),     intent(in) :: fieldname !< field to override
  real, dimension(:,:), intent(inout) :: data !< data returned by this call
  type(time_type),      intent(in) :: time !<  model time
  logical, intent(out), optional :: override !< true if the field has been overriden succesfully
  !local vars
  real, dimension(:,:,:), allocatable ::  data_SG
  real, dimension(:,:),   allocatable ::  data_UG
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

  nlevel = size(data,2)
  nlevel_max = nlevel
  call mpp_max(nlevel_max)

  call get_domainUG(gridname,UG_domain,comp_domain)
  allocate(data_SG(comp_domain(1):comp_domain(2),comp_domain(3):comp_domain(4),nlevel_max))
  allocate(data_UG(size(data,1), nlevel_max))
  data_SG = 0.0
  call data_override_3d(gridname,fieldname,data_SG,time,override)

  call mpp_pass_SG_to_UG(UG_domain, data_SG(:,:,:), data_UG(:,:))
  data(:,1:nlevel) = data_UG(:,1:nlevel)

  deallocate(data_SG, data_UG)

end subroutine data_override_UG_2d

end module data_override_mod
!> @}
! close documentation grouping
