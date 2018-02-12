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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module data_override_mod
!
! <CONTACT EMAIL="Zhi.Liang@noaa.gov">
! Z. Liang
! </CONTACT>
!
! <CONTACT EMAIL="Matthew.Harrison@noaa.gov">
!  M.J. Harrison 
! </CONTACT>
!
! <CONTACT EMAIL="Michael.Winton@noaa.gov">
! M. Winton
! </CONTACT>

!<OVERVIEW>
! Given a gridname, fieldname and model time this routine will get data in a file whose
! path is described in a user-provided data_table, do spatial and temporal interpolation if 
! necessary to convert data to model's grid and time.
!
! Before using data_override a data_table must be created with the following entries:
! gridname, fieldname_code, fieldname_file, file_name, ongrid, factor.
!
! More explainations about data_table entries can be found in the source code (defining data_type)
!
! If user wants to override fieldname_code with a const, set fieldname_file in data_table = ""
! and factor = const
!
! If user wants to override fieldname_code with data from a file, set fieldname_file = name in
! the netCDF data file, factor then will be for unit conversion (=1 if no conversion required)
!
! A field can be overriden globally (by default) or users can specify one or two regions in which
! data_override will take place, field values outside the region will not be affected. 
!</OVERVIEW>
#include <fms_platform.h>
use platform_mod, only: r8_kind
use constants_mod, only: PI
use mpp_io_mod, only: axistype,mpp_close,mpp_open,mpp_get_axis_data,MPP_RDONLY,MPP_ASCII
use mpp_mod, only : mpp_error,FATAL,WARNING,mpp_pe,stdout,stdlog,mpp_root_pe, NOTE, mpp_min, mpp_max, mpp_chksum
use mpp_mod, only : input_nml_file
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_new, horiz_interp_type, &
                             assignment(=), horiz_interp_del
use time_interp_external_mod, only:time_interp_external_init, time_interp_external, &
                                   init_external_field, get_external_field_size, &
                                   NO_REGION, INSIDE_REGION, OUTSIDE_REGION,     &
                                   set_override_region, reset_src_data_region
use fms_io_mod, only: field_size, read_data, fms_io_init,get_mosaic_tile_grid, get_mosaic_tile_file
use fms_mod, only: write_version_number, field_exist, lowercase, file_exist, open_namelist_file, check_nml_error, close_file
use axis_utils_mod, only: get_axis_bounds, nearest_index
use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, NULL_DOMAIN2D,operator(.NE.),operator(.EQ.)
use mpp_domains_mod, only : mpp_copy_domain, mpp_get_global_domain
use mpp_domains_mod, only : mpp_get_data_domain, mpp_set_compute_domain, mpp_set_data_domain
use mpp_domains_mod, only : mpp_set_global_domain, mpp_deallocate_domain
use mpp_domains_mod, only : domainUG, mpp_pass_SG_to_UG, mpp_get_UG_SG_domain, NULL_DOMAINUG

use time_manager_mod, only: time_type

implicit none
private

! Include variable "version" to be written to log file.
#include<file_version.h>

type data_type
   character(len=3)   :: gridname
   character(len=128) :: fieldname_code !fieldname used in user's code (model)
   character(len=128) :: fieldname_file ! fieldname used in the netcdf data file
   character(len=512) :: file_name   ! name of netCDF data file
   character(len=128) :: interpol_method   ! interpolation method (default "bilinear")
   real               :: factor ! For unit conversion, default=1, see OVERVIEW above
   real               :: lon_start, lon_end, lat_start, lat_end
   integer            :: region_type   
end type data_type


type override_type
   character(len=3)                 :: gridname  
   character(len=128)               :: fieldname
   integer                          :: t_index                 !index for time interp
   type(horiz_interp_type), pointer :: horz_interp(:) =>NULL() ! index for horizontal spatial interp
   integer                          :: dims(4)                 ! dimensions(x,y,z,t) of the field in filename
   integer                          :: comp_domain(4)          ! istart,iend,jstart,jend for compute domain
   integer                          :: numthreads
   real, _ALLOCATABLE               :: lon_in(:) _NULL
   real, _ALLOCATABLE               :: lat_in(:) _NULL
   logical, _ALLOCATABLE            :: need_compute(:) _NULL
   integer                          :: numwindows
   integer                          :: window_size(2)
   integer                          :: is_src, ie_src, js_src, je_src
end type override_type

 integer, parameter :: max_table=100, max_array=100
 integer            :: table_size ! actual size of data table
 integer, parameter :: ANNUAL=1, MONTHLY=2, DAILY=3, HOURLY=4, UNDEF=-1
 real, parameter    :: tpi=2*PI
 real               :: deg_to_radian, radian_to_deg 
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
integer:: num_fields = 0 ! number of fields in override_array already processed
type(data_type), dimension(max_table)           :: data_table ! user-provided data table
type(data_type)                                 :: default_table
type(override_type), dimension(max_array), save :: override_array ! to store processed fields
type(override_type), save                       :: default_array
logical                                         :: atm_on, ocn_on, lnd_on, ice_on
logical                                         :: lndUG_on
logical                                         :: debug_data_override
logical                                         :: grid_center_bug = .false.

namelist /data_override_nml/ debug_data_override, grid_center_bug

interface data_override
     module procedure data_override_0d
     module procedure data_override_2d
     module procedure data_override_3d
end interface

interface data_override_UG
     module procedure data_override_UG_1d
     module procedure data_override_UG_2d
end interface

public :: data_override_init, data_override, data_override_unset_domains
public :: data_override_UG

contains
!===============================================================================================
! <SUBROUTINE NAME="data_override_init">
!   <DESCRIPTION>
! Assign default values for default_table, get domain of component models,
! get global grids of component models.
! Users should call data_override_init before calling data_override
!   </DESCRIPTION>
!   <TEMPLATE>
! call data_override_init
!   </TEMPLATE>
subroutine data_override_init(Atm_domain_in, Ocean_domain_in, Ice_domain_in, Land_domain_in, Land_domainUG_in)
  type (domain2d), intent(in), optional :: Atm_domain_in
  type (domain2d), intent(in), optional :: Ocean_domain_in, Ice_domain_in 
  type (domain2d), intent(in), optional :: Land_domain_in
  type(domainUG) , intent(in), optional :: Land_domainUG_in

! <NOTE>
! This subroutine should be called in coupler_init after
! (ocean/atmos/land/ice)_model_init have been called.
!
! data_override_init can be called more than once, in one call the user can pass 
! up to 4 domains of component models, at least one domain must be present in
! any call
!
! Data_table is initialized here with default values. Users should provide "real" values
! that will override the default values. Real values can be given using data_table, each
! line of data_table contains one data_entry. Items of data_entry are comma separated.
!
! </NOTE>
  character(len=128)    :: grid_file = 'INPUT/grid_spec.nc'
  integer               :: is,ie,js,je,count
  integer               :: i, iunit, ntable, ntable_lima, ntable_new, unit,io_status, ierr
  character(len=256)    :: record
  logical               :: file_open
  logical               :: ongrid
  character(len=128)    :: region, region_type


  type(data_type)  :: data_entry

  debug_data_override = .false.

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, data_override_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'data_override_nml')
#else
  iunit = open_namelist_file ()
  ierr=1; do while (ierr /= 0)
  read  (iunit, nml=data_override_nml, iostat=io_status, end=10)
  ierr = check_nml_error(io_status, 'data_override_nml')
  enddo
10 call close_file (iunit)
#endif
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
    radian_to_deg = 180./PI
    deg_to_radian = PI/180.

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
    call mpp_open(iunit, 'data_table', action=MPP_RDONLY)
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
    call mpp_close(iunit)
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
 endif

 module_is_initialized = .TRUE.

 if ( .NOT. (atm_on .or. ocn_on .or. lnd_on .or. ice_on .or. lndUG_on)) return 
 call fms_io_init

! Test if grid_file is already opened
 inquire (file=trim(grid_file), opened=file_open)
 if(file_open) call mpp_error(FATAL, trim(grid_file)//' already opened')

 if(field_exist(grid_file, "x_T" ) .OR. field_exist(grid_file, "geolon_t" ) ) then
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
 else if(field_exist(grid_file, "ocn_mosaic_file" ) .OR. field_exist(grid_file, "gridfiles" ) ) then
    if(field_exist(grid_file, "gridfiles" ) ) then
       count = 0
       if (atm_on) count = count + 1 
       if (lnd_on) count = count + 1 
       if ( ocn_on .OR. ice_on ) count = count + 1 
       if(count .NE. 1) call mpp_error(FATAL, 'data_override_mod: the grid file is a solo mosaic, ' // &
            'one and only one of atm_on, lnd_on or ice_on/ocn_on should be true')
    endif
   if (atm_on .and. .not. allocated(lon_local_atm) ) then
       call mpp_get_compute_domain(atm_domain,is,ie,js,je) 
       allocate(lon_local_atm(is:ie,js:je), lat_local_atm(is:ie,js:je))
       call get_grid_version_2(grid_file, 'atm', atm_domain, is, ie, js, je, lon_local_atm, lat_local_atm, &
            min_glo_lon_atm, max_glo_lon_atm )
    endif

    if (ocn_on .and. .not. allocated(lon_local_ocn) ) then
       call mpp_get_compute_domain( ocn_domain,is,ie,js,je) 
       allocate(lon_local_ocn(is:ie,js:je), lat_local_ocn(is:ie,js:je))
       call get_grid_version_2(grid_file, 'ocn', ocn_domain, is, ie, js, je, lon_local_ocn, lat_local_ocn, &
            min_glo_lon_ocn, max_glo_lon_ocn )
    endif

    if (lnd_on .and. .not. allocated(lon_local_lnd) ) then
       call mpp_get_compute_domain( lnd_domain,is,ie,js,je) 
       allocate(lon_local_lnd(is:ie,js:je), lat_local_lnd(is:ie,js:je))
       call get_grid_version_2(grid_file, 'lnd', lnd_domain, is, ie, js, je, lon_local_lnd, lat_local_lnd, &
            min_glo_lon_lnd, max_glo_lon_lnd )
    endif

    if (ice_on .and. .not. allocated(lon_local_ice) ) then
       call mpp_get_compute_domain( ice_domain,is,ie,js,je) 
       allocate(lon_local_ice(is:ie,js:je), lat_local_ice(is:ie,js:je))
       call get_grid_version_2(grid_file, 'ocn', ice_domain, is, ie, js, je, lon_local_ice, lat_local_ice, &
            min_glo_lon_ice, max_glo_lon_ice )
    endif
 else
    call mpp_error(FATAL, 'data_override_mod: none of x_T, geolon_t, ocn_mosaic_file or gridfiles exist in '//trim(grid_file))
 end if

end subroutine data_override_init
! </SUBROUTINE>
!===============================================================================================

!===============================================================================================
! <SUBROUTINE NAME="data_override_unset_domain">
!   <DESCRIPTION>
! Unset domains that had previously been set for use by data_override.
!   </DESCRIPTION>
!   <TEMPLATE>
! call data_override_unset_domain
!   </TEMPLATE>
subroutine data_override_unset_domains(unset_Atm, unset_Ocean, &
                                      unset_Ice, unset_Land, must_be_set)
  logical, intent(in), optional :: unset_Atm, unset_Ocean, unset_Ice, unset_Land
  logical, intent(in), optional :: must_be_set

! <NOTE>
! This subroutine deallocates any data override domains that have been set.
! </NOTE>
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
! </SUBROUTINE>
!===============================================================================================


subroutine check_grid_sizes(domain_name, Domain, nlon, nlat)
character(len=12), intent(in) :: domain_name
type (domain2d),   intent(in) :: Domain
integer,           intent(in) :: nlon, nlat

character(len=184) :: error_message
integer            :: xsize, ysize

call mpp_get_global_domain(Domain, xsize=xsize, ysize=ysize)
if(nlon .NE. xsize .OR. nlat .NE. ysize) then
  error_message = 'Error in data_override_init. Size of grid as specified by '// &
                  '             does not conform to that specified by grid_spec.nc.'// &
                  '  From             :     by      From grid_spec.nc:     by    '
  error_message( 59: 70) = domain_name
  error_message(130:141) = domain_name
  write(error_message(143:146),'(i4)') xsize
  write(error_message(150:153),'(i4)') ysize
  write(error_message(174:177),'(i4)') nlon
  write(error_message(181:184),'(i4)') nlat
  call mpp_error(FATAL,error_message)
endif

end subroutine check_grid_sizes
!===============================================================================================
subroutine get_domain(gridname, domain, comp_domain)
! Given a gridname, this routine returns the working domain associated with this gridname
  character(len=3), intent(in) :: gridname
  type(domain2D), intent(inout) :: domain
  integer, intent(out), optional :: comp_domain(4) ! istart,iend,jstart,jend for compute domain

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

subroutine get_domainUG(gridname, UGdomain, comp_domain)
! Given a gridname, this routine returns the working domain associated with this gridname
  character(len=3), intent(in) :: gridname
  type(domainUG), intent(inout) :: UGdomain
  integer, intent(out), optional :: comp_domain(4) ! istart,iend,jstart,jend for compute domain
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

! <SUBROUTINE NAME="data_override_2d">
!   <DESCRIPTION>
! This routine performs data override for 2D fields; for usage, see data_override_3d.
!   </DESCRIPTION>
subroutine data_override_2d(gridname,fieldname,data_2D,time,override, is_in, ie_in, js_in, je_in)
  character(len=3), intent(in) :: gridname ! model grid ID
  character(len=*), intent(in) :: fieldname ! field to override
  logical, intent(out), optional :: override ! true if the field has been overriden succesfully
  type(time_type), intent(in) :: time !  model time
  real, dimension(:,:), intent(inout) :: data_2D !data returned by this call
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
! </SUBROUTINE>
!===============================================================================================

! <SUBROUTINE NAME="data_override_3d">
!   <DESCRIPTION>
! This routine performs data override for 3D fields
!   <TEMPLATE>
! call data_override(gridname,fieldname,data,time,override)
!   </TEMPLATE>
!   </DESCRIPTION>

!   <IN NAME="gridname"  TYPE="character" DIM="(*)">
! Grid name (Ocean, Ice, Atmosphere, Land)
!   </IN>
!   <IN NAME="fieldname_code" TYPE="character" DIM="(*)">
!    Field name as used in the code (may be different from the name in NetCDF data file)
!   </IN>
!   <OUT NAME="data" TYPE="real" DIM="(:,:,:)">
!    array containing output data
!   </OUT>
!   <IN NAME="time" TYPE="time_type">
!    model time
!   </IN>
!   <OUT NAME="override" TYPE="logical">
!    TRUE if the field is overriden, FALSE otherwise
!   </OUT>
!   <IN NAME="data_index" TYPE="integer">
!   </IN>
subroutine data_override_3d(gridname,fieldname_code,data,time,override,data_index, is_in, ie_in, js_in, je_in)
  character(len=3),             intent(in) :: gridname ! model grid ID
  character(len=*),             intent(in) :: fieldname_code ! field name as used in the model
  logical,           optional, intent(out) :: override ! true if the field has been overriden succesfully
  type(time_type),              intent(in) :: time !(target) model time
  integer,           optional,  intent(in) :: data_index
  real, dimension(:,:,:),    intent(inout) :: data !data returned by this call
  integer,           optional,  intent(in) :: is_in, ie_in, js_in, je_in
  logical, dimension(:,:,:),   allocatable :: mask_out

  character(len=512) :: filename, filename2 !file containing source data
  character(len=128) :: fieldname ! fieldname used in the data file
  integer            :: i,j
  integer            :: dims(4)
  integer            :: index1 ! field index in data_table
  integer            :: id_time !index for time interp in override array
  integer            :: axis_sizes(4)
  real, dimension(:,:), pointer :: lon_local =>NULL(), &
                                   lat_local =>NULL() !of output (target) grid cells
  real, dimension(:), allocatable :: lon_tmp, lat_tmp

  type(axistype) :: axis_centers(4), axis_bounds(4)
  logical :: data_file_is_2D = .false.  !data in netCDF file is 2D
  logical :: ongrid, use_comp_domain
  type(domain2D) :: domain
  integer :: curr_position ! position of the field currently processed in override_array
  real :: factor
  integer, dimension(4) :: comp_domain = 0  ! istart,iend,jstart,jend for compute domain
  integer :: nxd, nyd, nxc, nyc, nwindows
  integer :: nwindows_x, ipos, jpos, window_size(2)
  integer :: istart, iend, jstart, jend
  integer :: isw, iew, jsw, jew, n
  integer :: omp_get_num_threads, omp_get_thread_num, thread_id, window_id
  logical :: need_compute
  real    :: lat_min, lat_max
  integer :: is_src, ie_src, js_src, je_src
  logical :: exists

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
       if(mpp_pe() == mpp_root_pe() .and. debug_data_override) &
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
                                      use_comp_domain=use_comp_domain, nwindows=nwindows)
        dims = get_external_field_size(id_time)
        override_array(curr_position)%dims = dims
        if(id_time<0) call mpp_error(FATAL,'data_override:field not found in init_external_field 1') 
        override_array(curr_position)%t_index = id_time     
     else !ongrid=false
        id_time = init_external_field(filename,fieldname,domain=domain, axis_centers=axis_centers,&
             axis_sizes=axis_sizes, verbose=.false.,override=.true.,use_comp_domain=use_comp_domain, &
             nwindows = nwindows)  
        dims = get_external_field_size(id_time)
        override_array(curr_position)%dims = dims
        if(id_time<0) call mpp_error(FATAL,'data_override:field not found in init_external_field 2')
        override_array(curr_position)%t_index = id_time

        !  get lon and lat of the input (source) grid, assuming that axis%data contains
        !  lat and lon of the input grid (in degrees)
        call get_axis_bounds(axis_centers(1),axis_bounds(1), axis_centers)
        call get_axis_bounds(axis_centers(2),axis_bounds(2), axis_centers)

        allocate(override_array(curr_position)%horz_interp(nwindows))
        allocate(override_array(curr_position)%lon_in(axis_sizes(1)+1))
        allocate(override_array(curr_position)%lat_in(axis_sizes(2)+1))
        call mpp_get_axis_data(axis_bounds(1),override_array(curr_position)%lon_in)
        call mpp_get_axis_data(axis_bounds(2),override_array(curr_position)%lat_in)

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
           call mpp_get_axis_data(axis_centers(1), lon_tmp)
           call mpp_get_axis_data(axis_centers(2), lat_tmp)
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

  if(ongrid) then
!10 do time interp to get data in compute_domain
     if(data_file_is_2D) then        
        call time_interp_external(id_time,time,data(:,:,1),verbose=.false., &
                                  is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
        data(:,:,1) = data(:,:,1)*factor
        do i = 2, size(data,3)
           data(:,:,i) = data(:,:,1)
        enddo
     else
        call time_interp_external(id_time,time,data,verbose=.false., &
                                  is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
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
! </SUBROUTINE>

! <SUBROUTINE NAME="data_override_0d">
!   <DESCRIPTION>
! This routine performs data override for scalar fields
!   <TEMPLATE>
! call data_override(fieldname,data,time,override)
!   </TEMPLATE>
!   </DESCRIPTION>
!   <IN NAME="gridname"  TYPE="character" DIM="(*)">
! Grid name (Ocean, Ice, Atmosphere, Land)
!   </IN>
!   <IN NAME="fieldname_code" TYPE="character" DIM="(*)">
!    Field name as used in the code (may be different from the name in NetCDF data file)
!   </IN>
!   <OUT NAME="data" TYPE="real" DIM="(:,:,:)">
!    array containing output data
!   </OUT>
!   <IN NAME="time" TYPE="time_type">
!    model time
!   </IN>
!   <OUT NAME="override" TYPE="logical">
!    TRUE if the field is overriden, FALSE otherwise
!   </OUT>
!   <IN NAME="data_index" TYPE="integer">
!   </IN>
subroutine data_override_0d(gridname,fieldname_code,data,time,override,data_index)
  character(len=3), intent(in) :: gridname ! model grid ID
  character(len=*), intent(in) :: fieldname_code ! field name as used in the model
  logical, intent(out), optional :: override ! true if the field has been overriden succesfully
  type(time_type), intent(in) :: time !(target) model time
  real,             intent(out) :: data !data returned by this call
  integer, intent(in), optional :: data_index

  character(len=512) :: filename !file containing source data
  character(len=128) :: fieldname ! fieldname used in the data file
  integer :: index1 ! field index in data_table
  integer :: id_time !index for time interp in override array
  integer :: curr_position ! position of the field currently processed in override_array
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
       if(mpp_pe() == mpp_root_pe() .and. debug_data_override) &
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
! </SUBROUTINE>

subroutine data_override_UG_1d(gridname,fieldname,data,time,override)
  character(len=3),   intent(in) :: gridname ! model grid ID
  character(len=*),   intent(in) :: fieldname ! field to override
  real, dimension(:), intent(inout) :: data !data returned by this call
  type(time_type),    intent(in) :: time !  model time
  logical, intent(out), optional :: override ! true if the field has been overriden succesfully
  !local vars
  real, dimension(:,:), allocatable ::  data_SG
  type(domainUG) :: UG_domain
  integer       :: index1
  integer       :: i
  integer, dimension(4) :: comp_domain = 0  ! istart,iend,jstart,jend for compute domain

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

subroutine data_override_UG_2d(gridname,fieldname,data,time,override)
  character(len=3),     intent(in) :: gridname ! model grid ID
  character(len=*),     intent(in) :: fieldname ! field to override
  real, dimension(:,:), intent(inout) :: data !data returned by this call
  type(time_type),      intent(in) :: time !  model time
  logical, intent(out), optional :: override ! true if the field has been overriden succesfully
  !local vars
  real, dimension(:,:,:), allocatable ::  data_SG
  real, dimension(:,:),   allocatable ::  data_UG
  type(domainUG) :: UG_domain
  integer       :: index1
  integer       :: i, nlevel, nlevel_max
  integer, dimension(4) :: comp_domain = 0  ! istart,iend,jstart,jend for compute domain

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
  call data_override_3d(gridname,fieldname,data_SG(:,:,1:nlevel),time,override)    

  call mpp_pass_SG_to_UG(UG_domain, data_SG(:,:,:), data_UG(:,:))
  data(:,1:nlevel) = data_UG(:,1:nlevel)  

  deallocate(data_SG, data_UG)

end subroutine data_override_UG_2d

!===============================================================================================

! Get lon and lat of three model (target) grids from grid_spec.nc
subroutine get_grid_version_1(grid_file, mod_name, domain, isc, iec, jsc, jec, lon, lat, min_lon, max_lon)
  character(len=*),            intent(in) :: grid_file
  character(len=*),            intent(in) :: mod_name
  type(domain2d),              intent(in) :: domain
  integer,                     intent(in) :: isc, iec, jsc, jec
  real, dimension(isc:,jsc:), intent(out) :: lon, lat
  real,                       intent(out) :: min_lon, max_lon

  integer                                      :: i, j, siz(4)
  integer                                      :: nlon, nlat         ! size of global lon and lat
  real,          dimension(:,:,:), allocatable :: lon_vert, lat_vert !of OCN grid vertices
  real,          dimension(:),     allocatable :: glon, glat         ! lon and lat of 1-D grid of atm/lnd
  logical                                      :: is_new_grid
  integer                                      :: is, ie, js, je
  integer                                      :: isd, ied, jsd, jed
  integer                                      :: isg, ieg, jsg, jeg
  type(domain2d)                               :: domain2
  character(len=3)                             :: xname, yname

  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)

  select case(mod_name)
  case('ocn', 'ice')
    is_new_grid = .FALSE.
    if(field_exist(grid_file, 'x_T')) then
       is_new_grid = .true.
    else if(field_exist(grid_file, 'geolon_t')) then
       is_new_grid = .FALSE.
    else
       call mpp_error(FATAL,'data_override: both x_T and geolon_t is not in the grid file '//trim(grid_file) )
    endif

    if(is_new_grid) then
      call field_size(grid_file, 'x_T', siz)
      nlon = siz(1); nlat = siz(2)
      call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)
      allocate(lon_vert(isc:iec,jsc:jec,4), lat_vert(isc:iec,jsc:jec,4) )
      call read_data(trim(grid_file), 'x_vert_T', lon_vert, domain)
      call read_data(trim(grid_file), 'y_vert_T', lat_vert, domain)
      
!2 Global lon and lat of ocean grid cell centers are determined from adjacent vertices
      lon(:,:) = (lon_vert(:,:,1) + lon_vert(:,:,2) + lon_vert(:,:,3) + lon_vert(:,:,4))*0.25
      lat(:,:) = (lat_vert(:,:,1) + lat_vert(:,:,2) + lat_vert(:,:,3) + lat_vert(:,:,4))*0.25
    else      
      if(grid_center_bug) call mpp_error(NOTE, &
           'data_override: grid_center_bug is set to true, the grid center location may be incorrect')
      call field_size(grid_file, 'geolon_vert_t', siz)
      nlon = siz(1) - 1; nlat = siz(2) - 1;
      call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)
      call mpp_copy_domain(domain, domain2)
      call mpp_set_compute_domain(domain2, isc, iec+1, jsc, jec+1, iec-isc+2, jec-jsc+2 )
      call mpp_set_data_domain   (domain2, isd, ied+1, jsd, jed+1, ied-isd+2, jed-jsd+2 )   
      call mpp_set_global_domain (domain2, isg, ieg+1, jsg, jeg+1, ieg-isg+2, jeg-jsg+2 )    
      allocate(lon_vert(isc:iec+1,jsc:jec+1,1))
      allocate(lat_vert(isc:iec+1,jsc:jec+1,1))
      call read_data(trim(grid_file), 'geolon_vert_t', lon_vert, domain2)
      call read_data(trim(grid_file), 'geolat_vert_t', lat_vert, domain2)

      if(grid_center_bug) then
         do j = jsc, jec
            do i = isc, iec
               lon(i,j) = (lon_vert(i,j,1) + lon_vert(i+1,j,1))/2.
               lat(i,j) = (lat_vert(i,j,1) + lat_vert(i,j+1,1))/2.
            enddo
         enddo
      else
         do j = jsc, jec
            do i = isc, iec
               lon(i,j) = (lon_vert(i,j,1) + lon_vert(i+1,j,1) + &
                    lon_vert(i+1,j+1,1) + lon_vert(i,j+1,1))*0.25
               lat(i,j) = (lat_vert(i,j,1) + lat_vert(i+1,j,1) + &
                    lat_vert(i+1,j+1,1) + lat_vert(i,j+1,1))*0.25
            enddo
         enddo
      end if
      call mpp_deallocate_domain(domain2)
    endif
    deallocate(lon_vert)
    deallocate(lat_vert)
  case('atm', 'lnd')
     if(trim(mod_name) == 'atm') then
        xname = 'xta'; yname = 'yta'
     else
        xname = 'xtl'; yname = 'ytl'
     endif
     call field_size(grid_file, xname, siz)
     nlon = siz(1); allocate(glon(nlon))
     call read_data(grid_file, xname, glon, no_domain = .true.)

     call field_size(grid_file, yname, siz)
     nlat = siz(1); allocate(glat(nlat))
     call read_data(grid_file, yname, glat, no_domain = .true.)
     call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)

     is = isc - isg + 1; ie = iec - isg + 1
     js = jsc - jsg + 1; je = jec - jsg + 1
     do j = js, jec
        do i = is, ie
           lon(i,j) = glon(i)
           lat(i,j) = glat(j)
        enddo
     enddo
     deallocate(glon)
     deallocate(glat)
  case default
     call mpp_error(FATAL, "data_override_mod: mod_name should be 'atm', 'ocn', 'ice' or 'lnd' ")
  end select
 
  ! convert from degree to radian
  lon = lon * deg_to_radian
  lat = lat* deg_to_radian
  min_lon = minval(lon)
  max_lon = maxval(lon)
  call mpp_min(min_lon)
  call mpp_max(max_lon)


end subroutine get_grid_version_1

! Get global lon and lat of three model (target) grids from mosaic.nc
! z1l: currently we assume the refinement ratio is 2 and there is one tile on each pe.
subroutine get_grid_version_2(mosaic_file, mod_name, domain, isc, iec, jsc, jec, lon, lat, min_lon, max_lon)
  character(len=*),            intent(in) :: mosaic_file
  character(len=*),            intent(in) :: mod_name
  type(domain2d),              intent(in) :: domain
  integer,                     intent(in) :: isc, iec, jsc, jec
  real, dimension(isc:,jsc:), intent(out) :: lon, lat
  real,                       intent(out) :: min_lon, max_lon

  integer            :: i, j, siz(4)
  integer            :: nlon, nlat             ! size of global grid
  integer            :: nlon_super, nlat_super ! size of global supergrid.
  integer            :: isd, ied, jsd, jed
  integer            :: isg, ieg, jsg, jeg
  integer            :: isc2, iec2, jsc2, jec2
  character(len=256) :: solo_mosaic_file, grid_file
  real, allocatable  :: tmpx(:,:), tmpy(:,:)
  type(domain2d)     :: domain2

  if(trim(mod_name) .NE. 'atm' .AND. trim(mod_name) .NE. 'ocn' .AND. &
     trim(mod_name) .NE. 'ice' .AND. trim(mod_name) .NE. 'lnd' ) call mpp_error(FATAL, &
        "data_override_mod: mod_name should be 'atm', 'ocn', 'ice' or 'lnd' ")

  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)

  ! get the grid file to read
  if(field_exist(mosaic_file, trim(mod_name)//'_mosaic_file' )) then
     call read_data(mosaic_file, trim(mod_name)//'_mosaic_file', solo_mosaic_file) 
     solo_mosaic_file = 'INPUT/'//trim(solo_mosaic_file)
  else
     solo_mosaic_file = mosaic_file
  end if
  call get_mosaic_tile_grid(grid_file, solo_mosaic_file, domain)

  call field_size(grid_file, 'area', siz)
  nlon_super = siz(1); nlat_super = siz(2)
  if( mod(nlon_super,2) .NE. 0) call mpp_error(FATAL,  &
       'data_override_mod: '//trim(mod_name)//' supergrid longitude size can not be divided by 2')
  if( mod(nlat_super,2) .NE. 0) call mpp_error(FATAL,  &
       'data_override_mod: '//trim(mod_name)//' supergrid latitude size can not be divided by 2')
  nlon = nlon_super/2;
  nlat = nlat_super/2;     
  call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)

  !--- setup the domain for super grid.
  call mpp_copy_domain(domain, domain2)
  call mpp_set_compute_domain(domain2, 2*isc-1, 2*iec+1, 2*jsc-1, 2*jec+1, 2*iec-2*isc+3, 2*jec-2*jsc+3 )
  call mpp_set_data_domain   (domain2, 2*isd-1, 2*ied+1, 2*jsd-1, 2*jed+1, 2*ied-2*isd+3, 2*jed-2*jsd+3 )   
  call mpp_set_global_domain (domain2, 2*isg-1, 2*ieg+1, 2*jsg-1, 2*jeg+1, 2*ieg-2*isg+3, 2*jeg-2*jsg+3 )

  call mpp_get_compute_domain(domain2, isc2, iec2, jsc2, jec2)
  if(isc2 .NE. 2*isc-1 .OR. iec2 .NE. 2*iec+1 .OR. jsc2 .NE. 2*jsc-1 .OR. jec2 .NE. 2*jec+1) then
     call mpp_error(FATAL, 'data_override_mod: '//trim(mod_name)//' supergrid domain is not set properly')
  endif

  allocate(tmpx(isc2:iec2, jsc2:jec2), tmpy(isc2:iec2, jsc2:jec2) )
  call read_data( grid_file, 'x', tmpx, domain2)
  call read_data( grid_file, 'y', tmpy, domain2)     
  ! copy data onto model grid
  if(trim(mod_name) == 'ocn' .OR. trim(mod_name) == 'ice') then
     do j = jsc, jec
        do i = isc, iec
           lon(i,j) = (tmpx(i*2-1,j*2-1)+tmpx(i*2+1,j*2-1)+tmpx(i*2+1,j*2+1)+tmpx(i*2-1,j*2+1))*0.25
           lat(i,j) = (tmpy(i*2-1,j*2-1)+tmpy(i*2+1,j*2-1)+tmpy(i*2+1,j*2+1)+tmpy(i*2-1,j*2+1))*0.25
        end do
     end do
  else
     do j = jsc, jec
        do i = isc, iec
           lon(i,j) = tmpx(i*2,j*2)
           lat(i,j) = tmpy(i*2,j*2)
        end do
     end do
  endif

  ! convert to radian
  lon = lon * deg_to_radian
  lat = lat * deg_to_radian

  deallocate(tmpx, tmpy)
  min_lon = minval(lon)
  max_lon = maxval(lon)
  call mpp_min(min_lon)
  call mpp_max(max_lon)

  call mpp_deallocate_domain(domain2)

end subroutine get_grid_version_2

!===============================================================================================
end module data_override_mod

#ifdef test_data_override

 program test

 ! Input data and path_names file for this program is in:
 ! /archive/pjp/unit_tests/test_data_override/lima/exp1

 use           mpp_mod, only: input_nml_file, stdout, mpp_chksum
 use   mpp_domains_mod, only: domain2d, mpp_define_domains, mpp_get_compute_domain, mpp_define_layout
 use           fms_mod, only: fms_init, fms_end, mpp_npes, file_exist, open_namelist_file, check_nml_error, close_file
 use           fms_mod, only: error_mesg, FATAL, file_exist, field_exist, field_size
 use        fms_io_mod, only: read_data, fms_io_exit
 use     constants_mod, only: constants_init, pi
 use  time_manager_mod, only: time_type, set_calendar_type, set_date, NOLEAP, JULIAN, operator(+), set_time, print_time
 use  diag_manager_mod, only: diag_manager_init, diag_manager_end, register_static_field, register_diag_field
 use  diag_manager_mod, only: send_data, diag_axis_init
 use data_override_mod, only: data_override_init, data_override, data_override_UG
  use mpp_mod,         only : FATAL, WARNING, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_node, mpp_root_pe, mpp_error, mpp_set_warn_level
  use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync, mpp_sync_self
  use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod,         only : mpp_init, mpp_exit, mpp_chksum, stdout, stderr
  use mpp_mod,         only : input_nml_file
  use mpp_mod,         only : mpp_get_current_pelist, mpp_broadcast
  use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, CGRID_NE, DGRID_NE
  use mpp_domains_mod, only : FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE, FOLD_WEST_EDGE, FOLD_EAST_EDGE
  use mpp_domains_mod, only : MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
  use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D, BITWISE_EFP_SUM
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod, only : mpp_update_domains, mpp_check_field, mpp_redistribute, mpp_get_memory_domain
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains, mpp_modify_domain
  use mpp_domains_mod, only : mpp_get_neighbor_pe, mpp_define_mosaic, mpp_nullify_domain_list
  use mpp_domains_mod, only : NORTH, NORTH_EAST, EAST, SOUTH_EAST, CORNER, CENTER
  use mpp_domains_mod, only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST, mpp_define_mosaic_pelist
  use mpp_domains_mod, only : mpp_get_global_domain, ZERO, NINETY, MINUS_NINETY
  use mpp_domains_mod, only : mpp_get_boundary, mpp_start_update_domains, mpp_complete_update_domains
  use mpp_domains_mod, only : mpp_define_nest_domains, nest_domain_type
  use mpp_domains_mod, only : mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod, only : mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod, only : mpp_get_domain_shift, EDGEUPDATE, mpp_deallocate_domain
  use mpp_domains_mod, only : mpp_group_update_type, mpp_create_group_update
  use mpp_domains_mod, only : mpp_do_group_update, mpp_clear_group_update
  use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
  use mpp_domains_mod, only : WUPDATE, SUPDATE, mpp_get_compute_domains
  use mpp_domains_mod, only : domainUG, mpp_define_unstruct_domain, mpp_get_UG_domain_tile_id
  use mpp_domains_mod, only : mpp_get_UG_compute_domain, mpp_pass_SG_to_UG, mpp_pass_UG_to_SG
  use mpp_domains_mod, only : mpp_get_ug_global_domain, mpp_global_field_ug
  use mpp_memutils_mod, only : mpp_memuse_begin, mpp_memuse_end

 implicit none

 integer                           :: stdoutunit
 integer                           :: num_threads = 1
 integer                           :: omp_get_thread_num
 integer                           :: isw, iew, jsw, jew
 integer, allocatable              :: is_win(:), js_win(:)
 integer                           :: nx_dom, ny_dom, nx_win, ny_win
 type(domain2d)                    :: Domain
 integer                           :: nlon, nlat, siz(4)
 real, allocatable, dimension(:)   :: x, y
 real, allocatable, dimension(:,:) :: lon, lat
 real, allocatable, dimension(:,:) :: sst, ice
 integer                           :: id_x, id_y, id_lon, id_lat, id_sst, id_ice
 integer                           :: i, j, is, ie, js, je, unit, io, ierr, n
 real                              :: rad_to_deg
 character(len=36)                 :: message
 type(time_type)                   :: Time
 logical                           :: used
 logical, allocatable              :: ov_sst(:), ov_ice(:)
 integer, dimension(2)             :: layout = (/0,0/)
 character(len=256)                :: solo_mosaic_file, tile_file
 character(len=128)                :: grid_file   = "INPUT/grid_spec.nc"
 integer                           :: window(2) = (/1,1/)
 integer                           :: get_cpu_affinity, base_cpu
 integer                           :: nthreads=1
 integer                           :: nwindows

 namelist / test_data_override_nml / layout, window, nthreads

 call fms_init
 call constants_init
 call set_calendar_type(NOLEAP)
 call diag_manager_init

 rad_to_deg = 180./pi

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, test_data_override_nml, iostat=io)
      ierr = check_nml_error(io, 'test_data_override_nml')
#else
 if (file_exist('input.nml')) then
   unit = open_namelist_file ( )
   ierr=1
   do while (ierr /= 0)
     read(unit, nml=test_data_override_nml, iostat=io, end=10)
          ierr = check_nml_error(io, 'test_data_override_nml')
   enddo
10 call close_file (unit)
 endif
#endif

 if(field_exist(grid_file, "x_T" ) ) then
    call field_size(grid_file, 'x_T', siz)
    nlon = siz(1)
    nlat = siz(2)
 else if(field_exist(grid_file, "geolon_t" ) ) then
    call field_size(grid_file, 'geolon_t', siz)
    nlon = siz(1)
    nlat = siz(2)
 else if (field_exist(grid_file, "ocn_mosaic_file" )) then
    call read_data(grid_file, 'ocn_mosaic_file', solo_mosaic_file) 
    solo_mosaic_file = 'INPUT/'//trim(solo_mosaic_file)
    call field_size(solo_mosaic_file, 'gridfiles', siz)
    if( siz(2) .NE. 1) call error_mesg('test_data_override', 'only support single tile mosaic, contact developer', FATAL)
    call read_data(solo_mosaic_file, 'gridfiles', tile_file)
    tile_file = 'INPUT/'//trim(tile_file)
    call field_size(tile_file, 'area', siz)
    if(mod(siz(1),2) .NE. 0 .OR. mod(siz(2),2) .NE. 0 ) call error_mesg('test_data_override', &
        "test_data_override: supergrid size can not be divided by 2", FATAL)
    nlon = siz(1)/2
    nlat = siz(2)/2
 else
    call error_mesg('test_data_override', 'x_T, geolon_t and ocn_mosaic_file does not exist', FATAL)
 end if

 if(layout(1)*layout(2) .NE. mpp_npes() ) then
    call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes(), layout )
 end if


 call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, name='test_data_override')
 call data_override_init(Ice_domain_in=Domain, Ocean_domain_in=Domain)
 call data_override_init(Ice_domain_in=Domain, Ocean_domain_in=Domain)
 call mpp_get_compute_domain(Domain, is, ie, js, je)
 call get_grid

 allocate(x(nlon), y(nlat))

 do i=1,nlon
   x(i) = i
 enddo
 do j=1,nlat
   y(j) = j
 enddo

 Time = set_date(2,1,1,0,0,0)

 allocate(sst(is:ie,js:je), ice(is:ie,js:je))
 sst = 0
 ice = 0

 id_x  = diag_axis_init('x',  x,  'point_E', 'x', long_name='point_E', Domain2=Domain)
 id_y  = diag_axis_init('y',  y,  'point_N', 'y', long_name='point_N', Domain2=Domain)

 Time = Time + set_time(0,1)

 id_lon = register_static_field('test_data_override_mod', 'lon', (/id_x,id_y/), 'longitude', 'Degrees')
 id_lat = register_static_field('test_data_override_mod', 'lat', (/id_x,id_y/), 'longitude', 'Degrees')
 id_sst = register_diag_field  ('test_data_override_mod', 'sst', (/id_x,id_y/), Time, 'SST', 'K')
 id_ice = register_diag_field  ('test_data_override_mod', 'ice', (/id_x,id_y/), Time, 'ICE', ' ')
 used = send_data(id_lon, lon, Time)
 used = send_data(id_lat, lat, Time)

 sst = 0.
 ice = 0.

! get number of threads

nx_dom = ie - is + 1
ny_dom = je - js + 1
if( mod( nx_dom, window(1) ) .NE. 0 ) call error_mesg('test_data_override', &
        "nx_dom is not divisible by window(1)", FATAL)
if( mod( ny_dom, window(2) ) .NE. 0 ) call error_mesg('test_data_override', &
        "ny_dom is not divisible by window(2)", FATAL)

nwindows = window(1)*window(2)
!$ call omp_set_num_threads(nthreads)
!$ base_cpu = get_cpu_affinity()
!$OMP PARALLEL
!$ call set_cpu_affinity( base_cpu + omp_get_thread_num() )
!$OMP END PARALLEL

nx_win = nx_dom/window(1)
ny_win = ny_dom/window(2)
allocate(is_win(nwindows), js_win(nwindows))

i = 1
do jsw = js,je,ny_win
   do isw = is,ie,nx_win
      is_win(i) = isw
      js_win(i) = jsw
      i = i + 1
   enddo
enddo

allocate(ov_sst(nwindows), ov_ice(nwindows))
!$OMP parallel  do schedule(static) default(shared) private(isw, iew, jsw, jew)
do n = 1, nwindows
   isw = is_win(n)
   iew = isw + nx_win - 1
   jsw = js_win(n)
   jew = jsw + ny_win - 1
   call data_override('OCN','sst_obs',sst(isw:iew,jsw:jew),Time,override=ov_sst(n), &
                      is_in=isw-is+1, ie_in=iew-is+1, js_in=jsw-js+1, je_in=jew-js+1)
   call data_override('ICE', 'sic_obs', ice(isw:iew,jsw:jew), Time, override=ov_ice(n), &
                      is_in=isw-is+1, ie_in=iew-is+1, js_in=jsw-js+1, je_in=jew-js+1)
enddo

 if(ANY(.NOT. ov_sst) .or. ANY(.not.ov_ice)) then
   if(ANY(.NOT. ov_sst)) then
     message = 'override failed for ice'
   else if(ANY(.NOT. ov_ice)) then
     message = 'override failed for sst'
   else
     message = 'override failed for both sst and ice'
   endif
   call error_mesg('test_data_override', trim(message), FATAL)
 endif

 stdoutunit = stdout()
 write(stdoutunit,*)"===>NOTE from test_data_override: sst chksum = ", mpp_chksum(sst)
 write(stdoutunit,*)"===>NOTE from test_data_override: ice chksum = ", mpp_chksum(ice)


 if(id_sst > 0) used = send_data(id_sst, sst, Time)
 if(id_ice > 0) used = send_data(id_ice, ice, Time)

 call test_unstruct_grid( 'Cubic-Grid', Time )



!-------------------------------------------------------------------------------------------------------
! What follows is a test of calendar conversion
 
!Time = set_date(1980,2,27,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)
 
!Time = set_date(1980,2,28,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)
 
!Time = set_date(1980,2,29,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)
 
!Time = set_date(1980,3,1,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)
 
!Time = set_date(1980,3,2,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)
!-------------------------------------------------------------------------------------------------------

 call diag_manager_end(Time)
 call fms_io_exit
 call fms_end

contains 

!=================================================================================================================================
 subroutine get_grid
   real, allocatable, dimension(:,:,:) :: lon_vert_glo, lat_vert_glo
   real, allocatable, dimension(:,:)   :: lon_global, lat_global
   integer, dimension(4)  :: siz
   character(len=128) :: message


   if(field_exist(grid_file, 'x_T')) then
      call field_size(grid_file, 'x_T', siz)
      if(siz(1) /= nlon .or. siz(2) /= nlat) then
         write(message,'(a,2i4)') 'x_T is wrong shape. shape(x_T)=',siz(1:2)
         call error_mesg('test_data_override', trim(message), FATAL)
      endif
      allocate(lon_vert_glo(nlon,nlat,4), lat_vert_glo(nlon,nlat,4) )
      allocate(lon_global  (nlon,nlat  ), lat_global  (nlon,nlat  ) )
      call read_data(trim(grid_file), 'x_vert_T', lon_vert_glo, no_domain=.true.)
      call read_data(trim(grid_file), 'y_vert_T', lat_vert_glo, no_domain=.true.)
      lon_global(:,:)  = (lon_vert_glo(:,:,1) + lon_vert_glo(:,:,2) + lon_vert_glo(:,:,3) + lon_vert_glo(:,:,4))*0.25
      lat_global(:,:) =  (lat_vert_glo(:,:,1) + lat_vert_glo(:,:,2) + lat_vert_glo(:,:,3) + lat_vert_glo(:,:,4))*0.25
   else  if(field_exist(grid_file, "geolon_t" ) ) then
      call field_size(grid_file, 'geolon_vert_t', siz)
      if(siz(1) /= nlon+1 .or. siz(2) /= nlat+1) then
         write(message,'(a,2i4)') 'geolon_vert_t is wrong shape. shape(geolon_vert_t)=',siz(1:2)
         call error_mesg('test_data_override', trim(message), FATAL)
      endif
      allocate(lon_vert_glo(nlon+1,nlat+1,1), lat_vert_glo(nlon+1,nlat+1,1))
      allocate(lon_global  (nlon,  nlat    ), lat_global  (nlon,  nlat    ))
      call read_data(trim(grid_file), 'geolon_vert_t', lon_vert_glo, no_domain=.true.)
      call read_data(trim(grid_file), 'geolat_vert_t', lat_vert_glo, no_domain=.true.)

      do i = 1, nlon
         do j = 1, nlat
            lon_global(i,j) = (lon_vert_glo(i,j,1) + lon_vert_glo(i+1,j,1) + &
                 lon_vert_glo(i+1,j+1,1) + lon_vert_glo(i,j+1,1))*0.25
            lat_global(i,j) = (lat_vert_glo(i,j,1) + lat_vert_glo(i+1,j,1) + &
                 lat_vert_glo(i+1,j+1,1) + lat_vert_glo(i,j+1,1))*0.25
         enddo
      enddo
   else if( field_exist(grid_file, "ocn_mosaic_file") ) then ! reading from mosaic file
      call field_size(tile_file, 'area', siz)
      if(siz(1) /= nlon*2 .or. siz(2) /= nlat*2) then
         write(message,'(a,2i4)') 'area is wrong shape. shape(area)=',siz(1:2)
         call error_mesg('test_data_override', trim(message), FATAL)
      endif
      allocate(lon_vert_glo(siz(1)+1,siz(2)+1,1), lat_vert_glo(siz(1)+1,siz(2)+1,1))
      allocate(lon_global  (nlon,  nlat    ), lat_global  (nlon,  nlat    ))
      call read_data( tile_file, 'x', lon_vert_glo, no_domain=.true.)
      call read_data( tile_file, 'y', lat_vert_glo, no_domain=.true.)  
      do j = 1, nlat
         do i = 1, nlon
            lon_global(i,j) = lon_vert_glo(i*2,j*2,1)
            lat_global(i,j) = lat_vert_glo(i*2,j*2,1)
         end do
      end do
   end if

  allocate(lon(is:ie,js:je), lat(is:ie,js:je))
  lon = lon_global(is:ie,js:je)
  lat = lat_global(is:ie,js:je)

  deallocate(lon_vert_glo)
  deallocate(lat_vert_glo)
  deallocate(lon_global)
  deallocate(lat_global)

 end subroutine get_grid

  subroutine test_unstruct_grid( type, Time )

   character(len=*), intent(in) :: type
  type(time_type),              intent(in) :: Time !(target) model time

  integer :: pe, npes
  integer :: nx=128, ny=128, nz=40, stackmax=4000000
  integer :: unit=7
  integer :: stdunit = 6
  logical :: debug=.FALSE., opened
 
  integer :: mpes = 0
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
  integer :: x_cyclic_offset = 3   ! to be used in test_cyclic_offset
  integer :: y_cyclic_offset = -4  ! to be used in test_cyclic_offset
  character(len=32) :: warn_level = "fatal"
  integer :: wide_halo_x = 0, wide_halo_y = 0
  integer :: nx_cubic = 96, ny_cubic = 96
  logical :: test_performance = .false.
  logical :: test_interface = .true.
  logical :: test_nest_domain = .false.
  logical :: test_edge_update = .false.
  logical :: test_group = .false.
  logical :: test_cubic_grid_redistribute = .false.
  logical :: check_parallel = .FALSE.  ! when check_parallel set to false,
  logical :: test_get_nbr = .FALSE.
  logical :: test_boundary = .false.
  logical :: test_global_sum = .false.
  integer :: ensemble_size
  integer :: layout_cubic(2) = (/0,0/)
  integer :: layout_tripolar(2) = (/0,0/)
  integer :: layout_ensemble(2) = (/0,0/)
  logical :: do_sleep = .false.
  integer :: num_iter = 1
  integer :: num_fields = 4

  !--- namelist variable for nest domain
  integer :: tile_fine   = 1
  integer :: tile_coarse = 1
  integer :: istart_fine = 0, iend_fine = -1, jstart_fine = 0, jend_fine = -1
  integer :: istart_coarse = 0, iend_coarse = -1, jstart_coarse = 0, jend_coarse = -1
  integer :: npes_coarse = 0
  integer :: npes_fine   = 0
  integer :: extra_halo = 0
  logical :: mix_2D_3D = .false.
  logical :: test_subset = .false.
  logical :: test_unstruct = .false.
  integer :: nthreads = 1
  integer :: i, j, k, l, shift
  integer :: layout(2)
  integer :: id
  integer :: outunit, errunit, io_status
  integer :: get_cpu_affinity, base_cpu, omp_get_num_threads, omp_get_thread_num

    type(domain2D) :: SG_domain
    type(domainUG) :: UG_domain
    integer        :: num_contact, ntiles, npes_per_tile
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer        :: ism, iem, jsm, jem, lsg, leg

    integer, allocatable, dimension(:)       :: pe_start, pe_end, npts_tile, grid_index, ntiles_grid
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:)     :: x1, x2, g1, g2
    real,    allocatable, dimension(:,:,:)   :: a1, a2, gdata
    real,    allocatable, dimension(:,:)     :: rmask
    real,    allocatable, dimension(:)       :: frac_crit
    logical, allocatable, dimension(:,:,:)   :: lmask,msk
    integer, allocatable, dimension(:)       :: isl, iel, jsl, jel
    logical            :: cubic_grid
    character(len=3)   :: text
    integer            :: nx_save, ny_save, tile
    integer            :: ntotal_land, istart, iend, pos
    
  call mpp_memuse_begin()
  call mpp_init()
  npes = mpp_npes()
 
  outunit = stdout()
  errunit = stderr()
    cubic_grid         = .false.

    nx_save = nx
    ny_save = ny
    !--- check the type
    select case(type)
    case ( 'Cubic-Grid' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_unstruct_update: for Cubic_grid mosaic, nx_cubic is zero, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_unstruct_update: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       nx = nx_cubic
       ny = ny_cubic
       ntiles = 6
       num_contact = 12
       cubic_grid = .true.
       if( mod(npes, ntiles) == 0 ) then
          npes_per_tile = npes/ntiles
          write(outunit,*)'NOTE from test_unstruct_update ==> For Mosaic "', trim(type), &
               '", each tile will be distributed over ', npes_per_tile, ' processors.'
       else
          call mpp_error(NOTE,'test_unstruct_update: npes should be multiple of ntiles No test is done for '//trim(type))
          return
       endif
       if(layout_cubic(1)*layout_cubic(2) == npes_per_tile) then
          layout = layout_cubic
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
       allocate(frac_crit(ntiles))
       frac_crit(1) = 0.3; frac_crit(2) = 0.1; frac_crit(3) = 0.6
       frac_crit(4) = 0.2; frac_crit(5) = 0.4; frac_crit(6) = 0.5

    case default
       call mpp_error(FATAL, 'test_group_update: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    do n = 1, ntiles
       pe_start(n) = (n-1)*npes_per_tile
       pe_end(n)   = n*npes_per_tile-1
    end do

    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

    !--- define domain
    if( cubic_grid ) then
       call define_cubic_mosaic(type, SG_domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
            global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( SG_domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( SG_domain, isd, ied, jsd, jed )

    allocate(lmask(nx,ny,ntiles))
    allocate(npts_tile(ntiles))
    lmask = .false.
    if(mpp_pe() == mpp_root_pe() ) then
       allocate(rmask(nx,ny))
       !--- construct gmask.
       do n = 1, ntiles
          call random_number(rmask)
          do j = 1, ny
             do i = 1, nx
                if(rmask(i,j) > frac_crit(n)) then
                   lmask(i,j,n) = .true.
                endif
             enddo
          enddo
          npts_tile(n) = count(lmask(:,:,n))
       enddo
       ntotal_land = sum(npts_tile)
       allocate(grid_index(ntotal_land))
       l = 0
       allocate(isl(0:mpp_npes()-1), iel(0:mpp_npes()-1))
       allocate(jsl(0:mpp_npes()-1), jel(0:mpp_npes()-1))
       call mpp_get_compute_domains(SG_domain,xbegin=isl,xend=iel,ybegin=jsl,yend=jel)

       do n = 1, ntiles
          do j = 1, ny
             do i = 1, nx
                if(lmask(i,j,n)) then
                   l = l + 1
                   grid_index(l) = (j-1)*nx+i
                endif
             enddo
          enddo
       enddo
       deallocate(rmask, isl, iel, jsl, jel)
    endif
    call mpp_broadcast(npts_tile, ntiles, mpp_root_pe())
    if(mpp_pe() .NE. mpp_root_pe()) then
       ntotal_land = sum(npts_tile)
       allocate(grid_index(ntotal_land))
    endif
    call mpp_broadcast(grid_index, ntotal_land, mpp_root_pe())
    
    allocate(ntiles_grid(ntotal_land))
    ntiles_grid = 1
   !--- define the unstructured grid domain
    call mpp_define_unstruct_domain(UG_domain, SG_domain, npts_tile, ntiles_grid, mpp_npes(), 1, grid_index, name="LAND unstruct")
    call mpp_get_UG_compute_domain(UG_domain, istart, iend)

    !--- figure out lmask according to grid_index
    pos = 0
    do n = 1, ntiles
       do l = 1, npts_tile(n)
          pos = pos + 1
          j = (grid_index(pos)-1)/nx + 1
          i = mod((grid_index(pos)-1),nx) + 1
          lmask(i,j,n) = .true.
       enddo
    enddo

    !--- set up data
!    allocate(gdata(nx,ny,ntiles))
!    gdata = -999
!    do n = 1, ntiles
!       do j = 1, ny
!          do i = 1, nx
!             if(lmask(i,j,n)) then
!                gdata(i,j,n) = n*1.e+3 + i + j*1.e-3
!             endif
!          end do
!       end do
!    end do

    !--- test the 2-D data is on computing domain
    allocate( a1(isc:iec, jsc:jec,1), a2(isc:iec,jsc:jec,1 ) )
    allocate(msk(isc:iec, jsc:jec,1)); msk = .false.

    tile = mpp_pe()/npes_per_tile + 1
    do j = jsc, jec
       do i = isc, iec
!          a1(i,j,1) = gdata(i,j,tile)
          msk(i,j,1) = lmask(i,j,tile)
       enddo
    enddo
    !First override the test SG data from file/field
    call data_override_init(Land_domain_in=SG_domain)
    call data_override('LND','sst_obs',a1(:,:,1),Time)

    !Create the test UG data
    a2 = -9999
    !For this test on non-Land points a2 must match a1
    do j = jsc, jec
       do i = isc, iec
          if(.NOT. msk(i,j,1)) a2(i,j,1)=a1(i,j,1)
       enddo
    enddo

    allocate(x1(istart:iend,1), x2(istart:iend,1))
    x1 = -99999
    x2 = -999999
    !--- fill the value of x2

    !Now override the test UG data from the same file/field
    call data_override_init(Land_domainUG_in=UG_domain)
    call data_override_UG('LND','sst_obs',x2(:,1),Time)

    !Ensure you get the same UG data from the SG data
    call mpp_pass_SG_to_UG(UG_domain, a1(:,:,1), x1(:,1))
    call compare_checksums_2D(x1, x2, type//' SG2UG 2-D compute domain')

    !Ensure you get the same SG data from the UG data if you transform back
    call mpp_pass_UG_to_SG(UG_domain, x1(:,1), a2(:,:,1))
    call compare_checksums(a1(:,:,1:1),a2(:,:,1:1),type//' UG2SG 2-D compute domain')

    deallocate(a1,a2,x1,x2)
   
    !--- test the 3-D data is on computing domain
    allocate( a1(isc:iec, jsc:jec,nz), a2(isc:iec,jsc:jec,nz ) )
    
!    tile = mpp_pe()/npes_per_tile + 1
!    do k = 1, nz
!       do j = jsc, jec
!          do i = isc, iec
!             a1(i,j,k) = gdata(i,j,tile) 
!             if(a1(i,j,k) .NE. -999) a1(i,j,k) = a1(i,j,k) + k*1.e-6
!          enddo
!       enddo
!    enddo

    !First override the test SG data from file/field
    call data_override_init(Land_domain_in=SG_domain)
    call data_override('LND','sst_obs',a1,Time)

    a2 = -999
    !For this test on non-Land points a2 must match a1
    do k = 1, nz
    do j = jsc, jec
       do i = isc, iec
          if(.NOT. msk(i,j,1)) a2(i,j,k)=a1(i,j,k)
       enddo
    enddo
    enddo

    allocate(x1(istart:iend,nz), x2(istart:iend,nz))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
 !   tile = mpp_get_UG_domain_tile_id(UG_domain)
 !   pos = 0
 !   do n = 1, tile-1
 !      pos = pos + npts_tile(n)
 !   enddo
 !   do l = istart, iend
 !      i = mod((grid_index(pos+l)-1), nx) + 1
 !      j = (grid_index(pos+l)-1)/nx + 1
 !      do k = 1, nz
 !         x2(l,k) = gdata(i,j,tile) + k*1.e-6
 !      enddo     
 !   enddo

    !Now override the test UG data from the same file/field
    call data_override_init(Land_domainUG_in=UG_domain)
    call data_override_UG('LND','sst_obs',x2,Time)

    !Ensure you get the same UG data from the SG data
    call mpp_pass_SG_to_UG(UG_domain, a1, x1)
    call compare_checksums_2D(x1, x2, type//' SG2UG 3-D compute domain')
    !Ensure you get the same SG data from the UG data if you transform back
    call mpp_pass_UG_to_SG(UG_domain, x1, a2)
    call compare_checksums(a1,a2,type//' UG2SG 3-D compute domain')
    deallocate(a1,a2,x1,x2)


  end subroutine test_unstruct_grid

  subroutine compare_checksums( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(LONG_KIND) :: sum1, sum2
    integer :: i, j, k,pe

    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
    ! mpp_sync()
    call mpp_sync_self()
  pe = mpp_pe()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
         call mpp_error(FATAL,'compare_chksum: size of a and b does not match')

!print*, "pe,a(1,1,1),b(1,1,1)=", pe,a(1,1,1),b(1,1,1)

    do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,k) .ne. b(i,j,k)) then
            print*, "pe,i,j,k", pe,i,j,k
            print*, "a =", a(i,j,k)
            print*, "b =", b(i,j,k)
!                write(stdunit,'(a,i3,a,i3,a,i3,a,i3,a,f20.9,a,f20.9)')" at pe ", mpp_pe(), &
!                     ", at point (",i,", ", j, ", ", k, "), a = ", a(i,j,k), ", b = ", b(i,j,k)
                call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
             endif
          enddo
       enddo
    enddo

    sum1 = mpp_chksum( a, (/pe/) )
    sum2 = mpp_chksum( b, (/pe/) )

    if( sum1.EQ.sum2 )then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
        !--- in some case, even though checksum agree, the two arrays 
        !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
        !--- hence we need to check the value point by point.
    else
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums

  !###########################################################################
  subroutine compare_checksums_2D( a, b, string )
    real, intent(in), dimension(:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(LONG_KIND) :: sum1, sum2
    integer :: i, j,pe

    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
    ! mpp_sync()
    call mpp_sync_self()
  pe = mpp_pe()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) ) &
         call mpp_error(FATAL,'compare_chksum_2D: size of a and b does not match')

!print*, "a(1,1),b(1,1)=", a(1,1),b(1,1)

    do j = 1, size(a,2)
       do i = 1, size(a,1)
          if(a(i,j) .ne. b(i,j)) then
            print*, "i,j= ", i,j
            print*, "a =", a(i,j)
            print*, "b =", b(i,j)
!             write(stdunit,'(a,i3,a,i3,a,i3,a,f20.9,a,f20.9)')"at pe ", mpp_pe(), &
!                  ", at point (",i,", ", j, "),a=", a(i,j), ",b=", b(i,j)
             call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
          endif
       enddo
    enddo

    sum1 = mpp_chksum( a, (/pe/) )
    sum2 = mpp_chksum( b, (/pe/) )

    if( sum1.EQ.sum2 )then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
        !--- in some case, even though checksum agree, the two arrays 
        !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
        !--- hence we need to check the value point by point.
    else
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums_2D
  subroutine define_cubic_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end)
    character(len=*), intent(in)  :: type
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: global_indices(:,:), layout(:,:)
    integer,        intent(in)    :: ni(:), nj(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact, msize(2)
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2


    ntiles = 6
    num_contact = 12
    if(size(pe_start(:)) .NE. 6 .OR. size(pe_end(:)) .NE. 6 ) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of pe_start and pe_end should be 6")
    if(size(global_indices,1) .NE. 4) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of first dimension of global_indices should be 4")
    if(size(global_indices,2) .NE. 6) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of second dimension of global_indices should be 6")
    if(size(layout,1) .NE. 2) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of first dimension of layout should be 2")
    if(size(layout,2) .NE. 6) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of second dimension of layout should be 6")
    if(size(ni(:)) .NE. 6 .OR. size(nj(:)) .NE. 6) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of ni and nj should be 6")

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1; tile2(1) = 2
    istart1(1) = ni(1);  iend1(1) = ni(1);  jstart1(1) = 1;      jend1(1) = nj(1)
    istart2(1) = 1;      iend2(1) = 1;      jstart2(1) = 1;      jend2(1) = nj(2)
    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni(1);  jstart1(2) = nj(1);  jend1(2) = nj(1)
    istart2(2) = 1;      iend2(2) = 1;      jstart2(2) = nj(3);  jend2(2) = 1
    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1; tile2(3) = 5
    istart1(3) = 1;      iend1(3) = 1;      jstart1(3) = 1;      jend1(3) = nj(1)
    istart2(3) = ni(5);  iend2(3) = 1;      jstart2(3) = nj(5);  jend2(3) = nj(5)
    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni(1);  jstart1(4) = 1;      jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni(6);  jstart2(4) = nj(6);  jend2(4) = nj(6)       
    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2; tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni(2);  jstart1(5) = nj(2);  jend1(5) = nj(2)
    istart2(5) = 1;      iend2(5) = ni(3);  jstart2(5) = 1;      jend2(5) = 1
    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni(2);  iend1(6) = ni(2);  jstart1(6) = 1;      jend1(6) = nj(2)
    istart2(6) = ni(4);  iend2(6) = 1;      jstart2(6) = 1;      jend2(6) = 1
    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;      iend1(7) = ni(2);  jstart1(7) = 1;      jend1(7) = 1
    istart2(7) = ni(6);  iend2(7) = ni(6);  jstart2(7) = nj(6);  jend2(7) = 1
    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni(3);  iend1(8) = ni(3);  jstart1(8) = 1;      jend1(8) = nj(3)
    istart2(8) = 1;      iend2(8) = 1;      jstart2(8) = 1;      jend2(8) = nj(4)
    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni(3);  jstart1(9) = nj(3);  jend1(9) = nj(3)
    istart2(9) = 1;      iend2(9) = 1;      jstart2(9) = nj(5);  jend2(9) = 1
    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni(4); jstart1(10) = nj(4); jend1(10) = nj(4)
    istart2(10) = 1;     iend2(10) = ni(5); jstart2(10) = 1;     jend2(10) = 1
    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni(4); iend1(11) = ni(4); jstart1(11) = 1;     jend1(11) = nj(4)
    istart2(11) = ni(6); iend2(11) = 1;     jstart2(11) = 1;     jend2(11) = 1
    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni(5); iend1(12) = ni(5); jstart1(12) = 1;     jend1(12) = nj(5)
    istart2(12) = 1;     iend2(12) = 1;     jstart2(12) = 1;     jend2(12) = nj(6)
    msize(1) = maxval(ni(:)/layout(1,:)) + whalo + ehalo + 1 ! make sure memory domain size is no smaller than
    msize(2) = maxval(nj(:)/layout(2,:)) + shalo + nhalo + 1 ! data domain size       
    call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, tile2, &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
         pe_start, pe_end, symmetry = .true., whalo=whalo, ehalo=ehalo,   &
         shalo=shalo, nhalo=nhalo, name = trim(type), memory_size = msize  )  

    return 

  end subroutine define_cubic_mosaic

!=================================================================================================================================
 end program test
#endif
