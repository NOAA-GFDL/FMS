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
!> @defgroup interpolator_mod interpolator_mod
!> @ingroup interpolator
!> @brief A module to interpolate climatology data to model the grid.
!> @author William Cooke <William.Cooke@noaa.gov>

module interpolator_mod

use mpp_mod,           only : mpp_error, &
                              FATAL,     &
                              mpp_pe,    &
                              mpp_init,  &
                              mpp_exit,  &
                              mpp_npes,  &
                              WARNING,   &
                              NOTE,      &
                              input_nml_file
use mpp_domains_mod,   only : mpp_domains_init,      &
                              mpp_update_domains,    &
                              mpp_define_domains,    &
                              mpp_global_field,      &
                              domain2d,              &
                              mpp_define_layout,     &
                              mpp_get_compute_domain
use diag_manager_mod,  only : diag_manager_init, get_base_time, &
                              register_diag_field, send_data, &
                              diag_axis_init
use fms_mod,           only : lowercase, write_version_number, &
                              fms_init, &
                              mpp_root_pe, stdlog, &
                              check_nml_error
use fms2_io_mod,       only : FmsNetcdfFile_t, fms2_io_file_exist => file_exists, dimension_exists, &
                              open_file, fms2_io_read_data=>read_data,    &
                              variable_exists, get_variable_num_dimensions, &
                              get_num_variables, get_dimension_size,   &
                              get_variable_units, get_variable_names,  &
                              get_time_calendar, close_file,           &
                              get_variable_dimension_names, get_variable_sense

use horiz_interp_mod,  only : horiz_interp_type, &
                              horiz_interp_new,  &
                              horiz_interp_init, &
                              assignment(=), &
                              horiz_interp,      &
                              horiz_interp_del
use time_manager_mod,  only : time_type,   &
                              set_time,    &
                              set_date,    &
                              time_type_to_real, &
                              days_in_year, &
                              get_calendar_type, &
                              leap_year, &
                              JULIAN, NOLEAP, &
                              get_date, &
                              get_date_julian, set_date_no_leap, &
                              set_date_julian, get_date_no_leap, &
                              print_date, &
                              operator(+), &
                              operator(-), &
                              operator(*), &
                              operator(>), &
                              operator(<), &
                              assignment(=), &
                              decrement_time
use time_interp_mod,   only : time_interp, YEAR
use constants_mod,     only : grav, PI, SECONDS_PER_DAY
use platform_mod,      only : r4_kind, r8_kind, r16_kind

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!-------  interfaces --------

public interpolator_init, &
       interpolator,      &
       interpolate_type_eq, &
       obtain_interpolator_time_slices, &
       unset_interpolator_time_flag, &
       interpolator_end,  &
       init_clim_diag,    &
       query_interpolator,&
       read_data

!> Interpolates a field to a model grid
!!
!> Example usage:
!! ~~~~~~~~~~{.f90}
!! call interpolator (sulfate, model_time, p_half, model_data, name, is, js, clim_units)
!! call interpolator (o3, model_time, p_half, model_data, "ozone", is, js)
!! ~~~~~~~~~~
!!
!! The first option is used to generate sulfate models.
!!
!! The sulfate data is set by
!! ~~~~~~~~~~{.f90}
!! type(interpolate_type), intent(inout) :: sulfate
!! ~~~~~~~~~~
!! The name of the model is set by
!! ~~~~~~~~~~{.f90}
!! character(len=*), intent(in) :: name
!! ~~~~~~~~~~
!! The units used in this model are outputted to
!! ~~~~~~~~~~{.f90}
!! character(len=*), intent(out), optional :: clim_units
!! ~~~~~~~~~~
!!
!! The second option is generate ozone models.
!!
!! The ozone data is set by
!! ~~~~~~~~~~{.f90}
!! type(interpolate_type), intent(inout) :: o3
!! ~~~~~~~~~~
!!
!! Both of these options use the following variables in the model.
!!
!! The time used in the model is set by
!!
!! ~~~~~~~~~~{.f90}
!! type(time_type), intent(in) :: model_time
!! ~~~~~~~~~~
!! The model pressure field is set by
!! ~~~~~~~~~~{.f90}
!! real, intent(in), dimension(:,:,:) :: p_half
!! ~~~~~~~~~~
!!
!! @param [inout] <clim_type> The interpolate type previously defined by a call to interpolator_init
!! @param [in] <field_name> The name of a field that you wish to interpolate
!! @param [in] <Time> The model time that you wish to interpolate to
!! @param [in] <phalf> The half level model pressure field
!! @param [in] <is> Index for the physics window
!! @param [in] <js> Index for the physics window
!! @param [out] <interp_data> The model fields with the interpolated climatology data
!! @param [out] <clim_units> The units of field_name
!> @ingroup interpolator_mod
interface interpolator
   module procedure interpolator_4D_r4, interpolator_4D_r8
   module procedure interpolator_3D_r4, interpolator_3D_r8
   module procedure interpolator_2D_r4, interpolator_2D_r8
   module procedure interpolator_4D_no_time_axis_r4, interpolator_4D_no_time_axis_r8
   module procedure interpolator_3D_no_time_axis_r4, interpolator_3D_no_time_axis_r8
   module procedure interpolator_2D_no_time_axis_r4, interpolator_2D_no_time_axis_r8
end interface interpolator

!> Private assignment override interface for interpolate type
!> @ingroup interpolator_mod
interface assignment(=)
   module procedure interpolate_type_eq
end interface

interface interpolator_init
   module procedure interpolator_init_r4
   module procedure interpolator_init_r8
end interface interpolator_init

interface fms2io_interpolator_init
   module procedure fms2io_interpolator_init_r4
   module procedure fms2io_interpolator_init_r8
end interface fms2io_interpolator_init

interface get_axis_latlon_data
   module procedure get_axis_latlon_data_r4
   module procedure get_axis_latlon_data_r8
end interface get_axis_latlon_data

interface get_axis_level_data
   module procedure get_axis_level_data_r4
   module procedure get_axis_level_data_r8
end interface get_axis_level_data

interface cell_center2
   module procedure cell_center2_r4
   module procedure cell_center2_r8
end interface cell_center2

interface cart_to_latlon
   module procedure cart_to_latlon_r4
   module procedure cart_to_latlon_r8
end interface cart_to_latlon

interface latlon2xyz
   module procedure latlon2xyz_r4
   module procedure latlon2xyz_r8
end interface latlon2xyz

interface diag_read_data
  module procedure diag_read_data_r4
  module procedure diag_read_data_r8
end interface diag_read_data

interface read_data
   module procedure read_data_r4
   module procedure read_data_r8
end interface read_data

interface read_data_no_time_axis
   module procedure read_data_no_time_axis_r4
   module procedure read_data_no_time_axis_r8
end interface read_data_no_time_axis

interface interp_linear
   module procedure interp_linear_r4
   module procedure interp_linear_r8
end interface interp_linear

!> Private interface for weighted scalar interpolation
!!
!> Example usage:
!! ~~~~~~~~~~{.f90}
!! call interp_weighted_scalar (pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:,:),interp_data(ilon,j,:,:))
!! call interp_weighted_scalar (pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
!! ~~~~~~~~~~
!!
!! @param [in] <grdin> Input grid
!! @param [in] <grdout> Output grid
!! @param [in] <datin> Input data
!! @param [out] <datout> Output data
!> @ingroup interpolator_mod
interface interp_weighted_scalar
   module procedure interp_weighted_scalar_1D_r4, interp_weighted_scalar_1D_r8
   module procedure interp_weighted_scalar_2D_r4, interp_weighted_scalar_2d_r8
end interface interp_weighted_scalar

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

! Include variable "version" to be written to log file.
#include<file_version.h>

!> Redundant climatology data between fields
!> @ingroup interpolate_type

type, private :: interpolate_r4_type
logical :: is_allocated = .false.
real(r4_kind), allocatable :: lat(:)               !< No description
real(r4_kind), allocatable :: lon(:)               !< No description
real(r4_kind), allocatable :: latb(:)              !< No description
real(r4_kind), allocatable :: lonb(:)              !< No description
real(r4_kind), allocatable :: levs(:)              !< No description
real(r4_kind), allocatable :: halflevs(:)          !< No description
real(r4_kind), allocatable :: data5d(:,:,:,:,:)  !< (nlatmod,nlonmod,nlevclim,size(time_init,2),nfields)
real(r4_kind), allocatable :: pmon_pyear(:,:,:,:)           !< No description
real(r4_kind), allocatable :: pmon_nyear(:,:,:,:)           !< No description
real(r4_kind), allocatable :: nmon_nyear(:,:,:,:)           !< No description
real(r4_kind), allocatable :: nmon_pyear(:,:,:,:)           !< No description
real(r4_kind) :: tweight      !< No description
real(r4_kind) :: tweight1     !< The time weight between the climatology years
real(r4_kind) :: tweight2     !< No description
real(r4_kind) :: tweight3     !< The time weight between the month
end type interpolate_r4_type

type, private :: interpolate_r8_type
logical :: is_allocated = .false.
real(r8_kind), allocatable :: lat(:)               !< No description
real(r8_kind), allocatable :: lon(:)               !< No description
real(r8_kind), allocatable :: latb(:)              !< No description
real(r8_kind), allocatable :: lonb(:)              !< No description
real(r8_kind), allocatable :: levs(:)              !< No description
real(r8_kind), allocatable :: halflevs(:)          !< No description
real(r8_kind), allocatable :: data5d(:,:,:,:,:)  !< (nlatmod,nlonmod,nlevclim,size(time_init,2),nfields)
real(r8_kind), allocatable :: pmon_pyear(:,:,:,:)           !< No description
real(r8_kind), allocatable :: pmon_nyear(:,:,:,:)           !< No description
real(r8_kind), allocatable :: nmon_nyear(:,:,:,:)           !< No description
real(r8_kind), allocatable :: nmon_pyear(:,:,:,:)           !< No description
real(r8_kind) :: tweight          !< No description
real(r8_kind) :: tweight1     !< The time weight between the climatology years
real(r8_kind) :: tweight2     !< No description
real(r8_kind) :: tweight3     !< The time weight between the month
end type interpolate_r8_type

type, public  :: interpolate_type
private
!Redundant data between fields
!All climatology data
type(interpolate_r4_type) :: r4_type
type(interpolate_r8_type) :: r8_type
type(horiz_interp_type)  :: interph                         !< No description
type(time_type), allocatable :: time_slice(:) !< An array of the times within the climatology.
type(FmsNetcdfFile_t)    :: fileobj       ! object that stores opened file information
character(len=64)        :: file_name     !< Climatology filename
integer                  :: TIME_FLAG     !< Linear or seaonal interpolation?
integer                  :: level_type    !< Pressure or Sigma level
integer                  :: is,ie,js,je       !< No description
integer                  :: vertical_indices !< direction of vertical
                                              !! data axis
logical                  :: climatological_year !< Is data for year = 0000?

!Field specific data  for nfields
character(len=64), allocatable :: field_name(:)    !< name of this field
logical,           allocatable :: has_level(:)     !< indicate if the variable has level dimension
integer,           allocatable :: time_init(:,:)   !< second index is the number of time_slices being
                                                       !! kept. 2 or ntime.
integer,           allocatable :: mr(:)            !< Flag for conversion of climatology to mixing ratio.
integer,           allocatable :: out_of_bounds(:) !< Flag for when surface pressure is out of bounds.
!++lwh
integer,           allocatable :: vert_interp(:)   !< Flag for type of vertical interpolation.
!--lwh
!integer                    :: indexm, indexp, climatology
integer,dimension(:),  allocatable :: indexm        !< No description
integer,dimension(:),  allocatable :: indexp        !< No description
integer,dimension(:),  allocatable :: climatology   !< No description

type(time_type), allocatable :: clim_times(:,:) !< No description
logical :: separate_time_vary_calc              !< No description
integer :: itaum     !< No description
integer :: itaup     !< No description

end type interpolate_type

!> @addtogroup interpolator_mod
!> @{

logical            :: module_is_initialized = .false.
logical            :: clim_diag_initialized = .false.

integer :: ndim      !< No description
integer :: nvar      !< No description
integer :: ntime     !< No description
integer :: nlat      !< No description
integer :: nlatb     !< No description
integer :: nlon      !< No description
integer :: nlonb     !< No description
integer :: nlev      !< No description
integer :: nlevh     !< No description
integer ::          len, ntime_in, num_fields               !< No description

! pletzer real, allocatable :: time_in(:)
! sjs real, allocatable :: climdata(:,:,:), climdata2(:,:,:)

character(len=64) :: units                              !< No description
integer           :: sense                                        !< No description

integer, parameter :: max_diag_fields = 30                    !< No description

! flags to indicate direction of vertical axis in  data file
integer, parameter :: INCREASING_DOWNWARD = 1, INCREASING_UPWARD = -1          !< Flags to indicate direction
                                                                               !! of vertical axis in  data file
!++lwh
! Flags to indicate whether the time interpolation should be linear or some other scheme for seasonal data.
! NOTIME indicates that data file has no time axis.
integer, parameter :: LINEAR = 1, SEASONAL = 2, BILINEAR = 3, NOTIME = 4     !< Flags to indicate whether the time
                                                                             !! interpolation should be linear or some
                                                                             !! other scheme for seasonal data.
                                                                             !! NOTIME indicates
                                                                             !! that data file has no time axis.

! Flags to indicate where climatology pressure levels are pressure or sigma levels
integer, parameter :: PRESSURE = 1, SIGMA = 2          !< Flags to indicate where climatology pressure
                                                       !! levels are pressure or sigma levels

! Flags to indicate whether the climatology units are mixing ratio (kg/kg) or column integral (kg/m2).
! Vertical interpolation scheme requires mixing ratio at this time.
integer, parameter :: NO_CONV = 1, KG_M2 = 2          !< Flags to indicate whether the climatology units
                                                      !! are mixing ratio (kg/kg) or column integral (kg/m2).
                                                      !! Vertical interpolation scheme requires mixing ratio at
                                                      !! this time.

! Flags to indicate what to do when the model surface pressure exceeds the  climatology surface pressure level.
integer, parameter, public :: CONSTANT = 1, ZERO = 2  !< Flags to indicate what to do when the model surface
                                                      !! pressure exceeds the climatology surface pressure level.

! Flags to indicate the type of vertical interpolation
integer, parameter, public :: INTERP_WEIGHTED_P = 10, INTERP_LINEAR_P = 20, INTERP_LOG_P = 30 !< Flags to indicate
                                                                             !! the type of vertical interpolation
!--lwh

real(r8_kind), parameter :: TPI = (2.0_r8_kind*PI) ! 4.*acos(0.)
real(r8_kind), parameter :: DTR = TPI/360._r8_kind



integer :: num_clim_diag = 0                                            !< No description
character(len=64) :: climo_diag_name(max_diag_fields)                   !< No description
integer :: climo_diag_id(max_diag_fields), hinterp_id(max_diag_fields)  !< No description
real(r8_kind) ::  missing_value = -1.e10_r8_kind                        !< No description
! sjs integer :: itaum, itaup

#ifdef ENABLE_QUAD_PRECISION
! Higher precision (kind=16) for grid geometrical factors:
 integer, parameter:: f_p = r16_kind    !< Higher precision (kind=16) for grid geometrical factors
#else
! 64-bit precision (kind=8)
 integer, parameter:: f_p = r8_kind     !< 64-bit precision (kind=8)
#endif

logical :: read_all_on_init = .false.          !< No description
integer :: verbose = 0                         !< No description
logical :: conservative_interp = .true.        !< No description
logical :: retain_cm3_bug = .false.            !< No description
logical :: use_mpp_io = .false. !< Set to true to use mpp_io, otherwise fms2io is used

namelist /interpolator_nml/    &
                             read_all_on_init, verbose, conservative_interp, retain_cm3_bug, use_mpp_io

contains

!> @brief Assignment overload routine for interpolate_type, to be used
!! through the assignment interface
subroutine interpolate_type_eq (Out, In)

type(interpolate_type), intent(in) :: In
type(interpolate_type), intent(inout) :: Out

     if(In%r4_type%is_allocated) then
        if (allocated(In%r4_type%lat))      Out%r4_type%lat      =  In%r4_type%lat
        if (allocated(In%r4_type%lon))      Out%r4_type%lon      =  In%r4_type%lon
        if (allocated(In%r4_type%latb))     Out%r4_type%latb     =  In%r4_type%latb
        if (allocated(In%r4_type%lonb))     Out%r4_type%lonb     =  In%r4_type%lonb
        if (allocated(In%r4_type%levs))     Out%r4_type%levs     =  In%r4_type%levs
        if (allocated(In%r4_type%halflevs)) Out%r4_type%halflevs =  In%r4_type%halflevs
     else if(In%r8_type%is_allocated) then
        if (allocated(In%r8_type%lat))      Out%r8_type%lat      =  In%r8_type%lat
        if (allocated(In%r8_type%lon))      Out%r8_type%lon      =  In%r8_type%lon
        if (allocated(In%r8_type%latb))     Out%r8_type%latb     =  In%r8_type%latb
        if (allocated(In%r8_type%lonb))     Out%r8_type%lonb     =  In%r8_type%lonb
        if (allocated(In%r8_type%levs))     Out%r8_type%levs     =  In%r8_type%levs
        if (allocated(In%r8_type%halflevs)) Out%r8_type%halflevs =  In%r8_type%halflevs
     end if

     Out%interph = In%interph
     if (allocated(In%time_slice)) Out%time_slice =  In%time_slice
     Out%file_name  = In%file_name
     Out%time_flag  = In%time_flag
     Out%level_type = In%level_type
     Out%is = In%is
     Out%ie = In%ie
     Out%js = In%js
     Out%je = In%je
     Out%vertical_indices = In%vertical_indices
     Out%climatological_year = In%climatological_year
     if (allocated(In%has_level    )) Out%has_level     =  In%has_level
     if (allocated(In%field_name   )) Out%field_name    =  In%field_name
     if (allocated(In%time_init    )) Out%time_init     =  In%time_init
     if (allocated(In%mr           )) Out%mr            =  In%mr
     if (allocated(In%out_of_bounds)) Out%out_of_bounds =  In%out_of_bounds
     if (allocated(In%vert_interp  )) Out%vert_interp   =  In%vert_interp
     if(In%r4_type%is_allocated) then
        if (allocated(In%r4_type%data5d       )) Out%r4_type%data5d        =  In%r4_type%data5d
        if (allocated(In%r4_type%pmon_pyear   )) Out%r4_type%pmon_pyear    =  In%r4_type%pmon_pyear
        if (allocated(In%r4_type%pmon_nyear   )) Out%r4_type%pmon_nyear    =  In%r4_type%pmon_nyear
        if (allocated(In%r4_type%nmon_nyear   )) Out%r4_type%nmon_nyear    =  In%r4_type%nmon_nyear
        if (allocated(In%r4_type%nmon_pyear   )) Out%r4_type%nmon_pyear    =  In%r4_type%nmon_pyear
     else if(In%r8_type%is_allocated) then
        if (allocated(In%r8_type%data5d       )) Out%r8_type%data5d        =  In%r8_type%data5d
        if (allocated(In%r8_type%pmon_pyear   )) Out%r8_type%pmon_pyear    =  In%r8_type%pmon_pyear
        if (allocated(In%r8_type%pmon_nyear   )) Out%r8_type%pmon_nyear    =  In%r8_type%pmon_nyear
        if (allocated(In%r8_type%nmon_nyear   )) Out%r8_type%nmon_nyear    =  In%r8_type%nmon_nyear
        if (allocated(In%r8_type%nmon_pyear   )) Out%r8_type%nmon_pyear    =  In%r8_type%nmon_pyear
     end if
     if (allocated(In%indexm       )) Out%indexm        =  In%indexm
     if (allocated(In%indexp       )) Out%indexp        =  In%indexp
     if (allocated(In%climatology  )) Out%climatology   =  In%climatology
     if (allocated(In%clim_times   )) Out%clim_times    =  In%clim_times
     Out%separate_time_vary_calc = In%separate_time_vary_calc
     if(In%r4_type%is_allocated) then
        Out%r4_type%tweight  = In%r4_type%tweight
        Out%r4_type%tweight1 = In%r4_type%tweight1
        Out%r4_type%tweight2 = In%r4_type%tweight2
        Out%r4_type%tweight3 = In%r4_type%tweight3
     else if(In%r8_type%is_allocated) then
        Out%r8_type%tweight  = In%r8_type%tweight
        Out%r8_type%tweight1 = In%r8_type%tweight1
        Out%r8_type%tweight2 = In%r8_type%tweight2
        Out%r8_type%tweight3 = In%r8_type%tweight3
     end if
     Out%itaum = In%itaum
     Out%itaup = In%itaup

     Out%r4_type%is_allocated = In%r4_type%is_allocated
     Out%r8_type%is_allocated = Out%r8_type%is_allocated

end subroutine interpolate_type_eq




!#######################################################################
!
!---------------------------------------------------------------------
!> @brief check_climo_units checks the units that the climatology
!!        data is using. This is needed to allow for conversion of
!!        datasets to mixing ratios which is what the vertical
!!        interpolation scheme requires. The default is to assume no
!!        conversion is needed. If the units are those of a column
!!        burden (kg/m2) then conversion to mixing ratio is flagged.
!!
!! @param [in] <units> The units which you will be checking
function check_climo_units(units)
! Function to check the units that the climatology data is using.
! This is needed to allow for conversion of datasets to mixing ratios which is what the
! vertical interpolation scheme requires
! The default is to assume no conversion is needed.
! If the units are those of a column burden (kg/m2) then conversion to mixing ratio is flagged.
!
character(len=*), intent(in) :: units

integer :: check_climo_units

check_climo_units = NO_CONV
select case(chomp(units))
  case('kg/m2', 'kg/m^2', 'kg/m**2', 'kg m^-2', 'kg m**-2')
     check_climo_units = KG_M2
  case('molecules/cm2/s', 'molecule/cm2/s', 'molec/cm2/s')
     check_climo_units = KG_M2
  case('kg/m2/s')
     check_climo_units = KG_M2
end select

end function check_climo_units
!
!#######################################################################
!
!---------------------------------------------------------------------
!> @brief init_clim_diag is a routine to register diagnostic fields
!!        for the climatology file. This routine calculates the domain
!!        decompostion of the climatology fields for later export
!!        through send_data. The ids created here are for column
!!        burdens that will diagnose the vertical interpolation
!!        routine.
!!
!! @param [inout] <clim_type> The interpolate type containing the
!!                      names of the fields in the climatology file
!! @param [in] <mod_axes> The axes of the model
!! @param [in] <init_time> The model initialization time
!!
!! @throw FATAL, "init_clim_diag : You must call interpolator_init before calling init_clim_diag"
!! @throw FATAL, "init_clim_diag : Trying to set up too many diagnostic fields for the climatology data"
subroutine init_clim_diag(clim_type, mod_axes, init_time)
!
! Routine to register diagnostic fields for the climatology file.
! This routine calculates the domain decompostion of the climatology fields
! for later export through send_data.
! The ids created here are for column burdens that will diagnose the vertical interpolation routine.
! climo_diag_id : 'module_name = climo' is intended for use with the model vertical resolution.
! hinterp_id    : 'module_name = 'hinterp' is intended for use with the climatology vertical resolution.

! INTENT INOUT :
!    clim_type : The interpolate type containing the names of the fields in the climatology file.
!
! INTENT IN    :
!   mod_axes   : The axes of the model.
!   init_time  : The model initialization time.
!
type(interpolate_type), intent(inout)  :: clim_type
integer               , intent(in)     :: mod_axes(:)
type(time_type)       , intent(in)     :: init_time

integer :: axes(2),nxd,nyd,ndivs,i
type(domain2d) :: domain
integer :: domain_layout(2), iscomp, iecomp,jscomp,jecomp

if(clim_type%r4_type%is_allocated) call init_clim_diag_r4(clim_type, mod_axes, init_time)
if(clim_type%r8_type%is_allocated) call init_clim_diag_r8(clim_type, mod_axes, init_time)

end subroutine init_clim_diag



!
!---------------------------------------------------------------------
!> @brief obtain_interpolator_time_slices makes sure that the
!!        appropriate time slices are available for interpolation on
!!        this time step.
!!
!! @param [inout] <clim_type> The interpolate type previously defined
!!                      by a call to interpolator_init
!! @param [in] <Time> The model time that you wish to interpolate to
!!
!! @throw FATAL "interpolator_timeslice 1:  file="
!! @throw FATAL "interpolator_timeslice 2:  file="
!! @throw FATAL "interpolator_timeslice 3:  file="
!! @throw FATAL "interpolator_timeslice 4:  file="
!! @throw FATAL "interpolator_timeslice 5:  file="
!! @throw FATAL "interpolator_timeslice : No data from the previous
!!                    climatology time but we have the next time. How did
!!                    this happen?"
subroutine obtain_interpolator_time_slices (clim_type, Time)

!  Makes sure that appropriate time slices are available for interpolation
!  on this time step
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init
!
! INTENT IN
!   Time        : The model time that you wish to interpolate to.

type(interpolate_type), intent(inout)  :: clim_type
type(time_type)       , intent(in)  :: Time

integer :: taum, taup
integer :: modyear, modmonth, modday, modhour, modminute, modsecond
integer :: climyear, climmonth, climday, climhour, climminute, climsecond
integer :: year1, month1, day, hour, minute, second
integer :: climatology, m
type(time_type) :: t_prev, t_next
type(time_type), dimension(2) :: month
integer :: indexm, indexp, yearm, yearp
integer :: i, n
character(len=256) :: err_msg

if(clim_type%r4_type%is_allocated) call obtain_interpolator_time_slices_r4(clim_type, Time)
if(clim_type%r8_type%is_allocated) call obtain_interpolator_time_slices_r8(clim_type, Time)


end subroutine obtain_interpolator_time_slices


!#####################################################################
!
!---------------------------------------------------------------------
!> @brief unset_interpolator_time_flag sets a flag in clim_type to
!!        false.
!!
!! @param [inout] <clim_type> The interpolate type containing the names of the fields in the climatology file
subroutine unset_interpolator_time_flag (clim_type)

type(interpolate_type), intent(inout) :: clim_type


      clim_type%separate_time_vary_calc = .false.


end subroutine unset_interpolator_time_flag


!#####################################################################
!
!---------------------------------------------------------------------
!> @brief interpolator_end receives interpolate data as input
!!        and deallocates its memory.
!!
!! @param [inout] <clim_type> The interpolate type whose components will be deallocated
subroutine interpolator_end(clim_type)
! Subroutine to deallocate the interpolate type clim_type.
!
! INTENT INOUT
!  clim_type : allocate type whose components will be deallocated.
!
type(interpolate_type), intent(inout) :: clim_type
integer :: logunit

logunit=stdlog()
if ( mpp_pe() == mpp_root_pe() ) then
   write (logunit,'(/,(a))') 'Exiting interpolator, have a nice day ...'
end if

if(clim_type%r4_type%is_allocated) then
   if (allocated (clim_type%r4_type%lat     )) deallocate(clim_type%r4_type%lat)
   if (allocated (clim_type%r4_type%lon     )) deallocate(clim_type%r4_type%lon)
   if (allocated (clim_type%r4_type%latb    )) deallocate(clim_type%r4_type%latb)
   if (allocated (clim_type%r4_type%lonb    )) deallocate(clim_type%r4_type%lonb)
   if (allocated (clim_type%r4_type%levs    )) deallocate(clim_type%r4_type%levs)
   if (allocated (clim_type%r4_type%halflevs)) deallocate(clim_type%r4_type%halflevs)
   if (allocated (clim_type%r4_type%data5d  )) deallocate(clim_type%r4_type%data5d)
else if(clim_type%r8_type%is_allocated) then
   if (allocated (clim_type%r8_type%lat     )) deallocate(clim_type%r8_type%lat)
   if (allocated (clim_type%r8_type%lon     )) deallocate(clim_type%r8_type%lon)
   if (allocated (clim_type%r8_type%latb    )) deallocate(clim_type%r8_type%latb)
   if (allocated (clim_type%r8_type%lonb    )) deallocate(clim_type%r8_type%lonb)
   if (allocated (clim_type%r8_type%levs    )) deallocate(clim_type%r8_type%levs)
   if (allocated (clim_type%r8_type%halflevs)) deallocate(clim_type%r8_type%halflevs)
   if (allocated (clim_type%r8_type%data5d))   deallocate(clim_type%r8_type%data5d)
end if

if (allocated (clim_type%time_slice)) deallocate(clim_type%time_slice)
if (allocated (clim_type%field_name)) deallocate(clim_type%field_name)
if (allocated (clim_type%time_init )) deallocate(clim_type%time_init)
if (allocated (clim_type%has_level))  deallocate(clim_type%has_level)
if (allocated (clim_type%mr        )) deallocate(clim_type%mr)
if (allocated (clim_type%out_of_bounds )) deallocate(clim_type%out_of_bounds)
if (allocated (clim_type%vert_interp ))   deallocate(clim_type%vert_interp)
if (allocated(clim_type%indexm)) deallocate(clim_type%indexm)
if (allocated(clim_type%indexp)) deallocate(clim_type%indexp)
if (allocated(clim_type%clim_times))  deallocate(clim_type%clim_times)
if (allocated(clim_type%climatology)) deallocate(clim_type%climatology)

call horiz_interp_del(clim_type%interph)

if(clim_type%r4_type%is_allocated) then
   if (allocated(clim_type%r4_type%pmon_pyear)) then
      deallocate(clim_type%r4_type%pmon_pyear)
      deallocate(clim_type%r4_type%pmon_nyear)
      deallocate(clim_type%r4_type%nmon_nyear)
      deallocate(clim_type%r4_type%nmon_pyear)
   end if
else if(clim_type%r8_type%is_allocated) then
   if (allocated(clim_type%r8_type%pmon_pyear)) then
      deallocate(clim_type%r8_type%pmon_pyear)
      deallocate(clim_type%r8_type%pmon_nyear)
      deallocate(clim_type%r8_type%nmon_nyear)
      deallocate(clim_type%r8_type%nmon_pyear)
   end if
endif

clim_type%r4_type%is_allocated=.false.
clim_type%r8_type%is_allocated=.false.

!! RSH mod
if( .not.(clim_type%TIME_FLAG .eq. LINEAR  .and. read_all_on_init) &
   .and. (clim_type%TIME_FLAG.ne.NOTIME) ) then
!     read_all_on_init)) .or. clim_type%TIME_FLAG .eq. BILINEAR  ) then
     call close_file(clim_type%fileobj)
endif


module_is_initialized = .false.

end subroutine interpolator_end
!
!#######################################################################
!
!++lwh
!
!---------------------------------------------------------------------
!> @brief query_interpolator receives an interpolate type as input
!!        and returns the number of fields and field names.
!!
!! @param [in] <clim_type> The interpolate type which contains the data
!! @param [out] <nfields> OPTIONAL: No description
!! @param [out] <field_names> OPTIONAL: No description
subroutine query_interpolator( clim_type, nfields, field_names )
!
! Query an interpolate_type variable to find the number of fields and field names.
!
type(interpolate_type), intent(in)                    :: clim_type
integer, intent(out), optional                        :: nfields
character(len=*), dimension(:), intent(out), optional :: field_names

if( present( nfields ) )     nfields     = SIZE( clim_type%field_name(:) )
if( present( field_names ) ) field_names = clim_type%field_name

end subroutine query_interpolator
!--lwh
!
!#######################################################################
!
!---------------------------------------------------------------------
!> @brief chomp receives a string from NetCDF files and removes
!!        CHAR(0) from the end of this string.
!!
!! @param [in] <string> The string from the NetCDF file
function chomp(string)
!
! A function to remove CHAR(0) from the end of strings read from NetCDF files.
!
character(len=*), intent(in) :: string
character(len=64) :: chomp

integer :: len

len = len_trim(string)
if (string(len:len) == CHAR(0)) len = len -1

chomp = string(:len)

end function chomp
!
!########################################################################

#include "interpolator_r4.fh"
#include "interpolator_r8.fh"

end module interpolator_mod

!> @}
! close documentation grouping
!
!#######################################################################
