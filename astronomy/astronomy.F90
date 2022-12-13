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
!> @defgroup astronomy_mod astronomy_mod
!> @ingroup astronomy
!> @brief Provides astronomical variables for use
!!        by other modules within fms. The only currently used interface is
!!        for determination of astronomical values needed by the shortwave
!!        radiation packages.
!> @author Fei Liu

!> @addtogroup astronomy_mod
!> @{
module astronomy_mod


use fms_mod,           only: fms_init, mpp_pe, mpp_root_pe, stdlog, write_version_number, &
                             check_nml_error, error_mesg, FATAL, NOTE, WARNING
use time_manager_mod,  only: time_type, set_time, get_time, get_date_julian, set_date_julian, &
                              set_date, length_of_year, time_manager_init, &
                              operator(-), operator(+), operator( // ), operator(<)
use constants_mod,     only: constants_init, PI
use mpp_mod,           only: input_nml_file
use platform_mod,      only: r4_kind, r8_kind

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

! Include variable "version" to be written to log file.
#include<file_version.h>


!---------------------------------------------------------------------
!-----------------------------  interfaces ---------------------------

public :: astronomy_init, get_period, set_period
public :: set_orbital_parameters, get_orbital_parameters
public :: set_ref_date_of_ae, get_ref_date_of_ae
public :: diurnal_solar, daily_mean_solar, annual_mean_solar
public :: astronomy_end, universal_time, orbital_time

!> @}


interface diurnal_solar
   module procedure diurnal_solar_2d_r4, diurnal_solar_2d_r8
   module procedure diurnal_solar_1d_r4, diurnal_solar_1d_r8
   module procedure diurnal_solar_0d_r4, diurnal_solar_0d_r8
   module procedure diurnal_solar_cal_2d_r4, diurnal_solar_cal_2d_r8
   module procedure diurnal_solar_cal_1d_r4, diurnal_solar_cal_1d_r8
   module procedure diurnal_solar_cal_0d_r4, diurnal_solar_cal_0d_r8
end interface diurnal_solar


interface daily_mean_solar
   module procedure daily_mean_solar_2d_r4, daily_mean_solar_2d_r8
   module procedure daily_mean_solar_1d_r4, daily_mean_solar_1d_r8
   module procedure daily_mean_solar_2level_r4, daily_mean_solar_2level_r8 
   module procedure daily_mean_solar_0d_r4, daily_mean_solar_0d_r8
   module procedure daily_mean_solar_cal_2d_r4, daily_mean_solar_cal_2d_r8
   module procedure daily_mean_solar_cal_1d_r4, daily_mean_solar_cal_1d_r8
   module procedure daily_mean_solar_cal_2level_r4, daily_mean_solar_cal_2level_r8 
   module procedure daily_mean_solar_cal_0d_r4, daily_mean_solar_cal_0d_r8
end interface daily_mean_solar


interface annual_mean_solar
   module procedure annual_mean_solar_2d_r4, annual_mean_solar_2d_r8
   module procedure annual_mean_solar_1d_r4, annual_mean_solar_1d_r8
   module procedure annual_mean_solar_2level_r4, annual_mean_solar_2level_r8
end interface annual_mean_solar


interface get_period
   module procedure get_period_time_type_r4, get_period_time_type_r8
   module procedure get_period_integer_r4, get_period_integer_r8 
end interface get_period


interface set_period
   module procedure set_period_time_type_r4, set_period_time_type_r8
   module procedure set_period_integer_r4, set_period_integer_r8 
end interface set_period


private :: orbit         ! Called from astronomy_init and set_orbital_parameters
private :: r_inv_squared ! Called from diurnal_solar, daily_mean_solar and orbit
private :: angle,  declination, half_day ! called from  diurnal_solar and daily_mean_solar
!             half_day, orbital_time, & ! called from  diurnal_solar and daily_mean_solar
!             universal_time ! called from  diurnal_solar:


interface half_day
   module procedure half_day_2d_r4, half_day_2d_r8 
   module procedure half_day_0d_r4, half_day_0d_r8 
end interface half_day

!> @addtogroup astronomy_mod
!> @{

!---------------------------------------------------------------------
!-------- namelist  ---------

real(r8_kind)   :: ecc   = 0.01671   !< Eccentricity of Earth's orbit [dimensionless]
real(r8_kind)   :: obliq = 23.439    !< Obliquity [degrees]
real(r8_kind)   :: per   = 102.932   !< Longitude of perihelion with respect
                                     !! to autumnal equinox in NH [degrees]
integer         :: period = 0        !< Specified length of year [seconds];
                                     !! must be specified to override default
                                     !! value given by length_of_year in
                                     !! time_manager_mod
integer         :: day_ae    = 23    !< Day of specified autumnal equinox
integer         :: month_ae  = 9     !< Month of specified autumnal equinox
integer         :: year_ae   = 1998  !< Year of specified autumnal equinox
integer         :: hour_ae   = 5     !< Hour of specified autumnal equinox
integer         :: minute_ae = 37    !< Minute of specified autumnal equinox
integer         :: second_ae = 0     !< Second of specified autumnal equinox
integer         :: num_angles = 3600 !< Number of intervals into which the year
                                     !! is divided to compute orbital positions


namelist /astronomy_nml/ ecc, obliq, per, period, &
                         year_ae, month_ae,  day_ae,         &
                         hour_ae, minute_ae, second_ae, &
                         num_angles

!--------------------------------------------------------------------
!------   public data ----------


!--------------------------------------------------------------------
!------   private data ----------

type(time_type) :: autumnal_eq_ref  !< time_type variable containing
                                    !! specified time of reference
                                    !! NH autumnal equinox

type(time_type) :: period_time_type !< time_type variable containing
                                    !! period of one orbit

real(r8_kind), dimension(:), allocatable :: orb_angle !< table of orbital positions (0 to
                                             !! 2*pi) as a function of time  used
                                             !! to find actual orbital position
                                             !! via interpolation

real(r8_kind)    :: seconds_per_day=86400.   !< seconds in a day
real(r8_kind)    :: deg_to_rad               !< conversion from degrees to radians
real(r8_kind)    :: twopi                    !< 2 *PI
logical          :: module_is_initialized=.false. !< has the module been initialized ?

real(r8_kind), dimension(:,:), allocatable ::       &
                       cosz_ann, &  !< annual mean cos of zenith angle
                       solar_ann, & !< annual mean solar factor
                       fracday_ann  !< annual mean daylight fraction
real(r8_kind) :: rrsun_ann          !< annual mean earth-sun distance
logical       :: annual_mean_calculated=.false.  !< have the annual mean values been calculated?
integer       :: num_pts = 0              !< count of grid_boxes for which
                                          !! annual mean astronomy values have
                                          !! been calculated
integer       :: total_pts                !< number of grid boxes owned by the processor

!--------------------------------------------------------------------



contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!> @brief astronomy_init is the constructor for astronomy_mod.
!!
!! @throw FATAL, "astronomy_mod ecc must be between 0 and 0.99"
!! @throw FATAL, "astronomy_mod obliquity must be between -90 and 90 degrees"
!! @throw FATAL, "astronomy_mod perihelion must be between 0 and 360 degrees"
subroutine astronomy_init (latb, lonb)

real(r8_kind), dimension(:,:), intent(in), optional :: latb !< 2d array of model latitudes at cell corners [radians]
real(r8_kind), dimension(:,:), intent(in), optional :: lonb !< 2d array of model longitudes at cell corners [radians]

!-------------------------------------------------------------------
!  local variables:
!-------------------------------------------------------------------
integer :: unit, ierr, io, seconds, days, jd, id

!-------------------------------------------------------------------
!    if module has already been initialized, exit.
!-------------------------------------------------------------------
    if (module_is_initialized) return

!-------------------------------------------------------------------
!>    This routine will:
!>    Verify that modules used by this module have been initialized.
!-------------------------------------------------------------------
    call fms_init
    call time_manager_init
    call constants_init

!-----------------------------------------------------------------------
!>    Read namelist.
!-----------------------------------------------------------------------
    read (input_nml_file, astronomy_nml, iostat=io)
    ierr = check_nml_error(io,'astronomy_nml')
!---------------------------------------------------------------------
!>    Write version number and namelist to logfile.
!---------------------------------------------------------------------
    call write_version_number("ASTRONOMY_MOD", version)
    if (mpp_pe() == mpp_root_pe() ) then
       unit = stdlog()
       write (unit, nml=astronomy_nml)
    endif
!--------------------------------------------------------------------
!>    Be sure input values are within valid ranges.
!    QUESTION : ARE THESE THE RIGHT LIMITS ???
!---------------------------------------------------------------------
    if (ecc < real(0.0,kind=r8_kind) .or. ecc > real(0.99,kind=r8_kind)) &
       call error_mesg ('astronomy_mod', &
            'ecc must be between 0 and 0.99', FATAL)
    if (obliq < real(-90.0,kind=r8_kind) .or. obliq > real(90.0,kind=r8_kind)) &
        call error_mesg ('astronomy_mod', &
             'obliquity must be between -90 and 90 degrees', FATAL)
    if (per < real(0.0,kind=r8_kind) .or. per > real(360.0,kind=r8_kind)) &
        call error_mesg ('astronomy_mod', &
             'perihelion must be between 0 and 360 degrees', FATAL)

!----------------------------------------------------------------------
!>    Set up time-type variable defining specified time of autumnal equinox.
!----------------------------------------------------------------------
    autumnal_eq_ref = set_date (year_ae,month_ae,day_ae, &
                                hour_ae,minute_ae,second_ae)

!---------------------------------------------------------------------
!>    Set up time-type variable defining length of year.
!----------------------------------------------------------------------
    if (period == 0) then
        period_time_type = length_of_year()
        call get_time (period_time_type, seconds, days)
            period = int(seconds_per_day*days + seconds)
    else
        period_time_type = set_time(period,0)
    endif

!---------------------------------------------------------------------
!>    Define useful module variables.
!----------------------------------------------------------------------
    twopi = real(2.0,kind=r8_kind)*PI
    deg_to_rad = twopi/real(360.0,kind=r8_kind)

!---------------------------------------------------------------------
!>    Call orbit to define table of orbital angles as function of
!!    orbital time.
!----------------------------------------------------------------------
! wfc moved here from orbit
    allocate ( orb_angle(0:num_angles) )
    call orbit

!--------------------------------------------------------------------
!>    If annual mean radiation is desired, then latb will be present.
!!    allocate arrays to hold the needed astronomical factors. define
!!    the total number of points that the processor is responsible for.
!--------------------------------------------------------------------
    if (present(latb)) then
        jd = size(latb,2) - 1
        id = size(lonb,1) - 1
        allocate (cosz_ann(id, jd))
        allocate (solar_ann(id, jd))
        allocate (fracday_ann(id, jd))
        total_pts = jd*id
    endif

!---------------------------------------------------------------------
!>    Mark the module as initialized.
!---------------------------------------------------------------------
    module_is_initialized=.true.

!---------------------------------------------------------------------

end subroutine astronomy_init


!> @brief get_period_integer returns the length of the year as an
!!        integer number of seconds.
!!
!! @throw FATAL, "astronomy_mod module has not been initialized"
subroutine get_period_integer (period_out)

integer, intent(out) :: period_out !< Length of year [seconds]

!--------------------------------------------------------------------
!   local variables:

integer :: seconds, days

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) &
        call error_mesg ('astronomy_mod module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    define length of year in seconds.
!--------------------------------------------------------------------
    call get_time (period_time_type, seconds, days)
    period_out = int(seconds_per_day*days + seconds)


end subroutine get_period_integer

!> @brief get_period_time_type returns the length of the year as a time_type
!!        variable.
!!
!! @throw FATAL, "astronomy_mod module has not been initialized"
subroutine get_period_time_type (period_out)

type(time_type), intent(inout) :: period_out !< Length of year as time_type variable

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) &
       call error_mesg ('astronomy_mod module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    define length of year as a time_type variable.
!--------------------------------------------------------------------
    period_out = period_time_type

end subroutine get_period_time_type

!> @brief set_period_integer saves as the input length of the year (an
!!        integer) in a time_type module variable.
!!
!! @throw FATAL, "astronomy_mod module has not been initialized"
subroutine set_period_integer (period_in)

integer, intent(in) :: period_in !< Length of year as a time_type

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) &
        call error_mesg ('astronomy_mod module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    define time_type variable defining the length of year from the
!    input value (integer seconds).
!---------------------------------------------------------------------
    period_time_type = set_time(period_in, 0)

end subroutine set_period_integer


!> @brief Set_period_time_type saves the length of the year (input as a
!!        time_type variable) into a time_type module variable.
!!
!! @throw FATAL, "astronomy_mod module has not been initialized"
subroutine set_period_time_type(period_in)

type(time_type), intent(in) :: period_in

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) &
        call error_mesg ('astronomy_mod module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    define time_type variable defining the length of year from the
!    input value (time_type).
!---------------------------------------------------------------------
    period_time_type = period_in


end subroutine set_period_time_type

!> @brief set_ref_date_of_ae provides a means of specifying the reference
!!        date of the NH autumnal equinox for a particular year.
!!
!! @details set_ref_date_of_ae provides a means of specifying the reference
!!          date of the NH autumnal equinox for a particular year.  It is only
!!          used if calls are made to the calandar versions of the routines
!!          diurnal_solar and daily_mean_solar. If the NOLEAP calendar is
!!          used, then the date of autumnal equinox will be the same every
!!          year. If JULIAN is used, then the date of autumnal equinox will
!!          return to the same value every 4th year.
!!
!! @param [in] <day_in> Day of reference autumnal equinox
!! @param [in] <month_in> Month of reference autumnal equinox
!! @param [in] <year_in> Year of reference autumnal equinox
!! @param [out] <second_in> OPTIONAL: Second of reference autumnal equinox
!! @param [out] <minute_in> OPTIONAL: Minute of reference autumnal equinox
!! @param [out] <hour_in> OPTIONAL: Hour of reference autumnal equinox
!!
!! @throw FATAL, "astronomy_mod module has not been initialized"

subroutine set_ref_date_of_ae (day_in,month_in,year_in, &
                               second_in,minute_in,hour_in)

integer, intent(in)           :: day_in, month_in, year_in
integer, intent(in), optional :: second_in, minute_in, hour_in

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) &
        call error_mesg ('astronomy_mod module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    save the input time of ae specification into a time_type module
!    variable autumnal_eq_ref.
!--------------------------------------------------------------------
    day_ae =    day_in
    month_ae =  month_in
    year_ae =   year_in

    if (present(second_in)) then
        second_ae = second_in
        minute_ae = minute_in
        hour_ae =   hour_in
    else
        second_ae = 0
        minute_ae = 0
        hour_ae   = 0
    endif

    autumnal_eq_ref = set_date (year_ae,month_ae,day_ae, &
                                hour_ae,minute_ae,second_ae)

!---------------------------------------------------------------------


end subroutine set_ref_date_of_ae


!> @brief get_ref_date_of_ae retrieves the reference date of the autumnal
!!        equinox as integer variables.
!!
!! @throw FATAL, "astronomy_mod module has not been initialized"
!!
!! @param [out] <day_out> Day of reference autumnal equinox
!! @param [out] <month_out> Month of reference autumnal equinox
!! @param [out] <year_out> Year of reference autumnal equinox
!! @param [out] <second_out> Second of reference autumnal equinox
!! @param [out] <minute_out> Minute of reference autumnal equinox
!! @param [out] <hour_out> Hour of reference autumnal equinox
subroutine get_ref_date_of_ae (day_out,month_out,year_out,&
                               second_out,minute_out,hour_out)

integer, intent(out) :: day_out, month_out, year_out,  &
                        second_out, minute_out, hour_out

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) &
        call error_mesg ('astronomy_mod module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    fill the output fields with the proper module data.
!---------------------------------------------------------------------
    day_out    =  day_ae
    month_out  =  month_ae
    year_out   =  year_ae
    second_out =  second_ae
    minute_out =  minute_ae
    hour_out   =  hour_ae


end subroutine get_ref_date_of_ae

!> @brief astronomy_end is the destructor for astronomy_mod.
subroutine astronomy_end

    !----------------------------------------------------------------------
    !>    check if the module has been initialized.
    !----------------------------------------------------------------------
    if (.not. module_is_initialized)  return
    !                call error_mesg ( 'astronomy_mod',  &
    !                         ' module has not been initialized', FATAL)
    
    !----------------------------------------------------------------------
    !>    deallocate module variables.
    !----------------------------------------------------------------------
    deallocate (orb_angle)
    if (allocated(cosz_ann) ) then
        deallocate (cosz_ann)
        deallocate (fracday_ann)
        deallocate (solar_ann)
    endif
    
    !----------------------------------------------------------------------
    !>    mark the module as uninitialized.
    !----------------------------------------------------------------------
    module_is_initialized = .false.
    
end subroutine astronomy_end

#include "astronomy_r4.fh"
#include "astronomy_r8.fh"

end module astronomy_mod
!> @}
! close documentation grouping
