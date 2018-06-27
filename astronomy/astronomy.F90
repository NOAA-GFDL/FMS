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

!> \author Fei Liu <Fei.Liu@noaa.gov>
!! \link http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/ \endlink
!!
!! \brief astronomy_mod provides astronomical variables for use
!!        by other modules within fms. The only currently used interface is
!!        for determination of astronomical values needed by the shortwave
!!        radiation packages.
!!
!! Modules Included:
!!
!! <table>
!!   <tr>
!!     <th>Module Name</th>
!!     <th>Functions Included</th>
!!   </tr>
!!   <tr>
!!     <td>fms_mod</td>
!!     <td>open_namelist_file, fms_init, mpp_pe, mpp_root_pe, stdlog,
!!         file_exist, write_version_number, check_nml_error, error_mesg,
!!         FATAL, NOTE, WARNING, close_file</td>
!!   </tr>
!!   <tr>
!!     <td>time_manager_mod</td>
!!     <td>time_type, set_time, get_time, get_date_julian, set_date_julian,
!!         set_date, length_of_year, time_manager_init, operator(-),
!!         operator(+), operator( // ), operator(<)</td>
!!   </tr>
!!   <tr>
!!     <td>constants_mod</td>
!!     <td>constants_init, PI</td>
!!   </tr>
!!   <tr>
!!     <td>mpp_mod</td>
!!     <td>input_nml_file</td>
!!   </tr>
!! </table>
                      module astronomy_mod


use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, NOTE, WARNING, close_file
use time_manager_mod,  only: time_type, set_time, get_time, &
                             get_date_julian, set_date_julian, &
                             set_date, length_of_year, &
                             time_manager_init, &
                             operator(-), operator(+), &
                             operator( // ), operator(<)
use constants_mod,     only: constants_init, PI
use mpp_mod,           only: input_nml_file

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

! Include variable "version" to be written to log file.
#include<file_version.h>


!---------------------------------------------------------------------
!-------  interfaces --------

public       &
              astronomy_init, get_period, set_period, &
              set_orbital_parameters, get_orbital_parameters, &
              set_ref_date_of_ae, get_ref_date_of_ae,  &
              diurnal_solar, daily_mean_solar, annual_mean_solar,  &
              astronomy_end, universal_time, orbital_time

!> \page diurnal_solar diurnal_solar Interface
!!
!! ~~~~~~~~~~{.f90}
!! call diurnal_solar (lat, lon, time, cosz, fracday, rrsun, dt_time)
!! call diurnal_solar (lat, lon, gmt, time_since_ae, cosz, fracday, rrsun, dt)
!! ~~~~~~~~~~
!!
!! The first option (used in conjunction with time_manager_mod)
!! generates the real variables gmt and time_since_ae from the
!! time_type input, and then calls diurnal_solar with these real inputs.
!!
!! The time of day is set by
!! ~~~~~~~~~~{.f90}
!! real, intent(in) :: gmt
!! ~~~~~~~~~~
!! The time of year is set by
!! ~~~~~~~~~~{.f90}
!! real, intent(in) :: time_since_ae
!! ~~~~~~~~~~
!! with time_type input, both of these are extracted from
!! ~~~~~~~~~~{.f90}
!! type(time_type), intent(in) :: time
!! ~~~~~~~~~~
!!
!! Separate routines exist within this interface for scalar,
!! 1D or 2D input and output fields:
!!
!! ~~~~~~~~~~{.f90}
!! real, intent(in), dimension(:,:) :: lat, lon
!! real, intent(in), dimension(:)   :: lat, lon
!! real, intent(in)                 :: lat, lon
!!
!! real, intent(out), dimension(:,:) :: cosz, fracday
!! real, intent(out), dimension(:)   :: cosz, fracday
!! real, intent(out)                 :: cosz, fracday
!! ~~~~~~~~~~
!!
!! One may also average the output fields over the time interval
!! between gmt and gmt + dt by including the optional argument dt (or
!! dt_time). dt is measured in radians and must be less than pi
!! (1/2 day). This average is computed analytically, and should be
!! exact except for the fact that changes in earth-sun distance over
!! the time interval dt are ignored. In the context of a diurnal GCM,
!! this option should always be employed to insure that the total flux
!! at the top of the atmosphere is not modified by time truncation error.
!!
!! ~~~~~~~~~~{.f90}
!! real, intent(in), optional :: dt
!! type(time_type), optional :: dt_time
!! ~~~~~~~~~~
!! (see test.90 for examples of the use of these types)
!!
!! \param [in] <lat> Latitudes of model grid points [radians]
!! \param [in] <lon> Longitudes of model grid points [radians]
!! \param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi [radians]
!! \param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi [radians]
!! \param [in] <time> Time at which astronomical values are desired (time_type variable) [seconds, days]
!! \param [out] <cosz> Cosine of solar zenith angle, set to zero when entire period is in darkness [dimensionless]
!! \param [out] <fracday> Daylight fraction of time interval [dimensionless]
!! \param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2 [dimensionless]
!! \param [in] <dt> OPTIONAL: Time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous. [radians], (1 day = 2 * pi)
!! \param [in] <dt_time> OPTIONAL: Time interval after gmt over which the astronomical variables are to be
!!                       averaged. this produces averaged output rather than instantaneous. time_type, [days, seconds]
!! \param [in] <allow_negative_cosz> Allow negative values for cosz?
!! \param [out] <half_day_out> half_day_out
interface diurnal_solar
   module procedure diurnal_solar_2d
   module procedure diurnal_solar_1d
   module procedure diurnal_solar_0d
   module procedure diurnal_solar_cal_2d
   module procedure diurnal_solar_cal_1d
   module procedure diurnal_solar_cal_0d
end interface

!> \page daily_mean_solar daily_mean_solar Interface
!!
!! ~~~~~~~~~~{.f90}
!! call daily_mean_solar (lat, time, cosz, fracday, rrsun)
!! call daily_mean_solar (lat, time_since_ae, cosz, fracday, rrsun)
!! call daily_mean_solar (lat, time, cosz, solar)
!! call daily_mean_solar (lat, time_since_ae, cosz, solar)
!! ~~~~~~~~~~
!!
!! The first option (used in conjunction with time_manager_mod)
!! generates the real variable time_since_ae from the time_type
!! input time, and then calls daily_mean_solar with this real input
!! (option 2). The third and fourth options correspond to the first
!! and second and are used with then spectral 2-layer model, where
!! only cosz and solar are desired as output. These routines generate
!! dummy arguments and then call option 2, where the calculation is done.
!!
!! The time of year is set by
!! ~~~~~~~~~~{.f90}
!!    real, intent(in) :: time_since_ae
!! ~~~~~~~~~~
!! With time_type input, the time of year is extracted from
!! ~~~~~~~~~~{.f90}
!!    type(time_type), intent(in) :: time
!! ~~~~~~~~~~
!!
!! Separate routines exist within this interface for scalar,
!! 1D or 2D input and output fields:
!!
!! ~~~~~~~~~~{.f90}
!! real, intent(in), dimension(:,:) :: lat
!! real, intent(in), dimension(:)   :: lat
!! real, intent(in)                 :: lat
!!
!! real, intent(out), dimension(:,:) :: cosz, fracday
!! real, intent(out), dimension(:)   :: cosz, fracday
!! real, intent(out)                 :: cosz, fracday
!! ~~~~~~~~~~
!!
!! \param [in] <lat> Latitudes of model grid points [radians]
!! \param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi [radians]
!! \param [in] <time> Time at which astronomical values are desired (time_type variable) [seconds, days]
!! \param [out] <cosz> Cosine of solar zenith angle, set to zero when entire period is in darkness [dimensionless]
!! \param [out] <fracday> Daylight fraction of time interval [dimensionless]
!! \param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2 [dimensionless]
!! \param [out] <solar> shortwave flux factor: cosine of zenith angle * daylight fraction / (earth-sun distance squared) [dimensionless]
interface daily_mean_solar
   module procedure daily_mean_solar_2d
   module procedure daily_mean_solar_1d
   module procedure daily_mean_solar_2level
   module procedure daily_mean_solar_0d
   module procedure daily_mean_solar_cal_2d
   module procedure daily_mean_solar_cal_1d
   module procedure daily_mean_solar_cal_2level
   module procedure daily_mean_solar_cal_0d
end interface

!! \page annual_mean_solar annual_mean_solar Interface
!!
!! ~~~~~~~~~~{.f90}
!! call annual_mean_solar (js, je, lat, cosz, solar, fracday, rrsun)
!! call annual_mean_solar (lat, cosz, solar)
!! ~~~~~~~~~~
!!
!! The second interface above is used by the spectral 2-layer model,
!! which requires only cosz and solar as output arguments, and which
!! makes this call during the initialization phase of the model.
!! Separate routines exist within this interface for 1D or 2D input
!! and output fields:
!!
!! ~~~~~~~~~~{.f90}
!! real, intent(in), dimension(:,:) :: lat
!! real, intent(in), dimension(:)   :: lat
!!
!! real, intent(out), dimension(:,:) :: cosz, solar, fracday
!! real, intent(out), dimension(:)   :: cosz, solar, fracday
!! ~~~~~~~~~~
!!
!! \param [in] <jst> Starting subdomain j indices of data in the physics wiondow being integrated
!! \param [in] <jnd> Ending subdomain j indices of data in the physics wiondow being integrated
!! \param [in] <lat> Latitudes of model grid points [radians]
!! \param [out] <cosz> cosz is the average over the year of the cosine of an effective zenith angle
!!                     that would produce the correct daily solar flux if the sun were fixed at that
!!                     single position for the period of daylight on the given day. in this average,
!!                     the daily mean effective cosz is weighted by the daily mean solar flux. [dimensionless]
!! \param [out] <solar> Normalized solar flux, averaged over the year, equal to the product of
!!                      fracday*cosz*rrsun [dimensionless]
!! \param [out] <fracday> Daylight fraction calculated so as to make the average flux (solar) equal to the
!!                        product of the flux-weighted avg cosz * this fracday * assumed annual mean avg
!!                        Earth-Sun distance of 1.0. [dimensionless]
!! \param [out] <rrsun> Annual mean Earth-Sun distance (r) relative to semi-major axis of orbital ellipse
!!                      (a):(a/r)**2 [dimensionless]
interface annual_mean_solar
   module procedure annual_mean_solar_2d
   module procedure annual_mean_solar_1d
   module procedure annual_mean_solar_2level
end interface

!> \page get_period get_period Interface
!!
!! ~~~~~~~~~~{.f90}
!! call get_period (period)
!! ~~~~~~~~~~
!!
!! Separate routines exist within this interface for integer
!! and time_type output:
!!
!! ~~~~~~~~~~{.f90}
!! integer, intent(out)         :: period
!! type(time_type), intent(out) :: period
!! ~~~~~~~~~~
!!
!! \param [out] <period_out> Length of year for calendar in use
interface get_period
   module procedure get_period_time_type, get_period_integer
end interface

!> \page set_period set_period Interface
!!
!! ~~~~~~~~~~{.f90}
!! call set_period (period_in)
!! ~~~~~~~~~~
!!
!! Separate routines exist within this interface for integer
!! and time_type output:
!!
!! ~~~~~~~~~~{.f90}
!! integer, intent(out)         :: period_in
!! type(time_type), intent(out) :: period_in
!! ~~~~~~~~~~
!!
!! \param [in] <period_in> Length of year for calendar in use
interface set_period
   module procedure set_period_time_type, set_period_integer
end interface


private &
              orbit,  & ! Called from astronomy_init and set_orbital_parameters
              r_inv_squared, & ! Called from diurnal_solar, daily_mean_solar and orbit
              angle,  declination, half_day ! called from  diurnal_solar and daily_mean_solar
!             half_day, orbital_time, & ! called from  diurnal_solar and daily_mean_solar
!             universal_time ! called from  diurnal_solar:

!> \page half_day half_day Interface
!!
!! ~~~~~~~~~~{.f90}
!! half_day (latitude, dec) result (h)
!! ~~~~~~~~~~
!!
!! Separate routines exist within this interface for scalar,
!! or 2D input and output fields:
!!
!! ~~~~~~~~~~{.f90}
!! real, intent(in), dimension(:,:) :: latitude
!! real, intent(in)                 :: latitude
!!
!! real, dimension(size(latitude,1),size(latitude,2))  :: h
!! real                                                :: h
!! ~~~~~~~~~~
!!
!! \param [in] <latitude> Latitudes of model grid points [radians]
!! \param [in] <dec> Solar declination [radians]
!! \param [out] <h> Half of the length of daylight at the given latitude and orbital position (dec); value
!!                  ranges between 0 (all darkness) and pi (all daylight) [dimensionless]
interface half_day
   module procedure half_day_2d, half_day_0d
end interface


!---------------------------------------------------------------------
!-------- namelist  ---------

real    :: ecc   = 0.01671   !< Eccentricity of Earth's orbit [dimensionless]
real    :: obliq = 23.439    !< Obliquity [degrees]
real    :: per   = 102.932   !< Longitude of perihelion with respect
                             !! to autumnal equinox in NH [degrees]
integer :: period = 0        !< Specified length of year [seconds];
                             !! must be specified to override default
                             !! value given by length_of_year in
                             !! time_manager_mod
integer :: day_ae    = 23    !< Day of specified autumnal equinox
integer :: month_ae  = 9     !< Month of specified autumnal equinox
integer :: year_ae   = 1998  !< Year of specified autumnal equinox
integer :: hour_ae   = 5     !< Hour of specified autumnal equinox
integer :: minute_ae = 37    !< Minute of specified autumnal equinox
integer :: second_ae = 0     !< Second of specified autumnal equinox
integer :: num_angles = 3600 !< Number of intervals into which the year
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

real, dimension(:), allocatable :: orb_angle !< table of orbital positions (0 to
                                             !! 2*pi) as a function of time  used
                                             !! to find actual orbital position
                                             !! via interpolation

real    :: seconds_per_day=86400.   !< seconds in a day
real    :: deg_to_rad               !< conversion from degrees to radians
real    :: twopi                    !< 2 *PI
logical :: module_is_initialized=.false. !< has the module been initialized ?

real, dimension(:,:), allocatable ::       &
                       cosz_ann, &  !< annual mean cos of zenith angle
                       solar_ann, & !< annual mean solar factor
                       fracday_ann  !< annual mean daylight fraction
real    :: rrsun_ann                !< annual mean earth-sun distance
logical :: annual_mean_calculated=.false.  !< have the annual mean values been calculated?
integer :: num_pts = 0              !< count of grid_boxes for which
                                    !! annual mean astronomy values have
                                    !! been calculated
integer :: total_pts                !< number of grid boxes owned by the processor

!--------------------------------------------------------------------



                           contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!> \brief astronomy_init is the constructor for astronomy_mod.
!!
!! \throw FATAL, "astronomy_mod ecc must be between 0 and 0.99"
!! \throw FATAL, "astronomy_mod obliquity must be between -90 and 90 degrees"
!! \throw FATAL, "astronomy_mod perihelion must be between 0 and 360 degrees"
subroutine astronomy_init (latb, lonb)

real, dimension(:,:), intent(in), optional :: latb !< 2d array of model latitudes at cell corners [radians]
real, dimension(:,:), intent(in), optional :: lonb !< 2d array of model longitudes at cell corners [radians]

!-------------------------------------------------------------------
!  local variables:
!-------------------------------------------------------------------
integer :: unit, ierr, io, seconds, days, jd, id

!-------------------------------------------------------------------
!    if module has already been initialized, exit.
!-------------------------------------------------------------------
      if (module_is_initialized) return

!-------------------------------------------------------------------
!>    Verify that modules used by this module have been initialized.
!-------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call constants_init

!-----------------------------------------------------------------------
!>    Read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, astronomy_nml, iostat=io)
      ierr = check_nml_error(io,'astronomy_nml')
#else
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=astronomy_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'astronomy_nml')
        end do
10      call close_file (unit)
      endif
#endif
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
      if (ecc < 0.0 .or. ecc > 0.99) &
        call error_mesg ('astronomy_mod', &
                      'ecc must be between 0 and 0.99', FATAL)
      if (obliq < -90. .or. obliq > 90.) &
        call error_mesg ('astronomy_mod', &
                  'obliquity must be between -90 and 90 degrees', FATAL)
      if (per <  0.0 .or. per > 360.0) &
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
        period = seconds_per_day*days + seconds
      else
        period_time_type = set_time(period,0)
      endif

!---------------------------------------------------------------------
!>    Define useful module variables.
!----------------------------------------------------------------------
      twopi = 2.*PI
      deg_to_rad = twopi/360.

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


!> \brief get_period_integer returns the length of the year as an
!!        integer number of seconds.
!!
!! \throw FATAL, "astronomy_mod module has not been initialized"
subroutine get_period_integer (period_out)

integer, intent(out) :: period_out !< Length of year [seconds]

!--------------------------------------------------------------------
!   local variables:

      integer :: seconds, days

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    define length of year in seconds.
!--------------------------------------------------------------------
      call get_time (period_time_type, seconds, days)
      period_out = seconds_per_day*days + seconds


end subroutine get_period_integer

!> \brief get_period_time_type returns the length of the year as a time_type
!!        variable.
!!
!! \throw FATAL, "astronomy_mod module has not been initialized"
subroutine get_period_time_type (period_out)

type(time_type), intent(inout) :: period_out !< Length of year as time_type variable

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ( 'astronomy_mod',  &
              ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    define length of year as a time_type variable.
!--------------------------------------------------------------------
      period_out = period_time_type

end subroutine get_period_time_type

!> \brief set_period_integer saves as the input length of the year (an
!!        integer) in a time_type module variable.
!!
!! \throw FATAL, "astronomy_mod module has not been initialized"
subroutine set_period_integer (period_in)

integer, intent(in) :: period_in !< Length of year as a time_type

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    define time_type variable defining the length of year from the
!    input value (integer seconds).
!---------------------------------------------------------------------
      period_time_type = set_time(period_in, 0)

end subroutine set_period_integer


!> \brief Set_period_time_type saves the length of the year (input as a
!!        time_type variable) into a time_type module variable.
!!
!! \throw FATAL, "astronomy_mod module has not been initialized"
subroutine set_period_time_type(period_in)

type(time_type), intent(in) :: period_in

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    define time_type variable defining the length of year from the
!    input value (time_type).
!---------------------------------------------------------------------
      period_time_type = period_in


end subroutine set_period_time_type

!> \brief set_orbital_parameters saves the input values of eccentricity,
!!        obliquity and perihelion time as module variables for use by
!!        astronomy_mod.
!!
!! \throw FATAL, "astronomy_mod module has not been initialized"
!! \throw FATAL, "astronomy_mod ecc must be between 0 and 0.99"
!! \throw FATAL, "astronomy_mod obliquity must be between -90. and 90. degrees"
!! \throw FATAL, "astronomy_mod perihelion must be between 0.0 and 360. degrees"
subroutine set_orbital_parameters (ecc_in, obliq_in, per_in)

real, intent(in) :: ecc_in !< Eccentricity of orbital ellipse [dimensionless]
real, intent(in) :: obliq_in !< Obliquity [degrees]
real, intent(in) :: per_in !< Longitude of perihelion with respect to autumnal
                           !! equinox in northern hemisphere [degrees]

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    be sure input values are within valid ranges.
!    QUESTION : ARE THESE THE RIGHT LIMITS ???
!---------------------------------------------------------------------
      if (ecc_in < 0.0 .or. ecc_in > 0.99) &
        call error_mesg ('astronomy_mod', &
                      'ecc must be between 0 and 0.99', FATAL)
      if (obliq_in < -90.0 .or. obliq > 90.0) &
        call error_mesg ('astronomy_mod', &
                'obliquity must be between -90. and 90. degrees', FATAL)
      if (per_in < 0.0 .or. per_in > 360.0) &
        call error_mesg ('astronomy_mod', &
              'perihelion must be between 0.0 and 360. degrees', FATAL)

!---------------------------------------------------------------------
!    save input values into module variables.
!---------------------------------------------------------------------
      ecc   = ecc_in
      obliq = obliq_in
      per   = per_in

!---------------------------------------------------------------------
!    call orbit to define table of orbital angles as function of
!    orbital time using the input values of parameters just supplied.
!----------------------------------------------------------------------
      call orbit

!----------------------------------------------------------------------

end subroutine set_orbital_parameters

!> \brief get_orbital_parameters retrieves the orbital parameters for use
!!        by another module.
!!
!! \throw FATAL, "astronomy_mod module has not been initialized"
subroutine get_orbital_parameters (ecc_out, obliq_out, per_out)

!-------------------------------------------------------------------
!    get_orbital_parameters retrieves the orbital parameters for use
!    by another module.
!--------------------------------------------------------------------

real, intent(out) :: ecc_out !< Eccentricity of orbital ellipse [dimensionless]
real, intent(out) :: obliq_out !< Obliquity [degrees]
real, intent(out) :: per_out !< Longitude of perihelion with respect to autumnal
                             !! equinox in northern hemisphere [degrees]

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    fill the output arguments with the eccentricity, obliquity and
!    perihelion angle.
!--------------------------------------------------------------------
      ecc_out = ecc
      obliq_out = obliq
      per_out = per


end subroutine get_orbital_parameters


!> \brief set_ref_date_of_ae provides a means of specifying the reference
!!        date of the NH autumnal equinox for a particular year.
!!
!! \details set_ref_date_of_ae provides a means of specifying the reference
!!          date of the NH autumnal equinox for a particular year.  It is only
!!          used if calls are made to the calandar versions of the routines
!!          diurnal_solar and daily_mean_solar. If the NOLEAP calendar is
!!          used, then the date of autumnal equinox will be the same every
!!          year. If JULIAN is used, then the date of autumnal equinox will
!!          return to the same value every 4th year.
!!
!! \param [in] <day_in> Day of reference autumnal equinox
!! \param [in] <month_in> Month of reference autumnal equinox
!! \param [in] <year_in> Year of reference autumnal equinox
!! \param [out] <second_in> OPTIONAL: Second of reference autumnal equinox
!! \param [out] <minute_in> OPTIONAL: Minute of reference autumnal equinox
!! \param [out] <hour_in> OPTIONAL: Hour of reference autumnal equinox
!!
!! \throw FATAL, "astronomy_mod module has not been initialized"
subroutine set_ref_date_of_ae (day_in,month_in,year_in, &
                               second_in,minute_in,hour_in)

integer, intent(in)           :: day_in, month_in, year_in
integer, intent(in), optional :: second_in, minute_in, hour_in

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

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


!> \brief get_ref_date_of_ae retrieves the reference date of the autumnal
!!        equinox as integer variables.
!!
!! \throw FATAL, "astronomy_mod module has not been initialized"
!!
!! \param [out] <day_out> Day of reference autumnal equinox
!! \param [out] <month_out> Month of reference autumnal equinox
!! \param [out] <year_out> Year of reference autumnal equinox
!! \param [out] <second_out> Second of reference autumnal equinox
!! \param [out] <minute_out> Minute of reference autumnal equinox
!! \param [out] <hour_out> Hour of reference autumnal equinox
subroutine get_ref_date_of_ae (day_out,month_out,year_out,&
                               second_out,minute_out,hour_out)

integer, intent(out) :: day_out, month_out, year_out,  &
                        second_out, minute_out, hour_out

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

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


!> \brief diurnal_solar_2d returns 2d fields of cosine of zenith angle,
!!        daylight fraction and earth-sun distance at the specified latitudes,
!!        longitudes and time. These values may be instantaneous or averaged
!!        over a specified time interval.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <lon> Longitudes of model grid points
!! \param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! \param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! \param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! \param [in] <allow_negative_cosz> Allow negative values for cosz?
!! \param [out] <half_day_out> half_day_out
!!
!! \throw FATAL, "astronomy_mod time_since_ae not between 0 and 2pi"
!! \throw FATAL, "astronomy_mod gmt not between 0 and 2pi"
subroutine diurnal_solar_2d (lat, lon, gmt, time_since_ae, cosz, &
                             fracday, rrsun, dt, allow_negative_cosz, &
                             half_day_out)

real, dimension(:,:), intent(in)           :: lat, lon
real,                 intent(in)           :: gmt, time_since_ae
real, dimension(:,:), intent(out)          :: cosz, fracday
real,                 intent(out)          :: rrsun
real,                 intent(in), optional :: dt
logical,              intent(in), optional :: allow_negative_cosz
real, dimension(:,:), intent(out), optional :: half_day_out


!---------------------------------------------------------------------
!   local variables

      real, dimension(size(lat,1),size(lat,2)) :: t, tt, h, aa, bb,  &
                                                  st, stt, sh
      real                                     :: ang, dec
      logical :: Lallow_negative

!---------------------------------------------------------------------
!   local variables
!
!    t           time of day with respect to local noon (2 pi = 1 day)
!                [ radians ]
!    tt          end of averaging period [ radians ]
!    h           half of the daily period of daylight, centered at noon
!                [ radians, -pi --> pi ]
!    aa          sin(lat) * sin(declination)
!    bb          cos(lat) * cos(declination)
!    st          sine of time of day
!    stt         sine of time of day at end of averaging period
!    sh          sine of half-day period
!    ang         position of earth in its orbit wrt autumnal equinox
!                [ radians ]
!    dec         earth's declination [ radians ]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    be sure the time in the annual cycle is legitimate.
!---------------------------------------------------------------------
      if (time_since_ae < 0.0 .or. time_since_ae > twopi) &
          call error_mesg('astronomy_mod', &
                    'time_since_ae not between 0 and 2pi', FATAL)

!--------------------------------------------------------------------
!    be sure the time at longitude = 0.0 is legitimate.
!---------------------------------------------------------------------
      if (gmt < 0.0 .or. gmt > twopi) &
         call error_mesg('astronomy_mod', &
                    'gmt not between 0 and 2pi', FATAL)

!---------------------------------------------------------------------
!>    define the orbital angle (location in year), solar declination and
!!    earth sun distance factor. use functions contained in this module.
!---------------------------------------------------------------------
      ang = angle(time_since_ae)
      dec = declination(ang)
      rrsun  = r_inv_squared(ang)

!---------------------------------------------------------------------
!>    define terms needed in the cosine zenith angle equation.
!--------------------------------------------------------------------
      aa = sin(lat)*sin(dec)
      bb = cos(lat)*cos(dec)

!---------------------------------------------------------------------
!>    define local time. force it to be between -pi and pi.
!--------------------------------------------------------------------
      t = gmt + lon - PI
      where(t >= PI) t = t - twopi
      where(t < -PI) t = t + twopi

      Lallow_negative = .false.
      if (present(allow_negative_cosz)) then
         if (allow_negative_cosz) Lallow_negative = .true.
      end if

!---------------------------------------------------------------------
!>    perform a time integration to obtain cosz and fracday if desired.
!!    output is valid over the period from t to t + dt.
!--------------------------------------------------------------------
      h   = half_day   (lat,dec)

      if ( present(half_day_out) ) then
         half_day_out = h
      end if

      if ( present(dt) ) then   ! (perform time averaging)
        tt = t + dt
        st  = sin(t)
        stt = sin(tt)
        sh  = sin(h)
        cosz = 0.0

        if (.not. Lallow_negative) then
!-------------------------------------------------------------------
!    case 1: entire averaging period is before sunrise.
!-------------------------------------------------------------------
        where (t < -h .and. tt < -h) cosz = 0.0

!-------------------------------------------------------------------
!    case 2: averaging period begins before sunrise, ends after sunrise
!    but before sunset
!-------------------------------------------------------------------
        where ( (tt+h) /= 0.0 .and.   t < -h .and. abs(tt) <= h)   &
             cosz = aa + bb*(stt + sh)/ (tt + h)

!-------------------------------------------------------------------
!    case 3: averaging period begins before sunrise, ends after sunset,
!    but before the next sunrise. modify if averaging period extends
!    past the next day's sunrise, but if averaging period is less than
!    a half- day (pi) that circumstance will never occur.
!-------------------------------------------------------------------
        where (t < -h .and. h /= 0.0 .and. h < tt)    &
              cosz = aa + bb*( sh + sh)/(h+h)

!-------------------------------------------------------------------
!    case 4: averaging period begins after sunrise, ends before sunset.
!-------------------------------------------------------------------
        where ( abs(t) <= h .and. abs(tt) <= h)    &
             cosz = aa + bb*(stt - st)/ (tt - t)

!-------------------------------------------------------------------
!    case 5: averaging period begins after sunrise, ends after sunset.
!    modify when averaging period extends past the next day's sunrise.
!-------------------------------------------------------------------
        where ((h-t) /= 0.0 .and. abs(t) <= h .and.  h < tt)    &
              cosz = aa + bb*(sh - st)/(h-t)

!-------------------------------------------------------------------
!    case 6: averaging period begins after sunrise , ends after the
!    next day's sunrise. note that this includes the case when the
!    day length is one day (h = pi).
!-------------------------------------------------------------------
        where (twopi - h < tt .and. (tt+h-twopi) /= 0.0 .and. t <= h ) &
           cosz = (cosz*(h - t) + (aa*(tt + h - twopi) +     &
            bb*(stt + sh))) / ((h - t) + (tt + h - twopi))

!-------------------------------------------------------------------
!    case 7: averaging period begins after sunset and ends before the
!    next day's sunrise
!-------------------------------------------------------------------
        where(  h <  t .and. twopi - h >= tt  ) cosz = 0.0

!-------------------------------------------------------------------
!    case 8: averaging period begins after sunset and ends after the
!    next day's sunrise but before the next day's sunset. if the
!    averaging period is less than a half-day (pi) the latter
!    circumstance will never occur.
!-----------------------------------------------------------------
        where(  h <  t .and. twopi - h < tt  )
          cosz = aa + bb*(stt + sh) / (tt + h - twopi)
        end where

        else
           cosz = aa + bb*(stt - st)/ (tt - t)
        end if



!-------------------------------------------------------------------
!    day fraction is the fraction of the averaging period contained
!    within the (-h,h) period.
!-------------------------------------------------------------------
        where (t < -h .and.      tt < -h)      fracday = 0.0
        where (t < -h .and. abs(tt) <= h)      fracday = (tt + h )/dt
        where (t < -h .and.       h < tt)      fracday = ( h + h )/dt
        where (abs(t) <= h .and. abs(tt) <= h) fracday = (tt - t )/dt
        where (abs(t) <= h .and.       h < tt) fracday = ( h - t )/dt
        where (      h <  t                 )  fracday = 0.0
        where (twopi - h < tt)                 fracday = fracday +  &
                                                         (tt + h - &
                                                         twopi)/dt
!----------------------------------------------------------------------
!>    if instantaneous values are desired, define cosz at time t.
!----------------------------------------------------------------------
      else  ! (no time averaging)
        if (.not. Lallow_negative) then
           where (abs(t) < h)
             cosz = aa + bb*cos(t)
             fracday = 1.0
           elsewhere
             cosz = 0.0
             fracday = 0.0
           end where
        else
           cosz = aa + bb*cos(t)
           where (abs(t) < h)
             fracday = 1.0
           elsewhere
             fracday = 0.0
           end where
        end if
      end if

!----------------------------------------------------------------------
!>    Check that cosz is not negative, if desired.
!----------------------------------------------------------------------
      if (.not. Lallow_negative) then
         cosz = max(0.0, cosz)
      end if

end subroutine diurnal_solar_2d


!> \brief diurnal_solar_1d takes 1-d input fields, adds a second dimension
!!        and calls diurnal_solar_2d. on return, the 2d fields are returned
!!        to the original 1d fields.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <lon> Longitudes of model grid points
!! \param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! \param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! \param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! \param [in] <allow_negative_cosz> Allow negative values for cosz?
!! \param [out] <half_day_out> half_day_out
subroutine diurnal_solar_1d (lat, lon, gmt, time_since_ae, cosz, &
                             fracday, rrsun, dt, allow_negative_cosz, &
                             half_day_out)

!---------------------------------------------------------------------
real, dimension(:),  intent(in)           :: lat, lon
real,                intent(in)           :: gmt, time_since_ae
real, dimension(:),  intent(out)          :: cosz, fracday
real,                intent(out)          :: rrsun
real,                intent(in), optional :: dt
logical,             intent(in), optional :: allow_negative_cosz
real, dimension(:),  intent(out), optional :: half_day_out

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
      real, dimension(size(lat),1) :: lat_2d, lon_2d, cosz_2d,   &
                                      fracday_2d,halfday_2d

!--------------------------------------------------------------------
!>    define 2-d versions of input data arrays.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat
      lon_2d(:,1) = lon

!--------------------------------------------------------------------
!>    call diurnal_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
!     if (present(dt)) then
        call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae,&
                               cosz_2d, fracday_2d, rrsun, dt=dt, &
                               allow_negative_cosz=allow_negative_cosz, &
                               half_day_out=halfday_2d)
!     else
!       call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae, &
!                              cosz_2d, fracday_2d, rrsun)
!     endif

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(:,1)
      cosz  = cosz_2d (:,1)
      if (present(half_day_out)) then
         half_day_out = halfday_2d(:,1)
      end if

end subroutine diurnal_solar_1d


!> \brief diurnal_solar_0d takes scalar input fields, makes them into 2d
!!        arrays dimensioned (1,1), and calls diurnal_solar_2d. on return,
!!        the 2d fields are converted back to the desired scalar output.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <lon> Longitudes of model grid points
!! \param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! \param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! \param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! \param [in] <allow_negative_cosz> Allow negative values for cosz?
!! \param [out] <half_day_out> half_day_out
subroutine diurnal_solar_0d (lat, lon, gmt, time_since_ae, cosz,  &
                             fracday, rrsun, dt, allow_negative_cosz, &
                             half_day_out)

real, intent(in)           :: lat, lon, gmt, time_since_ae
real, intent(out)          :: cosz, fracday, rrsun
real, intent(in), optional :: dt
logical,intent(in),optional :: allow_negative_cosz
real, intent(out), optional :: half_day_out

!--------------------------------------------------------------------
!  local variables:
!--------------------------------------------------------------------
      real, dimension(1,1) :: lat_2d, lon_2d, cosz_2d, fracday_2d, halfday_2d

!---------------------------------------------------------------------
!>    create 2d arrays from the scalar input fields.
!---------------------------------------------------------------------
      lat_2d = lat
      lon_2d = lon

!--------------------------------------------------------------------
!>    call diurnal_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
!     if (present(dt)) then
        call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae,  &
                               cosz_2d, fracday_2d, rrsun, dt=dt, &
                               allow_negative_cosz=allow_negative_cosz, &
                               half_day_out=halfday_2d)
!     else
!       call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae, &
!                              cosz_2d, fracday_2d, rrsun)
!     end if

!-------------------------------------------------------------------
!>    place output fields into scalars for return to calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(1,1)
      cosz = cosz_2d(1,1)
      if (present(half_day_out)) then
         half_day_out = halfday_2d(1,1)
      end if

end subroutine diurnal_solar_0d


!> \brief diurnal_solar_cal_2d receives time_type inputs, converts
!!        them to real variables and then calls diurnal_solar_2d to
!!        compute desired astronomical variables.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <lon> Longitudes of model grid points
!! \param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! \param [in] <time> Time of year (time_type)
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! \param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! \param [in] <allow_negative_cosz> Allow negative values for cosz?
!! \param [out] <half_day_out> half_day_out
!!
!! \throw FATAL, "astronomy_mod radiation time step must be no longer than 12 hrs"
!! \throw FATAL, "astronomy_mod radiation time step must not be an integral number of days"
subroutine diurnal_solar_cal_2d (lat, lon, time, cosz, fracday,   &
                                 rrsun, dt_time, allow_negative_cosz, &
                                 half_day_out)

!-------------------------------------------------------------------
real, dimension(:,:), intent(in)            :: lat, lon
type(time_type),      intent(in)            :: time
real, dimension(:,:), intent(out)           :: cosz, fracday
real,                 intent(out)           :: rrsun
type(time_type),      intent(in), optional  :: dt_time
logical,              intent(in), optional  :: allow_negative_cosz
real, dimension(:,:), intent(out), optional  :: half_day_out
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
      real :: dt
      real :: gmt, time_since_ae

!---------------------------------------------------------------------
!>    Extract time of day (gmt) from time_type variable time with
!!    function universal_time.
!---------------------------------------------------------------------
      gmt = universal_time(time)

!---------------------------------------------------------------------
!>   Extract the time of year (time_since_ae) from time_type variable
!!   time using the function orbital_time.
!---------------------------------------------------------------------
      time_since_ae = orbital_time(time)

!---------------------------------------------------------------------
!>    Convert optional time_type variable dt_time (length of averaging
!!    period) to a real variable dt with the function universal_time.
!---------------------------------------------------------------------
      if (present(dt_time))  then
        dt = universal_time(dt_time)
        if (dt > PI) then
          call error_mesg ( 'astronomy_mod', &
             'radiation time step must be no longer than 12 hrs', &
                                                          FATAL)
        endif
        if (dt == 0.0) then
          call error_mesg ( 'astronomy_mod', &
              'radiation time step must not be an integral &
                                     &number of days', FATAL)
        endif

!--------------------------------------------------------------------
!>    Call diurnal_solar_2d to calculate astronomy fields, with or
!!    without the optional argument dt.
!--------------------------------------------------------------------
        call diurnal_solar_2d (lat, lon, gmt, time_since_ae, cosz, &
               fracday, rrsun, dt=dt, &
               allow_negative_cosz=allow_negative_cosz, &
               half_day_out=half_day_out)
      else
        call diurnal_solar_2d (lat, lon, gmt, time_since_ae, cosz, &
               fracday, rrsun, &
               allow_negative_cosz=allow_negative_cosz, &
               half_day_out=half_day_out)
      end if

end subroutine diurnal_solar_cal_2d


!> \brief diurnal_solar_cal_1d receives time_type inputs, converts
!!        them to real variables and then calls diurnal_solar_2d to
!!        compute desired astronomical variables.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <lon> Longitudes of model grid points
!! \param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! \param [in] <time> Time of year (time_type)
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! \param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! \param [in] <allow_negative_cosz> Allow negative values for cosz?
!! \param [out] <half_day_out> half_day_out
subroutine diurnal_solar_cal_1d (lat, lon, time, cosz, fracday,   &
                                 rrsun, dt_time, allow_negative_cosz, &
                                 half_day_out)

!--------------------------------------------------------------------
real, dimension(:), intent(in)           :: lat, lon
type(time_type),    intent(in)           :: time
real, dimension(:), intent(out)          :: cosz, fracday
real,               intent(out)          :: rrsun
type(time_type),    intent(in), optional :: dt_time
logical,            intent(in), optional :: allow_negative_cosz
real, dimension(:), intent(out), optional :: half_day_out
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables
!-------------------------------------------------------------------
      real, dimension(size(lat),1) :: lat_2d, lon_2d, cosz_2d, &
                                      fracday_2d, halfday_2d

!--------------------------------------------------------------------
!>    Define 2-d versions of input data arrays.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat
      lon_2d(:,1) = lon

!--------------------------------------------------------------------
!>    Call diurnal_solar_cal_2d to convert the time_types to reals and
!!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      if (present(dt_time)) then
        call diurnal_solar_cal_2d (lat_2d, lon_2d, time, cosz_2d,    &
           fracday_2d, rrsun, dt_time=dt_time, &
           allow_negative_cosz=allow_negative_cosz, &
           half_day_out=halfday_2d)
      else
        call diurnal_solar_cal_2d (lat_2d, lon_2d, time, cosz_2d,    &
           fracday_2d, rrsun, &
           allow_negative_cosz=allow_negative_cosz, &
           half_day_out=halfday_2d)
      end if

!-------------------------------------------------------------------
!>    Place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(:,1)
      cosz  = cosz_2d (:,1)
      if (present(half_day_out)) then
         half_day_out = halfday_2d(:,1)
      end if

end subroutine diurnal_solar_cal_1d


!> \brief diurnal_solar_cal_0d receives time_type inputs, converts them to real variables
!!        and then calls diurnal_solar_2d to compute desired astronomical variables.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <lon> Longitudes of model grid points
!! \param [in] <time> Time of year (time_type)
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a) : (a/r)**2
!! \param [out] <dt_time> OPTIONAL: time interval after gmt over which the astronomical variables are
!!                        to be averaged. this produces averaged output rather than instantaneous.
!! \param [in] <allow_negative_cosz> allow_negative_cosz
!! \param [out] <half_day_out> half_day_out
subroutine diurnal_solar_cal_0d (lat, lon, time, cosz, fracday,   &
                                 rrsun, dt_time, allow_negative_cosz, &
                                 half_day_out)

!---------------------------------------------------------------------
real,            intent(in)           :: lat, lon
type(time_type), intent(in)           :: time
real,            intent(out)          :: cosz, fracday, rrsun
type(time_type), intent(in), optional :: dt_time
logical,         intent(in), optional :: allow_negative_cosz
real,            intent(out), optional :: half_day_out
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
      real, dimension(1,1) :: lat_2d, lon_2d, cosz_2d, fracday_2d, halfday_2d

!--------------------------------------------------------------------
!>    Define 2-d versions of input data arrays.
!--------------------------------------------------------------------
      lat_2d = lat
      lon_2d = lon

!--------------------------------------------------------------------
!>    Call diurnal_solar_cal_2d to convert the time_types to reals and
!!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      if (present(dt_time)) then
        call diurnal_solar_cal_2d (lat_2d, lon_2d, time, cosz_2d,   &
           fracday_2d, rrsun, dt_time=dt_time, &
           allow_negative_cosz=allow_negative_cosz, &
           half_day_out=halfday_2d)
      else
        call diurnal_solar_cal_2d (lat_2d, lon_2d, time, cosz_2d,   &
           fracday_2d, rrsun, &
           allow_negative_cosz=allow_negative_cosz, &
           half_day_out=halfday_2d)
      end if

!-------------------------------------------------------------------
!>    Place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------
      fracday= fracday_2d(1,1)
      cosz = cosz_2d(1,1)
      if (present(half_day_out)) then
         half_day_out = halfday_2d(1,1)
      end if

end subroutine diurnal_solar_cal_0d


!> \brief daily_mean_solar_2d computes the daily mean astronomical parameters for
!!        the input points at latitude lat and time of year time_since_ae.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <h_out> 2-d array of half-day lengths at the latitudes
!! \param [out] <rr_out> the inverse of the square of the earth-sun distance relative
!!                       to the mean distance at angle ang in the earth's orbit.
!!
!! \throw FATAL, "astronomy_mod time_since_ae not between 0 and 2pi"
subroutine daily_mean_solar_2d (lat, time_since_ae, cosz, h_out, rr_out)

!----------------------------------------------------------------------
real, dimension(:,:), intent(in)   :: lat
real,                 intent(in)   :: time_since_ae
real, dimension(:,:), intent(out)  :: cosz, h_out
real,                 intent(out)  :: rr_out
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
      real, dimension(size(lat,1),size(lat,2)) :: h
      real :: ang, dec, rr

!--------------------------------------------------------------------
!    be sure the time in the annual cycle is legitimate.
!---------------------------------------------------------------------
      if (time_since_ae < 0.0 .or. time_since_ae > twopi) &
        call error_mesg('astronomy_mod', &
                        'time_since_ae not between 0 and 2pi', FATAL)

!---------------------------------------------------------------------
!>    Define the orbital angle (location in year), solar declination,
!!    half-day length and earth sun distance factor. Use functions
!!    contained in this module.
!---------------------------------------------------------------------
      ang = angle (time_since_ae)
      dec = declination(ang)
      h   = half_day    (lat, dec)
      rr  = r_inv_squared (ang)

!---------------------------------------------------------------------
!>    Where the entire day is dark, define cosz to be zero. otherwise
!!    use the standard formula. Define the daylight fraction and earth-
!!    sun distance.
!---------------------------------------------------------------------
      where (h == 0.0)
        cosz = 0.0
      elsewhere
        cosz = sin(lat)*sin(dec) + cos(lat)*cos(dec)*sin(h)/h
      end where
      h_out = h/PI
      rr_out = rr

end subroutine daily_mean_solar_2d


!> \brief daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!!        and calls daily_mean_solar_2d. on return, the 2d fields are
!!        returned to the original 1d fields.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <h_out> 2-d array of half-day lengths at the latitudes
!! \param [out] <rr_out> the inverse of the square of the earth-sun distance relative
!!                       to the mean distance at angle ang in the earth's orbit.
subroutine daily_mean_solar_1d (lat, time_since_ae, cosz, h_out, rr_out)

!----------------------------------------------------------------------
real, intent(in), dimension(:) :: lat
real, intent(in) :: time_since_ae
real, intent(out), dimension(size(lat(:))) ::        cosz
real, intent(out), dimension(size(lat(:)))           :: h_out
real, intent(out)           :: rr_out
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables
!----------------------------------------------------------------------
      real, dimension(size(lat),1) :: lat_2d, cosz_2d, hout_2d

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_2d (lat_2d, time_since_ae, cosz_2d,      &
                                hout_2d, rr_out)

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------
      h_out = hout_2d(:,1)
      cosz  = cosz_2d(:,1)

end subroutine daily_mean_solar_1d


!> \brief daily_mean_solar_2level takes 1-d input fields, adds a second
!!        dimension and calls daily_mean_solar_2d. on return, the 2d fields
!!        are returned to the original 1d fields.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <solar> Shortwave flux factor: cosine of zenith angle * daylight fraction / (earth-sun distance squared)
subroutine daily_mean_solar_2level (lat, time_since_ae, cosz, solar)

!----------------------------------------------------------------------
real, intent(in), dimension(:)          :: lat
real, intent(in)                        :: time_since_ae
real, intent(out), dimension(size(lat(:))) :: cosz, solar
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables
!----------------------------------------------------------------------
      real, dimension(size(lat),1) :: lat_2d, cosz_2d, hout_2d
      real                         :: rr_out

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_2d (lat_2d, time_since_ae, cosz_2d,      &
                                hout_2d, rr_out)

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------
      solar = cosz_2d(:,1)*hout_2d(:,1)*rr_out
      cosz  = cosz_2d(:,1)

end subroutine daily_mean_solar_2level


!> \brief daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!!        and calls daily_mean_solar_2d. on return, the 2d fields are
!!        returned to the original 1d fields.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <h_out> 2-d array of half-day lengths at the latitudes
!! \param [out] <rr_out> the inverse of the square of the earth-sun distance relative to
!!                       the mean distance at angle ang in the earth's orbit.
subroutine daily_mean_solar_0d (lat, time_since_ae, cosz, h_out, rr_out)

real, intent(in)         :: lat, time_since_ae
real, intent(out)        :: cosz, h_out, rr_out

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
      real, dimension(1,1) :: lat_2d, cosz_2d, hout_2d

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_2d (lat_2d, time_since_ae, cosz_2d,     &
                                hout_2d, rr_out)

!-------------------------------------------------------------------
!>    return output fields to scalars for return to calling routine.
!-------------------------------------------------------------------
      h_out = hout_2d(1,1)
      cosz  = cosz_2d(1,1)

end subroutine daily_mean_solar_0d


!> \brief daily_mean_solar_cal_2d receives time_type inputs, converts
!!        them to real variables and then calls daily_mean_solar_2d to
!!        compute desired astronomical variables.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <time> Time of year (time_type)
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!!
!! \throw FATAL, "astronomy_mod time_since_ae not between 0 and 2pi"
subroutine daily_mean_solar_cal_2d (lat, time, cosz, fracday, rrsun)

!-------------------------------------------------------------------
real, dimension(:,:), intent(in)  :: lat
type(time_type),      intent(in)  :: time
real, dimension(:,:), intent(out) :: cosz, fracday
real,                 intent(out) :: rrsun
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables
!-------------------------------------------------------------------
      real :: time_since_ae

!--------------------------------------------------------------------
!    be sure the time in the annual cycle is legitimate.
!---------------------------------------------------------------------
      time_since_ae = orbital_time(time)
      if (time_since_ae < 0.0 .or. time_since_ae > twopi) &
          call error_mesg ('astronomy_mod', &
                         'time_since_ae not between 0 and 2pi', FATAL)

!--------------------------------------------------------------------
!    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_2d (lat, time_since_ae, cosz,        &
                                fracday, rrsun)

end subroutine daily_mean_solar_cal_2d


!> \brief daily_mean_solar_cal_1d receives time_type inputs, converts
!!        them to real, 2d variables and then calls daily_mean_solar_2d to
!!        compute desired astronomical variables.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <time> Time of year (time_type)
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
subroutine daily_mean_solar_cal_1d (lat, time, cosz, fracday, rrsun)

real, dimension(:),  intent(in)   :: lat
type(time_type),     intent(in)   :: time
real, dimension(:),  intent(out)  :: cosz, fracday
real,                intent(out)  :: rrsun

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
      real, dimension(size(lat),1) :: lat_2d, cosz_2d, fracday_2d


!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_cal_2d to convert the time_types to reals and
!!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_cal_2d (lat_2d, time, cosz_2d,   &
                                    fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(:,1)
      cosz  = cosz_2d(:,1)

end subroutine daily_mean_solar_cal_1d


!> \brief daily_mean_solar_cal_2level receives 1d arrays and time_type input,
!!        converts them to real, 2d variables and then calls
!!        daily_mean_solar_2d to compute desired astronomical variables.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <time> Time of year (time_type)
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <solar> Shortwave flux factor: cosine of zenith angle * daylight fraction / (earth-sun distance squared)
subroutine daily_mean_solar_cal_2level (lat, time, cosz, solar)

real, dimension(:),  intent(in)   :: lat
type(time_type),     intent(in)   :: time
real, dimension(:),  intent(out)  :: cosz, solar

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
      real, dimension(size(lat),1) :: lat_2d, cosz_2d, fracday_2d
      real                         :: rrsun


!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_cal_2d to convert the time_types to reals and
!!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_cal_2d (lat_2d, time, cosz_2d,   &
                                    fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------
      solar = cosz_2d(:,1)*fracday_2d(:,1)*rrsun
      cosz  = cosz_2d(:,1)

end subroutine daily_mean_solar_cal_2level


!> \brief daily_mean_solar_cal_0d converts scalar input fields to real, 2d variables and
!!        then calls daily_mean_solar_2d to compute desired astronomical variables.
!!
!! \param [in] <lat> Latitudes of model grid points
!! \param [in] <time> Time of year (time_type)
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
subroutine daily_mean_solar_cal_0d (lat, time, cosz, fracday, rrsun)

!--------------------------------------------------------------------
real,             intent(in)  :: lat
type(time_type),  intent(in)  :: time
real,             intent(out) :: cosz, fracday, rrsun
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables
!--------------------------------------------------------------------
      real, dimension(1,1) :: lat_2d, cosz_2d, fracday_2d

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_cal_2d to convert the time_types to reals and
!!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_cal_2d (lat_2d, time, cosz_2d,           &
                                    fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into scalar arguments for return to
!!    calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(1,1)
      cosz  = cosz_2d(1,1)

end subroutine daily_mean_solar_cal_0d


!> \brief annual_mean_solar_2d returns 2d fields of annual mean values of the cosine of
!!        zenith angle, daylight fraction and earth-sun distance at the specified latitude.
!!
!! \param [in] <jst> Starting index of latitude window
!! \param [in] <jnd> Ending index of latitude window
!! \param [in] <lat> Latitudes of model grid points
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <solar> Shortwave flux factor: cosine of zenith angle * daylight fraction / (earth-sun distance squared)
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
subroutine annual_mean_solar_2d (js, je, lat, cosz, solar, fracday,  &
                                 rrsun)

!--------------------------------------------------------------------
integer,                 intent(in)    :: js, je
real, dimension(:,:),    intent(in)    :: lat
real, dimension(:,:),    intent(out)   :: solar, cosz, fracday
real,                    intent(out)   :: rrsun
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables
!--------------------------------------------------------------------
      real, dimension(size(lat,1),size(lat,2)) :: s,z
      real    :: t
      integer :: n, i

!--------------------------------------------------------------------
!    if the calculation has not yet been done, do it here.
!--------------------------------------------------------------------
      if (.not. annual_mean_calculated) then

!----------------------------------------------------------------------
!>    determine annual mean values of solar flux and product of cosz
!!    and solar flux by integrating the annual cycle in num_angles
!!    orbital increments.
!----------------------------------------------------------------------
        solar = 0.0
        cosz = 0.0
        do n =1, num_angles
          t = float(n-1)*twopi/float(num_angles)
          call daily_mean_solar(lat,t, z, fracday, rrsun)
          s = z*rrsun*fracday
          solar = solar + s
          cosz  = cosz  + z*s
        end do
        solar = solar/float(num_angles)
        cosz  = cosz/float(num_angles)

!--------------------------------------------------------------------
!>   define the flux-weighted annual mean cosine of the zenith angle.
!--------------------------------------------------------------------
        where(solar.eq.0.0)
          cosz = 0.0
        elsewhere
          cosz = cosz/solar
        end where

!-------------------------------------------------------------------
!>    define avg fracday such as to make the avg flux (solar) equal to
!!    the product of the avg cosz * avg fracday * assumed mean avg
!!    radius of 1.0. it is unlikely that these avg fracday and avg rr
!!    will ever be used.
!--------------------------------------------------------------------
        where(solar  .eq.0.0)
          fracday = 0.0
        elsewhere
          fracday = solar/cosz
        end where
        rrsun = 1.00

!---------------------------------------------------------------------
!>    save the values that have been calculated as module variables, if
!!    those variables are present; i.e., not the spectral 2-layer model.
!---------------------------------------------------------------------
        if (allocated (cosz_ann)) then
          cosz_ann    = cosz
          solar_ann   = solar
          fracday_ann = fracday
          rrsun_ann = rrsun

!--------------------------------------------------------------------
!>    increment the points computed counter. set flag to end execution
!!    once values have been calculated for all points owned by the
!!    processor.
!--------------------------------------------------------------------
          num_pts = num_pts + size(lat,1)*size(lat,2)
          if ( num_pts == total_pts)  annual_mean_calculated = .true.
        endif

!--------------------------------------------------------------------
!>    if the calculation has been done, return the appropriate module
!!    variables.
!--------------------------------------------------------------------
      else
        if (allocated (cosz_ann)) then
          cosz    = cosz_ann
          solar   = solar_ann
          fracday = fracday_ann
          rrsun = rrsun_ann
        endif
      endif

end subroutine annual_mean_solar_2d


!> \brief annual_mean_solar_1d creates 2-d input fields from 1-d input fields and then calls
!!        annual_mean_solar_2d to obtain 2-d output fields which are then stored into 1-d
!!        fields for return to the calling subroutine.
!!
!! \param [in] <jst> Starting index of latitude window
!! \param [in] <jnd> Ending index of latitude window
!! \param [in] <lat> Latitudes of model grid points
!! \param [out] <cosz> Cosine of solar zenith angle
!! \param [out] <solar> Shortwave flux factor: cosine of zenith angle * daylight fraction / (earth-sun distance squared)
!! \param [out] <fracday> Daylight fraction of time interval
!! \param [out] <rrsun_out> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
subroutine annual_mean_solar_1d (jst, jnd, lat, cosz, solar,  &
                                 fracday, rrsun_out)

!---------------------------------------------------------------------
integer,            intent(in)     :: jst, jnd
real, dimension(:), intent(in)     :: lat(:)
real, dimension(:), intent(out)    :: cosz, solar, fracday
real,               intent(out)    :: rrsun_out
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      real, dimension(size(lat),1) :: lat_2d, solar_2d, cosz_2d,   &
                                      fracday_2d
      real :: rrsun

!--------------------------------------------------------------------
!    if the calculation has not been done, do it here.
!--------------------------------------------------------------------
      if ( .not. annual_mean_calculated) then

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------
        lat_2d(:,1) = lat

!--------------------------------------------------------------------
!>    call annual_mean_solar_2d to calculate the astronomy fields.
!--------------------------------------------------------------------
        call annual_mean_solar_2d (jst, jnd, lat_2d, cosz_2d,   &
                                   solar_2d, fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into 1-D arrays for return to calling routine.
!-------------------------------------------------------------------
        fracday = fracday_2d(:,1)
        rrsun_out = rrsun
        solar = solar_2d(:,1)
        cosz  =  cosz_2d(:,1)

!--------------------------------------------------------------------
!>    if the calculation has been done, simply return the module
!!    variables contain the results at the desired latitudes.
!--------------------------------------------------------------------
      else
        cosz(:)    = cosz_ann(1,jst:jnd)
        solar(:)   = solar_ann(1,jst:jnd)
        fracday(:) = fracday_ann(1,jst:jnd)
        rrsun      = rrsun_ann
      endif

end subroutine annual_mean_solar_1d


!> \brief annual_mean_solar_2level creates 2-d input fields from 1-d input fields
!!        and then calls annual_mean_solar_2d to obtain 2-d output fields which are
!!        then stored into 1-d fields for return to the calling subroutine. This
!!        subroutine will be called during model initialization.
!!
!! \throw FATAL, "astronomy_mod annual_mean_solar_2level should be called only once"
subroutine annual_mean_solar_2level (lat, cosz, solar)

!---------------------------------------------------------------------
real, dimension(:), intent(in)     :: lat !< Latitudes of model grid points
real, dimension(:), intent(out)    :: cosz !< Cosine of solar zenith angle
real, dimension(:), intent(out)    :: solar !< shortwave flux factor: cosine of zenith angle *
                                            !! daylight fraction / (earth-sun distance squared)
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      real, dimension(size(lat),1) :: lat_2d, solar_2d, cosz_2d,   &
                                      fracday_2d
      integer :: jst, jnd
      real    :: rrsun

!--------------------------------------------------------------------
!    if the calculation has not been done, do it here.
!--------------------------------------------------------------------
      if ( .not. annual_mean_calculated) then

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------
        lat_2d(:,1) = lat
        jst = 1
        jnd = size(lat(:))

!--------------------------------------------------------------------
!>    call annual_mean_solar_2d to calculate the astronomy fields.
!--------------------------------------------------------------------
        call annual_mean_solar_2d (jst, jnd, lat_2d, cosz_2d,   &
                                   solar_2d, fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into 1-D arrays for return to calling routine.
!-------------------------------------------------------------------
        solar = solar_2d(:,1)
        cosz  =  cosz_2d(:,1)

!--------------------------------------------------------------------
!>    if the calculation has been done, print an error message since
!!    this subroutine should be called only once.
!--------------------------------------------------------------------
      else
        call error_mesg ('astronomy_mod', &
            'annual_mean_solar_2level should be called only once', &
                                                                 FATAL)
      endif
      annual_mean_calculated = .true.

end subroutine annual_mean_solar_2level


!> \brief astronomy_end is the destructor for astronomy_mod.
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


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!> \brief Orbit computes and stores a table of value of orbital angles as a
!!        function of orbital time (both the angle and time are zero at
!!        autumnal equinox in the NH, and range from 0 to 2*pi).
subroutine orbit

!---------------------------------------------------------------------
!   local variables

      integer :: n
      real    :: d1, d2, d3, d4, d5, dt, norm

!--------------------------------------------------------------------
!>    allocate the orbital angle array, sized by the namelist parameter
!!    num_angles, defining the annual cycle resolution of the earth's
!!    orbit. define some constants to be used.
!--------------------------------------------------------------------
! wfc moving to astronomy_init
!     allocate ( orb_angle(0:num_angles) )
      orb_angle(0) = 0.0
      dt = twopi/float(num_angles)
      norm = sqrt(1.0 - ecc**2)
      dt = dt*norm

!---------------------------------------------------------------------
!>    define the orbital angle at each of the num_angles locations in
!!    the orbit.
!---------------------------------------------------------------------
      do n = 1, num_angles
        d1 = dt*r_inv_squared(orb_angle(n-1))
        d2 = dt*r_inv_squared(orb_angle(n-1)+0.5*d1)
        d3 = dt*r_inv_squared(orb_angle(n-1)+0.5*d2)
        d4 = dt*r_inv_squared(orb_angle(n-1)+d3)
        d5 = d1/6.0 + d2/3.0 +d3/3.0 +d4/6.0
        orb_angle(n) = orb_angle(n-1) + d5
      end do

end subroutine orbit


!> \brief r_inv_squared returns the inverse of the square of the earth-sun
!!        distance relative to the mean distance at angle ang in the Earth's orbit.
function r_inv_squared (ang)

!--------------------------------------------------------------------
real, intent(in) :: ang !< angular position of earth in its orbit, relative to a
                        !! value of 0.0 at the NH autumnal equinox, value between
                        !! 0.0 and 2 * pi [radians]
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

real :: r_inv_squared !< The inverse of the square of the earth-sun distance relative
                      !! to the mean distance [dimensionless]
real :: r             !< Earth-Sun distance relative to mean distance [dimensionless]
real :: rad_per       !< Angular position of perihelion [radians]

!--------------------------------------------------------------------
!>    define the earth-sun distance (r) and then return the inverse of
!!    its square (r_inv_squared) to the calling routine.
!--------------------------------------------------------------------
      rad_per       = per*deg_to_rad
      r             = (1. - ecc**2)/(1. + ecc*cos(ang - rad_per))
      r_inv_squared = r**(-2)


end function r_inv_squared


!> \brief angle determines the position within the earth's orbit at time t
!!        in the year (t = 0 at NH autumnal equinox) by interpolating
!!        into the orbital position table.
function angle (t)

!--------------------------------------------------------------------
real, intent(in) :: t !< time of year (between 0 and 2*pi; t=0 at NH autumnal equinox
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
      real :: angle !< Orbital position relative to NH autumnal equinox [radians]
      real :: norm_time !< Index into orbital table corresponding to input time [dimensionless]
      real :: x !< Fractional distance between the orbital table entries bracketing the input time [dimensionless]
      integer :: int !< Table index which is lower than actual position, but closest to it [dimensionless]
      integer :: int_1 !< Next table index just larger than actual orbital position [dimensionless]

!--------------------------------------------------------------------
!>    Define orbital tables indices bracketing current orbital time
!!    (int and int_1). Define table index distance between the lower
!!    table value (int) and the actual orbital time (x). Define orbital
!!    position as being x of the way between int and int_1. Renormalize
!!    angle to be within the range 0 to 2*pi.
!--------------------------------------------------------------------
      norm_time = t*float(num_angles)/twopi
      int = floor(norm_time)
      int = modulo(int,num_angles)
      int_1 = int+1
      x = norm_time - floor(norm_time)
      angle = (1.0 -x)*orb_angle(int) + x*orb_angle(int_1)
      angle = modulo(angle, twopi)

end function angle


!> \brief Declination returns the solar declination angle at orbital
!!        position ang in earth's orbit.
function declination (ang)

!--------------------------------------------------------------------
real, intent(in) :: ang !< solar orbital position ang in earth's orbit
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables

      real :: declination !< Solar declination angle [radians]
      real :: rad_obliq !< Obliquity of the ecliptic [radians]
      real :: sin_dec !< Sine of the solar declination [dimensionless]

!---------------------------------------------------------------------
!    compute the solar declination.
!---------------------------------------------------------------------
      rad_obliq   =   obliq*deg_to_rad
      sin_dec     = - sin(rad_obliq)*sin(ang)
      declination =   asin(sin_dec)

end function declination


!> \brief half_day_2d returns a 2-d array of half-day lengths at the
!!        latitudes and declination provided.
!!
function half_day_2d (latitude, dec) result(h)

!---------------------------------------------------------------------
real, dimension(:,:), intent(in)                     :: latitude !< Latitutde of view point
real,                 intent(in)                     :: dec !< Solar declination angle at view point
real, dimension(size(latitude,1),size(latitude,2))   :: h
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
      real, dimension (size(latitude,1),size(latitude,2)):: &
                                                  cos_half_day, & !< Cosine of half-day length [dimensionless]
                                                                                                  lat !< Model latitude, adjusted so that it is never 0.5*pi or -0.5*pi
      real :: tan_dec !< tangent of solar declination [dimensionless]
      real :: eps = 1.0E-05 !< small increment

!--------------------------------------------------------------------
!>    define tangent of the declination.
!--------------------------------------------------------------------
      tan_dec = tan(dec)

!--------------------------------------------------------------------
!>    adjust latitude so that its tangent will be defined.
!--------------------------------------------------------------------
      lat = latitude
      where (latitude ==  0.5*PI) lat= latitude - eps
      where (latitude == -0.5*PI) lat= latitude + eps

!--------------------------------------------------------------------
!>    define the cosine of the half-day length. adjust for cases of
!!    all daylight or all night.
!--------------------------------------------------------------------
      cos_half_day = -tan(lat)*tan_dec
      where (cos_half_day <= -1.0)  h = PI
      where (cos_half_day >= +1.0)  h = 0.0
      where(cos_half_day > -1.0 .and. cos_half_day < 1.0) &
                                               h = acos(cos_half_day)

end function half_day_2d


!> \brief half_day_0d takes scalar input fields, makes them into 2-d fields
!!        dimensioned (1,1), and calls half_day_2d. On return, the 2-d
!!        fields are converted to the desired scalar output.
!!
!! \param [in] <latitude> Latitutde of view point
!! \param [in] <dec> Solar declination angle at view point
function half_day_0d(latitude, dec) result(h)

real, intent(in) :: latitude, dec
real             :: h

!----------------------------------------------------------------------
!  local variables
!----------------------------------------------------------------------
      real, dimension(1,1) :: lat_2d, h_2d

!---------------------------------------------------------------------
!    create 2d array from the input latitude field.
!---------------------------------------------------------------------
      lat_2d = latitude

!---------------------------------------------------------------------
!    call half_day with the 2d arguments to calculate half-day length.
!---------------------------------------------------------------------
      h_2d = half_day (lat_2d, dec)

!---------------------------------------------------------------------
!    create scalar from 2d array.
!---------------------------------------------------------------------
      h = h_2d(1,1)

end function half_day_0d


!> \brief Orbital time returns the time (1 year = 2*pi) since autumnal
!!        equinox
!!
!! \details Orbital time returns the time (1 year = 2*pi) since autumnal
!!          equinox; autumnal_eq_ref is a module variable of time_type and
!!          will have been defined by default or by a call to
!!          set_ref_date_of_ae; length_of_year is available through the time
!!          manager and is set at the value approriate for the calandar being used
function orbital_time(time) result(t)

type(time_type), intent(in) :: time !< time (1 year = 2*pi) since autumnal equinox
real                        :: t

      t = real ( (time - autumnal_eq_ref)//period_time_type)
      t = twopi*(t - floor(t))
      if (time < autumnal_eq_ref) t = twopi - t


end function orbital_time


!> \brief universal_time returns the time of day at longitude = 0.0
!!        (1 day = 2*pi)
function universal_time(time) result(t)

type(time_type), intent(in) :: time !< Time (1 year = 2*pi) since autumnal equinox
real                        :: t

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
      integer ::  seconds, days

      call get_time (time, seconds, days)
      t = twopi*real(seconds)/seconds_per_day

end function universal_time


                   end module astronomy_mod
