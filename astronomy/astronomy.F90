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


    use fms_mod,           only: fms_init, &
                                 mpp_pe, mpp_root_pe, stdlog, &
                                 write_version_number, &
                                 check_nml_error, error_mesg, &
                                 FATAL, NOTE, WARNING
    use time_manager_mod,  only: time_type, set_time, get_time, &
                                 get_date_julian, set_date_julian, &
                                 set_date, length_of_year, &
                                 time_manager_init, &
                                 operator(-), operator(+), &
                                 operator( // ), operator(<)
    use fms_io_utils_mod,       only: get_data_type_string
    use constants_mod,     only: constants_init, PI
    use mpp_mod,           only: input_nml_file
    use platform_mod,      only: r4_kind, r8_kind, i4_kind, i8_kind

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

    interface set_orbital_parameters
      module procedure set_orbital_parameters_r4, set_orbital_parameters_r8
    end interface set_orbital_parameters

    interface get_orbital_parameters
        module procedure get_orbital_parameters_r4, get_orbital_parameters_r8
    end interface get_orbital_parameters

    !> @}

    !> @brief Calculates solar information for the given location(lat & lon) and time
    !!
    !> ~~~~~~~~~~{.f90}
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
    !!
    !! @param [in] <lat> Latitudes of model grid points [radians]
    !! @param [in] <lon> Longitudes of model grid points [radians]
    !! @param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi [radians]
    !! @param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi [radians]
    !! @param [in] <time> Time at which astronomical values are desired (time_type variable) [seconds, days]
    !! @param [out] <cosz> Cosine of solar zenith angle, set to zero when entire period is in darkness [dimensionless]
    !! @param [out] <fracday> Daylight fraction of time interval [dimensionless]
    !! @param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse
    !! (a):(a/r)**2 [dimensionless]
    !! @param [in] <dt> OPTIONAL: Time interval after gmt over which the astronomical variables are to be
    !!                  averaged. this produces averaged output rather than instantaneous. [radians], (1 day = 2 * pi)
    !! @param [in] <dt_time> OPTIONAL: Time interval after gmt over which the astronomical variables are to be
    !!                       averaged. this produces averaged output rather than instantaneous. time_type,
    !!                       [days, seconds]
    !! @param [in] <allow_negative_cosz> Allow negative values for cosz?
    !! @param [out] <half_day_out> half_day_out
    !> @ingroup astronomy_mod
    interface diurnal_solar
      module procedure diurnal_solar_2d_r4, diurnal_solar_2d_r8
      module procedure diurnal_solar_1d_r4, diurnal_solar_1d_r8
      module procedure diurnal_solar_0d_r4, diurnal_solar_0d_r8
      module procedure diurnal_solar_cal_2d_r4, diurnal_solar_cal_2d_r8
      module procedure diurnal_solar_cal_1d_r4, diurnal_solar_cal_1d_r8
      module procedure diurnal_solar_cal_0d_r4, diurnal_solar_cal_0d_r8
    end interface diurnal_solar

    !> @brief Calculates the daily mean solar information for a given time and latitude.
    !!
    !> ~~~~~~~~~~{.f90}
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
    !! @param [in] <lat> Latitudes of model grid points [radians]
    !! @param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi [radians]
    !! @param [in] <time> Time at which astronomical values are desired (time_type variable) [seconds, days]
    !! @param [out] <cosz> Cosine of solar zenith angle, set to zero when entire period is in darkness [dimensionless]
    !! @param [out] <fracday> Daylight fraction of time interval [dimensionless]
    !! @param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse
    !! (a):(a/r)**2 [dimensionless]
    !! @param [out] <solar> shortwave flux factor: cosine of zenith angle * daylight fraction /
    !! (earth-sun distance squared) [dimensionless]
    !> @ingroup astronomy_mod
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

    !> Calculates the annual mean of solar information for a given latitude and time.
    !!
    !> ~~~~~~~~~~{.f90}
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
    !! @param [in] <jst> Starting subdomain j indices of data in the physics wiondow being integrated
    !! @param [in] <jnd> Ending subdomain j indices of data in the physics wiondow being integrated
    !! @param [in] <lat> Latitudes of model grid points [radians]
    !! @param [out] <cosz> cosz is the average over the year of the cosine of an effective zenith angle
    !!                     that would produce the correct daily solar flux if the sun were fixed at that
    !!                     single position for the period of daylight on the given day. in this average,
    !!                     the daily mean effective cosz is weighted by the daily mean solar flux. [dimensionless]
    !! @param [out] <solar> Normalized solar flux, averaged over the year, equal to the product of
    !!                      fracday*cosz*rrsun [dimensionless]
    !! @param [out] <fracday> Daylight fraction calculated so as to make the average flux (solar) equal to the
    !!                        product of the flux-weighted avg cosz * this fracday * assumed annual mean avg
    !!                        Earth-Sun distance of 1.0. [dimensionless]
    !! @param [out] <rrsun> Annual mean Earth-Sun distance (r) relative to semi-major axis of orbital ellipse
    !!                      (a):(a/r)**2 [dimensionless]
    !> @ingroup astronomy_mod
    interface annual_mean_solar
      module procedure annual_mean_solar_2d_r4, annual_mean_solar_2d_r8
      module procedure annual_mean_solar_1d_r4, annual_mean_solar_1d_r8
      module procedure annual_mean_solar_2level_r4, annual_mean_solar_2level_r8
    end interface annual_mean_solar

    !> Gets the length of year for current calendar
    !!
    !> ~~~~~~~~~~{.f90}
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
    !! @param [out] <period_out> Length of year for calendar in use
    !> @ingroup astronomy_mod
    interface get_period
       module procedure get_period_time_type, get_period_integer
    end interface

    !> Sets the length of a year for the calendar in use
    !!
    !> ~~~~~~~~~~{.f90}
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
    !! @param [in] <period_in> Length of year for calendar in use
    !> @ingroup astronomy_mod
    interface set_period
       module procedure set_period_time_type, set_period_integer
    end interface


    private &
                  orbit,  & ! Called from astronomy_init and set_orbital_parameters
                  r_inv_squared, & ! Called from diurnal_solar, daily_mean_solar and orbit
                  angle,  declination, half_day ! called from  diurnal_solar and daily_mean_solar
    !             half_day, orbital_time, & ! called from  diurnal_solar and daily_mean_solar
    !             universal_time ! called from  diurnal_solar:

    interface r_inv_squared
      module procedure r_inv_squared_r4, r_inv_squared_r8
    end interface r_inv_squared

    interface angle
      module procedure angle_r4, angle_r8
    end interface angle

    interface declination
      module procedure declination_r4, declination_r8
    end interface declination

    !> Private interface for internal use by dirunal_solar and daily_mean_solar.
    !!
    !> Example usage:
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
    !! @param [in] <latitude> Latitudes of model grid points [radians]
    !! @param [in] <dec> Solar declination [radians]
    !! @param [out] <h> Half of the length of daylight at the given latitude and orbital position (dec); value
    !!                  ranges between 0 (all darkness) and pi (all daylight) [dimensionless]
    !> @ingroup astronomy_mod
    interface half_day
      module procedure half_day_2d_r4, half_day_2d_r8
      module procedure half_day_0d_r4, half_day_0d_r8
    end interface half_day


!> @addtogroup astronomy_mod
!> @{

!---------------------------------------------------------------------
!-------- namelist  ---------

real(r8_kind)   :: ecc   = 0.01671_r8_kind  !< Eccentricity of Earth's orbit [dimensionless]
real(r8_kind)   :: obliq = 23.439_r8_kind   !< Obliquity [degrees]
real(r8_kind)   :: per   = 102.932_r8_kind  !< Longitude of perihelion with respect
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

real(r8_kind)    :: seconds_per_day = 86400.0_r8_kind !< seconds in a day
real(r8_kind)    :: deg_to_rad                        !< conversion from degrees to radians
real(r8_kind)    :: twopi                             !< 2 *PI
logical          :: module_is_initialized=.false.     !< has the module been initialized ?

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

class(*), dimension(:,:), intent(in), optional :: latb !< 2d array of model latitudes at cell corners [radians]
class(*), dimension(:,:), intent(in), optional :: lonb !< 2d array of model longitudes at cell corners [radians]
logical :: is_valid

!-------------------------------------------------------------------
!  local variables:
!-------------------------------------------------------------------
integer :: iunit, ierr, io, seconds, days, jd, id
character(len=17) :: err_str

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
       iunit = stdlog()
       write (iunit, nml=astronomy_nml)
    endif
!--------------------------------------------------------------------
!>    Be sure input values are within valid ranges.
!    QUESTION : ARE THESE THE RIGHT LIMITS ???
!---------------------------------------------------------------------
    if (ecc < 0.0_r8_kind .or. ecc > 0.99_r8_kind) &
       call error_mesg ('astronomy_mod', &
            'ecc must be between 0 and 0.99', FATAL)
    if (obliq < -90.0_r8_kind .or. obliq > 90.0_r8_kind) &
        call error_mesg ('astronomy_mod', &
             'obliquity must be between -90 and 90 degrees', FATAL)
    if (per < 0.0_r8_kind .or. per > 360.0_r8_kind) &
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
            period = int(seconds_per_day) * days + seconds
    else
        period_time_type = set_time(period,0)
    endif

!---------------------------------------------------------------------
!>    Define useful module variables.
!----------------------------------------------------------------------
    twopi = 2.0_r8_kind * PI
    deg_to_rad = twopi/360.0_r8_kind

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
    ! check that no invalid types (integers or characters) are given as optional arg

    is_valid = .false.
    if (present(latb) .and. present(lonb)) then
       select type (latb)
       type is (real(r4_kind))
          select type (lonb)
          type is (real(r4_kind))
             is_valid = .true.
          class default
             call get_data_type_string(lonb, err_str)
             call error_mesg('astronomy_mod', 'kind mismatch, argument latb is real(r4_kind) but lonb has type: '// &
                                               err_str, FATAL)
          end select
       type is (real(r8_kind))
          select type (lonb)
          type is (real(r8_kind))
             is_valid = .true.
          class default
             call get_data_type_string(lonb, err_str)
             call error_mesg('astronomy_mod', 'kind mismatch, argument latb is real(r8_kind) but lonb has type: '//&
                                               err_str, FATAL)
          end select
       end select
       if( is_valid ) then
          jd = size(latb,2) - 1
          id = size(lonb,1) - 1
          allocate (cosz_ann(id, jd))
          allocate (solar_ann(id, jd))
          allocate (fracday_ann(id, jd))
          total_pts = jd*id
       else
          call error_mesg('astronomy_mod', 'latb has unsupported kind size.' // &
                          'latb and lonb should both be real(r4_kind) or real(r8_kind)', FATAL)
       end if
    elseif ( (present(latb) .and. .not. present(lonb)) .or. (present(lonb) .and. .not. present(latb)) ) then
       call error_mesg ('astronomy_mod', 'lat and lon must both be present', FATAL)
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
    if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod', ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    define length of year in seconds.
!--------------------------------------------------------------------
    call get_time (period_time_type, seconds, days)
    period_out = int(seconds_per_day) * days + seconds


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
       call error_mesg ('astronomy_mod', 'module has not been initialized', FATAL)

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
        call error_mesg ('astronomy_mod', 'module has not been initialized', FATAL)

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
        call error_mesg ('astronomy_mod', 'module has not been initialized', FATAL)

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
        call error_mesg ('astronomy_mod', 'module has not been initialized', FATAL)

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
        call error_mesg ('astronomy_mod', 'module has not been initialized', FATAL)

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
    if (allocated(orb_angle)) deallocate (orb_angle)
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

!> @brief Orbit computes and stores a table of value of orbital angles as a
!!        function of orbital time (both the angle and time are zero at
!!        autumnal equinox in the NH, and range from 0 to 2*pi).
subroutine orbit
!---------------------------------------------------------------------
!   local variables

integer                  :: n
real(kind=r8_kind) :: d1, d2, d3, d4, d5, dt, norm

!--------------------------------------------------------------------
!>    allocate the orbital angle array, sized by the namelist parameter
!!    num_angles, defining the annual cycle resolution of the earth's
!!    orbit. define some constants to be used.
!--------------------------------------------------------------------
! wfc moving to astronomy_init
!     allocate ( orb_angle(0:num_angles) )
orb_angle(0) = 0.0_r8_kind
dt = twopi/real(num_angles, r8_kind)
norm = sqrt(1.0_r8_kind - ecc**2)
dt = dt*norm

!---------------------------------------------------------------------
!>    define the orbital angle at each of the num_angles locations in
!!    the orbit.
!---------------------------------------------------------------------
    do n = 1, num_angles
       d1 = dt*r_inv_squared(orb_angle(n-1))
       d2 = dt*r_inv_squared(orb_angle(n-1) + 0.5_r8_kind * d1)
       d3 = dt*r_inv_squared(orb_angle(n-1) + 0.5_r8_kind * d2)
       d4 = dt*r_inv_squared(orb_angle(n-1) + d3)
       d5 = d1/6.0_r8_kind + d2/3.0_r8_kind + d3/3.0_r8_kind + d4/6.0_r8_kind
        orb_angle(n) = orb_angle(n-1) + d5
    end do

end subroutine orbit


!> @brief Orbital time returns the time (1 year = 2*pi) since autumnal
!!        equinox
!!
!! @details Orbital time returns the time (1 year = 2*pi) since autumnal
!!          equinox; autumnal_eq_ref is a module variable of time_type and
!!          will have been defined by default or by a call to
!!          set_ref_date_of_ae; length_of_year is available through the time
!!          manager and is set at the value approriate for the calandar being used
function orbital_time(time) result(t)

type(time_type), intent(in) :: time !< time (1 year = 2*pi) since autumnal equinox
real(kind=r8_kind)    :: t

    t = (time - autumnal_eq_ref)//period_time_type
    t = twopi*(t - real(floor(t), r8_kind))
    if (time < autumnal_eq_ref) t = twopi - t

end function orbital_time


!> @brief universal_time returns the time of day at longitude = 0.0
!!       (1 day = 2*pi)
function universal_time(time) result(t)

type(time_type), intent(in) :: time !< Time (1 year = 2*pi) since autumnal equinox
real(kind=r8_kind)    :: t

    !--------------------------------------------------------------------
    !   local variables
    !--------------------------------------------------------------------
    integer ::  seconds, days

    call get_time(time, seconds, days)
        t = twopi * real(seconds, r8_kind)/real(seconds_per_day, r8_kind)

end function universal_time

#include "astronomy_r4.fh"
#include "astronomy_r8.fh"

end module astronomy_mod
!> @}
! close documentation grouping
