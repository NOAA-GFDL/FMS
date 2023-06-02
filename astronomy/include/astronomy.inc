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

!> @brief set_orbital_parameters saves the input values of eccentricity,
!!        obliquity and perihelion time as module variables for use by
!!        astronomy_mod.
!!
!! @throw FATAL, "astronomy_mod module has not been initialized"
!! @throw FATAL, "astronomy_mod ecc must be between 0 and 0.99"
!! @throw FATAL, "astronomy_mod obliquity must be between -90. and 90. degrees"
!! @throw FATAL, "astronomy_mod perihelion must be between 0.0 and 360. degrees"

subroutine SET_ORBITAL_PARAMETERS_(ecc_in, obliq_in, per_in)

real(kind=FMS_AST_KIND_), intent(in) :: ecc_in !< Eccentricity of orbital ellipse [dimensionless]
real(kind=FMS_AST_KIND_), intent(in) :: obliq_in !< Obliquity [degrees]
real(kind=FMS_AST_KIND_), intent(in) :: per_in !< Longitude of perihelion with respect to autumnal
                           !! equinox in northern hemisphere [degrees]

integer, parameter                   :: lkind = FMS_AST_KIND_

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) &
        call error_mesg('astronomy_mod', 'module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    be sure input values are within valid ranges.
!    QUESTION : ARE THESE THE RIGHT LIMITS ???
!---------------------------------------------------------------------

    if (ecc_in < 0.0_lkind .or. ecc_in > 0.99_lkind) &
        call error_mesg('astronomy_mod', 'ecc must be between 0 and 0.99', FATAL)
    if (obliq_in < -90.0_lkind .or. real(obliq, FMS_AST_KIND_) > 90.0_lkind) &
        call error_mesg('astronomy_mod', 'obliquity must be between -90. and 90. degrees', FATAL)
    if (per_in < 0.0_lkind .or. per_in > 360.0_lkind) &
        call error_mesg('astronomy_mod', 'perihelion must be between 0.0 and 360. degrees', FATAL)

!---------------------------------------------------------------------
!    save input values into module variables.
!---------------------------------------------------------------------

    ecc   = real(ecc_in,   FMS_AST_KIND_)
    obliq = real(obliq_in, FMS_AST_KIND_)
    per   = real(per_in,   FMS_AST_KIND_)

!---------------------------------------------------------------------
!    call orbit to define table of orbital angles as function of
!    orbital time using the input values of parameters just supplied.
!----------------------------------------------------------------------

    call orbit

!----------------------------------------------------------------------

end subroutine SET_ORBITAL_PARAMETERS_

!> @brief get_orbital_parameters retrieves the orbital parameters for use
!!        by another module.
!!
!! @throw FATAL, "astronomy_mod module has not been initialized"

subroutine GET_ORBITAL_PARAMETERS_(ecc_out, obliq_out, per_out)

!-------------------------------------------------------------------
!    get_orbital_parameters retrieves the orbital parameters for use
!    by another module.
!--------------------------------------------------------------------


real(kind=FMS_AST_KIND_), intent(out) :: ecc_out !< Eccentricity of orbital ellipse [dimensionless]
real(kind=FMS_AST_KIND_), intent(out) :: obliq_out !< Obliquity [degrees]
real(kind=FMS_AST_KIND_), intent(out) :: per_out !< Longitude of perihelion with respect to autumnal
                             !! equinox in northern hemisphere [degrees]

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------

    if (.not. module_is_initialized) &
        call error_mesg('astronomy_mod', 'module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    fill the output arguments with the eccentricity, obliquity and
!    perihelion angle.
!--------------------------------------------------------------------

    ecc_out   = real(ecc,   FMS_AST_KIND_)
    obliq_out = real(obliq, FMS_AST_KIND_)
    per_out   = real(per,   FMS_AST_KIND_)

end subroutine GET_ORBITAL_PARAMETERS_


!> @brief diurnal_solar_2d returns 2d fields of cosine of zenith angle,
!!        daylight fraction and earth-sun distance at the specified latitudes,
!!        longitudes and time. These values may be instantaneous or averaged
!!        over a specified time interval.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <lon> Longitudes of model grid points
!! @param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! @param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! @param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! @param [in] <allow_negative_cosz> Allow negative values for cosz?
!! @param [out] <half_day_out> half_day_out
!!
!! @throw FATAL, "astronomy_mod time_since_ae not between 0 and 2pi"
!! @throw FATAL, "astronomy_mod gmt not between 0 and 2pi"

subroutine DIURNAL_SOLAR_2D_(lat, lon, gmt, time_since_ae, cosz, &
    fracday, rrsun, dt, allow_negative_cosz, &
    half_day_out)

real(kind=FMS_AST_KIND_), dimension(:,:), intent(in)            :: lat, lon
real(kind=FMS_AST_KIND_),                 intent(in)            :: gmt, time_since_ae
real(kind=FMS_AST_KIND_), dimension(:,:), intent(out)           :: cosz, fracday
real(kind=FMS_AST_KIND_),                 intent(out)           :: rrsun
real(kind=FMS_AST_KIND_),                 intent(in),  optional :: dt
logical,                                  intent(in),  optional :: allow_negative_cosz
real(kind=FMS_AST_KIND_), dimension(:,:), intent(out), optional :: half_day_out

!---------------------------------------------------------------------
!   local variables


real(kind=FMS_AST_KIND_), dimension(size(lat,1),size(lat,2)) :: t, tt, h, aa, bb, st, stt, sh
real(kind=FMS_AST_KIND_)                                     :: ang, dec
logical                                                      :: Lallow_negative
integer, parameter                                           :: lkind = FMS_AST_KIND_


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

    if (time_since_ae < 0.0_lkind .or. time_since_ae > real(twopi, FMS_AST_KIND_)) &
        call error_mesg('astronomy_mod', 'time_since_ae not between 0 and 2pi', FATAL)

!--------------------------------------------------------------------
!    be sure the time at longitude = 0.0 is legitimate.
!---------------------------------------------------------------------

    if (gmt < 0.0_lkind .or. gmt > real(twopi, FMS_AST_KIND_)) &
        call error_mesg('astronomy_mod', 'gmt not between 0 and 2pi', FATAL)

!---------------------------------------------------------------------
!>    define the orbital angle (location in year), solar declination and
!!    earth sun distance factor. use functions contained in this module.
!---------------------------------------------------------------------

    ang   = angle(time_since_ae)
    dec   = declination(ang)
    rrsun = r_inv_squared(ang)

!---------------------------------------------------------------------
!>    define terms needed in the cosine zenith angle equation.
!--------------------------------------------------------------------

    aa = sin(lat)*sin(dec)
    bb = cos(lat)*cos(dec)

!---------------------------------------------------------------------
!>    define local time. force it to be between -pi and pi.
!--------------------------------------------------------------------

    t = gmt + lon - real(PI, FMS_AST_KIND_)
    where(t >= real(PI, FMS_AST_KIND_)) t = real(t, FMS_AST_KIND_) - real(twopi, FMS_AST_KIND_)
    where(t < real(-PI,FMS_AST_KIND_))  t = real(t, FMS_AST_KIND_) + real(twopi, FMS_AST_KIND_)

    Lallow_negative = .false.
    if (present(allow_negative_cosz)) then
        if (allow_negative_cosz) Lallow_negative = .true.
    end if


!---------------------------------------------------------------------
!>    perform a time integration to obtain cosz and fracday if desired.
!!    output is valid over the period from t to t + dt.
!--------------------------------------------------------------------

    h   = half_day(lat,dec)

    if ( present(half_day_out) ) then
        half_day_out = h
    end if

    if ( present(dt) ) then   ! (perform time averaging)
        tt   = t + dt
        st   = sin(t)
        stt  = sin(tt)
        sh   = sin(h)
        cosz = 0.0_lkind

        if (.not. Lallow_negative) then
!-------------------------------------------------------------------
!    case 1: entire averaging period is before sunrise.
!-------------------------------------------------------------------

        where (t < -h .and. tt < -h) cosz = 0.0_lkind

!-------------------------------------------------------------------
!    case 2: averaging period begins before sunrise, ends after sunrise
!    but before sunset
!-------------------------------------------------------------------

        where ( (tt+h) /= 0.0_lkind .and. t < -h .and. abs(tt) <= h)   &
            cosz = aa + bb*(stt + sh)/ (tt + h)

!-------------------------------------------------------------------
!    case 3: averaging period begins before sunrise, ends after sunset,
!    but before the next sunrise. modify if averaging period extends
!    past the next day's sunrise, but if averaging period is less than
!    a half- day (pi) that circumstance will never occur.
!-------------------------------------------------------------------

        where (t < -h .and. h /= 0.0_lkind .and. h < tt)    &
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

        where ((h-t) /= 0.0_lkind .and. abs(t) <= h .and. h < tt)    &
            cosz = aa + bb*(sh - st)/(h-t)

!-------------------------------------------------------------------
!    case 6: averaging period begins after sunrise , ends after the
!    next day's sunrise. note that this includes the case when the
!    day length is one day (h = pi).
!-------------------------------------------------------------------

        where (real(twopi, FMS_AST_KIND_) - h < tt .and. &
            (tt+h-real(twopi, FMS_AST_KIND_)) /= 0.0_lkind .and. t <= h ) &
            cosz = (real(cosz*(h - t), FMS_AST_KIND_) + (real(aa, FMS_AST_KIND_)*(real(tt, FMS_AST_KIND_) &
                   + real(h, FMS_AST_KIND_) - real(twopi, FMS_AST_KIND_)) + real(bb, FMS_AST_KIND_) &
                   * (real(stt, FMS_AST_KIND_) + real(sh, FMS_AST_KIND_)))) / ((real(h, FMS_AST_KIND_) &
                   - real(t, FMS_AST_KIND_)) + (real(tt, FMS_AST_KIND_) + real(h, FMS_AST_KIND_) &
                   - real(twopi, FMS_AST_KIND_)))

!-------------------------------------------------------------------
!    case 7: averaging period begins after sunset and ends before the
!    next day's sunrise
!-------------------------------------------------------------------

        where(h <  t .and. real(twopi, FMS_AST_KIND_) - h >= tt) &
            cosz = 0.0_lkind

!-------------------------------------------------------------------
!    case 8: averaging period begins after sunset and ends after the
!    next day's sunrise but before the next day's sunset. if the
!    averaging period is less than a half-day (pi) the latter
!    circumstance will never occur.
!-----------------------------------------------------------------

        where (h <  t .and. real(twopi, FMS_AST_KIND_) - h < tt)
            cosz = real(aa, FMS_AST_KIND_) + real(bb, FMS_AST_KIND_)&
                   *(real(stt, FMS_AST_KIND_) + real(sh,FMS_AST_KIND_)) &
                   / (real(tt, FMS_AST_KIND_) + real(h, FMS_AST_KIND_) - real(twopi, FMS_AST_KIND_))
        end where

        else
            cosz = aa + bb*(stt - st)/ (tt - t)
        end if



!-------------------------------------------------------------------
!    day fraction is the fraction of the averaging period contained
!    within the (-h,h) period.
!-------------------------------------------------------------------

    where (t < -h .and. tt < -h)           fracday = 0.0_lkind
    where (t < -h .and. abs(tt) <= h)      fracday = (tt + h )/dt
    where (t < -h .and. h < tt)            fracday = ( h + h )/dt
    where (abs(t) <= h .and. abs(tt) <= h) fracday = (tt - t )/dt
    where (abs(t) <= h .and. h < tt)       fracday = ( h - t )/dt
    where (h <  t)                         fracday = 0.0_lkind
    where (real(twopi, FMS_AST_KIND_) - h < tt)  &
                                           fracday = real(fracday, FMS_AST_KIND_) + &
                                                    (real(tt, FMS_AST_KIND_) + real(h, FMS_AST_KIND_) &
                                                    - real(twopi, FMS_AST_KIND_))/real(dt, FMS_AST_KIND_)
!----------------------------------------------------------------------
!>    if instantaneous values are desired, define cosz at time t.
!----------------------------------------------------------------------
    else  ! (no time averaging)
        if (.not. Lallow_negative) then
            where (abs(t) < h)
              cosz    = aa + bb*cos(t)
              fracday = 1.0_lkind
            elsewhere
              cosz    = 0.0_lkind
              fracday = 0.0_lkind
        end where
        else
            cosz = aa + bb*cos(t)
            where (abs(t) < h)
              fracday = 1.0_lkind
            elsewhere
              fracday = 0.0_lkind
            end where
        end if
    end if


!----------------------------------------------------------------------
!>    Check that cosz is not negative, if desired.
!----------------------------------------------------------------------

    if (.not. Lallow_negative) then
        cosz = max(0.0_lkind, cosz)
    end if

end subroutine DIURNAL_SOLAR_2D_


!> @brief diurnal_solar_1d takes 1-d input fields, adds a second dimension
!!        and calls diurnal_solar_2d. on return, the 2d fields are returned
!!        to the original 1d fields.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <lon> Longitudes of model grid points
!! @param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! @param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! @param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! @param [in] <allow_negative_cosz> Allow negative values for cosz?
!! @param [out] <half_day_out> half_day_out

subroutine DIURNAL_SOLAR_1D_(lat, lon, gmt, time_since_ae, cosz, &
                             fracday, rrsun, dt, allow_negative_cosz, &
                             half_day_out)

!---------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(:),  intent(in)            :: lat, lon
real(kind=FMS_AST_KIND_),                intent(in)            :: gmt, time_since_ae
real(kind=FMS_AST_KIND_), dimension(:),  intent(out)           :: cosz, fracday
real(kind=FMS_AST_KIND_),                intent(out)           :: rrsun
real(kind=FMS_AST_KIND_),                intent(in),  optional :: dt
logical,                                 intent(in),  optional :: allow_negative_cosz
real(kind=FMS_AST_KIND_), dimension(:),  intent(out), optional :: half_day_out

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(size(lat),1) :: lat_2d, lon_2d, cosz_2d, fracday_2d, halfday_2d

!--------------------------------------------------------------------
!>    define 2-d versions of input data arrays.
!--------------------------------------------------------------------

    lat_2d(:,1) = lat
    lon_2d(:,1) = lon


!--------------------------------------------------------------------
!>    call diurnal_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
!     if (present(dt)) then

    call diurnal_solar(lat_2d, lon_2d, gmt, time_since_ae, &
         cosz_2d, fracday_2d, rrsun, dt=dt, &
         allow_negative_cosz=allow_negative_cosz, half_day_out=halfday_2d)
!     else
!       call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae, &
!                              cosz_2d, fracday_2d, rrsun)
!     endif

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------

    fracday = fracday_2d(:,1)
    cosz    = cosz_2d (:,1)
    if (present(half_day_out)) then
        half_day_out = halfday_2d(:,1)
    end if

end subroutine DIURNAL_SOLAR_1D_



!> @brief diurnal_solar_0d takes scalar input fields, makes them into 2d
!!        arrays dimensioned (1,1), and calls diurnal_solar_2d. on return,
!!        the 2d fields are converted back to the desired scalar output.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <lon> Longitudes of model grid points
!! @param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! @param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! @param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! @param [in] <allow_negative_cosz> Allow negative values for cosz?
!! @param [out] <half_day_out> half_day_out

subroutine DIURNAL_SOLAR_0D_(lat, lon, gmt, time_since_ae, cosz,  &
                             fracday, rrsun, dt, allow_negative_cosz, &
                             half_day_out)

real(kind=FMS_AST_KIND_), intent(in)            :: lat, lon, gmt, time_since_ae
real(kind=FMS_AST_KIND_), intent(out)           :: cosz, fracday, rrsun
real(kind=FMS_AST_KIND_), intent(in),  optional :: dt
logical,                  intent(in),  optional :: allow_negative_cosz
real(kind=FMS_AST_KIND_), intent(out), optional :: half_day_out

!--------------------------------------------------------------------
!  local variables:
!--------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(1,1) :: lat_2d, lon_2d, cosz_2d, fracday_2d, halfday_2d

!---------------------------------------------------------------------
!>    create 2d arrays from the scalar input fields.
!---------------------------------------------------------------------

    lat_2d = lat
    lon_2d = lon

!--------------------------------------------------------------------
!>    call diurnal_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
!     if (present(dt)) then

    call diurnal_solar(lat_2d, lon_2d, gmt, time_since_ae,  &
         cosz_2d, fracday_2d, rrsun, dt=dt, &
         allow_negative_cosz=allow_negative_cosz, half_day_out=halfday_2d)
!     else
!       call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae, &
!                              cosz_2d, fracday_2d, rrsun)
!     end if

!-------------------------------------------------------------------
!>    place output fields into scalars for return to calling routine.
!-------------------------------------------------------------------

    fracday = fracday_2d(1,1)
    cosz    = cosz_2d(1,1)
    if (present(half_day_out)) then
        half_day_out = halfday_2d(1,1)
    end if

end subroutine DIURNAL_SOLAR_0D_


!> @brief diurnal_solar_cal_2d receives time_type inputs, converts
!!        them to real variables and then calls diurnal_solar_2d to
!!        compute desired astronomical variables.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <lon> Longitudes of model grid points
!! @param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! @param [in] <time> Time of year (time_type)
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! @param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! @param [in] <allow_negative_cosz> Allow negative values for cosz?
!! @param [out] <half_day_out> half_day_out
!!
!! @throw FATAL, "astronomy_mod radiation time step must be no longer than 12 hrs"
!! @throw FATAL, "astronomy_mod radiation time step must not be an integral number of days"

subroutine DIURNAL_SOLAR_CAL_2D_(lat, lon, time, cosz, fracday,   &
                                 rrsun, dt_time, allow_negative_cosz, &
                                 half_day_out)

!-------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(:,:), intent(in)             :: lat, lon
type(time_type),                          intent(in)             :: time
real(kind=FMS_AST_KIND_), dimension(:,:), intent(out)            :: cosz, fracday
real(kind=FMS_AST_KIND_),                 intent(out)            :: rrsun
type(time_type),                          intent(in),  optional  :: dt_time
logical,                                  intent(in),  optional  :: allow_negative_cosz
real(kind=FMS_AST_KIND_), dimension(:,:), intent(out), optional  :: half_day_out

!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------

real(kind=FMS_AST_KIND_) :: dt
real(kind=FMS_AST_KIND_) :: gmt, time_since_ae
integer, parameter       :: lkind = FMS_AST_KIND_

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
        if (dt > real(PI, FMS_AST_KIND_)) then
           call error_mesg('astronomy_mod', 'radiation time step must be no longer than 12 hrs', FATAL)
        endif
    if (dt == 0.0_lkind) then
        call error_mesg('astronomy_mod', 'radiation time step must not be an integral number of days', FATAL)
    endif

!--------------------------------------------------------------------
!>    Call diurnal_solar_2d to calculate astronomy fields, with or
!!    without the optional argument dt.
!--------------------------------------------------------------------

    call diurnal_solar(lat, lon, gmt, time_since_ae, cosz, &
        fracday, rrsun, dt=dt, &
        allow_negative_cosz=allow_negative_cosz, &
        half_day_out=half_day_out)
    else
        call diurnal_solar(lat, lon, gmt, time_since_ae, cosz, &
            fracday, rrsun, allow_negative_cosz=allow_negative_cosz, &
            half_day_out=half_day_out)
    end if

end subroutine DIURNAL_SOLAR_CAL_2D_


!> @brief diurnal_solar_cal_1d receives time_type inputs, converts
!!        them to real variables and then calls diurnal_solar_2d to
!!        compute desired astronomical variables.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <lon> Longitudes of model grid points
!! @param [in] <gmt> Time of day at longitude 0.0; midnight = 0.0, one day = 2 * pi
!! @param [in] <time> Time of year (time_type)
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> earth-sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!! @param [in] <dt> OPTIONAL: time interval after gmt over which the astronomical variables are to be
!!                  averaged. this produces averaged output rather than instantaneous.
!! @param [in] <allow_negative_cosz> Allow negative values for cosz?
!! @param [out] <half_day_out> half_day_out

subroutine DIURNAL_SOLAR_CAL_1D_(lat, lon, time, cosz, fracday,   &
                                 rrsun, dt_time, allow_negative_cosz, &
                                 half_day_out)

!--------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(:), intent(in)            :: lat, lon
type(time_type),                        intent(in)            :: time
real(kind=FMS_AST_KIND_), dimension(:), intent(out)           :: cosz, fracday
real(kind=FMS_AST_KIND_),               intent(out)           :: rrsun
type(time_type),                        intent(in),  optional :: dt_time
logical,                                intent(in),  optional :: allow_negative_cosz
real(kind=FMS_AST_KIND_), dimension(:), intent(out), optional :: half_day_out
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables
!-------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(size(lat),1) :: lat_2d, lon_2d, cosz_2d, &
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
        call diurnal_solar(lat_2d, lon_2d, time, cosz_2d,    &
             fracday_2d, rrsun, dt_time=dt_time, &
             allow_negative_cosz=allow_negative_cosz, &
             half_day_out=halfday_2d)
    else
        call diurnal_solar(lat_2d, lon_2d, time, cosz_2d,    &
             fracday_2d, rrsun, allow_negative_cosz=allow_negative_cosz, &
             half_day_out=halfday_2d)
    end if

!-------------------------------------------------------------------
!>    Place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------

    fracday = fracday_2d(:,1)
    cosz    = cosz_2d (:,1)
    if (present(half_day_out)) then
        half_day_out = halfday_2d(:,1)
    end if

end subroutine DIURNAL_SOLAR_CAL_1D_


!> @brief diurnal_solar_cal_0d receives time_type inputs, converts them to real variables
!!        and then calls diurnal_solar_2d to compute desired astronomical variables.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <lon> Longitudes of model grid points
!! @param [in] <time> Time of year (time_type)
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a) : (a/r)**2
!! @param [out] <dt_time> OPTIONAL: time interval after gmt over which the astronomical variables are
!!                        to be averaged. this produces averaged output rather than instantaneous.
!! @param [in] <allow_negative_cosz> allow_negative_cosz
!! @param [out] <half_day_out> half_day_out

subroutine DIURNAL_SOLAR_CAL_0D_(lat, lon, time, cosz, fracday,   &
                                 rrsun, dt_time, allow_negative_cosz, &
                                 half_day_out)

!---------------------------------------------------------------------

real(kind=FMS_AST_KIND_),  intent(in)            :: lat, lon
type(time_type),           intent(in)            :: time
real(kind=FMS_AST_KIND_),  intent(out)           :: cosz, fracday, rrsun
type(time_type),           intent(in),  optional :: dt_time
logical,                   intent(in),  optional :: allow_negative_cosz
real(kind=FMS_AST_KIND_),  intent(out), optional :: half_day_out

!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(1,1) :: lat_2d, lon_2d, cosz_2d, fracday_2d, halfday_2d

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
        call diurnal_solar(lat_2d, lon_2d, time, cosz_2d,   &
             fracday_2d, rrsun, dt_time=dt_time, &
             allow_negative_cosz=allow_negative_cosz, &
             half_day_out=halfday_2d)
    else
        call diurnal_solar(lat_2d, lon_2d, time, cosz_2d,   &
             fracday_2d, rrsun, allow_negative_cosz=allow_negative_cosz, &
             half_day_out=halfday_2d)
    end if

!-------------------------------------------------------------------
!>    Place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------

    fracday = fracday_2d(1,1)
    cosz    = cosz_2d(1,1)
    if (present(half_day_out)) then
          half_day_out = halfday_2d(1,1)
      end if

end subroutine DIURNAL_SOLAR_CAL_0D_


!> @brief daily_mean_solar_2d computes the daily mean astronomical parameters for
!!        the input points at latitude lat and time of year time_since_ae.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <h_out> 2-d array of half-day lengths at the latitudes
!! @param [out] <rr_out> the inverse of the square of the earth-sun distance relative
!!                       to the mean distance at angle ang in the earth's orbit.
!!
!! @throw FATAL, "astronomy_mod time_since_ae not between 0 and 2pi"

subroutine DAILY_MEAN_SOLAR_2D_(lat, time_since_ae, cosz, h_out, rr_out)

!----------------------------------------------------------------------
real(kind=FMS_AST_KIND_), dimension(:,:), intent(in)   :: lat
real(kind=FMS_AST_KIND_),                 intent(in)   :: time_since_ae
real(kind=FMS_AST_KIND_), dimension(:,:), intent(out)  :: cosz, h_out
real(kind=FMS_AST_KIND_),                 intent(out)  :: rr_out
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(size(lat,1),size(lat,2)) :: h
real(kind=FMS_AST_KIND_) :: ang, dec, rr
integer, parameter       :: lkind = FMS_AST_KIND_

!--------------------------------------------------------------------
!    be sure the time in the annual cycle is legitimate.
!---------------------------------------------------------------------

    if (time_since_ae < 0.0_lkind .or. time_since_ae > real(twopi, FMS_AST_KIND_)) &
        call error_mesg('astronomy_mod', 'time_since_ae not between 0 and 2pi', FATAL)

!---------------------------------------------------------------------
!>    Define the orbital angle (location in year), solar declination,
!!    half-day length and earth sun distance factor. Use functions
!!    contained in this module.
!---------------------------------------------------------------------

    ang = angle (time_since_ae)
    dec = declination (ang)
    h   = half_day (lat, dec)
    rr  = r_inv_squared (ang)

!---------------------------------------------------------------------
!>    Where the entire day is dark, define cosz to be zero. otherwise
!!    use the standard formula. Define the daylight fraction and earth-
!!    sun distance.
!---------------------------------------------------------------------

    where (h == 0.0_lkind)
          cosz = 0.0_lkind
    elsewhere
          cosz = sin(lat)*sin(dec) + cos(lat)*cos(dec)*sin(h)/h
    end where
      h_out    = h/real(PI, FMS_AST_KIND_)
      rr_out   = rr

end subroutine DAILY_MEAN_SOLAR_2D_


!> @brief daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!!        and calls daily_mean_solar_2d. on return, the 2d fields are
!!        returned to the original 1d fields.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <h_out> 2-d array of half-day lengths at the latitudes
!! @param [out] <rr_out> the inverse of the square of the earth-sun distance relative
!!                       to the mean distance at angle ang in the earth's orbit.

subroutine DAILY_MEAN_SOLAR_1D_(lat, time_since_ae, cosz, h_out, rr_out)

!----------------------------------------------------------------------
real(kind=FMS_AST_KIND_), intent(in), dimension(:)             :: lat
real(kind=FMS_AST_KIND_), intent(in)                           :: time_since_ae
real(kind=FMS_AST_KIND_), intent(out), dimension(size(lat(:))) :: cosz
real(kind=FMS_AST_KIND_), intent(out), dimension(size(lat(:))) :: h_out
real(kind=FMS_AST_KIND_), intent(out)                          :: rr_out

!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables
!----------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(size(lat),1) :: lat_2d, cosz_2d, hout_2d

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------

    lat_2d(:,1) = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------

    call daily_mean_solar(lat_2d, time_since_ae, cosz_2d,      &
                                hout_2d, rr_out)

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------

    h_out = hout_2d(:,1)
    cosz  = cosz_2d(:,1)

end subroutine DAILY_MEAN_SOLAR_1D_


!> @brief daily_mean_solar_2level takes 1-d input fields, adds a second
!!        dimension and calls daily_mean_solar_2d. on return, the 2d fields
!!        are returned to the original 1d fields.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <solar> Shortwave flux factor: cosine of zenith angle * daylight fraction /
!! (earth-sun distance squared)

subroutine DAILY_MEAN_SOLAR_2LEVEL_(lat, time_since_ae, cosz, solar)

!----------------------------------------------------------------------
real(kind=FMS_AST_KIND_), intent(in), dimension(:)             :: lat
real(kind=FMS_AST_KIND_), intent(in)                           :: time_since_ae
real(kind=FMS_AST_KIND_), intent(out), dimension(size(lat(:))) :: cosz, solar

!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables
!----------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(size(lat),1) :: lat_2d, cosz_2d, hout_2d
real(kind=FMS_AST_KIND_)                         :: rr_out

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------

    lat_2d(:,1) = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------

    call daily_mean_solar(lat_2d, time_since_ae, cosz_2d,      &
                                hout_2d, rr_out)

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------

    solar = cosz_2d(:,1)*hout_2d(:,1)*rr_out
    cosz  = cosz_2d(:,1)

end subroutine DAILY_MEAN_SOLAR_2LEVEL_


!> @brief daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!!        and calls daily_mean_solar_2d. on return, the 2d fields are
!!        returned to the original 1d fields.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <time_since_ae> Time of year; autumnal equinox = 0.0, one year = 2 * pi
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <h_out> 2-d array of half-day lengths at the latitudes
!! @param [out] <rr_out> the inverse of the square of the earth-sun distance relative to
!!                       the mean distance at angle ang in the earth's orbit.

subroutine DAILY_MEAN_SOLAR_0D_(lat, time_since_ae, cosz, h_out, rr_out)

real(kind=FMS_AST_KIND_), intent(in)     :: lat, time_since_ae
real(kind=FMS_AST_KIND_), intent(out)    :: cosz, h_out, rr_out

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(1,1) :: lat_2d, cosz_2d, hout_2d

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------

    lat_2d = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------

    call daily_mean_solar(lat_2d, time_since_ae, cosz_2d,     &
                                hout_2d, rr_out)

!-------------------------------------------------------------------
!>    return output fields to scalars for return to calling routine.
!-------------------------------------------------------------------

    h_out = hout_2d(1,1)
    cosz  = cosz_2d(1,1)

end subroutine DAILY_MEAN_SOLAR_0D_


!> @brief daily_mean_solar_cal_2d receives time_type inputs, converts
!!        them to real variables and then calls daily_mean_solar_2d to
!!        compute desired astronomical variables.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <time> Time of year (time_type)
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2
!!
!! @throw FATAL, "astronomy_mod time_since_ae not between 0 and 2pi"

subroutine DAILY_MEAN_SOLAR_CAL_2D_(lat, time, cosz, fracday, rrsun)

!-------------------------------------------------------------------
real(kind=FMS_AST_KIND_), dimension(:,:), intent(in)  :: lat
type(time_type),                          intent(in)  :: time
real(kind=FMS_AST_KIND_), dimension(:,:), intent(out) :: cosz, fracday
real(kind=FMS_AST_KIND_),                 intent(out) :: rrsun
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables
!-------------------------------------------------------------------

real(kind=FMS_AST_KIND_) :: time_since_ae
integer, parameter       :: lkind = FMS_AST_KIND_

!--------------------------------------------------------------------
!    be sure the time in the annual cycle is legitimate.
!---------------------------------------------------------------------

    time_since_ae = orbital_time(time)
    if (time_since_ae < 0.0_lkind .or. time_since_ae > real(twopi, FMS_AST_KIND_)) &
        call error_mesg('astronomy_mod', 'time_since_ae not between 0 and 2pi', FATAL)

!--------------------------------------------------------------------
!    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------

    call daily_mean_solar(lat, time_since_ae, cosz, fracday, rrsun)

end subroutine DAILY_MEAN_SOLAR_CAL_2D_


!> @brief daily_mean_solar_cal_1d receives time_type inputs, converts
!!        them to real, 2d variables and then calls daily_mean_solar_2d to
!!        compute desired astronomical variables.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <time> Time of year (time_type)
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2

subroutine DAILY_MEAN_SOLAR_CAL_1D_(lat, time, cosz, fracday, rrsun)

real(kind=FMS_AST_KIND_), dimension(:),  intent(in)   :: lat
type(time_type),                         intent(in)   :: time
real(kind=FMS_AST_KIND_), dimension(:),  intent(out)  :: cosz, fracday
real(kind=FMS_AST_KIND_),                intent(out)  :: rrsun

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(size(lat),1) :: lat_2d, cosz_2d, fracday_2d

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------

    lat_2d(:,1) = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_cal_2d to convert the time_types to reals and
!!    then calculate the astronomy fields.
!--------------------------------------------------------------------

    call daily_mean_solar(lat_2d, time, cosz_2d, fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------

    fracday = fracday_2d(:,1)
    cosz    = cosz_2d(:,1)

end subroutine DAILY_MEAN_SOLAR_CAL_1D_


!> @brief daily_mean_solar_cal_2level receives 1d arrays and time_type input,
!!        converts them to real, 2d variables and then calls
!!        daily_mean_solar_2d to compute desired astronomical variables.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <time> Time of year (time_type)
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <solar> Shortwave flux factor: cosine of zenith angle * daylight fraction /
!! (earth-sun distance squared)

subroutine DAILY_MEAN_SOLAR_CAL_2LEVEL_(lat, time, cosz, solar)

real(kind=FMS_AST_KIND_), dimension(:),  intent(in)   :: lat
type(time_type),                         intent(in)   :: time
real(kind=FMS_AST_KIND_), dimension(:),  intent(out)  :: cosz, solar
!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(size(lat),1) :: lat_2d, cosz_2d, fracday_2d
real(kind=FMS_AST_KIND_)                         :: rrsun


!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------

    lat_2d(:,1) = real(lat, FMS_AST_KIND_)

!--------------------------------------------------------------------
!>    call daily_mean_solar_cal_2d to convert the time_types to reals and
!!    then calculate the astronomy fields.
!--------------------------------------------------------------------

    call daily_mean_solar(lat_2d, time, cosz_2d, fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into 1-d arguments for return to
!!    calling routine.
!-------------------------------------------------------------------

    solar = cosz_2d(:,1)*fracday_2d(:,1)*rrsun
    cosz  = cosz_2d(:,1)

end subroutine DAILY_MEAN_SOLAR_CAL_2LEVEL_


!> @brief daily_mean_solar_cal_0d converts scalar input fields to real, 2d variables and
!!        then calls daily_mean_solar_2d to compute desired astronomical variables.
!!
!! @param [in] <lat> Latitudes of model grid points
!! @param [in] <time> Time of year (time_type)
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2

subroutine DAILY_MEAN_SOLAR_CAL_0D_(lat, time, cosz, fracday, rrsun)

!--------------------------------------------------------------------
real(kind=FMS_AST_KIND_), intent(in)  :: lat
type(time_type),          intent(in)  :: time
real(kind=FMS_AST_KIND_), intent(out) :: cosz, fracday, rrsun
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables
!--------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(1,1) :: lat_2d, cosz_2d, fracday_2d

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------

      lat_2d = lat

!--------------------------------------------------------------------
!>    call daily_mean_solar_cal_2d to convert the time_types to reals and
!!    then calculate the astronomy fields.
!--------------------------------------------------------------------

    call daily_mean_solar(lat_2d, time, cosz_2d, fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into scalar arguments for return to
!!    calling routine.
!-------------------------------------------------------------------

    fracday = fracday_2d(1,1)
      cosz  = cosz_2d(1,1)

end subroutine DAILY_MEAN_SOLAR_CAL_0D_


!> @brief annual_mean_solar_2d returns 2d fields of annual mean values of the cosine of
!!        zenith angle, daylight fraction and earth-sun distance at the specified latitude.
!!
!! @param [in] <jst> Starting index of latitude window
!! @param [in] <jnd> Ending index of latitude window
!! @param [in] <lat> Latitudes of model grid points
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <solar> Shortwave flux factor: cosine of zenith angle * daylight fraction /
!! (earth-sun distance squared)
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2

subroutine ANNUAL_MEAN_SOLAR_2D_(js, je, lat, cosz, solar, fracday,  &
                                 rrsun)

!--------------------------------------------------------------------
integer,                                     intent(in)    :: js, je
real(kind=FMS_AST_KIND_), dimension(:,:),    intent(in)    :: lat
real(kind=FMS_AST_KIND_), dimension(:,:),    intent(out)   :: solar, cosz, fracday
real(kind=FMS_AST_KIND_),                    intent(out)   :: rrsun
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables
!--------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(size(lat,1),size(lat,2)) :: s,z
real(kind=FMS_AST_KIND_)                                     :: t
integer                                                      :: n
integer, parameter                                           :: lkind = FMS_AST_KIND_

!--------------------------------------------------------------------
!    if the calculation has not yet been done, do it here.
!--------------------------------------------------------------------

    if (.not. annual_mean_calculated) then

!----------------------------------------------------------------------
!>    determine annual mean values of solar flux and product of cosz
!!    and solar flux by integrating the annual cycle in num_angles
!!    orbital increments.
!----------------------------------------------------------------------

        solar   = 0.0_lkind
        cosz    = 0.0_lkind
        do n =1, num_angles
          t     = real(float(n-1), FMS_AST_KIND_)*real(twopi, FMS_AST_KIND_)/real(float(num_angles), FMS_AST_KIND_)
          call daily_mean_solar(lat,t, z, fracday, rrsun)
          s     = z*rrsun*fracday
          solar = solar + s
          cosz  = cosz  + z*s
        end do
        solar   = solar/float(num_angles)
        cosz    = cosz/float(num_angles)

!--------------------------------------------------------------------
!>   define the flux-weighted annual mean cosine of the zenith angle.
!--------------------------------------------------------------------

        where(solar .eq. 0.0_lkind)
          cosz = 0.0_lkind
        elsewhere
          cosz = cosz/solar
        end where

!-------------------------------------------------------------------
!>    define avg fracday such as to make the avg flux (solar) equal to
!!    the product of the avg cosz * avg fracday * assumed mean avg
!!    radius of 1.0. it is unlikely that these avg fracday and avg rr
!!    will ever be used.
!--------------------------------------------------------------------

        where(solar .eq. 0.0_lkind)
          fracday = 0.0_lkind
        elsewhere
          fracday = solar/cosz
        end where
        rrsun = 1.00_lkind

!---------------------------------------------------------------------
!>    save the values that have been calculated as module variables, if
!!    those variables are present; i.e., not the spectral 2-layer model.
!---------------------------------------------------------------------
        if (allocated (cosz_ann)) then

          cosz_ann    = real(cosz,    FMS_AST_KIND_)
          solar_ann   = real(solar,   FMS_AST_KIND_)
          fracday_ann = real(fracday, FMS_AST_KIND_)
          rrsun_ann   = real(rrsun,   FMS_AST_KIND_)

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

          cosz    = real(cosz_ann,    FMS_AST_KIND_)
          solar   = real(solar_ann,   FMS_AST_KIND_)
          fracday = real(fracday_ann, FMS_AST_KIND_)
          rrsun   = real(rrsun_ann,   FMS_AST_KIND_)
        endif
    endif

end subroutine ANNUAL_MEAN_SOLAR_2D_


!> @brief annual_mean_solar_1d creates 2-d input fields from 1-d input fields and then calls
!!        annual_mean_solar_2d to obtain 2-d output fields which are then stored into 1-d
!!        fields for return to the calling subroutine.
!!
!! @param [in] <jst> Starting index of latitude window
!! @param [in] <jnd> Ending index of latitude window
!! @param [in] <lat> Latitudes of model grid points
!! @param [out] <cosz> Cosine of solar zenith angle
!! @param [out] <solar> Shortwave flux factor: cosine of zenith angle * daylight fraction /
!! (earth-sun distance squared)
!! @param [out] <fracday> Daylight fraction of time interval
!! @param [out] <rrsun_out> Earth-Sun distance (r) relative to semi-major axis of orbital ellipse (a):(a/r)**2

subroutine ANNUAL_MEAN_SOLAR_1D_(jst, jnd, lat, cosz, solar,  &
                                 fracday, rrsun_out)

!---------------------------------------------------------------------
integer,                                intent(in)     :: jst, jnd
real(kind=FMS_AST_KIND_), dimension(:), intent(in)     :: lat(:)
real(kind=FMS_AST_KIND_), dimension(:), intent(out)    :: cosz, solar, fracday
real(kind=FMS_AST_KIND_),               intent(out)    :: rrsun_out
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables


real(kind=FMS_AST_KIND_), dimension(size(lat),1) :: lat_2d, solar_2d, cosz_2d, fracday_2d
real(kind=FMS_AST_KIND_)                         :: rrsun

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

        call annual_mean_solar(jst, jnd, lat_2d, cosz_2d, solar_2d, fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into 1-D arrays for return to calling routine.
!-------------------------------------------------------------------

        fracday   = fracday_2d(:,1)
        rrsun_out = rrsun
        solar     = solar_2d(:,1)
        cosz      = cosz_2d(:,1)

!--------------------------------------------------------------------
!>    if the calculation has been done, simply return the module
!!    variables contain the results at the desired latitudes.
!--------------------------------------------------------------------

    else
        cosz(:)    = real(cosz_ann(1,jst:jnd),    FMS_AST_KIND_)
        solar(:)   = real(solar_ann(1,jst:jnd),   FMS_AST_KIND_)
        fracday(:) = real(fracday_ann(1,jst:jnd), FMS_AST_KIND_)
        rrsun      = real(rrsun_ann,              FMS_AST_KIND_)
    endif

end subroutine ANNUAL_MEAN_SOLAR_1D_


!> @brief annual_mean_solar_2level creates 2-d input fields from 1-d input fields
!!        and then calls annual_mean_solar_2d to obtain 2-d output fields which are
!!        then stored into 1-d fields for return to the calling subroutine. This
!!        subroutine will be called during model initialization.
!!
!! @throw FATAL, "astronomy_mod annual_mean_solar_2level should be called only once"

subroutine ANNUAL_MEAN_SOLAR_2LEVEL_(lat, cosz, solar)

!---------------------------------------------------------------------
real(kind=FMS_AST_KIND_), dimension(:), intent(in)     :: lat !< Latitudes of model grid points
real(kind=FMS_AST_KIND_), dimension(:), intent(out)    :: cosz !< Cosine of solar zenith angle
real(kind=FMS_AST_KIND_), dimension(:), intent(out)    :: solar !< shortwave flux factor: cosine of zenith angle *
                                            !! daylight fraction / (earth-sun distance squared)
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables


real(kind=FMS_AST_KIND_), dimension(size(lat),1) :: lat_2d, solar_2d, cosz_2d, fracday_2d
integer                                          :: jst, jnd
real(kind=FMS_AST_KIND_)                         :: rrsun

!--------------------------------------------------------------------
!    if the calculation has not been done, do it here.
!--------------------------------------------------------------------

    if ( .not. annual_mean_calculated) then

!--------------------------------------------------------------------
!>    define 2-d versions of input data array.
!--------------------------------------------------------------------

        lat_2d(:,1) = lat
        jst         = 1
        jnd         = size(lat(:))

!--------------------------------------------------------------------
!>    call annual_mean_solar_2d to calculate the astronomy fields.
!--------------------------------------------------------------------

        call annual_mean_solar(jst, jnd, lat_2d, cosz_2d, solar_2d, fracday_2d, rrsun)

!-------------------------------------------------------------------
!>    place output fields into 1-D arrays for return to calling routine.
!-------------------------------------------------------------------

        solar = solar_2d(:,1)
        cosz  = cosz_2d(:,1)

!--------------------------------------------------------------------
!>    if the calculation has been done, print an error message since
!!    this subroutine should be called only once.
!--------------------------------------------------------------------

    else
        call error_mesg('astronomy_mod', 'annual_mean_solar_2level should be called only once', FATAL)
    endif
    annual_mean_calculated = .true.

end subroutine ANNUAL_MEAN_SOLAR_2LEVEL_


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!> @brief r_inv_squared returns the inverse of the square of the earth-sun
!!        distance relative to the mean distance at angle ang in the Earth's orbit.
function R_INV_SQUARED_(ang)

!--------------------------------------------------------------------
real(kind=FMS_AST_KIND_), intent(in) :: ang !< angular position of earth in its orbit, relative to a
                        !! value of 0.0 at the NH autumnal equinox, value between
                        !! 0.0 and 2 * pi [radians]
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables


real(kind=FMS_AST_KIND_) :: R_INV_SQUARED_ !< The inverse of the square of the earth-sun distance relative
                      !! to the mean distance [dimensionless]
real(kind=FMS_AST_KIND_) :: r             !< Earth-Sun distance relative to mean distance [dimensionless]
real(kind=FMS_AST_KIND_) :: rad_per       !< Angular position of perihelion [radians]
integer, parameter       :: lkind = FMS_AST_KIND_

!--------------------------------------------------------------------
!>    define the earth-sun distance (r) and then return the inverse of
!!    its square (r_inv_squared) to the calling routine.
!--------------------------------------------------------------------

    rad_per        = real(per, FMS_AST_KIND_)*real(deg_to_rad, FMS_AST_KIND_)
    r              = (1.0_lkind - real(ecc, FMS_AST_KIND_)**2)/(1.0_lkind &
                    + real(ecc, FMS_AST_KIND_)*cos(real(ang, FMS_AST_KIND_) - real(rad_per, FMS_AST_KIND_)))
    R_INV_SQUARED_ = r**(-2)


end function R_INV_SQUARED_


!> @brief angle determines the position within the earth's orbit at time t
!!        in the year (t = 0 at NH autumnal equinox) by interpolating
!!        into the orbital position table.

function ANGLE_(t)

!--------------------------------------------------------------------
real(kind=FMS_AST_KIND_), intent(in) :: t !< time of year (between 0 and 2*pi; t=0 at NH autumnal equinox
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------

real(kind=FMS_AST_KIND_) :: ANGLE_ !< Orbital position relative to NH autumnal equinox [radians]
real(kind=FMS_AST_KIND_) :: norm_time !< Index into orbital table corresponding to input time [dimensionless]
real(kind=FMS_AST_KIND_) :: x !< Fractional distance between the orbital table entries bracketing
                              !! the input time [dimensionless]
integer                  :: int !< Table index which is lower than actual position, but closest to it [dimensionless]
integer                  :: int_1 !< Next table index just larger than actual orbital position [dimensionless]
integer, parameter       :: lkind = FMS_AST_KIND_

!--------------------------------------------------------------------
!>    Define orbital tables indices bracketing current orbital time
!!    (int and int_1). Define table index distance between the lower
!!    table value (int) and the actual orbital time (x). Define orbital
!!    position as being x of the way between int and int_1. Renormalize
!!    angle to be within the range 0 to 2*pi.
!--------------------------------------------------------------------

    norm_time = real(t, FMS_AST_KIND_)*real(float(num_angles), FMS_AST_KIND_)/real(twopi, FMS_AST_KIND_)
    int       = floor(norm_time)
    int       = modulo(int,num_angles)
    int_1     = int+1
    x         = norm_time - floor(norm_time)
    ANGLE_    = (1.0_lkind - real(x, FMS_AST_KIND_))*real(orb_angle(int), FMS_AST_KIND_) &
                + real(x, FMS_AST_KIND_)*real(orb_angle(int_1), FMS_AST_KIND_)
    ANGLE_    = modulo(real(ANGLE_, FMS_AST_KIND_), real(twopi, FMS_AST_KIND_))

end function ANGLE_


!> @brief Declination returns the solar declination angle at orbital
!!        position ang in earth's orbit.

function DECLINATION_(ang)

!--------------------------------------------------------------------
real(kind=FMS_AST_KIND_), intent(in) :: ang !< solar orbital position ang in earth's orbit

!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables


real(kind=FMS_AST_KIND_) :: DECLINATION_ !< Solar declination angle [radians]
real(kind=FMS_AST_KIND_) :: rad_obliq !< Obliquity of the ecliptic [radians]
real(kind=FMS_AST_KIND_) :: sin_dec !< Sine of the solar declination [dimensionless]

!---------------------------------------------------------------------
!    compute the solar declination.
!---------------------------------------------------------------------

    rad_obliq    = real(obliq, FMS_AST_KIND_)*real(deg_to_rad, FMS_AST_KIND_)
    sin_dec      = - sin(rad_obliq)*sin(ang)
    DECLINATION_ = asin(sin_dec)

end function DECLINATION_


!> @brief half_day_2d returns a 2-d array of half-day lengths at the
!!        latitudes and declination provided.
!!

function HALF_DAY_2D_(latitude, dec) result(h)

!---------------------------------------------------------------------
real(kind=FMS_AST_KIND_), dimension(:,:), intent(in)                   :: latitude !< Latitutde of view point
real(kind=FMS_AST_KIND_),                 intent(in)                   :: dec !< Solar declination angle at view point
real(kind=FMS_AST_KIND_), dimension(size(latitude,1),size(latitude,2)) :: h
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension (size(latitude,1),size(latitude,2)):: &
                                                  cos_half_day, & !< Cosine of half-day length [dimensionless]
                                                  lat !< Model latitude, adjusted so that it is never 0.5*pi or -0.5*pi
real(kind=FMS_AST_KIND_) :: tan_dec !< tangent of solar declination [dimensionless]
real(kind=FMS_AST_KIND_) :: eps = real(1.0E-05, FMS_AST_KIND_) !< small increment
integer, parameter       :: lkind = FMS_AST_KIND_

!--------------------------------------------------------------------
!>    define tangent of the declination.
!--------------------------------------------------------------------

    tan_dec = tan(dec)

!--------------------------------------------------------------------
!>    adjust latitude so that its tangent will be defined.
!--------------------------------------------------------------------

    lat = latitude
    where (latitude == 0.5_lkind*real(PI, FMS_AST_KIND_))  lat = latitude - eps
    where (latitude == -0.5_lkind*real(PI, FMS_AST_KIND_)) lat = latitude + eps

!--------------------------------------------------------------------
!>    define the cosine of the half-day length. adjust for cases of
!!    all daylight or all night.
!--------------------------------------------------------------------

    cos_half_day = -tan(lat)*tan_dec
    where (cos_half_day <= -1.0_lkind) h = real(PI, FMS_AST_KIND_)
    where (cos_half_day >=  1.0_lkind) h = 0.0_lkind
    where (cos_half_day >  -1.0_lkind .and.  &
           cos_half_day <   1.0_lkind) h = acos(cos_half_day)

end function HALF_DAY_2D_


!> @brief half_day_0d takes scalar input fields, makes them into 2-d fields
!!        dimensioned (1,1), and calls half_day_2d. On return, the 2-d
!!        fields are converted to the desired scalar output.
!!
!! @param [in] <latitude> Latitutde of view point
!! @param [in] <dec> Solar declination angle at view point

function HALF_DAY_0D_(latitude, dec) result(h)

real(kind=FMS_AST_KIND_), intent(in) :: latitude, dec
real(kind=FMS_AST_KIND_)             :: h

!----------------------------------------------------------------------
!  local variables
!----------------------------------------------------------------------

real(kind=FMS_AST_KIND_), dimension(1,1) :: lat_2d, h_2d

!---------------------------------------------------------------------
!    create 2d array from the input latitude field.
!---------------------------------------------------------------------

    lat_2d = latitude

!---------------------------------------------------------------------
!    call half_day with the 2d arguments to calculate half-day length.
!---------------------------------------------------------------------

    h_2d = half_day(lat_2d, dec)

!---------------------------------------------------------------------
!    create scalar from 2d array.
!---------------------------------------------------------------------

    h = h_2d(1,1)

end function HALF_DAY_0D_


!> @}
! close documentation grouping
