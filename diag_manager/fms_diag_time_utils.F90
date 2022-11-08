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
!> @defgroup fms_diag_time_utils_mod fms_diag_time_utils_mod
!> @ingroup diag_manager
!! @brief fms_diag_time_utils contains functions and subroutines necessary for the
!! <TT>diag_manager_mod</TT> related to time handling.
!! @author Uriel Ramirez

!> @addtogroup fms_diag_time_utils_mod
!> @{
module fms_diag_time_utils_mod

use time_manager_mod, only: time_type, increment_date, increment_time, get_calendar_type, NO_CALENDAR, leap_year, &
                            get_date, get_time,  operator(>), operator(<), operator(-)
use diag_data_mod,    only: END_OF_RUN, EVERY_TIME, DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, &
                            DIAG_YEARS
USE constants_mod,    ONLY: SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_MINUTE
use fms_mod,          only: fms_error_handler, error_mesg, FATAL

implicit none
private

public :: diag_time_inc
public :: get_time_string
public :: get_date_dif

contains

  !> @brief Return the next time data/file is to be written based on the frequency and units.
  TYPE(time_type) FUNCTION diag_time_inc(time, output_freq, output_units, err_msg)
    TYPE(time_type),  INTENT(in)            :: time         !< Current model time.
    INTEGER,          INTENT(in)            :: output_freq  !< Output frequency number value.
    INTEGER,          INTENT(in)            :: output_units !< Output frequency unit.
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg      !< Function error message.
                                                            !! An empty string indicates the next output
                                                            !! time was found successfully.

    CHARACTER(len=128) :: error_message_local !< Local variable to store the error_message

    IF ( PRESENT(err_msg) ) err_msg = ''
    error_message_local = ''

    ! special values for output frequency are -1 for output at end of run
    ! and 0 for every timestep.  Need to check for these here?
    ! Return zero time increment, hopefully this value is never used
    IF ( output_freq == END_OF_RUN .OR. output_freq == EVERY_TIME ) THEN
       diag_time_inc = time
       RETURN
    END IF

    ! Make sure calendar was not set after initialization
    IF ( output_units == DIAG_SECONDS ) THEN
       IF ( get_calendar_type() == NO_CALENDAR ) THEN
          diag_time_inc = increment_time(time, output_freq, 0, err_msg=error_message_local)
       ELSE
          diag_time_inc = increment_date(time, 0, 0, 0, 0, 0, output_freq, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_MINUTES ) THEN
       IF ( get_calendar_type() == NO_CALENDAR ) THEN
          diag_time_inc = increment_time(time, NINT(output_freq*SECONDS_PER_MINUTE), 0, &
               &err_msg=error_message_local)
       ELSE
          diag_time_inc = increment_date(time, 0, 0, 0, 0, output_freq, 0, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_HOURS ) THEN
       IF ( get_calendar_type() == NO_CALENDAR ) THEN
          diag_time_inc = increment_time(time, NINT(output_freq*SECONDS_PER_HOUR), 0, err_msg=error_message_local)
       ELSE
          diag_time_inc = increment_date(time, 0, 0, 0, output_freq, 0, 0, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_DAYS ) THEN
       IF (get_calendar_type() == NO_CALENDAR) THEN
          diag_time_inc = increment_time(time, 0, output_freq, err_msg=error_message_local)
       ELSE
          diag_time_inc = increment_date(time, 0, 0, output_freq, 0, 0, 0, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_MONTHS ) THEN
       IF (get_calendar_type() == NO_CALENDAR) THEN
          error_message_local = 'output units of months NOT allowed with no calendar'
       ELSE
          diag_time_inc = increment_date(time, 0, output_freq, 0, 0, 0, 0, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_YEARS ) THEN
       IF ( get_calendar_type() == NO_CALENDAR ) THEN
          error_message_local = 'output units of years NOT allowed with no calendar'
       ELSE
          diag_time_inc = increment_date(time, output_freq, 0, 0, 0, 0, 0, err_msg=error_message_local)
       END IF
    ELSE
       error_message_local = 'illegal output units'
    END IF

    IF ( error_message_local /= '' ) THEN
      IF ( fms_error_handler('diag_time_inc',error_message_local,err_msg) ) RETURN
    END IF
  END FUNCTION diag_time_inc

  !> @brief This function determines a string based on current time.
  !!     This string is used as suffix in output file name
  !! @return Character(len=128) get_time_string
  CHARACTER(len=128) FUNCTION get_time_string(filename, current_time)
    CHARACTER(len=128), INTENT(in) :: filename     !< File name.
    TYPE(time_type),    INTENT(in) :: current_time !< Current model time.

    INTEGER :: yr1 !< get from current time
    INTEGER :: mo1 !< get from current time
    INTEGER :: dy1 !< get from current time
    INTEGER :: hr1 !< get from current time
    INTEGER :: mi1 !< get from current time
    INTEGER :: sc1 !< get from current time
    INTEGER :: yr2 !< for computing next_level time unit
    INTEGER :: dy2 !< for computing next_level time unit
    INTEGER :: hr2 !< for computing next_level time unit
    INTEGER :: mi2 !< for computing next_level time unit
    INTEGER :: yr1_s !< actual values to write string
    INTEGER :: mo1_s !< actual values to write string
    INTEGER :: dy1_s !< actual values to write string
    INTEGER :: hr1_s !< actual values to write string
    INTEGER :: mi1_s !< actual values to write string
    INTEGER :: sc1_s !< actual values to write string
    INTEGER :: abs_day              !< component of current_time
    INTEGER :: abs_sec              !< component of current_time
    INTEGER :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    INTEGER :: julian_day, i, position, len, first_percent
    CHARACTER(len=1) :: width  !< width of the field in format write
    CHARACTER(len=10) :: format
    CHARACTER(len=20) :: yr !< string of current time (output)
    CHARACTER(len=20) :: mo !< string of current time (output)
    CHARACTER(len=20) :: dy !< string of current time (output)
    CHARACTER(len=20) :: hr !< string of current time (output)
    CHARACTER(len=20) :: mi !< string of current time (output)
    CHARACTER(len=20) :: sc !< string of current time (output)
    CHARACTER(len=128) :: filetail

    format = '("_",i*.*)'
    CALL get_date(current_time, yr1, mo1, dy1, hr1, mi1, sc1)
    len = LEN_TRIM(filename)
    first_percent = INDEX(filename, '%')
    filetail = filename(first_percent:len)
    ! compute year string
    position = INDEX(filetail, 'yr')
    IF ( position > 0 ) THEN
       width = filetail(position-1:position-1)
       yr1_s = yr1
       format(7:9) = width//'.'//width
       WRITE(yr, format) yr1_s
       yr2 = 0
    ELSE
       yr = ' '
       yr2 = yr1 - 1
    END IF
    ! compute month string
    position = INDEX(filetail, 'mo')
    IF ( position > 0 ) THEN
       width = filetail(position-1:position-1)
       mo1_s = yr2*12 + mo1
       format(7:9) = width//'.'//width
       WRITE(mo, format) mo1_s
    ELSE
       mo = ' '
    END IF
    ! compute day string
    IF ( LEN_TRIM(mo) > 0 ) THEN ! month present
       dy1_s = dy1
       dy2 = dy1_s - 1
    ELSE IF ( LEN_TRIM(yr) >0 )  THEN ! no month, year present
       ! compute julian day
       IF ( mo1 == 1 ) THEN
          dy1_s = dy1
       ELSE
          julian_day = 0
          DO i = 1, mo1-1
             julian_day = julian_day + days_per_month(i)
          END DO
          IF ( leap_year(current_time) .AND. mo1 > 2 ) julian_day = julian_day + 1
          julian_day = julian_day + dy1
          dy1_s = julian_day
       END IF
       dy2 = dy1_s - 1
    ELSE ! no month, no year
       CALL get_time(current_time, abs_sec, abs_day)
       dy1_s = abs_day
       dy2 = dy1_s
    END IF
    position = INDEX(filetail, 'dy')
    IF ( position > 0 ) THEN
       width = filetail(position-1:position-1)
       FORMAT(7:9) = width//'.'//width
       WRITE(dy, FORMAT) dy1_s
    ELSE
       dy = ' '
    END IF
    ! compute hour string
    IF ( LEN_TRIM(dy) > 0 ) THEN
       hr1_s = hr1
    ELSE
       hr1_s = dy2*24 + hr1
    END IF
    hr2 = hr1_s
    position = INDEX(filetail, 'hr')
    IF ( position > 0 ) THEN
       width = filetail(position-1:position-1)
       format(7:9) = width//'.'//width
       WRITE(hr, format) hr1_s
    ELSE
       hr = ' '
    END IF
    ! compute minute string
    IF ( LEN_TRIM(hr) > 0 ) THEN
       mi1_s = mi1
    ELSE
       mi1_s = hr2*60 + mi1
    END IF
    mi2 = mi1_s
    position = INDEX(filetail, 'mi')
    IF(position>0) THEN
       width = filetail(position-1:position-1)
       format(7:9) = width//'.'//width
       WRITE(mi, format) mi1_s
    ELSE
       mi = ' '
    END IF
    ! compute second string
    IF ( LEN_TRIM(mi) > 0 ) THEN
       sc1_s = sc1
    ELSE
       sc1_s = NINT(mi2*SECONDS_PER_MINUTE) + sc1
    END IF
    position = INDEX(filetail, 'sc')
    IF ( position > 0 ) THEN
       width = filetail(position-1:position-1)
       format(7:9) = width//'.'//width
       WRITE(sc, format) sc1_s
    ELSE
       sc = ' '
    ENDIF
    get_time_string = TRIM(yr)//TRIM(mo)//TRIM(dy)//TRIM(hr)//TRIM(mi)//TRIM(sc)
  END FUNCTION get_time_string

  !> @brief Return the difference between two times in units.
  !! @return Real get_data_dif
  REAL FUNCTION get_date_dif(t2, t1, units)
    TYPE(time_type), INTENT(in) :: t2 !< Most recent time.
    TYPE(time_type), INTENT(in) :: t1 !< Most distant time.
    INTEGER, INTENT(in) :: units !< Unit of return value.

    INTEGER :: dif_seconds, dif_days
    TYPE(time_type) :: dif_time

    ! Compute time axis label value
    ! <ERROR STATUS="FATAL">
    !   variable t2 is less than in variable t1
    ! </ERROR>
    IF ( t2 < t1 ) CALL error_mesg('diag_util_mod::get_date_dif', &
         & 'in variable t2 is less than in variable t1', FATAL)

    dif_time = t2 - t1

    CALL get_time(dif_time, dif_seconds, dif_days)

    IF ( units == DIAG_SECONDS ) THEN
       get_date_dif = dif_seconds + SECONDS_PER_DAY * dif_days
    ELSE IF ( units == DIAG_MINUTES ) THEN
       get_date_dif = 1440 * dif_days + dif_seconds / SECONDS_PER_MINUTE
    ELSE IF ( units == DIAG_HOURS ) THEN
       get_date_dif = 24 * dif_days + dif_seconds / SECONDS_PER_HOUR
    ELSE IF ( units == DIAG_DAYS ) THEN
       get_date_dif = dif_days + dif_seconds / SECONDS_PER_DAY
    ELSE IF ( units == DIAG_MONTHS ) THEN
       ! <ERROR STATUS="FATAL">
       !   months not supported as output units
       ! </ERROR>
       CALL error_mesg('diag_util_mod::get_date_dif', 'months not supported as output units', FATAL)
    ELSE IF ( units == DIAG_YEARS ) THEN
       ! <ERROR STATUS="FATAL">
       !   years not suppored as output units
       ! </ERROR>
       CALL error_mesg('diag_util_mod::get_date_dif', 'years not supported as output units', FATAL)
    ELSE
       ! <ERROR STATUS="FATAL">
       !   illegal time units
       ! </ERROR>
       CALL error_mesg('diag_util_mod::diag_date_dif', 'illegal time units', FATAL)
    END IF
  END FUNCTION get_date_dif
end module fms_diag_time_utils_mod
