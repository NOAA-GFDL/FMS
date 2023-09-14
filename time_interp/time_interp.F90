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
!> @defgroup time_interp_mod time_interp_mod
!> @ingroup time_interp
!> @brief Computes a weight and dates/indices for linearly interpolating between two dates.
!> @author Bruce Wyman
!!
!! A time type is converted into two consecutive dates plus
!! a fraction representing the distance between the dates.
!! This information can be used to interpolate between the dates.
!! The dates may be expressed as years, months, or days or
!! as indices in an array.

module time_interp_mod

use time_manager_mod, only: time_type, get_date, set_date, set_time, &
                            days_in_year, days_in_month, leap_year,  &
                            time_type_to_real, real_to_time_type,    &
                            get_calendar_type, JULIAN, GREGORIAN, NO_CALENDAR, &
                            operator(+), operator(-), operator(>),   &
                            operator(<), operator( // ), operator( / ),  &
                            operator(>=), operator(<=), operator( * ), &
                            operator(==), print_date, print_time,&
                            time_list_error, date_to_string

use          fms_mod, only: write_version_number, &
                            error_mesg, FATAL, stdout, stdlog, &
                            check_nml_error, &
                            fms_error_handler
use          mpp_mod, only: input_nml_file
use     platform_mod

implicit none
private

!-----------------------------------------------------------------------

public :: time_interp_init, time_interp, fraction_of_year

!> Returns a weight and dates or indices for interpolating between two dates. The
!! interface fraction_of_year is provided for backward compatibility with the
!! previous version.\n
!!
!! Returns weight by interpolating Time between Time1 and Time2.
!! i.e. weight = (Time-Time1)/(Time2-Time1)
!! Time1 and Time2 may be specified by any of several different ways,
!! which is the reason for multiple interfaces.\n
!!
!! - If Time1 and Time2 are the begining and end of the year in which
!!   Time falls, use first interface.\n
!!
!! - If Time1 and Time2 fall on year boundaries, use second interface.\n
!!
!! - If Time1 and Time2 fall on month boundaries, use third.\n
!!
!! - If Time1 and Time2 fall on day boundaries, use fourth.\n
!!
!! - If Time1 and Time2 are consecutive elements of an assending list, use fifth.
!!   The fifth also returns the indices of Timelist between which Time falls.\n
!!
!! - The sixth interface is for cyclical data. Time_beg and Time_end specify the
!!   begining and end of a repeating period. In this case:<br>
!!              weight = (Time_adjusted - Time1) / (Time2 - Time1)
!! <br>Where:
!! @code{.F90}
!!              Time1 = Timelist(index1)
!!              Time2 = Timelist(index2)
!!              Time_adjusted = Time - N*Period
!!              Period = Time_end-Time_beg
!! @endcode
!! N is between (Time-Time_end)/Period and (Time-Time_beg)/Period
!! That is, N is the integer that results in Time_adjusted that is between Time_beg and Time_end.
!!
!! <br>Example usages:
!! @code{.F90}
!!              call time_interp( Time, weight )
!!              call time_interp( Time, weight, year1, year2 )
!!              call time_interp( Time, weight, year1, year2, month1, month2 )
!!              call time_interp( Time, weight, year1, year2, month1, month2, day1, day2 )
!!              call time_interp( Time, Timelist, weight, index1, index2 [, modtime] )
!!              call time_interp( Time, Time_beg, Time_end, Timelist, weight, index1, index2
!!              [,correct_leap_year_inconsistency])
!! @endcode
!!
!!   For all routines in this module the calendar type in module
!!   time_manager must be set.
!!
!!     The following private parameters are set by this module:
!!
!!          seconds per minute = 60
!!          minutes per hour   = 60
!!          hours   per day    = 24
!!          months  per year   = 12
!!
!! @param Time The time at which the the weight is computed.
!! @param Time_beg For cyclical interpolation: Time_beg specifies the begining time of a cycle.
!! @param Time_end For cyclical interpolation: Time_end specifies the ending time of a cycle.
!! @param Timelist For cyclical interpolation: Timelist is an array of times between Time_beg and Time_end.
!!                 Must be monotonically increasing.
!! @param index1 Timelist(index1) = The largest value of Timelist which is less than mod(Time,Time_end-Time_beg)
!! @param index2 Timelist(index2) = The smallest value of Timelist which is greater than mod(Time,Time_end-Time_beg)
!! @param correct_leap_year_inconsistency Turns on a kluge for an inconsistency which may occur in a special case.
!!       When the modulo time period (i.e. Time_end - Time_beg) is a whole number of years
!!       and is not a multiple of 4, and the calendar in use has leap years, then it is
!!       likely that the interpolation will involve mapping a common year onto a leap year.
!!       In this case it is often desirable, but not absolutely necessary, to use data for
!!       Feb 28 of the leap year when it is mapped onto a common year.
!!       To turn this on, set correct_leap_year_inconsistency=.true.
!! @param weight weight = (mod(Time,Time_end-Time_beg) - Timelist(index1)) / (Timelist(index2) - Timelist(index1))
!> @ingroup time_interp_mod
interface time_interp
    module procedure time_interp_frac_r8,  time_interp_year_r8, &
                     time_interp_month_r8, time_interp_day_r8,  &
                     time_interp_list_r8,  time_interp_modulo_r8
    module procedure time_interp_frac_r4,  time_interp_year_r4, &
                     time_interp_month_r4, time_interp_day_r4,  &
                     time_interp_list_r4,  time_interp_modulo_r4
end interface

!> @addtogroup time_interp_mod
!> @{
integer, public, parameter :: NONE=0, YEAR=1, MONTH=2, DAY=3

!-----------------------------------------------------------------------

   integer, parameter ::  secmin = 60, minhour = 60, hourday = 24,  &
                          sechour = secmin*minhour,                  &
                          secday = secmin*minhour*hourday

   integer, parameter :: monyear = 12
   integer, parameter :: halfday = secday/2

   integer :: yrmod, momod, dymod
   logical :: mod_leapyear

! Include variable "version" to be written to log file.
#include<file_version.h>

   logical :: module_is_initialized=.FALSE.
   logical :: perthlike_behavior=.FALSE.

   namelist / time_interp_nml / perthlike_behavior

contains


 subroutine time_interp_init()
   integer :: ierr, io, logunit

   if ( module_is_initialized ) return

   read (input_nml_file, time_interp_nml, iostat=io)
   ierr = check_nml_error (io, 'time_interp_nml')

   call write_version_number("TIME_INTERP_MOD", version)
   logunit = stdlog()
   write(logunit,time_interp_nml)

   module_is_initialized = .TRUE.

 end subroutine time_interp_init


!> @brief Wrapper function to return the fractional time into the current year
!! Always returns an r8_kind, conversion to r4 will be done implicitly if needed
!> @param Time time to calculate fraction with
!> @return real(kind=8) fraction of time passed in current year
 function fraction_of_year (Time)
 type(time_type), intent(in)  :: Time
 real(r8_kind) :: fraction_of_year

  call time_interp ( Time, fraction_of_year )

 end function fraction_of_year

!#######################################################################
!> Given an array of times in ascending order and a specific time returns
!! values of index1 and index2 such that the Timelist(index1)<=Time and
!! Time<=Timelist(index2), and index2=index1+1
!! index1=0, index2=1 or index=n, index2=n+1 are returned to indicate that
!! the time is out of range
subroutine bisect(Timelist,Time,index1,index2)
  type(time_type)  , intent(in)  :: Timelist(:)
  type(time_type)  , intent(in)  :: Time
  integer, optional, intent(out) :: index1, index2

  integer :: i,il,iu,n,i1,i2

  n = size(Timelist(:))

  if (Time==Timelist(1)) then
     i1 = 1 ; i2 = 2
  else if (Time==Timelist(n)) then
     i1 = n ; i2 = n+1
  else
     il = 0; iu=n+1
     do while(iu-il > 1)
        i = (iu+il)/2
        if(Timelist(i) > Time) then
           iu = i
        else
           il = i
        endif
     enddo
     i1 = il ; i2 = il+1
  endif

  if(PRESENT(index1)) index1 = i1
  if(PRESENT(index2)) index2 = i2
end subroutine bisect

!#######################################################################
!  private routines
!#######################################################################

 function year_midpt (yr)

   integer, intent(in) :: yr
   type (time_type)    :: year_midpt, year_beg, year_end


   year_beg = set_date(yr  , 1, 1)
   year_end = set_date(yr+1, 1, 1)

   year_midpt = (year_beg + year_end) / 2

 end function year_midpt

 function month_midpt (yr, mo)

   integer, intent(in) :: yr, mo
   type (time_type)    :: month_midpt, month_beg, month_end

!  --- beginning of this month ---
   month_beg = set_date(yr, mo, 1)

!  --- start of next month ---
   if (mo < 12) then
      month_end = set_date(yr, mo+1, 1)
   else
      month_end = set_date(yr+1, 1, 1)
   endif

   month_midpt = (month_beg + month_end) / 2

 end function month_midpt

function set_modtime (Tin, modtime) result (Tout)
type(time_type), intent(in) :: Tin
integer, intent(in), optional :: modtime
type(time_type)             :: Tout
integer :: yr, mo, dy, hr, mn, se, mtime

  if(present(modtime)) then
    mtime = modtime
  else
    mtime = NONE
  endif

  select case (mtime)
    case (NONE)
       Tout = Tin
    case (YEAR)
       call get_date (Tin, yr, mo, dy, hr, mn, se)
       yr = yrmod
        ! correct leap year dates
          if (.not.mod_leapyear .and. mo == 2 .and. dy > 28) then
             mo = 3; dy = dy-28
          endif
       Tout = set_date (yr, mo, dy, hr, mn, se)
    case (MONTH)
       call get_date (Tin, yr, mo, dy, hr, mn, se)
       yr = yrmod; mo = momod
       Tout = set_date (yr, mo, dy, hr, mn, se)
    case (DAY)
       call get_date (Tin, yr, mo, dy, hr, mn, se)
       yr = yrmod; mo = momod; dy = dymod
       Tout = set_date (yr, mo, dy, hr, mn, se)
  end select

end function set_modtime

subroutine error_handler(string)
  character(len=*), intent(in) :: string

  call error_mesg ('time_interp_mod', trim(string), FATAL)

end subroutine error_handler


#include "time_interp_r4.fh"
#include "time_interp_r8.fh"

end module time_interp_mod

!> @}
! close documentation grouping
