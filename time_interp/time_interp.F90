
module time_interp_mod

! <CONTACT EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Computes a weight and dates/indices for linearly interpolating between two dates.
! </OVERVIEW>

! <DESCRIPTION>
!     A time type is converted into two consecutive dates plus
!     a fraction representing the distance between the dates.
!     This information can be used to interpolate between the dates.
!     The dates may be expressed as years, months, or days or
!     as indices in an array.
! </DESCRIPTION>

! <PUBLIC>
!   Description summarizing public interface.
! </PUBLIC>

!-----------------------------------------------------------------------

use time_manager_mod, only: time_type, get_date, set_date, set_time, &
                            days_in_year, days_in_month, leap_year,  &
                            operator(+), operator(-), operator(>),   &
                            operator(<), operator( // ), operator( / ),  &
                            operator(>=), operator(<=)

use          fms_mod, only: write_version_number, &
                            error_mesg, FATAL

implicit none
private

!-----------------------------------------------------------------------

public :: time_interp_init, time_interp, fraction_of_year

! <INTERFACE NAME="time_interp">

!   <OVERVIEW>
!      Returns a weight and dates or indices for interpolating between two dates. The
!      interface fraction_of_year is provided for backward compatibility with the
!      previous version. 
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns the fractional amount (between 0,1) that the input
!      Time is into the current year. The interface fraction_of_year
!      has the same functionality and is only provided for backward
!      compatibility with the version. 
!   </DESCRIPTION>
!   <TEMPLATE>
!      call time_interp( Time, weight )
!   </TEMPLATE>
!   <TEMPLATE>
!      call time_interp( Time, weight, year1, year2 )
!   </TEMPLATE>
!   <TEMPLATE>
!      call time_interp( Time, weight, year1, year2, month1, month2 )
!   </TEMPLATE>
!   <TEMPLATE>
!      call time_interp( Time, weight, year1, year2, month1, month2, day1, day2 )
!   </TEMPLATE>
!   <TEMPLATE>
!      call time_interp( Time, weightTime, Timelist, weight, index1, index2 [, modtime] )
!   </TEMPLATE>
!   <IN NAME="Time">
!      The time at which the the weight is computed.
!   </IN>
!   <IN NAME="Timelist">
!   </IN>
!   <IN NAME="modtime">
!   </IN>
!   <OUT NAME="weight">
!     The fractional amount (between 0,1) into the year as given by argument Time.
!   </OUT>
!   <OUT NAME="year1"> </OUT>
!   <OUT NAME="year2"> </OUT>
!   <OUT NAME="month1"> </OUT>
!   <OUT NAME="month2"> </OUT>
!   <OUT NAME="day1"> </OUT>
!   <OUT NAME="day2"> </OUT>
!   <OUT NAME="index1"> </OUT>
!   <OUT NAME="index2"> </OUT>
!   <ERROR MSG="input time list not ascending order" STATUS="ERROR">
!     The list of input time types must have ascending dates.
!   </ERROR>
!   <ERROR MSG="modulo months must have same length" STATUS="ERROR">
!     The length of the current month for input Time and Time_list
!     must be the same when using the modulo month option. The
!     modulo month option is available but not supported. 
!   </ERROR>
!   <ERROR MSG="invalid value for argument modtime" STATUS="ERROR">
!     The optional argument modtime must have a value set by one
!     of the public parameters: NONE, YEAR, MONTH, DAY. The
!     MONTH and DAY options are available but not supported. 
!   </ERROR>
!   <ERROR MSG="period of list exceeds modulo period" STATUS="ERROR">
!     The difference between the last and first values in the input
!     Time list/array exceeds the length of the modulo period.
!   </ERROR>
!   <ERROR MSG="time before range of list or time after range of list" STATUS="ERROR">
!     The difference between the last and first values in the input
!     These errors occur when you are not using a modulo axis and
!     the input Time occurs before the first value in the Time
!     list/array or after the last value in the Time list/array. 
!   </ERROR>
!   <NOTE>
!     Examples: 
!     <PRE>
!       Time: Jan 01 00z    weight = 0.0 
!       Time: Jul 01        weight ~ 0.5 
!       Time: Dec 31 23z    weight ~ 1.0
!     </PRE>
!   </NOTE>

interface time_interp
    module procedure time_interp_frac,  time_interp_year, &
                     time_interp_month, time_interp_day,  &
                     time_interp_list
end interface
! </INTERFACE>

integer, public, parameter :: NONE=0, YEAR=1, MONTH=2, DAY=3

!-----------------------------------------------------------------------

   integer, parameter ::  secmin = 60, minhour = 60, hourday = 24,  &
                         sechour = secmin*minhour,                  &
                          secday = secmin*minhour*hourday

   integer, parameter :: monyear = 12
   integer, parameter :: halfday = secday/2

   integer :: mtime
   integer :: yrmod, momod, dymod
   logical :: mod_leapyear

   character(len=128) :: version='$Id: time_interp.F90,v 1.4 2003/04/09 21:19:06 fms Exp $'
   character(len=128) :: tagname='$Name: inchon $'

   logical :: module_is_initialized=.FALSE.

contains


 subroutine time_interp_init()

   if ( module_is_initialized ) return

   call write_version_number( version, tagname )

   module_is_initialized = .TRUE.

 end subroutine time_interp_init

!#######################################################################

! <SUBROUTINE NAME="time_interp_frac" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
! </SUBROUTINE>
!  returns the fractional time into the current year

 subroutine time_interp_frac ( Time, weight )

   type(time_type), intent(in)  :: Time 
   real           , intent(out) :: weight

   integer         :: year, month, day, hour, minute, second
   type(time_type) :: Year_beg, Year_end


   if ( .not. module_is_initialized ) call time_interp_init

!  ---- compute fractional time of year -----

     call get_date (Time, year, month, day, hour, minute, second) 

     Year_beg = set_date(year  , 1, 1) 
     Year_end = set_date(year+1, 1, 1)

     weight = (Time - Year_beg) // (Year_end - Year_beg)

 end subroutine time_interp_frac

!#######################################################################
! <SUBROUTINE NAME="fraction_of_year">
! <OVERVIEW>
!  Wrapper for backward compatibility
! </OVERVIEW>
! </SUBROUTINE>

 function fraction_of_year (Time)
 type(time_type), intent(in)  :: Time
 real :: fraction_of_year

  call time_interp_frac ( Time, fraction_of_year )

 end function fraction_of_year

!#######################################################################
! <SUBROUTINE NAME="time_interp_year" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
!   <OUT NAME="year1" TYPE="integer"> </OUT>
!   <OUT NAME="year2" TYPE="integer"> </OUT>
! </SUBROUTINE>
!  returns fractional time between mid points of consecutive years

 subroutine time_interp_year ( Time, weight, year1, year2 )

   type(time_type), intent(in)  :: Time
   real           , intent(out) :: weight
   integer        , intent(out) :: year1, year2

   integer :: year, month, day, hour, minute, second
   type (time_type) :: Mid_year, Mid_year1, Mid_year2


   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, year, month, day, hour, minute, second)

    ! mid point of current year
      Mid_year = year_midpt(year)

      if ( Time >= Mid_year ) then
    ! current time is after mid point of current year
           year1  = year
           year2  = year+1
           Mid_year2 = year_midpt(year2)
           weight = (Time - Mid_year) // (Mid_year2 - Mid_year)
      else
    ! current time is before mid point of current year
           year2  = year
           year1  = year-1
           Mid_year1 = year_midpt(year1)
           weight = (Time - Mid_year1) // (Mid_year - Mid_year1)
      endif

 end subroutine time_interp_year

!#######################################################################
! <SUBROUTINE NAME="time_interp_month" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
!   <OUT NAME="year1" TYPE="integer"> </OUT>
!   <OUT NAME="year2" TYPE="integer"> </OUT>
!   <OUT NAME="month1" TYPE="integer"> </OUT>
!   <OUT NAME="month2" TYPE="integer"> </OUT>
! </SUBROUTINE>
!  returns fractional time between mid points of consecutive months

 subroutine time_interp_month ( Time, weight, year1, year2, month1, month2 )

   type(time_type), intent(in)  :: Time
   real           , intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2

   integer :: year, month, day, hour, minute, second,  &
              mid_month, cur_month, mid1, mid2

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, year, month, day, hour, minute, second)

    ! mid point of current month in seconds
      mid_month = days_in_month(Time) * halfday
    ! time into current month in seconds
      cur_month = second + secmin*minute + sechour*hour + secday*(day-1)

      if ( cur_month >= mid_month ) then
    ! current time is after mid point of current month
           year1  = year;  month1 = month
           year2  = year;  month2 = month+1
           if (month2 > monyear)  year2 = year2+1
           if (month2 > monyear) month2 = 1
           mid1 = mid_month
           mid2 = days_in_month(set_date(year2,month2,2)) * halfday
           weight = real(cur_month - mid1) / real(mid1+mid2)
      else
    ! current time is before mid point of current month
           year2  = year;  month2 = month
           year1  = year;  month1 = month-1
           if (month1 < 1)  year1 = year1-1
           if (month1 < 1) month1 = monyear
           mid1 = days_in_month(set_date(year1,month1,2)) * halfday
           mid2 = mid_month
           weight = real(cur_month + mid1) / real(mid1+mid2)
      endif

 end subroutine time_interp_month

!#######################################################################
! <SUBROUTINE NAME="time_interp_day" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
!   <OUT NAME="year1" TYPE="integer"> </OUT>
!   <OUT NAME="year2" TYPE="integer"> </OUT>
!   <OUT NAME="month1" TYPE="integer"> </OUT>
!   <OUT NAME="month2" TYPE="integer"> </OUT>
!   <OUT NAME="day1" TYPE="integer"> </OUT>
!   <OUT NAME="day2" TYPE="integer"> </OUT>
! </SUBROUTINE>
!  returns fractional time between mid points of consecutive days

 subroutine time_interp_day ( Time, weight, year1, year2, month1, month2, day1, day2 )

   type(time_type), intent(in)  :: Time
   real           , intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2, day1, day2

   integer :: year, month, day, hour, minute, second, sday

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, year, month, day, hour, minute, second)

    ! time into current day in seconds
      sday = second + secmin*minute + sechour*hour

      if ( sday >= halfday ) then
    ! current time is after mid point of day
           year1 = year;  month1 = month;  day1 = day
           year2 = year;  month2 = month;  day2 = day + 1
           weight  = real(sday - halfday) / real(secday)

           if (day2 > days_in_month(Time)) then
               month2 = month2 + 1
               day2 = 1
               if (month2 > monyear) then
                    month2 = 1;  year2 = year2+1
               endif
           endif
      else
    ! current time is before mid point of day
           year2 = year;  month2 = month;  day2 = day
           year1 = year;  month1 = month;  day1 = day - 1
           weight  = real(sday + halfday) / real(secday)

           if (day1 < 1) then
               month1 = month1 - 1
               if (month1 < 1) then
                   month1 = monyear;  year1 = year1-1
               endif
               day1 = days_in_month(set_date(year1,month1,2))
           endif
      endif

 end subroutine time_interp_day

!#######################################################################
! <SUBROUTINE NAME="time_interp_list" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <IN NAME="Timelist" TYPE="time_type" DIM="(:)"> </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
!   <OUT NAME="index1" TYPE="real"> </OUT>
!   <OUT NAME="index2" TYPE="real"> </OUT>
!   <IN NAME="modtime" TYPE="integer" > </IN>
! </SUBROUTINE>

subroutine time_interp_list ( Time, Timelist, weight, index1, index2, modtime )
type(time_type)  , intent(in)  :: Time, Timelist(:)
real             , intent(out) :: weight
integer          , intent(out) :: index1, index2
integer, optional, intent(in)  :: modtime

integer :: i, n, hr, mn, se
type(time_type) :: T, T1, T2, Ts, Te, Td, Period, Time_mod

  weight = 0.; index1 = 0; index2 = 0
  n = size(Timelist)

! check list for ascending order
  do i = 2, n
     if (Timelist(i) > Timelist(i-1)) cycle
     call error_handler ('input time list not ascending order')
  enddo

! setup modular time axis?
  mtime = NONE
  if (present(modtime)) then
     mtime = modtime
     Time_mod = (Timelist(1)+Timelist(n))/2
     call get_date (Time_mod, yrmod, momod, dymod, hr, mn, se)
     mod_leapyear = leap_year(Time_mod)
  endif

! starting and ending times from list
  Ts = set_modtime(Timelist(1))
  Te = set_modtime(Timelist(n))
  Td = Te-Ts
  T  = set_modtime(Time)

! set period for modulo axis
  select case (mtime)
     case (NONE)
       ! do nothing
     case (YEAR)
         Period = set_time(0,days_in_year(Time_mod))
     case (MONTH)
       ! month length must be equal
         if (days_in_month(Time_mod) /= days_in_month(Time)) &
         call error_handler ('modulo months must have same length')
         Period = set_time(0,days_in_month(Time_mod))
     case (DAY)
         Period = set_time(0,1)
     case default
         call error_handler ('invalid value for argument modtime')
  end select

! check length of modulo period
! list period cannot exceed modulo period
  if (mtime /= NONE) then
      if (Td >= Period) call error_handler &
                         ('period of list exceeds modulo period')
  endif

! time falls between start and end list values
  if ( T >= Ts .and. T <= Te ) then
     T1 = set_modtime(Timelist(1))
     do i = 2, n
       T2 = set_modtime(Timelist(i))
       if ( T >= T1 .and. T <= T2 ) then
          index1 = i-1
          index2 = i
          weight = (T-T1) // (T2-T1)
          exit
       endif
       T1 = T2
     enddo

! time falls before starting list value
  else if ( T < Ts ) then
     if (mtime == NONE) call error_handler ('time before range of list')
     Td = Te-Ts
     weight = 1. - ((Ts-T) // (Period-Td))
     index1 = n
     index2 = 1

! time falls after ending list value
  else if ( T > Te ) then
     if (mtime == NONE) call error_handler ('time after range of list')
     Td = Te-Ts
     weight = (T-Te) // (Period-Td)
     index1 = n
     index2 = 1
  endif

end subroutine time_interp_list

!#######################################################################
!  private routines
!#######################################################################

 function year_midpt (year)

   integer, intent(in) :: year
   type (time_type)    :: year_midpt, year_beg, year_end


   year_beg = set_date(year  , 1, 1)
   year_end = set_date(year+1, 1, 1)

   year_midpt = (year_beg + year_end) / 2
   
 end function year_midpt

!#######################################################################

 function month_midpt (year, month)

   integer, intent(in) :: year, month
   type (time_type)    :: month_midpt, month_beg, month_end

!  --- beginning of this month ---
   month_beg = set_date(year, month, 1)

!  --- start of next month ---
   if (month < 12) then
      month_end = set_date(year, month+1, 1)
   else
      month_end = set_date(year+1, 1, 1)
   endif

   month_midpt = (month_beg + month_end) / 2
   
 end function month_midpt

!#######################################################################

function set_modtime (Tin) result (Tout)
type(time_type), intent(in) :: Tin
type(time_type)             :: Tout
integer :: yr, mo, dy, hr, mn, se

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

!#######################################################################

subroutine error_handler (string)
character(len=*), intent(in) :: string

  call error_mesg ('time_interp_mod', trim(string), FATAL)

! write (*,'(a)') 'ERROR in time_interp: ' // trim(string)
! stop 111

end subroutine error_handler

!#######################################################################

end module time_interp_mod

! <INFO>

!   <ERROR MSG="input time list not ascending order" STATUS="">
!     The list of input time types must have ascending dates.
!   </ERROR> *
!   <ERROR MSG="modulo months must have same length" STATUS="">
!     The length of the current month for input Time and Time_list
!     must be the same when using the modulo month option.
!     The modulo month option is available but not supported.
!   </ERROR> *
!   <ERROR MSG="invalid value for argument modtime" STATUS="">
!     The optional argument modtime must have a value set by one
!     of the public parameters: NONE, YEAR, MONTH, DAY.
!     The MONTH and DAY options are available but not supported.
!   </ERROR> *
!   <ERROR MSG="period of list exceeds modulo period" STATUS="">
!     The difference between the last and first values in the
!     input Time list/array exceeds the length of the modulo period.
!   </ERROR> *
!   <ERROR MSG="time before range of list or time after range of list" STATUS="">
!     These errors occur when you are not using a modulo axis and the
!     input Time occurs before the first value in the Time list/array
!     or after the last value in the Time list/array.
!   </ERROR> *
!   <NOTE>
!   For all routines in this module the calendar type in module
!   time_manager must be set.
!   </NOTE>
!   <NOTE>
!     The following private parameters are set by this module:
! <PRE>
!           seconds per minute = 60
!           minutes per hour   = 60
!           hours   per day    = 24
!           months  per year   = 12
! </PRE>
!   </NOTE>

! </INFO>
