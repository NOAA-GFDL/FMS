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

!> @addtogroup time_interp_mod
!> @{

 !> @brief Calculates the fractional time into the current year
 subroutine TIME_INTERP_FRAC_ ( Time, weight )

   type(time_type),   intent(in)  :: Time
   real(FMS_TI_KIND_), intent(out) :: weight !< fractional time

   integer         :: yr, mo, dy, hour, minute, second
   type(time_type) :: Year_beg, Year_end


   if ( .not. module_is_initialized ) call time_interp_init

!  ---- compute fractional time of year -----

     call get_date (Time, yr, mo, dy, hour, minute, second)

     Year_beg = set_date(yr  , 1, 1)
     Year_end = set_date(yr+1, 1, 1)

     weight = real( (Time - Year_beg) // (Year_end - Year_beg) , kind=FMS_TI_KIND_)

 end subroutine TIME_INTERP_FRAC_


 !> @brief Calculates fractional time between mid points of consecutive years
 subroutine TIME_INTERP_YEAR_ ( Time, weight, year1, year2 )

   type(time_type),   intent(in)  :: Time
   real(FMS_TI_KIND_), intent(out) :: weight !< fractional time between midpoints of year1 and year2
   integer        ,   intent(out) :: year1, year2

   integer :: yr, mo, dy, hour, minute, second
   type (time_type) :: Mid_year, Mid_year1, Mid_year2


   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, yr, mo, dy, hour, minute, second)

    ! mid point of current year
      Mid_year = year_midpt(yr)

      if ( Time >= Mid_year ) then
    ! current time is after mid point of current year
           year1  = yr
           year2  = yr+1
           Mid_year2 = year_midpt(year2)
           weight = real( (Time - Mid_year) // (Mid_year2 - Mid_year) , kind=FMS_TI_KIND_ )
      else
    ! current time is before mid point of current year
           year2  = yr
           year1  = yr-1
           Mid_year1 = year_midpt(year1)
           weight = real( (Time - Mid_year1) // (Mid_year - Mid_year1), kind=FMS_TI_KIND_ )
      endif

 end subroutine TIME_INTERP_YEAR_

 !> @brief Calculates fractional time between mid points of consecutive months
 subroutine TIME_INTERP_MONTH_ ( Time, weight, year1, year2, month1, month2 )

   type(time_type), intent(in)  :: Time
   real(FMS_TI_KIND_)           , intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2

   integer :: yr, mo, dy, hour, minute, second,  &
              mid_month, cur_month, mid1, mid2

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, yr, mo, dy, hour, minute, second)

    ! mid point of current month in seconds
      mid_month = days_in_month(Time) * halfday
    ! time into current month in seconds
      cur_month = second + secmin*minute + sechour*hour + secday*(dy-1)

      if ( cur_month >= mid_month ) then
    ! current time is after mid point of current month
           year1  = yr;  month1 = mo
           year2  = yr;  month2 = mo+1
           if (month2 > monyear)  then
              year2 = year2+1;  month2 = 1
           endif
           mid1 = mid_month
           mid2 = days_in_month(set_date(year2,month2,2)) * halfday
           weight = real(cur_month - mid1, FMS_TI_KIND_) / real(mid1+mid2, FMS_TI_KIND_)
      else
    ! current time is before mid point of current month
           year2  = yr;  month2 = mo
           year1  = yr;  month1 = mo-1
           if (month1 < 1)  then
              year1 = year1-1;  month1 = monyear
           endif
           if (year1>0) then
              mid1 = days_in_month(set_date(year1,month1,2)) * halfday
           else
              ! this can happen if we are at the beginning of year 1. In this case
              ! use December 0001 to calculate the duration of December 0000.
              ! This should work for all calendars
              mid1 = days_in_month(set_date(1,month1,2)) * halfday
           endif
           mid2 = mid_month
           weight = real(cur_month + mid1, FMS_TI_KIND_) / real(mid1+mid2, FMS_TI_KIND_)
      endif

 end subroutine TIME_INTERP_MONTH_

 !> @brief Calculates fractional time between mid points of consecutive days
 subroutine TIME_INTERP_DAY_ ( Time, weight, year1, year2, month1, month2, day1, day2 )

   type(time_type), intent(in)  :: Time
   real(FMS_TI_KIND_), intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2, day1, day2

   integer :: yr, mo, dy, hour, minute, second, sday

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, yr, mo, dy, hour, minute, second)

    ! time into current day in seconds
      sday = second + secmin*minute + sechour*hour

      if ( sday >= halfday ) then
    ! current time is after mid point of day
           year1 = yr;  month1 = mo;  day1 = dy
           year2 = yr;  month2 = mo;  day2 = dy + 1
           weight  = real(sday - halfday, FMS_TI_KIND_) / real(secday, FMS_TI_KIND_)

           if (day2 > days_in_month(Time)) then
               month2 = month2 + 1
               day2 = 1
               if (month2 > monyear) then
                    month2 = 1;  year2 = year2+1
               endif
           endif
      else
    ! current time is before mid point of day
           year2 = yr;  month2 = mo   ;  day2 = dy
           year1 = yr;  month1 = mo;  day1 = dy - 1
           weight  = real(sday + halfday,FMS_TI_KIND_) / real(secday,FMS_TI_KIND_)

           if (day1 < 1) then
               month1 = month1 - 1
               if (month1 < 1) then
                   month1 = monyear;  year1 = year1-1
               endif
               day1 = days_in_month(set_date(year1,month1,2))
           endif
      endif

 end subroutine TIME_INTERP_DAY_

 !> Part of the time_interp interface, calculates for cyclical data
 !! Time_beg and Time_end mark a repeating period
 !!
 !! Finds mid points and fractional weight for a time perioid
subroutine TIME_INTERP_MODULO_(Time, Time_beg, Time_end, Timelist, weight, index1, index2, &
                              correct_leap_year_inconsistency, err_msg)
type(time_type), intent(in)  :: Time !< a specific time value
type(time_type), intent(in)  :: Time_beg !< begining of period to search with
type(time_type), intent(in)  :: Time_end !< end of period to search with
type(time_type), intent(in)  :: Timelist(:) !< ascending time values to search between
real(FMS_TI_KIND_)           , intent(out) :: weight
integer        , intent(out) :: index1, index2 !< indices of bounding time values within Timelist
logical, intent(in), optional :: correct_leap_year_inconsistency!< When true turns on a kluge for an
                                !! inconsistency which may occur in a special case.
                                !! When the modulo time period (i.e. Time_end - Time_beg) is a
                                !! whole number of years and is not a multiple of 4, and the calendar
                                !! in use has leap years, then it is likely that the interpolation
                                !! will involve mapping a common year onto a leap year. In this case
                                !! it is often desirable, but not absolutely necessary, to use data
                                !! for Feb 28 of the leap year when it is mapped onto a common year.
character(len=*), intent(out), optional :: err_msg

  type(time_type) :: Period, T
  integer :: is, ie,i1,i2
  integer :: ys,ms,ds,hs,mins,ss ! components of the starting date
  integer :: ye,me,de,he,mine,se ! components of the ending date
  integer :: yt,mt,dt,ht,mint,st ! components of the current date
  integer :: dt1                 ! temporary value for day
  integer :: n                   ! size of Timelist
  integer :: stdoutunit
  logical :: correct_lyr, calendar_has_leap_years, do_the_lyr_correction
  integer, parameter :: kindl = FMS_TI_KIND_

  if ( .not. module_is_initialized ) call time_interp_init
  if( present(err_msg) ) err_msg = ''

  stdoutunit = stdout()
  n = size(Timelist)

  if (Time_beg>=Time_end) then
     if(fms_error_handler('time_interp_modulo', &
     'end of the specified time loop interval must be later than its beginning',err_msg)) return
  endif

  calendar_has_leap_years = (get_calendar_type() == JULIAN .or. get_calendar_type() == GREGORIAN)

  Period = Time_end-Time_beg ! period of the time axis

  if(present(correct_leap_year_inconsistency)) then
    correct_lyr = correct_leap_year_inconsistency
  else
    correct_lyr = .false.
  endif

  ! bring the requested time inside the specified time period
  T = Time

  do_the_lyr_correction = .false.

  ! Determine if the leap year correction needs to be done.
  ! It never needs to be done unless 3 conditions are met:
  ! 1) We are using a calendar with leap years
  ! 2) optional argument correct_leap_year_inconsistency is present and equals .true.
  ! 3) The modulo time period is an integer number of years
  ! If all of these are true then set do_the_lyr_correction to .true.

  if(calendar_has_leap_years .and. correct_lyr) then
    call get_date(Time_beg,ys,ms,ds,hs,mins,ss)
    call get_date(Time_end,ye,me,de,he,mine,se)
    if(ms==me.and.ds==de.and.hs==he.and.mins==mine.and.ss==se) then
      ! whole number of years
      do_the_lyr_correction = .true.
    endif
  endif

  if(do_the_lyr_correction) then
     call get_date(T,yt,mt,dt,ht,mint,st)
     yt = ys+modulo(yt-ys,ye-ys)
     dt1 = dt
     ! If it is Feb 29, but we map into a common year, use Feb 28
     if(mt==2.and.dt==29.and..not.leap_year(set_date(yt,1,1))) dt1=28
     T = set_date(yt,mt,dt1,ht,mint,st)
     if (T < Time_beg) then
       ! the requested time is within the first year,
       ! but before the starting date. So we shift it to the last year.
       if(mt==2.and.dt==29.and..not.leap_year(set_date(ye,1,1))) dt=28
       T = set_date(ye,mt,dt,ht,mint,st)
     endif
  else
     do while ( T >= Time_end )
        T = T-Period
     enddo
     do while ( T < Time_beg )
        T = T+Period
     enddo
  endif

  ! find indices of the first and last records in the Timelist that are within
  ! the requested time period.
  if (Time_end<=Timelist(1).or.Time_beg>=Timelist(n)) then
     if(get_calendar_type() == NO_CALENDAR) then
       call print_time(Time_beg,    'Time_beg'    )
       call print_time(Time_end,    'Time_end'    )
       call print_time(Timelist(1), 'Timelist(1)' )
       call print_time(Timelist(n), 'Timelist(n)' )
     else
       call print_date(Time_beg,    'Time_beg'    )
       call print_date(Time_end,    'Time_end'    )
       call print_date(Timelist(1), 'Timelist(1)' )
       call print_date(Timelist(n), 'Timelist(n)' )
     endif
     write(stdoutunit,*)'where n = size(Timelist) =',n
     if(fms_error_handler('time_interp_modulo', &
     'the entire time list is outside the specified time loop interval',err_msg)) return
  endif

  call bisect(Timelist,Time_beg,index1=i1,index2=i2)
  if (i1 < 1) then
     is = 1 ! Time_beg before lower boundary
  else if (Time_beg == Timelist(i1)) then
     is = i1 ! Time_beg right on the lower boundary
  else
     is = i2 ! Time_beg inside the interval or on upper boundary
  endif
  call bisect(Timelist,Time_end,index1=i1,index2=i2)
  if (Time_end > Timelist(i1)) then
    ie = i1
  else if (Time_end == Timelist(i1)) then
    if(Time_beg == Timelist(is)) then
      ! Timelist includes time levels at both the lower and upper ends of the period.
      ! The endpoints of Timelist specify the same point in the cycle.
      ! This ambiguity is resolved by ignoring the last time level.
      ie = i1-1
    else
      ie = i1
    endif
  else
!   This should never happen because bisect does not return i1 such that Time_end < Timelist(i1)
  endif
  if (is>=ie) then
     if(get_calendar_type() == NO_CALENDAR) then
       call print_time(Time_beg,    'Time_beg   =')
       call print_time(Time_end,    'Time_end   =')
       call print_time(Timelist(1), 'Timelist(1)=')
       call print_time(Timelist(n), 'Timelist(n)=')
     else
       call print_date(Time_beg,    'Time_beg   =')
       call print_date(Time_end,    'Time_end   =')
       call print_date(Timelist(1), 'Timelist(1)=')
       call print_date(Timelist(n), 'Timelist(n)=')
     endif
     write(stdoutunit,*)'where n = size(Timelist) =',n
     write(stdoutunit,*)'is =',is,'ie =',ie
     if(fms_error_handler('time_interp_modulo', &
     'error in calculation of time list bounds within the specified time loop interval',err_msg)) return
  endif

  ! handle special cases:
  if( T>=Timelist(ie) ) then
     ! time is after the end of the portion of the time list within the requested period
     index1 = ie;   index2 = is
     weight = real((T-Timelist(ie))//(Period-(Timelist(ie)-Timelist(is))), FMS_TI_KIND_ )
  else if (T<Timelist(is)) then
     ! time is before the beginning of the portion of the time list within the requested period
     index1 = ie;   index2 = is
     weight = real( 1.0_kindl -((Timelist(is)-T)//(Period-(Timelist(ie)-Timelist(is)))), FMS_TI_KIND_ )
  else
     call bisect(Timelist,T,index1,index2)
     weight = real((T-Timelist(index1)) // (Timelist(index2)-Timelist(index1)), FMS_TI_KIND_ )
  endif

end subroutine TIME_INTERP_MODULO_

subroutine TIME_INTERP_LIST_ ( Time, Timelist, weight, index1, index2, modtime, err_msg )
type(time_type)  , intent(in)  :: Time, Timelist(:)
real(FMS_TI_KIND_)             , intent(out) :: weight
integer          , intent(out) :: index1, index2
integer, optional, intent(in)  :: modtime
character(len=*), intent(out), optional :: err_msg

integer :: n, hr, mn, se, mtime
type(time_type) :: T, Ts, Te, Td, Period, Time_mod
character(len=:),allocatable :: terr, tserr, teerr
integer, parameter :: kindl = FMS_TI_KIND_

  if ( .not. module_is_initialized ) call time_interp_init

  if( present(err_msg) ) err_msg = ''

  weight = 0.0_kindl; index1 = 0; index2 = 0
  n = size(Timelist(:))

! setup modular time axis?
  mtime = NONE
  if (present(modtime)) then
     mtime = modtime
     Time_mod = (Timelist(1)+Timelist(n))/2
     call get_date (Time_mod, yrmod, momod, dymod, hr, mn, se)
     mod_leapyear = leap_year(Time_mod)
  endif

! set period for modulo axis
  select case (mtime)
     case (NONE)
       ! do nothing
     case (YEAR)
         Period = set_time(0,days_in_year(Time_mod))
     case (MONTH)
       ! month length must be equal
         if (days_in_month(Time_mod) /= days_in_month(Time)) then
            if(fms_error_handler ('time_interp_list','modulo months must have same length',err_msg)) return
         endif
         Period = set_time(0,days_in_month(Time_mod))
     case (DAY)
         Period = set_time(0,1)
     case default
         if(fms_error_handler ('time_interp_list','invalid value for argument modtime',err_msg)) return
  end select

! If modulo time is in effect and Timelist spans a time interval exactly equal to
! the modulo period, then the endpoints of Timelist specify the same point in the cycle.
! This ambiguity is resolved by ignoring the last time level.
  if (mtime /= NONE .and. Timelist(size(Timelist))-Timelist(1) == Period) then
     n = size(Timelist) - 1
  else
     n = size(Timelist)
  endif

! starting and ending times from list
  Ts = Timelist(1)
  Te = Timelist(n)
  Td = Te-Ts
  T  = set_modtime(Time,mtime)

! Check that Timelist does not span a time interval greater than the modulo period
  if (mtime /= NONE) then
     if (Td > Period) then
        if(fms_error_handler ('time_interp_list','period of list exceeds modulo period',err_msg)) return
     endif
  endif

! time falls on start or between start and end list values
  if ( T >= Ts .and. T < Te ) then
     call bisect(Timelist(1:n),T,index1,index2)
     weight = real( (T-Timelist(index1)) // (Timelist(index2)-Timelist(index1)), FMS_TI_KIND_)

! time falls before starting list value
  else if ( T < Ts ) then
     if (mtime == NONE) then
        call time_list_error(T,terr)
        call time_list_error(Ts,tserr)
        call time_list_error(Te,teerr)
        if(fms_error_handler ('time_interp_list',&
           'time '//trim(terr)//' ('//date_to_string(T)//' is before range of list '//trim(tserr)//'-'//trim(teerr)//&
           '('//date_to_string(Ts)//' - '//date_to_string(Te)//')',&
           err_msg)) return
        deallocate(terr,tserr,teerr)
     endif
     Td = Te-Ts
     weight = real(1.0_kindl - ((Ts-T) // (Period-Td)), FMS_TI_KIND_ )
     index1 = n
     index2 = 1

! time falls on ending list value
  else if ( T == Te ) then
    if(perthlike_behavior) then
       weight = 1.0_kindl
       index1 = n-1
       index2 = n
    else
       weight = 0.0_kindl
       index1 = n
       if (mtime == NONE) then
         index2 = n
       else
         index2 = 1
       endif
    endif

! time falls after ending list value
  else if ( T > Te ) then
     if (mtime == NONE) then
        call time_list_error(T,terr)
        call time_list_error(Ts,tserr)
        call time_list_error(Te,teerr)
        if(fms_error_handler ('time_interp_list',&
           'time '//trim(terr)//' ('//date_to_string(T)//' is after range of list '//trim(tserr)//'-'//trim(teerr)//&
           '('//date_to_string(Ts)//' - '//date_to_string(Te)//')',&
           err_msg)) return
        deallocate(terr,tserr,teerr)
     endif
     Td = Te-Ts
     weight = real( (T-Te) // (Period-Td), FMS_TI_KIND_)
     index1 = n
     index2 = 1
  endif

end subroutine TIME_INTERP_LIST_
!> }