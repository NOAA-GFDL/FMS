
module time_interp_mod

!-----------------------------------------------------------------------

use time_manager_mod, only: time_type, get_date, set_date,  &
                            days_in_month, leap_year,       &
                            operator(+), operator(-), operator(//),  &
                            operator(/), operator(>=)

implicit none
private

public  time_interp, fraction_of_year
!-----------------------------------------------------------------------

interface time_interp
    module procedure time_interp_year, time_interp_month,  &
                     time_interp_day
end interface

!-----------------------------------------------------------------------

   integer, parameter ::  secmin = 60, minhour = 60, hourday = 24,  &
                         sechour = secmin*minhour,                  &
                          secday = secmin*minhour*hourday

   integer, parameter :: monyear = 12
   integer, parameter :: halfday = secday/2

contains

!#######################################################################

subroutine time_interp_year (time, frac, year1, year2)

   type (time_type), intent(in)  :: time
            real,    intent(out) :: frac
            integer, intent(out) :: year1, year2

   integer :: year, month, day, hour, minute, second
   type (time_type) :: mid_year, mid_year1, mid_year2

      call get_date (time, year, month, day, hour, minute, second)

      mid_year = year_midpt(year)

      if ( time >= mid_year ) then
           year1  = year
           year2  = year+1

           mid_year2 = year_midpt(year2)
           frac = (time - mid_year) // (mid_year2 - mid_year)
      else
           year2  = year
           year1  = year-1

           mid_year1 = year_midpt(year1)
           frac = (time - mid_year1) // (mid_year - mid_year1)
      endif

end subroutine time_interp_year

!#######################################################################

subroutine time_interp_month (time, frac, year1, year2, month1, month2)

   type (time_type), intent(in)  :: time
            real,    intent(out) :: frac
            integer, intent(out) :: year1, year2, month1, month2

   integer :: year, month, day, hour, minute, second,  &
              mid_month, cur_month, mid1, mid2


      call get_date (time, year, month, day, hour, minute, second)

      mid_month = days_in_month(time) * halfday
      cur_month = second + secmin*minute + sechour*hour + secday*(day-1)

      if ( cur_month >= mid_month ) then
           year1  = year;  month1 = month
           year2  = year;  month2 = month+1
           if (month2 > monyear)  year2 = year2+1
           if (month2 > monyear) month2 = 1
           mid1 = mid_month
           mid2 = days_in_month(set_date(year2,month2,2)) * halfday
           frac = float(cur_month - mid1) / float(mid1+mid2)
      else
           year2  = year;  month2 = month
           year1  = year;  month1 = month-1
           if (month1 < 1)  year1 = year1-1
           if (month1 < 1) month1 = monyear
           mid1 = days_in_month(set_date(year1,month1,2)) * halfday
           mid2 = mid_month
           frac = float(cur_month + mid1) / float(mid1+mid2)
      endif

end subroutine time_interp_month

!#######################################################################

subroutine time_interp_day (time, frac, year1, year2, month1, month2,  &
                                        day1,  day2)

   type (time_type), intent(in)  :: time
            real,    intent(out) :: frac
            integer, intent(out) :: year1, year2, month1, month2,  &
                                     day1,  day2

   integer :: year, month, day, hour, minute, second, sday


      call get_date (time, year, month, day, hour, minute, second)

      sday = second + secmin*minute + sechour*hour

      if ( sday >= halfday ) then
           year1 = year;  month1 = month;  day1 = day
           year2 = year;  month2 = month;  day2 = day + 1
           frac  = float(sday - halfday) / float(secday)

           if (day2 > days_in_month(time)) then
               month2 = month2 + 1
               day2 = 1
               if (month2 > monyear) then
                    month2 = 1;  year2 = year2+1
               endif
           endif
      else
           year2 = year;  month2 = month;  day2 = day
           year1 = year;  month1 = month;  day1 = day - 1
           frac  = float(sday + halfday) / float(secday)

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

function fraction_of_year (time)

   type (time_type), intent(in) :: time
   real                         :: fraction_of_year

   integer          :: year, month, day, hour, minute, second
   type (time_type) :: year_beg, year_end


!  ---- compute fractional time of year -----

    call get_date (time, year, month, day, hour, minute, second)

    year_beg = set_date(year  , 1, 1)
    year_end = set_date(year+1, 1, 1)

    fraction_of_year = (time - year_beg) // (year_end - year_beg)


end function fraction_of_year

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

end module time_interp_mod

