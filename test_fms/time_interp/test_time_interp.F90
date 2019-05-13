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

program test_time_interp

 use          fms_mod, only: fms_init, fms_end, stdout, stdlog, FATAL, mpp_error
 use time_manager_mod, only: get_date, set_time, set_date, time_manager_init, set_calendar_type, operator(+)
 use time_manager_mod, only: JULIAN, time_type, increment_time, NOLEAP, print_date
 use  time_interp_mod, only: time_interp_init, time_interp, NONE, YEAR, MONTH, DAY

 implicit none

 integer, parameter :: num_Time=6
 type(time_type) :: Time_beg, Time_end, Time(num_Time)
 type(time_type), allocatable, dimension(:) :: Timelist
 integer :: index1, index2, mo, yr, timelist_len, outunit, ntest, nline
 real :: weight

 integer :: nmin, nmax

 namelist / test_time_interp_nml / timelist_len

 call fms_init
 outunit = stdout()
 call set_calendar_type(JULIAN)
 call time_interp_init

 Time_beg = set_date(1, 1, 1)
 Time_end = set_date(2, 1, 1)
 Time(1) = Time_beg
 Time(2) = set_date(1, 1,16)
 Time(3) = set_date(1, 2, 1)
 Time(4) = set_date(1,12, 1)
 Time(5) = set_date(1,12,16)
 Time(6) = Time_end

! Tests with modulo time
 do nline=1,3
   if(nline == 1) then
     allocate(Timelist(12))
     do mo=1,12
       Timelist(mo) = set_date(1, mo, 1)
     enddo
   else if(nline == 2) then
     allocate(Timelist(13))
     do mo=1,12
       Timelist(mo) = set_date(1, mo, 1)
     enddo
     Timelist(13) = set_date(2, 1, 1)
   else if(nline == 3) then
     allocate(Timelist(12))
     do mo=2,12
       Timelist(mo-1) = set_date(1, mo, 1)
     enddo
     Timelist(12) = set_date(2, 1, 1)
   endif

   do ntest=1,num_Time
     call diagram(nline,ntest,modulo_time=.true.)
     call time_interp(Time(ntest), Time_beg, Time_end, Timelist, weight, index1, index2)
     write(outunit,*) 'time_interp_modulo:'
     write(outunit,'()')
     call print_date(Time(ntest),                'Time       =')
     call print_date(Time_beg,                   'Time_beg   =')
     call print_date(Time_end,                   'Time_end   =')
     call print_date(Timelist(1),                'Timelist(1)=')
     call print_date(Timelist(size(Timelist(:))),'Timelist(n)=')
     write(outunit,99) index1,index2,weight
     write(outunit,'()')

     call time_interp(Time(ntest), Timelist, weight, index1, index2, modtime=YEAR)
     write(outunit,*) 'time_interp_list with modtime=YEAR:'
     write(outunit,'()')
     call print_date(Time(ntest),                'Time       =')
     call print_date(Timelist(1),                'Timelist(1)=')
     call print_date(Timelist(size(Timelist(:))),'Timelist(n)=')
     write(outunit,99) index1,index2,weight
   enddo
   deallocate(Timelist)
 enddo

! Tests without modulo time
 do nline=1,3
   if(nline == 1) then
     allocate(Timelist(12))
     do mo=1,12
       Timelist(mo) = set_date(1, mo, 1)
     enddo
   else if(nline == 2) then
     allocate(Timelist(13))
     do mo=1,12
       Timelist(mo) = set_date(1, mo, 1)
     enddo
     Timelist(13) = set_date(2, 1, 1)
   else if(nline == 3) then
     allocate(Timelist(12))
     do mo=2,12
       Timelist(mo-1) = set_date(1, mo, 1)
     enddo
     Timelist(12) = set_date(2, 1, 1)
   endif

   if(nline == 1) then
     nmin = 1; nmax = 4
   else if(nline == 2) then
     nmin = 1; nmax = num_Time
   else if(nline == 3) then
     nmin = 3; nmax = num_Time
   endif
   do ntest=nmin,nmax
     call diagram(nline,ntest,modulo_time=.false.)
     call time_interp(Time(ntest), Timelist, weight, index1, index2, modtime=NONE)
     write(outunit,*) 'time_interp_list with modtime=NONE:'
     write(outunit,'()')
     call print_date(Time(ntest),                'Time       =')
     call print_date(Timelist(1),                'Timelist(1)=')
     call print_date(Timelist(size(Timelist(:))),'Timelist(n)=')
     write(outunit,99) index1,index2,weight
   enddo
   deallocate(Timelist)
 enddo

! More tests with modulo time
 Time_beg = set_date(1999, 1, 1)
 Time_end = set_date(2000, 1, 1)
 Time(1)  = set_date(1998, 1, 1)
 Time(2)  = set_date(1998, 2,28)
 Time(3)  = set_date(1998,12,16)
 Time(4)  = set_date(2000, 1, 1)
 Time(5)  = set_date(2000, 2,28)
 Time(6)  = set_date(2000, 2,29)

 allocate(Timelist(13))
 do mo=1,12
   Timelist(mo) = set_date(1999, mo, 1)
 enddo
 Timelist(13) = set_date(2000, 1, 1)

 write(outunit,'("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>",/)')
 write(outunit,'()')
 write(outunit,*) 'time_interp_modulo with correct_leap_year_inconsistency=.true.'
 write(outunit,'()')
 write(outunit,'(" Jan 1 1999                                     Jan 1 2000")')
 write(outunit,'("    |                                               |")')
 write(outunit,'("    v                                               v")')
 write(outunit,'("    x---x---x---x---x---x---x---x---x---x---x---x---x")')
 write(outunit,'("    ^                                               ^")')
 write(outunit,'("    |                                               |")')
 write(outunit,'(" Time_beg                                        Time_end ")')
 write(outunit,'()')

 do ntest=1,num_Time
   call time_interp(Time(ntest), Time_beg, Time_end, Timelist, weight, index1, index2, correct_leap_year_inconsistency=.true.)
   call print_date(Time(ntest),' Time =')
   write(outunit,99) index1,index2,weight
   write(outunit,'()')
 enddo
 deallocate(Timelist)

! Tests of modulo time and leap year inconsistency
 Time_beg = set_date(1978, 1, 1)
 Time_end = set_date(1981, 1, 1)
 Time(1)  = set_date(1976, 2,28)
 Time(2)  = set_date(1976, 2,29)
 Time(3)  = set_date(1976, 3, 1)
 Time(4)  = set_date(1983, 2,28)
 Time(5)  = set_date(1983, 3, 1)
 Time(6)  = set_date(1981, 1, 1)
 allocate(Timelist(37))
 do yr=1978,1980
   do mo=1,12
     Timelist(12*(yr-1978)+mo) = set_date(yr, mo, 1)
   enddo
 enddo
 Timelist(37) = set_date(1981, 1, 1)

 write(outunit,'("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")')
 write(outunit,'()')
 write(outunit,*) 'time_interp_modulo with correct_leap_year_inconsistency=.true.'
 write(outunit,'()')
 write(outunit,'(" Jan 1 1978              Jan 1 1979              Jan 1 1980              Jan 1 1981")')
 write(outunit,'("     |                       |                       | <---- leap year ----> |")')
 write(outunit,'("     v                       v                       v                       v")')
 write(outunit,'("     x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x")')
 write(outunit,'("     ^                                                                       ^")')
 write(outunit,'("     |                                                                       |")')
 write(outunit,'("  Time_beg                                                               Time_end")')
 write(outunit,'()')

 do ntest=1,num_Time
   call time_interp(Time(ntest), Time_beg, Time_end, Timelist, weight, index1, index2, correct_leap_year_inconsistency=.true.)
   call print_date(Time(ntest),' Time=')
   write(outunit,99) index1,index2,weight
   write(outunit,'()')
 enddo
 deallocate(Timelist)

 allocate(Timelist(12))
 Timelist( 1) = set_date(1,  1, 16, hour=12) ! Jan midmonth
 Timelist( 2) = set_date(1,  2, 15, hour= 0) ! Feb midmonth (common year)
 Timelist( 3) = set_date(1,  3, 16, hour=12) ! Mar midmonth
 Timelist( 4) = set_date(1,  4, 16, hour= 0) ! Apr midmonth
 Timelist( 5) = set_date(1,  5, 16, hour=12) ! May midmonth
 Timelist( 6) = set_date(1,  6, 16, hour= 0) ! Jun midmonth
 Timelist( 7) = set_date(1,  7, 16, hour=12) ! Jul midmonth
 Timelist( 8) = set_date(1,  8, 16, hour=12) ! Aug midmonth
 Timelist( 9) = set_date(1,  9, 16, hour= 0) ! Sep midmonth
 Timelist(10) = set_date(1, 10, 16, hour=12) ! Oct midmonth
 Timelist(11) = set_date(1, 11, 16, hour= 0) ! Nov midmonth
 Timelist(12) = set_date(1, 12, 16, hour=12) ! Dec midmonth
 Time_beg = set_date(1, 1, 1)
 Time_end = set_date(2, 1, 1)
 call diagram(nline=4, ntest=0, modulo_time=.true.)
 do ntest=0,73
   Time(1) = set_date(1996, 1, 1) + set_time(seconds=0, days=5*ntest)
   call print_date(Time(1),' Time=')
   call time_interp(Time(1), Timelist, weight, index1, index2, modtime=YEAR)
   write(outunit,89) 'time_interp_list with modtime=YEAR:   ', index1,index2,weight
   call time_interp(Time(1), Time_beg, Time_end, Timelist, weight, index1, index2, correct_leap_year_inconsistency=.true.)
   write(outunit,89) 'time_interp_modulo: ', index1,index2,weight
   write(outunit,'()')
 enddo

 99 format(' index1=',i3,'  index2=',i3,'  weight=',f18.15)
 89 format(a20,' index1=',i3,'  index2=',i3,'  weight=',f18.15)
 call fms_end

 contains

 subroutine diagram(nline,ntest,modulo_time)
 integer, intent(in) :: nline,ntest
 logical, intent(in) :: modulo_time

 write(outunit,'("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")')
 write(outunit,'()')
 if(modulo_time) then
   write(outunit,'(" Time_beg                                      Time_end")')
   write(outunit,'("  |                                               |")')
   write(outunit,'("  v                                               v")')
 endif

 if(nline == 1) then
   write(outunit,'("  x---x---x---x---x---x---x---x---x---x---x---x----")')
 else if(nline == 2) then
   write(outunit,'("  x---x---x---x---x---x---x---x---x---x---x---x---x")')
 else if(nline == 3) then
   write(outunit,'("  ----x---x---x---x---x---x---x---x---x---x---x---x")')
 else if(nline == 4) then
   write(outunit,'("  --x---x---x---x---x---x---x---x---x---x---x---x--")')
 endif

 if(ntest == 1) then
   write(outunit,'("  ^")  ')
   write(outunit,'("  |")  ')
   write(outunit,'(" Time")')
 else if(ntest == 2) then
   write(outunit,'("    ^")  ')
   write(outunit,'("    |")  ')
   write(outunit,'("   Time")')
 else if(ntest == 3) then
   write(outunit,'("      ^")  ')
   write(outunit,'("      |")  ')
   write(outunit,'("     Time")')
 else if(ntest == 4) then
   write(outunit,'("                                              ^")  ')
   write(outunit,'("                                              |")  ')
   write(outunit,'("                                             Time")')
 else if(ntest == 5) then
   write(outunit,'("                                                ^")  ')
   write(outunit,'("                                                |")  ')
   write(outunit,'("                                               Time")')
 else if(ntest == 6) then
   write(outunit,'("                                                  ^")  ')
   write(outunit,'("                                                  |")  ')
   write(outunit,'("                                                 Time")')
 endif
 write(outunit,'()')

 end subroutine diagram

 end program test_time_interp
