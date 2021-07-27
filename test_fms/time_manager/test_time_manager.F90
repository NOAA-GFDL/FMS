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

program test_time_manager

 use          mpp_mod, only: input_nml_file, mpp_error, NOTE, FATAL
 use          fms_mod, only: fms_init, fms_end, stderr
 use          fms_mod, only: check_nml_error
 use    constants_mod, only: constants_init, rseconds_per_day=>seconds_per_day
 use time_manager_mod, only: time_type, set_date, get_date, set_time, set_calendar_type, real_to_time_type
 use time_manager_mod, only: length_of_year, leap_year, days_in_month, days_in_year, print_time
 use time_manager_mod, only: set_ticks_per_second, get_ticks_per_second
 use time_manager_mod, only: decrement_date, increment_date, get_time, increment_time, decrement_time
 use time_manager_mod, only: JULIAN, GREGORIAN, THIRTY_DAY_MONTHS, NOLEAP
 use time_manager_mod, only: operator(-), operator(+),  operator(*),  operator(/),  &
                             operator(>), operator(>=), operator(==), operator(/=), &
                             operator(<), operator(<=), operator(//), assignment(=)

 implicit none

 type(time_type) :: Time, Time0, time1, time2
 real    :: xx
 integer :: yr, mo, day, hr, min, sec, ticks
 integer :: yr0, mo0, day0, hr0, min0, sec0, ticks0
 integer :: year, month, dday, days_this_month
 integer :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
 logical :: leap
 integer :: nr, icode, ierr, io, nn, errunit, outunit
 character(len=256) :: err_msg, char_date
 character(len=8),  allocatable, dimension(:) :: test_time
 character(len=23), allocatable, dimension(:) :: test_date
 character(len=8) :: test_name
 character(len=256) :: out_msg

 !: for testing set/get_date_gregorian
 integer, parameter :: days_in_400_year_period = 146097
 integer, dimension(days_in_400_year_period) :: coded_date
 integer, dimension(400,12,31) :: date_to_day

 logical :: test1 =.true.,test2 =.true.,test3 =.true.,test4 =.true.,test5 =.true.,test6 =.true.,test7 =.true.,test8 =.true.
 logical :: test9 =.true.,test10=.true.,test11=.true.,test12=.true.,test13=.true.,test14=.true.,test15=.true.,test16=.true.
 logical :: test17=.true.,test18=.true.,test19=.true.,test20=.true.

 namelist / test_nml / test1 ,test2 ,test3 ,test4 ,test5 ,test6 ,test7 ,test8,  &
                       test9 ,test10,test11,test12,test13,test14,test15,test16, &
                       test17,test18,test19,test20

 call fms_init
 call constants_init

 read (input_nml_file, test_nml, iostat=io)
 ierr = check_nml_error (io, 'test_nml')

 open(newunit = outunit, file='test_time_manager.out', status='replace', form='formatted')

 errunit = stderr()
 call set_ticks_per_second(10)

 !==============================================================================================
 ! Tests of set_time_i and get_time without ticks

 if(test1) then
   write(outunit,'(/,a)') '#################################  test1  #################################'
   Time = set_time(seconds=2, days=1)
   call get_time(Time, sec, day, ticks)
   write(outunit,'(a,i2,a,i8,a,i2)') ' test1.1: days=',day,' seconds=',sec,' ticks=',ticks
   call get_time(Time, sec, day)
   write(outunit,'(a,i2,a,i8)') ' test1.2: days=',day,' seconds=',sec
   call get_time(Time, sec)
   write(outunit,'(a,i8)') ' test1.2: seconds=',sec
 endif
 !==============================================================================================
 ! Tests of set_time_i and get_time with ticks

 if(test2) then
   write(outunit,'(/,a)') '#################################  test2  #################################'
   Time = set_time(seconds=2, days=1, ticks=5)
   call get_time(Time, sec, day, ticks)
   write(outunit,'(a,i2,a,i6,a,i2)') ' test2.1: days=',day,' seconds=',sec,' ticks=',ticks
   call get_time(Time, sec, ticks=ticks)
   write(outunit,'(a,i6,a,i2)') ' test2.2: seconds=',sec,' ticks=',ticks
   call get_time(Time, sec, day, err_msg=err_msg)
   if(err_msg /= '') then
     write(outunit,'(a)') ' test2.3 successful: '//trim(err_msg)
   else
     call mpp_error(FATAL, "ERROR: test2.3 fails, did not get expected error message")
   endif
   call get_time(Time, sec, err_msg=err_msg)
   if(err_msg /= '') then
     write(outunit,'(a)') ' test2.4 successful: '//trim(err_msg)
   else
     call mpp_error(FATAL, "ERROR: test2.4 fails, did not get expected error message")
   endif
 endif
 !==============================================================================================
 ! Tests of time operators
 ! Test of function scalar_time_mult is not necessary, it simply calls time_scalar_mult.
 ! Test of function time_ne is not necessary, it simply calls time_eq.
 ! Test of function time_ge is not necessary, it simply calls time_gt.
 ! Test of function time_le is not necessary, it simply calls time_lt and time_eq.
 ! Test of function time_ne is not necessary, it simply calls time_eq.

  if(test3) then
    write(outunit,'(/,a)') '#################################  test3  #################################'
 !  Test of function time_plus
    call print_time(set_time(seconds=0, days=2, ticks=5) + set_time(seconds=0, days=2, ticks=6), 'test3.1:', unit=outunit)

 !  Test of function time_minus
 !  The minus operator for time ensures a positive result. In effect is does this: abs(time1-time2)
    call print_time(set_time(seconds=0, days=2, ticks=5) - set_time(seconds=0, days=2, ticks=6), 'test3.2:', unit=outunit)

 !  Test of function time_scalar_mult.  Note that 25000*86399 is greater than huge = 2**31 - 1
    call print_time(2*set_time(seconds=0, days=2, ticks=6), 'test3.3:', unit=outunit)
    call print_time(25000*set_time(seconds=86399, days=0, ticks=0), 'test3.4:', unit=outunit)

 !  Test of function time_scalar_divide
    call print_time(set_time(seconds=0, days=60000, ticks=2)/2, 'test3.5:', unit=outunit)

 !  Test of function time_real_divide
    xx = set_time(seconds=0, days=60000, ticks=2)//set_time(seconds=86400)
    write(outunit,'("test3.6: xx=",f15.9)') xx

 !  Test of function time_divide
    nn = set_time(seconds=0, days=60000, ticks=2)//set_time(seconds=86400)
    write(outunit,'("test3.7: nn=",i6)') nn

 !  Test of function time_gt
    if(set_time(seconds=1, days=1, ticks=2) > set_time(seconds=1, days=1, ticks=1)) then
      write(outunit,'("test3.8 successful")')
    else
      call mpp_error(FATAL, "ERROR: test3.8 fails, did not get expected result")
    endif
    if(set_time(seconds=1, days=1, ticks=2) > set_time(seconds=1, days=1, ticks=2)) then
      call mpp_error(FATAL, "ERROR: test3.9 fails, did not get expected result")
    else
      write(outunit,'("test3.9 successful")')
    endif

 !  Test of function time_lt
    if(set_time(seconds=1, days=1, ticks=1) < set_time(seconds=1, days=1, ticks=2)) then
      write(outunit,'("test3.10 successful")')
    else
      call mpp_error(FATAL, "ERROR: test3.10 fails, did not get expected result")
    endif
    if(set_time(seconds=1, days=1, ticks=2) < set_time(seconds=1, days=1, ticks=2)) then
      call mpp_error(FATAL, "ERROR: test3.11 fails, did not get expected result")
    else
      write(outunit,'("test3.11 successful")')
    endif

 !  Test of function time_eq
    if(set_time(seconds=1, days=1, ticks=1) == set_time(seconds=1, days=1, ticks=1)) then
      write(outunit,'("test3.12 successful")')
    else
      call mpp_error(FATAL, "ERROR: test3.12 fails, did not get expected result")
    endif
    if(set_time(seconds=1, days=1, ticks=1) == set_time(seconds=1, days=1, ticks=2)) then
      call mpp_error(FATAL, "ERROR: test3.13 fails, did not get expected result")
    else
      write(outunit,'("test3.13 successful")')
    endif
  endif
 !==============================================================================================
 ! Tests of set_time_c

 if(test4) then
   write(outunit,'(/,a)') '#################################  test4  #################################'
   test_name = 'test4.  '
   allocate(test_time(15))
   test_time( 1: 6) = (/'1 10    ','1 10.   ','1 10.000','1  0.0  ','1   .000','1   .   '/)
   test_time( 7: 9) = (/'1 10.20 ','1 10.300','1  0.40 '/)
   test_time(10:15) = (/'1   .510','2 .50001','1.0 10.2','10.30000','10-0.40 ','10:1.510'/) ! invalid forms
   do nr=1,9
     write(test_name(7:8),'(i2.2)') nr
     Time = set_time(trim(test_time(nr)), err_msg=err_msg, allow_rounding=.false.)
     if(err_msg == '') then
       call print_time(Time, test_name//':', unit=outunit)
     else
       out_msg = test_name // ' fails: ' // trim(err_msg)
       call mpp_error(FATAL, out_msg)
     endif
   enddo

   test_time(1:6) = (/'1   .510','2 .50001','1.0 10.2','10.30000','10-0.40 ','10:1.510'/)
   do nr=10,15
     write(test_name(7:8),'(i2.2)') nr
     Time = set_time(trim(test_time(nr)), err_msg=err_msg, allow_rounding=.false.)
     if(err_msg /= '') then
       write(outunit,'(a)') test_name//' successful: '//trim(err_msg)
     else
       out_msg = test_name // ' fails: did not get expected error message'
       call mpp_error(FATAL, out_msg)
     endif
   enddo
 endif

 !==============================================================================================
 ! Tests of set_date_i

 if(test5) then
   write(outunit,'(/,a)') '#################################  test5  #################################'
   call set_calendar_type(JULIAN)
   call print_time(set_date(1980, 1, 1, 0, 0, 0),' test5.1:', unit=outunit)
   call print_time(set_date(1980, 1, 2, 3, 4, 5, 6),' test5.2:', unit=outunit)
   call print_time(set_date(1980, 1, 2, tick=6),' test5.3:', unit=outunit)
   Time = set_date(1980, 1, 2, tick=10, err_msg=err_msg)
   if(err_msg == '') then
     call mpp_error(FATAL, 'Test5.4 fails: did not get expected error')
   else
     write(outunit,'(a)') ' test5.4 successful: '//trim(err_msg)
   endif
 endif
 !==============================================================================================
 ! Tests of set_date_c

 if(test6) then
   write(outunit,'(/,a)') '#################################  test6  #################################'
   test_name = 'test6.  '
   call set_calendar_type(GREGORIAN)
   allocate(test_date(6))
   test_date(1:3) = (/' 1980-12-30 01:01:11   ',' 1980-12-30 01:01:11.50',' 1980-12-30 01:01:11.55'/)
   test_date(4:6) = (/' 1980-12-30 01:01:11.96','   1980-1-3 1:1:11     ','   1980-1-3 1:1:11.99  '/)
   do nr=1,6
     write(test_name(7:8),'(i2.2)') nr
     Time = set_date(trim(test_date(nr)), err_msg=err_msg, allow_rounding=.true., zero_year_warning=.true.)
     if(err_msg == '') then
       call print_time(Time,test_name//' successful:', unit=outunit)
     else
       out_msg = test_name // ' fails: ' // trim(err_msg)
       call mpp_error(FATAL, out_msg)
     endif
   enddo
   call set_calendar_type(THIRTY_DAY_MONTHS)
   call print_time(set_date('1900-02-30 00:00:00'),'test6.7:', unit=outunit)
   Time = set_date('1900-01-31 00:00:00', err_msg=err_msg)
   if(err_msg == '') then
     call mpp_error(FATAL, 'Test6.8 fails: did not get expected error')
   else
     write(outunit,'(a)') 'test6.8 successful '//trim(err_msg)
   endif
   call set_calendar_type(JULIAN)
   Time = set_date('1901-02-29 00:00:00', err_msg=err_msg)
   if(err_msg == '') then
     call mpp_error(FATAL, 'Test6.9 fails: did not get expected error')
   else
     write(outunit,'(a)') 'test6.9 successful '//trim(err_msg)
   endif
 endif
!==============================================================================================
! Tests of decrement_date and increment_date

 if(test7) then
   write(outunit,'(/,a)') '#################################  test7  #################################'
   char_date = '1904-01-01 00:00:00'
   write(outunit,'(a)') ' Initial date='//trim(char_date)//':00'

   do nr=1,4
     write(outunit,'("=================================================================")')
     if(nr == 1) then
       call set_calendar_type(THIRTY_DAY_MONTHS)
       write(outunit,'(" THIRTY_DAY_MONTHS")')
     endif
     if(nr == 2) then
       call set_calendar_type(NOLEAP)
       write(outunit,'(" NOLEAP")')
     endif
     if(nr == 3) then
       call set_calendar_type(JULIAN)
       write(outunit,'(" JULIAN")')
     endif
     if(nr == 4) then
       call set_calendar_type(GREGORIAN)
       write(outunit,'(" GREGORIAN")')
     endif
     time1 = set_date(trim(char_date))
     do year=-1,1
       do month=-1,1
         write(outunit,'(" test of decrement_date increments: year=",i2," month=",i2)') year,month
         time2 = decrement_date(time1, year, month, err_msg=err_msg)
         if(err_msg /= '') then
           write(outunit,'(a)') 'test of decrement_date fails '//trim(err_msg)
         else
           call get_date(time2, yr, mo, day, hr, min, sec, ticks)
           write(outunit,20) yr, mo, day, hr, min, sec, ticks
         endif
       enddo
     enddo
     time1 = set_date(1, 1, 2, 1, 1, 1, 1, err_msg)
     write(outunit,'(" Initial date = 01-01-02 01:01:01:01")')
     do icode=0,242
       day   = modulo(icode/81,3) - 1
       hr    = modulo(icode/27,3) - 1
       min   = modulo(icode/9, 3) - 1
       sec   = modulo(icode/3, 3) - 1
       ticks = modulo(icode   ,3) - 1
       write(outunit,11) day, hr, min, sec, ticks
       time2 = increment_date(time1, 0, 0, day, hr, min, sec, ticks, err_msg)
       call get_date(time2, yr, mo, day, hr, min, sec, ticks)
       write(outunit,20) yr, mo, day, hr, min, sec, ticks
     enddo
   enddo
 endif

  11 format(' test of increment_date increments: day=',i2,' hr=',i2,' min=',i2,' sec=',i2,' ticks=',i2)
  20 format(' time=',i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2, ':', i2.2)
 !==============================================================================================
 ! Tests involving Feb 29

  if(test8) then
    write(outunit,'(/,a)') '#################################  test8  #################################'
    call set_calendar_type(THIRTY_DAY_MONTHS)
    Time = set_date('1904-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
      call print_time(Time, 'test8.1 successful', unit=outunit)
    else
       out_msg = 'test8.1 fails: ' // trim(err_msg)
       call mpp_error(FATAL, out_msg)
    endif

    call set_calendar_type(NOLEAP)
    Time = set_date('1904-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test8.2 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test8.2 successful: '//trim(err_msg)
    endif

    call set_calendar_type(GREGORIAN)
    Time = set_date('1900-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test8.3 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test8.3 successful: '//trim(err_msg)
    endif
    Time = set_date('2000-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test8.4 successful'
    else
       out_msg = 'test8.4 fails: ' // trim(err_msg)
       call mpp_error(FATAL, out_msg)
    endif

    call set_calendar_type(JULIAN)
    Time = set_date('1900-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test8.5 successful'
    else
       out_msg = 'test8.5 fails: ' // trim(err_msg)
       call mpp_error(FATAL, out_msg)
    endif
    Time = set_date('1901-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test8.6 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test8.6 successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
 ! Tests of days_in_month

  if(test9) then
    write(outunit,'(/,a)') '#################################  test9  #################################'
    day = days_in_month(set_date('1901-02-28 00:00:00'))
    write(outunit,'(a,i4)') ' test9.1: day=',day
    day = days_in_month(set_date('1901-07-01 00:00:00'))
    write(outunit,'(a,i4)') ' test9.2: day=',day
  endif
 !==============================================================================================
 ! Tests of get_time error flag

  if(test10) then
    write(outunit,'(/,a)') '#################################  test10  #################################'
    Time = set_time(seconds=2, days=1, ticks=1)
    call get_time(Time, seconds=sec, days=day, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test10.1 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test10.1 successful: '//trim(err_msg)
    endif
    call set_calendar_type(GREGORIAN)
    Time = set_time(seconds=2, days=1, ticks=1)
    call get_date(Time, yr, mo, day, hr, min, sec, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test10.2 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test10.2 successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
 ! Tests of increment_time and decrement_time

  if(test11) then
    write(outunit,'(/,a)') '#################################  test11  #################################'
    call print_time(increment_time(set_time(seconds=0, days=2), seconds=0, days=1),'test11.1:', unit=outunit)
    call print_time(decrement_time(set_time(seconds=0, days=2), seconds=0, days=1),'test11.2:', unit=outunit)
    call print_time(increment_time(set_time(seconds=0, days=2, ticks=5), seconds=400, days=1, ticks=14),'test11.3:', unit=outunit)
    call print_time(decrement_time(set_time(seconds=0, days=2, ticks=5), seconds=400, days=1, ticks=14),'test11.4:', unit=outunit)
  endif
 !==============================================================================================
 !  Tests of negative increments in increment_time and decrement_time

  if(test12) then
    write(outunit,'(/,a)') '#################################  test12  #################################'
    call print_time(increment_time(set_time(seconds=0, days=2), seconds=0, days=-1),'test12.1:', unit=outunit)
    call print_time(decrement_time(set_time(seconds=0, days=2), seconds=0, days=-1),'test12.2:', unit=outunit)
    call print_time(increment_time(set_time(seconds=0, days=2, ticks=5),seconds=-400,days=-1,ticks=-14),'test12.3:',unit=outunit)
    call print_time(decrement_time(set_time(seconds=0, days=2, ticks=5),seconds=-400,days=-1,ticks=-14),'test12.4:',unit=outunit)
  endif
 !==============================================================================================
 !  Test of trap for negative time

  if(test13) then
    write(outunit,'(/,a)') '#################################  test13  #################################'
    Time = set_time(seconds= 2, days=0, ticks=-21, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test13.1 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test13.1 successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
 !  Tests of negative seconds and/or ticks

  if(test14) then
    write(outunit,'(/,a)') '#################################  test14  #################################'
    call print_time(set_time(seconds=-86399, days=2, ticks=-10),'test14.1:', unit=outunit)
    call print_time(set_time(seconds=-86390, days=2, ticks=-95),'test14.2:', unit=outunit)
    call print_time(set_time(seconds= 86400, days=2, ticks= 95),'test14.3:', unit=outunit)
  endif
 !==============================================================================================
 !  Tests of consistency of day numbering between calendars

  if(test15) then
    write(outunit,'(/,a)') '#################################  test15  #################################'
    call set_calendar_type(GREGORIAN)
    Time = set_date(1, 1, 1)
    call get_time(Time, sec, day)
    write(outunit,10) 'GREGORIAN',day

    call set_calendar_type(JULIAN)
    Time = set_date(1, 1, 1)
    call get_time(Time, sec, day)
    write(outunit,10) 'JULIAN',day

    call set_calendar_type(THIRTY_DAY_MONTHS)
    Time = set_date(1, 1, 1)
    call get_time(Time, sec, day)
    write(outunit,10) 'THIRTY_DAY_MONTHS',day

    call set_calendar_type(NOLEAP)
    Time = set_date(1, 1, 1)
    call get_time(Time, sec, day)
    write(outunit,10) 'NOLEAP',day
  endif

  10 format(a17,' Jan 1 year 1 is day=',i6)

 !==============================================================================================
 ! Tests of error message for invalid dates

  if(test16) then
    write(outunit,'(/,a)') '#################################  test16  #################################'
    call set_calendar_type(GREGORIAN)
    Time = set_date(1900, 1, 32, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test16.1 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test16.1 successful: '//trim(err_msg)
    endif

    Time = set_date(1900, 4, 31, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test16.2 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test16.2 successful: '//trim(err_msg)
    endif

    Time = set_date(1900, 2, 29, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test16.3 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test16.3 successful: '//trim(err_msg)
    endif

    call set_calendar_type(JULIAN)
    Time = set_date(1900, 1, 0, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test16.4 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test16.4 successful: '//trim(err_msg)
    endif

    call set_calendar_type(NOLEAP)
    Time = set_date(1900, 0, 1, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test16.5 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test16.5 successful: '//trim(err_msg)
    endif

    Time = set_date(1900, 1, 1, tick=11, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test16.6 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test16.6 successful: '//trim(err_msg)
    endif

    call set_calendar_type(THIRTY_DAY_MONTHS)
    Time = set_date(1900, 13, 1, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test16.7 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test16.7 successful: '//trim(err_msg)
    endif

    Time = set_date(1900, 12, 31, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test16.8 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test16.8 successful: '//trim(err_msg)
    endif

    call set_calendar_type(JULIAN)
    Time = set_date(1900, 4, 31, err_msg=err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test16.9 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test16.9 successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
 !  Tests of Gregorian calendar
 !  This test loops through every day of an 400 year period and writes a line to the output file for each day.

  if(test17) then
    write(outunit,'(/,a)') '#################################  test17  #################################'
    write(errunit,'(/,a)') ' ====================================================='
    write(errunit,'(a)')   '  Warning: test17 produces voluminous output.'
    write(errunit,'(a)')   '  It can be turned off with: &test_nml test17=.false./'
    write(errunit,'(a,/)') ' ====================================================='
    call set_calendar_type(GREGORIAN)
    do year=1801,2200
      leap = mod(year,4) == 0
      leap = leap .and. .not.mod(year,100) == 0
      leap = leap .or. mod(year,400) == 0
      do month=1,12
        days_this_month = days_per_month(month)
        if(leap .and. month == 2) days_this_month = 29
        do dday=1,days_this_month
          Time = set_date(year, month, dday, 0, 0, 0)
          call get_date(Time, yr, mo, day, hr, min, sec)
          write(outunit,100) yr, mo, day, leap_year(Time), days_in_month(Time), days_in_year(Time)
        enddo
      enddo
    enddo
  endif
  100 format('yr=',i4,' mo=',i2,' day=',i2,' leap=',L1,' days_in_month=',i2,' days_in_year=',i3)
 !==============================================================================================
 !  Tests of length_of_year

  if(test18) then
    write(outunit,'(/,a)') '#################################  test18  #################################'
    call set_calendar_type(THIRTY_DAY_MONTHS)
    call print_time(length_of_year(), 'length_of_year for THIRTY_DAY_MONTHS:', unit=outunit)
    call set_calendar_type(NOLEAP)
    call print_time(length_of_year(), 'length_of_year for NOLEAP:', unit=outunit)
    call set_calendar_type(JULIAN)
    call print_time(length_of_year(), 'length_of_year for JULIAN:', unit=outunit)
    call set_calendar_type(GREGORIAN)
    call print_time(length_of_year(), 'length_of_year for GREGORIAN:', unit=outunit)
  endif
 !==============================================================================================
 !  Tests of real_to_time_type

  if(test19) then
    write(outunit,'(/,a)') '#################################  test19  #################################'
    call print_time(real_to_time_type(86401.1), 'real_to_time_type(86401.1):', unit=outunit)
    Time = real_to_time_type(-1.0, err_msg)
    if(err_msg == '') then
       call mpp_error(FATAL, 'test19.3 fails: did not get the expected error message')
    else
      write(outunit,'(a)') 'test successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
  write(outunit,'(/,a)') '############################################################################'
  write(outunit,'(a,i6)') ' ticks_per_second=',get_ticks_per_second()

 !==============================================================================================
 !  Tests the new set/get_date_gregorian by comparing against the old set/get_date_gregorian copied over to this test program
 !  This test loops through every day up to year 3200

  if(test20) then
    write(outunit,'(/,a)') '#################################  test20  #################################'
    write(outunit,'(/,a)') ' ====================================================='
    write(outunit,'(a)')   '  Test get/set_date_gregorian with get/set_date_gregorian_old'
    write(outunit,'(a,/)') ' ====================================================='
    call set_calendar_type(GREGORIAN)
    call get_coded_date( coded_date, date_to_day ) ! assign coded_date and date_to_day used by get/set_date_gregorian_old

    ! Check that the get/set_date_gregorian_old here are the same as in time_manager
    ! This part of the test will be deleted when the old methods are removed from time_manager
    do year=1, 3200
      leap = mod(year,4) == 0
      leap = leap .and. .not.mod(year,100) == 0
      leap = leap .or. mod(year,400) == 0
      do month=1,12
        days_this_month = days_per_month(month)
        if(leap .and. month == 2) days_this_month = 29
        do dday=1,days_this_month
          ! test set_date_gregorian
          Time  = set_date(year, month, dday, 0, 0, 0, old_method=.true.)
          Time0 = set_date_gregorian_old(year, month, dday, 0, 0, 0, 0, date_to_day)
          if( .not. (Time==Time0) ) then
            write(outunit,'("ERROR with year",i5,"mo",i5,"dday",i5)') year, month, dday
            call mpp_error(FATAL, 'ERROR testing set_date_gregorian_old:  Time!=Time0')
          end if
          ! test #1 get_date
          call get_date(Time0, yr, mo, day, hr, min, sec, old_method=.true.)
          call get_date_gregorian_old(Time0, coded_date, yr0, mo0, day0, hr0, min0, sec0, ticks0)
          if( yr0.ne.yr .or. mo0.ne.mo .or. day0.ne.day ) then
            write(outunit,"('expected year ',i5,'but got year ',i5)") yr0, yr
            write(outunit,"('expected month',i5,'but got month',i5)") mo0, mo
            write(outunit,"('expected day  ',i5,'but got day  ',i5)") day0, day
            call mpp_error(FATAl,'Error testing get_date_gregorian_old 1')
          end if
          ! test #2 get_date
          call get_date(Time, yr, mo, day, hr, min, sec, old_method=.true.)
          call get_date_gregorian_old(Time, coded_date, yr0, mo0, day0, hr0, min0, sec0, ticks0)
          if( yr0.ne.yr .or. mo0.ne.mo .or. day0.ne.day ) then
            write(outunit,"('expected year ',i5,'but got year ',i5)") yr0, yr
            write(outunit,"('expected month',i5,'but got month',i5)") mo0, mo
            write(outunit,"('expected day  ',i5,'but got day  ',i5)") day0, day
            call mpp_error(FATAl,'Error testing get_date_gregorian 2')
          end if
          ! test #3 get_date
          call get_date(Time, yr, mo, day, hr, min, sec, old_method=.true.)
          call get_date_gregorian_old(Time0, coded_date, yr0, mo0, day0, hr0, min0, sec0, ticks0)
          if( yr0.ne.yr .or. mo0.ne.mo .or. day0.ne.day ) then
            write(outunit,"('expected year ',i5,'but got year ',i5)") yr0, yr
            write(outunit,"('expected month',i5,'but got month',i5)") mo0, mo
            write(outunit,"('expected day  ',i5,'but got day  ',i5)") day0, day
            call mpp_error(FATAl,'Error testing get_date_gregorian 3')
          end if
          ! test #4 get_date
          call get_date(Time0, yr, mo, day, hr, min, sec, old_method=.true.)
          call get_date_gregorian_old(Time, coded_date, yr0, mo0, day0, hr0, min0, sec0, ticks0)
          if( yr0.ne.yr .or. mo0.ne.mo .or. day0.ne.day ) then
            write(outunit,"('expected year ',i5,'but got year ',i5)") yr0, yr
            write(outunit,"('expected month',i5,'but got month',i5)") mo0, mo
            write(outunit,"('expected day  ',i5,'but got day  ',i5)") day0, day
            call mpp_error(FATAl,'Error testing get_date_gregorian 4')
          end if
        enddo
      enddo
    enddo
    write(outunit,'(a)') 'set_date_gregorian_old and get_date_gregorian_old tests successful'

    ! test the new Gregorian methods and compare with the old methods
    do year=1, 3200
      leap = mod(year,4) == 0
      leap = leap .and. .not.mod(year,100) == 0
      leap = leap .or. mod(year,400) == 0
      do month=1,12
        days_this_month = days_per_month(month)
        if(leap .and. month == 2) days_this_month = 29
        do dday=1,days_this_month
          ! test new set_date_gregorian
          Time  = set_date(year, month, dday, 0, 0, 0)
          Time0 = set_date_gregorian_old(year, month, dday, 0, 0, 0, 0, date_to_day)
          if( .not. (Time==Time0) ) then
             write(outunit,'("ERROR with year",i5,"mo",i5,"dday",i5)') year, month, dday
             call mpp_error(FATAL, 'ERROR testing set_date_gregorian:  Time!=Time0')
          end if
          ! test #1 get_date
          call get_date(Time0, yr, mo, day, hr, min, sec)
          call get_date_gregorian_old(Time0, coded_date, yr0, mo0, day0, hr0, min0, sec0, ticks0)
          if( yr0.ne.yr .or. mo0.ne.mo .or. day0.ne.day ) then
            write(outunit,"('expected year ',i5,'but got year ',i5)") yr0, yr
            write(outunit,"('expected month',i5,'but got month',i5)") mo0, mo
            write(outunit,"('expected day  ',i5,'but got day  ',i5)") day0, day
            call mpp_error(FATAl,'Error testing get_date_gregorian 1')
          end if
          ! test #2 get_date
          call get_date(Time, yr, mo, day, hr, min, sec)
          call get_date_gregorian_old(Time, coded_date, yr0, mo0, day0, hr0, min0, sec0, ticks0)
          if( yr0.ne.yr .or. mo0.ne.mo .or. day0.ne.day ) then
            write(outunit,"('expected year ',i5,'but got year ',i5)") yr0, yr
            write(outunit,"('expected month',i5,'but got month',i5)") mo0, mo
            write(outunit,"('expected day  ',i5,'but got day  ',i5)") day0, day
            call mpp_error(FATAl,'Error testing get_date_gregorian 2')
          end if
          ! test #3 get_date
          call get_date(Time, yr, mo, day, hr, min, sec)
          call get_date_gregorian_old(Time0, coded_date, yr0, mo0, day0, hr0, min0, sec0, ticks0)
          if( yr0.ne.yr .or. mo0.ne.mo .or. day0.ne.day ) then
            write(outunit,"('expected year ',i5,'but got year ',i5)") yr0, yr
            write(outunit,"('expected month',i5,'but got month',i5)") mo0, mo
            write(outunit,"('expected day  ',i5,'but got day  ',i5)") day0, day
            call mpp_error(FATAl,'Error testing get_date_gregorian 3')
          end if
          ! test #4 get_date
          call get_date(Time0, yr, mo, day, hr, min, sec)
          call get_date_gregorian_old(Time, coded_date, yr0, mo0, day0, hr0, min0, sec0, ticks0)
          if( yr0.ne.yr .or. mo0.ne.mo .or. day0.ne.day ) then
            write(outunit,"('expected year ',i5,'but got year ',i5)") yr0, yr
            write(outunit,"('expected month',i5,'but got month',i5)") mo0, mo
            write(outunit,"('expected day  ',i5,'but got day  ',i5)") day0, day
            call mpp_error(FATAl,'Error testing get_date_gregorian 4')
          end if
        enddo
      enddo
    enddo
    write(outunit,'(a)') 'set_date_gregorian and get_date_gregorian tests successful'
 endif

  call fms_end

contains

  ! get_coded_date:  copied from subroutine set_calendar_type in time_manager and slightly modified
  ! to work in this test program.
  subroutine get_coded_date(coded_date_old, date_to_day_old)

    implicit none

    integer, intent(out), dimension(146097) :: coded_date_old
    integer, intent(out), dimension(400,12,31) :: date_to_day_old

    integer :: iday, days_this_month, year, month, day
    logical :: leap

    iday = 0
     date_to_day = -1 ! invalid_date = -1 in time_manager
     do year=1,400
       leap = mod(year,4) == 0
       leap = leap .and. .not.mod(year,100) == 0
       leap = leap .or. mod(year,400) == 0
       do month=1,12
         days_this_month = days_per_month(month)
         if(leap .and. month ==2) days_this_month = 29
         do day=1,days_this_month
           date_to_day_old(year,month,day) = iday
           iday = iday+1
           coded_date_old(iday) = day + 32*(month + 16*year)
         enddo ! do day
       enddo ! do month
     enddo ! do year

  end subroutine get_coded_date

  ! get_date_gregorian_old:  original get_date_gregorian subroutine in time_manager that has been slightly
  ! modified to work in this test program
  subroutine get_date_gregorian_old(time, coded_date, year, month, day, hour, minute, second, tick)

    use time_manager_mod, only : set_time

    integer, parameter :: days_in_400_year_period = 146097

    type(time_type), intent(in) :: time
    integer, intent(in), dimension(days_in_400_year_period) :: coded_date
    integer, intent(out) :: year, month, day, hour, minute, second
    integer, intent(out) :: tick

    integer :: iday, isec, time_days, time_seconds, time_ticks

    ! set time_days=Time%days and time_seconds=Time%seconds, time_ticks=Time%ticks
    call get_time(Time, seconds=time_seconds, days=time_days, ticks=time_ticks)

    iday = mod(time_days+1, days_in_400_year_period)
    if(iday == 0) iday = days_in_400_year_period

    year = coded_date(iday)/512
    day = mod(coded_date(iday),32)
    month = coded_date(iday)/32 - 16*year

    year = year + 400*(time_days/days_in_400_year_period)

    hour   = time_seconds / 3600
    isec   = time_seconds - 3600*hour
    minute = isec / 60
    second = isec - 60*minute
    tick   = time_ticks

  end subroutine get_date_gregorian_old

  ! set_date_gregorian_old: original set_date_gregorian function in time_manager that has been slightly
  ! modified to work in this test program
  function set_date_gregorian_old(year, month, day, hour, minute, second, tick, date_to_day)

    use time_manager_mod, only: set_time

    type(time_type) :: set_date_gregorian_old

    integer, parameter :: days_in_400_year_period = 146097

    integer, intent(in)  :: year, month, day, hour, minute, second, tick
    integer, intent(in),  dimension(400,12,31) :: date_to_day

    integer :: yr1, day1, second1

    second1 = second + 60*(minute + 60*hour)

    yr1 = mod(year,400)
    if(yr1 == 0) yr1 = 400
    day1 = date_to_day(yr1,month,day)

    day1 = day1 + days_in_400_year_period*((year-1)/400)

    set_date_gregorian_old = set_time(seconds=second1, days=day1, ticks=tick)
    
  end function set_date_gregorian_old

end program test_time_manager
