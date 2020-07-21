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

!> @file
!! @brief A collection of tests for mpp_clock_begin, mpp_clock_end, mpp_clock_id
!! @author Christopher Dupuis
!! @email gfdl.climate.model.info@noaa.gov
program test_mpp_clock_begin_end_id

  use mpp_mod, only : mpp_init, mpp_exit
  use mpp_mod, only : mpp_error, FATAL
  use mpp_mod, only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod, only : mpp_record_time_start, mpp_record_time_end
  use mpp_mod, only : CLOCK_LOOP
  use mpp_parameter_mod, only : MAX_CLOCKS

  integer :: test_number

  namelist / test_mpp_clock_begin_end_id_nml / test_number

  open(3, file="clock.nml", form="formatted", status="old")
  read(3, nml = test_mpp_clock_begin_end_id_nml)
  close(3)

  select case(test_number)
    case(1)
      call test1()
    case(2)
      call test2()
    case(3)
      call test3()
    case(4)
      call test4()
    case(5)
      call test5()
    case(6)
      continue!call test6()
    case(7)
      call test7()
    !case(8)
    !  call test8()
    case(9)
      call test9()
    case(10)
      call test10()
    case(11)
      call test11()
    case(12)
      call test12()
    case(13)
      call test13()
    case(14)
      call test14()
    case(15)
      call test15()
    case(16)
      call test16()
    case(17)
      call test17()
    case default
      call mpp_init() !LEVEL 0
      call mpp_error(FATAL, "ERROR: No test was selected.")
  end select
  contains

    subroutine segfault()

      logical,pointer :: segfault_array => NULL()
      segfault_array = .true.

    end subroutine segfault


    subroutine test1()
      integer :: clock_id
      integer :: bp
      bp = 1

      call mpp_init()!test_level=bp)
      clock_id = mpp_clock_id("Ultraclock")
      call mpp_exit()

    end subroutine test1

    subroutine test2()
      integer :: clock_id
      integer :: bp
      bp = 1

      call mpp_init()!test_level=bp)
      clock_id = mpp_clock_id("Ultraclock", flags=1)
      call mpp_exit()

    end subroutine test2

    subroutine test3()
      integer :: clock_id
      integer :: bp
      bp = 1

      call mpp_init()!test_level=bp)
      clock_id = mpp_clock_id("Ultraclock", flags=2)
      call mpp_exit()

    end subroutine test3

    subroutine test4()
      integer :: clock_id
      integer :: bp
      bp = 1

      call mpp_init()!test_level=bp)
      clock_id = mpp_clock_id("Ultraclock", flags=3)
      call mpp_exit()

    end subroutine test4

    subroutine test5()
      integer :: clock_id
      integer :: bp
      bp = 1

      clock_id = mpp_clock_id("Ultraclock")
      call mpp_init()!test_level=bp)
      call mpp_exit()

    end subroutine test5

    subroutine test6()
      integer :: clock_id
      integer :: bp
      bp = 1

      call mpp_init()!test_level=bp)
      call mpp_exit()
      clock_id = mpp_clock_id("Ultraclock")

    end subroutine test6

    subroutine test7()
      integer :: clock_id
      integer :: new_grain
      integer :: bp
      bp = 1

      new_grain = CLOCK_LOOP !clock_grain + 1

      call mpp_init()!test_level=bp)
      clock_id = mpp_clock_id("Ultraclock", grain = new_grain)
      call mpp_exit()

      if(clock_id .ne. 0) then
          call segfault()
      end if

    end subroutine test7

    !subroutine test8()
    !  integer :: clock_id
    !  integer :: clock_num_tmp
    !  integer :: bp
    !  bp = 1

    !  call mpp_init()!test_level=bp)
    !  clock_num_tmp = clock_num
    !  clock_num = MAX_CLOCKS + 1
    !  clock_id = mpp_clock_id("Ultraclock")
    !  clock_num = clock_num_tmp
    !  call mpp_exit()

    !end subroutine test8
    
    subroutine test9()
      integer :: clock_id
      integer :: bp
      bp = 6

      clock_id = 1
      call mpp_clock_begin(clock_id)
      call mpp_init()!test_level=bp)
      call mpp_exit()

    end subroutine test9

    subroutine test10()
      integer :: clock_id
      integer :: bp
      bp = 6

      clock_id = 1
      call mpp_init()!test_level=bp)
      call mpp_exit()
      call mpp_clock_begin(clock_id)

    end subroutine test10

    subroutine test11()
      integer :: clock_id
      integer :: bp
      bp = 6

      call mpp_init()!test_level=bp)
      clock_id = mpp_clock_id("Ultraclock")
      call mpp_clock_begin(clock_id)
      call mpp_exit()

    end subroutine test11

    subroutine test12()
      integer :: clock_id
      integer :: bp
      bp = 6

      call mpp_init()!test_level=bp)
      clock_id = 0
      call mpp_clock_begin(clock_id)
      call mpp_exit()

    end subroutine test12

    subroutine test13()
      integer :: clock_id
      integer :: bp
      bp = 6

      call mpp_init()!test_level=bp)
      clock_id = MAX_CLOCKS + 10
      call mpp_clock_begin(clock_id)
      call mpp_exit()

    end subroutine test13

    subroutine test14()
      integer :: clock_id
      integer :: bp
      bp = 6

      call mpp_init()!test_level=bp)
      clock_id = mpp_clock_id("Ultraclock")
      call mpp_clock_begin(clock_id)
      call mpp_clock_end(clock_id)
      call mpp_exit()

    end subroutine test14

    subroutine test15()
      integer :: clock_id
      integer :: bp
      bp = 6

      call mpp_init()!test_level=bp)
      clock_id = mpp_clock_id("Ultraclock")
      call mpp_clock_end(clock_id)
      call mpp_exit()

    end subroutine test15

    subroutine test16()
      integer :: clock_id
      integer :: bp
      bp = 6

      call mpp_init()!test_level=bp)
      clock_id = mpp_clock_id("Ultraclock")
      call mpp_clock_begin(clock_id)
      call mpp_clock_begin(clock_id)
      call mpp_exit()

    end subroutine test16

    subroutine test17()
      integer :: clock_id
      integer :: bp
      bp = 6

      call mpp_init()!test_level=bp)
      clock_id = MAX_CLOCKS + 10
      call mpp_clock_end(clock_id)
      call mpp_exit()

    end subroutine test17

end program test_mpp_clock_begin_end_id
