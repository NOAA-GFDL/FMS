!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @file
!! @brief Tests the clock_init subroutine
!! @author Colin Gladue
!! @email gfdl.climate.model.info@noaa.gov

program test_clock_init

  use mpp_mod, only : mpp_init, mpp_init_test_init_true_only
  use mpp_mod, only : mpp_error, FATAL
  use mpp_mod, only : mpp_clock_id, MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED

  integer :: ierr

  call mpp_init(test_level=mpp_init_test_init_true_only)
  ! We cannot directly access the clock_init subroutine,
  ! but we can use the mpp_clock_id function to call clock_init.
  write(*,*) "Testing simple clock"
  call create_and_check_clock("name1", 1) ! Simple name input clock init
  write(*,*) "Testing clock with empty string name"
  call create_and_check_clock("", 2)      ! Simple clock init with empty name
  ! All possible combinations of specified flags
  write(*,*) "Testing clock with MPP_CLOCK_SYNC flag"
  call create_and_check_clock("name3", 3, flagsIn=MPP_CLOCK_SYNC )
  write(*,*) "Testing clock with MPP_CLOCK_DETAILED flag"
  call create_and_check_clock("name4", 4, flagsIn=MPP_CLOCK_DETAILED)
  write(*,*) "Testing clock with both flags"
  call create_and_check_clock("name5", 5, flagsIn=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED)
  call MPI_FINALIZE(ierr)

  contains
    !> @brief Helper subroutine to test different scenarios when testing mpp's clock_init subroutine
    subroutine create_and_check_clock(nameIn, expected, flagsIn)
      character(len=*), intent(in)  :: nameIn !< Name of clock
      integer, intent(in)           :: expected !< Expected clock id
      integer, intent(in), optional :: flagsIn !< Clock creation flags
      integer                       :: dummy !< Dummy variable to aid in calling mpp_clock_id function
      integer                       :: id !< Id of a given clock
      if (PRESENT(flagsIN)) then
        dummy = mpp_clock_id(nameIn, flags=flagsIn)
      else
        dummy = mpp_clock_id(nameIn)
      end if
      ! At this point we have not encountered any errors after initializing the clock
      ! Let's make sure this clock can be accessed and has the expected ID.
      id = mpp_clock_id(nameIn)
      if (id.ne.expected) then
        call mpp_error(FATAL, "Clock not initialized correctly, ID != expected")
      end if
    end subroutine create_and_check_clock

end program test_clock_init
