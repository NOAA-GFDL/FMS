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
!! @brief unit test for the mpp_root_pe() function
!! @author MiKyung Lee
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program tests the function mpp_root_pe().
!! The value of root_pe is set to 0 upon declaration.
!! Thus, the test passes if the return value of mpp_root_pe is 0 and fails if otherwise.


program test_mpp_root_pe


  use mpp_mod, only:  mpp_init, FATAL, mpp_error, mpp_root_pe, mpp_init_test_init_true_only

  implicit none
  integer :: my_root_pe, test_root_pe, ierr


  !> call mpp_init at the lowest level
  call mpp_init( test_level=mpp_init_test_init_true_only )

  !> assign my_root_pe=root_pe=0
  my_root_pe = 0
  !> call mpp_root_pe()
  test_root_pe = mpp_root_pe()
  !> call mpp_error if mpp_root_pe did not return value eq to my_root_pe
  if( test_root_pe .ne. my_root_pe ) call mpp_error(FATAL, "mpp_root_pe does not equal root_pe")

  !> end
  call MPI_FINALIZE(ierr)


end program test_mpp_root_pe
