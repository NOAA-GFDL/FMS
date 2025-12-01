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
!! @brief Tests the mpp_init_logfile interface.
!! The files to be replaced are created and checked by the calling
!! shell program test_mpp_init_logfile.sh.
!!
!! @author Miguel Zuniga
!! @email gfdl.climate.model.info@noaa.gov

program test_mpp_init_logfile

  use mpp_mod, only : mpp_init
  use mpp_mod, only : mpp_init_test_logfile_init

  IMPLICIT NONE

  integer :: err_no

  ! Initialize mpp but have mpp_init() exit right after (or soon after) it is called.
  ! Note mpp_init may write to te root log file once.
  call mpp_init( test_level =  mpp_init_test_logfile_init)

  ! With the unifinished initialization, mpp_exit() may cause a crash. Use MPI_FINALIZE:
  call MPI_FINALIZE(err_no)

end program test_mpp_init_logfile
