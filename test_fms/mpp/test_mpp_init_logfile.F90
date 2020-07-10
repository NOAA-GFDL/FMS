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
!! @brief Tests the mpp_init_logfile interface
!! @author Miguel Zuniga
!! @email gfdl.climate.model.info@noaa.gov

program test_mpp_init_logfile
#include <fms_platform.h>

  use mpp_mod, only : mpp_init, mpp_init_logfile
  use mpp_mod, only : mpp_init_test_logfile_init

  IMPLICIT NONE

  integer :: my_mpp_init_level
  integer :: err_no

  !Initialize mpp but have mpp_init() exit before log files are initalized.
  my_mpp_init_level = mpp_init_test_logfile_init - 1
  call mpp_init( test_level = my_mpp_init_level )
  
  call mpp_init_logfile()
  !! With the unifinished initialization, mpp_exit() is causing crash. Use MPI_FINALIZE:
  call MPI_FINALIZE(err_no)

  end program test_mpp_init_logfile
