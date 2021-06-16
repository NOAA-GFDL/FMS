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
!> @defgroup mpp_data_mod mpp_data_mod
!> @ingroup mpp
!> @brief Module to hold pointer and stack data for use in @ref mpp modules.
!!
!> Makes stack and pointer data publicly available from @ref mpp_data_mpi.inc or @ref
!! mpp_data_nocomm.inc for use in @ref mpp modules. This module is mainly
!! for internal use within @ref mpp_mod and @ref mpp_domains_mod .

!> @file
!> @brief File for @ref mpp_data_mod

!> @addtogroup mpp_data_mod
!> @{
module mpp_data_mod

#if defined(use_libMPI)
  use mpi
#endif

  use mpp_parameter_mod, only : MAXPES
  use platform_mod

  implicit none
  private

! Include variable "version" to be written to log file.
#include<file_version.h>
  public version

  !> public data used by mpp_mod
  public :: stat, mpp_stack, ptr_stack, status, ptr_status, sync, ptr_sync
  public :: mpp_from_pe, ptr_from, remote_data_loc, ptr_remote

  !--- All othere modules should import these parameters from mpp_domains_mod.
  !> public data which is used by mpp_domains_mod.
  public :: mpp_domains_stack, ptr_domains_stack
  public :: mpp_domains_stack_nonblock, ptr_domains_stack_nonblock

  !-------------------------------------------------------------------------------!
  ! The following data included in the .inc file are diffrent for sma or mpi case !
  !-------------------------------------------------------------------------------!

#ifdef use_libMPI
#include <mpp_data_mpi.inc>
#else
#include <mpp_data_nocomm.inc>
#endif

end module mpp_data_mod
!> @}
! close documentation grouping
