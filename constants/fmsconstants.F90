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
!> @defgroup fmsconstants FMSConstants
!> @ingroup libfms
!> @brief Defines useful constants for Earth. Constants are defined as real
!!
!>    FMSconstants have been declared as REAL(kind=sizeof(rvar)), PARAMETER.
!!
!!    The value of a constant defined and used from here cannot be changed
!!    in a users program. New constants can be defined in terms of values
!!    from the FMSconstants module and their includes using a parameter
!!    statement.<br><br>
!!
!!    The currently support contant systems are:
!!       GFDL constants (gfdl_constants.fh)
!!       GEOS constants (geos_constants.fh)
!!       GFS  constants (gfs_constants.fh)
!!       <br><br>
!!
!!    The name given to a particular constant may be changed.<br><br>
!!
!!    Constants can only be used on the right side on an assignment statement
!!    (their value can not be reassigned).
!!
!!    Example:
!!
!! @verbatim
!!    use FMSConstants, only:  TFREEZE, grav_new => GRAV
!!    real, parameter :: grav_inv = 1.0 / grav_new
!!    tempc(:,:,:) = tempk(:,:,:) - TFREEZE
!!    geopotential(:,:) = height(:,:) * grav_new
!! @endverbatim
!> @file
!> @brief File for @ref FMSconstants_mod

!> @addtogroup FMSconstants_mod
!> @{
module FMSconstants

  use platform_mod, only: r4_kind, r8_kind

  !--- default scoping
  implicit none

  !--- needed with implicit none
  real :: dum  !< dummy real variable

#define RKIND sizeof(dum)

!--- set a default for the FMSConstants
#if !defined(GFDL_CONSTANTS) && !defined(GFS_CONSTANTS) && !defined(GEOS_CONSTANTS)
#define GFDL_CONSTANTS
#endif

!--- perform error checking and include the correct system of constants
#if defined(GFDL_CONSTANTS) && !defined(GFS_CONSTANTS) && !defined(GEOS_CONSTANTS)
#warning "Using GFDL constants"
#include <gfdl_constants.fh>
#elif !defined(GFDL_CONSTANTS) && defined(GFS_CONSTANTS) && !defined(GEOS_CONSTANTS)
#warning "Using GFS constants"
#include <gfs_constants.fh>
#elif !defined(GFDL_CONSTANTS) && !defined(GFS_CONSTANTS) && defined(GEOS_CONSTANTS)
#warning "Using GEOS constants"
#include <geos_constants.fh>
#else
#error FATAL FMSConstants error -  multiple constants macros are defined for FMS
#endif

  !--- public interfaces
  public :: FMSConstants_init

  contains

    !> @brief FMSconstants init routine
    subroutine FMSconstants_init
      use mpp_mod, only: stdlog
      integer :: logunit
      logunit = stdlog()

      write (logunit,'(/,80("="),/(a))') trim(constants_version)

    end subroutine FMSconstants_init

end module FMSconstants
!> @}
! close documentation grouping
