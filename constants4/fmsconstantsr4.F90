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
!> @defgroup fmsconstantsR4 FMSConstantsR4
!> @ingroup libfms
!> @brief Defines useful constants for Earth. Constants are defined as real
!!
!>    FMSconstantsR4 have been declared as r4_kind or r8_kind PARAMETER.
!!
!!    The value of a constant defined and used from here cannot be changed
!!    in a users program. New constants can be defined in terms of values
!!    from the FMSconstants module and their includes using a parameter
!!    statement.<br><br>
!!
!!    The currently support contant systems are:
!!       GFDL constants (gfdl_constantsR4.fh)
!!       GEOS constants (geos_constantsR4.fh)
!!       GFS  constants (gfs_constantsR4.fh)
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
!!    use FMSConstantsR4, only:  TFREEZE, grav_new => GRAV
!!    real, parameter :: grav_inv = 1.0 / grav_new
!!    tempc(:,:,:) = tempk(:,:,:) - TFREEZE
!!    geopotential(:,:) = height(:,:) * grav_new
!! @endverbatim
!> @file
!> @brief File for @ref FMSconstantsR4_mod

!> @addtogroup FMSconstantsR4_mod
!> @{
module FMSconstantsR4

  use platform_mod, only: r4_kind, r8_kind

  !--- default scoping
  implicit none

#define RKIND r4_kind

!--- set a default for the FMSConstantsR4
#if !defined(GFDL_CONSTANTS) && !defined(GFS_CONSTANTS) && !defined(GEOS_CONSTANTS)
#define GFDL_CONSTANTS
#endif

!--- perform error checking and include the correct system of constants
#if defined(GFDL_CONSTANTS) && !defined(GFS_CONSTANTS) && !defined(GEOS_CONSTANTS)
#warning "Using GFDL constantsR4"
#include <gfdl_constantsR4.fh>
#elif !defined(GFDL_CONSTANTS) && defined(GFS_CONSTANTS) && !defined(GEOS_CONSTANTS)
#warning "Using GFS constantsR4"
#include <gfs_constantsR4.fh>
#elif !defined(GFDL_CONSTANTS) && !defined(GFS_CONSTANTS) && defined(GEOS_CONSTANTS)
#warning "Using GEOS constantsR4"
#include <geos_constantsR4.fh>
#else
#error FATAL FMSConstantsR4 error -  multiple constants macros are defined for FMS
#endif

  !--- public interfaces
  public :: FMSConstantsR4_init

  contains

    !> @brief FMSconstantsR4 init routine
    subroutine FMSconstantsR4_init
      use mpp_mod, only: stdlog
      integer :: logunit
      logunit = stdlog()

      write (logunit,'(/,80("="),/(a))') trim(constantsR4_version)

    end subroutine FMSconstantsR4_init

end module FMSconstantsR4
!> @}
! close documentation grouping
