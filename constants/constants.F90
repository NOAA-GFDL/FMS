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
!> @brief File for @ref constants_mod

!> @defgroup constants_mod constants_mod
!> @ingroup constants
!> @brief compatibility module as we transition to an FMSConstants module
!!
!>    Constants have been declared as type REAL, PARAMETER.
!!
!!    The value a constant can not be changed in a users program.
!!    New constants can be defined in terms of values from the
!!    constants module using a parameter statement.<br><br>
!!
!!    The name given to a particular constant may be changed.<br><br>
!!
!!    Constants can be used on the right side on an assignment statement
!!    (their value can not be reassigned).
!!
!!    Example:
!!
!! @verbatim
!!    use constants_mod, only:  TFREEZE, grav_new => GRAV
!!    real, parameter :: grav_inv = 1.0 / grav_new
!!    tempc(:,:,:) = tempk(:,:,:) - TFREEZE
!!    geopotential(:,:) = height(:,:) * grav_new
!! @endverbatim

!> @addtogroup constants_mod
!> @{
module constants_mod

!> rename to not conflict with any other version vars
use FMSConstants, version => constants_version, constants_init => FMSconstants_init

contains

end module constants_mod
