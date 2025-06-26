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
!> @defgroup constantsR4_mod constantsR4_mod
!> @ingroup constantsR4
!> @brief compatibility module as we transition to an FMSConstantsR4 module
!!
!> @file
!> @brief File for @ref constantsR4_mod

module constantsR4_mod

!> rename to not conflict with any other version vars
use FMSConstantsR4, version => constantsR4_version, constants_init => FMSconstantsR4_init

contains

end module constantsR4_mod
