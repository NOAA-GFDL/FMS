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
