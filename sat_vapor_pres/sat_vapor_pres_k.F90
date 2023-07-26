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
!> @defgroup sat_vapor_pres_k_mod sat_vapor_pres_k_mod
!> @ingroup sat_vapor_pres
!> @brief Kernel module to be used by @ref sat_vapor_pres_mod for
!! table lookups and calculations

 module sat_vapor_pres_k_mod

! This module is what I (pjp) think a kernel should be.
! There have been many proposals as to what a kernel should look like.
! If fact, so many different ideas have been expressed that the lack
! of agreement has greatly hampered progress.
! The only way to move forward is to limit the requirments for a kernel
! to only what is widely agreeded upon.
! I believe that there are only two things widely agreeded upon.

! 1) A kernel should be independent of the rest of FMS so that it can
!    easily be ported into another programming system.
!    This requires that a kernel does not access anything by use association.
!    The one exception is this kernel, because it is not practical for physics
!    modules to avoid using a module that computes the saturation vapor
!    pressure of water vapor.

! 2) For the sake of thread safety, module globals should be written only at initialization.
!    In this case, the module globals are the tables and a handful of scalars.

! 3) A kernel should not read from an external file.

! One of the things that was not widely agreeded upon is that a kernel should
! not be a fortran module. This complicates things greatly for questionable
! benefit and could be done as a second step anyway, if necessary.

 use platform_mod, only : r4_kind, r8_kind

 implicit none
 private

! Include variable "version" to be written to log file.
#include<file_version.h>

 public :: sat_vapor_pres_init_k
 public :: lookup_es_k
 public :: lookup_des_k
 public :: lookup_es_des_k
 public :: lookup_es2_k
 public :: lookup_des2_k
 public :: lookup_es2_des2_k
 public :: lookup_es3_k
 public :: lookup_des3_k
 public :: lookup_es3_des3_k
 public :: compute_qs_k
 public :: compute_mrs_k

 !> @ingroup sat_vapor_pres_k_mod
 interface sat_vapor_pres_init_k
    module procedure sat_vapor_pres_init_k_r4
    module procedure sat_vapor_pres_init_k_r8
 end interface sat_vapor_pres_init_k

 !> @ingroup sat_vapor_pres_k_mod
 interface compute_es_k
    module procedure compute_es_k_r4
    module procedure compute_es_k_r8
 end interface compute_es_k

 interface compute_es_liq_k
    module procedure compute_es_liq_k_r4
    module procedure compute_es_liq_k_r8
 end interface compute_es_liq_k

 interface compute_es_liq_ice_k
    module procedure compute_es_liq_ice_k_r4
    module procedure compute_es_liq_ice_k_r8
 end interface compute_es_liq_ice_k

 !> @ingroup sat_vapor_pres_k_mod
 interface lookup_es_k
   module procedure lookup_es_k_0d_r4
   module procedure lookup_es_k_0d_r8
   module procedure lookup_es_k_1d_r4
   module procedure lookup_es_k_1d_r8
   module procedure lookup_es_k_2d_r4
   module procedure lookup_es_k_2d_r8
   module procedure lookup_es_k_3d_r4
   module procedure lookup_es_k_3d_r8
 end interface

 !> @ingroup sat_vapor_pres_k_mod
 interface lookup_des_k
   module procedure lookup_des_k_0d_r4
   module procedure lookup_des_k_0d_r8
   module procedure lookup_des_k_1d_r4
   module procedure lookup_des_k_1d_r8
   module procedure lookup_des_k_2d_r4
   module procedure lookup_des_k_2d_r8
   module procedure lookup_des_k_3d_r4
   module procedure lookup_des_k_3d_r8
 end interface

 !> @ingroup sat_vapor_pres_k_mod
 interface lookup_es_des_k
   module procedure lookup_es_des_k_0d_r4
   module procedure lookup_es_des_k_0d_r8
   module procedure lookup_es_des_k_1d_r4
   module procedure lookup_es_des_k_1d_r8
   module procedure lookup_es_des_k_2d_r4
   module procedure lookup_es_des_k_2d_r8
   module procedure lookup_es_des_k_3d_r4
   module procedure lookup_es_des_k_3d_r8
 end interface

 !> @ingroup sat_vapor_pres_k_mod
 interface lookup_es2_k
   module procedure lookup_es2_k_0d_r4
   module procedure lookup_es2_k_0d_r8
   module procedure lookup_es2_k_1d_r4
   module procedure lookup_es2_k_1d_r8
   module procedure lookup_es2_k_2d_r4
   module procedure lookup_es2_k_2d_r8
   module procedure lookup_es2_k_3d_r4
   module procedure lookup_es2_k_3d_r8
 end interface

 !> @ingroup sat_vapor_pres_k_mod
 interface lookup_des2_k
   module procedure lookup_des2_k_0d_r4
   module procedure lookup_des2_k_0d_r8
   module procedure lookup_des2_k_1d_r4
   module procedure lookup_des2_k_1d_r8
   module procedure lookup_des2_k_2d_r4
   module procedure lookup_des2_k_2d_r8
   module procedure lookup_des2_k_3d_r4
   module procedure lookup_des2_k_3d_r8
 end interface

 !> @ingroup sat_vapor_pres_k_mod
 interface lookup_es2_des2_k
   module procedure lookup_es2_des2_k_0d_r4
   module procedure lookup_es2_des2_k_0d_r8
   module procedure lookup_es2_des2_k_1d_r4
   module procedure lookup_es2_des2_k_1d_r8
   module procedure lookup_es2_des2_k_2d_r4
   module procedure lookup_es2_des2_k_2d_r8
   module procedure lookup_es2_des2_k_3d_r4
   module procedure lookup_es2_des2_k_3d_r8
 end interface

 !> @ingroup sat_vapor_pres_k_mod
 interface lookup_es3_k
   module procedure lookup_es3_k_0d_r4
   module procedure lookup_es3_k_0d_r8
   module procedure lookup_es3_k_1d_r4
   module procedure lookup_es3_k_1d_r8
   module procedure lookup_es3_k_2d_r4
   module procedure lookup_es3_k_2d_r8
   module procedure lookup_es3_k_3d_r4
   module procedure lookup_es3_k_3d_r8
 end interface

 !> @ingroup sat_vapor_pres_k_mod
 interface lookup_des3_k
   module procedure lookup_des3_k_0d_r4
   module procedure lookup_des3_k_0d_r8
   module procedure lookup_des3_k_1d_r4
   module procedure lookup_des3_k_1d_r8
   module procedure lookup_des3_k_2d_r4
   module procedure lookup_des3_k_2d_r8
   module procedure lookup_des3_k_3d_r4
   module procedure lookup_des3_k_3d_r8
 end interface

 !> @ingroup sat_vapor_pres_k_mod
 interface lookup_es3_des3_k
   module procedure lookup_es3_des3_k_0d_r4
   module procedure lookup_es3_des3_k_0d_r8
   module procedure lookup_es3_des3_k_1d_r4
   module procedure lookup_es3_des3_k_1d_r8
   module procedure lookup_es3_des3_k_2d_r4
   module procedure lookup_es3_des3_k_2d_r8
   module procedure lookup_es3_des3_k_3d_r4
   module procedure lookup_es3_des3_k_3d_r8
 end interface

 !> @ingroup sat_vapor_pres_k_mod
 interface compute_qs_k
   module procedure compute_qs_k_0d_r4
   module procedure compute_qs_k_0d_r8
   module procedure compute_qs_k_1d_r4
   module procedure compute_qs_k_1d_r8
   module procedure compute_qs_k_2d_r4
   module procedure compute_qs_k_2d_r8
   module procedure compute_qs_k_3d_r4
   module procedure compute_qs_k_3d_r8
 end interface
 !> @ingroup sat_vapor_pres_k_mod
 interface compute_mrs_k
   module procedure compute_mrs_k_0d_r4
   module procedure compute_mrs_k_0d_r8
   module procedure compute_mrs_k_1d_r4
   module procedure compute_mrs_k_1d_r8
   module procedure compute_mrs_k_2d_r4
   module procedure compute_mrs_k_2d_r8
   module procedure compute_mrs_k_3d_r4
   module procedure compute_mrs_k_3d_r8
 end interface compute_mrs_k

!> @addtogroup sat_vapor_pres_k_mod
!> @{

 real(kind=r8_kind) :: dtres, tepsl, tminl, dtinvl
 integer :: table_siz
 real(kind=r8_kind), dimension(:), allocatable :: TABLE   !  sat vapor pres (es)
 real(kind=r8_kind), dimension(:), allocatable :: DTABLE  !  first derivative of es
 real(kind=r8_kind), dimension(:), allocatable :: D2TABLE ! second derivative of es
 real(kind=r8_kind), dimension(:), allocatable :: TABLE2  !  sat vapor pres (es)
 real(kind=r8_kind), dimension(:), allocatable :: DTABLE2 !  first derivative of es
 real(kind=r8_kind), dimension(:), allocatable :: D2TABLE2 ! second derivative of es
 real(kind=r8_kind), dimension(:), allocatable :: TABLE3  !  sat vapor pres (es)
 real(kind=r8_kind), dimension(:), allocatable :: DTABLE3 !  first derivative of es
 real(kind=r8_kind), dimension(:), allocatable :: D2TABLE3 ! second derivative of es

 logical  :: use_exact_qs
 logical  :: module_is_initialized = .false.

 contains

!#######################################################################
!#######################################################################

#include "sat_vapor_pres_k_r4.fh"
#include "sat_vapor_pres_k_r8.fh"

 end module sat_vapor_pres_k_mod
!> @}
! close documentation grouping
