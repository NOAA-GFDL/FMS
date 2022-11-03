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
!> @defgroup horiz_interp_bilinear_mod horiz_interp_bilinear_mod
!> @ingroup horiz_interp
!> @brief Performs spatial interpolation between grids using bilinear interpolation
!!
!> @author Zhi Liang <Zhi.Liang@noaa.gov>
!> This module can interpolate data from regular rectangular grid
!! to rectangular/tripolar grid. The interpolation scheme is bilinear interpolation.
!! There is an optional mask field for missing input data.
!! An optional output mask field may be used in conjunction with
!! the input mask to show where output data exists.

module horiz_interp_bilinear_mod

  use mpp_mod,               only: mpp_error, FATAL, stdout, mpp_pe, mpp_root_pe
  use fms_mod,               only: write_version_number
  use constants_mod,         only: PI
  use horiz_interp_type_mod, only: horiz_interp_type, stats

  implicit none
  private


  public :: horiz_interp_bilinear_new, horiz_interp_bilinear, horiz_interp_bilinear_del
  public :: horiz_interp_bilinear_init

  !> Creates a @ref horiz_interp_type for bilinear interpolation.
  !> @ingroup horiz_interp_bilinear_mod
  interface horiz_interp_bilinear_new
    module procedure horiz_interp_bilinear_new_1d_r4
    module procedure horiz_interp_bilinear_new_1d_r8
    module procedure horiz_interp_bilinear_new_2d_r4
    module procedure horiz_interp_bilinear_new_2d_r8
  end interface

  interface horiz_interp_bilinear
    module procedure horiz_interp_bilinear_r4
    module procedure horiz_interp_bilinear_r8
  end interface

!> @addtogroup horiz_interp_bilinear_mod
!> @{

  real, parameter :: epsln=1.e-10
  integer, parameter :: DUMMY = -999

  !-----------------------------------------------------------------------
! Include variable "version" to be written to log file.
#include<file_version.h>
  logical            :: module_is_initialized = .FALSE.

contains

  !> Initialize this module and writes version number to logfile.
  subroutine horiz_interp_bilinear_init

    if(module_is_initialized) return
    call write_version_number("HORIZ_INTERP_BILINEAR_MOD", version)
    module_is_initialized = .true.

  end subroutine horiz_interp_bilinear_init

  !> @brief Deallocates memory used by "horiz_interp_type" variables.
  !!
  !> Must be called before reinitializing with horiz_interp_bilinear_new.
  subroutine horiz_interp_bilinear_del( Interp )

    type (horiz_interp_type), intent(inout) :: Interp!< A derived-type variable returned by previous
                                   !! call to horiz_interp_bilinear_new. The input variable must
                                   !! have allocated arrays. The returned variable will contain
                                   !! deallocated arrays

    if( allocated(Interp%kind4_reals)) then
      if(associated(Interp%kind4_reals%wti))   deallocate(Interp%kind4_reals%wti)
      if(associated(Interp%kind4_reals%wtj))   deallocate(Interp%kind4_reals%wtj)
    else if (allocated(Interp%kind8_reals)) then
      if(associated(Interp%kind8_reals%wti))   deallocate(Interp%kind8_reals%wti)
      if(associated(Interp%kind8_reals%wtj))   deallocate(Interp%kind8_reals%wtj)
    endif
    if(associated(Interp%i_lon)) deallocate(Interp%i_lon)
    if(associated(Interp%j_lat)) deallocate(Interp%j_lat)

  end subroutine horiz_interp_bilinear_del

#undef FMS_HI_KIND
#define FMS_HI_KIND 4
#undef HI_KIND_TYPE
#define HI_KIND_TYPE kind4_reals 
#undef BILINEAR_NEW_1D
#define BILINEAR_NEW_1D horiz_interp_bilinear_new_1d_r4
#undef BILINEAR_NEW_2D
#define BILINEAR_NEW_2D horiz_interp_bilinear_new_2d_r4
#undef DO_BILINEAR
#define DO_BILINEAR horiz_interp_bilinear_r4
#undef FIND_NEIGHBOR
#define FIND_NEIGHBOR find_neighbor_r4
#undef FIND_NEIGHBOR_NEW
#define FIND_NEIGHBOR_NEW find_neighbor_new_r4
#undef inside_polygon
#define inside_polygon inside_polygon_r4
#undef intersect
#define intersect intersect_r4
#undef indp
#define indp indp_r4
#include <horiz_interp_bilinear.inc>
#undef FMS_HI_KIND
#define FMS_HI_KIND 8
#undef HI_KIND_TYPE
#define HI_KIND_TYPE kind8_reals 
#undef BILINEAR_NEW_1D
#define BILINEAR_NEW_1D horiz_interp_bilinear_new_1d_r8
#undef BILINEAR_NEW_2D
#define BILINEAR_NEW_2D horiz_interp_bilinear_new_2d_r8
#undef DO_BILINEAR
#define DO_BILINEAR horiz_interp_bilinear_r8
#undef FIND_NEIGHBOR
#define FIND_NEIGHBOR find_neighbor_r8
#undef FIND_NEIGHBOR_NEW
#define FIND_NEIGHBOR_NEW find_neighbor_new_r8
#undef inside_polygon 
#define inside_polygon inside_polygon_r8
#undef intersect
#define intersect intersect_r8
#undef indp
#define indp indp_r8
#include <horiz_interp_bilinear.inc>


end module horiz_interp_bilinear_mod
!> @}
! close documentation grouping
