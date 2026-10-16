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
!> @defgroup horiz_interp_bilinear_mod horiz_interp_bilinear_mod
!! @ingroup horiz_interp
!! @{
!! @author Zhi Liang <Zhi.Liang@noaa.gov>
!! @parblock
!! Horiz_interp_bilinear_mod contains methods called from horiz_interp_mod to
!! interpolate data on a regular rectangular grid to a rectangular/tripolar grid.
!! Users are recommened to use the top-level module horiz_interp_mod with "interp"
!! set to "bilinear" for bilinear interpolation.
!! @endparblock

module horiz_interp_bilinear_mod

  use mpp_mod,               only: mpp_error, FATAL, stdout, mpp_pe, mpp_root_pe
  use fms_mod,               only: write_version_number
  use constants_mod,         only: PI
  use horiz_interp_type_mod, only: horiz_interp_type, stats, BILINEAR
  use platform_mod,          only: r4_kind, r8_kind
  use axis_utils2_mod,       only: nearest_index
  use fms2_io_mod,           only: open_file, close_file, read_data, FmsNetcdfFile_t, get_dimension_size
  use fms_string_utils_mod,  only: string

  implicit none
  private


  public :: horiz_interp_bilinear_new, horiz_interp_bilinear, horiz_interp_bilinear_del
  public :: horiz_interp_bilinear_init, horiz_interp_read_weights_bilinear

  !> Generic interface to compute interpolation weights and mapping indices.
  !! Member subroutines and their main arguments are:
  !! horiz_interp_bilinear_new_1d_r4:
  !!   input grids provided as 1D arrays, output grids as 2D arrays, both in 32-bit precision.
  !! horiz_interp_bilinear_new_1d_r8:
  !!   input grids provided as 1D arrays, output grids as 2D arrays, both in 64-bit precision.
  !! horiz_interp_bilinear_new_2d_r4:
  !!   input and output grids provided as 2D arrays in 32-bit precision.
  !! horiz_interp_bilinear__new_2d_r8:
  !!   input and output grids provided as 2D arrays in 64-bit precision.
  interface horiz_interp_bilinear_new
    module procedure horiz_interp_bilinear_new_1d_r4
    module procedure horiz_interp_bilinear_new_1d_r8
    module procedure horiz_interp_bilinear_new_2d_r4
    module procedure horiz_interp_bilinear_new_2d_r8
  end interface

  !> Generic interface to populate Interp from a weight file generated from fregrid.
  !! Calls horiz_interp_read_weights_bilinear_r4 if grid arguments are provided as 32-bit precision.
  !! Calls horiz_interp_read_weights_bilienar_r8 if grid arguments are provided as 64-bit precision.
  interface horiz_interp_read_weights_bilinear
    module procedure horiz_interp_read_weights_bilinear_r4
    module procedure horiz_interp_read_weights_bilinear_r8
  end interface

  !> Generic interface to interpolate data from a source grid to target grid using the weights
  !! stored in Interp.  Calls horiz_interp_bilinear_r4 if input and output data are 2D arrays
  !! in 32-bit precision.  Calls horiz_interp_bilinear_r8 if input and output data are 2D arrays
  !! in 64-bit precision.
  interface horiz_interp_bilinear
    module procedure horiz_interp_bilinear_r4
    module procedure horiz_interp_bilinear_r8
  end interface

  real(r8_kind), parameter :: epsln=1.e-10_r8_kind !< is a really small number
  real(r4_kind), parameter :: epsln_r4=1.e-4_r4_kind !< is 0.0001
  integer, parameter :: DUMMY = -999 !< is -999

  !> Generic interface used internally when finding nearest neighboring input cells
  !! for an output cell.  Calls intersect_r4 when grid points are in 32-bit precision.
  !! Calls intersect_r8 when grid points are in 64-bit precision.
  interface intersect
    module procedure intersect_r4
    module procedure intersect_r8
  end interface

#include<file_version.h>
  logical :: module_is_initialized = .FALSE. !< is a flag to prevent re-initialization

contains

  !> @parblock
  !! Initializes horiz_interp_bilinear_mod.  Called from horiz_interp_init in horiz_interp_mod.
  !! @endparblock
  subroutine horiz_interp_bilinear_init

    if(module_is_initialized) return
    call write_version_number("HORIZ_INTERP_BILINEAR_MOD", version)
    module_is_initialized = .true.

  end subroutine horiz_interp_bilinear_init

  !> @brief Deallocates memory used by "horiz_interp_type" variables.
  !!
  !> Must be called before reinitializing with horiz_interp_bilinear_new.

  !> @parblock
  !! Deallocates arrays holding bilinear interpolation weights and mapping indices in Interp.
  !! Resets %is_allocated to .false.  Called from horiz_interp_del in horiz_interp_mod.
  !! @endparblock
  subroutine horiz_interp_bilinear_del( Interp )

    type (horiz_interp_type), intent(inout) :: Interp !< will be reset with deallocated memory.

    if( Interp%horizInterpReals4_type%is_allocated) then
      if(allocated(Interp%horizInterpReals4_type%wti))   deallocate(Interp%horizInterpReals4_type%wti)
      if(allocated(Interp%horizInterpReals4_type%wtj))   deallocate(Interp%horizInterpReals4_type%wtj)
    else if (Interp%horizInterpReals8_type%is_allocated) then
      if(allocated(Interp%horizInterpReals8_type%wti))   deallocate(Interp%horizInterpReals8_type%wti)
      if(allocated(Interp%horizInterpReals8_type%wtj))   deallocate(Interp%horizInterpReals8_type%wtj)
    endif
    if(allocated(Interp%i_lon)) deallocate(Interp%i_lon)
    if(allocated(Interp%j_lat)) deallocate(Interp%j_lat)

    Interp%horizInterpReals4_type%is_allocated = .false.
    Interp%horizInterpReals8_type%is_allocated = .false.

  end subroutine horiz_interp_bilinear_del

#include "horiz_interp_bilinear_r4.fh"
#include "horiz_interp_bilinear_r8.fh"

end module horiz_interp_bilinear_mod
!> @}
