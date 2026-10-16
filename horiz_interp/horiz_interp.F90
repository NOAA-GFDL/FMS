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
!> @defgroup horiz_interp_mod horiz_interp_mod
!! @ingroup horiz_interp
!! @{
!! @author Zhi Liang, Bruce Wyman
!!
!! @parblock
!! Horiz_interp_mod contains subroutines and derived types to interpolate
!! data from  any logically rectangular grid to any logically rectangular grid with
!! the following interpolation schemes:  conservative in horiz_interp_conserve_mod,
!! bilinear in horiz_interp_bilinear_mod, bicubic in horiz_interp_bicubic_mod, and
!! inverse of square distance weighted in horiz_interp_spherical.mod.
!! @endparblock

module horiz_interp_mod

!-----------------------------------------------------------------------
!
!        Performs spatial interpolation between grids.
!
!-----------------------------------------------------------------------

use fms_mod,                    only: write_version_number, fms_error_handler
use fms_mod,                    only: check_nml_error
use mpp_mod,                    only: mpp_error, FATAL, stdout, stdlog, mpp_min
use mpp_mod,                    only: input_nml_file, WARNING, mpp_pe, mpp_root_pe
use constants_mod,              only: pi
use horiz_interp_type_mod,      only: horiz_interp_type, assignment(=)
use horiz_interp_type_mod,      only: CONSERVE, BILINEAR, SPHERICAL, BICUBIC
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_init, horiz_interp_conserve
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_new, horiz_interp_conserve_del
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_init, horiz_interp_bilinear
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_new, horiz_interp_bilinear_del
use horiz_interp_bilinear_mod,  only: horiz_interp_read_weights_bilinear
use horiz_interp_bicubic_mod,   only: horiz_interp_bicubic_init, horiz_interp_bicubic
use horiz_interp_bicubic_mod,   only: horiz_interp_bicubic_new, horiz_interp_bicubic_del
use horiz_interp_spherical_mod, only: horiz_interp_spherical_init, horiz_interp_spherical
use horiz_interp_spherical_mod, only: horiz_interp_spherical_new, horiz_interp_spherical_del
use platform_mod,               only: r4_kind, r8_kind

 implicit none
 private

!---- interfaces ----

 public   horiz_interp_type, horiz_interp, horiz_interp_new, horiz_interp_del, &
          horiz_interp_init, horiz_interp_end, assignment(=), horiz_interp_read_weights

 !> @parblock
 !! Generic interface to horiz_interp_new procedures to compute interpolation
 !! weights and mapping indices. Following captures the main input arguemnts for each
 !! member subroutine:
 !! horiz_interp_new_1d_r4:
 !!   input and output grids are provided as 1D arrays in 32-bit floating point precision.
 !! horiz_interp_new_1d_r8:
 !!   input and output grids are provided as 1D arrays in 64-bit floating point precision.
 !! horiz_interp_new_1d_src_r4:
 !!   input grid is provided as 1D arrays, output grid as 2D arrays, both in 32-bit precision.
 !! horiz_interp_new_1d_src_r8:
 !!   input grid is provided as 1D arrays, output grid as 2D arrays, both in 64-bit precision.
 !! horiz_interp_new_2d_r4:
 !!   input and output grids are provided as 2D arrays in 32-bit precision.
 !! horiz_interp_new_2d_r8:
 !!   input and output grids are provided as 2D arrays in 64-bit precision.
 !! horiz_interp_new_1d_dst_r4:
 !!   input grid is provided as 2D arrays, output grid as 1D arrays, both in 32-bit precision.
 !! horiz_interp_new_1d_dst_r8:
 !!   input grid is provided as 2D arrays, output grid as 1D arrays, both in 64-bit precision.
 !! @endparblock
 interface horiz_interp_new
    ! Source grid is 1d, destination grid is 1d
    module procedure horiz_interp_new_1d_r4
    module procedure horiz_interp_new_1d_r8
    ! Source grid is 1d, destination grid is 2d
    module procedure horiz_interp_new_1d_src_r4
    module procedure horiz_interp_new_1d_src_r8
    ! Source grid is 2d, destination grid is 2d
    module procedure horiz_interp_new_2d_r4
    module procedure horiz_interp_new_2d_r8
    ! Source grid is 2d, destination grid is 1d
    module procedure horiz_interp_new_1d_dst_r4
    module procedure horiz_interp_new_1d_dst_r8
 end interface

 !> Generic interface to horiz_interp_read_weights_r4 and horiz_interp_read_weights_r8
 !! to read in weight files generated from fregrid for bilinear interpolation.
 !! Calls horiz_interp_read_weights_r4 if lat_out, lon_out, lat_in, and lon_in are
 !! 32-bit floating point representation.  Calls horiz_interp_read_weights_r8 if
 !! the gridpoints are in 64-bit floating point representation.
 interface horiz_interp_read_weights
   module procedure horiz_interp_read_weights_r4
   module procedure horiz_interp_read_weights_r8
 end interface horiz_interp_read_weights

 !> @parblock
 !! Generic interface to interpolate data on input (source) grid to output (target) grid.
 !! Calls horiz_interp_base_2d_r4, horiz_interp_base_2d_r8, horiz_interp_base_3d_r4, or
 !! horiz_interp_base_3d_r8 are if Interp (populated from horiz_interp_new) is the
 !! first argument.  Else, if Interp is not provided, calls horiz_interp_solo_* blackbox methods:
 !! Following captures the main input arguments for each subroutine:
 !! horiz_interp_base_2d_r4:
 !!   2d input data to 2d output data in 32-bit floating point precision.  Takes Interp as first argument.
 !! horiz_interp_base_2d_r8:
 !!   2d input data to 2d output data in 64-bit floating point precision.  Takes Interp as first argument.
 !! horiz_interp_base_3d_r4:
 !!   3d input data to 3d output data in 32-bit floating point precision.  Takes Interp as first argument.
 !! horiz_interp_base_3d_r8:
 !!   3d input data to 3d output data in 64-bit floating point precision.  Takes Interp as first argument.
 !! horiz_interp_solo_1d_r4:
 !!   input and output grids provided as 1D arrays to interpolate 32-bit 2D data.
 !!   Does not take Interp as argument.
 !! horiz_interp_solo_1d_r8:
 !!   input and output grids provided as 1D arrays to interpolate 64-bit 2D data.
 !!   Does not take Interp as argument.
 !! horiz_interp_solo_1d_src_r4:
 !!   input grid provided as 32-bit 1D arrays, output as 2D arrays to interpolate 32-bit 2D data.
 !!   Does not take Inerp as argument.
 !! horiz_interp_solo_1d_src_r8:
 !!   input grid provided as 64-bit 1D arrays, output as 2D arrays to interpolate  32-bit 2D data.
 !!   Does not take Interp as argument.
 !! horiz_interp_solo_2d_r4:
 !!   input and output grids provided as 32-bit 2D arrays to interpolate 32-bit 2d data.
 !!   Does not take Interp as argument.
 !! horiz_interp_solo_2d_r8:
 !!   input and output grids provided as 64-bit 2D arrays to interpolate 64-bit 2d data.
 !!   Does not take Interp as argument.
 !! horiz_interp_solo_1d_dst_r4:
 !!   input grid provided as 32-bit 2D arrays, input as 1D arrays to interpolate 32-bit 2D data.
 !!   Does not take Interp as argument.
 !! horiz_interp_solo_1d_dst_r8
 !!   input grid provided as 64-bit 2D arrays, input as 1D arrays to interpolate 64-bit 2D data.
 !!   Does not take Interp as argument.
 !! @endparblock
 interface horiz_interp
    module procedure horiz_interp_base_2d_r4
    module procedure horiz_interp_base_2d_r8
    module procedure horiz_interp_base_3d_r4
    module procedure horiz_interp_base_3d_r8
    module procedure horiz_interp_solo_1d_r4
    module procedure horiz_interp_solo_1d_r8
    module procedure horiz_interp_solo_1d_src_r4
    module procedure horiz_interp_solo_1d_src_r8
    module procedure horiz_interp_solo_2d_r4
    module procedure horiz_interp_solo_2d_r8
    module procedure horiz_interp_solo_1d_dst_r4
    module procedure horiz_interp_solo_1d_dst_r8
    module procedure horiz_interp_solo_old_r4
    module procedure horiz_interp_solo_old_r8
 end interface


 !> Internally used subroutine to determine if the grid is a lat-lon grid.
 !! Calls is_lat_lon_r4 if grids are in 32-bit floating point precision.
 !! Calls is_lat_lon_r8 if grids are in 64-bit floating point precision.
 interface is_lat_lon
    module procedure is_lat_lon_r4
    module procedure is_lat_lon_r8
  end interface is_lat_lon

 !> Old interfaces, not recommended.
 interface horiz_interp_solo_1d
   module procedure horiz_interp_solo_1d_r4
   module procedure horiz_interp_solo_1d_r8
 end interface horiz_interp_solo_1d


 logical :: reproduce_siena = .false.
   !< is a namelist flag to reproduces siena results if set to .true.  Else, defaults
   !! to false to decrease truncation error in function poly_area in file mosaic_util.c.
   !! The truncation error of second order conservative remapping might be big for high resolution grid.
   !! This namelist flag will be removed soon.

 namelist /horiz_interp_nml/ reproduce_siena

!-----------------------------------------------------------------------
! Include variable "version" to be written to log file.
#include<file_version.h>
 logical            :: module_is_initialized = .FALSE.
    !< is an internally used flag to prevent module re-initialization
!-----------------------------------------------------------------------

contains

!#######################################################################

  !> @parblock
  !! Initializes horiz_interp_mod and writes version number to the logfile
  !! @endparblock
  subroutine horiz_interp_init
  integer :: iunit, ierr, io

  if(module_is_initialized) return
  call write_version_number("HORIZ_INTERP_MOD", version)

  read (input_nml_file, horiz_interp_nml, iostat=io)
  ierr = check_nml_error(io,'horiz_interp_nml')
  if (mpp_pe() == mpp_root_pe() ) then
     iunit = stdlog()
     write (iunit, nml=horiz_interp_nml)
  endif

  if (reproduce_siena) then
    call mpp_error(FATAL, "horiz_interp_mod: You have overridden the default value of " // &
       "reproduce_siena and set it to .true. in horiz_interp_nml. This was a temporary workaround to " // &
       "allow for consistency in continuing experiments and is no longer supported. " // &
       "Please remove this namelist.")
  endif

  call horiz_interp_conserve_init
  call horiz_interp_bilinear_init
  call horiz_interp_bicubic_init
  call horiz_interp_spherical_init

  module_is_initialized = .true.

  end subroutine horiz_interp_init

  !> @parblock
  !! Deallocates memory in Interp holding the weights and mapping indices for the interpolation
  !! method in Interp%interp_method.
  !! @endparblock
  subroutine horiz_interp_del ( Interp )

   type (horiz_interp_type), intent(inout) :: Interp
     !< will be reset after arrays holding interpolation weights and mapping indices are
     !! deallocated and %is_allocated is set to .false.

!-----------------------------------------------------------------------
!  releases space used by horiz_interp_type variables
!  must be called before re-initializing the same variable
!-----------------------------------------------------------------------
   select case(Interp % interp_method)
   case (CONSERVE)
      call horiz_interp_conserve_del(Interp )
   case (BILINEAR)
      call horiz_interp_bilinear_del(Interp )
   case (BICUBIC)
      call horiz_interp_bicubic_del(Interp )
   case (SPHERICAL)
      call horiz_interp_spherical_del(Interp )
   end select

   Interp%I_am_initialized = .false.
!-----------------------------------------------------------------------

 end subroutine horiz_interp_del

 !#####################################################################

 !> @parblock
 !! Dummy routine
 !! @endparblock
 subroutine horiz_interp_end
 return
 end subroutine horiz_interp_end

#include "horiz_interp_r4.fh"
#include "horiz_interp_r8.fh"

end module horiz_interp_mod
!> @}
! close documentation grouping
