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
!> @defgroup horiz_interp_bicubic_mod horiz_interp_bicubic_mod
!! @ingroup horiz_interp
!! @{
!!
!! @author martin.schmidt@io-warnemuende.de (2004)
!!
!! @parblock
!! Horiz_interp_bicubic_mod contains methods called from horiz_interp_mod for bicubic interpolation
!! from a coarse regular grid on a fine regular grids.  Users are recommended  to use the top-level
!! module horiz_interp_mod with "interp" set to "bicubic" for bicubic interpolation.
!! Internally used subroutines bcuint and bcuof, used internally, have been implemented following
!!
!!       W. H. Press, S. A. Teukolski, W. T. Vetterling and B. P. Flannery,
!!       Numerical Recipies in FORTRAN, The Art of Scientific Computing.
!!       Cambridge University Press, 1992
!!
!! All public interfaces and subroutines in this module can be accessed from public interfaces
!! and subroutines in horiz_interp_mod.  For simplicity, users are recommended to call
!! all procedures throug horiz_interp_mod.
!! Note from the author:  The module is thought to interact with MOM-4.
!! Alle benotigten Felder werden extern von MOM verwaltet, da sie
!! nicht fur alle interpolierten Daten die gleiche Dimension haben mussen.
!! @endparblock
module horiz_interp_bicubic_mod

  use mpp_mod,               only: mpp_error, FATAL, stdout, mpp_pe, mpp_root_pe
  use fms_mod,               only: write_version_number
  use horiz_interp_type_mod, only: horiz_interp_type, BICUBIC
  use constants_mod,         only: PI
  use platform_mod,          only: r4_kind, r8_kind


 implicit none

   private

   public  :: horiz_interp_bicubic, horiz_interp_bicubic_new, horiz_interp_bicubic_del, fill_xy
   public  :: horiz_interp_bicubic_init

   !> Generic interface to horiz_interp_bicubic_new_* to compute interpolation weights and,
   !! mapping indices. The member subroutines and their main arguments are:
   !! horiz_interp_bicubic_new_1d_r8:
   !!   input and output grids provided as 1D arrays in 64-bit precision.
   !! horiz_interp_bicubic_new_1d_r4:
   !!   input and output grids provided as 1D arrays in 32-bit precision.
   !! horiz_interp_bicubic_new_1d_s_r8:
   !!   input grid provided as 1D arrays, output grids as 2D arrays, both in 64-bit precision.
   !! horiz_interp_bicubic_new_1d_s_r4:
   !!   input grid provided as 1D arrays, output grids as 2D arrays, both in 32-bit precision.
   interface horiz_interp_bicubic_new
     module procedure horiz_interp_bicubic_new_1d_r8
     module procedure horiz_interp_bicubic_new_1d_s_r8
     module procedure horiz_interp_bicubic_new_1d_r4
     module procedure horiz_interp_bicubic_new_1d_s_r4
   end interface horiz_interp_bicubic_new

   !> Generic interface to horiz_interp_bicubic_* to interpolate data using the weights
   !! stored in Interp.  Calls horiz_interp_bicubic_r4 if input and output data are 2D arrays
   !! in 32-bit precision.  Calls horiz_interp_bicubic_r8 if input and output data are 2D
   !! arrays in 64-bit precision.
   interface horiz_interp_bicubic
     module procedure horiz_interp_bicubic_r4
     module procedure horiz_interp_bicubic_r8
   end interface horiz_interp_bicubic

!> Include variable "version" to be written to log file.
#include<file_version.h>
   logical :: module_is_initialized = .FALSE.  !< is a flag to prevent module re-initialization.
   integer :: verbose_bicubic = 0 !< sets the default verbosity level.

   real(r8_kind) :: tpi  !< is 2*PI.

   !> Unused interface.
   interface fill_xy
      module procedure fill_xy_r4
      module procedure fill_xy_r8
   end interface

   !> Generic interface used internally when interpolating data from source grid to
   !! target grid.  Calls bcuint_r4 when input and output data are in 32-bit precision.
   !! Calls bcuint_r8 when input and output data are in 64-bit precision.
   interface bcuint
      module procedure bcuint_r4
      module procedure bcuint_r8
   end interface

   !> Generic interface used internally in bcuint to compute bicubic coefficients when
   !! interpolating data from source grid to target grid.  Calls bcuof_r4 when all arguments
   !! are in 32-bit precision.  Calls bcuof_r8 when all arguemnts are in 64-bit precision.
   interface bcucof
      module procedure bcucof_r4
      module procedure bcucof_r8
   end interface

   !> Generic interface used internally in finding neighbors
   interface indl
      module procedure indl_r4
      module procedure indl_r8
   end interface

   !> Generic interface used internally in finding neighbors
   interface indu
      module procedure indu_r4
      module procedure indu_r8
   end interface

 contains

   !> @parblock
   !! Initializes horiz_interp_bicubic_mod including setting tpi.
   !! Called from horiz_interp_init in horiz_interp_mod.
   !! @endparblock
   subroutine horiz_interp_bicubic_init

     if(module_is_initialized) return
     call write_version_number("HORIZ_INTERP_BICUBIC_MOD", version)
     module_is_initialized = .true.
     tpi = real(2.0_r8_kind*PI, R8_KIND)

  end subroutine horiz_interp_bicubic_init

  !> @parblock
  !! Deallocates arrays holding bicubic interpolation weights and mapping indices in Interp.
  !! Resets %is_allocated to .false.  Called from horiz_interp_del in horiz_interp_mod.
  !! @endparblock
  subroutine horiz_interp_bicubic_del( Interp )
    type(horiz_interp_type), intent(inout) :: Interp !< will be reset with deallocated memory.

    if(Interp%horizInterpReals8_type%is_allocated) then
      if(allocated(Interp%horizInterpReals8_type%rat_x))  deallocate ( Interp%horizInterpReals8_type%rat_x )
      if(allocated(Interp%horizInterpReals8_type%rat_y))  deallocate ( Interp%horizInterpReals8_type%rat_y )
      if(allocated(Interp%horizInterpReals8_type%lon_in)) deallocate ( Interp%horizInterpReals8_type%lon_in )
      if(allocated(Interp%horizInterpReals8_type%lat_in)) deallocate ( Interp%horizInterpReals8_type%lat_in )
      if(allocated(Interp%horizInterpReals8_type%wti))    deallocate ( Interp%horizInterpReals8_type%wti )
    else if(Interp%horizInterpReals4_type%is_allocated) then
      if(allocated(Interp%horizInterpReals4_type%rat_x))  deallocate ( Interp%horizInterpReals4_type%rat_x )
      if(allocated(Interp%horizInterpReals4_type%rat_y))  deallocate ( Interp%horizInterpReals4_type%rat_y )
      if(allocated(Interp%horizInterpReals4_type%lon_in)) deallocate ( Interp%horizInterpReals4_type%lon_in )
      if(allocated(Interp%horizInterpReals4_type%lat_in)) deallocate ( Interp%horizInterpReals4_type%lat_in )
      if(allocated(Interp%horizInterpReals4_type%wti))    deallocate ( Interp%horizInterpReals4_type%wti )
    endif
    if( allocated(Interp%i_lon) ) deallocate( Interp%i_lon )
    if( allocated(Interp%j_lat) ) deallocate( Interp%j_lat )

    Interp%horizInterpReals8_type%is_allocated = .false.
    Interp%horizInterpReals4_type%is_allocated = .false.

  end subroutine horiz_interp_bicubic_del

#include "horiz_interp_bicubic_r4.fh"
#include "horiz_interp_bicubic_r8.fh"

end module horiz_interp_bicubic_mod
!> @}
! close documentation
