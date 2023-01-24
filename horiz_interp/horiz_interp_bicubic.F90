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
!> @defgroup horiz_interp_bicubic_mod horiz_interp_bicubic_mod
!> @ingroup horiz_interp
!> @brief Delivers methods for bicubic interpolation from a coarse regular grid
!! on a fine regular grid
!!
!> This module delivers methods for bicubic interpolation from a
!! coarse regular grid on a fine regular grid.
!! Subroutines
!!
!! - @ref bcuint
!! - @ref bcucof
!!
!! are methods taken from
!!
!!       W. H. Press, S. A. Teukolski, W. T. Vetterling and B. P. Flannery,
!!       Numerical Recipies in FORTRAN, The Art of Scientific Computing.
!!       Cambridge University Press, 1992
!!
!! written by
!!       martin.schmidt@io-warnemuende.de (2004)
!! revised by
!!       martin.schmidt@io-warnemuende.de (2004)
!!
!! Version 1.0.0.2005-07-06
!! The module is thought to interact with MOM-4.
!! Alle benotigten Felder werden extern von MOM verwaltet, da sie
!! nicht fur alle interpolierten Daten die gleiche Dimension haben mussen.
module horiz_interp_bicubic_mod

  use mpp_mod,               only: mpp_error, FATAL, stdout, mpp_pe, mpp_root_pe
  use fms_mod,               only: write_version_number
  use horiz_interp_type_mod, only: horiz_interp_type
  use constants_mod,         only: PI
  use platform_mod,          only: r4_kind, r8_kind


 implicit none

   private

   public  :: horiz_interp_bicubic, horiz_interp_bicubic_new, horiz_interp_bicubic_del, fill_xy
   public  :: horiz_interp_bicubic_init

  !> Creates a new @ref horiz_interp_type for bicubic interpolation.
  !> @ingroup horiz_interp_bicubic_mod
  interface horiz_interp_bicubic_new
    module procedure horiz_interp_bicubic_new_1d_r8
    module procedure horiz_interp_bicubic_new_1d_s_r8
    module procedure horiz_interp_bicubic_new_1d_r4
    module procedure horiz_interp_bicubic_new_1d_s_r4
  end interface

  interface horiz_interp_bicubic
    module procedure horiz_interp_bicubic_r4
    module procedure horiz_interp_bicubic_r8
  end interface

!> @addtogroup horiz_interp_bicubic_mod
!> @{

! Include variable "version" to be written to log file.
#include<file_version.h>
   logical            :: module_is_initialized = .FALSE.
   integer            :: verbose_bicubic = 0

!     Grid variables
!     xc, yc : co-ordinates of the coarse grid
!     xf, yf : co-ordinates of the fine grid
!     fc     : variable to be interpolated at the coarse grid
!     dfc_x  : x-derivative of fc at the coarse grid
!     dfc_y  : y-derivative of fc at the coarse grid
!     dfc_xy : x-y-derivative of fc at the coarse grid
!     ff     : variable to be interpolated at the fine grid
!     dff_x  : x-derivative of fc at the fine grid
!     dff_y  : y-derivative of fc at the fine grid
!     dff_xy : x-y-derivative of fc at the fine grid


   real(r8_kind)               :: tpi

   interface fill_xy
      module procedure fill_xy_r4
      module procedure fill_xy_r8
   end interface


   contains

  !> @brief Initializes module and writes version number to logfile.out
  subroutine horiz_interp_bicubic_init

     if(module_is_initialized) return
     call write_version_number("HORIZ_INTERP_BICUBIC_MOD", version)
     module_is_initialized = .true.
     tpi = 2.0*PI

  end subroutine horiz_interp_bicubic_init

  !#######################################################################


  subroutine horiz_interp_bicubic_del( Interp )
    type(horiz_interp_type), intent(inout) :: Interp

    if(allocated(Interp%kind8_reals)) then
      if(associated(Interp%kind8_reals%rat_x))  deallocate ( Interp%kind8_reals%rat_x )
      if(associated(Interp%kind8_reals%rat_y))  deallocate ( Interp%kind8_reals%rat_y )
      if(associated(Interp%kind8_reals%lon_in)) deallocate ( Interp%kind8_reals%lon_in )
      if(associated(Interp%kind8_reals%lat_in)) deallocate ( Interp%kind8_reals%lat_in )
      if(associated(Interp%kind8_reals%wti))    deallocate ( Interp%kind8_reals%wti )
      deallocate(Interp%kind8_reals)
    else if(allocated(Interp%kind4_reals)) then
      if(associated(Interp%kind4_reals%rat_x))  deallocate ( Interp%kind4_reals%rat_x )
      if(associated(Interp%kind4_reals%rat_y))  deallocate ( Interp%kind4_reals%rat_y )
      if(associated(Interp%kind4_reals%lon_in)) deallocate ( Interp%kind4_reals%lon_in )
      if(associated(Interp%kind4_reals%lat_in)) deallocate ( Interp%kind4_reals%lat_in )
      if(associated(Interp%kind4_reals%wti))    deallocate ( Interp%kind4_reals%wti )
      deallocate(Interp%kind4_reals)
    endif

    if(allocated(Interp%rat_x))  deallocate ( Interp%rat_x )
    if(allocated(Interp%rat_y))  deallocate ( Interp%rat_y )
    if(allocated(Interp%lon_in)) deallocate ( Interp%lon_in )
    if(allocated(Interp%lat_in)) deallocate ( Interp%lat_in )
    if(allocated(Interp%i_lon))  deallocate ( Interp%i_lon )
    if(allocated(Interp%j_lat))  deallocate ( Interp%j_lat )
    if(allocated(Interp%wti))    deallocate ( Interp%wti )

  end subroutine horiz_interp_bicubic_del

#include "horiz_interp_bicubic_r4.fh"
#include "horiz_interp_bicubic_r8.fh"

end module horiz_interp_bicubic_mod
!> @}
! close documentation
