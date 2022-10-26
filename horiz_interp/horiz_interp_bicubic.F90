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
  use horiz_interp_type_mod, only: horiz_interp_type, horiz_interp_r4_type, horiz_interp_r8_type
  use constants_mod,         only: PI


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


   real               :: tpi

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

    class(horiz_interp_type), intent(inout) :: Interp
    select type(Interp)
      type is(horiz_interp_r4_type)
        if(associated(Interp%rat_x))  deallocate ( Interp%rat_x )
        if(associated(Interp%rat_y))  deallocate ( Interp%rat_y )
        if(associated(Interp%lon_in)) deallocate ( Interp%lon_in )
        if(associated(Interp%lat_in)) deallocate ( Interp%lat_in )
        if(associated(Interp%i_lon))  deallocate ( Interp%i_lon )
        if(associated(Interp%j_lat))  deallocate ( Interp%j_lat )
        if(associated(Interp%wti))    deallocate ( Interp%wti )
      type is(horiz_interp_r8_type)
        if(associated(Interp%rat_x))  deallocate ( Interp%rat_x )
        if(associated(Interp%rat_y))  deallocate ( Interp%rat_y )
        if(associated(Interp%lon_in)) deallocate ( Interp%lon_in )
        if(associated(Interp%lat_in)) deallocate ( Interp%lat_in )
        if(associated(Interp%i_lon))  deallocate ( Interp%i_lon )
        if(associated(Interp%j_lat))  deallocate ( Interp%j_lat )
        if(associated(Interp%wti))    deallocate ( Interp%wti )
    end select

  end subroutine horiz_interp_bicubic_del

  
#undef FMS_HI_KIND
#define FMS_HI_KIND 4
#undef BICUBIC_NEW_1D_S 
#define BICUBIC_NEW_1D_S horiz_interp_bicubic_new_1d_s_r4 
#undef BICUBIC_NEW_1D
#define BICUBIC_NEW_1D horiz_interp_bicubic_new_1d_r4 
#undef BICUBIC_NEW
#define BICUBIC_NEW horiz_interp_bicubic_r4
#undef BCUINT
#define BCUINT bcuint_r4
#undef BCUCOF
#define BCUCOF bcucof_r4 
#undef HI_KIND_TYPE
#define HI_KIND_TYPE horiz_interp_r4_type
#undef INDL
#define INDL indl_r4
#undef INDU
#define INDU indu_r4
#undef FILL_XY
#define FILL_XY fill_xy_r4
#include <horiz_interp_bicubic.inc> 

#undef FMS_HI_KIND
#define FMS_HI_KIND 8
#undef BICUBIC_NEW_1D_S 
#define BICUBIC_NEW_1D_S horiz_interp_bicubic_new_1d_s_r8 
#undef BICUBIC_NEW_1D
#define BICUBIC_NEW_1D horiz_interp_bicubic_new_1d_r8 
#undef BICUBIC_NEW
#define BICUBIC_NEW horiz_interp_bicubic_r8
#undef BCUINT
#define BCUINT bcuint_r8
#undef BCUCOF
#define BCUCOF bcucof_r8 
#undef HI_KIND_TYPE
#define HI_KIND_TYPE horiz_interp_r8_type
#undef INDL
#define INDL indl_r8
#undef INDU
#define INDU indu_r8
#undef FILL_XY
#define FILL_XY fill_xy_r8
#include <horiz_interp_bicubic.inc> 

end module horiz_interp_bicubic_mod
!> @}
! close documentation
