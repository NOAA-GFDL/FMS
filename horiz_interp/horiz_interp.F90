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
!> @defgroup horiz_interp_mod horiz_interp_mod
!> @ingroup horiz_interp
!> @brief Performs spatial interpolation between grids.
!!
!> @author Zhi Liang, Bruce Wyman
!!
!! This module can interpolate data from any logically rectangular grid
!! to any logically rectangular grid. Four interpolation schems are used here:
!! conservative, bilinear, bicubic and inverse of square distance weighted.
!! The four interpolation schemes are implemented seperately in
!! horiz_interp_conserver_mod, horiz_interp_blinear_mod, horiz_interp_bicubic_mod
!! and horiz_interp_spherical_mod. bicubic interpolation requires the source grid
!! is regular lon/lat grid. User can choose the interpolation method in the
!! public interface horiz_interp_new through optional argument interp_method,
!! with acceptable value "conservative", "bilinear", "bicubic" and "spherical".
!! The default value is "conservative". There is an optional mask field for
!! missing input data. An optional output mask field may be used in conjunction with
!! the input mask to show where output data exists.

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
use horiz_interp_type_mod,      only: CONSERVE, BILINEAR, SPHERICA, BICUBIC
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_init, horiz_interp_conserve
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_new, horiz_interp_conserve_del
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_init, horiz_interp_bilinear
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_new, horiz_interp_bilinear_del
use horiz_interp_bicubic_mod,   only: horiz_interp_bicubic_init, horiz_interp_bicubic
use horiz_interp_bicubic_mod,   only: horiz_interp_bicubic_new, horiz_interp_bicubic_del
use horiz_interp_spherical_mod, only: horiz_interp_spherical_init, horiz_interp_spherical
use horiz_interp_spherical_mod, only: horiz_interp_spherical_new, horiz_interp_spherical_del
use platform_mod,               only: r4_kind, r8_kind

 implicit none
 private

!---- interfaces ----

 public   horiz_interp_type, horiz_interp, horiz_interp_new, horiz_interp_del, &
          horiz_interp_init, horiz_interp_end, assignment(=)

!> Allocates space and initializes a derived-type variable
!! that contains pre-computed interpolation indices and weights.
!!
!> Allocates space and initializes a derived-type variable
!! that contains pre-computed interpolation indices and weights
!! for improved performance of multiple interpolations between
!! the same grids. This routine does not need to be called if you
!! are doing a single grid-to-grid interpolation.
!!
!! @param lon_in
!!      Longitude (in radians) for source data grid. You can pass 1-D lon_in to
!!      represent the geographical longitude of regular lon/lat grid, or just
!!      pass geographical longitude(lon_in is 2-D). The grid location may be
!!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!!
!! @param lat_in
!!      Latitude (in radians) for source data grid. You can pass 1-D lat_in to
!!      represent the geographical latitude of regular lon/lat grid, or just
!!      pass geographical latitude(lat_in is 2-D). The grid location may be
!!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!!
!! @param lon_out
!!      Longitude (in radians) for destination data grid. You can pass 1-D lon_out to
!!      represent the geographical longitude of regular lon/lat grid, or just
!!      pass geographical longitude(lon_out is 2-D). The grid location may be
!!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!!
!! @param lat_out
!!      Latitude (in radians) for destination data grid. You can pass 1-D lat_out to
!!      represent the geographical latitude of regular lon/lat grid, or just
!!      pass geographical latitude(lat_out is 2-D). The grid location may be
!!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!!
!! @param verbose
!!      Integer flag that controls the amount of printed output.
!!      verbose = 0, no output; = 1, min,max,means; = 2, still more
!!
!! @param interp_method
!!      interpolation method, = "conservative", using conservation scheme,
!!      = "bilinear", using bilinear interpolation, = "spherical",using spherical regrid.
!!      = "bicubic", using bicubic interpolation. The default value is "convervative".
!!
!! @param src_modulo
!!      Indicate the source data grid is cyclic or not.
!!
!! @param grid_at_center
!!      Indicate the data is on the center of grid box or the edge of grid box.
!!      When true, the data is on the center of grid box. default vaule is false.
!!      This option is only available when interp_method = "bilinear" or "bicubic".
!!
!! @param Interp
!!      A derived-type variable containing indices and weights used for subsequent
!!      interpolations. To reinitialize this variable for a different grid-to-grid
!!      interpolation you must first use the "horiz_interp_del" interface.
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


!> Subroutine for performing the horizontal interpolation between two grids.
!!
!> Subroutine for performing the horizontal interpolation between
!! two grids. There are two forms of this interface.
!! Form A requires first calling horiz_interp_new, while Form B
!! requires no initialization.
!!
!! @param Interp
!!     Derived-type variable containing interpolation indices and weights.
!!     Returned by a previous call to horiz_interp_new.
!!
!! @param data_in
!!      Input data on source grid.
!!
!! @param verbose
!!      flag for the amount of print output.
!!               verbose = 0, no output; = 1, min,max,means; = 2, still more
!!
!! @param mask_in
!!      Input mask, must be the same size as the input data. The real value of
!!      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points
!!      that should not be used or have missing data. It is Not needed for
!!      spherical regrid.
!!
!! @param missing_value
!!      Use the missing_value to indicate missing data.
!!
!! @param missing_permit
!!      numbers of points allowed to miss for the bilinear interpolation. The value
!!      should be between 0 and 3.
!!
!! @param lon_in, lat_in
!!      longitude and latitude (in radians) of source grid. More explanation can
!!      be found in the documentation of horiz_interp_new.
!!
!! @param lon_out, lat_out
!!      longitude and latitude (in radians) of destination grid. More explanation can
!!      be found in the documentation of horiz_interp_new.
!!
!! @param data_out
!!      Output data on destination grid.
!!
!! @param mask_out
!!      Output mask that specifies whether data was computed.
!!
!!
!! @throws FATAL, size of input array incorrect
!!      The input data array does not match the size of the input grid edges
!!      specified. If you are using the initialization interface make sure you
!!      have the correct grid size.
!!
!! @throws FATAL, size of output array incorrect
!!      The output data array does not match the size of the input grid
!!      edges specified. If you are using the initialization interface make
!!      sure you have the correct grid size.
!> @ingroup horiz_interp_mod
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

!> Private helper routines
interface is_lat_lon
    module procedure is_lat_lon_r4
    module procedure is_lat_lon_r8
end interface

interface horiz_interp_solo_1d
  module procedure horiz_interp_solo_1d_r4
  module procedure horiz_interp_solo_1d_r8
end interface


!> @addtogroup horiz_interp_mod
!> @{

 logical :: reproduce_siena = .false. !< Set reproduce_siena = .true. to reproduce siena results.
                 !! Set reproduce_siena = .false. to decrease truncation error
                 !! in routine poly_area in file mosaic_util.c. The truncation error of
                 !! second order conservative remapping might be big for high resolution
                 !! grid.

 namelist /horiz_interp_nml/ reproduce_siena

!-----------------------------------------------------------------------
! Include variable "version" to be written to log file.
#include<file_version.h>
 logical            :: module_is_initialized = .FALSE.
!-----------------------------------------------------------------------

contains

!#######################################################################

  !> Initialize module and writes version number to logfile.out
  subroutine horiz_interp_init
  integer :: unit, ierr, io

  if(module_is_initialized) return
  call write_version_number("HORIZ_INTERP_MOD", version)

  read (input_nml_file, horiz_interp_nml, iostat=io)
  ierr = check_nml_error(io,'horiz_interp_nml')
  if (mpp_pe() == mpp_root_pe() ) then
     unit = stdlog()
     write (unit, nml=horiz_interp_nml)
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

!> Deallocates memory used by "horiz_interp_type" variables.
!! Must be called before reinitializing with horiz_interp_new.
 subroutine horiz_interp_del ( Interp )

   type (horiz_interp_type), intent(inout) :: Interp !< A derived-type variable returned by previous
                                           !! call to horiz_interp_new. The input variable must have
                                           !! allocated arrays. The returned variable will contain
                                           !! deallocated arrays

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
   case (SPHERICA)
      call horiz_interp_spherical_del(Interp )
   end select

   Interp%I_am_initialized = .false.
!-----------------------------------------------------------------------

 end subroutine horiz_interp_del

 !#####################################################################

 !> Dummy routine
 subroutine horiz_interp_end
 return
 end subroutine horiz_interp_end

#include "horiz_interp_r4.fh"
#include "horiz_interp_r8.fh"

end module horiz_interp_mod
!> @}
! close documentation grouping
