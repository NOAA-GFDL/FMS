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
!> @defgroup horiz_interp_conserve_mod horiz_interp_conserve_mod
!> @ingroup horiz_interp
!> @brief Performs spatial interpolation between grids using conservative interpolation
!!
!> @author Bruce Wyman, Zhi Liang
!!
!> This module can conservatively interpolate data from any logically rectangular grid
!! to any rectangular grid. The interpolation scheme is area-averaging
!! conservative scheme. There is an optional mask field for missing input data in both
!! horiz_interp__conserveinit and horiz_interp_conserve. For efficiency purpose, mask should only be
!! kept in horiz_interp_init (will remove the mask in horiz_interp in the future).
!! There are 1-D and 2-D version of horiz_interp_conserve_init for 1-D and 2-D grid.
!! There is a optional argument mask in horiz_interp_conserve_init_2d and no mask should
!! to passed into horiz_interp_conserv. optional argument mask will not be passed into
!! horiz_interp_conserve_init_1d and optional argument mask may be passed into
!! horiz_interp_conserve (For the purpose of reproduce Memphis??? results).
!! An optional output mask field may be used in conjunction with the input mask to show
!! where output data exists.

module horiz_interp_conserve_mod

  use platform_mod,          only: r4_kind, r8_kind
  use fms2_io_mod,           only: FmsNetcdfFile_t, open_file, read_data, close_file, get_dimension_size
  use mpp_mod,               only: mpp_send, mpp_recv, mpp_pe, mpp_root_pe, mpp_npes
  use mpp_mod,               only: mpp_error, FATAL,  mpp_sync_self
  use mpp_mod,               only: COMM_TAG_1, COMM_TAG_2
  use fms_mod,               only: write_version_number
  use grid2_mod,             only: get_great_circle_algorithm
  use constants_mod,         only: PI
  use horiz_interp_type_mod, only: horiz_interp_type, CONSERVE


  implicit none
  private

  ! public interface


  !> @brief Allocates space and initializes a derived-type variable
  !! that contains pre-computed interpolation indices and weights.
  !!
  !> Allocates space and initializes a derived-type variable
  !! that contains pre-computed interpolation indices and weights
  !! for improved performance of multiple interpolations between
  !! the same grids.
  !! @param lon_in
  !!      Longitude (in radians) for source data grid.
  !!
  !! @param lat_in
  !!      Latitude (in radians) for source data grid.
  !!
  !! @param lon_out
  !!      Longitude (in radians) for destination data grid.
  !!
  !! @param lat_out
  !!      Latitude (in radians) for destination data grid.
  !!
  !! @param verbose
  !!      flag for the amount of print output.
  !!
  !! @param mask_in
  !!      Input mask.  must be the size (size(lon_in)-1, size(lon. The real value of
  !!      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points
  !!      that should not be used or have missing data.
  !!
  !! @param mask_out
  !!      Output mask that specifies whether data was computed.
  !!
  !! @param Interp
  !!      A derived-type variable containing indices and weights used for subsequent
  !!      interpolations. To reinitialize this variable for a different grid-to-grid
  !!      interpolation you must first use the "horiz_interp_del" interface.
  !!
  !> @ingroup horiz_interp_conserve_mod
  interface horiz_interp_conserve_new
     module procedure horiz_interp_conserve_new_1dx1d_r4
     module procedure horiz_interp_conserve_new_1dx2d_r4
     module procedure horiz_interp_conserve_new_2dx1d_r4
     module procedure horiz_interp_conserve_new_2dx2d_r4
     module procedure horiz_interp_conserve_new_1dx1d_r8
     module procedure horiz_interp_conserve_new_1dx2d_r8
     module procedure horiz_interp_conserve_new_2dx1d_r8
     module procedure horiz_interp_conserve_new_2dx2d_r8
  end interface

  interface horiz_interp_conserve
    module procedure horiz_interp_conserve_r4
    module procedure horiz_interp_conserve_r8
  end interface

!> private helper routines
  interface data_sum
    module procedure data_sum_r4
    module procedure data_sum_r8
  end interface

  interface stats
    module procedure stats_r4
    module procedure stats_r8
  end interface

  interface horiz_interp_conserve_version1
    module procedure horiz_interp_conserve_version1_r8
    module procedure horiz_interp_conserve_version1_r4
  end interface

  interface horiz_interp_conserve_version2
    module procedure horiz_interp_conserve_version2_r8
    module procedure horiz_interp_conserve_version2_r4
  end interface



  !> @addtogroup horiz_interp_conserve_mod
  !> @{
  public :: horiz_interp_conserve_init
  public :: horiz_interp_read_weights_conserve
  public :: horiz_interp_conserve_new, horiz_interp_conserve, horiz_interp_conserve_del

  integer :: pe, root_pe
  !-----------------------------------------------------------------------
  ! Include variable "version" to be written to log file.
#include<file_version.h>
  logical            :: module_is_initialized = .FALSE.

  logical         :: great_circle_algorithm = .false.

contains

  !> Writes version number to logfile.
  subroutine horiz_interp_conserve_init

    if(module_is_initialized) return
    call write_version_number("HORIZ_INTERP_CONSERVE_MOD", version)

    great_circle_algorithm = get_great_circle_algorithm()

    module_is_initialized = .true.

  end subroutine horiz_interp_conserve_init

  !> Deallocates memory used by "HI_KIND_TYPE" variables.
  !! Must be called before reinitializing with horiz_interp_new.
  subroutine horiz_interp_conserve_del ( Interp )

    type (horiz_interp_type), intent(inout) :: Interp !< A derived-type variable returned by
                         !! previous call to horiz_interp_new. The input variable must have
                         !! allocated arrays. The returned variable will contain deallocated arrays.

    select case(Interp%version)
    case (1)
      if( Interp%horizInterpReals8_type%is_allocated) then
        if(allocated(Interp%horizInterpReals8_type%area_src)) deallocate(Interp%horizInterpReals8_type%area_src)
        if(allocated(Interp%horizInterpReals8_type%area_dst)) deallocate(Interp%horizInterpReals8_type%area_dst)
        if(allocated(Interp%horizInterpReals8_type%facj))     deallocate(Interp%horizInterpReals8_type%facj)
        if(allocated(Interp%jlat))                 deallocate(Interp%jlat)
        if(allocated(Interp%horizInterpReals8_type%faci))     deallocate(Interp%horizInterpReals8_type%faci)
        if(allocated(Interp%ilon))                 deallocate(Interp%ilon)
      else if( Interp%horizInterpReals4_type%is_allocated) then
        if(allocated(Interp%horizInterpReals4_type%area_src)) deallocate(Interp%horizInterpReals4_type%area_src)
        if(allocated(Interp%horizInterpReals4_type%area_dst)) deallocate(Interp%horizInterpReals4_type%area_dst)
        if(allocated(Interp%horizInterpReals4_type%facj))     deallocate(Interp%horizInterpReals4_type%facj)
        if(allocated(Interp%jlat))                 deallocate(Interp%jlat)
        if(allocated(Interp%horizInterpReals4_type%faci))     deallocate(Interp%horizInterpReals4_type%faci)
        if(allocated(Interp%ilon))                 deallocate(Interp%ilon)
      endif
    case (2)
      if( Interp%horizInterpReals8_type%is_allocated) then
        if(allocated(Interp%i_src)) deallocate(Interp%i_src)
        if(allocated(Interp%j_src)) deallocate(Interp%j_src)
        if(allocated(Interp%i_dst)) deallocate(Interp%i_dst)
        if(allocated(Interp%j_dst)) deallocate(Interp%j_dst)
        if(allocated(Interp%horizInterpReals8_type%area_frac_dst)) &
            deallocate(Interp%horizInterpReals8_type%area_frac_dst)
      else if( Interp%horizInterpReals4_type%is_allocated ) then
        if(allocated(Interp%i_src)) deallocate(Interp%i_src)
        if(allocated(Interp%j_src)) deallocate(Interp%j_src)
        if(allocated(Interp%i_dst)) deallocate(Interp%i_dst)
        if(allocated(Interp%j_dst)) deallocate(Interp%j_dst)
        if(allocated(Interp%horizInterpReals4_type%area_frac_dst)) &
            deallocate(Interp%horizInterpReals4_type%area_frac_dst)
       endif
    end select
    if(allocated(Interp%xgrid_area)) deallocate(Interp%xgrid_area)
    Interp%horizInterpReals4_type%is_allocated = .false.
    Interp%horizInterpReals8_type%is_allocated = .false.

  end subroutine horiz_interp_conserve_del

  subroutine horiz_interp_read_weights_conserve(Interp, weight_filename, weight_file_source, &
    nlon_src, nlat_src, nlon_dst, nlat_dst, isw, iew, jsw, jew, src_tile)

    type(horiz_interp_type), intent(inout) :: Interp
    character(len=*), intent(in) :: weight_filename
    character(len=*), intent(in) :: weight_file_source
    integer, intent(in) :: nlon_src, nlat_src, nlon_dst, nlat_dst
    integer, intent(in) :: isw, iew, jsw, jew
    integer, intent(in), optional :: src_tile

    integer :: i, j, ncells, domain_ncells
    integer :: istart, iend, i_dst, j_dst, index
    real(8), allocatable :: dst_area2(:,:), dst_area1(:), read1(:), xarea(:)
    integer, allocatable :: tile1(:), read2(:,:)
    logical, allocatable :: mask(:)

    type(FmsNetcdfFile_t) :: weight_fileobj !< FMS2io fileob for the weight file

    ! check if weight_file was generated from fregrid
    if(trim(weight_file_source) /= "fregrid") then
      call mpp_error(FATAL, trim(weight_file_source)//&
        &" is not a supported weight file source. fregrid is the only supported weight file source.")
    end if

    if(open_file(weight_fileobj, trim(weight_filename), "read")) then

      ! get ncells
      call get_dimension_size(weight_fileobj, "ncells", ncells)

      istart = 1
      iend = ncells

      !get section of xgrid on src_tile
      if(present(src_tile)) then
        allocate(tile1(ncells))
        call read_data(weight_fileobj, "tile1", tile1)
        !find istart
        do i=1, ncells
          if(tile1(i) == src_tile) then
            istart = i
            exit
          end if
        end do
        !find iend
        do i=istart, ncells
          if(tile1(i) /= src_tile) then
            iend = i - 1
            exit
          end if
        end do
        ncells = iend - istart + 1
        deallocate(tile1)
      end if

      ! allocate arrays for reading data
      allocate(read2(2, ncells))

      ! get section of xgrid for the specified window (compute domain) on the tgt grid
      call read_data(weight_fileobj, "tile2_cell", read2, corner=[1,istart], edge_lengths=[2,ncells])

      ! get xgrid indices, used copilot for this section of code
      allocate(mask(ncells))
      mask = (read2(1,:) >= isw .and. read2(1,:) <= iew .and. read2(2,:) >= jsw .and. read2(2,:) <= jew)

      domain_ncells = count(mask)

      write(*,*) istart, iend, ncells, "and", isw, iew, jsw, jew, domain_ncells

      ! allocate data to store xgrid
      allocate(Interp%i_src(domain_ncells))
      allocate(Interp%j_src(domain_ncells))
      allocate(Interp%i_dst(domain_ncells))
      allocate(Interp%j_dst(domain_ncells))
      allocate(Interp%horizInterpReals8_type%area_frac_dst(domain_ncells))

      Interp%i_dst = pack(read2(1,:), mask) - isw + 1
      Interp%j_dst = pack(read2(2,:), mask) - jsw + 1

      !save src parent cell indices
      call read_data(weight_fileobj, "tile1_cell", read2, corner=[1,istart], edge_lengths=[2,ncells])
      Interp%i_src = pack(read2(1,:), mask)
      Interp%j_src = pack(read2(2,:), mask)

      deallocate(read2)

      ! allocate arrays to compute weights
      allocate(read1(ncells), dst_area1(domain_ncells), dst_area2(nlon_dst, nlat_dst), xarea(domain_ncells))

      ! read xgrid area
      call read_data(weight_fileobj, "xgrid_area", read1, corner=[istart], edge_lengths=[ncells])

      xarea = pack(read1, mask)

      !sum over xgrid area to get destination grid area
      dst_area2 = 0.0
      do i = 1, domain_ncells
        i_dst = Interp%i_dst(i)
        j_dst = Interp%j_dst(i)
        dst_area2(i_dst, j_dst) = dst_area2(i_dst, j_dst) + xarea(i)
        dst_area1(i) = dst_area2(i_dst, j_dst)
      end do

      Interp%horizInterpReals8_type%area_frac_dst = xarea/dst_area1

      deallocate(read1)
      deallocate(dst_area1)
      deallocate(dst_area2)
      deallocate(xarea)

      call close_file(weight_fileobj)

    else
      call mpp_error(FATAL, "cannot open weight file")
    end if

    Interp%nxgrid = domain_ncells
    Interp%nlon_src = nlon_src
    Interp%nlat_src = nlat_src
    Interp%nlon_dst = nlon_dst
    Interp%nlat_dst = nlat_dst
    Interp%horizInterpReals8_type%is_allocated = .true.
    Interp%interp_method = CONSERVE
    Interp%version = 2
    Interp%I_am_initialized = .true.

  end subroutine horiz_interp_read_weights_conserve

#include "horiz_interp_conserve_r4.fh"
#include "horiz_interp_conserve_r8.fh"

end module horiz_interp_conserve_mod
!> @}
! close documentation grouping
