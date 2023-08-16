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
!> @file
!> @author Caitlyn McAllister
!> @brief Unit test for astronomy/test_topography
!> @email gfdl.climate.model.info@noaa.gov
!> @description FILL THIS OUT

program test_top

  use gaussian_topog_mod, only: gaussian_topog_init, get_gaussian_topog
  use topography_mod,     only: topography_init, get_topog_mean, get_topog_stdev, &
                                get_ocean_frac, get_ocean_mask, get_water_frac, &
                                get_water_mask
  use fms_mod,            only: fms_init, fms_end
  use fms2_io_mod,        only: fms2_io_init, FmsNetcdfFile_t, open_file, close_file, register_axis, register_field, &
                                register_variable_attribute, write_data, read_data, unlimited
  use mpp_mod,            only: mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_sync, input_nml_file
  use horiz_interp_mod,   only: horiz_interp_type, horiz_interp_new, &
                                horiz_interp, horiz_interp_del
  use constants_mod,      only: pi
  use platform_mod,       only: r4_kind, r8_kind
  
  implicit none

  type(FmsNetcdfFile_t)     :: top_fileobj                 ! fileobj for fms2_io
  character(len=128)        :: topog_file, water_file      ! filenames needed for topography_mod
  real(kind=TEST_TOP_KIND_) :: ipts(2), jpts(2)            ! axis for files
  real(kind=TEST_TOP_KIND_) :: iptsp1(3), jptsp1(3)        ! sub axis, size +1 from other axis
  real(kind=TEST_TOP_KIND_) :: xdat(3), ydat(3), zdat(2,2) ! specifc data topog_mod looks for
  integer                   :: i, j                        ! counters for loops
  integer                   :: io_status                   ! namelist
  integer, parameter        :: lkind = TEST_TOP_KIND_      ! kind parameter for mixed precision

  call fms_init
  call topography_init

  ! name files
  topog_file = "topography.data.nc"
  water_file = "water.data.nc"

  ! create data for both topog and water files
  ipts = (/-180.0_lkind, 180.0_lkind/)
  jpts = (/-90.0_lkind, 90.0_lkind/)

  iptsp1 = (/-180.0_lkind, 0.0_lkind, 180.0_lkind/)
  jptsp1 = (/-90.0_lkind, 0.0_lkind, 90.0_lkind/)

  xdat = (/1.0_lkind, 2.0_lkind, 3.0_lkind/)
  ydat = (/1.0_lkind, 2.0_lkind, 3.0_lkind/)

  !do i = 1, 2
  !  do j = 1, 2
  !    zdat(j,i) = (real(i, TEST_TOP_KIND_) + real(j, TEST_TOP_KIND_)) - 1.0_lkind
  !  end do
  !end do
  zdat = 24.0
  ! write topog file
  if (open_file(top_fileobj, topog_file, "overwrite")) then
    call register_axis(top_fileobj, "ipts",   2)
    call register_axis(top_fileobj, "jpts",   2)
    call register_axis(top_fileobj, "iptsp1", 3)
    call register_axis(top_fileobj, "jptsp1", 3)

    call register_field(top_fileobj, "ipts",   "double", dimensions=(/"ipts"/))
    call register_field(top_fileobj, "jpts",   "double", dimensions=(/"jpts"/))
    call register_field(top_fileobj, "iptsp1", "double", dimensions=(/"iptsp1"/))
    call register_field(top_fileobj, "jptsp1", "double", dimensions=(/"jptsp1"/))
    call register_field(top_fileobj, "xdat",   "double", dimensions=(/"iptsp1"/))
    call register_field(top_fileobj, "ydat",   "double", dimensions=(/"jptsp1"/))
    call register_field(top_fileobj, "zdat",   "double", dimensions=(/"ipts", "jpts"/))

    call write_data(top_fileobj, "ipts",   ipts)
    call write_data(top_fileobj, "jpts",   jpts)
    call write_data(top_fileobj, "iptsp1", iptsp1)
    call write_data(top_fileobj, "jptsp1", jptsp1)
    call write_data(top_fileobj, "xdat",   xdat)
    call write_data(top_fileobj, "ydat",   ydat)
    call write_data(top_fileobj, "zdat",   zdat)

    call close_file(top_fileobj)

  else
    call mpp_error(FATAL, "test_topography: error opening topog_file")
  end if

  call test_topog_mean

  call fms_end

  contains

  subroutine test_topog_mean

    implicit none
    real(kind=TEST_TOP_KIND_), dimension(2,2)              :: lon2d, lat2d
    real(kind=TEST_TOP_KIND_), dimension(:,:), allocatable :: zmean
    logical                                                :: get_topog_mean2d
    integer                                                :: ix, iy

    integer :: js, je
    type (horiz_interp_type) :: Interp

    real(kind=TEST_TOP_KIND_)              :: ybeg, yend
    integer                               :: j

    lon2d = 0.0_lkind ; lat2d = 0.0_lkind
    ix = size(lon2d,1) - 1 ; iy = size(lat2d,2) - 1
    allocate (zmean(ix, iy))


    get_topog_mean2d = get_topog_mean(lon2d, lat2d, zmean)
    !print *, ref_1d_zmean
    !deallocate (ref_1d_zmean)

    ! input_data(TOPOG_INDEX, xdat, ydat, zdat)
    if (open_file(top_fileobj, topog_file, "read")) then
    call read_data(top_fileobj, 'xdat', xdat)
    call read_data(top_fileobj, 'ydat', ydat)
    call read_data(top_fileobj, 'zdat', zdat)
    end if

    ! expected output
    print *, "input_data", xdat
    print *, "input_data", ydat
    print *, "input_data", zdat

    ! find_indices(minval(blat), maxval(blat), ydat, js, je)
    js = 1
    do j = 1, size(ydat(:))-1
       if (ybeg >= ydat(j) .and. ybeg <= ydat(real(j,TEST_TOP_KIND_)+1.0_lkind)) then
          js = j
          exit
       endif
    enddo
 
    je = size(ydat(:))-1
    do j = js, size(ydat(:))-1
       if (yend >= ydat(j) .and. yend <= ydat(real(j,TEST_TOP_KIND_)+1.0_lkind)) then
          je = j
          exit
       endif
    enddo

    ! expected output
    print *, "find_indices", js
    print *, "find_indices", je

    call horiz_interp_new (Interp, xdat, ydat(js:je+1), lon2d, lat2d)
    call horiz_interp     (Interp, zdat(:,js:je), zmean)

    print *, "zmean", zmean


  end subroutine test_topog_mean



end program test_top