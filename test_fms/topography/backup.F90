



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
    real(kind=TEST_TOP_KIND_), dimension(3)                :: lon1d, lat1d
    real(kind=TEST_TOP_KIND_), dimension(:,:), allocatable :: ref_1d_zmean
    logical                                                :: mean_answer1d
    integer                                                :: ix, iy

    lon1d = 0.0_lkind ; lat1d = 0.0_lkind
    ix = size(lon1d) - 1 ; iy = size(lat1d) - 1
    allocate (ref_1d_zmean(ix, iy))

    mean_answer1d = get_topog_mean(lon1d, lat1d, ref_1d_zmean)
    !print *, ref_1d_zmean
    deallocate (ref_1d_zmean)


  end subroutine test_topog_mean



end program test_top

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

  type(FmsNetcdfFile_t)     :: fileobj
  character(len=128)        :: topog_file, water_file
  real(kind=TEST_TOP_KIND_) :: height(4)
  real(kind=TEST_TOP_KIND_) :: olon(4)
  real(kind=TEST_TOP_KIND_) :: olat(4)
  real(kind=TEST_TOP_KIND_) :: wlon(4)
  real(kind=TEST_TOP_KIND_) :: wlat(4)
  real(kind=TEST_TOP_KIND_) :: ipts(2), jpts(2), iptsp1(3), jptsp1(3)
  real(kind=TEST_TOP_KIND_) :: xdat(3), ydat(3), zdat(2,2)
  integer                   :: io_status, i, j
  real(kind=TEST_TOP_KIND_), dimension(2, 2) :: topog_data, water_data

  namelist /test_topog_nml/ height, olon, olat, wlon, wlat

  call fms_init
  call topography_init

  ! create data files
  topog_file = "DATA/navy_topography.data.nc"
  water_file = "DATA/navy_pctwater.data.nc"
  ipts = (/-180.0, 180.0/)
  jpts = (/-90.0, 90.0/)
  iptsp1 = (/-180.0, 0.0, 180.0/)
  jptsp1 = (/-90.0, 0.0, 90.0/)

  xdat = (/1.0, 2.0, 3.0/)
  ydat = (/1.0, 2.0, 3.0/)

  do i = 1, 2
    do j = 1, 2
      zdat(j,i) = (real(i, TEST_TOP_KIND_) + real(j, TEST_TOP_KIND_)) - 1.0
    end do
  end do

  do i = 1, 2
    do j = 1, 2
      water_data(j,i) = (real(i, TEST_TOP_KIND_) + real(j, TEST_TOP_KIND_)) - 1.0
    end do
  end do

  ! topog file
  if (open_file(fileobj, topog_file, "overwrite")) then
    call register_axis(fileobj, "ipts", 2)
    call register_axis(fileobj, "jpts", 2)
    call register_axis(fileobj, "iptsp1", 3)
    call register_axis(fileobj, "jptsp1", 3)

    call register_field(fileobj, "ipts", "double", dimensions=(/"ipts"/))
    call register_field(fileobj, "jpts", "double", dimensions=(/"jpts"/))
    call register_field(fileobj, "iptsp1", "double", dimensions=(/"iptsp1"/))
    call register_field(fileobj, "jptsp1", "double", dimensions=(/"jptsp1"/))
    call register_field(fileobj, "xdat", "double", dimensions=(/"iptsp1"/))
    call register_field(fileobj, "ydat", "double", dimensions=(/"jptsp1"/))
    call register_field(fileobj, "zdat", "double", dimensions=(/"ipts", "jpts"/))

    call write_data(fileobj, "ipts", ipts)
    call write_data(fileobj, "jpts", jpts)
    call write_data(fileobj, "iptsp1", iptsp1)
    call write_data(fileobj, "jptsp1", jptsp1)
    call write_data(fileobj, "xdat", xdat)
    call write_data(fileobj, "ydat", ydat)
    call write_data(fileobj, "zdat", zdat)

    call close_file(fileobj)
  else 
    call mpp_error(FATAL, "test_topography: error opening topog_file")
  end if

  ! water file
  if (open_file(fileobj, water_file, "overwrite")) then
    call register_axis(fileobj, "nlon", 2)
    call register_axis(fileobj, "nlat", 2)

    call register_field(fileobj, "nlon", "int", dimensions=(/"nlon"/))
    call register_field(fileobj, "nlat", "int", dimensions=(/"nlat"/))
    call register_field(fileobj, "blon", "double", dimensions=(/"nlon"/))
    call register_field(fileobj, "blat", "double", dimensions=(/"nlat"/))
    call register_field(fileobj, "water_data", "double", dimensions=(/"nlon", "nlat"/))

    call write_data(fileobj, "nlon", nlon)
    call write_data(fileobj, "nlat", nlat)
    call write_data(fileobj, "blon", nlon+1)
    call write_data(fileobj, "blat", nlat+1)
    call write_data(fileobj, "water_data", water_data)

    call close_file(fileobj)
  else 
    call mpp_error(FATAL, "test_topography: error opening water_file")
  end if

  read (input_nml_file, test_topog_nml, iostat=io_status)
  if (io_status > 0) then
     call mpp_error(FATAL,'test_topog: Error reading input.nml')
  endif

  call try_this

  call fms_end

  contains

  subroutine try_this

    implicit none
    real(kind=TEST_TOP_KIND_), dimension(:,:), allocatable :: ref_1d_zmean, ref_2d_zmean
    logical                                                :: mean_fnc1d, mean_fnc2d
    integer, parameter                                     :: lkind = TEST_TOP_KIND_
    integer                                                :: ix1d, iy1d

    if (open_file(fileobj, topog_file, "read")) then

      allocate(ref_1d_zmean(ix1d,iy1d))

      mean_fnc1d = get_topog_mean(blon, blat, ref_1d_zmean)
    end if 

  end subroutine try_this

end program test_top

program test_top

  use gaussian_topog_mod, only: gaussian_topog_init, get_gaussian_topog
  use topography_mod,     only: topography_init, get_topog_mean, get_topog_stdev, &
                                get_ocean_frac, get_ocean_mask, get_water_frac, &
                                get_water_mask
  use fms_mod,            only: fms_init, fms_end
  use fms2_io_mod,        only: fms2_io_init, FmsNetcdfFile_t, open_file, close_file, register_axis, register_field, &
                                register_variable_attribute, write_data, unlimited
  use mpp_mod,            only: mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_sync
  use horiz_interp_mod,   only: horiz_interp_type, horiz_interp_new, &
                                horiz_interp, horiz_interp_del
  use constants_mod,      only: pi
  use platform_mod,       only: r4_kind, r8_kind

  namelist /test_topog_nml/ height, olon, olat, wlon, wlat

  implicit none

  character(len=128)                                     :: filename = 'INPUT/test_topog.nc'
  type(FmsNetcdfFile_t)                                  :: fileobj


  call fms_init
  call topography_init

  call create_file

  call test_topog_mean ; call test_topog_stdev
  call test_ocean_frac ; call test_ocean_mask
  call test_water_frac ; call test_water_mask

  call fms_end

  contains

  subroutine test_topog_mean
    real(kind=TEST_TOP_KIND_), dimension(:,:), allocatable :: ref_1d_zmean, ref_2d_zmean
    logical                                                :: mean_fnc1d, mean_fnc2d

    call create_file
      ! -------------------------- test get_topog_mean_1d -------------------------- !

      allocate(ref_1d_zmean(ix1d,iy1d) ) ; ref_1d_zmean = 0.0_lkind

      mean_fnc1d = get_topog_mean(lon1d, lat1d, ref_1d_zmean)
      print *, mean_fnc1d

      ! -------------------------- test get_topog_mean_2d -------------------------- !

      allocate(ref_2d_zmean(ix2d,iy2d) ) ; ref_2d_zmean = 0.0

      mean_fnc2d = get_topog_mean(lon2d, lat2d, ref_2d_zmean)

      ! --------------------------------- end tests -------------------------------- !

      deallocate(ref_1d_zmean, ref_2d_zmean)
    else 
      call mpp_error(FATAL, "Error reading test_topography.nc")
    end if
  
  end subroutine test_topog_mean

  subroutine test_topog_stdev
    real(kind=TEST_TOP_KIND_), dimension(:,:), allocatable :: ref_1d_stdev, ref_2d_stdev
    logical                                                :: stdev_fnc1d, stdev_fnc2d

    ! -------------------------- test get_topog_stdev_1d -------------------------- !

    allocate(ref_1d_stdev(ix1d,iy1d) ) ; ref_1d_stdev = 0.0

    stdev_fnc1d = get_topog_stdev(lon1d, lat1d, ref_1d_stdev)

    ! -------------------------- test get_topog_stdev_2d -------------------------- !

    allocate(ref_2d_stdev(ix2d,iy2d) ) ; ref_2d_stdev = 0.0

    stdev_fnc2d = get_topog_stdev(lon2d, lat2d, ref_2d_stdev)

    ! --------------------------------- end tests --------------------------------- !

    deallocate(ref_1d_stdev, ref_2d_stdev)

  end subroutine test_topog_stdev

  subroutine test_ocean_frac
    real(kind=TEST_TOP_KIND_), dimension(:,:), allocatable :: ref_1d_ocnfrac, ref_2d_ocnfrac
    logical                                                :: ocn_frac_fnc1d, ocn_frac_fnc2d

    ! -------------------------- test get_topog_ocean_frac_1d -------------------------- !

    allocate(ref_1d_ocnfrac(ix1d,iy1d) ) ; ref_1d_ocnfrac = 0.0

    ocn_frac_fnc1d = get_ocean_frac(lon1d, lat1d, ref_1d_ocnfrac)

    ! -------------------------- test get_topog_ocean_frac_2d -------------------------- !

    allocate(ref_2d_ocnfrac(ix2d,iy2d) ) ; ref_2d_ocnfrac = 0.0

    ocn_frac_fnc2d = get_ocean_frac(lon2d, lat2d, ref_2d_ocnfrac)

    ! ------------------------------------ end tests ------------------------------------ !

    deallocate(ref_1d_ocnfrac, ref_2d_ocnfrac)

  end subroutine test_ocean_frac

  subroutine test_ocean_mask
    logical, dimension(:,:), allocatable :: ref_1d_ocnmask, ref_2d_ocnmask
    logical                              :: ocn_mask_fnc1d, ocn_mask_fnc2d

    ! -------------------------- test get_topog_ocean_mask_1d -------------------------- !

    allocate(ref_1d_ocnmask(ix1d,iy1d) )

    ocn_mask_fnc1d = get_ocean_mask(lon1d, lat1d, ref_1d_ocnmask)

    ! -------------------------- test get_topog_ocean_mask_2d -------------------------- !

    allocate(ref_2d_ocnmask(ix2d,iy2d) )

    ocn_mask_fnc2d = get_ocean_mask(lon2d, lat2d, ref_2d_ocnmask)

    ! ------------------------------------ end tests ------------------------------------ !

    deallocate(ref_1d_ocnmask, ref_2d_ocnmask)

  end subroutine test_ocean_mask

  subroutine test_water_frac
    real(kind=TEST_TOP_KIND_), dimension(:,:), allocatable :: ref_1d_wtrfrac, ref_2d_wtrfrac
    logical                                                :: wtr_frac_fnc1d, wtr_frac_fnc2d

    ! -------------------------- test get_topog_water_frac_1d -------------------------- !

    allocate(ref_1d_wtrfrac(ix1d,iy1d) ) ; ref_1d_wtrfrac = 0.0

    wtr_frac_fnc1d = get_water_frac(lon1d, lat1d, ref_1d_wtrfrac)

    ! -------------------------- test get_topog_water_frac_2d -------------------------- !

    allocate(ref_2d_wtrfrac(ix2d,iy2d) ) ; ref_2d_wtrfrac = 0.0

    wtr_frac_fnc2d = get_water_frac(lon2d, lat2d, ref_2d_wtrfrac)

    ! ------------------------------------ end tests ------------------------------------ !

    deallocate(ref_1d_wtrfrac, ref_2d_wtrfrac)

  end subroutine test_water_frac

  subroutine test_water_mask
    logical, dimension(:,:), allocatable :: ref_1d_wtrmask, ref_2d_wtrmask
    logical                              :: wtr_mask_fnc1d, wtr_mask_fnc2d

    ! -------------------------- test get_topog_water_mask_1d -------------------------- !

    allocate(ref_1d_wtrmask(ix1d,iy1d) )

    wtr_mask_fnc1d = get_water_mask(lon1d, lat1d, ref_1d_wtrmask)

    ! -------------------------- test get_topog_ocean_mask_2d -------------------------- !

    allocate(ref_2d_wtrmask(ix2d,iy2d) )

    wtr_mask_fnc2d = get_water_mask(lon2d, lat2d, ref_2d_wtrmask)

    ! ------------------------------------ end tests ------------------------------------ !

    deallocate(ref_1d_wtrmask, ref_2d_wtrmask)

  end subroutine test_water_mask

  subroutine create_file

    implicit none
    real(kind=TEST_TOP_KIND_), allocatable, dimension(:,:) :: mtn_data
    character(len=128)                                     :: fieldname = "mtn_data"
    integer                                                :: i, ix, iy

    if (open_file(fileobj, filename, "overwrite")) then
      call register_axis(fileobj, "lon", 5)
      call register_axis(fileobj, "lat", 5)
      call register_axis(fileobj, "time", unlimited)

      call register_field(fileobj, "lon", "double", dimensions=(/"lon"/))
      call register_field(fileobj, "lat", "double", dimensions=(/"lat"/))
      call register_field(fileobj, "time", "double", dimensions=(/"time"/))
      call register_field(fileobj, fieldname, "double", dimensions=(/"lon", "lat"/))

      call register_variable_attribute(fileobj, "lon", "cartesian_axis", "X", str_len=1)
      call register_variable_attribute(fileobj, "lat", "cartesian_axis", "Y", str_len=1)
    
      call write_data(fileobj, "lon", (/(0.0+i*2.0,i=1,5)/))
      call write_data(fileobj, "lat", (/(0.0+i*2.0,i=1,5)/))
      call write_data(fileobj, "time", (/(1+i*2, i=0,2)/))

      allocate (mtn_data(5,5))
      mtn_data = 2.0_lkind
      call write_data(fileobj, "mtn_data", mtn_data)

      call close_file(fileobj)
    else
      call mpp_error(FATAL, "Error opening test_topography.nc to write")
    end if

    deallocate (mtn_data)
  
  end subroutine create_file

end program test_top