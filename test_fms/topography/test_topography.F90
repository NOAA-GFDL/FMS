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
  real(kind=TEST_TOP_KIND_) :: xdat(3), ydat(3), zdat(2,2) ! specifc data topog_mod looks for
  integer                   :: ipts, jpts                  ! axis for files
  integer                   :: iptsp1, jptsp1              ! serves as a counter for data
  integer                   :: ipts_r, jpts_r              ! sub axis
  integer, parameter        :: lkind = TEST_TOP_KIND_      ! kind parameter for mixed precision

  real(kind=TEST_TOP_KIND_), parameter        :: deg2rad = real(pi, TEST_TOP_KIND_)/180.0_lkind
  real(kind=TEST_TOP_KIND_), dimension(2,2)   :: lon2d, lat2d   ! in radians
  real(kind=TEST_TOP_KIND_), dimension(4)     :: lon1d, lat1d   ! in radians

  call fms_init
  call topography_init

  ! define blat and blon, in this test they'll be referred to as lat2d/lon2d, lat1d/lon1d
  lon2d(1,1) = 1.5_lkind*deg2rad ; lat2d(1,1) = 1.5_lkind*deg2rad
  lon2d(2,1) = 2.5_lkind*deg2rad ; lat2d(2,1) = 1.5_lkind*deg2rad
  lon2d(1,2) = 1.5_lkind*deg2rad ; lat2d(1,2) = 2.5_lkind*deg2rad
  lon2d(2,2) = 2.5_lkind*deg2rad ; lat2d(2,2) = 2.5_lkind*deg2rad

  lon1d(1) = 1.5_lkind*deg2rad ; lat1d(1) = 1.5_lkind*deg2rad
  lon1d(2) = 2.5_lkind*deg2rad ; lat1d(2) = 1.5_lkind*deg2rad
  lon1d(3) = 1.5_lkind*deg2rad ; lat1d(3) = 2.5_lkind*deg2rad
  lon1d(4) = 2.5_lkind*deg2rad ; lat1d(4) = 2.5_lkind*deg2rad

  ! name files
  topog_file = "topography.data.nc"
  water_file = "water.data.nc"

  ! create data for both topog and water files
  ipts_r = 1 ; ipts = 2  ; iptsp1 = 3
  jpts_r = 1 ; jpts = 2  ; jptsp1 = 3


  xdat = (/1.0_lkind*deg2rad, 2.0_lkind*deg2rad, 3.0_lkind*deg2rad/)   !size of iptsp1, in radians
  ydat = (/1.0_lkind*deg2rad, 2.0_lkind*deg2rad, 3.0_lkind*deg2rad/)   !size of jptsp1, in radians

  zdat(1,1) = 0.0_lkind ; zdat(1,2) = 1.0_lkind
  zdat(2,1) = 1.0_lkind ; zdat(2,2) = 0.0_lkind   !size of (ipts, jpts)



  ! write topog file
  if (open_file(top_fileobj, topog_file, "overwrite")) then
    call register_axis(top_fileobj, "i_zdat", ipts)   !first index dimension in zdat
    call register_axis(top_fileobj, "j_zdat", jpts)   !second index dimension in zdat
    call register_axis(top_fileobj, "ipts_r", ipts_r)
    call register_axis(top_fileobj, "jpts_r", jpts_r)
    call register_axis(top_fileobj, "i_xdat", iptsp1) !# of points in xdat variable
    call register_axis(top_fileobj, "j_ydat", jptsp1) !# of point in ydat variable

    call register_field(top_fileobj, "ipts",   "double", dimensions=(/"ipts_r"/))
    call register_field(top_fileobj, "jpts",   "double", dimensions=(/"jpts_r"/))
    call register_field(top_fileobj, "xdat",   "double", dimensions=(/"i_xdat"/))
    call register_field(top_fileobj, "ydat",   "double", dimensions=(/"j_ydat"/))
    call register_field(top_fileobj, "zdat",   "double", dimensions=(/"i_zdat", "j_zdat"/))

    call write_data(top_fileobj, "ipts",   real(ipts, TEST_TOP_KIND_))
    call write_data(top_fileobj, "jpts",   real(jpts, TEST_TOP_KIND_))
    call write_data(top_fileobj, "xdat",   xdat)
    call write_data(top_fileobj, "ydat",   ydat)
    call write_data(top_fileobj, "zdat",   zdat)

    call close_file(top_fileobj)

  else
    call mpp_error(FATAL, "test_topography: error opening topog_file")
  end if

  ! write water file
  if (open_file(top_fileobj, water_file, "overwrite")) then
    call register_axis(top_fileobj, "i_zdat", ipts)   !first index dimension in zdat
    call register_axis(top_fileobj, "j_zdat", jpts)   !second index dimension in zdat
    call register_axis(top_fileobj, "ipts_r", ipts_r)
    call register_axis(top_fileobj, "jpts_r", jpts_r)
    call register_axis(top_fileobj, "i_xdat", iptsp1) !# of points in xdat variable
    call register_axis(top_fileobj, "j_ydat", jptsp1) !# of point in ydat variable

    call register_field(top_fileobj, "ipts",   "double", dimensions=(/"ipts_r"/))
    call register_field(top_fileobj, "jpts",   "double", dimensions=(/"jpts_r"/))
    call register_field(top_fileobj, "xdat",   "double", dimensions=(/"i_xdat"/))
    call register_field(top_fileobj, "ydat",   "double", dimensions=(/"j_ydat"/))
    call register_field(top_fileobj, "zdat",   "double", dimensions=(/"i_zdat", "j_zdat"/))

    call write_data(top_fileobj, "ipts",   real(ipts, TEST_TOP_KIND_))
    call write_data(top_fileobj, "jpts",   real(jpts, TEST_TOP_KIND_))
    call write_data(top_fileobj, "xdat",   xdat)
    call write_data(top_fileobj, "ydat",   ydat)
    call write_data(top_fileobj, "zdat",   zdat)

    call close_file(top_fileobj)

  else
    call mpp_error(FATAL, "test_topography: error opening water_file")
  end if

  call test_topog_mean
  call test_topog_stdev
  call test_get_ocean_frac
  call test_get_ocean_mask
  call test_get_water_frac
  call test_get_water_mask

  call fms_end

  contains

  subroutine test_topog_mean() 

    implicit none
    real(kind=TEST_TOP_KIND_), dimension(size(lon2d,1)-1,size(lat2d,2)-1) :: zmean2d
    real(kind=TEST_TOP_KIND_), dimension(size(lon1d)-1,size(lat1d)-1)     :: zmean1d
    logical                                                               :: get_mean_answer

    !---------------------------------------- test topog mean 2d ---------------------------------------------!

    get_mean_answer = get_topog_mean(lon2d, lat2d, zmean2d)

    if (get_mean_answer .neqv. .true.) call mpp_error(FATAL, "topog field not read correctly")
    call check_answers(zmean2d(1,1), 0.5_lkind, "Error in test_topog_mean 2d")

    !---------------------------------------- test topog mean 1d ---------------------------------------------!

    get_mean_answer = get_topog_mean(lon1d, lat1d, zmean1d)

    if (get_mean_answer .neqv. .true.) call mpp_error(FATAL, "topog field not read correctly")
    call check_answers(zmean1d(1,1), 0.5_lkind, "Error in test_topog_mean 1d")

  end subroutine test_topog_mean

  subroutine test_topog_stdev

    implicit none
    real(kind=TEST_TOP_KIND_), dimension(size(lon2d,1)-1,size(lat2d,2)-1) :: stdev2d
    real(kind=TEST_TOP_KIND_), dimension(size(lon1d)-1,size(lat1d)-1)     :: stdev1d
    logical                                                               :: get_stdev_answer

    !---------------------------------------- test topog stdev 2d ---------------------------------------------!

    get_stdev_answer = get_topog_stdev(lon2d, lat2d, stdev2d)

    if (get_stdev_answer .neqv. .true.) call mpp_error(FATAL, "topog field not read correctly")
    call check_answers(stdev2d(1,1), 0.5_lkind, "Error in test_topog_stdev 2d")

    !---------------------------------------- test topog stdev 2d ---------------------------------------------!

    get_stdev_answer = get_topog_stdev(lon1d, lat1d, stdev1d)

    if (get_stdev_answer .neqv. .true.) call mpp_error(FATAL, "topog field not read correctly")
    call check_answers(stdev1d(1,1), 0.5_lkind, "Error in test_topog_stdev 1d")

  end subroutine test_topog_stdev

  subroutine test_get_ocean_frac

    implicit none
    real(kind=TEST_TOP_KIND_), dimension(size(lon2d,1)-1,size(lat2d,2)-1) :: ocean_frac2d
    real(kind=TEST_TOP_KIND_), dimension(size(lon1d)-1,size(lat1d)-1)     :: ocean_frac1d
    logical                                                               :: get_ocean_frac_answer

    !---------------------------------------- test get_ocean_frac 2d ---------------------------------------------!

    get_ocean_frac_answer = get_ocean_frac(lon2d, lat2d, ocean_frac2d)

    if (get_ocean_frac_answer .neqv. .true.) call mpp_error(FATAL, "ocean field not read correctly")
    call check_answers(ocean_frac2d(1,1), 0.5_lkind, "Error in test_get_ocean_frac 2d")

    !---------------------------------------- test get_ocean_frac 1d ---------------------------------------------!

    get_ocean_frac_answer = get_ocean_frac(lon1d, lat1d, ocean_frac1d)

    if (get_ocean_frac_answer .neqv. .true.) call mpp_error(FATAL, "ocean field not read correctly")
    call check_answers(ocean_frac1d(1,1), 0.5_lkind, "Error in test_get_ocean_frac 1d")

  end subroutine test_get_ocean_frac

  subroutine test_get_ocean_mask

    implicit none
    logical, dimension(size(lon2d,1)-1,size(lat2d,2)-1) :: ocean_mask2d
    logical, dimension(size(lon1d)-1,size(lat1d)-1)     :: ocean_mask1d
    logical                                             :: get_ocean_mask_answer

    !---------------------------------------- test get_ocean_mask 2d ---------------------------------------------!

    get_ocean_mask_answer = get_ocean_mask(lon2d, lat2d, ocean_mask2d)


    if (get_ocean_mask_answer .neqv. .true.) call mpp_error(FATAL, "ocean field not read correctly")
    if (ocean_mask2d(1,1) .neqv. .false.) call mpp_error(FATAL, "test_get_ocean_mask 2d: ocean mask should be false")

    !---------------------------------------- test get_ocean_mask 1d ---------------------------------------------!

    get_ocean_mask_answer = get_ocean_mask(lon1d, lat1d, ocean_mask1d)

    if (get_ocean_mask_answer .neqv. .true.) call mpp_error(FATAL, "ocean field not read correctly")
    if (ocean_mask1d(1,1) .neqv. .false.) call mpp_error(FATAL, "test_get_ocean_mask 1d: ocean mask should be false")

  end subroutine test_get_ocean_mask

  subroutine test_get_water_frac

    implicit none
    real(kind=TEST_TOP_KIND_), dimension(size(lon2d,1)-1,size(lat2d,2)-1) :: water_frac2d
    real(kind=TEST_TOP_KIND_), dimension(size(lon1d)-1,size(lat1d)-1)     :: water_frac1d
    logical                                                               :: get_water_frac_answer

    !---------------------------------------- test get_water_frac 2d ---------------------------------------------!

    get_water_frac_answer = get_water_frac(lon2d, lat2d, water_frac2d)

    if (get_water_frac_answer .neqv. .true.) call mpp_error(FATAL, "ocean field not read correctly")
    call check_answers(water_frac2d(1,1), 0.5_lkind, "Error in test_get_water_frac 2d")

    !---------------------------------------- test get_water_frac 1d ---------------------------------------------!

    get_water_frac_answer = get_water_frac(lon1d, lat1d, water_frac1d)

    if (get_water_frac_answer .neqv. .true.) call mpp_error(FATAL, "ocean field not read correctly")
    call check_answers(water_frac1d(1,1), 0.5_lkind, "Error in test_get_ocean_frac 1d")

  end subroutine test_get_water_frac

  subroutine test_get_water_mask

    implicit none
    logical, dimension(size(lon2d,1)-1,size(lat2d,2)-1) :: water_mask2d
    logical, dimension(size(lon1d)-1,size(lat1d)-1)     :: water_mask1d
    logical                                             :: get_water_mask_answer

    !---------------------------------------- test get_water_mask 2d ---------------------------------------------!

    get_water_mask_answer = get_water_mask(lon2d, lat2d, water_mask2d)

    if (get_water_mask_answer .neqv. .true.) call mpp_error(FATAL, "ocean field not read correctly")
    if (water_mask2d(1,1) .neqv. .false.) call mpp_error(FATAL, "test_get_water_mask 2d: ocean mask should be false")

    !---------------------------------------- test get_water_mask 1d ---------------------------------------------!

    get_water_mask_answer = get_ocean_mask(lon1d, lat1d, water_mask1d)

    if (get_water_mask_answer .neqv. .true.) call mpp_error(FATAL, "ocean field not read correctly")
    if (water_mask1d(1,1) .neqv. .false.) call mpp_error(FATAL, "test_get_ocean_mask 1d: ocean mask should be false")

  end subroutine test_get_water_mask

  subroutine check_answers(result, answer, what_error)
    
    implicit none
    real(kind=TEST_TOP_KIND_) :: result       ! value calculated from script
    real(kind=TEST_TOP_KIND_) :: answer       ! expected answer
    character(*)              :: what_error   ! error message to print

    if (result .ne. answer) then
      write(*,*) 'Expected ', answer, ' but computed ', result
      call mpp_error(FATAL, trim(what_error))
    end if

  end subroutine check_answers



end program test_top