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
!> @brief Unit test for astronomy/diurnal_solar interfaces
!> @email gfdl.climate.model.info@noaa.gov
!> @description Performs calculations done in astronomy_mod using the diurnal_solar
!> interfaces using 32 and 64 bit reals
!! TODO: Add more comprehensive testing with optional arguments for any dirunal
!! solar interface including dt, allow_negative_cosz, and half_day_out
!! TODO: Add more comprehensive testing changing angle and declination values
!! TODO: XFAILS for out of bounds argument values

program test_diurnal_solar

  use astronomy_mod
  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only: mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use time_manager_mod,  only: JULIAN, set_calendar_type, time_type, set_date
  use constants_mod,     only: PI
  use platform_mod,      only: r4_kind, r8_kind

  implicit none

  real(kind=r8_kind), parameter :: twopi = 2.0_r8_kind * real(PI,r8_kind)
  real(kind=r8_kind), parameter :: deg_to_rad = twopi / 360.0_r8_kind


  real(kind=r8_kind)   :: ecc   = 0.01671_r8_kind  !< Eccentricity of Earth's orbit [dimensionless]
  real(kind=r8_kind)   :: obliq = 23.439_r8_kind   !< Obliquity [degrees]
  real(kind=r8_kind)   :: per   = 102.932_r8_kind  !< Longitude of perihelion with respect


  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init

  call fms_end()

  call test_diurnal_solar_allow_neg_2d
  call test_diurnal_solar_allow_neg_1d
  call test_diurnal_solar_allow_neg_0d

  call test_diurnal_solar_no_neg_2d
  call test_diurnal_solar_no_neg_1d
  call test_diurnal_solar_no_neg_0d

  contains

  ! test for diurnal_solar inteface below, reference values calcualted underneath tests
  subroutine test_diurnal_solar_allow_neg_2d
    !! for this test, only the routine input variables lat and lon
    !! will be changed since these are not dependent on private function calculations,
    !! a more comprehensive test that includes different values for the 
    !! time_since_ae, orbit angle and declination should be created
    !! this test will allow negative cosz (allow_negative_cosz=.true.)
    implicit none
    integer, parameter                        :: n = 4
    real(kind=TEST_AST_KIND_), dimension(n,n) :: lat_2d, lon_2d, h
    real(kind=TEST_AST_KIND_)                 :: gmt, time_since_ae
    real(kind=TEST_AST_KIND_), dimension(n,n) :: ref_cosz_2d, cosz_2d, ref_fracday_2d, fracday_2d
    real(kind=TEST_AST_KIND_)                 :: ref_rrsun, rrsun
    real(kind=TEST_AST_KIND_)                 :: ang, dec
    integer                                   :: i, j, counter
    integer, parameter                        :: lkind = TEST_AST_KIND_

    ! test only changes in lat
    lon_2d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lat_2d = 0.0_lkind
    counter = 1
    do i = 1, 4
      do j = 1, 4
        lat_2d(j,i) = real(counter, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
        counter = counter + 1
      end do
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, dec)
    call create_ref_cosz_allow_neg_2d(lat_2d, lon_2d, gmt, h, dec, ref_cosz_2d, ref_fracday_2d)

    call diurnal_solar(lat_2d, lon_2d, gmt, time_since_ae, cosz_2d, fracday_2d, rrsun, allow_negative_cosz=.true.)

    do i = 1, 4
      do j = 1, 4

        if (ref_cosz_2d(i,j) .ne. cosz_2d(i,j)) then
        print *, ref_cosz_2d(i,j), " @ ", i,j, " does not equal ", cosz_2d(i,j)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 2d cosz value, with only changes in lat, does not match reference value")
        end if

        if (ref_fracday_2d(i,j) .ne. fracday_2d(i,j)) then
          print *, ref_fracday_2d(i,j), " @ ", i,j, " does not equal ", fracday_2d(i,j)
          call mpp_error(FATAL, &
          "test_diurnal_solar: 2d fracday value, with only changes in lat, does not match reference value")
        end if

      end do
    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " @ ", i,j, " does not equal ", rrsun
      call mpp_error(FATAL, &
      "test_diurnal_solar: 2d rrsun value, with only changes in lat, does not match reference value")
    end if

    ! test only changes in lon
    lat_2d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lon_2d = 0.0_lkind
    counter = 1
    do i = 1, 4
      do j = 1, 4
        lon_2d(j,i) = real(counter, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
        counter = counter + 1
      end do
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, dec)
    call create_ref_cosz_allow_neg_2d(lat_2d, lon_2d, gmt, h, dec, ref_cosz_2d, ref_fracday_2d)

    call diurnal_solar(lat_2d, lon_2d, gmt, time_since_ae, cosz_2d, fracday_2d, rrsun, allow_negative_cosz=.true.)

    do i = 1, 4
      do j = 1, 4

        if (ref_cosz_2d(i,j) .ne. cosz_2d(i,j)) then
          print *, ref_cosz_2d(i,j), " @ ", i,j, " does not equal ", cosz_2d(i,j)
          call mpp_error(FATAL, &
          "test_diurnal_solar: 2d cosz value, with only changes in lon, does not match reference value")
        end if

        if (ref_fracday_2d(i,j) .ne. fracday_2d(i,j)) then
          print *, ref_fracday_2d(i,j), " @ ", i, " does not equal ", fracday_2d(i,j)
          call mpp_error(FATAL, &
          "test_diurnal_solar: 2d fracday value, with only changes in lon, does not match reference value")
        end if

      end do
    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " @ ", i,j, " does not equal ", rrsun
      call mpp_error(FATAL, &
      "test_diurnal_solar: 2d rrsun value, with only changes in lon, does not match reference value")
    end if

  end subroutine test_diurnal_solar_allow_neg_2d

  subroutine test_diurnal_solar_allow_neg_1d
    !! for this tests, only the routine input variables lat, lon, and gmt
    !! will be changed since these are not dependent on private function calculations,
    !! a more comprehensive test that includes different values for the 
    !! time_since_ae, orbit angle and declination should be created
    !! this test will allow negative cosz (allow_negative_cosz=.true.)
    implicit none
    integer, parameter                      :: n = 16
    real(kind=TEST_AST_KIND_), dimension(n) :: lat_1d, lon_1d, h
    real(kind=TEST_AST_KIND_)               :: gmt, time_since_ae
    real(kind=TEST_AST_KIND_), dimension(n) :: ref_cosz_1d, cosz_1d, ref_fracday_1d, fracday_1d
    real(kind=TEST_AST_KIND_)               :: ref_rrsun, rrsun
    real(kind=TEST_AST_KIND_)               :: ang, dec
    integer                                 :: i
    integer, parameter                      :: lkind = TEST_AST_KIND_

    ! test only changes in lat
    lon_1d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lat_1d = 0.0_lkind
    do i = 1, 16
      lat_1d(i) = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, dec)
    call create_ref_cosz_allow_neg_1d(lat_1d, lon_1d, gmt, h, dec, ref_cosz_1d, ref_fracday_1d)

    call diurnal_solar(lat_1d, lon_1d, gmt, time_since_ae, cosz_1d, fracday_1d, rrsun, allow_negative_cosz=.true.)

    do i = 1, 16

      if (ref_cosz_1d(i) .ne. cosz_1d(i)) then
        print *, ref_cosz_1d(i), " @ ", i, " does not equal ", cosz_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d cosz value, with only changes in lat, does not match reference value")
      end if

      if (ref_fracday_1d(i) .ne. fracday_1d(i)) then
        print *, ref_fracday_1d(i), " @ ", i, " does not equal ", fracday_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d fracday value, with only changes in lat, does not match reference value")
      end if

    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
      call mpp_error(FATAL, &
      "test_diurnal_solar: 1d rrsun value, with only changes in lat, does not match reference value")
    end if

    ! test only changes in lon
    lat_1d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lon_1d = 0.0_lkind
    do i = 1, 16
      lon_1d(i) = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, dec)
    call create_ref_cosz_allow_neg_1d(lat_1d, lon_1d, gmt, h, dec, ref_cosz_1d, ref_fracday_1d)

    call diurnal_solar(lat_1d, lon_1d, gmt, time_since_ae, cosz_1d, fracday_1d, rrsun, allow_negative_cosz=.true.)

    do i = 1, 16

      if (ref_cosz_1d(i) .ne. cosz_1d(i)) then
        print *, ref_cosz_1d(i), " @ ", i, " does not equal ", cosz_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d cosz value, with only changes in lon, does not match reference value")
      end if

      if (ref_fracday_1d(i) .ne. fracday_1d(i)) then
        print *, ref_fracday_1d(i), " @ ", i, " does not equal ", fracday_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d fracday value, with only changes in lon, does not match reference value")
      end if

    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
      call mpp_error(FATAL, &
      "test_diurnal_solar: 1d rrsun value, with only changes in lon, does not match reference value")
    end if

    ! test only changes in time of day (gmt)
    lat_1d = 0.0_lkind ; lon_1d = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    gmt = 0.0_lkind
    do i = 1, 16
      gmt = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind

      call create_ref_rrsun(ang, ref_rrsun)
      call create_ref_declination(ang, dec)
      call create_ref_cosz_allow_neg_1d(lat_1d, lon_1d, gmt, h, dec, ref_cosz_1d, ref_fracday_1d)

      call diurnal_solar(lat_1d, lon_1d, gmt, time_since_ae, cosz_1d, fracday_1d, rrsun, allow_negative_cosz=.true.)

      if (ref_cosz_1d(i) .ne. cosz_1d(i)) then
        print *, ref_cosz_1d(i), " @ ", i, " does not equal ", cosz_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d cosz value, with only changes in gmt, does not match reference value")
      end if

      if (ref_fracday_1d(i) .ne. fracday_1d(i)) then
        print *, ref_fracday_1d(i), " @ ", i, " does not equal ", fracday_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d fracday value, with only changes in gmt, does not match reference value")
      end if

      if (ref_rrsun .ne. rrsun) then
        print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d rrsun value, with only changes in gmt, does not match reference value")
      end if
    end do

  end subroutine test_diurnal_solar_allow_neg_1d

  subroutine test_diurnal_solar_allow_neg_0d
    !! for this test, only the routine input variables lat, lon, and gmt
    !! will be changed since these are not dependent on private function calculations,
    !! a more comprehensive test that includes different values for the 
    !! time_since_ae, orbit angle and declination should be created
    !! this test will allow negative cosz (allow_negative_cosz=.true.)
    implicit none
    real(kind=TEST_AST_KIND_) :: lat_0d, lon_0d, h
    real(kind=TEST_AST_KIND_) :: gmt, time_since_ae
    real(kind=TEST_AST_KIND_) :: ref_cosz_0d, cosz_0d, ref_fracday_0d, fracday_0d
    real(kind=TEST_AST_KIND_) :: ref_rrsun, rrsun
    real(kind=TEST_AST_KIND_) :: ang, dec
    integer                   :: i
    integer, parameter        :: lkind = TEST_AST_KIND_

    ! test only changes in latitude 
    lon_0d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lat_0d = 0.0_lkind
    do i = 1, 16
      lat_0d = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind

      call create_ref_rrsun(ang, ref_rrsun)
      call create_ref_declination(ang, dec)
      call create_ref_cosz_allow_neg_0d(lat_0d, lon_0d, gmt, h, dec, ref_cosz_0d, ref_fracday_0d)

      call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, fracday_0d, rrsun, allow_negative_cosz=.true.)

      if (ref_cosz_0d .ne. cosz_0d) then
        print *, ref_cosz_0d, " @ ", i, " does not equal ", cosz_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d cosz value, with changes in lat, does not match reference value")
      end if

      if (ref_fracday_0d .ne. fracday_0d) then
        print *, ref_fracday_0d, " @ ", i, " does not equal ", fracday_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d fracday value, with changes in lat, does not match reference value")
      end if

      if (ref_rrsun .ne. rrsun) then
        print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d rrsun value, with changes in lat, does not match reference value")
      end if
    end do

    ! test only changes in longitude 
    lat_0d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lon_0d = 0.0_lkind
    do i = 1, 16
      lon_0d = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind

      call create_ref_rrsun(ang, ref_rrsun)
      call create_ref_declination(ang, dec)
      call create_ref_cosz_allow_neg_0d(lat_0d, lon_0d, gmt, h, dec, ref_cosz_0d, ref_fracday_0d)

      call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, fracday_0d, rrsun, allow_negative_cosz=.true.)

      if (ref_cosz_0d .ne. cosz_0d) then
        print *, ref_cosz_0d, " @ ", i, " does not equal ", cosz_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d cosz value, with changes in lon, does not match reference value")
      end if

      if (ref_fracday_0d .ne. fracday_0d) then
        print *, ref_fracday_0d, " @ ", i, " does not equal ", fracday_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d fracday value, with changes in lon, does not match reference value")
      end if

      if (ref_rrsun .ne. rrsun) then
        print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d rrsun value, with changes in lon, does not match reference value")
      end if
    end do

    ! test only changes in time of day (gmt)
    lat_0d = 0.0_lkind ; lon_0d = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    gmt = 0.0_lkind
    do i = 1, 16
      gmt = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind

      call create_ref_rrsun(ang, ref_rrsun)
      call create_ref_declination(ang, dec)
      call create_ref_cosz_allow_neg_0d(lat_0d, lon_0d, gmt, h, dec, ref_cosz_0d, ref_fracday_0d)

      call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, fracday_0d, rrsun, allow_negative_cosz=.true.)

      if (ref_cosz_0d .ne. cosz_0d) then
        print *, ref_cosz_0d, " @ ", i, " does not equal ", cosz_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d cosz value, with changes in gmt, does not match reference value")
      end if

      if (ref_fracday_0d .ne. fracday_0d) then
        print *, ref_fracday_0d, " @ ", i, " does not equal ", fracday_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d fracday value, with changes in gmt, does not match reference value")
      end if

      if (ref_rrsun .ne. rrsun) then
        print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d rrsun value, with changes in gmt, does not match reference value")
      end if
    end do

  end subroutine test_diurnal_solar_allow_neg_0d

  subroutine test_diurnal_solar_no_neg_2d
    !! for this test, only the routine input variables lat, lon, and gmt
    !! will be changed since these are not dependent on private function calculations,
    !! a more comprehensive test that includes different values for the 
    !! time_since_ae, orbit angle and declination should be created
    !! this test will not allow negative cosz and will explicitly have
    !! allow_negative_cosz=.false. although omitting the optional argument 
    !! should produce the same effect
    implicit none
    integer, parameter                        :: n = 4
    real(kind=TEST_AST_KIND_), dimension(n,n) :: lat_2d, lon_2d, h
    real(kind=TEST_AST_KIND_)                 :: gmt, time_since_ae
    real(kind=TEST_AST_KIND_), dimension(n,n) :: ref_cosz_2d, cosz_2d, ref_fracday_2d, fracday_2d
    real(kind=TEST_AST_KIND_)                 :: ref_rrsun, rrsun
    real(kind=TEST_AST_KIND_)                 :: ang, dec
    integer                                   :: i, j, counter
    integer, parameter                        :: lkind = TEST_AST_KIND_

    ! test only changes in lat
    lon_2d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lat_2d = 0.0_lkind
    counter = 1
    do i = 1, 4
      do j = 1, 4
        lat_2d(j,i) = real(counter, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
        counter = counter + 1
      end do
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, dec)
    call create_ref_cosz_no_neg_2d(lat_2d, lon_2d, gmt, h, dec, ref_cosz_2d, ref_fracday_2d)

    call diurnal_solar(lat_2d, lon_2d, gmt, time_since_ae, cosz_2d, fracday_2d, rrsun, allow_negative_cosz=.false.)

    do i = 1, 4
      do j = 1, 4

        if (ref_cosz_2d(i,j) .ne. cosz_2d(i,j)) then
        print *, ref_cosz_2d(i,j), " @ ", i,j, " does not equal ", cosz_2d(i,j)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 2d cosz value, with only changes in lat, does not match reference value")
        end if

        if (ref_fracday_2d(i,j) .ne. fracday_2d(i,j)) then
          print *, ref_fracday_2d(i,j), " @ ", i,j, " does not equal ", fracday_2d(i,j)
          call mpp_error(FATAL, &
          "test_diurnal_solar: 2d fracday value, with only changes in lat, does not match reference value")
        end if

      end do
    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " @ ", i,j, " does not equal ", rrsun
      call mpp_error(FATAL, &
      "test_diurnal_solar: 2d rrsun value, with only changes in lat, does not match reference value")
    end if

    ! test only changes in lon
    lat_2d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lon_2d = 0.0_lkind
    counter = 1
    do i = 1, 4
      do j = 1, 4
        lon_2d(j,i) = real(counter, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
        counter = counter + 1
      end do
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, dec)
    call create_ref_cosz_no_neg_2d(lat_2d, lon_2d, gmt, h, dec, ref_cosz_2d, ref_fracday_2d)

    call diurnal_solar(lat_2d, lon_2d, gmt, time_since_ae, cosz_2d, fracday_2d, rrsun, allow_negative_cosz=.false.)

    do i = 1, 4
      do j = 1, 4

        if (ref_cosz_2d(i,j) .ne. cosz_2d(i,j)) then
          print *, ref_cosz_2d(i,j), " @ ", i,j, " does not equal ", cosz_2d(i,j)
          call mpp_error(FATAL, &
          "test_diurnal_solar: 2d cosz value, with only changes in lon, does not match reference value")
        end if

        if (ref_fracday_2d(i,j) .ne. fracday_2d(i,j)) then
          print *, ref_fracday_2d(i,j), " @ ", i, " does not equal ", fracday_2d(i,j)
          call mpp_error(FATAL, &
          "test_diurnal_solar: 2d fracday value, with only changes in lon, does not match reference value")
        end if

      end do
    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " @ ", i,j, " does not equal ", rrsun
      call mpp_error(FATAL, &
      "test_diurnal_solar: 2d rrsun value, with only changes in lon, does not match reference value")
    end if

  end subroutine test_diurnal_solar_no_neg_2d

  subroutine test_diurnal_solar_no_neg_1d
    !! for this tests, only the routine input variables lat, lon, and gmt
    !! will be changed since these are not dependent on private function calculations,
    !! a more comprehensive test that includes different values for the 
    !! time_since_ae, orbit angle and declination should be created
    !! this test will not allow negative cosz and will explicitly have
    !! allow_negative_cosz=.false. although omitting the optional argument 
    !! should produce the same effect

    implicit none
    integer, parameter                      :: n = 16
    real(kind=TEST_AST_KIND_), dimension(n) :: lat_1d, lon_1d, h
    real(kind=TEST_AST_KIND_)               :: gmt, time_since_ae
    real(kind=TEST_AST_KIND_), dimension(n) :: ref_cosz_1d, cosz_1d, ref_fracday_1d, fracday_1d
    real(kind=TEST_AST_KIND_)               :: ref_rrsun, rrsun
    real(kind=TEST_AST_KIND_)               :: ang, dec
    integer                                 :: i
    integer, parameter                      :: lkind = TEST_AST_KIND_

    ! test only changes in lat
    lon_1d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lat_1d = 0.0_lkind
    do i = 1, 16
      lat_1d(i) = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, dec)
    call create_ref_cosz_no_neg_1d(lat_1d, lon_1d, gmt, h, dec, ref_cosz_1d, ref_fracday_1d)

    call diurnal_solar(lat_1d, lon_1d, gmt, time_since_ae, cosz_1d, fracday_1d, rrsun, allow_negative_cosz=.false.)

    do i = 1, 16

      if (ref_cosz_1d(i) .ne. cosz_1d(i)) then
        print *, ref_cosz_1d(i), " @ ", i, " does not equal ", cosz_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d cosz value, with only changes in lat, does not match reference value")
      end if

      if (ref_fracday_1d(i) .ne. fracday_1d(i)) then
        print *, ref_fracday_1d(i), " @ ", i, " does not equal ", fracday_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d fracday value, with only changes in lat, does not match reference value")
      end if

    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
      call mpp_error(FATAL, &
      "test_diurnal_solar: 1d rrsun value, with only changes in lat, does not match reference value")
    end if

    ! test only changes in lon
    lat_1d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lon_1d = 0.0_lkind
    do i = 1, 16
      lon_1d(i) = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind
    end do

    call create_ref_rrsun(ang, ref_rrsun)
    call create_ref_declination(ang, dec)
    call create_ref_cosz_no_neg_1d(lat_1d, lon_1d, gmt, h, dec, ref_cosz_1d, ref_fracday_1d)

    call diurnal_solar(lat_1d, lon_1d, gmt, time_since_ae, cosz_1d, fracday_1d, rrsun, allow_negative_cosz=.false.)

    do i = 1, 16

      if (ref_cosz_1d(i) .ne. cosz_1d(i)) then
        print *, ref_cosz_1d(i), " @ ", i, " does not equal ", cosz_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d cosz value, with only changes in lon, does not match reference value")
      end if

      if (ref_fracday_1d(i) .ne. fracday_1d(i)) then
        print *, ref_fracday_1d(i), " @ ", i, " does not equal ", fracday_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d fracday value, with only changes in lon, does not match reference value")
      end if

    end do

    if (ref_rrsun .ne. rrsun) then
      print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
      call mpp_error(FATAL, &
      "test_diurnal_solar: 1d rrsun value, with only changes in lon, does not match reference value")
    end if

    ! test only changes in time of day (gmt)
    lat_1d = 0.0_lkind ; lon_1d = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    gmt = 0.0_lkind
    do i = 1, 16
      gmt = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind

      call create_ref_rrsun(ang, ref_rrsun)
      call create_ref_declination(ang, dec)
      call create_ref_cosz_no_neg_1d(lat_1d, lon_1d, gmt, h, dec, ref_cosz_1d, ref_fracday_1d)

      call diurnal_solar(lat_1d, lon_1d, gmt, time_since_ae, cosz_1d, fracday_1d, rrsun, allow_negative_cosz=.false.)

      if (ref_cosz_1d(i) .ne. cosz_1d(i)) then
        print *, ref_cosz_1d(i), " @ ", i, " does not equal ", cosz_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d cosz value, with only changes in gmt, does not match reference value")
      end if

      if (ref_fracday_1d(i) .ne. fracday_1d(i)) then
        print *, ref_fracday_1d(i), " @ ", i, " does not equal ", fracday_1d(i)
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d fracday value, with only changes in gmt, does not match reference value")
      end if

      if (ref_rrsun .ne. rrsun) then
        print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
        call mpp_error(FATAL, &
        "test_diurnal_solar: 1d rrsun value, with only changes in gmt, does not match reference value")
      end if
    end do

  end subroutine test_diurnal_solar_no_neg_1d

  subroutine test_diurnal_solar_no_neg_0d
    !! for this test, only the routine input variables lat, lon, and gmt
    !! will be changed since these are not dependent on private function calculations,
    !! a more comprehensive test that includes different values for the 
    !! time_since_ae, orbit angle and declination should be created
    !! this test will not allow negative cosz and will explicitly have
    !! allow_negative_cosz=.false. although omitting the optional argument 
    !! should produce the same effect
    implicit none
    real(kind=TEST_AST_KIND_) :: lat_0d, lon_0d, h
    real(kind=TEST_AST_KIND_) :: gmt, time_since_ae
    real(kind=TEST_AST_KIND_) :: ref_cosz_0d, cosz_0d, ref_fracday_0d, fracday_0d
    real(kind=TEST_AST_KIND_) :: ref_rrsun, rrsun
    real(kind=TEST_AST_KIND_) :: ang, dec
    integer                   :: i
    integer, parameter        :: lkind = TEST_AST_KIND_

    ! test only changes in latitude 
    lon_0d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lat_0d = 0.0_lkind
    do i = 1, 16
      lat_0d = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind

      call create_ref_rrsun(ang, ref_rrsun)
      call create_ref_declination(ang, dec)
      call create_ref_cosz_no_neg_0d(lat_0d, lon_0d, gmt, h, dec, ref_cosz_0d, ref_fracday_0d)

      call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, fracday_0d, rrsun, allow_negative_cosz=.false.)

      if (ref_cosz_0d .ne. cosz_0d) then
        print *, ref_cosz_0d, " @ ", i, " does not equal ", cosz_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d cosz value, with changes in lat, does not match reference value")
      end if

      if (ref_fracday_0d .ne. fracday_0d) then
        print *, ref_fracday_0d, " @ ", i, " does not equal ", fracday_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d fracday value, with changes in lat, does not match reference value")
      end if

      if (ref_rrsun .ne. rrsun) then
        print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d rrsun value, with changes in lat, does not match reference value")
      end if
    end do

    ! test only changes in longitude 
    lat_0d = 0.0_lkind ; gmt = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    lon_0d = 0.0_lkind
    do i = 1, 16
      lon_0d = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind

      call create_ref_rrsun(ang, ref_rrsun)
      call create_ref_declination(ang, dec)
      call create_ref_cosz_no_neg_0d(lat_0d, lon_0d, gmt, h, dec, ref_cosz_0d, ref_fracday_0d)

      call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, fracday_0d, rrsun, allow_negative_cosz=.false.)

      if (ref_cosz_0d .ne. cosz_0d) then
        print *, ref_cosz_0d, " @ ", i, " does not equal ", cosz_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d cosz value, with changes in lon, does not match reference value")
      end if

      if (ref_fracday_0d .ne. fracday_0d) then
        print *, ref_fracday_0d, " @ ", i, " does not equal ", fracday_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d fracday value, with changes in lon, does not match reference value")
      end if

      if (ref_rrsun .ne. rrsun) then
        print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d rrsun value, with changes in lon, does not match reference value")
      end if
    end do

    ! test only changes in time of day (gmt)
    lat_0d = 0.0_lkind ; lon_0d = 0.0_lkind ; time_since_ae = 0.0_lkind
    ang = 0.0_lkind ; dec = 0.0_lkind ; h = real(PI, TEST_AST_KIND_)/2.0_lkind

    gmt = 0.0_lkind
    do i = 1, 16
      gmt = real(i, TEST_AST_KIND_) * real(PI,TEST_AST_KIND_)/8.0_lkind

      call create_ref_rrsun(ang, ref_rrsun)
      call create_ref_declination(ang, dec)
      call create_ref_cosz_no_neg_0d(lat_0d, lon_0d, gmt, h, dec, ref_cosz_0d, ref_fracday_0d)

      call diurnal_solar(lat_0d, lon_0d, gmt, time_since_ae, cosz_0d, fracday_0d, rrsun, allow_negative_cosz=.false.)

      if (ref_cosz_0d .ne. cosz_0d) then
        print *, ref_cosz_0d, " @ ", i, " does not equal ", cosz_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d cosz value, with changes in gmt, does not match reference value")
      end if

      if (ref_fracday_0d .ne. fracday_0d) then
        print *, ref_fracday_0d, " @ ", i, " does not equal ", fracday_0d
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d fracday value, with changes in gmt, does not match reference value")
      end if

      if (ref_rrsun .ne. rrsun) then
        print *, ref_rrsun, " @ ", i, " does not equal ", rrsun
        call mpp_error(FATAL, &
        "test_diurnal_solar: 0d rrsun value, with changes in gmt, does not match reference value")
      end if
    end do
  end subroutine test_diurnal_solar_no_neg_0d


  ! create all reference values
  subroutine create_ref_rrsun(ang, rrsun)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)  :: ang
    real(kind=TEST_AST_KIND_), intent(out) :: rrsun
    real(kind=TEST_AST_KIND_)              :: rad_per, rr
    integer, parameter                     :: lkind = TEST_AST_KIND_
  
    rad_per = real(per, TEST_AST_KIND_) * real(deg_to_rad, TEST_AST_KIND_)
    rr      = (1.0_lkind - real((ecc**2),TEST_AST_KIND_))/(1.0_lkind + real(ecc, TEST_AST_KIND_) * cos(ang - rad_per))
    rrsun  = rr**(-2)
  
  end subroutine create_ref_rrsun
  
  subroutine create_ref_declination(ang, dec)
  
    implicit none
    real(kind=TEST_AST_KIND_), intent(in)  :: ang
    real(kind=TEST_AST_KIND_), intent(out) :: dec
    real(kind=TEST_AST_KIND_)              :: rad_obliq, sin_dec
  
    rad_obliq    = real(obliq, TEST_AST_KIND_) * real(deg_to_rad, TEST_AST_KIND_)
    sin_dec      = - sin(rad_obliq)*sin(ang)
    dec          = asin(sin_dec)
  
  end subroutine create_ref_declination

  subroutine create_ref_cosz_allow_neg_2d(lat_2d, lon_2d, gmt, h, dec, cosz_2d, fracday_2d)

    implicit none
    integer, parameter :: n = 4
    real(kind=TEST_AST_KIND_), intent(in), dimension(n,n)  :: lat_2d, lon_2d, h
    real(kind=TEST_AST_KIND_), intent(in)                  :: gmt, dec
    real(kind=TEST_AST_KIND_), intent(out), dimension(n,n) :: cosz_2d, fracday_2d
    real(kind=TEST_AST_KIND_), dimension(n,n)              :: aa, bb, t
    integer, parameter                                     :: lkind = TEST_AST_KIND_


    !! these reference values will allow negative cosz values,
    !! i.e. Lallow_nagative = .true.
    !! they will also not include time averaging, i.e. .not. present(dt)
    !! fracday values calculated here

    aa = sin(lat_2d)*sin(dec)
    bb = cos(lat_2d)*cos(dec)

    t = gmt + lon_2d - real(PI, TEST_AST_KIND_)
    where (t >= real(PI, TEST_AST_KIND_))  t = t - real(twopi, TEST_AST_KIND_)
    where (t < real(-PI, TEST_AST_KIND_))  t = t + real(twopi, TEST_AST_KIND_)

    cosz_2d = aa + bb * cos(t)

    where (abs(t) < h)
      fracday_2d = 1.0_lkind
    else where
      fracday_2d = 0.0_lkind
    end where

  end subroutine create_ref_cosz_allow_neg_2d

  subroutine create_ref_cosz_allow_neg_1d(lat_1d, lon_1d, gmt, h, dec, cosz_1d, fracday_1d)

    implicit none
    integer, parameter :: n = 16
    real(kind=TEST_AST_KIND_), intent(in), dimension(n)  :: lat_1d, lon_1d, h
    real(kind=TEST_AST_KIND_), intent(in)                :: gmt, dec
    real(kind=TEST_AST_KIND_), intent(out), dimension(n) :: cosz_1d, fracday_1d
    real(kind=TEST_AST_KIND_), dimension(n)              :: aa, bb, t
    integer, parameter                                   :: lkind = TEST_AST_KIND_


    !! these reference values will allow negative cosz values,
    !! i.e. Lallow_nagative = .true.
    !! they will also not include time averaging, i.e. .not. present(dt)
    !! fracday values calculated here

    aa = sin(lat_1d)*sin(dec)
    bb = cos(lat_1d)*cos(dec)

    t = gmt + lon_1d - real(PI, TEST_AST_KIND_)
    where (t >= real(PI, TEST_AST_KIND_))  t = t - real(twopi, TEST_AST_KIND_)
    where (t < real(-PI, TEST_AST_KIND_))  t = t + real(twopi, TEST_AST_KIND_)

    cosz_1d = aa + bb * cos(t)

    where (abs(t) < h)
      fracday_1d = 1.0_lkind
    else where
      fracday_1d = 0.0_lkind
    end where

  end subroutine create_ref_cosz_allow_neg_1d

  subroutine create_ref_cosz_allow_neg_0d(lat_0d, lon_0d, gmt, h, dec, cosz_0d, fracday_0d)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)  :: lat_0d, lon_0d, gmt, h, dec
    real(kind=TEST_AST_KIND_), intent(out) :: cosz_0d, fracday_0d
    real(kind=TEST_AST_KIND_)              :: aa, bb, t
    integer, parameter                     :: lkind = TEST_AST_KIND_

    !! these reference values will allow negative cosz values,
    !! i.e. Lallow_nagative = .true.
    !! they will also not include time averaging, i.e. .not. present(dt)
    !! fracday values calculated here

    aa = sin(lat_0d)*sin(dec)
    bb = cos(lat_0d)*cos(dec)

    t = gmt + lon_0d - real(PI, TEST_AST_KIND_)
    if (t >= real(PI, TEST_AST_KIND_))  t = t - real(twopi, TEST_AST_KIND_)
    if (t < real(-PI, TEST_AST_KIND_))  t = t + real(twopi, TEST_AST_KIND_)

    cosz_0d = aa + bb * cos(t)

    if (abs(t) < h) then
      fracday_0d = 1.0_lkind
    else
      fracday_0d = 0.0_lkind
    end if

  end subroutine create_ref_cosz_allow_neg_0d

  subroutine create_ref_cosz_no_neg_2d(lat_2d, lon_2d, gmt, h, dec, cosz_2d, fracday_2d)

    implicit none
    integer, parameter :: n = 4
    real(kind=TEST_AST_KIND_), intent(in), dimension(n,n)  :: lat_2d, lon_2d, h
    real(kind=TEST_AST_KIND_), intent(in)                  :: gmt, dec
    real(kind=TEST_AST_KIND_), intent(out), dimension(n,n) :: cosz_2d, fracday_2d
    real(kind=TEST_AST_KIND_), dimension(n,n)              :: aa, bb, t
    integer, parameter                                     :: lkind = TEST_AST_KIND_

    !! these reference values will not allow negative cosz values,
    !! i.e. Lallow_nagative = .false.
    !! they will also not include time averaging, i.e. .not. present(dt)
    !! fracday values calculated here

    aa = sin(lat_2d)*sin(dec)
    bb = cos(lat_2d)*cos(dec)

    t = gmt + lon_2d - real(PI, TEST_AST_KIND_)
    where (t >= real(PI, TEST_AST_KIND_))  t = t - real(twopi, TEST_AST_KIND_)
    where (t < real(-PI, TEST_AST_KIND_))  t = t + real(twopi, TEST_AST_KIND_)

    where (abs(t) < h)
      cosz_2d    = aa + bb * cos(t)
      fracday_2d = 1.0_lkind
    else where
      cosz_2d    = 0.0_lkind
      fracday_2d = 0.0_lkind
    end where

  end subroutine create_ref_cosz_no_neg_2d

  subroutine create_ref_cosz_no_neg_1d(lat_1d, lon_1d, gmt, h, dec, cosz_1d, fracday_1d)

    implicit none
    integer, parameter :: n = 16
    real(kind=TEST_AST_KIND_), intent(in), dimension(n)  :: lat_1d, lon_1d, h
    real(kind=TEST_AST_KIND_), intent(in)                :: gmt, dec
    real(kind=TEST_AST_KIND_), intent(out), dimension(n) :: cosz_1d, fracday_1d
    real(kind=TEST_AST_KIND_), dimension(n)              :: aa, bb, t
    integer, parameter                                   :: lkind = TEST_AST_KIND_

    !! these reference values will not allow negative cosz values,
    !! i.e. Lallow_nagative = .false.
    !! they will also not include time averaging, i.e. .not. present(dt)
    !! fracday values calculated here

    aa = sin(lat_1d)*sin(dec)
    bb = cos(lat_1d)*cos(dec)

    t = gmt + lon_1d - real(PI, TEST_AST_KIND_)
    where (t >= real(PI, TEST_AST_KIND_))  t = t - real(twopi, TEST_AST_KIND_)
    where (t < real(-PI, TEST_AST_KIND_))  t = t + real(twopi, TEST_AST_KIND_)

    where (abs(t) < h)
      cosz_1d    = aa + bb * cos(t)
      fracday_1d = 1.0_lkind
    else where
      cosz_1d    = 0.0_lkind
      fracday_1d = 0.0_lkind
    end where

  end subroutine create_ref_cosz_no_neg_1d

  subroutine create_ref_cosz_no_neg_0d(lat_0d, lon_0d, gmt, h, dec, cosz_0d, fracday_0d)

    implicit none
    real(kind=TEST_AST_KIND_), intent(in)  :: lat_0d, lon_0d, gmt, h, dec
    real(kind=TEST_AST_KIND_), intent(out) :: cosz_0d, fracday_0d
    real(kind=TEST_AST_KIND_)              :: aa, bb, t
    integer, parameter                     :: lkind = TEST_AST_KIND_


    !! these reference values will not allow negative cosz values,
    !! i.e. Lallow_nagative = .false.
    !! they will also not include time averaging, i.e. .not. present(dt)
    !! fracday values calculated here

    aa = sin(lat_0d)*sin(dec)
    bb = cos(lat_0d)*cos(dec)

    t = gmt + lon_0d - real(PI, TEST_AST_KIND_)
    if (t >= real(PI, TEST_AST_KIND_))  t = t - real(twopi, TEST_AST_KIND_)
    if (t < real(-PI, TEST_AST_KIND_))  t = t + real(twopi, TEST_AST_KIND_)

    if (abs(t) < h) then
      cosz_0d    = aa + bb * cos(t)
      fracday_0d = 1.0_lkind
    else
      cosz_0d    = 0.0_lkind
      fracday_0d = 0.0_lkind
    end if

  end subroutine create_ref_cosz_no_neg_0d

end program test_diurnal_solar