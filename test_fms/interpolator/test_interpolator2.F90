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
!! @brief unit tests for interpolator_mod
!! @author MiKyung Lee
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program tests the various subroutines in interpolator_mod.  The code in
!! in test_interpolator.F90 shows how to use interpolator_mod.  This program contains the
!! actual testing.
!! TODO:

program test_interpolator2

  use fms2_io_mod, only: FmsNetcdfFile_t, UNLIMITED,                  &
                         register_field, register_variable_attribute, &
                         register_axis,                               &
                         write_data, close_file, open_file
  use mpp_mod,          only: mpp_error, FATAL
  use time_manager_mod, only: time_type, set_date, set_calendar_type, time_manager_init
  use fms_mod,          only: fms_init
  use constants_mod,    only: PI
  use platform_mod,     only: r4_kind, r8_kind

  use interpolator_mod

  implicit none

  character(100), parameter :: ncfile='immadeup.o3.climatology.nc' !< fake climatology file.
  integer, parameter :: calendar_type=2  !< JULIAN calendar
  integer, parameter :: lkind=TEST_INTP_KIND_
  real(r8_kind), parameter :: tol=1.e-5_r8_kind !< the interpolation methods are not perfect.
                                                !! Will not get perfectly agreeing answers

  !> climatology related variables and arrays (made up data)
  integer :: nlonlat       !< number of latitude and longitudinal center coordinates
  integer :: nlonlatb      !< number of latitude and longitudinal boundary coordinates
  integer :: ntime         !< number of time slices
  integer :: npfull        !< number of p levels
  integer :: nphalf        !< number of half p levels
  real(TEST_INTP_KIND_), allocatable :: lat(:)   !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: lon(:)   !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: latb(:)  !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: lonb(:)  !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: clim_time (:) !< climatology time
  real(TEST_INTP_KIND_), allocatable :: pfull(:) !< climatology p level
  real(TEST_INTP_KIND_), allocatable :: phalf(:) !< climatology p half level
  real(TEST_INTP_KIND_), allocatable :: ozone(:,:,:,:) !< climatology ozone data

  !> model related variables and arrays
  integer :: nlonlat_mod, nlonlatb_mod !< number of latitude and longitude coordinates in the model
  real(TEST_INTP_KIND_), allocatable :: lat_mod(:,:)  !< model coordinates
  real(TEST_INTP_KIND_), allocatable :: lon_mod(:,:)  !< model coordinates
  real(TEST_INTP_KIND_), allocatable :: latb_mod(:,:) !< model coordinates
  real(TEST_INTP_KIND_), allocatable :: lonb_mod(:,:) !< model coordinates

  type(interpolate_type) :: o3 !< recyclable interpolate_type

  call fms_init
  call time_manager_init
  call set_calendar_type(calendar_type)

  !> set data
  call set_write_data(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=5, npfull_in=3)

  !> test interpolator_init
  write(*,*) '===== test_interpolator_init ====='
  call test_interpolator_init(o3)

  !> test interpolator 2D-4D
  write(*,*) '===== test_intepolator ======='
  call test_interpolator(o3)

  !> test interpolate_type_eq
  !! This test has been commented out and will be included
  !! in the testing suite once fileobj cp is added into
  !! test_interpolate_type_eq
  !write(*,*) '===== test_interpolate_type_eq ====='
  !call test_interpolate_type_eq()

  !> test query_interpolator
  write(*,*) '===== test_query_interpolator ====='
  call test_query_interpolator()

  !> test interpolator end
  write(*,*) '===== test_interpolator_end ====='
  call test_interpolator_end(o3)

  !> deallocate all arrays used to write the .nc file and used for model coordinates
  call deallocate_arrays()

  !> test interpolator_no_time_axis
  !! Write out new set of data that will have a time axis, but will have "0" time points
  !! because that's how interpolator_init is set up.
  call set_write_data(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=0, npfull_in=3)
  call test_interpolator_init(o3)
  write(*,*) '===== test_intepolator_no_time_axis ======='
  call test_interpolator_no_time_axis(o3)

contains
  !===============================================!
  subroutine test_interpolator_init(clim_type)

    !> interopolator_init initializes the interpolate_type after reading in the
    !! climatology data.  The required arguments are clim_type, file_name, lonb_mod, latb_mod
    !! where lonb_mod and latb_mod contain the model longitude and latitude values on the grid.

    implicit none
    type(interpolate_type), intent(inout) :: clim_type
    integer, dimension(1) :: data_out_of_bounds

    data_out_of_bounds(1)=CONSTANT
    call interpolator_init(clim_type,trim(ncfile), lonb_mod, latb_mod, data_out_of_bounds=data_out_of_bounds)

  end subroutine test_interpolator_init
  !===============================================!
  subroutine test_interpolator(clim_type)

    !> call the variants of interpolator (4D-2d) that interpolates data at a given time-point
    !! The tests here do not test the "no_axis" interpolator routines
    !! This subroutine also tests obtain_interpolator_time_slices for the 2D case.

    implicit none

    type(interpolate_type), intent(inout) :: clim_type
    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,npfull,1) :: interp_data !< last column, there is only one field
    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,nphalf) :: phalf
    type(time_type) :: model_time
    integer :: itime, i, j, k, l

    phalf(:,:,1)=0.0000_lkind
    phalf(:,:,2)=0.0002_lkind
    phalf(:,:,3)=0.0004_lkind
    phalf(:,:,4)=0.0005_lkind

    do itime=1, ntime

       model_time=set_date(1849,1,1+int(clim_time(itime)))

       !> test interpolator_4D_r4/8
       call interpolator(clim_type, model_time, phalf, interp_data, 'ozone')
       do i=1, npfull
          do j=1, nlonlat
             do k=1, nlonlat
                call check_answers(interp_data(k,j,i,1), ozone(k,j,i,itime), tol, 'test interpolator_4D')
             end do
          end do
       end do

       !> test interpolator_3_r4/8
       call interpolator(clim_type, model_time, phalf, interp_data(:,:,:,1), 'ozone')
       do i=1, npfull
          do j=1, nlonlat
             do k=1, nlonlat
                call check_answers(interp_data(k,j,i,1), ozone(k,j,i,itime), tol, 'test interpolator_3D')
             end do
          end do
       end do

       !> test interpolator_2D_r4/8
       call interpolator(clim_type, model_time, interp_data(:,:,1,1), 'ozone')
       do j=1, nlonlat_mod
          do k=1, nlonlat_mod
             call check_answers(interp_data(k,j,1,1), ozone(k,j,1,itime), tol, 'test interpolator_2D')
          end do
       end do

       !> Test obtain_interpolator_time_slices
       call obtain_interpolator_time_slices(clim_type,model_time)
       call interpolator(clim_type, model_time, interp_data(:,:,1,1), 'ozone')
       call unset_interpolator_time_flag(clim_type)
       do j=1, nlonlat_mod
          do k=1, nlonlat_mod
             call check_answers(interp_data(k,j,1,1), ozone(k,j,1,itime), tol, 'test interpolator_2D')
          end do
       end do

    end do

  end subroutine test_interpolator
  !===============================================!
  subroutine test_interpolator_end(clim_type)

    !> This subroutine tests interpolator_end

    implicit none

    type(interpolate_type) :: clim_type

    call interpolator_end(clim_type)

  end subroutine test_interpolator_end
  !===============================================!
  subroutine test_interpolator_no_time_axis(clim_type)

    !> This subroutine tests the variants (42-2D) of interpolator_no_time_axis

    implicit none

    type(interpolate_type) :: clim_type

    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,nphalf-1,1) :: interp_data !< last column, there is only one field
    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,nphalf) :: phalf
    integer :: i, j, k

    phalf(:,:,1)=0.0000_lkind
    phalf(:,:,2)=0.0002_lkind
    phalf(:,:,3)=0.0004_lkind
    phalf(:,:,4)=0.0005_lkind

    !> test interpolator_4D_no_time_axis_r4/8
    call interpolator(clim_type, phalf, interp_data, 'ozone')
    do i=1, nphalf-1
       do j=1, nlonlat
          do k=1, nlonlat
             call check_answers(interp_data(k,j,i,1), ozone(k,j,i,1), tol, 'test interpolator_4D_no_time_axis')
          end do
       end do
    end do

    !> test interpolator_3D_no_time_axis_r4/8
    call interpolator(clim_type, phalf, interp_data(:,:,:,1), 'ozone')
    do i=1, nphalf-1
       do j=1, nlonlat
          do k=1, nlonlat
             call check_answers(interp_data(k,j,i,1), ozone(k,j,i,1), tol, 'test interpolator_3D_no_time_axis')
          end do
       end do
    end do

    !> test interpolator_2D_no_time_axis_r4/8
    call interpolator(clim_type, interp_data(:,:,1,1), 'ozone')
    do j=1, nlonlat
       do k=1, nlonlat
          call check_answers(interp_data(k,j,1,1), ozone(k,j,1,1), tol, 'test interpolator_2D_no_time_axis')
       end do
    end do

  end subroutine test_interpolator_no_time_axis
  !===============================================!
  subroutine test_interpolate_type_eq

    !> This subroutine tests interpolaote_type_eq (assignment = operator)
    !! The success of "=" is insured by checking to see if interpolation with o3_copy succeds.

    implicit none

    type(interpolate_type) :: o3_copy

    o3_copy = o3
    call test_interpolator(o3_copy)

  end subroutine test_interpolate_type_eq
  !===============================================!
  subroutine test_query_interpolator

    !> This subroutne tests query_interpolator

    implicit none

    character(100), parameter :: answer_field_name='ozone'
    integer, parameter :: answer_nfields=1

    character(100), allocatable :: field_names(:)
    integer :: nfields

    call query_interpolator(o3,nfields)
    if( nfields .ne. answer_nfields) call mpp_error(FATAL, '')

    allocate(field_names(nfields))
    call query_interpolator(o3,nfields,field_names)
    if(trim(answer_field_name).ne.trim(field_names(1))) call mpp_error(FATAL,'')

    deallocate(field_names)

  end subroutine test_query_interpolator
  !===============================================!
  subroutine check_answers(results, answers, tol, whoami)

    implicit none
    real(TEST_INTP_KIND_), intent(in) :: results, answers
    real(r8_kind), intent(in) :: tol
    character(*) :: whoami

    if (real(abs(results-answers),r8_kind).gt.tol) then
       write(*,*) '      EXPECTED ', answers, ' but computed ', results
       call mpp_error(FATAL, trim(whoami))
    end if

  end subroutine check_answers
  !===============================================!
#include "test_interpolator_write_climatology.inc"

end program test_interpolator2
