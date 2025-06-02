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
  use mpp_mod,          only: mpp_error, FATAL, WARNING
  use time_manager_mod, only: time_type, set_calendar_type, time_manager_init, &
                              get_date_no_leap, get_date_julian, set_date, set_time, &
                              set_date_no_leap, set_date_julian, operator(/), &
                              operator(+), operator(-), time_type_to_real, increment_time, &
                              leap_year, days_in_month, print_date, print_time
  use fms_mod,          only: fms_init
  use constants_mod,    only: PI
  use platform_mod,     only: r4_kind, r8_kind

  use interpolator_mod

  implicit none

  character(100), parameter :: ncfile='immadeup.o3.climatology.nc' !< fake climatology file
  integer, parameter :: lkind=TEST_INTP_KIND_
  !> the interpolation methods are not perfect.Will not get perfectly agreeing answers
  real(r8_kind), parameter :: tol=0.1_lkind
  integer :: calendar_type

  !> climatology related variables and arrays (made up data)
  integer :: nlonlat       !< number of latitude and longitudinal center coordinates in file
  integer :: nlonlatb      !< number of latitude and longitudinal boundary coordinates in file
  integer :: ntime         !< number of time slices
  integer :: npfull        !< number of p levels
  integer :: nphalf        !< number of half p levels
  real(TEST_INTP_KIND_), allocatable :: lat(:)   !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: lon(:)   !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: latb(:)  !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: lonb(:)  !< climatology coordinates
  real(r8_kind), allocatable :: clim_time (:) !< climatology time
  real(TEST_INTP_KIND_), allocatable :: pfull(:) !< climatology p level
  real(TEST_INTP_KIND_), allocatable :: phalf(:) !< climatology p half level
  real(TEST_INTP_KIND_), allocatable :: ozone(:,:,:,:) !< climatology ozone data

  !> model related variables and arrays
  integer :: nlonlat_mod, nlonlatb_mod !< number of latitude and longitude coordinates in the model
  real(TEST_INTP_KIND_), allocatable :: lat_mod(:,:)  !< model coordinates
  real(TEST_INTP_KIND_), allocatable :: lon_mod(:,:)  !< model coordinates
  real(TEST_INTP_KIND_), allocatable :: latb_mod(:,:) !< model coordinates
  real(TEST_INTP_KIND_), allocatable :: lonb_mod(:,:) !< model coordinates

  !> array holding model times
  type(time_type), allocatable :: model_time_julian(:), model_time_noleap(:)

  type(interpolate_type) :: o3 !< recyclable interpolate_type

  !> whether the file input is yearly, daily data
  logical :: yearly, daily
  logical :: noleap

  logical :: test_file_daily_julian=.true.,   test_file_daily_noleap=.false.
  logical :: test_file_yearly_noleap=.false., test_file_yearly_julian=.false.
  logical :: test_file_no_time=.false.
  integer :: nml_unit_var=99
  character(*), parameter :: nml_file='test_interpolator.nml'
  NAMELIST / test_interpolator_nml / test_file_daily_noleap, test_file_daily_julian, &
                                     test_file_yearly_noleap, test_file_yearly_julian, test_file_no_time

  open(unit=nml_unit_var, file=nml_file)
  read(unit=nml_unit_var, nml=test_interpolator_nml)
  close(nml_unit_var)

  call fms_init
  call time_manager_init
  call write_header

  !> set data
  call set_parameters_wrapper
  call set_and_write_data

  if(.not.test_file_no_time) then
     !> test interpolator when  model calendar is JULIAN
     calendar_type=2
     call set_calendar_type(calendar_type)
     call run_test_set
     !---------------------------------------------------------
     !> test interpolator when model calendar is NOLEAP
     calendar_type=4
     call set_calendar_type(calendar_type)
     call run_test_set
  end if

  !> test interpolator_no_time_axis
  if(test_file_no_time) then
     write(*,*) 'test_intepolator_no_time_axis'
     calendar_type=2  !< still need to set model calendar
     call set_calendar_type(calendar_type)
     call test_interpolator_init(o3)
     call test_interpolator_no_time_axis(o3)
     call test_interpolator_end(o3)
  end if

contains

#include "test_interpolator_write_climatology.inc"

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
  subroutine test_interpolator(clim_type, model_time)

    !> call the variants of interpolator (4D-2d) that interpolates data at a given time-point
    !! The tests here do not test the "no_axis" interpolator routines
    !! This subroutine also tests obtain_interpolator_time_slices for the 2D case.

    implicit none

    type(interpolate_type), intent(inout) :: clim_type
    type(time_type), dimension(ntime), intent(in) :: model_time
    type(time_type) :: test_time
    real(TEST_INTP_KIND_), dimension(nlonlat_mod,nlonlat_mod,npfull,1) :: interp_data !<only 1 field
    real(TEST_INTP_KIND_), dimension(nlonlat_mod,nlonlat_mod,nphalf) :: phalf_in
    integer :: itime, i, j, k, l

    real(TEST_INTP_KIND_) :: answer


    do i=1, nphalf
       phalf_in(:,:,i)=phalf(i)
    end do

    do itime=2, ntime-1

       answer=0.5_lkind*ozone(1,1,1,itime-1)+0.5_lkind*ozone(1,1,1,itime)
       test_time=model_time(itime-1) + (model_time(itime)-model_time(itime-1))/2

       !> test interpolator_4D_r4/8
       call interpolator(clim_type, test_time, phalf_in, interp_data, 'ozone')
       do i=1, npfull
          do j=1, nlonlat_mod
             do k=1, nlonlat_mod
                call check_answers(interp_data(k,j,i,1), answer, tol, 'test interpolator_4D')
             end do
          end do
       end do

       !> test interpolator_3_r4/8
       call interpolator(clim_type, test_time, phalf_in, interp_data(:,:,:,1), 'ozone')
       do i=1, npfull
          do j=1, nlonlat_mod
             do k=1, nlonlat_mod
                call check_answers(interp_data(k,j,i,1), answer, tol, 'test interpolator_3D')
             end do
          end do
       end do

       !> test interpolator_2D_r4/8
       call interpolator(clim_type, test_time, interp_data(:,:,1,1), 'ozone')
       do j=1, nlonlat_mod
          do k=1, nlonlat_mod
             call check_answers(interp_data(k,j,1,1), answer, tol, 'test interpolator_2D')
          end do
       end do

       !> Test obtain_interpolator_time_slices
       call obtain_interpolator_time_slices(clim_type,test_time)
       call interpolator(clim_type, test_time, interp_data(:,:,1,1), 'ozone')
       call unset_interpolator_time_flag(clim_type)
       do j=1, nlonlat_mod
          do k=1, nlonlat_mod
             call check_answers(interp_data(k,j,1,1), answer, tol, 'test interpolator_2D')
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

    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,npfull,1) :: interp_data !< last column, there is only one field
    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,nphalf) :: phalf_in
    integer :: i, j, k

    do i=1, nphalf
       phalf_in(:,:,i)=phalf(i)
    end do

    !> test interpolator_4D_no_time_axis_r4/8
    call interpolator(clim_type, phalf_in, interp_data, 'ozone')
    do i=1, npfull
       do j=1, nlonlat
          do k=1, nlonlat
             call check_answers(interp_data(k,j,i,1), ozone(k,j,i,1), tol, 'test interpolator_4D_no_time_axis')
          end do
       end do
    end do

    !> test interpolator_3D_no_time_axis_r4/8
    call interpolator(clim_type, phalf_in, interp_data(:,:,:,1), 'ozone')
    do i=1, npfull
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
    !! The success of "=" is insured by checking to see if interpolation with o3_copy succeeds.

    implicit none

    type(interpolate_type) :: o3_copy

    o3_copy = o3
    if(calendar_type==2) call test_interpolator(o3_copy, model_time_julian)
    if(calendar_type==4) call test_interpolator(o3_copy, model_time_noleap)

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
  subroutine run_test_set

    implicit none
    integer :: i

    if(calendar_type==2) write(*,*) "** MODEL CALENDAR JULIAN ** MODEL CALENDAR JULIAN ** MODEL CALENDAR JULIAN **"
    if(calendar_type==4) write(*,*) "** MODEL CALENDAR NOLEAP ** MODEL CALENDAR NOLEAP ** MODEL CALENDAR NOLEAP **"

    write(*,*) '1.  interpolator_init'
    call test_interpolator_init(o3)

    !> test interpolator 2D-4D
    write(*,*) '2.  interpolator4D-2D'
    if(calendar_type==2) call test_interpolator(o3, model_time_julian)
    if(calendar_type==4) call test_interpolator(o3, model_time_noleap)

    !> test interpolate_type_eq
    !! This test has been commented out and will be included in the testing suite once fileobj cp is added into
    write(*,*) '3.  skipping interpolate_type_eq until bug in interpolator is fixed'
    !call test_interpolate_type_eq()

    !> test query_interpolator
    write(*,*) '4.  query_interpolator'
    call test_query_interpolator()

    !> test interpolator end
    write(*,*) '5.  interpolator_end'
    call test_interpolator_end(o3)

  end subroutine run_test_set
  !===============================================!
  subroutine check_answers(results, answers, tol, whoami)

    implicit none
    real(TEST_INTP_KIND_), intent(in) :: results, answers
    real(r8_kind), intent(in) :: tol
    character(*) :: whoami

    if (real(abs(results-answers),r8_kind).gt.tol) then
    !if (results.ne.answers) then
       write(*,*) '      EXPECTED ', answers, ' but computed ', results
       call mpp_error(FATAL, trim(whoami))
    end if

  end subroutine check_answers
  !===============================================!
  subroutine write_header


    if(test_file_daily_noleap) &
         write(*,"(////10x,a,i0////)") &
         " ** DAILY FILE CAL NOLEAP ** DAILY FILE CAL NOLEAP ** DAILY FILE CAL NOLEAP ** ", lkind
    if(test_file_daily_julian)  &
         write(*,"(////10x,a,i0/////)") &
         ' ** DAILY FILE CAL JULIAN ** DAILY FILE CAL JULIAN ** DAILY FILE CAL JULIAN ** ', lkind
    if(test_file_yearly_noleap) &
         write(*,"(////10x,a,i0/////)") &
         ' ** YEARLY FILE CAL NOLEAP ** YEARLY FILE CAL NOLEAP ** YEARLY FILE CAL NOLEAP ** ', lkind
    if(test_file_yearly_julian) &
         write(*,"(////10x,a,i0/////)") &
         ' ** YEARLY FILE CAL JULIAN ** YEARLY FILE CAL JULIAN ** YEARLY FILE CAL JULIAN ** ', lkind
    if(test_file_no_time) &
         write(*,"(////10x,a/////)") " ** NO TIME AXIS ** NO TIME AXIS ** NO TIME AXIS **"

  end subroutine write_header
  !===============================================!
  subroutine set_parameters_wrapper

    if(test_file_daily_noleap)  then
       call set_parameters(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=240, npfull_in=3, &
            daily_in=.true., yearly_in=.false., noleap_in=.true.)
    else if(test_file_daily_julian) then
       call set_parameters(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=240, npfull_in=3, &
            daily_in=.true., yearly_in=.false., noleap_in=.false.)
    else if(test_file_yearly_noleap) then
       call set_parameters(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=240, npfull_in=3, &
            daily_in=.false., yearly_in=.true., noleap_in=.true.)
    else if(test_file_yearly_julian) then
       call set_parameters(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=240, npfull_in=3, &
            daily_in=.false., yearly_in=.true., noleap_in=.false.)
    else if(test_file_no_time) then
       call set_parameters(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=0, npfull_in=3, &
            daily_in=.true., yearly_in=.false., noleap_in=.false.)
    end if

  end subroutine set_parameters_wrapper
  !===============================================!
end program test_interpolator2
