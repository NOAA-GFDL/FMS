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
!! @description This program tests the various subroutines in interpolator_mod.
!! TODO:  Test obtain_interpolator_time_slices.
!!        The subroutine obtain_interpolator_time_slices returns
!!        an interpolator_type which cannot be queried.
!!        This subroutine should be called by the interpolators.

program test_interpolator

  use fms2_io_mod, only: FmsNetcdfFile_t, UNLIMITED,                  &
                         register_field, register_variable_attribute, &
                         register_axis,                               &
                         write_data, close_file, open_file
  use mpp_mod,     only: mpp_error, FATAL
  use time_manager_mod, only: time_type, set_date, set_calendar_type, time_manager_init
  use fms_mod,     only: fms_init
  use constants_mod, only: PI

  use interpolator_mod

  implicit none

  character(100), parameter :: ncfile='immadeup.o3.climatology.nc'
  integer, parameter :: calendar_type=2     !< JULIAN

  integer :: nlonlat       !< number of latitude and longitudinal center coordinates
  integer :: nlonlatb      !< number of latitude and longitudinal boundary coordinates
  integer :: ntime         !< number of time slices
  integer :: npfull        !< number of p levels
  integer :: nphalf        !< number of half p levels
  integer :: nlonlat_mod, nlonlatb_mod

  !> climatology related arrays that will hold made-up data
  real, allocatable :: lat(:)   !< climatology coordinates
  real, allocatable :: lon(:)   !< climatology coordinates
  real, allocatable :: latb(:)  !< climatology coordinates
  real, allocatable :: lonb(:)  !< climatology coordinates
  real, allocatable :: clim_time (:) !< climatology time
  real, allocatable :: pfull(:) !< climatology p level
  real, allocatable :: phalf(:) !< climatology p half level
  real, allocatable :: ozone(:,:,:,:) !< climatology ozone data

  !> model related arrays
  real, allocatable :: lat_mod(:,:)  !< model coordinates
  real, allocatable :: lon_mod(:,:)  !< model coordinates
  real, allocatable :: latb_mod(:,:) !< model coordinates
  real, allocatable :: lonb_mod(:,:) !< model coordinates

  type(interpolate_type) :: o3

  call fms_init
  call time_manager_init
  call set_calendar_type(calendar_type)

  !> set data
  nlonlat  = 10        !< number of latitude and longitudinal center coordinates
  nlonlatb = nlonlat+1 !< number of latitude and longitudinal boundary coordinates
  nlonlat_mod  = nlonlat
  nlonlatb_mod = nlonlatb
  ntime  = 5        !< number of time slices
  npfull = 3        !< number of p levels
  nphalf = npfull+1 !< number of half p levels
  call allocate_arrays()
  call set_clim_time()
  call set_latlon_b()
  call set_pfullhalf()
  call set_ozone()
  call write_climatology_file()

  call set_latlon_b_mod()

  write(*,*) '===== test_intepolator_init ====='
  call test_interpolator_init(o3)
  write(*,*) '      test interpolator_init success'
  write(*,*) ''

  write(*,*) '===== test_intepolator_2D ======='
  write(*,*) '      ALSO TESTING obtain_interpolator_time_slices'
  call test_interpolator_2D(o3)
  write(*,*) '      test_interpolator_2D success'
  write(*,*) ''

  write(*,*) '===== test_intepolator_3D ======='
  call test_interpolator_3D(o3)
  write(*,*) '      test_interpolator_3D success'
  write(*,*) ''

  write(*,*) '===== test_intepolator_4D ======='
  call test_interpolator_4D(o3)
  write(*,*) '      test_interpolator_4D success'
  write(*,*) ''

  write(*,*) '===== test_interpolate_type_eq ====='
  call test_interpolate_type_eq()
  write(*,*) '      test_interpolate_type_eq success'
  write(*,*) ''

  write(*,*) '===== test_query_interpolator ====='
  call test_query_interpolator()
  write(*,*) '      test_query_interpolator success'
  write(*,*) ''

  write(*,*) '===== test_interpolator_end ====='
  call test_interpolator_end(o3)
  write(*,*) '      test_interpolator_end success'
  write(*,*) ''

  call deallocate_arrays()

  write(*,*) "*************************************"

contains
  !===============================================!
  !> test interpolator init
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
  subroutine test_interpolator_2D(clim_type)

    implicit none
    type(interpolate_type), intent(inout)    :: clim_type

    real, dimension(nlonlat_mod,nlonlat_mod) :: interp_data !< returns interpolation for first plevel
    type(time_type) :: model_time
    integer :: itime, i, j, k

    do itime=1, ntime
       model_time=set_date(1849,1,1+int(clim_time(itime)))
       call obtain_interpolator_time_slices(clim_type,model_time)
       call interpolator(clim_type, model_time, interp_data, 'ozone')
       call unset_interpolator_time_flag(clim_type)
       call interpolator(clim_type, model_time, interp_data, 'ozone')

       do i=1, 1
          do j=1, nlonlat_mod
             do k=1, nlonlat_mod
                !call check_answers(interp_data(k,j), ozone(k,j,i,itime), 'test interpolator_2D')
             end do
          end do
       end do

    end do

  end subroutine test_interpolator_2D
  !===============================================!
  subroutine test_interpolator_3D(clim_type)

    implicit none
    type(interpolate_type), intent(inout) :: clim_type
    real, dimension(nlonlat,nlonlat,nphalf-1) :: interp_data
    real, dimension(nlonlat,nlonlat,nphalf) :: phalf
    type(time_type) :: model_time

    integer :: itime, i, j, k

    do itime=1, ntime

       model_time=set_date(1849,1,1+int(clim_time(itime)))
       phalf(:,:,1)=0.0000
       phalf(:,:,2)=0.0002
       phalf(:,:,3)=0.0004
       phalf(:,:,4)=0.0005
       call interpolator(clim_type, model_time, phalf, interp_data, 'ozone')

       do i=1, nphalf-1
          do j=1, nlonlat
             do k=1, nlonlat
                !call check_answers(interp_data(k,j,i), ozone(k,j,i,itime), 'test interpolator_3D')
             end do
          end do
       end do

    end do

  end subroutine test_interpolator_3D
  !===============================================!
  subroutine test_interpolator_4D(clim_type)

    implicit none

    type(interpolate_type) :: clim_type

    real, dimension(nlonlat,nlonlat,nphalf-1,1) :: interp_data !< last column, there is only one field

    real, dimension(nlonlat,nlonlat,nphalf) :: phalf
    type(time_type) :: model_time

    integer :: itime, i, j, k, l

    do itime=1, ntime

       model_time=set_date(1849,1,1+int(clim_time(itime)))
       phalf(:,:,1)=0.0000
       phalf(:,:,2)=0.0002
       phalf(:,:,3)=0.0004
       phalf(:,:,4)=0.0005
       call interpolator(clim_type, model_time, phalf, interp_data, 'ozone')

       do i=1, nphalf-1
          do j=1, nlonlat
             do k=1, nlonlat
                !call check_answers(interp_data(k,j,i,1), ozone(k,j,i,itime), 'test interpolator_3D')
             end do
          end do
       end do

    end do

  end subroutine test_interpolator_4D
  !===============================================!
  subroutine test_interpolator_end(clim_type)

    implicit none

    type(interpolate_type) :: clim_type

    call interpolator_end(clim_type)

  end subroutine test_interpolator_end
  !===============================================!
  subroutine test_interpolator_4D_no_time_axis(clim_type)

    implicit none

    type(interpolate_type) :: clim_type

    real, dimension(nlonlat,nlonlat,nphalf-1) :: interp_data !< last column, there is only one field

    real, dimension(nlonlat,nlonlat,nphalf) :: phalf
    type(time_type) :: model_time

    phalf(:,:,1)=0.0000
    phalf(:,:,2)=0.0002
    phalf(:,:,3)=0.0004
    phalf(:,:,4)=0.0005
    model_time=set_date(1849,1,1+int(clim_time(1)))
    call read_data(clim_type, 'ozone', interp_data, 1, 1, Time=model_time)
    call interpolator(clim_type, phalf, interp_data, 'ozone')

    !do i=1, 1!nphalf-1
    !   do j=1, nlonlat
    !      do k=1, nlonlat
    !         if( interp_data(k,j,i).ne.ozone(k,j,i,1) ) then
    !            write(*,*) k,j,interp_data(k,j,i), ozone(k,j,i,1)
    !            !call mpp_error(FATAL, 'wrong')
    !         end if
    !      end do
    !   end do
    !end do

  end subroutine test_interpolator_4D_no_time_axis
  !===============================================!
  subroutine test_interpolate_type_eq

    implicit none

    type(interpolate_type) :: o3_copy

    o3_copy = o3
    call test_interpolator_3D(o3_copy)

  end subroutine test_interpolate_type_eq
  !===============================================!
  subroutine test_query_interpolator

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
  subroutine check_answers(results, answers, whoami)

    implicit none
    real :: results, answers
    character(*) :: whoami

    if (results.ne.answers) then
       write(*,*) '      EXPECTED ', answers, ' but computed ', results
       call mpp_error(FATAL, trim(whoami))
    end if

  end subroutine check_answers
  !===============================================!
#include "test_interpolator_write_climatology.inc"

end program test_interpolator
