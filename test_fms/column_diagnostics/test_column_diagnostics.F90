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
!! @brief unit test for column_diagnostics_mod
!! @author MiKyung Lee
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program mainly tests initialize_diagnostics_columns.
!! TODO:  The current test only tests with 1 processor.  A test that uses
!! domain decomposition is needed.
program test_column_diagnostics

  use column_diagnostics_mod
  use fms_mod, only: fms_init
  use mpp_mod, only: FATAL, mpp_error
  use time_manager_mod, only: time_manager_init, time_type, set_time, set_calendar_type
  use constants_mod, only : PI, DEG_TO_RAD
  use platform_mod, only: r4_kind, r8_kind

  implicit none

  character(13), parameter :: mod_name='pemberley_mod' !< made up module name; Mr. Darcy's estate
  integer, parameter :: num_diag_pts_latlon=2 !< number of diagnostics column described in terms of latlon coordinates
  integer, parameter :: num_diag_pts_ij=2     !< number of diagnostics column describes in terms of i/j indices
  integer :: global_i(num_diag_pts_ij) ! global i coordinates of the diagnostic column
  integer :: global_j(num_diag_pts_ij) ! global j coordinates of the diagnostic column
  real(TEST_FMS_KIND_) :: global_lat_latlon(num_diag_pts_latlon)!< latitude value for the diagnostic column
  real(TEST_FMS_KIND_) :: global_lon_latlon(num_diag_pts_latlon)!< longitude value for the diagnostic columns

  integer, parameter :: nlatlon=6 !< number of latlon grid points
  real(TEST_FMS_KIND_) :: lonb_in(nlatlon,nlatlon) !< model longitude grid point
  real(TEST_FMS_KIND_) :: latb_in(nlatlon,nlatlon) !< model latitude point
  logical :: do_column_diagnostics(nlatlon,nlatlon) !< out

  integer, parameter :: num_diag_pts=num_diag_pts_latlon + num_diag_pts_ij !< total number of diagnostics column
  integer :: diag_i(num_diag_pts) !< out
  integer :: diag_j(num_diag_pts) !< out
  real(TEST_FMS_KIND_) :: diag_lat(num_diag_pts) !< out
  real(TEST_FMS_KIND_) :: diag_lon(num_diag_pts) !< out
  integer :: diag_units(num_diag_pts)

  integer, parameter :: lkind=TEST_FMS_KIND_ !< local kind; either r4_kind or r8_kind

  call fms_init()
  call time_manager_init()
  call initialize_variables(0.0_lkind) !< set up input arrays
  call column_diagnostics_init() !< initialize diagnostics column
  call initialize_variables(0.01_lkind) !< set up input arrays;
  call test_initialize_diagnostic_columns !< initialize diagnostics column
  call test_column_diagnostics_header

contains
  !------------------------------------------!
  subroutine initialize_variables(dlatlon)

    !> This subroutine initializes all the input arrays for initialize_diagnostic_columns

    implicit none

    real(lkind), intent(in) :: dlatlon !< in degrees; displace lat/lon grid by dlatlon
    real(lkind) :: dlat, dlon
    integer :: i

    !> lat lon coordinates in degrees; made up to match the diagnostic column coordinates +/- dlatlon
    !! see initialize_diagnostic_columns.  A-Grid coordinates
    dlat=15.0_lkind !< randomly chosen value
    dlon=15.0_lkind !< randomly chosen value
    do i=1, nlatlon
       lonb_in(i,:)=real(i,lkind)*dlat - 0.5_lkind*dlat
       latb_in(:,i)=-90._lkind + real(i,lkind)*dlon -0.5_lkind*dlat
    end do

    !> initialize_diagnostic_columns coordinates expects these values to be in degrees
    global_lon_latlon(1)=lonb_in(2,1)
    global_lon_latlon(2)=lonb_in(3,1)
    global_lat_latlon(1)=latb_in(1,2)
    global_lat_latlon(2)=latb_in(1,3)
    global_i(1)=4 ; global_i(2)=5
    global_j(1)=4 ; global_j(2)=5

    !> initialize_diagnostic_columns expects these values to be in radians
    lonb_in=(lonb_in+dlatlon)*DEG_TO_RAD
    latb_in=(latb_in+dlatlon)*DEG_TO_RAD


  end subroutine initialize_variables
 !------------------------------------------!
  subroutine test_initialize_diagnostic_columns

    !> this subroutine tests initialize_diagnostics_columns

    implicit none
    integer :: i

    integer :: i_answers(num_diag_pts), j_answers(num_diag_pts)
    real(TEST_FMS_KIND_) :: lon_answers(num_diag_pts), lat_answers(num_diag_pts)

    call initialize_diagnostic_columns(mod_name, num_diag_pts_latlon, num_diag_pts_ij, &
                                       global_i, global_j, global_lat_latlon, global_lon_latlon, &
                                       lonb_in, latb_in, do_column_diagnostics, &
                                       diag_lon, diag_lat, diag_i, diag_j, diag_units)

    !> the edge points do not count
    i_answers=(/2,3,4,5/)
    j_answers=(/2,3,4,5/)
    lon_answers=lonb_in(2:5,1)/DEG_TO_RAD
    lat_answers=latb_in(1,2:5)/DEG_TO_RAD

    do i=1, num_diag_pts
       call check_answers(i_answers(i), diag_i(i), 'test_initialize_diagnostics_column diag_i')
       call check_answers(j_answers(i), diag_j(i), 'test_initialize_diagnostics_column diag_j')
       call check_answers(lon_answers(i), diag_lon(i), 'test_initialize_diagnostics_column diag_lon')
       call check_answers(lat_answers(i), diag_lat(i), 'test_initialize_diagnostics_column diag_lon')
    end do

  end subroutine test_initialize_diagnostic_columns
  !------------------------------------------!
  subroutine test_column_diagnostics_header

    !> This subroutine only tests that column_diagnostics_header works

    implicit none
    integer :: nn, diag_unit
    type(time_type) :: Time

    diag_unit=45  !< will produce fort.45 file
    call set_calendar_type(2)
    Time=set_time(12,14,1)
    do nn=1, num_diag_pts
       call column_diagnostics_header(mod_name, diag_unit, Time, nn, diag_lon, diag_lat, diag_i, diag_j)
    end do

  end subroutine test_column_diagnostics_header
  !------------------------------------------!
  subroutine check_answers(answer, myvalue, whoami)

    implicit none
    class(*) :: answer
    class(*) :: myvalue
    character(*) :: whoami

    select type(answer)
    type is ( integer )
       select type(myvalue)
       type is( integer )
          if( answer .ne. myvalue ) then
             write(*,*) '*************************************'
             write(*,*) 'EXPECTED ', answer, 'but got ', myvalue
             call mpp_error( FATAL,'failed '//trim(whoami) )
          end if
       end select
    type is( real(r4_kind) )
       select type( myvalue)
       type is(real(r4_kind) )
          if( answer .ne. myvalue ) then
             write(*,*) '*************************************'
             write(*,*) 'EXPECTED ', answer, 'but got ', myvalue
             write(*,*) 'difference of', abs(answer-myvalue)
             call mpp_error( FATAL,'failed '//trim(whoami) )
          end if
       end select
    type is( real(r8_kind) )
       select type( myvalue)
       type is(real(r4_kind) )
          if( answer .ne. myvalue ) then
             write(*,*) '*************************************'
             write(*,*) 'EXPECTED ', answer, 'but got ', myvalue
             write(*,*) 'difference of', abs(answer-myvalue)
             call mpp_error( FATAL,'failed '//trim(whoami) )
          end if
       end select
    end select

  end subroutine check_answers
  !------------------------------------------!
end program test_column_diagnostics
