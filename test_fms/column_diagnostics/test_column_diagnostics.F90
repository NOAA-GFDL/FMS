program test_column_diagnostics

  use column_diagnostics_mod
  use fms_mod, only: fms_init
  use mpp_mod, only: FATAL, mpp_error
  use time_manager_mod, only: time_manager_init, time_type, set_time, set_calendar_type
  use constants_mod, only : PI, DEG_TO_RAD
  use platform_mod, only: r4_kind, r8_kind

  implicit none

  character(13), parameter :: mod_name='pemberley_mod'
  integer, parameter :: num_diag_pts_latlon=2 !< number of diagnostics column
  integer, parameter :: num_diag_pts_ij=2     !< number of diagnostics column
  integer :: global_i(num_diag_pts_ij) !<in; global i coordinates for ij diagnostic column
  integer :: global_j(num_diag_pts_ij) !<in; global j coordinates for ij diagnostic column
  real(TEST_CD_KIND_) :: global_lat_latlon(num_diag_pts_latlon)!<in; latitudes for lat-lon diagnostic columns
  real(TEST_CD_KIND_) :: global_lon_latlon(num_diag_pts_latlon)!<in; longitude for lat-lon diagnostic column

  integer, parameter :: nlatlon=6
  real(TEST_CD_KIND_) :: lonb_in(nlatlon,nlatlon) !< in
  real(TEST_CD_KIND_) :: latb_in(nlatlon,nlatlon) !< in
  logical :: do_column_diagnostics(nlatlon,nlatlon) !< out

  integer, parameter :: num_diag_pts=num_diag_pts_latlon + num_diag_pts_ij
  integer :: diag_i(num_diag_pts) !< out;
  integer :: diag_j(num_diag_pts) !< out
  real(TEST_CD_KIND_) :: diag_lat(num_diag_pts)   !< out
  real(TEST_CD_KIND_) :: diag_lon(num_diag_pts)   !< out
  integer :: diag_units(num_diag_pts)

  integer, parameter :: lkind=TEST_CD_KIND_

  call fms_init()
  call time_manager_init()
  call initialize_variables
  call column_diagnostics_init()
  call test_initialize_diagnostic_columns
  call test_column_diagnostics_header

contains
  !------------------------------------------!
  subroutine initialize_variables

    implicit none

    real(lkind) :: dlat, dlon
    integer :: i

    !> lat lon coordinates in degrees; made up.
    dlat=15.0_lkind
    dlon=15.0_lkind
    do i=1, nlatlon
       lonb_in(i,:)=real(i,lkind)*dlat - 0.5_lkind*dlat
       latb_in(:,i)=-90._lkind + real(i,lkind)*dlon -0.5_lkind*dlat
    end do

    !> initialize_diagnostic_columns expects these values to be in degrees
    global_lon_latlon(1)=lonb_in(2,1) ; global_lon_latlon(2)=lonb_in(3,1)
    global_lat_latlon(1)=latb_in(1,2) ; global_lat_latlon(2)=latb_in(1,3)
    global_i(1)=4 ; global_i(2)=5
    global_j(1)=4 ; global_j(2)=5

    !> intialize_diagnostic_columns expects these values to be in radians
    lonb_in=lonb_in*DEG_TO_RAD
    latb_in=latb_in*DEG_TO_RAD


  end subroutine initialize_variables
 !------------------------------------------!
  subroutine test_initialize_diagnostic_columns

    implicit none
    integer :: i

    integer :: i_answers(num_diag_pts), j_answers(num_diag_pts)
    real(TEST_CD_KIND_) :: lon_answers(num_diag_pts), lat_answers(num_diag_pts)

    call initialize_diagnostic_columns(mod_name, num_diag_pts_latlon, num_diag_pts_ij, &
                                       global_i, global_j, global_lat_latlon, global_lon_latlon, &
                                       lonb_in, latb_in, do_column_diagnostics, &
                                       diag_lon, diag_lat, diag_i, diag_j, diag_units)

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
