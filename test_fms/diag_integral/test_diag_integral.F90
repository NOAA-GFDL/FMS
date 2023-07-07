program test_diag_integral

  use diag_integral_mod
  use fms_mod, only : fms_init, fms_end
  use fms2_io_mod, only: fms2_io_init, open_file, close_file, FmsNetcdfFile_t
  use fms2_io_mod, only: register_axis, register_field, write_data
  use mpp_mod, only: mpp_init, mpp_sync, mpp_npes, mpp_get_current_pelist, mpp_error, FATAL
  use time_manager_mod, only:time_type, set_time
  use constants_mod, only: PI
  use platform_mod, only : r4_kind, r8_kind

  implicit none

  integer, parameter :: nxy=20
  integer, parameter :: nxyp=nxy+1
  real(r8_kind) :: lat(nxyp,nxyp), lon(nxyp,nxyp), area(nxy,nxy)
  !> test arrays
  real(TEST_DI_KIND_) :: immadeup2(nxy,nxy), immadeup3(nxy,nxy,nxy)
  real(TEST_DI_KIND_) :: immadeuph(nxy,nxy), weight(nxy,nxy,nxy)

  type(FmsNetcdfFile_t):: fileobj        !< Fileobj for the files written by the test

  character(100), parameter :: ncfile='INPUT/sample.tile1.nc'
  character(9), parameter :: field_name2='immadeup2'
  character(9), parameter :: field_name3='immadeup3'
  character(9), parameter :: field_nameh='immadeuph'
  character(8), parameter :: std_digits   = 'e23.15e3' !'e13.6e2'

  type(time_type) :: Time_init, Time

  !testing and generating answers
  integer :: i, j, k
  real(r8_kind) :: area_sum, itime, field_avg2, field_avg3, field_wght_avg3, field_avgh


  integer, parameter :: lkind=TEST_DI_KIND_

  call fms_init
  call write_netcdf_file !< initializes arrays

  !> set time; made up time
  Time_init=set_time(0,1,0)
  Time=set_time(0,2,0)

  call test_diag_integral_init
  call test_diag_integral_field_init

  call test_call_diag_integral_field
  call test_call_sum_diag_integral_field
  call diag_integral_end(Time)
  call fms_end
  !> diag_integral.out file is printed out by diag_integral_end

  !> compute total area
  area_sum=0.0_r8_kind
  do j=1, nxy
     do i=1, nxy
        area_sum = area_sum + area(i,j)
     end do
  end do

  call read_diag_integral_file
  call test_sum_diag_integral_field_2d
  call test_sum_diag_integral_field_3d
  call test_sum_field_wght_3d

contains
  !-------------------------------------
  !-------------------------------------
  subroutine test_diag_integral_init

    !> this test subroutine ony checks that diag_integral_init
    !! can be called successfully

    implicit none

    call diag_integral_init(Time_init=Time_init, Time=Time, blon=lon, blat=lat, area_in=area)

  end subroutine test_diag_integral_init
  !-------------------------------------
  !-------------------------------------
  subroutine test_diag_integral_field_init

    !> this test subroutine ony checks that diag_integral_init
    !! can be called successfully

    implicit none

    call diag_integral_field_init(field_name2, std_digits)
    call diag_integral_field_init(field_name3, std_digits)
    call diag_integral_field_init(field_nameh, std_digits)

  end subroutine test_diag_integral_field_init
  !-------------------------------------
  !-------------------------------------
  subroutine test_call_diag_integral_field

    implicit none
    integer :: is, ie, js, je

    is=1 ; ie=nxy/2
    js=1 ; je=nxy !je is not used

    Time=set_time(0,3,0)
    call sum_diag_integral_field(field_name2, immadeup2)
    call sum_diag_integral_field(field_name3, immadeup3)
    do js=1, nxy
       call sum_diag_integral_field(field_nameh, immadeuph, is, ie, js, je)
    end do
    call diag_integral_output(Time)

  end subroutine test_call_diag_integral_field
  !-------------------------------------
  !-------------------------------------
  subroutine test_call_sum_diag_integral_field

    implicit none
    integer :: is, ie, js, je

    is=1 ; ie=nxy/2
    js=1 ; je=nxy

    Time=set_time(0,4,0)
    call sum_diag_integral_field(field_name2, immadeup2)
    call sum_diag_integral_field(field_name3, immadeup3, weight)
    do js=1, nxy
       call sum_diag_integral_field(field_nameh, immadeuph, is, ie, js, je)
    end do
    call diag_integral_output(Time)

  end subroutine test_call_sum_diag_integral_field
  !-------------------------------------
  !-------------------------------------
  subroutine test_sum_diag_integral_field_2d

    implicit none

    real(r8_kind) :: answer2

    !> compute answer for 2d
    answer2=0.0_r8_kind
    do j=1, nxy
       do i=1, nxy
          answer2 = answer2 + real(immadeup2(i,j),r8_kind)*area(i,j)
       end do
    end do

    answer2=answer2/area_sum

    call check_answers(answer2, field_avg2, 'sum_diag_integral_field failed for 2d')
    call check_answers(answer2/2.0,field_avgh, 'sum_diag_integral_field_hemi failed')

  end subroutine test_sum_diag_integral_field_2d
  !-------------------------------------
  !-------------------------------------
  subroutine test_sum_diag_integral_field_3d

    implicit none

    real(r8_kind) :: answer3, tmp

    !> compute answer for 3d
    answer3=0.0_r8_kind
    do j=1, nxy
       do i=1, nxy
          tmp=0.0_r8_kind
          do k=1, nxy
             tmp=tmp+real(immadeup3(i,j,k),r8_kind)
          end do
          answer3 = answer3 + tmp*area(i,j)
       end do
    end do

    answer3=answer3/area_sum
    call check_answers(answer3,field_avg3,'sum_diag_integral_field failed for 3d')

  end subroutine test_sum_diag_integral_field_3d
  !-------------------------------------
  !-------------------------------------
  subroutine test_sum_field_wght_3d

    implicit none

    real(r8_kind) :: tmp3(nxy,nxy,nxy)
    real(r8_kind) :: answer3, tmp, wght

    !> adding onto it with weights
    tmp3=0.0_r8_kind
    do k=1, nxy
       do j=1, nxy
          do i=1, nxy
             tmp3(i,j,k)=real(immadeup3(i,j,k),r8_kind)*real(weight(i,j,k),r8_kind)
          end do
       end do
    end do

    answer3=0.0_r8_kind
    do j=1, nxy
       do i=1, nxy
          tmp=0.0_r8_kind
          wght=0.0_r8_kind
          do k=1, nxy
             wght=wght+real(weight(i,j,k),r8_kind)
             tmp=tmp+tmp3(i,j,k)
          end do
          answer3 = answer3 + tmp*area(i,j)/wght
       end do
    end do

    answer3 = answer3/area_sum
    call check_answers(answer3,field_wght_avg3,'sum_diag_field_wght failed for 3d')

  end subroutine test_sum_field_wght_3d
  !-------------------------------------
  !-------------------------------------
  subroutine read_diag_integral_file

    character(17), parameter :: di_file='diag_integral.out'
    integer, parameter  :: iunit=100

    character(100) :: cline1, cline2, cline3, cline4, cline5

    !> read in computed values
    open(unit=iunit,file=trim(di_file))
    read(iunit,*) cline1, cline2, cline3, cline4, cline5
    read(iunit,*) itime, field_avg2, field_avg3, field_avgh
    read(iunit,*) cline1, cline2, cline3, cline4, cline5
    !> only the answers for the last time step will be checked
    read(iunit,*) cline1, field_avg2, field_wght_avg3, field_avgh
    close(iunit)

  end subroutine read_diag_integral_file
  !-------------------------------------
  !-------------------------------------
  subroutine check_answers(answer, outresult, whoami)

    implicit none
    real(r8_kind), intent(in) :: answer, outresult
    character(len=*), intent(in) :: whoami

    real(r8_kind), parameter :: tol=1.e-3_r8_kind

    if( abs(answer-outresult)>tol ) then
       write(*,*) '*******************************************'
       write(*,*) 'expected', answer, 'but computed',  outresult
       call mpp_error(FATAL,'ERROR: '//trim(whoami))
    end if

  end subroutine check_answers
  !-------------------------------------
  !-------------------------------------
  subroutine write_netcdf_file

    !> made up numbers

    implicit none

    integer :: i,j,k
    integer, allocatable :: pes(:)
    real(r8_kind), parameter :: dxy=0.1_r8_kind

    do i=1, nxyp
       lon(i,:)=real(i,r8_kind)*dxy
    end do

    do j=1, nxyp
       lat(:,j)=15.0_r8_kind*real(j,r8_kind)*dxy
    end do

    do j=1, nxy
       do i=1, nxy
          area(i,j)=0.12_r8_kind*real(i,r8_kind)+real(j,r8_kind)
       end do
    end do

    do j=1, nxy
       do i=1, nxy
          immadeup2(i,j)=100.0_lkind*real(j-1,lkind)+real(i,lkind)*real(PI,lkind)
       end do
    end do

    immadeuph=immadeup2

    do k=1, nxy
       do j=1, nxy
          do i=1, nxy
             immadeup3(i,j,k)=real(k-1,lkind) + 150_lkind*real(j-1,lkind) + real(i,lkind)*real(PI,lkind)
             weight(i,j,k)=real(i*k, lkind)/real(j,lkind)*real(PI,lkind)
          end do
       end do
    end do

    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)

    if( open_file(fileobj, ncfile, 'overwrite', pelist=pes) ) then
       call register_axis(fileobj, 'nx', nxy)
       call register_axis(fileobj, 'ny', nxy)
       call register_axis(fileobj, 'nz', nxy)
       call register_axis(fileobj, 'nxp', nxyp)
       call register_axis(fileobj, 'nyp', nxyp)

       call register_field(fileobj, 'x', 'double', dimensions=(/'nxp', 'nyp'/))
       call register_field(fileobj, 'y', 'double', dimensions=(/'nxp', 'nyp'/))
       call register_field(fileobj, 'area', 'double', dimensions=(/'nx','ny'/))
       call register_field(fileobj, 'immadeup2', 'double', dimensions=(/'nx','ny'/))
       call register_field(fileobj, 'immadeup3', 'double', dimensions=(/'nx','ny','nz'/))
       call register_field(fileobj, 'immadeuph', 'double', dimensions=(/'nx','ny'/))

       call write_data(fileobj, 'x', lon)
       call write_data(fileobj, 'y', lat)
       call write_data(fileobj, 'area', area)
       call write_data(fileobj, 'immadeup2', immadeup2)
       call write_data(fileobj, 'immadeup3', immadeup3)
       call write_data(fileobj, 'immadeuph', immadeuph)

    end if
    call close_file(fileobj)

  end subroutine write_netcdf_file

end program test_diag_integral
