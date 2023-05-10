program test_diag_integral

  use diag_integral_mod
  use fms_mod, only : fms_init
  use fms2_io_mod, only: fms2_io_init, open_file, close_file, FmsNetcdfFile_t
  use fms2_io_mod, only: register_axis, register_field, write_data
  use mpp_mod, only: mpp_init, mpp_sync, mpp_npes, mpp_get_current_pelist
  use time_manager_mod, only:time_type, set_time
  use constants_mod, only: PI

  implicit none

  integer, parameter :: nxy=20 !< supergrid
  integer, parameter :: nxyp=nxy+1
  real :: lat(nxyp,nxyp), lon(nxyp,nxyp), area(nxy,nxy)
  real :: immadeup2(nxy,nxy), immadeup3(nxy,nxy,nxy)

  type(FmsNetcdfFile_t):: fileobj        !< Fileobj for the files written by the test

  character(100), parameter :: ncfile='INPUT/sample.tile1.nc'
  character(9), parameter :: field_name2='immadeup2'
  character(9), parameter :: field_name3='immadeup3'
  character(8), parameter :: std_digits   = 'f20.10'

  type(time_type) :: Time_init, Time !< current time

  call fms_init
  call write_netcdf_file

  !> set time; made up time
  Time_init=set_time(0,1,0)
  Time=set_time(0,3,0)

  call test_diag_integral_init
  call test_diag_integral_field_init
  call test_sum_diag_integral_field

contains
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

  end subroutine test_diag_integral_field_init
  !-------------------------------------
  !-------------------------------------
  subroutine test_sum_diag_integral_field

    implicit none

    call sum_diag_integral_field(field_name2, immadeup2)
    call sum_diag_integral_field(field_name3, immadeup3)
    call diag_integral_output(Time)

  end subroutine test_sum_diag_integral_field
  !-------------------------------------
  !-------------------------------------
  subroutine write_netcdf_file

    !> made up numbers

    implicit none

    integer :: i,j,k
    integer, allocatable :: pes(:)
    real, parameter :: dxy=0.1

    do i=1, nxyp
       lon(i,:)=real(i)*dxy
    end do

    do j=1, nxyp
       lat(:,j)=real(j)*dxy*15.0
    end do

    do j=1, nxy
       do i=1, nxy
          area(i,j)=0.12*real(i)+real(j)
       end do
    end do

    do j=1, nxy
       do i=1, nxy
          immadeup2(i,j)=real(j-1)*100.0+real(i)*PI
       end do
    end do

    do k=1, nxy
       do j=1, nxy
          do i=1, nxy
             immadeup3(i,j,k)=real(k-1)*1000. + real(j-1)*100. + real(i)*PI
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

       call write_data(fileobj, 'x', lon)
       call write_data(fileobj, 'y', lat)
       call write_data(fileobj, 'area', area)
       call write_data(fileobj, 'immadeup2', immadeup2)
       call write_data(fileobj, 'immadeup3', immadeup3)

    end if
    call close_file(fileobj)

  end subroutine write_netcdf_file

end program test_diag_integral
