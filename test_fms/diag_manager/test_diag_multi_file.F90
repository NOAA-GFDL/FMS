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
program test_multi_file
  use fms_mod, only: fms_init, fms_end
  use diag_manager_mod, only: diag_manager_init, diag_axis_init, register_diag_field, send_data, &
                              diag_send_complete, diag_manager_set_time_end, diag_manager_end
  use mpp_mod, only: mpp_error, mpp_pe, mpp_root_pe, FATAL, input_nml_file
  use time_manager_mod, only: time_type, set_calendar_type, JULIAN, set_time, set_date, operator(+)
  use fms2_io_mod

  implicit none

  type(time_type)                    :: Time             !< Time of the simulation
  type(time_type)                    :: Time_step        !< Time_step of the simulation
  integer                            :: nx               !< Number of x points
  integer                            :: ny               !< Number of y points
  integer                            :: id_x             !< Axis id for the x dimension
  integer                            :: id_y             !< Axis id for the y dimension
  integer                            :: id_var1          !< Field id for 1st variable
  integer                            :: id_var2          !< Field id for 2nd variable
  logical                            :: used             !< Dummy argument to send_data
  real,                  allocatable :: x(:)             !< X axis data
  real,                  allocatable :: y(:)             !< Y axis_data
  real,                  allocatable :: var1_data(:,:)   !< Data for variable 1
  integer                            :: i                !< For do loops
  integer                            :: ntimes           !< Number of times to run the simulation for
  integer                            :: io_status        !< Status after reading the namelist

  integer :: test_case = 1 !< Test case to do:
                           !! 1 Test new_file_freq
                           !! 2 Test file start_time and file_duration with new_file_freq
                           !! 3 Test file start_time and file_duration without new_file_freq
                           !! 4 Flexible output timings
  namelist / test_multi_file_nml / test_case

  call fms_init()
  call set_calendar_type(JULIAN)
  call diag_manager_init()

  read (input_nml_file, test_multi_file_nml, iostat=io_status)
  if (io_status > 0) call mpp_error(FATAL,'=>test_multi_file: Error reading input.nml')

  nx = 10
  ny = 15

  ntimes = 24

  allocate(x(nx), y(ny))
  allocate(var1_data(nx,ny))
  do i=1,nx
    x(i) = i
  enddo
  do i=1,ny
    y(i) = -91 + i
  enddo

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time (3600,0) !< 1 hour

  id_x  = diag_axis_init('x',  x,  'point_E', 'x', long_name='point_E')
  id_y  = diag_axis_init('y',  y,  'point_N', 'y', long_name='point_N')

  id_var1 = register_diag_field  ('atmos', 'ua', (/id_x, id_y/), Time)

  call diag_manager_set_time_end(set_date(2,1,2,0,0,0))
  do i = 1, ntimes
    Time = Time + Time_step
    var1_data = real(i)
    used = send_data(id_var1, var1_data, Time)
    call diag_send_complete(Time_step)
  enddo

  call diag_manager_end(Time)

  call check_answers()
  call fms_end()

  contains

  subroutine check_file(fileobj, expected_ntimes, start, filename, expected_ans)
    type(FmsNetcdfFile_t), intent(in) :: fileobj
    integer,               intent(in) :: start
    integer,               intent(in) :: expected_ntimes
    character(len=*),      intent(in) :: filename
    real, optional,     intent(in) :: expected_ans

    integer :: j
    real :: ans_var
    real :: vardata_out(nx, ny)
    integer :: ntimes

    call get_dimension_size(fileobj, "time", ntimes)
    if (ntimes .ne. expected_ntimes) call mpp_error(FATAL, "The time dimension for "//trim(filename)//&
      " is not the correct size!")
    do j = 1, ntimes
      ans_var = real((start + 2*(j-1) + 1) + (start + 2*(j-1) + 2) ) / real(2.)
      if (present(expected_ans) .and. j .eq. 1) ans_var = expected_ans
      call read_data(fileobj, "ua", vardata_out, unlim_dim_level = j)
      if (any(vardata_out .ne. ans_var)) &
        call mpp_error(FATAL, "The data in "//trim(filename)//" is not the expected result!")
    enddo
  end subroutine

  subroutine check_answers_case_1()
    type(FmsNetcdfFile_t) :: fileobj

    if (.not. open_file(fileobj, "test_multi_file_0002_01_01_12.nc", "read")) &
      call mpp_error(FATAL, "Unable to open the file: test_multi_file_0002_01_01_12.nc")
    call check_file(fileobj, 6, 0, "test_multi_file_0002_01_01_12.nc")
    call close_file(fileobj)

    if (.not. open_file(fileobj, "test_multi_file_0002_01_02_00.nc", "read")) &
      call mpp_error(FATAL, "Unable to open the file: test_multi_file_0002_01_02_00.nc")
    call check_file(fileobj, 6, 12, "test_multi_file_0002_01_02_00.nc")

    if (file_exists("test_multi_file_0002_01_02_12.nc")) &
      call mpp_error(FATAL, "The file test_multi_file_0002_01_02_12.nc should not exist!")

    call close_file(fileobj)
  end subroutine check_answers_case_1

  subroutine check_answers_case_2()
    type(FmsNetcdfFile_t) :: fileobj

    if (file_exists("test_multi_file_0002_01_01_04.nc")) &
      call mpp_error(FATAL, "The file test_multi_file_0002_01_01_04 should not exist!")

    if (file_exists("test_multi_file_0002_01_01_08.nc")) &
      call mpp_error(FATAL, "The file test_multi_file_0002_01_01_08 should not exist!")

    if (file_exists("test_multi_file_0002_01_01_12.nc")) &
      call mpp_error(FATAL, "The file test_multi_file_0002_01_01_12 should not exist!")

    if (file_exists("test_multi_file_0002_01_02_00.nc")) &
      call mpp_error(FATAL, "The file test_multi_file_0002_01_02_00 should not exist!")

    if (.not. open_file(fileobj, "test_multi_file_0002_01_01_16.nc", "read")) &
      call mpp_error(FATAL, "Unable to open the file: test_multi_file_0002_01_01_16.nc")
    call check_file(fileobj, 2, 12, "test_multi_file_0002_01_01_16.nc", expected_ans = real(13))
    call close_file(fileobj)

    if (.not. open_file(fileobj, "test_multi_file_0002_01_01_20.nc", "read")) &
      call mpp_error(FATAL, "Unable to open the file: test_multi_file_0002_01_01_20.nc")
    call check_file(fileobj, 2, 16, "test_multi_file_0002_01_01_20.nc")
    call close_file(fileobj)

  end subroutine check_answers_case_2

  subroutine check_answers_case_3()
    type(FmsNetcdfFile_t) :: fileobj

    if (.not. open_file(fileobj, "test_multi_file.nc", "read")) &
      call mpp_error(FATAL, "Unable to open the file: test_multi_file.nc")
    call check_file(fileobj, 4, 12, "test_multi_file.nc", expected_ans = real(13))
    call close_file(fileobj)
  end subroutine check_answers_case_3

  subroutine check_answers_case_4()
    type(FmsNetcdfFile_t) :: fileobj

    if (.not. open_file(fileobj, "test_multi_file_0002_01_01_04.nc", "read")) &
      call mpp_error(FATAL, "Unable to open the file: test_multi_file_0002_01_01_04.nc")
    call check_file(fileobj, 2, 0, "test_multi_file_0002_01_01_04.nc")
    call close_file(fileobj)

    if (.not. open_file(fileobj, "test_multi_file_0002_01_01_08.nc", "read")) &
      call mpp_error(FATAL, "Unable to open the file: test_multi_file_0002_01_01_08.nc")
    call check_file(fileobj, 2, 4, "test_multi_file_0002_01_01_08.nc")
    call close_file(fileobj)

    if (.not. open_file(fileobj, "test_multi_file_0002_01_01_10.nc", "read")) &
      call mpp_error(FATAL, "Unable to open the file: test_multi_file_0002_01_01_10.nc")
    call check_file(fileobj, 1, 8, "test_multi_file_0002_01_01_10.nc")
    call close_file(fileobj)

    if (.not. open_file(fileobj, "test_multi_file_0002_01_01_22.nc", "read")) &
      call mpp_error(FATAL, "Unable to open the file: test_multi_file_0002_01_01_22.nc")
    call check_file(fileobj, 6, 10, "test_multi_file_0002_01_01_22.nc")
    call close_file(fileobj)

  end subroutine check_answers_case_4

  subroutine check_answers()
    select case (test_case)
    case (1)
      call check_answers_case_1()
    case (2)
      call check_answers_case_2()
    case (3)
      call check_answers_case_3()
    case (4)
      call check_answers_case_4()
    end select
  end subroutine check_answers
end program test_multi_file