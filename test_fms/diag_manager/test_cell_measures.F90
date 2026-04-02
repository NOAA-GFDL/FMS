!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @brief This program tests the diag_manager with fields with cell measures (area, volume)
program test_cell_measures
  use fms_mod, only: fms_init, fms_end
  use diag_manager_mod, only: diag_axis_init, register_static_field, diag_send_complete, send_data
  use diag_manager_mod, only : register_diag_field
  use diag_manager_mod, only: diag_manager_init, diag_manager_end, diag_manager_set_time_end
  use time_manager_mod, only: time_type, set_calendar_type, set_date, JULIAN, set_time, OPERATOR(+)
  use mpp_mod, only: mpp_error, FATAL
  use fms2_io_mod
  use platform_mod, only: r4_kind

  implicit none

  type(time_type)                 :: Time             !< Time of the simulation
  type(time_type)                 :: Time_step        !< Time_step of the simulation
  integer                         :: i                !< For looping
  integer                         :: id_axis1         !< Id of axis1
  integer                         :: naxis1           !< Size of axis1
  real(kind=r4_kind), allocatable :: axis1_data(:)    !< Data for axis1
  integer                         :: id_var1          !< Id of var1
  integer                         :: id_var2b         !< Id of var2b
  integer                         :: id_var2          !< Id of var2
  integer                         :: id_var3          !< Id of var3
  real(kind=r4_kind), allocatable :: var1_data(:)     !< Data for "var1"
  real(kind=r4_kind), allocatable :: area_data(:)     !< Data for the "area"
  integer                         :: id_area1         !< Id of the "area" field for var 1
  integer                         :: id_area2         !< Id of the "area" field for var 2
  integer                         :: id_area3         !< Id of the "area" field for var 3
  logical                         :: used             !< Used for send_data call

  naxis1 = 10
  call fms_init()
  call set_calendar_type(JULIAN)
  call diag_manager_init()

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time (3600,0)
  call diag_manager_set_time_end(set_date(2,1,2,0,0,0))

  allocate(axis1_data(naxis1))
  allocate(var1_data(naxis1))
  allocate(area_data(naxis1))
  do i = 1, naxis1
    axis1_data = real(i, kind=r4_kind)
    area_data = real(i/100, kind=r4_kind)
    var1_data = real(i*10, kind=r4_kind)
  enddo

  id_axis1  = diag_axis_init('axis1',  axis1_data,  'axis1', 'x')
  id_area1 = register_static_field ('fun_mod', 'var1_area', (/id_axis1/))
  id_area2 = register_static_field ('fun_mod', 'var2_area', (/id_axis1/))
  id_area3 = register_static_field ('fun_mod', 'var3_area', (/id_axis1/))
  id_var1 = register_diag_field  ('fun_mod', 'var1', (/id_axis1/), init_time=Time, area=id_area1)
  id_var2 = register_diag_field  ('fun_mod', 'var2', (/id_axis1/), init_time=Time, area=id_area2)
  id_var2b = register_diag_field  ('fun_mod', 'var2b', (/id_axis1/), init_time=Time, area=id_area1)
  id_var3 = register_diag_field  ('fun_mod', 'var3', (/id_axis1/), init_time=Time, area=id_area3)

  used = send_data(id_area1, area_data)
  used = send_data(id_area2, area_data)
  used = send_data(id_area3, area_data)

  do i = 1, 6
    Time = Time + Time_step
    call diag_send_complete(Time_step)
    used = send_data(id_var1, var1_data, Time)
    used = send_data(id_var2, var1_data, Time)
    used = send_data(id_var2b, var1_data, Time)
    used = send_data(id_var3, var1_data, Time)
  enddo
  call diag_manager_end(Time)

  call check_output()
  call fms_end()

  contains
    subroutine check_output()
      type(FmsNetcdfFile_t) :: fileobj !< FMS2io fileobj
      character(len=256) :: buffer !< Buffer to read stuff into
      character(len=256) :: associated_files !< Expected associated files attribute
      character(len=256) :: cell_measures !< Expected cell measures attribute

      ! Check that the static_file.nc was created and it contains the area variables
      ! defined in the diag_table.yaml
      if (.not. open_file(fileobj, "static_file.nc", "read")) &
        call mpp_error(FATAL, "static_file.nc was not created by the diag manager!")
      if (.not. variable_exists(fileobj, "area1")) &
        call mpp_error(FATAL, "area2 is not in static_file.nc")
      if (.not. variable_exists(fileobj, "area2")) &
        call mpp_error(FATAL, "area2 is not in static_file.nc")
      call close_file(fileobj)

      ! Check that file1.nc exists, that it contains the associated files attribute and it is correct,
      ! that the variables exists and it contains the cell_measures attributes
      if (.not. open_file(fileobj, "file1.nc", "read")) &
        call mpp_error(FATAL, "file1.nc was not created by the diag manager!")

      associated_files = "area1: static_file.nc area2: static_file.nc"
      call get_global_attribute(fileobj, "associated_files", buffer)
      call compare_answers("associated_files", associated_files, buffer)

      cell_measures = "area: area1"
      call get_variable_attribute(fileobj, "var1", "cell_measures", buffer)
      call compare_answers("cell_measures", cell_measures, buffer)

      cell_measures = "area: area2"
      call get_variable_attribute(fileobj, "var2", "cell_measures", buffer)
      call compare_answers("cell_measures", cell_measures, buffer)

      cell_measures = "area: area1"
      call get_variable_attribute(fileobj, "var2b", "cell_measures", buffer)
      call compare_answers("cell_measures", cell_measures, buffer)

      call close_file(fileobj)

      if (.not. open_file(fileobj, "file2.nc", "read")) &
        call mpp_error(FATAL, "file2.nc was not created by the diag manager!")
      cell_measures = "area: area3"
      call get_variable_attribute(fileobj, "var3", "cell_measures", buffer)
      call compare_answers("cell_measures", cell_measures, buffer)
      call close_file(fileobj)
    end subroutine check_output

    subroutine compare_answers(label, expected_answer, answer)
      character(len=*), intent(in) :: label
      character(len=*), intent(in) :: expected_answer
      character(len=*), intent(in) :: answer

      if (trim(answer) .ne. trim(expected_answer)) then
        call mpp_error(FATAL, "The "//trim(label)//" attribute is not the expected result! "//&
                              trim(answer)//" does not equal "//trim(expected_answer))
      endif
    end subroutine
end program
