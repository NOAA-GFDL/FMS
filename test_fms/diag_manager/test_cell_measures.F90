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
  real(kind=r4_kind), allocatable :: var1_data(:)     !< Data for "var1"
  real(kind=r4_kind), allocatable :: area_data(:)     !< Data for the "area"
  integer                         :: id_area          !< Id of the "area" field
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
  id_area = register_static_field ('fun_mod', 'area', (/id_axis1/))
  id_var1 = register_diag_field  ('fun_mod', 'var1', (/id_axis1/), init_time=Time, area=id_area)

  used = send_data(id_area, area_data)

  do i = 1, 6
    Time = Time + Time_step
    call diag_send_complete(Time_step)
    used = send_data(id_var1, var1_data, Time)
  enddo
  call diag_manager_end(Time)

  call check_output()
  call fms_end()

  contains
    subroutine check_output()
      type(FmsNetcdfFile_t) :: fileobj !< FMS2io fileobj
      character(len=256) :: buffer !< Buffer to read stuff into

      ! Check that the static_file.nc was created and it contains the area attribute
      if (.not. open_file(fileobj, "static_file.nc", "read")) &
        call mpp_error(FATAL, "static_file.nc was not created by the diag manager!")
      if (.not. variable_exists(fileobj, "area")) &
        call mpp_error(FATAL, "area is not in static_file.nc")
      call close_file(fileobj)

      ! Check that file1.nc exists, that it contains the associated files attribute and it is correct,
      ! that the var1 exists and it contains the cell_measures attributes
      if (.not. open_file(fileobj, "file1.nc", "read")) &
        call mpp_error(FATAL, "file1.nc was not created by the diag manager!")

      call get_global_attribute(fileobj, "associated_files", buffer)
      if (trim(buffer) .ne. "area: static_file.nc") &
        call mpp_error(FATAL, "The associated_files global attribute is not the expected result! "//trim(buffer))

      call get_variable_attribute(fileobj, "var1", "cell_measures", buffer)
      if (trim(buffer) .ne. "area: area") &
        call mpp_error(FATAL, "The cell_measures attribute is not the expected result! "//trim(buffer))
      call close_file(fileobj)
    end subroutine check_output
end program
