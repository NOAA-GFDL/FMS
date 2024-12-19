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

program test_diag_attribute_add
  use platform_mod, only: r4_kind, r8_kind
  use mpp_mod, only: FATAL, mpp_error
  use fms_mod, only: fms_init, fms_end
  use diag_manager_mod, only: diag_axis_init, register_static_field, diag_send_complete, send_data
  use diag_manager_mod, only: register_diag_field, diag_field_add_attribute
  use diag_manager_mod, only: diag_manager_init, diag_manager_end, diag_manager_set_time_end
  use time_manager_mod, only: time_type, set_calendar_type, set_date, JULIAN, set_time, OPERATOR(+)
  use fms2_io_mod

  implicit none

  integer :: id_potatoes
  integer :: i
  type(time_type) :: Time
  type(time_type) :: Time_step
  logical :: used
  real(kind=r4_kind) :: fbuffer(2) = (/ 13., 14./)
  real(kind=r8_kind) :: dbuffer(2) = (/ 23., 24./)
  integer :: ibuffer(2) = (/ 551, 552/)
  character(len=20) :: cbuffer = "Hello World"

  call fms_init()
  call set_calendar_type(JULIAN)
  call diag_manager_init()

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time (3600*4,0)
  call diag_manager_set_time_end(set_date(2,1,2,0,0,0))

  id_potatoes = register_diag_field ('food_mod', 'potatoes', init_time=Time)
  call diag_field_add_attribute(id_potatoes, "real_32", fbuffer(1))
  call diag_field_add_attribute(id_potatoes, "real_32_1d", fbuffer)
  call diag_field_add_attribute(id_potatoes, "real_64", dbuffer(1))
  call diag_field_add_attribute(id_potatoes, "real_64_1d", dbuffer )
  call diag_field_add_attribute(id_potatoes, "integer", ibuffer(1))
  call diag_field_add_attribute(id_potatoes, "integer_1d", ibuffer)
  call diag_field_add_attribute(id_potatoes, "some_string", cbuffer)

  do i = 1, 6
    Time = Time + Time_step
    used = send_data(id_potatoes, real(103.201), Time)
    call diag_send_complete(Time_step)
  enddo

  call diag_manager_end(Time)

  call check_output()
  call fms_end()

  contains

  subroutine check_output()
    type(FmsNetcdfFile_t) :: fileobj !< FMS2io fileobj
    character(len=256) :: cbuffer_out !< Buffer to read stuff into
    integer :: ibuffer_out(2)
    real(kind=r4_kind) :: fbuffer_out(2)
    real(kind=r8_kind) :: dbuffer_out(2)

    if (.not. open_file(fileobj, "food_file.nc", "read")) &
      call mpp_error(FATAL, "food_file.nc was not created by the diag manager!")
    if (.not. variable_exists(fileobj, "potatoes")) &
      call mpp_error(FATAL, "potatoes is not in food_file.nc")

    !! Checking the string attributes
    call get_variable_attribute(fileobj, "potatoes", "some_string", cbuffer_out)
    if (trim(cbuffer_out) .ne. trim(cbuffer)) call mpp_error(FATAL, "some_string is not the expected attribute")

    !! Checking the integer attributes
    ibuffer_out = -999
    call get_variable_attribute(fileobj, "potatoes", "integer", ibuffer_out(1))
    if (ibuffer(1) .ne. ibuffer_out(1)) call mpp_error(FATAL, "integer is not the expected attribute")

    ibuffer_out = -999
    call get_variable_attribute(fileobj, "potatoes", "integer_1d", ibuffer_out)
    if (ibuffer(1) .ne. ibuffer_out(1)) call mpp_error(FATAL, "integer_1d is not the expected attribute")
    if (ibuffer(2) .ne. ibuffer_out(2)) call mpp_error(FATAL, "integer_1d is not the expected attribute")

    !! Checking the double attributes
    dbuffer_out = -999
    call get_variable_attribute(fileobj, "potatoes", "real_64", dbuffer_out(1))
    if (dbuffer(1) .ne. dbuffer_out(1)) call mpp_error(FATAL, "real_64 is not the expected attribute")

    dbuffer_out = -999
    call get_variable_attribute(fileobj, "potatoes", "real_64_1d", dbuffer_out)
    if (dbuffer(1) .ne. dbuffer_out(1)) call mpp_error(FATAL, "real_64_1d is not the expected attribute")
    if (dbuffer(2) .ne. dbuffer_out(2)) call mpp_error(FATAL, "real_64_1d is not the expected attribute")

    !! Checking the float attributes
    fbuffer_out = -999
    call get_variable_attribute(fileobj, "potatoes", "real_32", fbuffer_out(1))
    if (fbuffer(1) .ne. fbuffer_out(1)) call mpp_error(FATAL, "real_32 is not the expected attribute")

    fbuffer_out = -999
    call get_variable_attribute(fileobj, "potatoes", "real_32_1d", fbuffer_out)
    if (fbuffer(1) .ne. fbuffer_out(1)) call mpp_error(FATAL, "real_32_1d is not the expected attribute")
    if (fbuffer(2) .ne. fbuffer_out(2)) call mpp_error(FATAL, "real_32_1d is not the expected attribute")

    call close_file(fileobj)
  end subroutine check_output
end program test_diag_attribute_add