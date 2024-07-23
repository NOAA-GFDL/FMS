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

!> @brief  This programs tests diag manager when the init date is prepended to the file name
program test_prepend_date

  use fms_mod,          only: fms_init, fms_end, string
  use diag_manager_mod, only: diag_axis_init, send_data, diag_send_complete, diag_manager_set_time_end, &
                              register_diag_field, diag_manager_init, diag_manager_end, register_static_field, &
                              diag_axis_init
  use time_manager_mod, only: time_type, operator(+), JULIAN, set_time, set_calendar_type, set_date
  use mpp_mod,          only: FATAL, mpp_error, input_nml_file
  use fms2_io_mod,      only: FmsNetcdfFile_t, open_file, close_file, read_data, get_dimension_size
  use platform_mod,     only: r4_kind

  implicit none

  integer         :: id_var0, id_var2, id_var1 !< diag field ids
  integer         :: id_axis1         !< Id for axis
  logical         :: used             !< for send_data calls
  integer         :: ntimes = 48      !< Number of time steps
  real            :: vdata            !< Buffer to store the data
  type(time_type) :: Time             !< "Model" time
  type(time_type) :: Time_step        !< Time step for the "simulation"
  integer         :: i                !< For do loops
  logical         :: pass_diag_time = .True.   !< .True. if passing the time to diag_manager_init

  integer :: io_status !< Status when reading the namelist

  namelist / test_prepend_date_nml / pass_diag_time

  call fms_init

  read (input_nml_file, test_prepend_date_nml, iostat=io_status)
  if (io_status > 0) call mpp_error(FATAL,'=>test_prepend_date: Error reading input.nml')

  call set_calendar_type(JULIAN)

  ! This is going to be different from the base_date
  if (pass_diag_time) then
    call diag_manager_init(time_init=(/2, 1, 1, 0, 0, 0/))
  else
    call diag_manager_init()
  endif

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time (3600,0) !< 1 hour
  call diag_manager_set_time_end(set_date(2,1,3,0,0,0))

  id_axis1 = diag_axis_init('dummy_axis', (/real(1.)/), "mullions", "X")
  id_var0 = register_diag_field  ('ocn_mod', 'var0', Time)
  id_var2 = register_static_field ('ocn_mod', 'var2', (/id_axis1/))

  ! This is a different start_time, should lead to a crash if the variable is in the diag table yaml
  id_var1 = register_diag_field  ('ocn_mod', 'var1', set_date(2,1,6,0,0,0))

  used = send_data(id_var2, real(123.456))
  do i = 1, ntimes
    Time = Time + Time_step
    vdata = real(i)

    used = send_data(id_var0, vdata, Time) !< Sending data every hour!

    call diag_send_complete(Time_step)
  enddo

  call diag_manager_end(Time)

  call check_output()
  call fms_end

  contains

  !< @brief Check the diag manager output
  subroutine check_output()
    type(FmsNetcdfFile_t) :: fileobj     !< Fms2io fileobj
    integer               :: var_size    !< Size of the variable reading
    real(kind=r4_kind), allocatable     :: var_data(:) !< Buffer to read variable data to
    integer               :: j           !< For looping

    if (.not. open_file(fileobj, "00020101.test_non_static.nc", "read")) &
      call mpp_error(FATAL, "Error opening file:00020101.test_non_static.nc to read")

    call get_dimension_size(fileobj, "time", var_size)
    if (var_size .ne. 48) call mpp_error(FATAL, "The dimension of time in the file:test_0days is not the "//&
                                                 "correct size!")
    allocate(var_data(var_size))
    var_data = -999.99

    call read_data(fileobj, "var0", var_data)
    do j = 1, var_size
      if (var_data(j) .ne. real(j, kind=r4_kind)) call mpp_error(FATAL, "The variable data for var1 at time level:"//&
                                                          string(j)//" is not the correct value!")
    enddo

    call close_file(fileobj)

    if (.not. open_file(fileobj, "00020101.test_static.nc", "read")) &
      call mpp_error(FATAL, "Error opening file:00020101.test_static.nc to read")

    call read_data(fileobj, "var2", var_data(1))
    if (var_data(1) .ne. real(123.456, kind=r4_kind)) call mpp_error(FATAL, &
      "The variable data for var2 is not the correct value!")

    call close_file(fileobj)

  end subroutine check_output
end program test_prepend_date
