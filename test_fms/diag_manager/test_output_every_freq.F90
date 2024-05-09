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

!> @brief  This programs tests diag manager when the file frequency is set to 0 days
program test_output_every_freq

  use fms_mod,          only: fms_init, fms_end, string
  use diag_manager_mod, only: diag_axis_init, send_data, diag_send_complete, diag_manager_set_time_end, &
                              register_diag_field, diag_manager_init, diag_manager_end
  use time_manager_mod, only: time_type, operator(+), JULIAN, set_time, set_calendar_type, set_date
  use mpp_mod,          only: FATAL, mpp_error
  use fms2_io_mod,      only: FmsNetcdfFile_t, open_file, close_file, read_data, get_dimension_size

  implicit none

  integer         :: id_var0, id_var1 !< diag field ids
  logical         :: used             !< for send_data calls
  integer         :: ntimes = 48      !< Number of time steps
  real            :: vdata            !< Buffer to store the data
  type(time_type) :: Time             !< "Model" time
  type(time_type) :: Time_step        !< Time step for the "simulation"
  integer         :: i                !< For do loops

  call fms_init
  call set_calendar_type(JULIAN)
  call diag_manager_init

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time (3600,0) !< 1 hour
  call diag_manager_set_time_end(set_date(2,1,3,0,0,0))

  id_var0 = register_diag_field  ('ocn_mod', 'var0', Time)
  id_var1 = register_diag_field  ('ocn_mod', 'var1', Time)

  do i = 1, ntimes
    Time = Time + Time_step
    vdata = real(i)

    used = send_data(id_var0, vdata, Time) !< Sending data every hour!
    if (mod(i,2) .eq. 0) used = send_data(id_var1, vdata, Time) !< Sending data every two hours!

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
    real, allocatable     :: var_data(:) !< Buffer to read variable data to
    integer               :: j           !< For looping

    if (.not. open_file(fileobj, "test_0days.nc", "read")) &
      call mpp_error(FATAL, "Error opening file:test_0days.nc to read")

    call get_dimension_size(fileobj, "time", var_size)
    if (var_size .ne. 48) call mpp_error(FATAL, "The dimension of time in the file:test_0days is not the "//&
                                                 "correct size!")
    allocate(var_data(var_size))
    var_data = -999.99

    call read_data(fileobj, "var0", var_data)
    do j = 1, var_size
      if (var_data(j) .ne. real(j)) call mpp_error(FATAL, "The variable data for var1 at time level:"//&
                                                          string(j)//" is not the correct value!")
    enddo

    call close_file(fileobj)
  end subroutine check_output
end program test_output_every_freq
