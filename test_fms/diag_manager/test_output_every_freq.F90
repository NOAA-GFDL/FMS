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
                              register_diag_field, diag_manager_init, diag_manager_end, register_static_field, &
                              diag_axis_init
  use time_manager_mod, only: time_type, operator(+), JULIAN, set_time, set_calendar_type, set_date
  use mpp_mod,          only: FATAL, mpp_error
  use fms2_io_mod,      only: FmsNetcdfFile_t, open_file, close_file, read_data, get_dimension_size
  use block_control_mod,only: define_blocks, block_control_type

  implicit none

  integer           :: id_var0, id_var01, id_var1, id_var2 !< diag field ids
  integer           :: id_axis1, id_axis2                  !< Id for axis
  logical           :: used                                !< for send_data calls
  integer           :: ntimes = 48                         !< Number of time steps
  real, allocatable :: vdata(:,:)                          !< Buffer to store the data
  type(time_type)   :: Time                                !< "Model" time
  type(time_type)   :: Time_step                           !< Time step for the "simulation"
  integer           :: i, iblock                           !< For do loops
  integer           :: npoints = 20                        !< Number of points in the x and y direction
  real, allocatable :: dummy_axis_data(:)                  !< Variable to store axis "data"
  integer           :: is, ie, js, je                      !< Starting and ending index of each block
  logical           :: message                             !< Flag for outputting debug message
  type(block_control_type) :: my_block                     !< Variable storing the indices of each openmp block

  call fms_init
  call set_calendar_type(JULIAN)
  call diag_manager_init

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time (3600,0) !< 1 hour
  call diag_manager_set_time_end(set_date(2,1,3,0,0,0))

  allocate(dummy_axis_data(npoints))
  do i = 1, npoints
    dummy_axis_data(i) = i
  enddo

  id_axis1 = diag_axis_init('x_axis', dummy_axis_data, "mullions", "X")
  id_axis2 = diag_axis_init('y_axis', dummy_axis_data, "mullions", "Y")

  allocate(vdata(npoints, npoints))
  vdata = -999.99

  ! These are the same variable but one will called from an openmp region and one will not
  id_var0 = register_diag_field  ('ocn_mod', 'var0', (/id_axis1, id_axis2 /), Time)
  id_var01 = register_diag_field  ('ocn_mod', 'var0_openmp', (/id_axis1, id_axis2 /), Time)

  id_var1 = register_diag_field  ('ocn_mod', 'var1', (/id_axis1, id_axis2 /), Time)
  id_var2 = register_static_field ('ocn_mod', 'var2', (/id_axis1, id_axis2/))

  ! The npoints x npoints data are going to divided in blocks
  call define_blocks ('testing_model', my_block, 1, npoints, 1, npoints, kpts=0, &
                         nx_block=1, ny_block=4, message=message)

  vdata = real(123.456)
  used = send_data(id_var2, vdata)
  do i = 1, ntimes
    Time = Time + Time_step

!$OMP parallel do default(shared) private(iblock, is, ie, js, je)
    do iblock=1, 4
      ! These are already in 1-based arrays since there is only 1 PE
      is = my_block%ibs(iblock)
      ie = my_block%ibe(iblock)
      js = my_block%jbs(iblock)
      je = my_block%jbe(iblock)

      vdata(is:ie, js:je) = real(i) + real(iblock)/100.

      used = send_data(id_var01, vdata(is:ie, js:je), Time, is_in=is, js_in=js)
    enddo

    !< Sending data every hour!
    used = send_data(id_var0, vdata, Time)

    !< Sending data every two hours!
    !! Because data is sent at different frequencies var0 and var1 cannot be in the same file,
    !! if writting at every time step and a failure is expected
    if (mod(i,2) .eq. 0) used = send_data(id_var1, vdata, Time)

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
    real, allocatable     :: var0_data(:,:,:) !< Buffer to read variable data to
    real, allocatable     :: var01_data(:,:,:) !< Buffer to read variable data to

    if (.not. open_file(fileobj, "test_0days.nc", "read")) &
      call mpp_error(FATAL, "Error opening file:test_0days.nc to read")

    call get_dimension_size(fileobj, "time", var_size)
    if (var_size .ne. 48) call mpp_error(FATAL, "The dimension of time in the file:test_0days is not the "//&
                                                 "correct size!")
    allocate(var0_data(npoints, npoints, var_size))
    var0_data = -999.99
    allocate(var01_data(npoints, npoints, var_size))
    var01_data = -999.99

    call read_data(fileobj, "var0", var0_data)
    call read_data(fileobj, "var0_openmp", var01_data)

    if (sum(var01_data) .ne. sum(var0_data)) call mpp_error(FATAL, "Data is not the expected result")

    call close_file(fileobj)
  end subroutine check_output
end program test_output_every_freq
