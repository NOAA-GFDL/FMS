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

!> @brief  This programs tests when flush_nc_file=.true.
program test_flush_nc_file
  use fms_mod, only: fms_init, fms_end
  use diag_manager_mod
  use mpp_mod
  use mpp_domains_mod
  use platform_mod,     only: r8_kind, r4_kind
  use time_manager_mod, only: time_type, set_calendar_type, set_date, JULIAN, set_time, OPERATOR(+)
  use fms2_io_mod
  use fms_diag_yaml_mod

  implicit none

  type(time_type)                    :: Time             !< Time of the simulation
  type(time_type)                    :: Time_step        !< Time_step of the simulation
  integer                            :: nx               !< Number of x points
  integer                            :: ny               !< Number of y points
  integer                            :: nz               !< Number of z points
  integer                            :: id_x             !< Axis id for the x dimension
  integer                            :: id_y             !< Axis id for the y dimension
  integer                            :: id_var1          !< Field id for 1st variable
  logical                            :: used             !< Dummy argument to send_data
  real,                  allocatable :: x(:)             !< X axis data
  real,                  allocatable :: y(:)             !< Y axis_data
  real,                  allocatable :: var1_data(:,:)   !< Data for variable 1
  integer                            :: i                !< For do loops

  call fms_init
  call set_calendar_type(JULIAN)
  call diag_manager_init

  nx = 360
  ny = 180

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
  do i = 1, 24
    Time = Time + Time_step
    var1_data = real(i)
    used = send_data(id_var1, var1_data, Time)
    if (mpp_pe() .eq. mpp_root_pe()) print *, "Calling send_data::", i
    call diag_send_complete(Time_step)
    call mpp_sync()
    ! The file should have been flushed by now
    call check_answers(i)
  enddo

  call diag_manager_end(Time)
  !call check_answers(i)
  call fms_end

  contains
    subroutine check_answers(time_level)
      integer, intent(in) :: time_level
      type(FmsNetcdfFile_t) :: fileobj
      integer                            :: ntimes
      real,   allocatable                :: vardata(:,:)

      if (mpp_pe() .ne. mpp_root_pe()) return
      if (.not. open_file(fileobj, "test_flush.nc", "read")) &
        call mpp_error(FATAL, "unable to open test_flush.nc for reading")

      call get_dimension_size(fileobj, "time", ntimes)
      if (ntimes .ne. time_level) call mpp_error(FATAL, "time is not the correct size::", time_level)

      allocate(vardata(nx,ny))

      call read_data(fileobj, "ua", vardata, unlim_dim_level=i)
      if (any(vardata .ne. time_level)) &
        call mpp_error(FATAL, "ua is not the expected result")

      call close_file(fileobj)
    end subroutine check_answers
end program test_flush_nc_file