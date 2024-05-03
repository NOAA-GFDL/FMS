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

!> @brief  This programs tests the flexible timing capability in the modern diag_manager
program test_flexible_time
use   fms_mod,          only: fms_init, fms_end
use   time_manager_mod, only: set_date, time_type, increment_date, set_calendar_type, &
                              JULIAN, set_time, operator(+)
use   diag_manager_mod, only: diag_manager_init, diag_axis_init, register_diag_field, &
                              diag_manager_set_time_end, diag_send_complete, diag_manager_end, &
                              send_data
use   mpp_mod,          only: FATAL, mpp_error
use   platform_mod,     only: r8_kind

implicit none

real(kind=r8_kind)                :: var_data(2)      !< Dummy data
logical                           :: used             !< .True. if send_data was sucessful
type(time_type)                   :: Time             !< Time of the simulation
type(time_type)                   :: Time_step        !< Start time of the simulation
type(time_type)                   :: End_Time         !< End Time of the simulation
integer :: i
integer :: id_z, id_var

call fms_init()
call set_calendar_type(JULIAN)
call diag_manager_init

!< Starting time of the simulation
Time = set_date(2,1,1,3,0,0) !02/01/01 hour 3

!< Set up a dummy variable
id_z  = diag_axis_init('z',  (/1. ,2. /),  'point_Z', 'z', long_name='point_Z')
id_var = register_diag_field  ('atm_mod', 'var1', (/id_z/), Time, 'Var not domain decomposed', 'mullions')

!< Set up the end of the simulation (i.e 2 days long)
End_Time = set_date(2,1,3,3,0,0)
call diag_manager_set_time_end(End_Time)

!< Set up the simulation
Time_step = set_time (3600,0) !< 1 hour
do i=1,48
  var_data = real(i, kind=r8_kind)
  Time = Time + Time_step
  used = send_data(id_var, var_data, Time)
  call diag_send_complete(set_time(3600,0))
enddo

call diag_manager_end(End_Time)

call fms_end()

end program test_flexible_time
