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
!> @defgroup time_manager_mod time_manager_mod
!> @ingroup time_manager


!> @brief Convert a real number of seconds into a time_type variable.
!! @return A filled time type variable, and an error message if an
!!         error occurs.
function REAL_TO_TIME_TYPE(x,err_msg) result(t)
  real(FMS_TM_KIND),intent(in) :: x !< Number of seconds.
  character(len=*),intent(out),optional :: err_msg !< Error message.
  type(time_type) :: t
  integer :: days
  integer :: seconds
  integer :: ticks
  character(len=128) :: err_msg_local
  real(FMS_TM_KIND),parameter :: spd = real(86400)
  real(FMS_TM_KIND) :: tps
  real(FMS_TM_KIND) :: a
  tps = real(ticks_per_second, FMS_TM_KIND)
  a = x/spd
  days = safe_rtoi(a,do_floor)
  a = x - real(days, FMS_TM_KIND)*spd
  seconds = safe_rtoi(a,do_floor)
  a = (a - real(seconds, FMS_TM_KIND))*tps
  ticks = safe_rtoi(a,do_nearest)
  if (.not. set_time_private(seconds,days,ticks,t,err_msg_local)) then
    if (error_handler('function real_to_time_type',err_msg_local,err_msg)) then
      return
    endif
  endif
end function REAL_TO_TIME_TYPE

!> @brief Convert a floating point value to an integer value.
!! @return The integer value, using the input rounding mode.
function SAFE_RTOI(rval,mode) result(ival)
  real(FMS_TM_KIND),intent(in) :: rval !< A floating point value.
  integer,intent(in) :: mode !< A rouding mode (either "do_floor" or
                             !! "do_nearest")
  integer :: ival
  real(FMS_TM_KIND) :: big
  big = real(huge(ival), FMS_TM_KIND)
  if (rval .le. big .and. real(-1., FMS_TM_KIND)*rval .ge. real(-1., FMS_TM_KIND)*big) then
    if (mode .eq. do_floor) then
      ival = floor(rval)
    elseif (mode .eq. do_nearest) then
      ival = nint(rval)
    else
      call error_mesg("safe_rtoi","mode must be either do_floor" &
                      //" or do_nearest.",FATAL)
    endif
  else
    call error_mesg("safe_rtoi","input value cannot be safely" &
                   //" converted to a 32-bit integer.",FATAL)
  endif
end function SAFE_RTOI



!> @}
! close documentation grouping
