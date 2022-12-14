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
!> @defgroup get_cal_time_mod get_cal_time_mod
!> @ingroup time_manager
!> @brief Given a time increment as a real number, and base time and calendar
!!  as a character strings, returns time as a time_type variable.

!> @addtogroup get_cal_time_mod
!> @{
module get_cal_time_mod

use          fms_mod, only: error_mesg, FATAL, write_version_number, lowercase, &
                            check_nml_error, stdlog, &
                            mpp_pe, mpp_root_pe

use time_manager_mod, only: time_type, operator(+), operator(-), set_time, get_time, &
                            NO_CALENDAR, THIRTY_DAY_MONTHS, NOLEAP, JULIAN, GREGORIAN, &
                            set_calendar_type, get_calendar_type, set_date, &
                            get_date, days_in_month, valid_calendar_types
use mpp_mod,          only: input_nml_file
use platform_mod,     only: r4_kind, r8_kind

implicit none
private

logical :: module_is_initialized=.false. !> This module is initialized on
                                         !! the first call to get_cal_time
                                         !! because there is no constructor.
! <NAMELIST NAME="get_cal_time_nml">
! <DATA NAME="allow_calendar_conversion" TYPE="logical"  DEFAULT=".true.">
!   This sets the default value of the optional argument named "permit_calendar_conversion" of get_cal_time.
!   This namelist is deprecated as of the memphis release.
!   If calendar conversion is not desired, then it is recommended that permit_calendar_conversion
!   be present in the call to get_cal_time and that it be set to .false.
! </DATA>

logical :: allow_calendar_conversion=.true.

namelist / get_cal_time_nml / allow_calendar_conversion
! </NAMELIST>

interface get_cal_time
  module procedure get_cal_time_r4
  module procedure get_cal_time_r8
end interface

public :: get_cal_time

! Include variable "version" to be written to log file.
#include<file_version.h>

contains

function cut0(string)
character(len=256) :: cut0
character(len=*), intent(in) :: string
integer :: i

cut0 = string

do i=1,len(string)
  if(ichar(string(i:i)) == 0 ) then
    cut0(i:i) = ' '
  endif
enddo

return
end function cut0

#include "get_cal_time_r4.fh"
#include "get_cal_time_r8.fh"

end module get_cal_time_mod
!> @}
! close documentation grouping
