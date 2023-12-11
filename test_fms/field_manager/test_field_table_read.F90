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

!********************* Sample field table required: *********************
! "TRACER", "ocean_mod", "biotic1"
!           "diff_horiz", "linear", "slope=ok"
!           "longname", "biotic one" /
! "TRACER", "ocean_mod", "age_ctl" /
! "TRACER", "atmos_mod","radon"
!           "longname","radon-222"
!           "units","VMR*1E21"
!           "profile_type","fixed","surface_value=0.0E+00"
!           "convection","all"/
! "TRACER", "land_mod", "sphum"
!           "longname",     "specific humidity"
!            "units",        "kg/kg" /
!***********************************************************************

program test_field_table_read

use field_manager_mod, only: field_manager_init
use fms_mod,           only: fms_init, fms_end
use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_error, NOTE, FATAL

implicit none

integer :: nfields

call fms_init
call field_manager_init(nfields)
if (nfields .ne. 4) &
  call mpp_error(FATAL, "test_field_table_read:: The number fields returned is not the expected result")
call fms_end
end program test_field_table_read
