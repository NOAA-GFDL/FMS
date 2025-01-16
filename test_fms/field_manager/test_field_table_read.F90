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
use fms2_io_mod,       only: set_filename_appendix
use ensemble_manager_mod, only: get_ensemble_size, ensemble_manager_init
use mpp_mod,           only : mpp_pe, mpp_root_pe, mpp_error, NOTE, FATAL, input_nml_file, mpp_npes, &
                              mpp_set_current_pelist, mpp_get_current_pelist

implicit none

integer :: nfields
integer :: nfields_expected
integer :: io_status
integer :: npes
integer, allocatable :: pelist(:)
integer :: ens_siz(6)
integer :: ensemble_id
character(len=10) :: text
integer, parameter :: default_test = 0
integer, parameter :: ensemble_test = 1
integer, parameter :: ensemble_same_yaml_test = 2

! namelist parameters
integer :: test_case = default_test

namelist / test_field_table_read_nml / test_case

call fms_init
read (input_nml_file, test_field_table_read_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>test_field_table_read: Error reading input.nml')

npes = mpp_npes()
allocate(pelist(npes))
call mpp_get_current_pelist(pelist)

nfields_expected = 4
select case (test_case)
case (ensemble_test, ensemble_same_yaml_test)
  if (npes .ne. 2) &
    call mpp_error(FATAL, "test_field_table_read:: this test requires 2 PEs!")

  call ensemble_manager_init
  ens_siz = get_ensemble_size()
  if (ens_siz(1) .ne. 2) &
    call mpp_error(FATAL, "This test requires 2 ensembles")

  nfields_expected = 3
  if (mpp_pe() .eq. 0) then
    !PEs 0 is the first ensemble
    ensemble_id = 1
    call mpp_set_current_pelist((/0/))
  else
    !PEs 1 is the second ensemble
    ensemble_id = 2
    call mpp_set_current_pelist((/1/))
    if (test_case .eq. ensemble_test) nfields_expected = 4
  endif

  write( text,'(a,i2.2)' ) 'ens_', ensemble_id
  call set_filename_appendix(trim(text))

end select

call field_manager_init(nfields)
print *, nfields
if (nfields .ne. nfields_expected) &
  call mpp_error(FATAL, "test_field_table_read:: The number fields returned is not the expected result")

select case (test_case)
case (ensemble_test)
  call mpp_set_current_pelist(pelist)
end select

call fms_end
end program test_field_table_read
