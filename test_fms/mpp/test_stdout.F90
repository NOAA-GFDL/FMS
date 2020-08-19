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

!> @file
!! @brief Unit test for the STDOUT function
!! @author Colin Gladue
!! @email gfdl.climate.model.info@noaa.gov

program test_stdout
  use mpp_mod, only : mpp_init, mpp_init_test_peset_allocated, stdlog
  use mpp_mod, only : stdout, mpp_pe, mpp_root_pe
  use iso_fortran_env, only : OUTPUT_UNIT

  integer :: out_unit !< Stores the returned standard output unit number
  integer :: log_unit !< Stores the returned standard log unit number
  integer :: pe !< pe value
  integer :: root_pe !< root pe value
  integer :: ierr !< Error code

  call mpp_init(test_level=mpp_init_test_peset_allocated)

  out_unit = stdout()
  pe = mpp_pe()
  root_pe = mpp_root_pe()

  if (root_pe.eq.pe) then
    if (OUTPUT_UNIT.ne.out_unit) then
      write(OUTPUT_UNIT,*) "ERROR: stdout() did not return correct out_unit"
      write(OUTPUT_UNIT,*) "root_pe.eq.pe case"
    end if
  else
    log_unit = stdlog()
    if (log_unit.ne.out_unit) then
      write(OUTPUT_UNIT,*) "ERROR: stdout() did not return correct value"
      write(OUTPUT_UNIT,*) "root_pe.ne.pe case. stdout() should equal stdlog()."
    end if
  end if
  
  call MPI_FINALIZE(ierr)

end program test_stdout
