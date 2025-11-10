!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
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
