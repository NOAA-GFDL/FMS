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
!! @brief Tests the read_input_nml subroutine in four difference scenarios
!! @author Colin Gladue
!! @email gfdl.climate.model.info@noaa.gov

program test_read_input_nml

  use mpp_mod, only : mpp_init, mpp_init_test_peset_allocated
  use mpp_mod, only : mpp_error, FATAL, NOTE
  use mpp_mod, only : read_input_nml, mpp_get_current_pelist_name
  use mpp_mod, only : input_nml_file
#include<file_version.h>

character(len=200) :: line !< Storage location of lines read from the input nml
character(len=128) :: filename !< Name of input nml file to be read
integer :: stat !< IOSTAT output integer
integer :: n, m !< Loop counting variable
integer :: current_pelist_name_len_plus1 !< Current pelist name length plus 1
integer :: ierr !< used by MPI_FINALIZE
character(len=:), allocatable :: toobig !< String passed as argument into read_input_nml that is
                                        !!larger than pelis_name and should raise an error

namelist /test_read_input_nml_nml/ test_numb

open(10, file="test_numb.nml", form="formatted", status="old")
read(10, nml = test_read_input_nml_nml)
close(10)

call mpp_init(test_level=mpp_init_test_peset_allocated)

if (test_numb == 1 .or. test_numb == 2 .or. test_numb == 4) then
  ! Test 1: Tests the subroutine on a valid input nml full of data,
  ! with no arguments passed to read_input_nml()
  ! Test 2: Tests the subroutine on a valid input nml full of data,
  ! with a string passed to read_input_nml() in order to read a different nml
  ! Test 4: Tests the subroutine on a valid empty input nml,
  ! with no arguments passed to read_input_nml()
  if (test_numb == 1) then
    filename = "input.nml"
    call read_input_nml()
  else if (test_numb == 4) then
    filename = "input_blank.nml"
    call read_input_nml("blank")
  else if (test_numb == 2) then
    filename = "input_alternative.nml"
    call read_input_nml("alternative")
  end if
  open(1, file=filename, iostat=stat) ! Open input nml or alternative
  n = 1
  do
    read(1, '(A)', iostat=stat) line
    if (stat.eq.-1) then
      exit
    end if
    if (input_nml_file(n).ne.line) then
      call mpp_error(FATAL, "data resident in ./input.nml does not match&
                             & that read into input_nml_file by read_input_nml")
    end if
    n = n + 1
  end do
  close(1)

else if (test_numb.eq.3) then
  ! Test 3: Tests with an invalid pelist_name_in pass as an argument. An invalid
  ! pelist_name_in would be one who's size is greater than local pelist_name
  current_pelist_name_len_plus1 = LEN(mpp_get_current_pelist_name())
  allocate(character(len=current_pelist_name_len_plus1) :: toobig)
  call read_input_nml(pelist_name_in=toobig)
                                                          ! Call read_input_nml
                                                          ! with the local
                                                          ! pelist_name plus an
                                                          ! extra character "e"
  deallocate(toobig)
end if

call MPI_FINALIZE(ierr)

end program test_read_input_nml
