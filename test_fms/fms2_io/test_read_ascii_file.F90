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
!! @brief Tests the read_ascii_file subroutine
!! @author Colin Gladue
!! @email gfdl.climate.model.info@noaa.gov

program test_read_ascii_file

  use mpp_mod, only : mpp_init
  use mpp_mod, only : mpp_error, FATAL, NOTE
  use fms2_io_mod, only : fms2_io_init, ascii_read

  character(len=:), dimension(:), allocatable :: test_array !< Content array
  character(len=256) :: filename !< Name of ascii file to be read
  character(len=256) :: filename2 !< Name of alternative ascii file to be read
  character(len=256) :: line !< Content of a line of the read ascii file
  integer :: num_lines !< Number of lines in the ascii file
  integer, dimension(2) :: stat !< IOSTATUS from the read method
  integer, allocatable :: cur_pelist(:) !< PELIST is read into this variable
  integer :: ierr !< used by MPI_FINALIZE

  call mpp_init()
  call fms2_io_init()
  filename = "ascii_test1"
  call ascii_read(filename, test_array)
  read(test_array(1), *) stat
  if (stat(1)*6 - (stat(2)+3) /= 13) call mpp_error(FATAL, "test_read_ascii: failed to read integers")
  read(test_array(2), *) num_lines
  if (num_lines-11 /= 12) call mpp_error(FATAL, "test_read_ascii: failed to read integer")
  read(test_array(3), *) line
  if (trim(line)//"wut" /= "forlendulawut") call mpp_error(FATAL, "test_read_ascii: failed to read string")
  call MPI_FINALIZE(ierr)
end program test_read_ascii_file
