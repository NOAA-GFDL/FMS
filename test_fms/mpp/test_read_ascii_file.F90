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
!! @brief Tests the read_ascii_file subroutine
!! @author Colin Gladue
!! @email gfdl.climate.model.info@noaa.gov

program test_read_ascii_file

  use mpp_mod, only : mpp_init, mpp_init_test_peset_allocated
  use mpp_mod, only : mpp_error, FATAL, NOTE
  use mpp_mod, only : read_ascii_file, INPUT_STR_LENGTH, get_ascii_file_num_lines
  use mpp_mod, only : mpp_get_current_pelist, mpp_npes

  character(len=256), dimension(:), allocatable :: test_array !< Content array
  character(len=256) :: filename !< Name of ascii file to be read
  character(len=256) :: filename2 !< Name of alternative ascii file to be read
  character(len=256) :: line !< Content of a line of the read ascii file
  integer :: num_lines !< Number of lines in the ascii file
  integer :: stat !< IOSTATUS from the read method
  integer, allocatable :: cur_pelist(:) !< PELIST is read into this variable
  integer :: ierr !< used by MPI_FINALIZE

  namelist /test_read_ascii_file_nml/ test_numb

  open(20, file="test_numb_ascii.nml", form="formatted", status="old")
  read(20, nml=test_read_ascii_file_nml)
  close(20)

  ! Tests not meant to raise errors
  call mpp_init(test_level=mpp_init_test_peset_allocated)
  if (test_numb == 1 .or. test_numb == 7 .or. test_numb == 8) then
    if (test_numb == 1) then
      filename = "input.nml"
      num_lines = get_ascii_file_num_lines(filename, INPUT_STR_LENGTH)
      allocate(test_array(num_lines))
      call read_ascii_file(filename, INPUT_STR_LENGTH, test_array)
    else if (test_numb == 7) then
      filename = "input.nml"
      num_lines = get_ascii_file_num_lines(filename, INPUT_STR_LENGTH)
      allocate(test_array(num_lines))
      allocate(cur_pelist(0:mpp_npes()-1))
      call mpp_get_current_pelist(cur_pelist)
      call read_ascii_file(filename, INPUT_STR_LENGTH, test_array, PELIST=cur_pelist)
    else if (test_numb == 8) then
      filename = "empty.nml"
      num_lines = get_ascii_file_num_lines(filename, INPUT_STR_LENGTH)
      allocate(test_array(num_lines))
      call read_ascii_file(filename, INPUT_STR_LENGTH, test_array)
    end if
    ! Content check
    open(2, file=filename, iostat=stat)
    do i=1, num_lines-1
      read(2, '(A)', iostat=stat) line
      if (stat.eq.-1) then
        call mpp_error(FATAL, "Problem reading the ascii file")
      end if
      if (test_array(i).ne.line) then
        call mpp_error(FATAL, "Content array variable does not&
                                   & match the ascii file content")
      end if
    end do
  ! Tests meant to raise errors
  else
    if (test_numb == 2) then
      filename = "input.nml"
      allocate(test_array(20))
      call read_ascii_file(filename, INPUT_STR_LENGTH, test_array)
    else if (test_numb == 3) then
      filename = "doesnotexist.txt"
      ! Need to pass in an exist file name below to avoid raising error on 
      ! get_ascii_file_num_lines call in order to get to the error in read_ascii_file
      filename2 = "input.nml"
      num_lines = get_ascii_file_num_lines(filename2, INPUT_STR_LENGTH)
      allocate(test_array(num_lines))
      call read_ascii_file(filename, INPUT_STR_LENGTH, test_array)
    else if (test_numb == 4) then
      filename = "input.nml"
      filename2 = "empty.nml"
      num_lines = get_ascii_file_num_lines(filename2, INPUT_STR_LENGTH)
      allocate(test_array(num_lines))
      call read_ascii_file(filename, INPUT_STR_LENGTH, test_array)
    else if (test_numb == 5) then
      filename = "input.nml"
      num_lines = get_ascii_file_num_lines(filename, INPUT_STR_LENGTH)
      allocate(test_array(num_lines))
      call read_ascii_file(filename, 0, test_array)
    else if (test_numb == 6) then
      filename = "input.nml"
      num_lines = get_ascii_file_num_lines(filename, INPUT_STR_LENGTH)
      allocate(test_array(num_lines-1))
      call read_ascii_file(filename, INPUT_STR_LENGTH, test_array)
    end if
  end if
  call MPI_FINALIZE(ierr)
end program test_read_ascii_file
