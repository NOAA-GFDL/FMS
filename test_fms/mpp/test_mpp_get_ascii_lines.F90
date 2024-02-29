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

program test_get_ascii_lines
  use mpp_mod, only: mpp_init, mpp_init_test_logfile_init, get_ascii_file_num_lines, read_ascii_file
  use mpp_mod, only: input_nml_file
  use, intrinsic ::  iso_fortran_env, only: INT8

  implicit none

  interface assertEquals
     procedure assertEquals_int_int
  end interface assertEquals

  integer, parameter :: str_length = 256
  integer, parameter :: num_tests = 5
  integer(KIND=INT8) :: test_number
  integer, dimension(num_tests) :: my_num_lines=(/5,25,0,5,5/)
  integer :: f_num_lines, io_status
  integer :: ierr
  character(len=10), dimension(num_tests) :: file_name=(/"ascii_5   ",&
                                                         "ascii_25  ",&
                                                         "ascii_0   ",&
                                                         "ascii_skip",&
                                                         "ascii_long"/)
  character(len=15), dimension(num_tests) :: test_name=(/"5 line test    ",&
                                                         "25 line test   ",&
                                                         "0 line test    ",&
                                                         "blank line test",&
                                                         "long line test "/)
  character(len=str_length), allocatable :: file_contents(:)
  namelist /test_mpp_get_ascii_lines_nml/ test_number

  open(30, file="test_numb2.nml", form="formatted", status="old")
  read(30, nml = test_mpp_get_ascii_lines_nml)
  close(30)

  call mpp_init(test_level=mpp_init_test_logfile_init)
  my_num_lines(test_number) = my_num_lines(test_number)+1 !!!!! Please See Note At End of File
  f_num_lines = get_ascii_file_num_lines(trim(file_name(test_number)), str_length)
  call assertEquals(f_num_lines, my_num_lines(test_number), trim(test_name(test_number)))
  call MPI_FINALIZE(ierr)

contains

  subroutine assertEquals_int_int(tstval, expval, test_name)
    use, intrinsic ::  iso_fortran_env, only: ERROR_UNIT
    integer, intent(in) :: tstval
    integer, intent(in) :: expval
    character(len=*), intent(in) :: test_name



    if (tstval .eq. expval) then
       write (ERROR_UNIT, '(A," ",I0,": ",A," - ",A)') "Test number", test_number, test_name, "success"
    else
       write (ERROR_UNIT, '(A," ",I0,": ",A," - ",A)') "Test number", test_number, test_name, "failure"
       write (ERROR_UNIT, '("Expected """,I0,""" got """,I0,""".""")') expval, tstval
    end if
  end subroutine assertEquals_int_int
end program test_get_ascii_lines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPORTANT NOTE REGARDING LINE COUNT !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The current implmenetation of mpp (as of release 2020.02) requires that the
! number of lines returned from a file be one more than the actual number of
! lines contained within the file. This is set up so that there are no attempts
! to read nml information from an array of dimension zero, as that results in a
! segmentation fault. This workaround may be addressed in future versions of
! mpp, in which case this test must be adapted.
