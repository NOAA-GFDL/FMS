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
!! @brief Unit test for the stdlog and checking warning log functionality
!! @author Ryan Mulhall
!! @email gfdl.climate.model.info@noaa.gov
program test_stdlog
  use mpp_mod, only : mpp_init, mpp_init_test_peset_allocated, stdlog
  use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_error, FATAL, WARNING, NOTE
  use fms_mod, only : input_nml_file, check_nml_error

  integer :: log_unit !< Stores the returned standard log unit number
  integer :: warn_unit
  integer :: pe !< pe value
  integer :: root_pe !< root pe value
  integer :: ierr !< Error code

  integer :: test_num = 1
  namelist / test_stdlog_nml / test_num

  call mpp_init()

  read(input_nml_file, nml=test_stdlog_nml, iostat=io)
  ierr = check_nml_error(io, 'test_stdlog_nml')

  pe = mpp_pe()
  root_pe = mpp_root_pe()
  log_unit = stdlog()

  print * , "running test num: ", test_num

  select case(test_num)
  case(1)
    call test_write(.false.)
  case(2)
    call test_write(.true.)
  case(3)
    call check_write()
  end select

  call MPI_FINALIZE(ierr)

  contains

  subroutine test_write(do_error_test)
    logical, intent(in) :: do_error_test !< causes a fatal error to check output if true

    write(log_unit, *) "asdf"
    call mpp_error(NOTE, "test note output")
    call mpp_error(WARNING, "test warning output")
    if(do_error_test) call mpp_error(FATAL, "test fatal output")
  end subroutine test_write

  subroutine check_write()
    integer :: i, ref_num, u_num_warn
    character(len=128) :: line
    character(len=23), parameter :: warn_fname = 'warnfile.000000.out.old'
    character(len=128) :: ref_line(4)

    ref_line(1) = "NOTE from PE     0: MPP_DOMAINS_SET_STACK_SIZE: stack size set to    32768."
    ref_line(2) = "NOTE from PE     0: test note output"
    ref_line(3) = "WARNING from PE     0: test warning output"
    ref_line(4) = "FATAL from PE     0: test fatal output"
    open(newunit=u_num_warn, file=warn_fname, status="old", action="read")
    ref_num = 1
    do i=1, 7
      read(u_num_warn, '(A)') line
      if (trim(line) == '') cycle
      !! if we're testing with the old io enabled, we'll have some additional output we can skip
      if (trim(line) == 'NOTE from PE     0: MPP_IO_SET_STACK_SIZE: stack size set to     131072.') cycle
      if(trim(line) .ne. trim(ref_line(ref_num))) call mpp_error(FATAL, "warnfile output does not match reference data"&
                                                                //"reference line:"//ref_line(ref_num) &
                                                                //"output line:"//line)
      ref_num = ref_num + 1
    enddo
    close(u_num_warn)
  end subroutine check_write

end program test_stdlog

