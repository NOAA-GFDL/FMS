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
!> @author Ryan Mulhall
!> @brief Simple test program for diag manager output yaml file.
!! Just checks output from previous test
program test_diag_out_yaml

use fms_mod,          only: fms_init, fms_end
use time_manager_mod, only: set_calendar_type, JULIAN, time_type
use mpp_mod,          only: mpp_root_pe, mpp_pe, mpp_error, FATAL

implicit none

type(time_type) :: time

call fms_init
call check_output_yaml
call fms_end

contains

!> checks output and reference file are equivalent
subroutine check_output_yaml
  integer :: i, un_out, un_ref
  integer, parameter :: yaml_len = 402
  character(len=128) :: out_yaml_line, ref_yaml_line
  character(len=17), parameter :: ref_fname = 'diag_out_ref.yaml'
  character(len=13), parameter :: out_fname = 'diag_out.yaml'
  if( mpp_root_pe() .ne. mpp_pe()) return
  open(newunit=un_out, file=out_fname, status="old", action="read")
  open(newunit=un_ref, file=ref_fname, status="old", action="read")
  do i=1, yaml_len
    read(un_out, '(A)') out_yaml_line
    read(un_ref, '(A)') ref_yaml_line
    if(out_yaml_line .ne. ref_yaml_line) call mpp_error(FATAL, 'diag_out.yaml does not match reference file.' &
                                                               //'reference line:'//ref_yaml_line &
                                                               //'output line:'//out_yaml_line)
  enddo
  close(un_out)
  close(un_ref)

end subroutine


end program