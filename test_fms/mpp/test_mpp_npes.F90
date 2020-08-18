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
program test_mpp_npes
!> @author Tom Robinson
!> @description This test first reads the environment variable NUM_PES which
!! must be set before running the program and should be set to the number
!! of PEs that are being used to run the program.  MPP_NPES is then called
!! and the value is compared to NUM_PES.  If the numbers are the same, then
!! the test is successful.
 use mpp_mod, only : mpp_init, mpp_init_test_peset_allocated, mpp_pe, mpp_npes, &
                     stderr, stdout, mpp_error, FATAL
  implicit none
 integer :: npes !< The total number of PEs returned from mpp_npes
 character (len=20) :: string_npes !< npes converted to a string
 integer :: test_npes !< The number of PEs used to run the program
 character(len=:),allocatable :: env_pes !< test_npes as a string
 integer :: len_env_var !< The length of env_pes
 integer :: ierr !< MPI error return
!> Initialize MPI
  call mpp_init(test_level=mpp_init_test_peset_allocated)
!> Find the number of PEs that were used to run the program
 CALL GET_ENVIRONMENT_VARIABLE("NUM_PES",length=len_env_var) !> Get the length
 allocate(character(len=len_env_var) :: env_pes) !> Allocate the string
 CALL GET_ENVIRONMENT_VARIABLE("NUM_PES",value=env_pes) !> Get the number
 read(env_pes,*) test_npes !> Convert to an integer
!> Get the total number of PEs
 npes = mpp_npes()
!> Check that the number of PEs matches the number used to run the program
 if (npes .ne. test_npes) then 
     write(string_npes,*) npes
     call mpp_error(FATAL, "The number of PEs used to run the program is "&
     //env_pes//" but mpp_npes returned "//trim(string_npes) )
 endif
!> Finalize MPI
 call MPI_FINALIZE (ierr)
end program test_mpp_npes
