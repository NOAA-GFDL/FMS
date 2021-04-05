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
!> @author Tom Robinson
!> @description This test program is for testing the mpp_pe routine.  The
!! test initializes mpp, then calls the integer function MPP_PE.  MPP_PE
!! returns a positive integer value.  If the value returned by MPP_PE is
!! less than 0, then an error is thrown.
program test_mpp_p5

 use mpp_mod, only : mpp_init, mpp_init_test_peset_allocated, mpp_pe, mpp_npes, &
                     stderr, stdout, mpp_error, FATAL

 implicit none
 integer :: ierr
 integer :: total_pes !< The total number of PEs returned from mpp_npes
 integer :: my_mpp_pe !< The unique PE identifier


!> Initialize MPI to do what mpp_init would do

  call mpp_init(test_level=mpp_init_test_peset_allocated)

!> Get the total number of PEs
 total_pes = mpp_npes()
!> Get the PE number
 my_mpp_pe = mpp_pe()
!> Check that the PE number is between 0 and npes-1
 if (my_mpp_pe < 0) then
     call mpp_error(FATAL, "The PE number is less than 0")
 elseif (my_mpp_pe > total_pes-1) then
     call mpp_error(FATAL, "The PE number is greater than npes-1")
 endif
!> Finalize mpp
 call MPI_FINALIZE (ierr)
end program test_mpp_p5
