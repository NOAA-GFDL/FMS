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

program test_mpp_root_pe

  
  use mpp_mod, only:  mpp_init, FATAL, mpp_error, mpp_root_pe
  
  implicit none
  
  integer :: my_root_pe, test_root_pe, ierr


  call mpp_init( test_level=0 )
  
  
  !: check default root_pe=0
  my_root_pe = 0
  test_root_pe = mpp_root_pe()
  if( test_root_pe .ne. my_root_pe ) call mpp_error(FATAL, "mpp_root_pe does not equal root_pe")
  
  
  call MPI_FINALIZE(ierr)

  
end program test_mpp_root_pe
