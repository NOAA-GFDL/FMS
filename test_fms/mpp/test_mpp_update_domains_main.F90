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
!> @author Jessica Liptak
!> @brief run mpp_domains tests on arrays of integers and real numbers
!! using different layouts and data precision
!> @note This test calls extensions of the routine test_halo_upate in test_mpp_domains.
program test_mpp_update_domains_main

  use test_mpp_update_domains_real, only : test_halo_update_r8, test_halo_update_r4
  use test_mpp_update_domains_real, only : test_subset_update_r8, test_subset_update_r4
  use test_mpp_update_domains_int , only : test_halo_update_i8, test_halo_update_i4
  use test_mpp_update_domains_int, only : test_subset_update_i8, test_subset_update_i4
  use mpp_mod, only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe
  use mpp_mod, only : mpp_set_stack_size
  use mpp_mod, only : mpp_init_test_requests_allocated
  use mpp_domains_mod, only : MPP_DOMAIN_TIME, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit
  use mpp_io_mod, only: mpp_io_init
  use platform_mod

  implicit none

  integer :: ierr
  integer :: stackmax=10000000
  !> Initialize mpp and mpp IO modules
  call mpp_init(test_level=mpp_init_test_requests_allocated)
  call mpp_domains_init(MPP_DOMAIN_TIME)
  call mpp_io_init()
  call mpp_domains_set_stack_size(stackmax)
  !> run the tests
  !> run the tests
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Calling real test_halo_update <-------------------'

  call test_halo_update_r4( 'Simple' ) !includes global field, global sum tests
  call test_halo_update_r4( 'Cyclic' )
  call test_halo_update_r4( 'Folded-north' ) !includes vector field test
  call test_halo_update_r4( 'Masked' ) !includes vector field test
  call test_halo_update_r4( 'Folded xy_halo' ) !
  call test_halo_update_r4( 'Simple symmetry' ) !includes global field, global sum tests
  call test_halo_update_r4( 'Cyclic symmetry' )
  call test_halo_update_r4( 'Folded-north symmetry' ) !includes vector field test
  call test_halo_update_r4( 'Folded-south symmetry' ) !includes vector field test
  call test_halo_update_r4( 'Folded-west symmetry' ) !includes vector field test
  call test_halo_update_r4( 'Folded-east symmetry' ) !includes vector field test

  call test_halo_update_r8( 'Simple' ) !includes global field, global sum tests
  call test_halo_update_r8( 'Cyclic' )
  call test_halo_update_r8( 'Folded-north' ) !includes vector field test
  call test_halo_update_r8( 'Masked' ) !includes vector field test
  call test_halo_update_r8( 'Folded xy_halo' ) !
  call test_halo_update_r8( 'Simple symmetry' ) !includes global field, global sum tests
  call test_halo_update_r8( 'Cyclic symmetry' )
  call test_halo_update_r8( 'Folded-north symmetry' ) !includes vector field test
  call test_halo_update_r8( 'Folded-south symmetry' ) !includes vector field test
  call test_halo_update_r8( 'Folded-west symmetry' ) !includes vector field test
  call test_halo_update_r8( 'Folded-east symmetry' ) !includes vector field test
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Calling integer test_halo_update <-------------------'
  !> @note mpp_update_domains vector interfaces are not defined for integer data types
  call test_halo_update_i4( 'Simple' ) !includes global field, global sum tests
  call test_halo_update_i4( 'Cyclic' )
  call test_halo_update_i4( 'Simple symmetry' ) !includes global field, global sum tests
  call test_halo_update_i4( 'Cyclic symmetry' )
  call test_halo_update_i8( 'Simple' ) !includes global field, global sum tests
  call test_halo_update_i8( 'Cyclic' )
  call test_halo_update_i8( 'Simple symmetry' ) !includes global field, global sum tests
  call test_halo_update_i8( 'Cyclic symmetry' )
  ! pe subset test
  !> @todo resolve issue. This test triggers an error in mpp_clock_begin called by mpp_update_domains
  !! cannot change pelist context of a clock.
  if (mpp_npes() .GE. 16) then
    if (mpp_pe() == mpp_root_pe()) &
      print *, '--------------------> Calling test_subset_update <-------------------'
    call test_subset_update_r8
    call test_subset_update_r4
    call test_subset_update_i8
    call test_subset_update_i4
  endif
  call mpp_domains_exit()
  !> Finalize mpp
  call MPI_FINALIZE(ierr)
end program test_mpp_update_domains_main