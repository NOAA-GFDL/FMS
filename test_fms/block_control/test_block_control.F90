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

program test_block_control
  use fms_mod,              only: fms_init, fms_end
  use mpp_domains_mod,      only: domain2d, mpp_define_domains, mpp_get_compute_domain
  use block_control_mod,    only: block_control_type, define_blocks
  use mpp_mod,              only: mpp_pe, mpp_root_pe, mpp_error, FATAL
  use fms_string_utils_mod, only: string

  implicit none

  integer, parameter :: nx=96                !< Size of the x grid
  integer, parameter :: ny=96                !< Size of the y grid
  type(domain2d)     :: Domain               !< 2D domain
  integer            :: layout(2) = (/2, 3/) !< Layout of the domain
  type(block_control_type) :: my_block       !< Block control type
  integer            :: isc, iec, jsc, jec   !< Starting and ending index for the commute domain
  integer            :: expected_startingy   !< Expected starting y index for the current block
  integer            :: expected_endingy     !< Expected ending y index for the current block
  integer            :: ncy(3)               !< Size of the y for each block
  logical            :: message              !< Set to .True., to output the warning message
  integer            :: i                    !< For do loops

  call fms_init()
  message = .True. !< Needs to be .true. so that the error message can be printed
  call mpp_define_domains( (/1,nx,1,ny/), layout, Domain)
  call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
  call define_blocks ('testing_model', my_block, isc, iec, jsc, jec, kpts=0, &
                         nx_block=1, ny_block=3, message=message)

  !< Message will be set to .false. if the blocks are not uniform
  if (message) &
    call mpp_error(FATAL, "test_block_control::define_blocks did not output the warning message"//&
                          " about uneven blocks")

  !Expected size of each block for every PE
  ncy = (/11, 10, 11/)
  expected_endingy = jsc-1
  do i = 1, 3
    ! Check the starting and ending "x" indices for each block
    if (my_block%ibs(i) .ne. isc .or. my_block%ibe(i) .ne. iec) &
      call mpp_error(FATAL, "The starting and ending 'x' index for the "//string(i)//" block is not expected value!")

    ! Check the starting and ending "y" indices for each block
    expected_startingy = expected_endingy + 1
    expected_endingy = expected_startingy + ncy(i) - 1
    if (my_block%jbs(i) .ne. expected_startingy .or. my_block%jbe(i) .ne. expected_endingy) &
    call mpp_error(FATAL, "The starting and ending 'y' index for the "//string(i)//" block is not expected value!")
  enddo

  call fms_end()
end program
