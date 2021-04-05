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

!> @brief Initialize a cubed-sphere atomsphere domain.
subroutine create_atmosphere_domain(ni, nj, global_indices, layout, pe_start, pe_end, &
                                    io_layout, domain)

  integer, dimension(:), intent(in) :: ni
  integer, dimension(:), intent(in) :: nj
  integer, dimension(:,:), intent(in) :: global_indices
  integer, dimension(:,:), intent(in) :: layout
  integer, dimension(:), intent(in) :: pe_start
  integer, dimension(:), intent(in) :: pe_end
  integer, dimension(2), intent(in) :: io_layout
  type(domain2d), intent(inout) :: domain !< A cubed-sphere domain.

  integer, dimension(12) :: tile1
  integer, dimension(12) :: tile2
  integer, dimension(12) :: istart1
  integer, dimension(12) :: iend1
  integer, dimension(12) :: jstart1
  integer, dimension(12) :: jend1
  integer, dimension(12) :: istart2
  integer, dimension(12) :: iend2
  integer, dimension(12) :: jstart2
  integer, dimension(12) :: jend2
  integer :: ntiles
  integer :: num_contact
  integer, dimension(2) :: msize
  integer :: whalo
  integer :: ehalo
  integer :: shalo
  integer :: nhalo

  ntiles = 6
  num_contact = 12
  whalo = 2
  ehalo = whalo
  shalo = whalo
  nhalo = whalo
  if (size(pe_start) .ne. 6 .or. size(pe_end) .ne. 6 ) then
    call mpp_error(FATAL, "size of pe_start and pe_end should be 6.")
  endif
  if (size(global_indices,1) .ne. 4) then
    call mpp_error(FATAL, "size of first dimension of global_indices should be 4.")
  endif
  if (size(global_indices,2) .ne. 6) then
    call mpp_error(FATAL, "size of second dimension of global_indices should be 6.")
  endif
  if (size(layout,1) .ne. 2) then
    call mpp_error(FATAL, "size of first dimension of layout should be 2.")
  endif
  if (size(layout,2) .ne. 6) then
    call mpp_error(FATAL, "size of second dimension of layout should be 6.")
  endif
  if (size(ni) .ne. 6 .or. size(nj) .ne. 6) then
    call mpp_error(FATAL, "size of ni and nj should be 6.")
  endif

  !Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
  tile1(1) = 1
  tile2(1) = 2
  istart1(1) = ni(1)
  iend1(1) = ni(1)
  jstart1(1) = 1
  jend1(1) = nj(1)
  istart2(1) = 1
  iend2(1) = 1
  jstart2(1) = 1
  jend2(1) = nj(2)

  !Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
  tile1(2) = 1
  tile2(2) = 3
  istart1(2) = 1
  iend1(2) = ni(1)
  jstart1(2) = nj(1)
  jend1(2) = nj(1)
  istart2(2) = 1
  iend2(2) = 1
  jstart2(2) = nj(3)
  jend2(2) = 1

  !Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
  tile1(3) = 1
  tile2(3) = 5
  istart1(3) = 1
  iend1(3) = 1
  jstart1(3) = 1
  jend1(3) = nj(1)
  istart2(3) = ni(5)
  iend2(3) = 1
  jstart2(3) = nj(5)
  jend2(3) = nj(5)

  !Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
  tile1(4) = 1
  tile2(4) = 6
  istart1(4) = 1
  iend1(4) = ni(1)
  jstart1(4) = 1
  jend1(4) = 1
  istart2(4) = 1
  iend2(4) = ni(6)
  jstart2(4) = nj(6)
  jend2(4) = nj(6)

  !Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
  tile1(5) = 2
  tile2(5) = 3
  istart1(5) = 1
  iend1(5) = ni(2)
  jstart1(5) = nj(2)
  jend1(5) = nj(2)
  istart2(5) = 1
  iend2(5) = ni(3)
  jstart2(5) = 1
  jend2(5) = 1

  !Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
  tile1(6) = 2
  tile2(6) = 4
  istart1(6) = ni(2)
  iend1(6) = ni(2)
  jstart1(6) = 1
  jend1(6) = nj(2)
  istart2(6) = ni(4)
  iend2(6) = 1
  jstart2(6) = 1
  jend2(6) = 1

  !Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
  tile1(7) = 2
  tile2(7) = 6
  istart1(7) = 1
  iend1(7) = ni(2)
  jstart1(7) = 1
  jend1(7) = 1
  istart2(7) = ni(6)
  iend2(7) = ni(6)
  jstart2(7) = nj(6)
  jend2(7) = 1

  !Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
  tile1(8) = 3
  tile2(8) = 4
  istart1(8) = ni(3)
  iend1(8) = ni(3)
  jstart1(8) = 1
  jend1(8) = nj(3)
  istart2(8) = 1
  iend2(8) = 1
  jstart2(8) = 1
  jend2(8) = nj(4)

  !Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
  tile1(9) = 3
  tile2(9) = 5
  istart1(9) = 1
  iend1(9) = ni(3)
  jstart1(9) = nj(3)
  jend1(9) = nj(3)
  istart2(9) = 1
  iend2(9) = 1
  jstart2(9) = nj(5)
  jend2(9) = 1

  !Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
  tile1(10) = 4
  tile2(10) = 5
  istart1(10) = 1
  iend1(10) = ni(4)
  jstart1(10) = nj(4)
  jend1(10) = nj(4)
  istart2(10) = 1
  iend2(10) = ni(5)
  jstart2(10) = 1
  jend2(10) = 1

  !Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
  tile1(11) = 4
  tile2(11) = 6
  istart1(11) = ni(4)
  iend1(11) = ni(4)
  jstart1(11) = 1
  jend1(11) = nj(4)
  istart2(11) = ni(6)
  iend2(11) = 1
  jstart2(11) = 1
  jend2(11) = 1

  !Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
  tile1(12) = 5
  tile2(12) = 6
  istart1(12) = ni(5)
  iend1(12) = ni(5)
  jstart1(12) = 1
  jend1(12) = nj(5)
  istart2(12) = 1
  iend2(12) = 1
  jstart2(12) = 1
  jend2(12) = nj(6)
  msize(1) = maxval(ni(:)/layout(1,:)) + whalo + ehalo + 1
  msize(2) = maxval(nj(:)/layout(2,:)) + shalo + nhalo + 1
  call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, &
                         tile2, istart1, iend1, jstart1, jend1, istart2, iend2, &
                         jstart2, jend2, pe_start, pe_end, symmetry = .true., &
                         whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                         name=trim("Cubed-sphere"), memory_size=msize)
  call mpp_define_io_domain(domain, io_layout)
end subroutine create_atmosphere_domain
