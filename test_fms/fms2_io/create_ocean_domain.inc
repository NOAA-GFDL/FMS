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

!> @brief Initialize an ocean style tripolar domain.
subroutine create_ocean_domain(nx, ny, npes, domain, layout, io_layout)

  integer, intent(in) :: nx !< Number of points in the x-direction.
  integer, intent(in) :: ny !< Number of points in the y-direction.
  integer, intent(in) :: npes !< Number of MPI ranks.
  type(domain2d), intent(inout) :: domain !< A tripolar domain.
  integer, dimension(2), intent(in) :: layout
  integer, dimension(2), intent(in) :: io_layout !< I/O layout.

  integer :: whalo
  integer :: ehalo
  integer :: shalo
  integer :: nhalo

  call mpp_domains_init()
! call mpp_define_layout((/1, nx, 1, ny/), npes, layout)
  whalo = 2
  ehalo = whalo
  shalo = whalo
  nhalo = whalo
  call mpp_define_domains((/1,nx,1,ny/), layout, domain, xflags=CYCLIC_GLOBAL_DOMAIN, &
                          yflags=FOLD_NORTH_EDGE, whalo=whalo, ehalo=ehalo, shalo=shalo, &
                          nhalo=nhalo, symmetry=.true., name="Tripolar Folded North Symmetry")
  call mpp_define_io_domain(domain, io_layout)
end subroutine create_ocean_domain
