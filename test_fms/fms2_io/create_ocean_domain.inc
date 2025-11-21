!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
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
