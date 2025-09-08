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

!> @brief Create an unstructured land-like domain.
subroutine create_land_domain(atmos_domain, nx, ny, ntiles, land_domain, npes_group)

  type(domain2d), intent(in) :: atmos_domain
  integer, intent(in) :: nx
  integer, intent(in) :: ny
  integer, intent(in) :: ntiles
  type(domainug), intent(inout) :: land_domain
  integer, intent(in) :: npes_group

  integer :: isc
  integer :: iec
  integer :: jsc
  integer :: jec
  integer :: isd
  integer :: ied
  integer :: jsd
  integer :: jed
  logical, dimension(:,:,:), allocatable :: lmask
  integer, dimension(:), allocatable :: npts_tile
  real :: rmask
  integer :: ntotal_land
  integer, dimension(:), allocatable :: grid_index
  integer, dimension(:), allocatable :: isl
  integer, dimension(:), allocatable :: iel
  integer, dimension(:), allocatable :: jsl
  integer, dimension(:), allocatable :: jel
  integer, dimension(:), allocatable :: ntiles_grid
  integer :: i
  integer :: j
  integer :: l
  integer :: n

  call mpp_get_compute_domain(atmos_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(atmos_domain, isd, ied, jsd, jed)
  allocate(npts_tile(ntiles))
  if (mpp_pe() .eq. mpp_root_pe()) then
    allocate(lmask(nx, ny, ntiles))
    lmask = .false.
    do n = 1, ntiles
      do j = 1, ny
        do i = 1, nx
          call random_number(rmask)
          if (rmask .gt. 0.5) then
            lmask(i, j, n) = .true.
          endif
        enddo
      enddo
      npts_tile(n) = count(lmask(:, :, n))
    enddo
    ntotal_land = sum(npts_tile)
    allocate(grid_index(ntotal_land))
    l = 0
    allocate(isl(0:mpp_npes()-1))
    allocate(iel(0:mpp_npes()-1))
    allocate(jsl(0:mpp_npes()-1))
    allocate(jel(0:mpp_npes()-1))
    call mpp_get_compute_domains(atmos_domain, xbegin=isl, xend=iel, ybegin=jsl, yend=jel)
    do n = 1, ntiles
      do j = 1, ny
        do i = 1, nx
          if (lmask(i, j, n)) then
            l = l + 1
            grid_index(l) = (j-1)*nx+i
          endif
        enddo
      enddo
    enddo
    deallocate(lmask)
    deallocate(isl)
    deallocate(iel)
    deallocate(jsl)
    deallocate(jel)
  endif
  call mpp_broadcast(npts_tile, ntiles, mpp_root_pe())
  if (mpp_pe() .ne. mpp_root_pe()) then
    ntotal_land = sum(npts_tile)
    allocate(grid_index(ntotal_land))
  endif
  call mpp_broadcast(grid_index, ntotal_land, mpp_root_pe())
  allocate(ntiles_grid(ntotal_land))
  ntiles_grid = 1
  call mpp_define_unstruct_domain(land_domain, atmos_domain, npts_tile, ntiles_grid, &
                                  mpp_npes(), npes_group, grid_index, name="Unstructured domain")
  deallocate(npts_tile)
  deallocate(grid_index)
  deallocate(ntiles_grid)
end subroutine create_land_domain
