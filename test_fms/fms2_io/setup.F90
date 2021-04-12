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

module setup
use, intrinsic :: iso_fortran_env, only : error_unit
use mpi
use argparse
use mpp_mod
use mpp_domains_mod
use platform_mod

implicit none
private


character(len=8), parameter, public :: green = achar(27)//"[1;32m"
character(len=8), parameter, public :: pink = achar(27)//"[1;95m"
character(len=8), parameter, public :: color_end = achar(27)//"[0m"


!> @brief Test parameters.
type :: Params
  integer :: nx !< Number of grid points in the x-direction.
  integer :: ny !< Number of grid points in the y-direction.
  integer :: nz !< Number of grid points in the z-direction.
  integer, dimension(:,:), allocatable :: global_indices
  integer, dimension(:,:), allocatable :: layout
  integer, dimension(2) :: io_layout
  integer :: npes_group
  integer, dimension(:), allocatable :: pe_start
  integer, dimension(:), allocatable :: pe_end
  logical :: debug !< Flag controlling output verbosity.
endtype Params


public :: Params
public :: mpi_check
public :: init
public :: cleanup
public :: create_cubed_sphere_domain
public :: create_unstructured_domain
public :: create_tripolar_domain
public :: create_data


interface create_data
  module procedure create_data_double_1d
  module procedure create_data_double_2d
  module procedure create_data_int_2d
  module procedure create_data_double_3d
end interface create_data


contains


subroutine mpi_check(err)

  integer, intent(in) :: err

  if (err .ne. mpi_success) then
    call mpp_error(fatal, "Error: MPI returned error code: ", err)
  endif
end subroutine mpi_check


!> @brief Initial setup for tests.
subroutine init(test_params, ntiles)

  type(Params), intent(out) :: test_params
  integer, intent(in) :: ntiles

  integer :: my_rank
  integer :: npes
  type(Parser_t) :: parser

  character(len=32) :: buf
  integer :: err
  integer :: i

  !Allocate arrays.
  allocate(test_params%global_indices(4,ntiles))
  allocate(test_params%layout(2,ntiles))
  allocate(test_params%pe_start(ntiles))
  allocate(test_params%pe_end(ntiles))

  !Initialize mpp.
  call mpp_init()

  !Define command line arguments.
  parser = get_parser()
  call add_argument(parser, "-x", 'Number of points (per domain tile) in the "x" direction.', &
                    requires_val=.true.)
  call add_argument(parser, "-y", 'Number of points (per domain tile) in the "y" direction.', &
                    requires_val=.true.)
  call add_argument(parser, "-z", 'Number of points in the "z" direction.', &
                    requires_val=.true.)
  call add_argument(parser, "-d", "Have each rank occupy its own I/O domain tile.", &
                    requires_val=.false.)
  call add_argument(parser, "-v", "Increase verbosity.", requires_val=.false., &
                    long_name="--verbose")
  call parse_args(parser)

  my_rank = mpp_pe()
  npes = mpp_npes()
  if (my_rank .eq. 0) then
    write(error_unit,*) "Using ", npes, " MPI ranks."
  endif
  call mpi_barrier(mpi_comm_world, err)
  call mpi_check(err)

  !Set defaults.
  test_params%nx = 96
  test_params%ny = 96
  test_params%nz = 30
  test_params%io_layout(:) = 1
  test_params%npes_group = 1
  test_params%debug = .false.

  !Parse command line arguments.
  call get_argument(parser, "-x", buf)
  if (trim(buf) .ne. "not present") then
    read(buf,*) test_params%nx
  endif
  call get_argument(parser, "-y", buf)
  if (trim(buf) .ne. "not present") then
    read(buf,*) test_params%ny
  endif
  call get_argument(parser, "-z", buf)
  if (trim(buf) .ne. "not present") then
    read(buf,*) test_params%nz
  endif
  call get_argument(parser, "-d", buf)
  if (trim(buf) .eq. "present") then
    test_params%io_layout(2) = npes/ntiles
  else
    test_params%npes_group = npes/ntiles
  endif
  call get_argument(parser, "-v", buf)
  if (trim(buf) .eq. "present") then
    test_params%debug = .true.
  endif

  !Prepare for domains creation.
  call mpp_domains_init()
  do i = 1, ntiles
    test_params%global_indices(:, i) = (/1, test_params%nx, 1, test_params%ny/)
    test_params%layout(:, i) = (/1, npes/ntiles/)
    test_params%pe_start(i) = (i-1)*(npes/ntiles)
    test_params%pe_end(i) = i*(npes/ntiles) - 1
  enddo
end subroutine init


!> @brief Cleanup for tests.
subroutine cleanup(test_params)

  type(Params), intent(inout) :: test_params

  !Deallocate arrays.
  deallocate(test_params%global_indices)
  deallocate(test_params%layout)
  deallocate(test_params%pe_start)
  deallocate(test_params%pe_end)
  call mpp_domains_exit()
  call mpp_exit()
end subroutine cleanup


!> @brief Initialize a cubed-sphere domain.
subroutine create_cubed_sphere_domain(test_params, domain, io_layout)

  type(Params), intent(in) :: test_params !< Test parameters.
  type(domain2d), intent(inout) :: domain !< A cubed-sphere domain.
  integer, dimension(2), intent(in) :: io_layout

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
  integer, parameter :: ntiles = 6
  integer, parameter :: num_contact = 12
  integer, dimension(2) :: msize
  integer :: whalo
  integer :: ehalo
  integer :: shalo
  integer :: nhalo

  whalo = 2
  ehalo = whalo
  shalo = whalo
  nhalo = whalo
  if (size(test_params%pe_start) .ne. ntiles .or. size(test_params%pe_end) .ne. ntiles) then
    call mpp_error(fatal, "size of pe_start and pe_end should be 6.")
  endif
  if (size(test_params%global_indices, 1) .ne. 4) then
    call mpp_error(fatal, "size of first dimension of global_indices should be 4.")
  endif
  if (size(test_params%global_indices, 2) .ne. ntiles) then
    call mpp_error(fatal, "size of second dimension of global_indices should be 6.")
  endif
  if (size(test_params%layout, 1) .ne. 2) then
    call mpp_error(fatal, "size of first dimension of layout should be 2.")
  endif
  if (size(test_params%layout, 2) .ne. ntiles) then
    call mpp_error(fatal, "size of second dimension of layout should be 6.")
  endif

  !Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
  tile1(1) = 1
  tile2(1) = 2
  istart1(1) = test_params%nx
  iend1(1) = test_params%nx
  jstart1(1) = 1
  jend1(1) = test_params%ny
  istart2(1) = 1
  iend2(1) = 1
  jstart2(1) = 1
  jend2(1) = test_params%ny

  !Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
  tile1(2) = 1
  tile2(2) = 3
  istart1(2) = 1
  iend1(2) = test_params%nx
  jstart1(2) = test_params%ny
  jend1(2) = test_params%ny
  istart2(2) = 1
  iend2(2) = 1
  jstart2(2) = test_params%ny
  jend2(2) = 1

  !Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
  tile1(3) = 1
  tile2(3) = 5
  istart1(3) = 1
  iend1(3) = 1
  jstart1(3) = 1
  jend1(3) = test_params%ny
  istart2(3) = test_params%nx
  iend2(3) = 1
  jstart2(3) = test_params%ny
  jend2(3) = test_params%ny

  !Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
  tile1(4) = 1
  tile2(4) = 6
  istart1(4) = 1
  iend1(4) = test_params%nx
  jstart1(4) = 1
  jend1(4) = 1
  istart2(4) = 1
  iend2(4) = test_params%nx
  jstart2(4) = test_params%ny
  jend2(4) = test_params%ny

  !Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
  tile1(5) = 2
  tile2(5) = 3
  istart1(5) = 1
  iend1(5) = test_params%nx
  jstart1(5) = test_params%ny
  jend1(5) = test_params%ny
  istart2(5) = 1
  iend2(5) = test_params%nx
  jstart2(5) = 1
  jend2(5) = 1

  !Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
  tile1(6) = 2
  tile2(6) = 4
  istart1(6) = test_params%nx
  iend1(6) = test_params%nx
  jstart1(6) = 1
  jend1(6) = test_params%ny
  istart2(6) = test_params%nx
  iend2(6) = 1
  jstart2(6) = 1
  jend2(6) = 1

  !Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
  tile1(7) = 2
  tile2(7) = 6
  istart1(7) = 1
  iend1(7) = test_params%nx
  jstart1(7) = 1
  jend1(7) = 1
  istart2(7) = test_params%nx
  iend2(7) = test_params%nx
  jstart2(7) = test_params%ny
  jend2(7) = 1

  !Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
  tile1(8) = 3
  tile2(8) = 4
  istart1(8) = test_params%nx
  iend1(8) = test_params%nx
  jstart1(8) = 1
  jend1(8) = test_params%ny
  istart2(8) = 1
  iend2(8) = 1
  jstart2(8) = 1
  jend2(8) = test_params%ny

  !Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
  tile1(9) = 3
  tile2(9) = 5
  istart1(9) = 1
  iend1(9) = test_params%nx
  jstart1(9) = test_params%ny
  jend1(9) = test_params%ny
  istart2(9) = 1
  iend2(9) = 1
  jstart2(9) = test_params%ny
  jend2(9) = 1

  !Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
  tile1(10) = 4
  tile2(10) = 5
  istart1(10) = 1
  iend1(10) = test_params%nx
  jstart1(10) = test_params%ny
  jend1(10) = test_params%ny
  istart2(10) = 1
  iend2(10) = test_params%nx
  jstart2(10) = 1
  jend2(10) = 1

  !Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
  tile1(11) = 4
  tile2(11) = 6
  istart1(11) = test_params%nx
  iend1(11) = test_params%nx
  jstart1(11) = 1
  jend1(11) = test_params%ny
  istart2(11) = test_params%nx
  iend2(11) = 1
  jstart2(11) = 1
  jend2(11) = 1

  !Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
  tile1(12) = 5
  tile2(12) = 6
  istart1(12) = test_params%nx
  iend1(12) = test_params%nx
  jstart1(12) = 1
  jend1(12) = test_params%ny
  istart2(12) = 1
  iend2(12) = 1
  jstart2(12) = 1
  jend2(12) = test_params%ny
  msize(1) = test_params%nx/test_params%layout(1,1) + whalo + ehalo + 1
  msize(2) = test_params%ny/test_params%layout(2,1) + shalo + nhalo + 1
  call mpp_define_mosaic(test_params%global_indices, test_params%layout, domain, ntiles, &
                         num_contact, tile1, tile2, istart1, iend1, jstart1, jend1, &
                         istart2, iend2, jstart2, jend2, test_params%pe_start, &
                         test_params%pe_end, symmetry=.true., whalo=whalo, ehalo=ehalo, &
                         shalo=shalo, nhalo=nhalo, name=trim("Cubed-sphere"), &
                         memory_size=msize)
  call mpp_define_io_domain(domain, io_layout)
end subroutine create_cubed_sphere_domain


!> @brief Create an unstructured land-like domain.
subroutine create_unstructured_domain(cubed_sphere_domain, test_params, unstructured_domain)

  type(domain2d), intent(in) :: cubed_sphere_domain
  type(Params), intent(in) :: test_params
  type(domainug), intent(inout) :: unstructured_domain

  integer, parameter :: ntiles = 6
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

  call mpp_get_compute_domain(cubed_sphere_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(cubed_sphere_domain, isd, ied, jsd, jed)
  allocate(npts_tile(ntiles))
  if (mpp_pe() .eq. mpp_root_pe()) then
    allocate(lmask(test_params%nx, test_params%ny, ntiles))
    lmask = .false.
    do n = 1, ntiles
      do j = 1, test_params%ny
        do i = 1, test_params%nx
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
    call mpp_get_compute_domains(cubed_sphere_domain, xbegin=isl, xend=iel, ybegin=jsl, yend=jel)
    do n = 1, ntiles
      do j = 1, test_params%ny
        do i = 1, test_params%nx
          if (lmask(i, j, n)) then
            l = l + 1
            grid_index(l) = (j-1)*test_params%nx+i
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
  call mpp_define_unstruct_domain(unstructured_domain, cubed_sphere_domain, npts_tile, ntiles_grid, &
                                  mpp_npes(), test_params%npes_group, grid_index, &
                                  name="Unstructured domain")
  deallocate(npts_tile)
  deallocate(grid_index)
  deallocate(ntiles_grid)
end subroutine create_unstructured_domain


!> @brief Initialize an ocean style tripolar domain.
subroutine create_tripolar_domain(test_params, domain)

  type(Params), intent(in) :: test_params
  type(domain2d), intent(inout) :: domain !< A tripolar domain.

  integer :: whalo
  integer :: ehalo
  integer :: shalo
  integer :: nhalo

  whalo = 2
  ehalo = whalo
  shalo = whalo
  nhalo = whalo
  call mpp_define_domains((/1,test_params%nx,1,test_params%ny/), test_params%layout(:,1), &
                          domain, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          symmetry=.true., name="Tripolar Folded North Symmetry")
  call mpp_define_io_domain(domain, test_params%io_layout)
end subroutine create_tripolar_domain


subroutine create_data_double_1d(array, xsize)

  real(kind=r8_kind), dimension(:), allocatable, intent(inout) :: array
  integer, intent(in), optional :: xsize

  integer :: i

  if (present(xsize)) then
    if (allocated(array)) then
      deallocate(array)
    endif
    allocate(array(xsize))
  endif
  do i = 1, size(array)
    call random_number(array(i))
    array(i) = array(i) + real(mpp_pe(), kind=r8_kind)
  enddo
end subroutine create_data_double_1d


subroutine create_data_double_2d(array, sizes)

  real(kind=r8_kind), dimension(:,:), allocatable, intent(inout) :: array
  integer, dimension(2), intent(in), optional :: sizes

  integer :: i
  integer :: j

  if (present(sizes)) then
    if (allocated(array)) then
      deallocate(array)
    endif
    allocate(array(sizes(1), sizes(2)))
  endif
  do j = 1, size(array, 2)
    do i = 1, size(array, 1)
      call random_number(array(i,j))
      array(i,j) = array(i,j) + real(mpp_pe(), kind=r8_kind)
    enddo
  enddo
end subroutine create_data_double_2d


subroutine create_data_int_2d(array, sizes)

  integer(kind=i4_kind), dimension(:,:), allocatable, intent(inout) :: array
  integer, dimension(2), intent(in), optional :: sizes

  integer :: i
  integer :: j
  real(kind=r4_kind) :: r

  if (present(sizes)) then
    if (allocated(array)) then
      deallocate(array)
    endif
    allocate(array(sizes(1), sizes(2)))
  endif
  do j = 1, size(array, 2)
    do i = 1, size(array, 1)
      call random_number(r)
      array(i,j) = int(10.*r, kind=i4_kind) + mpp_pe()
    enddo
  enddo
end subroutine create_data_int_2d


subroutine create_data_double_3d(array, sizes)

  real(kind=r8_kind), dimension(:,:,:), allocatable, intent(inout) :: array
  integer, dimension(3), intent(in), optional :: sizes

  integer :: i
  integer :: j
  integer :: k

  if (present(sizes)) then
    if (allocated(array)) then
      deallocate(array)
    endif
    allocate(array(sizes(1), sizes(2), sizes(3)))
  endif
  do k = 1, size(array, 3)
    do j = 1, size(array, 2)
      do i = 1, size(array, 1)
        call random_number(array(i,j,k))
        array(i,j,k) = array(i,j,k) + real(mpp_pe(), kind=r8_kind)
      enddo
    enddo
  enddo
end subroutine create_data_double_3d


end module setup
