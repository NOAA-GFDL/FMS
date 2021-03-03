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

program main
use, intrinsic :: iso_fortran_env
use argparse
use mpi
use mpp_mod
use mpp_domains_mod
use fms2_io_mod
use platform_mod
implicit none

character(len=8), parameter :: green = achar(27)//"[1;32m"
character(len=8), parameter :: pink = achar(27)//"[1;95m"
character(len=8) :: color_end = achar(27)//"[0m"
integer, parameter :: atmos = 1
integer, parameter :: land = 2
integer, parameter :: ocean = 3
integer, parameter :: ntiles = 6

type(Parser_t) :: parser
logical, dimension(3) :: tests
integer :: nx
integer :: ny
integer :: nz
logical :: debug
character(len=32) :: buf
integer :: err
integer :: my_rank
integer :: npes
integer, dimension(4,ntiles) :: global_indices
integer, dimension(2,ntiles) :: layout
integer, dimension(2) :: ocn_layout
integer, dimension(2) :: io_layout
integer, dimension(2) :: ocn_io_layout
integer :: npes_group
integer, dimension(ntiles) :: pe_start
integer, dimension(ntiles) :: pe_end
type(domain2d) :: ocean_domain
type(domain2d) :: atmosphere_domain
type(domainug) :: land_domain
integer :: i

!Initialize mpp.
call mpp_init()
my_rank = mpp_pe()
npes = mpp_npes()
if (my_rank .eq. 0) then
  write(error_unit,*) "Using ", npes, " MPI ranks."
endif
call mpi_barrier(mpi_comm_world, err)
call mpi_check(err)

!Define command line arguments.
parser = get_parser()
call add_argument(parser, "-t", "Test to run (atmos,land,ocean).", requires_val=.true., &
                  long_name="--test")
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

!Set defaults.
tests(:) = .true.
nx = 96
ny = 96
nz = 30
io_layout(:) = 1
ocn_io_layout(:) = 1
npes_group = 1
debug = .false.

!Parse command line arguments.
call get_argument(parser, "-t",  buf)
if (trim(buf) .ne. "not present") then
  if (trim(buf) .eq. "atmos") then
    tests(land) = .false.
    tests(ocean) = .false.
  elseif (trim(buf) .eq. "land") then
    tests(atmos) = .false.
    tests(ocean) = .false.
  elseif (trim(buf) .eq. "ocean") then
    tests(atmos) = .false.
    tests(land) = .false.
  else
    write(error_unit, *) "unsupported value for -t.  Please use -h."
    stop 1
  endif
endif
call get_argument(parser, "-x", buf)
if (trim(buf) .ne. "not present") then
  read(buf,*) nx
endif
call get_argument(parser, "-y", buf)
if (trim(buf) .ne. "not present") then
  read(buf,*) ny
endif
call get_argument(parser, "-z", buf)
if (trim(buf) .ne. "not present") then
  read(buf,*) nz
endif
call get_argument(parser, "-d", buf)
if (trim(buf) .eq. "present") then
  io_layout(2) = npes/ntiles
  ocn_io_layout(2) = npes
else
  npes_group = npes/ntiles
endif
call get_argument(parser, "-v", buf)
if (trim(buf) .eq. "present") then
  debug = .true.
endif

!Prepare for domains creation.
call mpp_domains_init()
do i = 1,ntiles
  global_indices(:, i) = (/1, nx, 1, ny/)
  layout(:, i) = (/1, npes/ntiles/)
  pe_start(i) = (i-1)*(npes/ntiles)
  pe_end(i) = i*(npes/ntiles) - 1
enddo
ocn_layout = (/1, npes/)

call fms2_io_init()
!Run tests.
if (tests(atmos)) then
  if (mod(npes,ntiles) .ne. 0) then
    call mpp_error(FATAL, "the number of ranks must be divisible by ntiles.")
  endif
  call create_atmosphere_domain((/nx, nx, nx, nx, nx, nx/), &
                                (/ny, ny, ny, ny, ny, ny/), &
                                global_indices, layout, pe_start, pe_end, &
                                io_layout, atmosphere_domain)
  call atmosphere_restart_file(atmosphere_domain, nz, 3, debug)
endif
if (tests(land)) then
  if (.not. tests(atmos)) then
    if (mod(npes, ntiles) .ne. 0) then
      call mpp_error(FATAL, "the number of ranks must be divisible by ntiles.")
    endif
    call create_atmosphere_domain((/nx, nx, nx, nx, nx, nx/), &
                                  (/ny, ny, ny, ny, ny, ny/), &
                                  global_indices, layout, pe_start, pe_end, &
                                  io_layout, atmosphere_domain)
  endif
  call create_land_domain(atmosphere_domain, nx, ny, ntiles, land_domain, npes_group)
  call land_unstructured_restart_file(land_domain, nz, 2, debug)
  call land_compressed_restart_file(nz, 4, debug)
endif
if (tests(ocean)) then
  call create_ocean_domain(ntiles*nx, ntiles*ny, npes, ocean_domain, ocn_layout, &
                           ocn_io_layout)
  call ocean_restart_file(ocean_domain, nz, 5, debug)
endif


!Clean up.
call mpp_domains_exit()
call mpp_exit()
call destroy_parser(parser)


contains


include "create_atmosphere_domain.inc"
include "create_land_domain.inc"
include "create_ocean_domain.inc"
include "atmosphere_restart_file_test.inc"
include "land_unstructured_restart_file_test.inc"
include "land_compressed_restart_file_test.inc"
include "ocean_restart_file_test.inc"


subroutine mpi_check(err)

  integer, intent(in) :: err

  integer :: e

  if (err .ne. mpi_success) then
    call mpp_error(fatal, "Error: MPI returned error code: ", err)
  endif
end subroutine mpi_check

subroutine chksum_match(out_chksum, in_chksum, var_name, debug)

  integer(kind=i8_kind), intent(in) :: out_chksum
  integer(kind=i8_kind), intent(in) :: in_chksum
  character(len=*), intent(in) :: var_name
  logical, intent(in) :: debug

  character(len=32) :: res
  integer :: err

  if (in_chksum .eq. out_chksum) then
    res = trim(green)
  else
    res = trim(pink)
  endif
  if (debug) then
    write(error_unit, '(a,1x,i3.3,1x,a,1x,i20,1x,a,1x,i20,1x,a)') &
      trim(res)//"Rank", mpp_pe(), "- chksum("//trim(var_name)//"):", &
      out_chksum,"(written)", in_chksum,"(read) "//trim(color_end)
    call mpi_barrier(mpi_comm_world, err)
    call mpi_check(err)
  endif
  if (in_chksum .ne. out_chksum) then
    call mpp_error(fatal, "data that is read in does not match the" &
                   //" data that was written out.")
  endif
end subroutine chksum_match

end program main
