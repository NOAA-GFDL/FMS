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

!> @brief  This programs tests test_get_mosaic_tile_grid.
program test_get_mosaic_tile_grid
#ifndef use_mpp_io
use   mpp_mod,         only: mpp_npes, mpp_pe, mpp_root_pe, mpp_error, FATAL
use   mpp_domains_mod, only: domain2d, mpp_define_mosaic, mpp_define_io_domain
use   fms2_io_mod,     only: get_mosaic_tile_grid
use   fms_mod,         only: fms_init, fms_end
use   netcdf,          only: nf90_create, nf90_def_dim, nf90_clobber, nf90_64bit_offset, &
                             nf90_def_var, nf90_enddef, nf90_put_var, nf90_close, nf90_char, &
                             NF90_NOERR
use   mpi,             only: mpi_barrier, mpi_comm_world

implicit none

type(domain2d)               :: domain             !< 6 tile cube sphere domain
integer, parameter           :: ntiles = 6         !< Number of tiles
integer                      :: nx=60              !< Number of grid points in x
integer                      :: ny=60              !< Number of grid points in y
integer, dimension(4,ntiles) :: global_indices     !< Global indices for each tile
integer, dimension(2,ntiles) :: layout             !< Domain layout for each title
integer, dimension(2)        :: io_layout          !< Io layout for each tile
integer, dimension(ntiles)   :: pe_start           !< Starting PE for the given domain
integer, dimension(ntiles)   :: pe_end             !< Ending PE for the given domain
integer                      :: npes               !< Number of PEs
character(len=256)           :: expected_grid_file !< Expected grid_file name
character(len=256)           :: grid_file          !< Grid file read with the get_mosaic_tile_grid call
character(len=256)           :: mosaic_file        !< Filename where the grid_file is read from
integer                      :: i                  !< Index
integer :: err, ncid, dim1d, dim2d, varid
character (len = 255), dimension (7) :: char_in = (/'C48_grid.tile1.nc', &
                                               'C48_grid.tile2.nc', &
                                               'C48_grid.tile3.nc', &
                                               'C48_grid.tile4.nc', &
                                               'C48_grid.tile5.nc', &
                                               'C48_grid.tile6.nc', &
                                               'C48_grid.tile7.nc'  /)
call fms_init

!< Create a sample grid_spec file with the needed information
if (mpp_pe() .eq. mpp_root_pe()) then
   err = nf90_create('grid_spec.nc', ior(nf90_clobber, nf90_64bit_offset), ncid)
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_get_mosaic_tile_grid: Error creating file")

   err = nf90_def_dim(ncid, 'str', 255, dim1d)
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_get_mosaic_tile_grid: Error creating dimension")

   err = nf90_def_dim(ncid, 'ntiles', 7, dim2d)
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_get_mosaic_tile_grid: Error creating dimension")

   err = nf90_def_var(ncid, 'gridfiles', nf90_char, (/dim1d, dim2d/), varid)
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_get_mosaic_tile_grid: Error defining variable")

   err = nf90_enddef(ncid)
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_get_mosaic_tile_grid: Error leaving define mode")

   err = nf90_put_var(ncid, varid, char_in, start=(/ 1, 1/), count=(/ 255, 7 /))
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_get_mosaic_tile_grid: Error adding variable")

   err = nf90_close(ncid)
   if (err .ne. NF90_NOERR) call mpp_error(FATAL, "test_get_mosaic_tile_grid: Error closing file")
endif

call mpi_barrier(mpi_comm_world, err)

npes = mpp_npes()

!< Create a 6 tile cube sphere domain
do i = 1,ntiles
  global_indices(:, i) = (/1, nx, 1, ny/)
  layout(:, i) = (/1, npes/ntiles/)
  pe_start(i) = (i-1)*(npes/ntiles)
  pe_end(i) = i*(npes/ntiles) - 1
enddo
io_layout(:) = 1

call create_atmosphere_domain((/nx, nx, nx, nx, nx, nx/), &
                              (/ny, ny, ny, ny, ny, ny/), &
                              global_indices, layout, pe_start, pe_end, &
                              io_layout, domain)

!< Try to read in the grid_spec.nc file to get the name of the grid file for whatever tile
!! you are in
mosaic_file = 'grid_spec.nc'
call get_mosaic_tile_grid(grid_file,mosaic_file, domain)

!> Check that the grid_file was read correctly

!> Because this domain has 6 tiles, and you are running on 6 pes. Each pe will be on its
!! own tile. Pe 0 will be in tile1, pe 1 will be in tile 2, etc ..

i = mpp_pe()+1
write(expected_grid_file,'(A,I1,A)') "INPUT/C48_grid.tile", i, ".nc"
if (trim(grid_file) .ne. trim(expected_grid_file)) then
  print *, "grid_file: ", trim(grid_file), " expected: ", trim(expected_grid_file)
  call mpp_error(FATAL, "test_get_mosaic_tile_grid: the grid_file read is not correct")
endif

call fms_end

contains

include "create_atmosphere_domain.inc"
#endif
end program
