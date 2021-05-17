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

!> @brief  This programs tests calls to get_mosaic_ntiles, get_mosaic_ncontacts,
!! get_mosaic_grid_sizes, get_mosaic_contact
program test_mosaic

use mosaic2_mod, only : get_mosaic_ntiles, get_mosaic_ncontacts
use mosaic2_mod, only : get_mosaic_grid_sizes, get_mosaic_contact
use mpp_mod,     only : mpp_init, mpp_error, FATAL, mpp_sync, mpp_npes, mpp_get_current_pelist
use fms2_io_mod, only : open_file, close_file, FmsNetcdfFile_t
use fms2_io_mod, only : register_axis, register_field, write_data
use fms_mod,     only : fms_init, fms_end

implicit none

integer              :: ntiles         !< Number of tiles
integer              :: ncontacts      !< Number of contacts
integer              :: n              !< For do loops
integer, allocatable :: tile1(:)       !< tile number for first contact
integer, allocatable :: tile2(:)       !< tile number of the second contact
integer, allocatable :: nx(:), ny(:)   !< Number of x/y points for each tile
integer, allocatable :: istart1(:), iend1(:), jstart1(:), jend1(:) !< Indexes of first contact point
integer, allocatable :: istart2(:), iend2(:), jstart2(:), jend2(:) !< Indexes of second contact point
character(len=128)   :: mosaic_file    !< Mosaic filename
type(FmsNetcdfFile_t):: mosaic_fileobj !< Fileobj for the file read by the test
integer              :: answers(2, 8)  !< Expected results
integer, allocatable :: pes(:)         !< List of pes in the current pelist

call mpp_init()
call fms_init()

mosaic_file = "INPUT/ocean_mosaic.nc"
answers(1,:) = (/1440, 1440, 1, 1080, 1, 1, 1, 1080 /)
answers(2,:) = (/1, 720, 1080, 1080, 1440, 721, 1080, 1080 /)

allocate(pes(mpp_npes()))
call mpp_get_current_pelist(pes)

call create_files(pes)

!< Open the mosaic file
if(.not. open_file(mosaic_fileobj, mosaic_file, 'read', pelist=pes)) then
  call mpp_error(FATAL, 'test_mosaic: error in opening file '//trim(mosaic_file))
endif

ntiles = get_mosaic_ntiles(mosaic_fileobj)
ncontacts = get_mosaic_ncontacts(mosaic_fileobj)
allocate(nx(ntiles), ny(ntiles))
allocate(tile1(ncontacts), tile2(ncontacts) )
allocate(istart1(ncontacts), iend1(ncontacts), jstart1(ncontacts), jend1(ncontacts) )
allocate(istart2(ncontacts), iend2(ncontacts), jstart2(ncontacts), jend2(ncontacts) )

call get_mosaic_grid_sizes(mosaic_fileobj, nx, ny )
call get_mosaic_contact(mosaic_fileobj, tile1, tile2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2)

!< Compare with expected results:
if (ntiles .ne. 1) call mpp_error(FATAL, "ntiles is not equal to 1")

do n = 1, ntiles
   if (nx(n) .ne. 2880/2)  call mpp_error(FATAL, "nx is not the expected result")
   if (ny(n) .ne. 2160/2)  call mpp_error(FATAL, "ny is not the expected result")
end do

if (ncontacts .ne. 2) call mpp_error(FATAL, "ncontacts is not the expected result")
do n = 1, ncontacts
   if (istart1(n) .ne. answers(n,1)) call mpp_error(FATAL, "istart1 is not the expected result")
   if (iend1(n)   .ne. answers(n,2)) call mpp_error(FATAL, "iend1 is not the expected result")

   if (jstart1(n) .ne. answers(n,3)) call mpp_error(FATAL, "jstart1 is not the expected result")
   if (jend1(n)   .ne. answers(n,4)) call mpp_error(FATAL, "jend1 is not the expected result")

   if (istart2(n) .ne. answers(n,5)) call mpp_error(FATAL, "istart2 is not the expected result")
   if (iend2(n)   .ne. answers(n,6)) call mpp_error(FATAL, "iend2 is not the expected result")

   if (jstart2(n) .ne. answers(n,7)) call mpp_error(FATAL, "jstart2 is not the expected result")
   if (jend2(n)   .ne. answers(n,8)) call mpp_error(FATAL, "jend2 is not the expected result")
end do

deallocate(tile1, tile2, nx, ny)
deallocate(istart1, iend1, jstart1, jend1)
deallocate(istart2, iend2, jstart2, jend2)

call close_file(mosaic_fileobj)
call fms_end()

contains

subroutine create_files(pes)
   integer, intent(in)  :: pes(:)         !< List of pes

   type(FmsNetcdfFile_t):: fileobj        !< Fileobj for the files written by the test
   character(len=255)     :: str_array(2)   !< Array of strings because GNU

   if( open_file(fileobj, mosaic_file, 'overwrite', pelist=pes)) then
      call register_axis(fileobj, "ntiles", 1)
      call register_axis(fileobj, "ncontact", 2)
      call register_axis(fileobj, "string", 255)

      str_array(1) = "string"
      str_array(2) = "ncontact"
      call register_field(fileobj, "contacts", "char",  dimensions=str_array)
      call register_field(fileobj, "contact_index", "char",  dimensions=str_array)
      call register_field(fileobj, "gridfiles", "char", dimensions=(/"string", "ntiles"/))
      call register_field(fileobj, "gridtiles", "char", dimensions=(/"string", "ntiles"/))

      call write_data(fileobj, "gridfiles", (/"ocean_hgrid.nc"/))
      call write_data(fileobj, "gridtiles", (/"tile1"/))

      str_array(1) = "2880:2880,1:2160::1:1,1:2160"
      str_array(2) = "1:1440,2160:2160::2880:1441,2160:2160"
      call write_data(fileobj, "contact_index", str_array)
      call write_data(fileobj, "contacts", &
         & (/"ocean_mosaic:tile1::ocean_mosaic:tile1", "ocean_mosaic:tile1::ocean_mosaic:tile1" /))

      call close_file(fileobj)
   endif
   call mpp_sync()

   if( open_file(fileobj, "INPUT/ocean_hgrid.nc", "overwrite", pelist=pes)) then
      call register_axis(fileobj, "nx", 2880)
      call register_axis(fileobj, "ny", 2160)

      call close_file(fileobj)
   endif
   call mpp_sync()
end subroutine create_files

end program test_mosaic
