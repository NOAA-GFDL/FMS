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

program test_global_att
#ifndef use_mpp_io
use   fms2_io_mod
use   mpp_mod
use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64

type(FmsNetcdfFile_t) :: fileobj  !< fms2io netcd file obj
real(kind=real64) :: buf_real64 !< real64 buffer
real(kind=real64), dimension(2) :: buf_real64_1d !< real64 1D buffer
real(kind=real32) :: buf_real32 !< real32 buffer
real(kind=real32), dimension(2)  :: buf_real32_1d !< real32 1D buffer
integer(kind=int32)   :: buf_int32 !< int32 buffer
integer(kind=int32), dimension(2)    :: buf_int32_1d !< int32 1D buffer
integer(kind=int64)   :: buf_int64 !< int64 buffer
integer(kind=int64), dimension(2)    :: buf_int64_1d !< int64 1D buffer

character(len=20) :: buf_str !< character buffer

call fms2_io_init
call mpp_init

!> Write out the different possible global attributes to a netcdf file
if (open_file(fileobj, "test_global_att.nc", "overwrite")) then
   call register_global_attribute(fileobj, "buf_real64", real(7., kind=real64))
   call register_global_attribute(fileobj, "buf_real64_1d", (/ real(7., kind=real64), real(9., kind=real64) /))

   call register_global_attribute(fileobj, "buf_real32", real(4., kind=real32))
   call register_global_attribute(fileobj, "buf_real32_1d", (/ real(4., kind=real32), real(6., kind=real32)/) )

   call register_global_attribute(fileobj, "buf_int32", int(3, kind=int32))
   call register_global_attribute(fileobj, "buf_int32_1d", (/ int(3, kind=int32), int(5, kind=int32) /) )

   call register_global_attribute(fileobj, "buf_int64", int(2, kind=int64))
   call register_global_attribute(fileobj, "buf_int64_1d", (/ int(2, kind=int64), int(4, kind=int64) /) )

   call register_global_attribute(fileobj, "buf_str", "some text"//char(0), str_len=10)

   call close_file(fileobj)
else
   call mpp_error(FATAL, "test_global_att: error opening the file for writting")
endif

!> Read the global attributes from the netcdf file
if (open_file(fileobj, "test_global_att.nc", "read")) then
   call get_global_attribute(fileobj, "buf_real64", buf_real64)
   call get_global_attribute(fileobj, "buf_real64_1d", buf_real64_1d)

   call get_global_attribute(fileobj, "buf_real32", buf_real32)
   call get_global_attribute(fileobj, "buf_real32_1d", buf_real32_1d)

   call get_global_attribute(fileobj, "buf_int32", buf_int32)
   call get_global_attribute(fileobj, "buf_int32_1d", buf_int32_1d)

   call get_global_attribute(fileobj, "buf_int64", buf_int64)
   call get_global_attribute(fileobj, "buf_int64_1d", buf_int64_1d)

   call get_global_attribute(fileobj, "buf_str", buf_str)

   call close_file(fileobj)
else
   call mpp_error(FATAL, "test_global_att: error opening the file for reading")
endif

!> Compares the values read with the expected values
if (buf_real64 /= real(7., kind=real64)) call mpp_error(FATAL, "test_global_att: error reading buf_real64")
if (buf_real64_1d(1) /= real(7., kind=real64) .or. buf_real64_1d(2) /= real(9., kind=real64)) &
     call mpp_error(FATAL, "test_global_att: error reading buf_real64_1d")

if (buf_real32 /= real(4., kind=real32)) call mpp_error(FATAL, "test_global_att: error reading buf_real32")
if (buf_real32_1d(1) /= real(4., kind=real32) .or. buf_real32_1d(2) /= real(6., kind=real32)) &
    call mpp_error(FATAL, "test_global_att: error reading buf_real32_1d")

if (buf_int32 /= int(3, kind=int32)) call mpp_error(FATAL, "test_global_att: error reading buf_int32")
if (buf_int32_1d(1) /= int(3, kind=int32) .or. buf_int32_1d(2) /= int(5, kind=int32)) &
    call mpp_error(FATAL, "test_global_att: error reading buf_int32_1d")

if (buf_int64 /= int(2, kind=int64)) call mpp_error(FATAL, "test_global_att: error reading buf_int64")
if (buf_int64_1d(1) /= int(2, kind=int64) .or. buf_int64_1d(2) /= int(4, kind=int64)) &
    call mpp_error(FATAL, "test_global_att: error reading buf_int64_1d")

if (trim(buf_str) /= "some text") then
    print *, "buf_str read in = ", trim(buf_str)
    call mpp_error(FATAL, "test_global_att: error reading buf_str")
endif
call mpp_exit()
#endif
end program
