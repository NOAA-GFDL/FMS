program test_global_att

use   fms2_io_mod
use   mpp_mod
use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64

type(FmsNetcdfFile_t) :: fileobj
real(kind=real64) :: buf_real64
real(kind=real64), dimension(2) :: buf_real64_1d
real(kind=real32) :: buf_real32
real(kind=real32), dimension(2)  :: buf_real32_1d
integer(kind=int32)   :: buf_int32
integer(kind=int32), dimension(2)    :: buf_int32_1d
integer(kind=int64)   :: buf_int64
integer(kind=int64), dimension(2)    :: buf_int64_1d

character(len=20) :: buf_str

call fms2_io_init

!> Write out the different possible global attributes
if (open_file(fileobj, "test_global_att.nc", "overwrite")) then
   call register_global_attribute(fileobj, "buf_real64", real(7., kind=real64))
   call register_global_attribute(fileobj, "buf_real64_1d", (/ real(7., kind=real64), real(9., kind=real64) /))

   call register_global_attribute(fileobj, "buf_real32", real(4., kind=real32))
   call register_global_attribute(fileobj, "buf_real32_1d", (/ real(4., kind=real32), real(6., kind=real32)/) )

   call register_global_attribute(fileobj, "buf_int32", int(3, kind=int32))
   call register_global_attribute(fileobj, "buf_int32_1d", (/ int(3, kind=int32), int(5, kind=int32) /) )

   call register_global_attribute(fileobj, "buf_int64", int(2, kind=int64))
   call register_global_attribute(fileobj, "buf_int64_1d", (/ int(2, kind=int64), int(4, kind=int64) /) )

   call register_global_attribute(fileobj, "buf_str", "some text")

   call close_file(fileobj)
else
   call mpp_error(FATAL, "test_global_att: error opening the file for writting")
endif

!> Read the global attributes back
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

if (trim(buf_str) /= "some text") call mpp_error(FATAL, "test_global_att: error reading buf_str")

end program
