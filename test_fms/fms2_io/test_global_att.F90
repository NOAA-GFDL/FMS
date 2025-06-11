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

!> @brief  This programs tests the "register_global_attribute" and
!! "get_global_attribute" interfaces
program test_global_att
use   fms2_io_mod
use   mpp_mod
use   platform_mod

type(FmsNetcdfFile_t) :: fileobj  !< fms2io netcd file obj
real(kind=r8_kind) :: buf_r8_kind !< r8_kind buffer
real(kind=r8_kind), dimension(2) :: buf_r8_kind_1d !< r8_kind 1D buffer
real(kind=r4_kind) :: buf_r4_kind !< r4_kind buffer
real(kind=r4_kind), dimension(2)  :: buf_r4_kind_1d !< r4_kind 1D buffer
integer(kind=i4_kind)   :: buf_i4_kind !< i4_kind buffer
integer(kind=i4_kind), dimension(2)    :: buf_i4_kind_1d !< i4_kind 1D buffer
integer(kind=i8_kind)   :: buf_i8_kind !< i8_kind buffer
integer(kind=i8_kind), dimension(2)    :: buf_i8_kind_1d !< i8_kind 1D buffer
character (len = 120), dimension(3) :: my_format !< Array of formats to try.

character(len=20) :: buf_str !< character buffer
integer :: i !< For Do loop

call fms2_io_init
call mpp_init

my_format(1) = '64bit'
my_format(2) = 'classic'
my_format(3) = 'netcdf4'

do i = 1, size(my_format)
!< Write out the different possible global attributes to a netcdf file
if (open_file(fileobj, "test_global_att_"//trim(my_format(i))//".nc", "overwrite", nc_format=my_format(i))) then
   call register_global_attribute(fileobj, "buf_r8_kind", real(7., kind=r8_kind))
   call register_global_attribute(fileobj, "buf_r8_kind_1d", (/ real(7., kind=r8_kind), real(9., kind=r8_kind) /))

   call register_global_attribute(fileobj, "buf_r4_kind", real(4., kind=r4_kind))
   call register_global_attribute(fileobj, "buf_r4_kind_1d", (/ real(4., kind=r4_kind), real(6., kind=r4_kind)/) )

   call register_global_attribute(fileobj, "buf_i4_kind", int(3, kind=i4_kind))
   call register_global_attribute(fileobj, "buf_i4_kind_1d", (/ int(3, kind=i4_kind), int(5, kind=i4_kind) /) )

   !< int8 is only supported with the "netcdf4" type
   if(i .eq. 3) then
     call register_global_attribute(fileobj, "buf_i8_kind", int(2, kind=i8_kind))
     call register_global_attribute(fileobj, "buf_i8_kind_1d", (/ int(2, kind=i8_kind), int(4, kind=i8_kind) /) )
   endif

   call register_global_attribute(fileobj, "buf_str", "some text"//char(0), str_len=10)

   call close_file(fileobj)
else
   call mpp_error(FATAL, "test_global_att: error opening the file for writing")
endif

!< Read the global attributes from the netcdf file
if (open_file(fileobj, "test_global_att_"//trim(my_format(i))//".nc", "read", nc_format=my_format(i))) then
   call get_global_attribute(fileobj, "buf_r8_kind", buf_r8_kind)
   call get_global_attribute(fileobj, "buf_r8_kind_1d", buf_r8_kind_1d)

   call get_global_attribute(fileobj, "buf_r4_kind", buf_r4_kind)
   call get_global_attribute(fileobj, "buf_r4_kind_1d", buf_r4_kind_1d)

   call get_global_attribute(fileobj, "buf_i4_kind", buf_i4_kind)
   call get_global_attribute(fileobj, "buf_i4_kind_1d", buf_i4_kind_1d)

   !< int8 is only supported with the "netcdf4" type
   if(i .eq. 3) then
     call get_global_attribute(fileobj, "buf_i8_kind", buf_i8_kind)
     call get_global_attribute(fileobj, "buf_i8_kind_1d", buf_i8_kind_1d)
   endif

   call get_global_attribute(fileobj, "buf_str", buf_str)

   call close_file(fileobj)
else
   call mpp_error(FATAL, "test_global_att: error opening the file for reading")
endif

!< Compares the values read with the expected values
if (buf_r8_kind /= real(7., kind=r8_kind)) call mpp_error(FATAL, "test_global_att-"// &
    & trim(my_format(i))//": error reading buf_r8_kind")
if (buf_r8_kind_1d(1) /= real(7., kind=r8_kind) .or. buf_r8_kind_1d(2) /= real(9., kind=r8_kind)) &
     call mpp_error(FATAL, "test_global_att-"//trim(my_format(i))//": error reading buf_r8_kind_1d")

if (buf_r4_kind /= real(4., kind=r4_kind)) call mpp_error(FATAL, "test_global_att-"// &
    & trim(my_format(i))//": error reading buf_r4_kind")
if (buf_r4_kind_1d(1) /= real(4., kind=r4_kind) .or. buf_r4_kind_1d(2) /= real(6., kind=r4_kind)) &
    call mpp_error(FATAL, "test_global_att-"//trim(my_format(i))//": error reading buf_r4_kind_1d")

if (buf_i4_kind /= int(3, kind=i4_kind)) call mpp_error(FATAL, "test_global_att-"// &
    & trim(my_format(i))//": error reading buf_i4_kind")
if (buf_i4_kind_1d(1) /= int(3, kind=i4_kind) .or. buf_i4_kind_1d(2) /= int(5, kind=i4_kind)) &
    call mpp_error(FATAL, "test_global_att-"//trim(my_format(i))//": error reading buf_i4_kind_1d")

!< int8 is only supported with the "netcdf4" type
if(i .eq. 3) then
   if (buf_i8_kind /= int(2, kind=i8_kind)) call mpp_error(FATAL, "test_global_att-"// &
       & trim(my_format(i))//": error reading buf_i8_kind")
   if (buf_i8_kind_1d(1) /= int(2, kind=i8_kind) .or. buf_i8_kind_1d(2) /= int(4, kind=i8_kind)) &
       & call mpp_error(FATAL, "test_global_att-"//trim(my_format(i))//": error reading buf_i8_kind_1d")
endif

if (trim(buf_str) /= "some text") then
    print *, "buf_str read in = ", trim(buf_str)
    call mpp_error(FATAL, "test_global_att-"//trim(my_format(i))//": error reading buf_str")
endif

enddo !< do i = 1, size(my_format)

call mpp_exit()
end program
