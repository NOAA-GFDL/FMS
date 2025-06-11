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

!> @brief  This programs tests calls to
!! set_filename_appendix(string_in) -> sets string_in as a module variable, filename_appendix, in fms2io
!! nullify_filename_appendix -> Resets the module variable, filename_appendix, in fms2io
!! get_instance_filename(name_in,name_out) -> Adds the filename_appendix to name_in and sets it as name_out
!! get_filename_appendix(string_out) -> gets the filename_appendix and sets it as string_out

program test_file_appendix

use   fms2_io_mod  , only: fms2_io_init, set_filename_appendix, get_filename_appendix, &
                           get_instance_filename, open_file, FmsNetcdfFile_t, close_file, &
                           nullify_filename_appendix
use   mpp_mod      , only: mpp_init, mpp_exit, mpp_error, FATAL

implicit none

logical :: test_passed                    !> Flag indicating if the test_passed
type(FmsNetcdfFile_t) :: fileobj          !> Fms2io file obj
character(len=50) :: buf                  !> String buffer

test_passed = .false.

call mpp_init()
call fms2_io_init()

call set_filename_appendix("nest01")

!< Call get_filename_appendix to check if the filename appendix was set up
buf = ""
call get_filename_appendix(buf)
if (trim(buf) .ne. "nest01") call mpp_error(FATAL, "get_filename_appendix is not set up")

!< Call get_instance_filename to check if it adds the filename_appendix to the original string
call get_instance_filename("nestfile.nc", buf)
if (trim(buf) .ne. "nestfile.nest01.nc") &
    call mpp_error(FATAL, "get_instance_filename does not add the filename appendix")

!< Call get instance_filename where the string send in has tile in the name
call get_instance_filename("nestfile.tile1.nc", buf)
if (trim(buf) .ne. "nestfile.nest01.tile1.nc") &
    call mpp_error(FATAL, "get_instance_filename does not add the filename appendix before the .tile")

!< Call get_instance_filename where the string send in has the pe number in the end
call get_instance_filename("nestfile.nc.0001", buf)
if (trim(buf) .ne. "nestfile.nest01.nc.0001") &
    call mpp_error(FATAL, trim(buf)//": get_instance_filename does not add the filename appendix before the .nc")

!< Call get_instance_filename where the string send in has the pe number in the
!end
call get_instance_filename("nestfile.tile1.nc.0001", buf)
if (trim(buf) .ne. "nestfile.nest01.tile1.nc.0001") &
    call mpp_error(FATAL, trim(buf)//": get_instance_filename does not add the filename appendix before the tile")

!< Call open_file to check if the filename_appendix was appended to the file created
if (open_file(fileobj, "nestfile.nc", "overwrite", is_restart=.true.)) then
    call close_file(fileobj)
else
    call mpp_error(FATAL, "Error opening nestfile.nc")
end if

inquire(file="nestfile.res.nest01.nc", exist=test_passed)
if (.not. test_passed) call mpp_error(FATAL, "Error creating nestfile.res.nest01.nc")

!< Call nullify_filename_appendix and check if the filename_appendix was removed
call nullify_filename_appendix()
call get_filename_appendix(buf)
if (trim(buf) .ne. "") call mpp_error(FATAL, "get_filename_appendix was not removed")

!< Call open file again and check if the expected file: "nestfile.res.nc" was created
if (open_file(fileobj, "nestfile.nc", "overwrite", is_restart=.true.)) then
    call close_file(fileobj)
else
    call mpp_error(FATAL, "Error opening nestfile.nc")
end if

test_passed = .false.
inquire(file="nestfile.res.nc", exist=test_passed)
if (.not. test_passed) call mpp_error(FATAL, "Error creating nestfile.res.nc")

call mpp_exit()

end program test_file_appendix
