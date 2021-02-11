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

!> @brief  This programs tests calls to fms2_io get_valid and is_valid

program test_get_is_valid

use   fms2_io_mod, only: fms2_io_init, open_file, FmsNetcdfFile_t, valid_t, &
                         close_file, is_valid, get_valid
use   mpp_mod    , only: mpp_init, mpp_exit, mpp_npes, mpp_get_current_pelist, &
                         mpp_root_pe, mpp_pe, FATAL, mpp_error, mpp_sync
use   netcdf     , only: nf90_create, nf90_def_var, nf90_put_att, nf90_enddef, &
                         nf90_close, nf90_clobber, nf90_64bit_offset, nf90_double

use, intrinsic :: iso_fortran_env
use platform_mod

implicit none

type(valid_t) :: valid_type               !> Fms2io valid object type
type(FmsNetcdfFile_t) :: fileobj          !> Fms2io file obj
real(kind=r8_kind) :: sst(10,10)           !> Data
integer :: mask_in(10,10)                 !> Array defining if sst(:,:) is valid (1) or not (0)
integer :: ncid, varid, err               !> Netcdf integers needed
integer, dimension(:), allocatable :: pes !> Current pelist

call mpp_init()
call fms2_io_init()

allocate(pes(mpp_npes()))
call mpp_get_current_pelist(pes)

if (mpp_root_pe() .eq. mpp_pe()) then
!> Create your own file with a variable to test with:
   err = nf90_create('test_file.nc', ior(nf90_clobber, nf90_64bit_offset), ncid)
   err = nf90_def_var(ncid, 'sst', nf90_double,varid)
   err = nf90_put_att(ncid, varid, "_FillValue", real(999,kind=r8_kind))
   err = nf90_put_att(ncid, varid, "missing_value", real(999,kind=r8_kind))
   err = nf90_enddef(ncid)
   err = nf90_close(ncid)
endif

!> Wait for the root pe to catch up
call mpp_sync()

!> Open the file and set up the valid type
if (open_file(fileobj, "test_file.nc", "read", pelist=pes)) then
   valid_type = get_valid(fileobj, "sst")
   call close_file(fileobj)
else
   call mpp_error(FATAL, "test_get_valid: Error opening file")
endif

deallocate(pes)

!> Error checking:
if (.not. valid_type%has_fill) call mpp_error(FATAL, "test_get_valid: the fill valid_type%has_fill is not .true.")
if (valid_type%fill_val .ne. real(999,kind=r8_kind)) call mpp_error(FATAL, &
        "test_get_valid: the fill value is not correct")

if (.not. valid_type%has_missing) call mpp_error(FATAL, "test_get_valid: the fill valid_type%has_fill is not .true.")
if (valid_type%missing_val .ne. real(999,kind=r8_kind)) call mpp_error(FATAL, &
        "test_get_valid: the fill value is not correct")

if (.not. valid_type%has_range) call mpp_error(FATAL, "test_get_valid: the valid_type%has_range is not .true.")
if (.not. valid_type%has_max) call mpp_error(FATAL, "test_get_valid: the valid_type%has_max is not .true.")
if (valid_type%has_min) call mpp_error(FATAL, "test_get_valid: the valid_type%has_min is .true.")

!> Now test is_valid

!> Create some fake data
mask_in = 1
sst = real(0, kind=r8_kind)
sst(1,1) = real(999, kind=r8_kind)

!> Check where the data is valid
!> Set the mask_in to 0 wherever sst is valid
where(is_valid(sst,valid_type)) mask_in = 0

!> Error checking:
!> The sum of mask_in should be 1 because only one value is equal to the fill_value
if (sum(mask_in) .ne. 1) call mpp_error(FATAL, "test_is_valid: the mask should have a 1")

call mpp_exit()

end program test_get_is_valid
