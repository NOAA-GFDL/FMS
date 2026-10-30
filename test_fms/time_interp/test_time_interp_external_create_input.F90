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

program test_time_interp_external_create_input

use constants_mod,   only : constants_init
use mpp_mod,         only : mpp_init, mpp_exit, mpp_pe, mpp_root_pe, mpp_sync
use mpp_domains_mod, only : mpp_domains_init
use fms2_io_mod,     only : FmsNetcdfFile_t, fms2_io_init, open_file, close_file, write_data, register_axis
use fms2_io_mod,     only : register_field, unlimited, register_variable_attribute
use platform_mod,    only : r8_kind, r4_kind

implicit none

integer            :: i
character(len=128) :: filename='INPUT/aerosol.climatology.nc'
character(len=128) :: filename_1d='INPUT/solar_constant.nc'
character(len=128) :: fieldname='so4_anthro'
character(len=128) :: fieldname_1d_band='ssi_band'
character(len=128) :: fieldname_1d_band2='ssi_band2'
integer, parameter :: kindl = TEST_FMS_KIND_
type(FmsNetcdfFile_t) :: fileobj
real(TEST_FMS_KIND_), allocatable  :: data_in(:,:)
real(TEST_FMS_KIND_), allocatable  :: data_in_1d(:)

call mpp_init
call constants_init
call fms2_io_init
call mpp_domains_init

call create_input_files
call create_input_files_1d

call mpp_exit

contains

    !> Writes netcdf input files with fields to interpolate
    subroutine create_input_files
        if (mpp_pe() .eq. mpp_root_pe()) then
            if (open_file(fileobj, filename, "overwrite")) then
                call register_axis(fileobj, "lon", 179)
                call register_axis(fileobj, "lat", 89)
                call register_axis(fileobj, "time", unlimited)

                call register_field(fileobj, "lon", "double", dimensions=(/"lon"/))
                call register_field(fileobj, "lat", "double", dimensions=(/"lat"/))
                call register_field(fileobj, "time", "double", dimensions=(/"time"/))

                call register_field(fileobj, fieldname, "double", dimensions=(/"lon ", "lat ", "time"/))
                call register_field(fileobj, trim(fieldname)//"_random", "double", dimensions=(/"lon ","lat ","time"/))

                call register_variable_attribute(fileobj, "lon", "cartesian_axis", "X", str_len=1)
                call register_variable_attribute(fileobj, "lat", "cartesian_axis", "Y", str_len=1)

                call register_variable_attribute(fileobj, "time", "cartesian_axis", "T", str_len=1)
                call register_variable_attribute(fileobj, "time", "units", "days since 1800-01-01 00:00:00",str_len=30)
                call register_variable_attribute(fileobj, "time", "calendar", "julian", str_len=6)

                call write_data(fileobj, "lat", (/(-90.0_kindl+i*2.0_kindl,i=1,89)/))
                call write_data(fileobj, "lon", (/(-180.0_kindl+i*2.0_kindl,i=1,179)/))
                call write_data(fileobj, "time", (/(1+i*2, i=0,2)/))

                allocate(data_in(179, 89))
                do i=0, 2
                    data_in = real((1+i*2), TEST_FMS_KIND_)
                    call write_data(fileobj, fieldname, data_in, unlim_dim_level=i+1)
                    call random_number(data_in)
                    call write_data(fileobj, trim(fieldname)//"_random", data_in, unlim_dim_level=i+1)
                enddo
                call close_file(fileobj)
                deallocate(data_in)
            endif
        endif

        call mpp_sync()
    end subroutine create_input_files

    !> Writes netcdf input files with fields to interpolate
    subroutine create_input_files_1d
        if (mpp_pe() .eq. mpp_root_pe()) then
            if (open_file(fileobj, filename_1d, "overwrite")) then
                call register_axis(fileobj, "lon", 179)
                call register_axis(fileobj, "time", unlimited)

                call register_field(fileobj, "lon", "double", dimensions=(/"lon"/))
                call register_field(fileobj, "time", "double", dimensions=(/"time"/))

                call register_field(fileobj, fieldname_1d_band, "double", dimensions=(/"lon ", "time"/))
                call register_field(fileobj, trim(fieldname_1d_band2), "double", dimensions=(/"lon ","time"/))

                call register_variable_attribute(fileobj, "lon", "cartesian_axis", "X", str_len=1)

                call register_variable_attribute(fileobj, "time", "cartesian_axis", "T", str_len=1)
                call register_variable_attribute(fileobj, "time", "units", "days since 1800-01-01 00:00:00",str_len=30)
                call register_variable_attribute(fileobj, "time", "calendar", "julian", str_len=6)

                call write_data(fileobj, "lon", (/(-180.0_kindl+i*2.0_kindl,i=1,179)/))
                call write_data(fileobj, "time", (/(1+i*2, i=0,2)/))

                allocate(data_in_1d(179))
                do i=0, 2
                    data_in_1d = real((1+i*2), TEST_FMS_KIND_)
                    call write_data(fileobj, fieldname_1d_band, data_in_1d, unlim_dim_level=i+1)
                    data_in_1d = - data_in_1d
                    call write_data(fileobj, trim(fieldname_1d_band2), data_in_1d, unlim_dim_level=i+1)
                enddo
                call close_file(fileobj)
                deallocate(data_in_1d)
            endif
        endif

        call mpp_sync()
    end subroutine create_input_files_1d

end program test_time_interp_external_create_input
