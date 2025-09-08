
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

program test_time_interp_external

use constants_mod,             only : constants_init
use fms_mod,                   only : check_nml_error
use mpp_mod,                   only : mpp_init, mpp_exit, mpp_npes, stdout, stdlog, FATAL, mpp_error
use mpp_mod,                   only : input_nml_file, mpp_pe, mpp_root_pe, mpp_sync
use mpp_domains_mod,           only : mpp_domains_init, domain2d, mpp_define_layout, mpp_define_domains
use mpp_domains_mod, only : mpp_global_sum, mpp_global_max, mpp_global_min, BITWISE_EXACT_SUM, mpp_get_compute_domain
use mpp_domains_mod,           only : mpp_domains_set_stack_size
use time_interp_external2_mod, only : time_interp_external, time_interp_external_init, get_external_fileobj
use time_interp_external2_mod, only : time_interp_external_exit, time_interp_external, init_external_field, &
                                   &  get_external_field_size
use time_manager_mod,          only : get_date, set_date, time_manager_init, set_calendar_type, JULIAN, &
                                   &  time_type, increment_time
use time_manager_mod,          only : NOLEAP
use horiz_interp_mod,          only : horiz_interp, horiz_interp_init, horiz_interp_new, horiz_interp_del, &
                                    & horiz_interp_type
use axis_utils2_mod,           only : axis_edges
use fms2_io_mod,               only : FmsNetcdfFile_t, fms2_io_init, open_file, close_file, write_data, register_axis
use fms2_io_mod,               only : register_field, unlimited, register_variable_attribute
use platform_mod

implicit none

integer            :: id !< Time_interp_external id
integer            :: i  !< Index for loops
integer            :: ierr !< io return status
character(len=128) :: filename='INPUT/aerosol.climatology.nc'
character(len=128) :: filename_1d='INPUT/solar_constant.nc'
character(len=128) :: fieldname='so4_anthro'
character(len=128) :: fieldname_1d_band='ssi_band'
character(len=128) :: fieldname_1d_band2='ssi_band2'
type(time_type)    :: time !< "model" time
integer, parameter :: kindl = TEST_FMS_KIND_
real(TEST_FMS_KIND_) :: data_d_0d = 1.0_kindl !< interpolated data in compute domain
real(TEST_FMS_KIND_) :: data_g_0d = 1.0_kindl !< interpolated global data
type(domain2d)     :: domain !<Domain to interpolate data to
integer            :: layout(2) !< Domain layout
integer            :: fld_size(4) !< Size of fieldname
integer            :: isc, iec, jsc, jec !< starting (s)/ending (e) x/y indexes of the compute domain
integer            :: isd, ied, jsd, jed !< starting (s)/ending (e) x/y indexes of the data domain
real(TEST_FMS_KIND_)               :: sm !< Sum of a data array
real(TEST_FMS_KIND_)               :: mx !< Max value of a data array
real(TEST_FMS_KIND_)               :: mn !< Min value of a data_array
character(len=12)  :: cal_type="julian" !< Calendar type
type(horiz_interp_type) :: Hinterp !< horix interp type
real(TEST_FMS_KIND_)               :: lon_out(180,89) !< lat grid to interpolate to (global)
real(TEST_FMS_KIND_)               :: lat_out(180,89) !< lon grid to interpolate to (global)
integer            :: outunit !< stdout unit number
type(FmsNetcdfFile_t) :: fileobj !< fileobj
character(len=12)  :: axis_names(4) = (/ "lon ", "lat ", "time", "none"/) !< axis_names
real(TEST_FMS_KIND_), allocatable  :: data_in(:,:)  !< data added to file
real(TEST_FMS_KIND_), allocatable  :: data_in_1d(:)  !< data added to file

namelist /test_time_interp_external_nml/ cal_type

call mpp_init
call constants_init
call fms2_io_init
call mpp_domains_init
call time_interp_external_init
call time_manager_init
call horiz_interp_init

read (input_nml_file, test_time_interp_external_nml, iostat=ierr)
ierr = check_nml_error (ierr, 'test_time_interp_external_nml')

call create_input_files
call create_input_files_1d

select case (trim(cal_type))
case ('julian')
   call set_calendar_type(JULIAN)
case ('no_leap')
   call set_calendar_type(NOLEAP)
case default
   call mpp_error(FATAL,'invalid calendar type')
end select

call time_interpolate_3d_data
call time_interpolate_2d_data
call time_interpolate_1d_data
call time_interpolate_scalar_data

call time_interp_external_exit
call mpp_exit


 contains

    !> Writes netcdf input files with fields to interpolate
    subroutine create_input_files
        !< Create a file to test with:
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

                !call register_global_attribute(fileobj, "global_scalar", 1024.0_r8_kind)

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
        !< Wait for the root pe to catch up
        call mpp_sync()
    end subroutine


    !> Writes netcdf input files with fields to interpolate
    subroutine create_input_files_1d
        !< Create a file to test with:
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

                !call register_global_attribute(fileobj, "global_scalar", 1024.0_r8_kind)

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
        !< Wait for the root pe to catch up
        call mpp_sync()
    end subroutine create_input_files_1d

    subroutine time_interpolate_3d_data
        real(TEST_FMS_KIND_), allocatable  :: data_d_3d(:,:,:) !< interpolated data in compute domain
        real(TEST_FMS_KIND_), allocatable  :: data_g_3d(:,:,:) !< interpolated global data
        real(TEST_FMS_KIND_), allocatable  :: lon_local_out(:,:) !< lat grid to interpolate to (compute)
        real(TEST_FMS_KIND_), allocatable  :: lat_local_out(:,:) !< lon grid to interpolate to (compute)
        real(TEST_FMS_KIND_), allocatable  :: lon_in(:) !< lat grid in file
        real(TEST_FMS_KIND_), allocatable  :: lat_in(:) !< lat grid in file

        outunit = stdout()
        write(outunit,*) 'INTERPOLATING NON DECOMPOSED FIELDS'
        write(outunit,*) '======================================'
        !< Here every rank is on the same x/y grid

        id = init_external_field(filename,fieldname,verbose=.false.)
        fld_size = get_external_field_size(id)
        allocate(data_g_3d(fld_size(1),fld_size(2),fld_size(3)))
        data_g_3d = 0.0_kindl

        !< Corresponds to time=1 in the file, so no time_interpolation
        time = set_date(1800,1,2,0,0,0)
        call time_interp_external(id,time,data_g_3d,verbose=.true.)
        sm = sum(data_g_3d)
        mn = minval(data_g_3d)
        mx = maxval(data_g_3d)
        write(outunit,*) 'sum= ', sm
        write(outunit,*) 'max= ', mx
        write(outunit,*) 'min= ', mn

        !< Corresponds to time=2, so there will be time_interpolation between time=1 and
        !time =3 in the file
        time = set_date(1800,1,3,0,0,0)
        call time_interp_external(id,time,data_g_3d,verbose=.true.)
        sm = sum(data_g_3d)
        mn = minval(data_g_3d)
        mx = maxval(data_g_3d)
        write(outunit,*) 'sum= ', sm
        write(outunit,*) 'max= ', mx
        write(outunit,*) 'min= ', mn

        write(outunit,*) 'INTERPOLATING DOMAIN DECOMPOSED FIELDS'
        write(outunit,*) '======================================'
        !< Here every rank has its own section

        call mpp_define_layout((/1,fld_size(1),1,fld_size(2)/),mpp_npes(),layout)
        call mpp_define_domains((/1,fld_size(1),1,fld_size(2)/),layout,domain)
        call mpp_get_compute_domain(domain,isc,iec,jsc,jec)
        call mpp_get_compute_domain(domain,isd,ied,jsd,jed)

        call mpp_domains_set_stack_size(fld_size(1)*fld_size(2)*min(fld_size(3),1)*2)
        allocate(data_d_3d(isd:ied,jsd:jed,fld_size(3)))
        data_d_3d = 0

        id = init_external_field(filename,fieldname,domain=domain, verbose=.false.)

        time = set_date(1800,1,2,0,0,0)
        call time_interp_external(id,time,data_d_3d,verbose=.true.)
        sm = mpp_global_sum(domain,data_d_3d,flags=BITWISE_EXACT_SUM)
        mx = mpp_global_max(domain,data_d_3d)
        mn = mpp_global_min(domain,data_d_3d)
        write(outunit,*) 'global sum= ', sm
        write(outunit,*) 'global max= ', mx
        write(outunit,*) 'global min= ', mn

        !< Corresponds to time=2, so there will be time_interpolation between time=1 and
        !time =3 in the file
        time = set_date(1800,1,3,0,0,0)
        call time_interp_external(id,time,data_d_3d,verbose=.true.)
        sm = mpp_global_sum(domain,data_d_3d,flags=BITWISE_EXACT_SUM)
        mx = mpp_global_max(domain,data_d_3d)
        mn = mpp_global_min(domain,data_d_3d)
        write(outunit,*) 'global sum= ', sm
        write(outunit,*) 'global max= ', mx
        write(outunit,*) 'global min= ', mn

        write(outunit,*) 'INTERPOLATING DOMAIN DECOMPOSED FIELDS USING HORIZ INTERP'
        write(outunit,*) '======================================'
        !< Here there will be interpolation in x/y/time

        ! define a global 2 degree output grid
        do i=1,180
        lon_out(i,:) = 2.0_kindl*real(i,kindl)*atan(1.0_kindl)/45.0_kindl
        enddo

        do i=1,89
        lat_out(:,i) = (i-45)*2.0_kindl*atan(1.0_kindl)/45.0_kindl
        enddo

        id = init_external_field(filename,trim(fieldname)//"_random",domain=domain,axis_names=axis_names,&
            verbose=.false., override=.true.)

        if (.not. get_external_fileobj(filename, fileobj)) call mpp_error(FATAL, "Couldn't find the fileobj for: "// &
            & trim(filename))

        allocate (lon_local_out(isc:iec,jsc:jec))
        allocate (lat_local_out(isc:iec,jsc:jec))

        lon_local_out(isc:iec,jsc:jec) = lon_out(isc:iec,jsc:jec)
        lat_local_out(isc:iec,jsc:jec) = lat_out(isc:iec,jsc:jec)

        !< Read in the axis data
        allocate(lon_in(fld_size(1)+1))
        allocate(lat_in(fld_size(2)+1))

        call axis_edges(fileobj, axis_names(1), lon_in)
        call axis_edges(fileobj, axis_names(2), lat_in)

        lon_in = lon_in*atan(1.0_kindl)/45.0_kindl
        lat_in = lat_in*atan(1.0_kindl)/45.0_kindl

        call horiz_interp_new(Hinterp,lon_in,lat_in, lon_local_out, lat_local_out, &
            interp_method='bilinear')

        deallocate(data_d_3d)
        allocate(data_d_3d(isc:iec,jsc:jec,fld_size(3)))

        time = set_date(1800,1,3,0,0,0)
        data_d_3d = 0.0_kindl
        call time_interp_external(id,time,data_d_3d,verbose=.true.,horz_interp=Hinterp)
        sm = mpp_global_sum(domain,data_d_3d,flags=BITWISE_EXACT_SUM)
        mx = mpp_global_max(domain,data_d_3d)
        mn = mpp_global_min(domain,data_d_3d)
        write(outunit,*) 'global sum= ', sm
        write(outunit,*) 'global max= ', mx
        write(outunit,*) 'global min= ', mn

        sm = mpp_global_sum(domain,data_d_3d,flags=BITWISE_EXACT_SUM)
        write(outunit,*) 'n valid points= ', sm

        call horiz_interp_del(Hinterp)

        deallocate(data_d_3d, data_g_3d, lon_local_out, lat_local_out, lon_in, lat_in)
    end subroutine time_interpolate_3d_data

    subroutine time_interpolate_2d_data
        real(TEST_FMS_KIND_), allocatable  :: data_d_2d(:,:) !< interpolated data in compute domain
        real(TEST_FMS_KIND_), allocatable  :: data_g_2d(:,:) !< interpolated global data
        real(TEST_FMS_KIND_), allocatable  :: lon_local_out(:,:) !< lat grid to interpolate to (compute)
        real(TEST_FMS_KIND_), allocatable  :: lat_local_out(:,:) !< lon grid to interpolate to (compute)
        real(TEST_FMS_KIND_), allocatable  :: lon_in(:) !< lat grid in file
        real(TEST_FMS_KIND_), allocatable  :: lat_in(:) !< lat grid in file

        outunit = stdout()
        write(outunit,*) 'INTERPOLATING NON DECOMPOSED FIELDS'
        write(outunit,*) '======================================'
        !< Here every rank is on the same x/y grid

        id = init_external_field(filename,fieldname,verbose=.false.)
        fld_size = get_external_field_size(id)
        allocate(data_g_2d(fld_size(1),fld_size(2)))
        data_g_2d = 0.0_kindl

        !< Corresponds to time=1 in the file, so no time_interpolation
        time = set_date(1800,1,2,0,0,0)
        call time_interp_external(id,time,data_g_2d,verbose=.true.)
        sm = sum(data_g_2d)
        mn = minval(data_g_2d)
        mx = maxval(data_g_2d)
        write(outunit,*) 'sum= ', sm
        write(outunit,*) 'max= ', mx
        write(outunit,*) 'min= ', mn

        !< Corresponds to time=2, so there will be time_interpolation between time=1 and
        !time =3 in the file
        time = set_date(1800,1,3,0,0,0)
        call time_interp_external(id,time,data_g_2d,verbose=.true.)
        sm = sum(data_g_2d)
        mn = minval(data_g_2d)
        mx = maxval(data_g_2d)
        write(outunit,*) 'sum= ', sm
        write(outunit,*) 'max= ', mx
        write(outunit,*) 'min= ', mn

        write(outunit,*) 'INTERPOLATING DOMAIN DECOMPOSED FIELDS'
        write(outunit,*) '======================================'
        !< Here every rank has its own section

        call mpp_define_layout((/1,fld_size(1),1,fld_size(2)/),mpp_npes(),layout)
        call mpp_define_domains((/1,fld_size(1),1,fld_size(2)/),layout,domain)
        call mpp_get_compute_domain(domain,isc,iec,jsc,jec)
        call mpp_get_compute_domain(domain,isd,ied,jsd,jed)

        call mpp_domains_set_stack_size(fld_size(1)*fld_size(2)*min(fld_size(3),1)*2)
        allocate(data_d_2d(isd:ied,jsd:jed))
        data_d_2d = 0

        id = init_external_field(filename,fieldname,domain=domain, verbose=.false.)

        time = set_date(1800,1,2,0,0,0)
        call time_interp_external(id,time,data_d_2d,verbose=.true.)
        sm = mpp_global_sum(domain,data_d_2d,flags=BITWISE_EXACT_SUM)
        mx = mpp_global_max(domain,data_d_2d)
        mn = mpp_global_min(domain,data_d_2d)
        write(outunit,*) 'global sum= ', sm
        write(outunit,*) 'global max= ', mx
        write(outunit,*) 'global min= ', mn

        !< Corresponds to time=2, so there will be time_interpolation between time=1 and
        !time =3 in the file
        time = set_date(1800,1,3,0,0,0)
        call time_interp_external(id,time,data_d_2d,verbose=.true.)
        sm = mpp_global_sum(domain,data_d_2d,flags=BITWISE_EXACT_SUM)
        mx = mpp_global_max(domain,data_d_2d)
        mn = mpp_global_min(domain,data_d_2d)
        write(outunit,*) 'global sum= ', sm
        write(outunit,*) 'global max= ', mx
        write(outunit,*) 'global min= ', mn

        write(outunit,*) 'INTERPOLATING DOMAIN DECOMPOSED FIELDS USING HORIZ INTERP'
        write(outunit,*) '======================================'
        !< Here there will be interpolation in x/y/time

        ! define a global 2 degree output grid
        do i=1,180
        lon_out(i,:) = 2.0_kindl*real(i,kindl)*atan(1.0_kindl)/45.0_kindl
        enddo

        do i=1,89
        lat_out(:,i) = (i-45)*2.0_kindl*atan(1.0_kindl)/45.0_kindl
        enddo

        id = init_external_field(filename,trim(fieldname)//"_random",domain=domain,axis_names=axis_names,&
            verbose=.false., override=.true.)

        if (.not. get_external_fileobj(filename, fileobj)) call mpp_error(FATAL, "Couldn't find the fileobj for: "// &
            & trim(filename))

        allocate (lon_local_out(isc:iec,jsc:jec))
        allocate (lat_local_out(isc:iec,jsc:jec))

        lon_local_out(isc:iec,jsc:jec) = lon_out(isc:iec,jsc:jec)
        lat_local_out(isc:iec,jsc:jec) = lat_out(isc:iec,jsc:jec)

        !< Read in the axis data
        allocate(lon_in(fld_size(1)+1))
        allocate(lat_in(fld_size(2)+1))

        call axis_edges(fileobj, axis_names(1), lon_in)
        call axis_edges(fileobj, axis_names(2), lat_in)

        lon_in = lon_in*atan(1.0_kindl)/45.0_kindl
        lat_in = lat_in*atan(1.0_kindl)/45.0_kindl

        call horiz_interp_new(Hinterp,lon_in,lat_in, lon_local_out, lat_local_out, &
            interp_method='bilinear')

        deallocate(data_d_2d)
        allocate(data_d_2d(isc:iec,jsc:jec))

        time = set_date(1800,1,3,0,0,0)
        data_d_2d = 0.0_kindl
        call time_interp_external(id,time,data_d_2d,verbose=.true.,horz_interp=Hinterp)
        sm = mpp_global_sum(domain,data_d_2d,flags=BITWISE_EXACT_SUM)
        mx = mpp_global_max(domain,data_d_2d)
        mn = mpp_global_min(domain,data_d_2d)
        write(outunit,*) 'global sum= ', sm
        write(outunit,*) 'global max= ', mx
        write(outunit,*) 'global min= ', mn

        sm = mpp_global_sum(domain,data_d_2d,flags=BITWISE_EXACT_SUM)
        write(outunit,*) 'n valid points= ', sm

        call horiz_interp_del(Hinterp)

        deallocate(data_d_2d, data_g_2d, lon_local_out, lat_local_out, lon_in, lat_in)
    end subroutine time_interpolate_2d_data

    subroutine time_interpolate_1d_data
        real(TEST_FMS_KIND_), allocatable  :: data_d_1d(:) !< interpolated data in compute domain
        real(TEST_FMS_KIND_), allocatable  :: data_g_1d(:) !< interpolated global data
        real(TEST_FMS_KIND_), allocatable  :: lon_local_out(:) !< lat grid to interpolate to (compute)
        real(TEST_FMS_KIND_), allocatable  :: lat_local_out(:) !< lon grid to interpolate to (compute)
        real(TEST_FMS_KIND_), allocatable  :: lon_in(:) !< lat grid in file
        real(TEST_FMS_KIND_), allocatable  :: lat_in(:) !< lat grid in file

        outunit = stdout()
        write(outunit,*) 'INTERPOLATING NON DECOMPOSED FIELDS 1D'
        write(outunit,*) '======================================'
        !< Here every rank is on the same x/y grid

        id = init_external_field(filename_1d,fieldname_1d_band,verbose=.false.)
        fld_size = get_external_field_size(id)

        ! check field size
        write(*,*) 'FIELD SIZE', fld_size

        if(fld_size(1) /= 179) call mpp_error(FATAL, "incorrect X field size for time_interp_external_1d")
        if(fld_size(2) /= 1) call mpp_error(FATAL, "incorrect Y field size for time_interp_external_1d")
        if(fld_size(3) /= 1) call mpp_error(FATAL, "incorrect Z field size for time_interp_external_1d")
        if(fld_size(4) /= 3) call mpp_error(FATAL, "incorrect T field size for time_interp_external_1d")

        allocate(data_g_1d(fld_size(1)))
        data_g_1d = 0.0_kindl

        !< Corresponds to time=1 in the file, so no time_interpolation
        time = set_date(1800,1,2,0,0,0)
        call time_interp_external(id,time,data_g_1d,verbose=.true.)
        !check answers
        do i=1, fld_size(1)
          if(data_g_1d(i) /= 1.0) call mpp_error(FATAL, "incorrect interpolation for t=1 in time_interp_external_1d")
        end do

        !< Corresponds to time=2, so there will be time_interpolation between time=1 and
        !time =3 in the file
        time = set_date(1800,1,3,0,0,0)
        call time_interp_external(id,time,data_g_1d,verbose=.true.)
        !check answers
        do i=1, fld_size(1)
          if(data_g_1d(i) /= 2.0) call mpp_error(FATAL, "incorrect interpolation for t=3 in time_interp_external_1d")
        end do

        !< Here every rank has its own section

        call mpp_define_layout((/1,fld_size(1),1,fld_size(2)/),mpp_npes(),layout)
        call mpp_define_domains((/1,fld_size(1),1,fld_size(2)/),layout,domain)
        call mpp_get_compute_domain(domain,isc,iec,jsc,jec)
        call mpp_get_compute_domain(domain,isd,ied,jsd,jed)

        call mpp_domains_set_stack_size(fld_size(1)*fld_size(2)*min(fld_size(3),1)*2)
        allocate(data_d_1d(isd:ied))
        data_d_1d = 0

        id = init_external_field(filename_1d,fieldname_1d_band2,domain=domain, verbose=.false.)

        time = set_date(1800,1,2,0,0,0)
        call time_interp_external(id,time,data_d_1d,verbose=.true.)
        ! answer check
        do i=isd, ied
          if(data_d_1d(i) /= -1.0) call mpp_error(FATAL, "incorrect interpolation for t=1 in time_interp_external_1d")
        end do

        !< Corresponds to time=2, so there will be time_interpolation between time=1 and
        !time =3 in the file
        time = set_date(1800,1,3,0,0,0)
        call time_interp_external(id,time,data_d_1d,verbose=.true.)
        ! answer check
        do i=isd, ied
          if(data_d_1d(i) /= -2.0) call mpp_error(FATAL, "incorrect interpolation for t=1 in time_interp_external_1d")
        end do

        deallocate(data_d_1d, data_g_1d)
    end subroutine time_interpolate_1d_data

    subroutine time_interpolate_scalar_data
        real(TEST_FMS_KIND_) :: data_scalar !< scalar to interpolate

        id = init_external_field(filename, fieldname,verbose=.false.)

        !< Corresponds to time=1 in the file, so no time_interpolation
        time = set_date(1800,1,2,0,0,0)
        call time_interp_external(id,time,data_scalar,verbose=.true.)
        print *, 'data= ', data_scalar

        !< Corresponds to midpoint, so interpolation between 1 and 2 periods
        time = set_date(1800,1,2,12,0,0)
        call time_interp_external(id,time,data_scalar,verbose=.true.)
        print *, 'data= ', data_scalar

        !< Corresponds to time=2, so there will be time_interpolation between time=1 and
        !time =3 in the file
        time = set_date(1800,1,3,0,0,0)
        call time_interp_external(id,time,data_scalar,verbose=.true.)
        print *, 'data= ', data_scalar

    end subroutine time_interpolate_scalar_data


end program test_time_interp_external
