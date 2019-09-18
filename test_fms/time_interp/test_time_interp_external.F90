
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

program test_time_interp_external

use constants_mod, only: constants_init
use fms_mod,       only: open_namelist_file, check_nml_error, close_file
use mpp_mod, only : mpp_init, mpp_exit, mpp_npes, stdout, stdlog, FATAL, mpp_error
use mpp_mod, only : input_nml_file
use mpp_io_mod, only : mpp_io_init, mpp_io_exit, mpp_open, MPP_RDONLY, MPP_ASCII, mpp_close, &
                       axistype, mpp_get_axis_data
use mpp_domains_mod, only : mpp_domains_init, domain2d, mpp_define_layout, mpp_define_domains,&
     mpp_global_sum, mpp_global_max, mpp_global_min, BITWISE_EXACT_SUM, mpp_get_compute_domain, &
     mpp_domains_set_stack_size
use time_interp_external_mod, only : time_interp_external, time_interp_external_init,&
     time_interp_external_exit, time_interp_external, init_external_field, get_external_field_size
use time_manager_mod, only : get_date, set_date, time_manager_init, set_calendar_type, JULIAN, time_type, increment_time,&
                             NOLEAP
use horiz_interp_mod, only: horiz_interp, horiz_interp_init, horiz_interp_new, horiz_interp_del, horiz_interp_type
use axis_utils_mod, only: get_axis_bounds
implicit none



integer :: id, i, io_status, unit, ierr
character(len=128) :: filename='INPUT/aerosol.climatology.nc', fieldname='so4_anthro'
type(time_type) :: time
real, allocatable, dimension(:,:,:) :: data_d, data_g
logical, allocatable, dimension(:,:,:) :: mask_d
type(domain2d) :: domain, domain_out
integer :: layout(2), fld_size(4)
integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
integer :: yy, mm, dd, hh, ss
real :: sm,mx,mn
character(len=12) :: cal_type
integer :: ntime=12,year0=1990,month0=1,day0=1,days_inc=31
type(horiz_interp_type) :: Hinterp
type(axistype) :: Axis_centers(4), Axis_bounds(4)
real :: lon_out(180,89), lat_out(180,89)
real, allocatable, dimension(:,:) :: lon_local_out, lat_local_out
real, allocatable, dimension(:) :: lon_in, lat_in
integer :: isc_o, iec_o, jsc_o, jec_o, outunit

namelist /test_time_interp_ext_nml/ filename, fieldname,ntime,year0,month0,&
     day0,days_inc, cal_type

call constants_init
call mpp_init
call mpp_io_init
call mpp_domains_init
call time_interp_external_init
call time_manager_init
call horiz_interp_init

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, test_time_interp_ext_nml, iostat=io_status)
      ierr = check_nml_error(io_status, 'test_time_interp_ext_nml')
#else
      unit = open_namelist_file ()
      ierr=1; do while (ierr /= 0)
      read  (unit, nml=test_time_interp_ext_nml, iostat=io_status, end=10)
      ierr = check_nml_error(io_status, 'test_time_interp_ext_nml')
      enddo
10    call close_file (unit)
#endif

outunit = stdlog()
write(outunit,test_time_interp_ext_nml)

select case (trim(cal_type))
case ('julian')
   call set_calendar_type(JULIAN)
case ('no_leap')
   call set_calendar_type(NOLEAP)
case default
   call mpp_error(FATAL,'invalid calendar type')
end select

outunit = stdout()
write(outunit,*) 'INTERPOLATING NON DECOMPOSED FIELDS'
write(outunit,*) '======================================'

call time_interp_external_init

id = init_external_field(filename,fieldname,verbose=.true.)

fld_size = get_external_field_size(id)

allocate(data_g(fld_size(1),fld_size(2),fld_size(3)))
data_g = 0

time = set_date(year0,month0,day0,0,0,0)

do i=1,ntime
   call time_interp_external(id,time,data_g,verbose=.true.)
   sm = sum(data_g)
   mn = minval(data_g)
   mx = maxval(data_g)
   write(outunit,*) 'sum= ', sm
   write(outunit,*) 'max= ', mx
   write(outunit,*) 'min= ', mn
   time = increment_time(time,0,days_inc)
enddo

call mpp_define_layout((/1,fld_size(1),1,fld_size(2)/),mpp_npes(),layout)
call mpp_define_domains((/1,fld_size(1),1,fld_size(2)/),layout,domain)
call mpp_get_compute_domain(domain,isc,iec,jsc,jec)
call mpp_get_compute_domain(domain,isd,ied,jsd,jed)

call mpp_domains_set_stack_size(fld_size(1)*fld_size(2)*min(fld_size(3),1)*2)
allocate(data_d(isd:ied,jsd:jed,fld_size(3)))
data_d = 0

write(outunit,*) 'INTERPOLATING DOMAIN DECOMPOSED FIELDS'
write(outunit,*) '======================================'

id = init_external_field(filename,fieldname,domain=domain, verbose=.true.)

time = set_date(year0,month0,day0)

do i=1,ntime
   call time_interp_external(id,time,data_d,verbose=.true.)
   sm = mpp_global_sum(domain,data_d,flags=BITWISE_EXACT_SUM)
   mx = mpp_global_max(domain,data_d)
   mn = mpp_global_min(domain,data_d)
   write(outunit,*) 'global sum= ', sm
   write(outunit,*) 'global max= ', mx
   write(outunit,*) 'global min= ', mn
   time = increment_time(time,0,days_inc)
enddo

write(outunit,*) 'INTERPOLATING DOMAIN DECOMPOSED FIELDS USING HORIZ INTERP'
write(outunit,*) '======================================'


! define a global 2 degree output grid

do i=1,180
   lon_out(i,:) = 2.0*i*atan(1.0)/45.0
enddo

do i=1,89
   lat_out(:,i) = (i-45)*2.0*atan(1.0)/45.0
enddo

call mpp_define_layout((/1,180,1,89/),mpp_npes(),layout)
call mpp_define_domains((/1,180,1,89/),layout,domain_out)
call mpp_get_compute_domain(domain_out,isc_o,iec_o,jsc_o,jec_o)

id = init_external_field(filename,fieldname,domain=domain_out,axis_centers=axis_centers,&
      verbose=.true., override=.true.)

allocate (lon_local_out(isc_o:iec_o,jsc_o:jec_o))
allocate (lat_local_out(isc_o:iec_o,jsc_o:jec_o))

lon_local_out(isc_o:iec_o,jsc_o:jec_o) = lon_out(isc_o:iec_o,jsc_o:jec_o)
lat_local_out(isc_o:iec_o,jsc_o:jec_o) = lat_out(isc_o:iec_o,jsc_o:jec_o)

call get_axis_bounds(axis_centers(1), axis_bounds(1), axis_centers)
call get_axis_bounds(axis_centers(2), axis_bounds(2), axis_centers)

allocate(lon_in(fld_size(1)+1))
allocate(lat_in(fld_size(2)+1))

call mpp_get_axis_data(axis_bounds(1), lon_in) ; lon_in = lon_in*atan(1.0)/45
call mpp_get_axis_data(axis_bounds(2), lat_in) ; lat_in = lat_in*atan(1.0)/45

call horiz_interp_new(Hinterp,lon_in,lat_in, lon_local_out, lat_local_out, &
     interp_method='bilinear')

time = set_date(year0,month0,day0)

deallocate(data_d)
allocate(data_d(isc_o:iec_o,jsc_o:jec_o,fld_size(3)))
allocate(mask_d(isc_o:iec_o,jsc_o:jec_o,fld_size(3)))
do i=1,ntime
   data_d = 0
   call time_interp_external(id,time,data_d,verbose=.true.,horz_interp=Hinterp, mask_out=mask_d)
   sm = mpp_global_sum(domain_out,data_d,flags=BITWISE_EXACT_SUM)
   mx = mpp_global_max(domain_out,data_d)
   mn = mpp_global_min(domain_out,data_d)
   write(outunit,*) 'global sum= ', sm
   write(outunit,*) 'global max= ', mx
   write(outunit,*) 'global min= ', mn

   where(mask_d)
      data_d = 1.0
   elsewhere
      data_d = 0.0
   endwhere
   sm = mpp_global_sum(domain_out,data_d,flags=BITWISE_EXACT_SUM)
   write(outunit,*) 'n valid points= ', sm

   time = increment_time(time,0,days_inc)
enddo

call horiz_interp_del(Hinterp)


call time_interp_external_exit


call mpp_io_exit
call mpp_exit
stop

end program test_time_interp_external
