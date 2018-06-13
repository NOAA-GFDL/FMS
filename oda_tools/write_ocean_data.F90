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
module write_ocean_data_mod

 use mpp_io_mod, only : fieldtype, axistype, mpp_open,&
      MPP_OVERWR, MPP_NETCDF, MPP_MULTI, MPP_SINGLE,&
      mpp_write_meta, mpp_write, mpp_close
 use mpp_mod, only : mpp_error, FATAL
 use oda_types_mod, only : missing_value
 use oda_types_mod, only : ocean_profile_type, max_levels_file
 use time_manager_mod, only : time_type, get_time, set_date, operator ( - )

 implicit none

 private

 type(fieldtype), save :: lon_field, lat_field, time_field, data_t_field, data_s_field, &
      project_field,probe_field,ref_inst_field, fix_depth_field, database_id_field,&
      profile_flag_field, profile_flag_s_field, temp_err_field, salt_err_field, &
      flag_t_field, flag_s_field, ocn_vehicle_field,&
      depth_field, nvar_field, lon_index_field, lat_index_field, &
      yyyy_field, mmdd_field, link_field

 integer, parameter :: ref_yr=1900, ref_mon=1, ref_day=1,&
                       ref_hr=0, ref_min=0, ref_sec=0,max_files=1000

 integer :: ref_seconds, ref_days, chid, wmo_id

 integer,save :: nvar_out

 integer, save :: sta_num(max_files), unit_num(max_files), nfiles

 type(time_type) :: ref_time, time

 logical :: module_is_initialized=.false.

 public :: open_profile_file, write_profile, close_profile_file, &
      write_ocean_data_init

#include <netcdf.inc>

contains

function open_profile_file(name, nvar, grid_lon, grid_lat,thread,fset)

  character(len=*), intent(in) :: name
  integer, intent(in), optional :: nvar
  real, dimension(:), optional, intent(in) :: grid_lon, grid_lat
  integer, intent(in), optional :: thread, fset

  integer :: i, open_profile_file, unit
  integer :: threading, fileset
  character(len=128) :: units, time_units
  real, dimension(max_levels_file) :: array

type(axistype) :: depth_axis, station_axis, lon_axis, lat_axis

threading=MPP_MULTI
fileset=MPP_SINGLE

if (PRESENT(thread)) threading=thread
if (PRESENT(fset)) fileset=fset

ref_time = set_date(ref_yr, ref_mon, ref_day, ref_hr, ref_min, ref_sec)
call get_time(ref_time, ref_seconds, ref_days)
call mpp_open(unit, trim(name), action=MPP_OVERWR, form=MPP_NETCDF,&
              threading=threading, fileset=fileset)

open_profile_file = unit

nfiles=nfiles+1

if (nfiles > max_files) call mpp_error(FATAL,'max number of profiles exceeded&
     &in module write_ocean_data, increase param : max_files')

unit_num(nfiles) = unit


nvar_out = 2
if (PRESENT(nvar)) nvar_out = nvar

if (PRESENT(grid_lon) .and. PRESENT(grid_lat)) then
   call mpp_write_meta(unit, lon_axis, 'grid_longitude','degrees_E',&
        'observational grid longitude',cartesian='X',sense=1,data=grid_lon)

   call mpp_write_meta(unit, lat_axis, 'grid_latitude','degrees_N',&
        'observational grid latitude', cartesian='Y',sense=1,data=grid_lat)
endif

!call mpp_write_meta(unit,depth_axis,'depth_index','none','depth index',&
!                  cartesian='Z',sense=-1)!,data=(/(float(i),i=1,max_levels_file)/))
!pgf90 complains about the above. This is a compiler bug. Workaround:
array = (/(float(i),i=1,max_levels_file)/)
call mpp_write_meta(unit,depth_axis,'depth_index','none','depth index',&
                    cartesian='Z',sense=-1,data=array)

call mpp_write_meta(unit,station_axis,'station_index','none',&
                    'station index', cartesian='T',sense=1)

if (PRESENT(grid_lon) .and. PRESENT(grid_lat)) then
   call mpp_write_meta(unit, lon_index_field, (/station_axis/),&
        'longitude_index','none','longitude_index', missing=missing_value)
   call mpp_write_meta(unit, lat_index_field, (/station_axis/),&
        'latitude_index','none','latitude_index',missing=missing_value)
endif

call mpp_write_meta(unit,nvar_field,(/station_axis/),&
     'nvar','none','temp (1) or temp and salt (2)')

call mpp_write_meta(unit,lon_field,(/station_axis/),&
                   'longitude','degrees_E','longitude',&
                    min=-1.0,max=361.0)

call mpp_write_meta(unit,lat_field,(/station_axis/),&
                   'latitude','degrees_N','latitude',&
                    min=-91.0,max=91.0)

call mpp_write_meta(unit,profile_flag_field,(/station_axis/),&
                   'profile_flag','none','profile_flag',&
                   min=0.0,max=10.0,missing=missing_value)


if (nvar_out .eq. 2) call mpp_write_meta(unit,profile_flag_s_field,(/station_axis/),&
                   'profile_flag_s','none','profile_flag for salt',&
                    min=0.0,max=10.0,missing=missing_value)


write(time_units,'(a,i4.4,a,i2.2,a,i2.2,a)')  'days since ',ref_yr,'-',ref_mon,'-',ref_day,' 00:00:00'

call mpp_write_meta(unit,time_field,(/station_axis/),&
                   'time',trim(time_units),'time')

call mpp_write_meta(unit,yyyy_field,(/station_axis/),&
     'yyyy','none','yyyy')

call mpp_write_meta(unit,mmdd_field,(/station_axis/),&
                   'mmdd','none','mmdd')



units='deg_C'
call mpp_write_meta(unit,temp_err_field,(/station_axis/),&
                   'temp_error',trim(units),'measurement error of temperature',missing=missing_value)

units='g/kg'
if (nvar_out .eq. 2) call mpp_write_meta(unit,salt_err_field,(/station_axis/),&
                   'salt_error',trim(units),'measurement error of salinity',missing=missing_value)

call mpp_write_meta(unit,project_field,(/station_axis/),&
     'project','none','see NODC codes')

call mpp_write_meta(unit,probe_field,(/station_axis/),&
     'probe','none','see NODC codes')

call mpp_write_meta(unit,ref_inst_field,(/station_axis/),&
     'ref_inst','none','see NODC codes')

call mpp_write_meta(unit,fix_depth_field,(/station_axis/),&
     'fix_depth','none','see NODC codes')

call mpp_write_meta(unit,database_id_field,(/station_axis/),&
     'database_id','none','see NODC codes')

call mpp_write_meta(unit,ocn_vehicle_field,(/station_axis/),&
     'ocn_vehicle','none','see NODC codes')

call mpp_write_meta(unit,link_field,(/station_axis/),&
     'link','none','partial_profile flag')


 units='degrees_C'
 call mpp_write_meta(unit,data_t_field,(/depth_axis,station_axis/),&
              'temp',trim(units),'in-situ temperature',&
                min=-10.0,max=50.0,missing=missing_value)

 units='g/kg'
 if (nvar_out .eq. 2) call mpp_write_meta(unit,data_s_field,(/depth_axis,station_axis/),&
                   'salt',trim(units),'salinity',&
                    min=0.0,max=50.0,missing=missing_value)

call mpp_write_meta(unit,depth_field,(/depth_axis,station_axis/),&
                   'depth','meters','depth of obs',&
                    min=0.0,max=7000.0,missing=missing_value)



call mpp_write_meta(unit,flag_t_field,(/depth_axis,station_axis/),&
     'temp_flag','none','temperature level flag (see NODC codes)',missing=missing_value)

if (nvar_out .eq. 2) call mpp_write_meta(unit,flag_s_field,(/depth_axis,station_axis/),&
     'salt_flag','none','salinity level flag (see NODC codes)',missing=missing_value)



call mpp_write(unit, depth_axis)

if (PRESENT(grid_lon).and.PRESENT(grid_lat)) then
   call mpp_write(unit, lon_axis)
   call mpp_write(unit, lat_axis)
endif

end function open_profile_file


subroutine write_profile(unit,profile)

use mpp_domains_mod, only : domain2d,mpp_get_compute_domain, &
                            mpp_get_data_domain
use mpp_mod, only : mpp_pe

integer, intent(in) :: unit
type(ocean_profile_type), intent(in) :: profile

real, dimension(max_levels_file) :: data_t, data_s, depth
integer :: levels, secs, days, i, j, nlinks
real :: profile_flag, profile_flag_s, days_since, error, nvar, station
real :: tmp_s
real, dimension(max_levels_file) :: flag_t, flag_s
logical :: grid_ptr = .false.
integer :: findex
integer :: isc,iec,jsc,jec,isd,ied,jsd,jed
logical :: debug=.false.

! find file index from file unit list

findex=-1
do i=1,nfiles
   if (unit_num(i) .eq. unit) then
       findex=i
       exit
   endif
enddo

if (findex .eq. -1) call mpp_error(FATAL,'Attempt write to unopened file in&
     &write_ocean_data_mod:write_profile_data')


sta_num(findex)=sta_num(findex)+1

station=sta_num(findex)

levels = min(profile%levels,max_levels_file)
data_t=missing_value;data_s=missing_value;depth=missing_value
flag_t=missing_value;flag_s=missing_value
data_t(1:levels)=profile%data_t(1:levels)
flag_t(1:levels)=profile%flag_t(1:levels)

if (ASSOCIATED(profile%Model_Grid)) grid_ptr = .true.


if (grid_ptr) then
    call mpp_get_compute_domain(profile%Model_Grid%Dom, isc, iec, jsc, jec)
    if (floor(profile%i_index) .lt. isc .or. floor(profile%i_index) .gt. iec) return
    if (floor(profile%j_index) .lt. jsc .or. floor(profile%j_index) .gt. jec) return
endif

if (profile%nvar == 2) then
    data_s(1:levels)   = profile%data_s(1:levels)
    flag_s(1:levels)=profile%flag_s(1:levels)
endif

depth(1:levels)=profile%depth(1:levels)
time = profile%time - ref_time
call get_time(time, secs, days)
days_since = days + secs/86400.


 nvar = profile%nvar
 call mpp_write(unit,nvar_field,nvar,station)
 call mpp_write(unit,data_t_field,data_t,station)
 if (nvar_out .eq. 2) call mpp_write(unit,data_s_field,data_s,station)
 call mpp_write(unit,depth_field,depth,station)
 call mpp_write(unit,project_field,profile%project,station)
 call mpp_write(unit,probe_field,profile%probe,station)
 call mpp_write(unit,ref_inst_field,profile%ref_inst,station)
 call mpp_write(unit,fix_depth_field,profile%fix_depth,station)
 call mpp_write(unit,ocn_vehicle_field,profile%ocn_vehicle,station)
 call mpp_write(unit,database_id_field,profile%database_id,station)
 profile_flag = profile%profile_flag
 call mpp_write(unit,profile_flag_field,profile_flag,station)
 profile_flag = profile%profile_flag_s
 if (nvar_out .eq. 2) call mpp_write(unit,profile_flag_s_field,profile_flag,station)
 call mpp_write(unit,lon_field,profile%lon,station)
 call mpp_write(unit,lat_field,profile%lat,station)
 call mpp_write(unit,time_field,days_since,station)
 tmp_s = real(profile%yyyy)
 call mpp_write(unit,yyyy_field,tmp_s,station)
 tmp_s = real(profile%mmdd)
 call mpp_write(unit,mmdd_field,tmp_s,station)
 call mpp_write(unit,temp_err_field,profile%temp_err,station)
 if (nvar_out .eq. 2) call mpp_write(unit,salt_err_field,profile%salt_err,station)
 nlinks = 0

 if (profile%levels .gt. max_levels_file) then
     nlinks = ceiling(float(profile%levels)/float(max_levels_file)) - 1
 endif

 if (nlinks .gt. 0) then
     call mpp_write(unit,link_field,1.,station)
 else
     call mpp_write(unit,link_field,0.,station)
 endif

if (profile%i_index .ne. -1.0 .and. profile%j_index .ne. -1.0) then
   call mpp_write(unit, lon_index_field,profile%i_index)
   call mpp_write(unit, lat_index_field,profile%j_index)
endif

do i = 1, nlinks
   sta_num(findex)=sta_num(findex)+1
   station=sta_num(findex)
   if (i.eq.nlinks) then
       levels = mod(profile%levels,max_levels_file)
       if (levels .eq. 0) levels = max_levels_file
   else
       levels = max_levels_file
   endif
   data_t = missing_value; data_s = missing_value; depth = missing_value
   flag_t=missing_value;flag_s=missing_value

   data_t(1:levels)=profile%data_t((max_levels_file*i)+1:(max_levels_file*i)+levels)
   flag_t(1:levels)=profile%flag_t((max_levels_file*i)+1:(max_levels_file*i)+levels)

   if (profile%nvar == 2) then
       data_s(1:levels)   = profile%data_s((max_levels_file*i)+1:(max_levels_file*i)+levels)
       flag_s(1:levels)= profile%flag_s((max_levels_file*i)+1:(max_levels_file*i)+levels)

   endif


   depth(1:levels)=profile%depth((max_levels_file*i)+1:(max_levels_file*i)+levels)

   call mpp_write(unit,nvar_field,nvar,station)
   call mpp_write(unit,data_t_field,data_t,station)
   if (nvar_out .eq. 2) call mpp_write(unit,data_s_field,data_s,station)
   call mpp_write(unit,depth_field,depth,station)

   call mpp_write(unit,project_field,profile%project,station)
   call mpp_write(unit,probe_field,profile%probe,station)
   call mpp_write(unit,ref_inst_field,profile%ref_inst,station)
   call mpp_write(unit,fix_depth_field,profile%fix_depth,station)
   call mpp_write(unit,ocn_vehicle_field,profile%ocn_vehicle,station)
   call mpp_write(unit,database_id_field,profile%database_id,station)
   profile_flag = profile%profile_flag
   call mpp_write(unit,profile_flag_field,profile_flag,station)
   profile_flag = profile%profile_flag_s
   if (nvar_out .eq. 2)   call mpp_write(unit,profile_flag_s_field,profile_flag,station)
   call mpp_write(unit,lon_field,profile%lon,station)
   call mpp_write(unit,lat_field,profile%lat,station)
   call mpp_write(unit,time_field,days_since,station)
   tmp_s = real(profile%yyyy)
   call mpp_write(unit,yyyy_field,tmp_s,station)
   tmp_s = real(profile%mmdd)
   call mpp_write(unit,mmdd_field,tmp_s,station)
   call mpp_write(unit,temp_err_field,profile%temp_err,station)
   if (nvar_out .eq. 2)   call mpp_write(unit,salt_err_field,profile%salt_err,station)

   if (profile%i_index .ne. -1.0 .and. profile%j_index .ne. -1.0) then
       call mpp_write(unit, lon_index_field,profile%i_index)
       call mpp_write(unit, lat_index_field,profile%j_index)
   endif

   if (i .lt. nlinks) then
       call mpp_write(unit,link_field,1.,station)
   else
       call mpp_write(unit,link_field,0.,station)
   endif

enddo

end subroutine write_profile

subroutine close_profile_file(unit)

  integer, intent(in) :: unit

  call mpp_close(unit)

end subroutine close_profile_file

subroutine write_ocean_data_init()

  module_is_initialized=.true.

  sta_num=0;unit_num=0;nfiles=0

  return

end subroutine write_ocean_data_init

end module write_ocean_data_mod
