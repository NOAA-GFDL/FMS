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

!> This file contains subroutines that initialize data and write out the climatology file.

subroutine write_climatology_file

  implicit none

  type(FmsNetcdfFile_t) :: fileobj

  ! write netcdf file
  if(open_file(fileobj, 'INPUT/'//trim(ncfile), 'overwrite')) then
     call register_axis(fileobj, 'time', UNLIMITED)
     call register_axis(fileobj, 'lon',   nlonlat)
     call register_axis(fileobj, 'lat',   nlonlat)
     call register_axis(fileobj, 'latb',  nlonlatb)
     call register_axis(fileobj, 'lonb',  nlonlatb)
     call register_axis(fileobj, 'pfull', npfull)
     call register_axis(fileobj, 'phalf', nphalf)
     call register_axis(fileobj, 'longitude', nlonlat)
     call register_axis(fileobj, 'latitude',  nlonlat)

     !time axes
     !if(ntime /=0 ) then
     call register_field(fileobj, 'time', 'double', dimensions=(/'time'/))
     call register_string_attribute(fileobj, 'time', 'axis', 'T')
     if(noleap) call register_string_attribute(fileobj, 'time', 'calendar', 'noleap')
     if(.not.noleap) call register_string_attribute(fileobj, 'time', 'calendar', 'julian')
     if(daily)  call register_string_attribute(fileobj, 'time', 'units', 'days since 1849-01-01 01:01:01')
     if(yearly) call register_string_attribute(fileobj, 'time', 'units', 'years since 1849-01-01 01:01:01')
     !end if

     !lon axes
     call register_field(fileobj, 'lon', 'double', dimensions=(/'longitude'/))
     call register_string_attribute(fileobj, 'lon', 'edges', 'lonb')
     call register_string_attribute(fileobj, 'lon', 'cartesian_axis', 'X')
     call register_string_attribute(fileobj, 'lon', 'units', 'degrees_east')

     !lonb  axes longitude edges
     call register_field(fileobj, 'lonb', 'double', dimensions=(/'lonb'/))
     call register_string_attribute(fileobj, 'lonb', 'cartesian_axis','X')
     call register_string_attribute(fileobj, 'lonb', 'units', 'degrees_edges')

     !lat
     call register_field(fileobj, 'lat', 'double', dimensions=(/'latitude'/))
     call register_string_attribute(fileobj, 'lat', 'edges', 'latb')
     call register_string_attribute(fileobj, 'lat', 'cartesian_axis', 'Y')
     call register_string_attribute(fileobj, 'lat', 'units', 'degrees_north')

     !latb latitude_edges
     call register_field(fileobj, 'latb', 'double', dimensions=(/'latb'/))
     call register_string_attribute(fileobj, 'latb', 'cartesian_axis', 'Y')
     call register_string_attribute(fileobj, 'latb', 'units', 'degrees_N')

     !pfull layer-mean pressure
     call register_field(fileobj, 'pfull', 'double', dimensions=(/'pfull'/))
     call register_string_attribute(fileobj, 'pfull', 'positive', 'down')
     call register_string_attribute(fileobj, 'pfull', 'cartesian_axis', 'Z')
     call register_string_attribute(fileobj, 'pfull', 'units', 'hPa')

     !phalf half pressure level
     call register_field(fileobj, 'phalf', 'double', dimensions=(/'phalf'/))
     call register_string_attribute(fileobj, 'phalf', 'positive', 'down')
     call register_string_attribute(fileobj, 'phalf', 'cartesian_axis', 'Z')
     call register_string_attribute(fileobj, 'phalf', 'units', 'hPa')

     !ozone ozone mass mixing ratio
     if( ntime==0 ) then
        call register_field(fileobj, 'ozone', 'double', dimensions=(/'lon  ','lat  ','pfull'/))
     else
        call register_field(fileobj, 'ozone', 'double', dimensions=(/'lon  ','lat  ','pfull','time '/))
     end if
     call register_string_attribute(fileobj, 'ozone', 'units', 'kg/kg')

     if(ntime /= 0 ) call write_data(fileobj, 'time',  clim_time)
     call write_data(fileobj, 'lon',   lon)
     call write_data(fileobj, 'lonb',  lonb)
     call write_data(fileobj, 'lat',   lat)
     call write_data(fileobj, 'latb',  latb)
     call write_data(fileobj, 'pfull', pfull)
     call write_data(fileobj, 'phalf', phalf)
     call write_data(fileobj, 'ozone', ozone)

     call close_file(fileobj)

  end if

end subroutine write_climatology_file
!===============================================!
subroutine register_string_attribute(fileobj, variable, attribute, att_value)

  implicit none

  type(FmsNetcdfFile_t) :: fileobj
  character(*) :: variable, attribute, att_value
  integer :: str_len

  str_len=len(trim(att_value))
  call register_variable_attribute(fileobj, trim(variable), trim(attribute), trim(att_value), str_len=str_len)

end subroutine register_string_attribute
!===============================================!
subroutine deallocate_arrays()

  implicit none
  deallocate( lat,  &
              lon,  &
              latb, &
              lonb, &
              clim_time, &
              model_time_noleap, &
              model_time_julian, &
              pfull, &
              phalf, &
              lat_mod,  &
              lon_mod,  &
              latb_mod, &
              lonb_mod, &
              ozone )

end subroutine deallocate_arrays
!===============================================!
subroutine set_parameters(nlonlat_in, nlonlat_mod_in, ntime_in, npfull_in, daily_in, yearly_in, noleap_in)

  implicit none
  integer, intent(in) :: nlonlat_in, nlonlat_mod_in, ntime_in, npfull_in
  logical, intent(in) :: daily_in, yearly_in, noleap_in

  nlonlat  = nlonlat_in   !< number of latitude and longitudinal center coordinates
  nlonlatb = nlonlat_in+1 !< number of latitude and longitudinal boundary coordinates
  nlonlat_mod  = nlonlat_mod_in
  nlonlatb_mod = nlonlat_mod_in+1
  ntime  = ntime_in    !< number of time slices
  npfull = npfull_in   !< number of p levels
  nphalf = npfull_in+1 !< number of half p levels

  daily  = daily_in
  yearly = yearly_in
  noleap = noleap_in

end subroutine set_parameters
!===============================================!
subroutine set_and_write_data

  implicit none

  if( ntime /= 0) call set_clim_time()
  call set_latlon_b()
  call set_latlon_b_mod()
  call set_pfullhalf()
  call set_ozone()
  call write_climatology_file

end subroutine set_and_write_data
!===============================================!
subroutine set_clim_time()

  implicit none
  integer :: i

  type(time_type) :: base_time
  integer :: l, ii, yr, mo, dy, hr, mn, sc
  integer :: ntime2

  if(allocated(clim_time)) deallocate(clim_time)
  if(allocated(model_time_julian)) deallocate(model_time_julian)
  if(allocated(model_time_noleap)) deallocate(model_time_noleap)
  allocate(clim_time(ntime))
  allocate(model_time_julian(ntime), model_time_noleap(ntime))

  !write(*,*) ' -- SETTING TIME -- '

  hr = 1 ; mn = 1 ; sc = 1

  if(test_file_daily_noleap) then
     ! base_time must match that in the time attribute
     base_time=set_date_no_leap(1849,1,1,hr,mn,sc)
     yr = 1849 ; mo = 1 ; dy = 15
     do i=1, ntime
        yr = yr + 1
        mo = mo + 1 ; if( mo > 12 ) mo=1
        dy = 15
        model_time_julian(i)=set_date_julian(yr, mo, dy, hr, mn, sc)
        model_time_noleap(i)=set_date_no_leap(yr, mo, dy, hr, mn, sc)
        clim_time(i)=time_type_to_real(model_time_noleap(i)-base_time)/86400._r8_kind
        !call print_time(model_time_noleap(i))
     end do
  else if(test_file_daily_julian) then
     ! base_time must match that in the time attribute
     base_time=set_date_julian(1849,1,1,hr,mn,sc)
     yr = 1849 ; mo = 1 ; dy = 15
     do i=1, ntime
        yr = yr + 1
        mo = mo + 1 ; if( mo > 12 ) mo=1
        dy = 15
        model_time_julian(i)=set_date_julian(yr, mo, dy, hr, mn, sc)
        model_time_noleap(i)=set_date_no_leap(yr, mo, dy, hr, mn, sc)
        clim_time(i)=time_type_to_real(model_time_julian(i)-base_time)/86400._r8_kind
        !call print_time(model_time_julian(i))
     end do
  else if(test_file_yearly_noleap) then
     ! base_time must match that in the time attribute
     base_time=set_date_no_leap(1849,1,1,hr,mn,sc)
     yr = 1849 ; mo = 1 ; dy = 15
     do i=1, ntime
        yr = yr + 1
        mo = mo + 1 ; if( mo > 12 ) mo=1
        dy = 15
        model_time_julian(i)=set_date_julian(yr, mo, dy, hr, mn, sc)
        model_time_noleap(i)=set_date_no_leap(yr, mo, dy, hr, mn, sc)
        clim_time(i)=real(yr-1849,r8_kind)
        do ii=1, mo-1
           select case(ii)
           case(1,3,5,7,8,10,12) ; clim_time(i) = clim_time(i) + 31.0_r8_kind/365._r8_kind
           case(2) ; clim_time(i) = clim_time(i) + 28._r8_kind/365._r8_kind
           case(4,6,9,11) ; clim_time(i) = clim_time(i) + 30.0_r8_kind/365._r8_kind
           end select
        end do
        clim_time(i) = clim_time(i) + real(dy-1,r8_kind)/365._r8_kind
        !call print_time(model_time_noleap(i))
     end do
  else if(test_file_yearly_julian) then
     ! base_time must match that in the time attribute
     base_time=set_date_julian(1849,1,1,hr,mn,sc)
     yr = 1849 ; mo = 1 ; dy = 15
     call set_calendar_type(2)
     do i=1, ntime
        yr = yr + 1
        mo = mo + 1 ; if( mo > 12 ) mo=1
        dy = 15
        model_time_julian(i)=set_date_julian(yr,mo,dy,hr,mn,sc)
        model_time_noleap(i)=set_date_no_leap(yr,mo,dy,hr,mn,sc)
        clim_time(i)=real(yr-1849,lkind)
        l=0 ; if(leap_year(model_time_julian(i))) l=l+1
        do ii=1, mo-1
           select case(ii)
           case(1,3,5,7,8,10,12) ; clim_time(i) = clim_time(i) + 31.0_r8_kind/real(365+l,r8_kind)
           case(2) ; clim_time(i) = clim_time(i) + real(28+l,r8_kind)/real(365+l,r8_kind)
           case(4,6,9,11) ; clim_time(i) = clim_time(i) + 30.0_r8_kind/real(365+l,r8_kind)
           end select
        end do
        clim_time(i) = clim_time(i) + real(dy-1,r8_kind)/real(365+l,r8_kind)
        !call print_time(model_time_julian(i))
     end do
     call set_calendar_type(0)
  end if

end subroutine set_clim_time
!===============================================!
subroutine set_latlon_b()

  implicit none
  integer :: i

  if(allocated(lat))  deallocate(lat)
  if(allocated(lon))  deallocate(lon)
  if(allocated(latb)) deallocate(latb)
  if(allocated(lonb)) deallocate(lonb)

  allocate(lat(nlonlat), lon(nlonlat))
  allocate(latb(nlonlat+1), lonb(nlonlat+1))

  do i=1, nlonlat
     lat(i)= real(2*i-1,TEST_INTP_KIND_)
     lon(i)= real(2*i-1,TEST_INTP_KIND_)
  end do

  lonb(1)=1.0_lkind
  latb(1)=1.0_lkind
  do i=2,nlonlat
     lonb(i)=0.5_lkind*(lon(i-1)+lon(i))
     latb(i)=0.5_lkind*(lat(i-1)+lat(i))
  enddo
  lonb(nlonlat+1)=real(2*nlonlat-1)
  latb(nlonlat+1)=real(2*nlonlat-1)

end subroutine set_latlon_b
!===============================================!
!> create model lat(b) and lon(b) arrays
subroutine set_latlon_b_mod

  implicit none
  integer :: i, j

  if(allocated(lat_mod))  deallocate(lat_mod)
  if(allocated(lon_mod))  deallocate(lon_mod)
  if(allocated(latb_mod)) deallocate(latb_mod)
  if(allocated(lonb_mod)) deallocate(lonb_mod)

  allocate(lat_mod(nlonlat_mod,nlonlat_mod), lon_mod(nlonlat_mod,nlonlat_mod))
  allocate(latb_mod(nlonlatb_mod,nlonlatb_mod), lonb_mod(nlonlatb_mod,nlonlatb_mod))

  !> nlonlat = nlonlat_mod
  ! the model coordinates are the same as the file coordinates
  do i=1, nlonlat_mod
     lat_mod(i,:) = lat(i)
     lon_mod(:,i) = lon(i)
  end do
  !> nlonlatb_mod = nlonlatb
  !! the model coordinates are the same as the file coordinates
  do i=1, nlonlatb_mod
     latb_mod(i,:) = latb(i)
     lonb_mod(:,i) = lonb(i)
  end do

  !> convert from degrees to radians
  lon_mod  = lon_mod*real(PI,TEST_INTP_KIND_)/180.0_lkind
  lat_mod  = lat_mod*real(PI,TEST_INTP_KIND_)/180.0_lkind
  lonb_mod = lonb_mod*real(PI,TEST_INTP_KIND_)/180.0_lkind
  latb_mod = latb_mod*real(PI,TEST_INTP_KIND_)/180.0_lkind


end subroutine set_latlon_b_mod
!===============================================!
subroutine set_pfullhalf()

  implicit none
  integer :: i

  if(allocated(pfull)) deallocate(pfull)
  if(allocated(phalf)) deallocate(phalf)
  allocate(pfull(npfull), phalf(nphalf))

  do i=1, npfull
     pfull(i) = 0.0001_lkind * real(i-1,TEST_INTP_KIND_)
  end do

  phalf(1)=0.0_lkind
  do i=2, npfull
     phalf(i)=0.5_lkind*( pfull(i)+pfull(i-1) )
  end do
  phalf(npfull+1)=pfull(npfull)

end subroutine set_pfullhalf
!===============================================!
subroutine set_ozone()

  implicit none
  integer :: i, j, k, l
  integer :: ntime2

  ntime2=ntime
  if(ntime == 0 ) ntime2=1
  if(allocated(ozone)) deallocate(ozone)
  allocate(ozone(nlonlat, nlonlat, npfull, ntime2))
  do i=1, ntime2
     do j=1, npfull
        do k=1, nlonlat
           do l=1, nlonlat
              ozone(l,k,j,i)= real(i,TEST_INTP_KIND_)
           end do
        end do
     end do
  end do

end subroutine set_ozone
!===============================================!
