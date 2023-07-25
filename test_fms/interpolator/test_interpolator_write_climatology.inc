subroutine write_climatology_file

  implicit none

  type(FmsNetcdfFile_t) :: fileobj
  integer :: i,j,k,l, itime
  integer :: stringlen
  character(100) :: mystring

  real :: dx, dy

  ! write netcdf file
  if(open_file(fileobj, 'INPUT'//trim(ncfile), 'overwrite')) then
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
     call register_field(fileobj, 'time', 'double', dimensions=(/'time'/))
     call register_string_attribute(fileobj, 'time', 'axis', 'T')
     call register_string_attribute(fileobj, 'time', 'calendar', 'noleap')
     call register_string_attribute(fileobj, 'time', 'climatology', '1979-01-01 00:00:00, 1998-01-01 00:00:00')
     call register_string_attribute(fileobj, 'time', 'units', 'days since 1849-01-01 00:00:00')

     !lon axes
     call register_field(fileobj, 'lon', 'float', dimensions=(/'longitude'/))
     call register_string_attribute(fileobj, 'lon', 'edges', 'lonb')
     call register_string_attribute(fileobj, 'lon', 'cartesian_axis', 'X')
     call register_string_attribute(fileobj, 'lon', 'units', 'degrees_east')

     !lonb  axes longitude edges
     call register_field(fileobj, 'lonb', 'float', dimensions=(/'lonb'/))
     call register_string_attribute(fileobj, 'lonb', 'cartesian_axis','X')
     call register_string_attribute(fileobj, 'lonb', 'units', 'degrees_edges')

     !lat
     call register_field(fileobj, 'lat', 'float', dimensions=(/'latitude'/))
     call register_string_attribute(fileobj, 'lat', 'edges', 'latb')
     call register_string_attribute(fileobj, 'lat', 'cartesian_axis', 'Y')
     call register_string_attribute(fileobj, 'lat', 'units', 'degrees_north')

     !latb latitude_edges
     call register_field(fileobj, 'latb', 'float', dimensions=(/'latb'/))
     call register_string_attribute(fileobj, 'latb', 'cartesian_axis', 'Y')
     call register_string_attribute(fileobj, 'latb', 'units', 'degrees_N')

     !pfull layer-mean pressure
     call register_field(fileobj, 'pfull', 'float', dimensions=(/'pfull'/))
     call register_string_attribute(fileobj, 'pfull', 'positive', 'down')
     call register_string_attribute(fileobj, 'pfull', 'cartesian_axis', 'Z')
     call register_string_attribute(fileobj, 'pfull', 'units', 'hPa')

     !phalf half pressure level
     call register_field(fileobj, 'phalf', 'float', dimensions=(/'phalf'/))
     call register_string_attribute(fileobj, 'phalf', 'positive', 'down')
     call register_string_attribute(fileobj, 'phalf', 'cartesian_axis', 'Z')
     call register_string_attribute(fileobj, 'phalf', 'units', 'hPa')

     !ozone ozone mass mixing ratio
     call register_field(fileobj, 'ozone', 'float', dimensions=(/'lon  ','lat  ','pfull','time '/))
     call register_string_attribute(fileobj, 'ozone', 'units', 'kg/kg')

     call write_data(fileobj, 'time',  clim_time)
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
subroutine allocate_arrays()

  implicit none
  allocate( lat(nlonlat),   &
            lon(nlonlat),   &
            latb(nlonlatb), &
            lonb(nlonlatb), &
            clim_time(ntime),    &
            pfull(npfull),  &
            phalf(nphalf),  &
            ozone(nlonlat, nlonlat, npfull, ntime) )

  allocate( lat_mod(nlonlat_mod, nlonlat_mod),    &
            lon_mod(nlonlat_mod, nlonlat_mod),    &
            latb_mod(nlonlatb_mod, nlonlatb_mod), &
            lonb_mod(nlonlatb_mod, nlonlatb_mod) )

end subroutine allocate_arrays
!===============================================!
subroutine deallocate_arrays()

  implicit none
  deallocate( lat,   &
              lon,   &
              latb,  &
              lonb,  &
              clim_time,  &
              pfull, &
              phalf, &
              lat_mod,   &
              lon_mod,   &
              latb_mod,  &
              lonb_mod,  &
              ozone )

end subroutine deallocate_arrays
!===============================================!
subroutine set_clim_time()

  implicit none
  integer :: i

  do i=1, ntime
     clim_time(i)=2*i
  end do

end subroutine set_clim_time
!===============================================!
subroutine set_latlon_b()

  implicit none
  integer :: i

  do i=1, nlonlat
     lat(i)= real(2*i-1)
     lon(i)= real(2*i-1)
  end do

  lonb(1)=0.0
  latb(1)=0.0
  do i=2,nlonlat
     lonb(i)=0.5*(lon(i-1)+lon(i))
     latb(i)=0.5*(lat(i-1)+lat(i))
  enddo
  lonb(nlonlat+1)=real(2*nlonlat-1)
  latb(nlonlat+1)=real(2*nlonlat-1)

end subroutine set_latlon_b
!===============================================!
!> create model lat(b) and lon(b) arrays
subroutine set_latlon_b_mod

  implicit none
  integer :: i

  !> the model grid and the climatology grid are the same
  do i=1, nlonlat_mod
     lon_mod(:,i) = real(2*i-1)
     lat_mod(i,:) = real(2*i-1)
  end do

  lonb_mod(1:nlonlat,1)=0.0
  latb_mod(1,1:nlonlat)=0.0
  do i=2, nlonlatb_mod
     lonb_mod(1:nlonlat,i)=0.5*(lon_mod(:,i-1)+lon_mod(:,i))
     latb_mod(i,1:nlonlat)=0.5*(lat_mod(i-1,:)+lat_mod(i,:))
  end do
  lonb_mod(1:nlonlat,nlonlat+1)=real(2*nlonlat-1)
  latb_mod(nlonlat+1,1:nlonlat)=real(2*nlonlat-1)

  !> convert from degrees to radians
  lon_mod  = lon_mod*PI/180.0
  lat_mod  = lat_mod*PI/180.0
  lonb_mod = lonb_mod*PI/180.0
  latb_mod = latb_mod*PI/180.0

end subroutine set_latlon_b_mod
!===============================================!
subroutine set_pfullhalf()

  implicit none
  integer :: i

  do i=1, npfull
     pfull(i) = 0.0001 * real(i-1)
  end do

  phalf(1)=0.0
  do i=2, npfull
     phalf(i)=0.5*( pfull(i)+pfull(i-1) )
  end do
  phalf(npfull+1)=pfull(npfull)

end subroutine set_pfullhalf
!===============================================!
subroutine set_ozone()

  implicit none
  integer :: i, j, k, l

  do i=1, ntime
     do j=1, npfull
        do k=1, nlonlat
           do l=1, nlonlat
              ozone(l,k,j,i)=real(i)
           end do
        end do
     end do
  end do

end subroutine set_ozone
!===============================================!
