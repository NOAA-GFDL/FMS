module get_grid_version_mpp_mod
use constants_mod, only: PI
use mpp_mod, only : mpp_error,FATAL,WARNING,NOTE, mpp_min, mpp_max
use fms_io_mod, only: field_size, read_data, get_mosaic_tile_grid
use fms_mod, only: field_exist
use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, operator(.NE.),operator(.EQ.)
use mpp_domains_mod, only : mpp_copy_domain, mpp_get_global_domain
use mpp_domains_mod, only : mpp_get_data_domain, mpp_set_compute_domain, mpp_set_data_domain
use mpp_domains_mod, only : mpp_set_global_domain, mpp_deallocate_domain

implicit none

real, parameter    :: deg_to_radian=PI/180.
contains
! Get lon and lat of three model (target) grids from grid_spec.nc

subroutine check_grid_sizes(domain_name, Domain, nlon, nlat)
character(len=12), intent(in) :: domain_name
type (domain2d),   intent(in) :: Domain
integer,           intent(in) :: nlon, nlat

character(len=184) :: error_message
integer            :: xsize, ysize

call mpp_get_global_domain(Domain, xsize=xsize, ysize=ysize)
if(nlon .NE. xsize .OR. nlat .NE. ysize) then
  error_message = 'Error in data_override_init. Size of grid as specified by '// &
                  '             does not conform to that specified by grid_spec.nc.'// &
                  '  From             :     by      From grid_spec.nc:     by    '
  error_message( 59: 70) = domain_name
  error_message(130:141) = domain_name
  write(error_message(143:146),'(i4)') xsize
  write(error_message(150:153),'(i4)') ysize
  write(error_message(174:177),'(i4)') nlon
  write(error_message(181:184),'(i4)') nlat
  call mpp_error(FATAL,error_message)
endif
end subroutine check_grid_sizes

subroutine get_grid_version_classic_1(grid_file, mod_name, domain, isc, iec, jsc, jec, lon, lat, min_lon, max_lon, grid_center_bug)
  character(len=*),            intent(in) :: grid_file
  character(len=*),            intent(in) :: mod_name
  type(domain2d),              intent(in) :: domain
  integer,                     intent(in) :: isc, iec, jsc, jec
  real, dimension(isc:,jsc:), intent(out) :: lon, lat
  real,                       intent(out) :: min_lon, max_lon
  logical,           intent(in), optional :: grid_center_bug

  integer                                      :: i, j, siz(4)
  integer                                      :: nlon, nlat         ! size of global lon and lat
  real,          dimension(:,:,:), allocatable :: lon_vert, lat_vert !of OCN grid vertices
  real,          dimension(:),     allocatable :: glon, glat         ! lon and lat of 1-D grid of atm/lnd
  logical                                      :: is_new_grid
  integer                                      :: is, ie, js, je
  integer                                      :: isd, ied, jsd, jed
  integer                                      :: isg, ieg, jsg, jeg
  type(domain2d)                               :: domain2
  character(len=3)                             :: xname, yname

  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)

  select case(mod_name)
  case('ocn', 'ice')
    is_new_grid = .FALSE.
    if(field_exist(grid_file, 'x_T')) then
       is_new_grid = .true.
    else if(field_exist(grid_file, 'geolon_t')) then
       is_new_grid = .FALSE.
    else
       call mpp_error(FATAL,'data_override: both x_T and geolon_t is not in the grid file '//trim(grid_file) )
    endif

    if(is_new_grid) then
      call field_size(grid_file, 'x_T', siz)
      nlon = siz(1); nlat = siz(2)
      call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)
      allocate(lon_vert(isc:iec,jsc:jec,4), lat_vert(isc:iec,jsc:jec,4) )
      call read_data(trim(grid_file), 'x_vert_T', lon_vert, domain)
      call read_data(trim(grid_file), 'y_vert_T', lat_vert, domain)

!2 Global lon and lat of ocean grid cell centers are determined from adjacent vertices
      lon(:,:) = (lon_vert(:,:,1) + lon_vert(:,:,2) + lon_vert(:,:,3) + lon_vert(:,:,4))*0.25
      lat(:,:) = (lat_vert(:,:,1) + lat_vert(:,:,2) + lat_vert(:,:,3) + lat_vert(:,:,4))*0.25
    else
      if(grid_center_bug) call mpp_error(NOTE, &
           'data_override: grid_center_bug is set to true, the grid center location may be incorrect')
      call field_size(grid_file, 'geolon_vert_t', siz)
      nlon = siz(1) - 1; nlat = siz(2) - 1;
      call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)
      call mpp_copy_domain(domain, domain2)
      call mpp_set_compute_domain(domain2, isc, iec+1, jsc, jec+1, iec-isc+2, jec-jsc+2 )
      call mpp_set_data_domain   (domain2, isd, ied+1, jsd, jed+1, ied-isd+2, jed-jsd+2 )
      call mpp_set_global_domain (domain2, isg, ieg+1, jsg, jeg+1, ieg-isg+2, jeg-jsg+2 )
      allocate(lon_vert(isc:iec+1,jsc:jec+1,1))
      allocate(lat_vert(isc:iec+1,jsc:jec+1,1))
      call read_data(trim(grid_file), 'geolon_vert_t', lon_vert, domain2)
      call read_data(trim(grid_file), 'geolat_vert_t', lat_vert, domain2)

      if(grid_center_bug) then
         do j = jsc, jec
            do i = isc, iec
               lon(i,j) = (lon_vert(i,j,1) + lon_vert(i+1,j,1))/2.
               lat(i,j) = (lat_vert(i,j,1) + lat_vert(i,j+1,1))/2.
            enddo
         enddo
      else
         do j = jsc, jec
            do i = isc, iec
               lon(i,j) = (lon_vert(i,j,1) + lon_vert(i+1,j,1) + &
                    lon_vert(i+1,j+1,1) + lon_vert(i,j+1,1))*0.25
               lat(i,j) = (lat_vert(i,j,1) + lat_vert(i+1,j,1) + &
                    lat_vert(i+1,j+1,1) + lat_vert(i,j+1,1))*0.25
            enddo
         enddo
      end if
      call mpp_deallocate_domain(domain2)
    endif
    deallocate(lon_vert)
    deallocate(lat_vert)
  case('atm', 'lnd')
     if(trim(mod_name) == 'atm') then
        xname = 'xta'; yname = 'yta'
     else
        xname = 'xtl'; yname = 'ytl'
     endif
     call field_size(grid_file, xname, siz)
     nlon = siz(1); allocate(glon(nlon))
     call read_data(grid_file, xname, glon, no_domain = .true.)

     call field_size(grid_file, yname, siz)
     nlat = siz(1); allocate(glat(nlat))
     call read_data(grid_file, yname, glat, no_domain = .true.)
     call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)

     is = isc - isg + 1; ie = iec - isg + 1
     js = jsc - jsg + 1; je = jec - jsg + 1
     do j = js, jec
        do i = is, ie
           lon(i,j) = glon(i)
           lat(i,j) = glat(j)
        enddo
     enddo
     deallocate(glon)
     deallocate(glat)
  case default
     call mpp_error(FATAL, "data_override_mod: mod_name should be 'atm', 'ocn', 'ice' or 'lnd' ")
  end select

  ! convert from degree to radian
  lon = lon * deg_to_radian
  lat = lat * deg_to_radian
  min_lon = minval(lon)
  max_lon = maxval(lon)
  call mpp_min(min_lon)
  call mpp_max(max_lon)


end subroutine get_grid_version_classic_1

! Get global lon and lat of three model (target) grids from mosaic.nc
! z1l: currently we assume the refinement ratio is 2 and there is one tile on each pe.
subroutine get_grid_version_classic_2(mosaic_file, mod_name, domain, isc, iec, jsc, jec, lon, lat, min_lon, max_lon)
  character(len=*),            intent(in) :: mosaic_file
  character(len=*),            intent(in) :: mod_name
  type(domain2d),              intent(in) :: domain
  integer,                     intent(in) :: isc, iec, jsc, jec
  real, dimension(isc:,jsc:), intent(out) :: lon, lat
  real,                       intent(out) :: min_lon, max_lon

  integer            :: i, j, siz(4)
  integer            :: nlon, nlat             ! size of global grid
  integer            :: nlon_super, nlat_super ! size of global supergrid.
  integer            :: isd, ied, jsd, jed
  integer            :: isg, ieg, jsg, jeg
  integer            :: isc2, iec2, jsc2, jec2
  character(len=256) :: solo_mosaic_file, grid_file
  real, allocatable  :: tmpx(:,:), tmpy(:,:)
  type(domain2d)     :: domain2

  if(trim(mod_name) .NE. 'atm' .AND. trim(mod_name) .NE. 'ocn' .AND. &
     trim(mod_name) .NE. 'ice' .AND. trim(mod_name) .NE. 'lnd' ) call mpp_error(FATAL, &
        "data_override_mod: mod_name should be 'atm', 'ocn', 'ice' or 'lnd' ")

  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)

  ! get the grid file to read
  if(field_exist(mosaic_file, trim(mod_name)//'_mosaic_file' )) then
     call read_data(mosaic_file, trim(mod_name)//'_mosaic_file', solo_mosaic_file)
     solo_mosaic_file = 'INPUT/'//trim(solo_mosaic_file)
  else
     solo_mosaic_file = mosaic_file
  end if
  call get_mosaic_tile_grid(grid_file, solo_mosaic_file, domain)

  call field_size(grid_file, 'area', siz)
  nlon_super = siz(1); nlat_super = siz(2)
  if( mod(nlon_super,2) .NE. 0) call mpp_error(FATAL,  &
       'data_override_mod: '//trim(mod_name)//' supergrid longitude size can not be divided by 2')
  if( mod(nlat_super,2) .NE. 0) call mpp_error(FATAL,  &
       'data_override_mod: '//trim(mod_name)//' supergrid latitude size can not be divided by 2')
  nlon = nlon_super/2;
  nlat = nlat_super/2;
  call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)

  !--- setup the domain for super grid.
  call mpp_copy_domain(domain, domain2)
  call mpp_set_compute_domain(domain2, 2*isc-1, 2*iec+1, 2*jsc-1, 2*jec+1, 2*iec-2*isc+3, 2*jec-2*jsc+3 )
  call mpp_set_data_domain   (domain2, 2*isd-1, 2*ied+1, 2*jsd-1, 2*jed+1, 2*ied-2*isd+3, 2*jed-2*jsd+3 )
  call mpp_set_global_domain (domain2, 2*isg-1, 2*ieg+1, 2*jsg-1, 2*jeg+1, 2*ieg-2*isg+3, 2*jeg-2*jsg+3 )

  call mpp_get_compute_domain(domain2, isc2, iec2, jsc2, jec2)
  if(isc2 .NE. 2*isc-1 .OR. iec2 .NE. 2*iec+1 .OR. jsc2 .NE. 2*jsc-1 .OR. jec2 .NE. 2*jec+1) then
     call mpp_error(FATAL, 'data_override_mod: '//trim(mod_name)//' supergrid domain is not set properly')
  endif

  allocate(tmpx(isc2:iec2, jsc2:jec2), tmpy(isc2:iec2, jsc2:jec2) )
  call read_data( grid_file, 'x', tmpx, domain2)
  call read_data( grid_file, 'y', tmpy, domain2)
  ! copy data onto model grid
  if(trim(mod_name) == 'ocn' .OR. trim(mod_name) == 'ice') then
     do j = jsc, jec
        do i = isc, iec
           lon(i,j) = (tmpx(i*2-1,j*2-1)+tmpx(i*2+1,j*2-1)+tmpx(i*2+1,j*2+1)+tmpx(i*2-1,j*2+1))*0.25
           lat(i,j) = (tmpy(i*2-1,j*2-1)+tmpy(i*2+1,j*2-1)+tmpy(i*2+1,j*2+1)+tmpy(i*2-1,j*2+1))*0.25
        end do
     end do
  else
     do j = jsc, jec
        do i = isc, iec
           lon(i,j) = tmpx(i*2,j*2)
           lat(i,j) = tmpy(i*2,j*2)
        end do
     end do
  endif

  ! convert to radian
  lon = lon * deg_to_radian
  lat = lat * deg_to_radian

  deallocate(tmpx, tmpy)
  min_lon = minval(lon)
  max_lon = maxval(lon)
  call mpp_min(min_lon)
  call mpp_max(max_lon)

  call mpp_deallocate_domain(domain2)

end subroutine get_grid_version_classic_2

end module get_grid_version_mpp_mod

