
module gaussian_topog_mod

use  fms_mod, only: file_exist, open_namelist_file,  &
                    check_nml_error, close_file,     &
                    stdlog, write_version_number,    &
                    mpp_pe, mpp_root_pe,             &
                    error_mesg, FATAL

use constants_mod, only: pi

implicit none
private

public :: gaussian_topog_init, get_gaussian_topog

!-----------------------------------------------------------------------
!----- namelist ------
!----- multiple mountains can be generated ------
!
!      height = height in meters
!      olon, olat = longitude,latitude origin              (degrees)
!      rlon, rlat = longitude,latitude half-width of ridge (degrees)
!      wlon, wlat = longitude,latitude half-width of tail  (degrees)
!
!      Note: For the standard gaussian mountain
!            set rlon = rlat = 0 .
!
!-----------------------------------------------------------------------

   integer, parameter :: maxmts = 10

   real, dimension(maxmts) :: height = 0.
   real, dimension(maxmts) ::  olon  = 0.
   real, dimension(maxmts) ::  olat  = 0.
   real, dimension(maxmts) ::  wlon  = 0.
   real, dimension(maxmts) ::  wlat  = 0.
   real, dimension(maxmts) ::  rlon  = 0.
   real, dimension(maxmts) ::  rlat  = 0.

   namelist /gaussian_topog_nml/ height, olon, olat, wlon, wlat, rlon, rlat

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: gaussian_topog.F90,v 1.2 2002/07/16 22:57:15 fms Exp $'
character(len=128) :: tagname = '$Name: havana $'

logical :: do_nml = .true.
logical :: module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine gaussian_topog_init ( lon, lat, zsurf )

real, intent(in)  :: lon(:), lat(:)
real, intent(out) :: zsurf(:,:)

!  lon, lat = mean grid box longitude and latitude. (radians)
!  zsurf    = surface height (meters)

integer :: n

  if (.not.module_is_initialized) then
     call write_version_number( version, tagname )
  endif

  if(any(shape(zsurf) /= (/size(lon),size(lat)/))) then
    call error_mesg ('get_gaussian_topog in topography_mod', &
     'shape(zsurf) is not equal to (/size(lon),size(lat)/)', FATAL)
  endif

  if (do_nml) call read_namelist

! compute sum of all non-zero mountains
  zsurf(:,:) = 0.
  do n = 1, maxmts
    if ( height(n) == 0. ) cycle
    zsurf = zsurf + get_gaussian_topog ( lon, lat, height(n), &
                olon(n), olat(n), wlon(n), wlat(n), rlon(n), rlat(n))
  enddo
 module_is_initialized = .TRUE.                    

end subroutine gaussian_topog_init

!#######################################################################

function get_gaussian_topog ( lon, lat, height,                          &
                              olond, olatd, wlond, wlatd, rlond, rlatd ) &
                     result ( zsurf )

real, intent(in)  :: lon(:), lat(:)
real, intent(in)  :: height
real, intent(in), optional :: olond, olatd, wlond, wlatd, rlond, rlatd
real :: zsurf(size(lon),size(lat))

!  lon, lat    = mean grid box longitude and latitude (in radians)
!  height      = mountain height in meters
!  olond,olatd = lon/lat of mountain position/origin (in degrees)
!  wlond,wlatd = lon/lat of mountain half-width      (in degrees)
!  rlond,rlatd = lon/lat of ridge half-width         (in degrees)
!  zsurf       = mountain height (meters)

integer :: i, j
real    :: olon, olat, wlon, wlat, rlon, rlat
real    :: tpi, dtr, dx, dy, xx, yy

  if (do_nml) call read_namelist

! no need to compute mountain if height=0
  if ( height == 0. ) then
       zsurf(:,:) = 0.
       return
  endif

  tpi = 2.0*pi
  dtr = tpi/360.

! defaults and convert degrees to radians (dtr)
  olon = 90.*dtr;  if (present(olond)) olon=olond*dtr
  olat = 45.*dtr;  if (present(olatd)) olat=olatd*dtr
  wlon = 15.*dtr;  if (present(wlond)) wlon=wlond*dtr
  wlat = 15.*dtr;  if (present(wlatd)) wlat=wlatd*dtr
  rlon =  0.    ;  if (present(rlond)) rlon=rlond*dtr
  rlat =  0.    ;  if (present(rlatd)) rlat=rlatd*dtr

! compute gaussian-shaped mountain
    do j=1,size(lat)
      dy = abs(lat(j) - olat)   ! dist from y origin
      yy = max(0., dy-rlat)/wlat
      do i=1,size(lon)
        dx = abs(lon(i) - olon) ! dist from x origin
        dx = min(dx, abs(dx-tpi))  ! To ensure that: -pi <= dx <= pi
        xx = max(0., dx-rlon)/wlon
        zsurf(i,j) = height*exp(-xx**2 - yy**2)
      enddo
    enddo

end function get_gaussian_topog

!#######################################################################

subroutine read_namelist

   integer :: unit, ierr, io
   real    :: dtr

!  read namelist

   if ( file_exist('input.nml')) then
      unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=gaussian_topog_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'gaussian_topog_nml')
      enddo
 10   call close_file (unit)
   endif

!  write version and namelist to log file

   if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=gaussian_topog_nml)

   do_nml = .false.

end subroutine read_namelist

!#######################################################################

end module gaussian_topog_mod

