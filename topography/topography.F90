
module topography_mod

use horiz_interp_mod, only: horiz_interp

use    utilities_mod, only: file_exist, open_file, check_nml_error, &
                            get_my_pe, close_file, error_mesg, FATAL

use    constants_mod, only: grav

implicit none
private

public :: simple_mountain, compute_real_topog, compute_stdev, compute_ocean_mask

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
!
!      height -->   ___________________________
!                  /                           \
!                 /              |              \
!   gaussian     /               |               \
!     sides --> /                |                \
!              /               olon                \
!        _____/                olat                 \______
!                    
!             |    |             |
!             |<-->|<----------->|
!             |wlon|    rlon     |
!              wlat     rlat      
!
!-----------------------------------------------------------------------

   integer, parameter :: maxmts = 100

   real, dimension(maxmts) :: height = 0.
   real, dimension(maxmts) ::  olon  = 0.
   real, dimension(maxmts) ::  olat  = 0.
   real, dimension(maxmts) ::  wlon  = 0.
   real, dimension(maxmts) ::  wlat  = 0.
   real, dimension(maxmts) ::  rlon  = 0.
   real, dimension(maxmts) ::  rlat  = 0.

   namelist /topography_nml/ height, olon, olat, wlon, wlat, rlon, rlat

!-----------------------------------------------------------------------
! --- resolution of the topography data set ---

  integer, parameter :: ipts = 2160, jpts = 1080

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: topography.F90,v 1.2 2001/07/05 17:02:33 fms Exp $'
character(len=128) :: tag = '$Name: fez $'

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine simple_mountain ( lon, lat, surf_geo )

real, intent(in)  :: lon(:), lat(:)
real, intent(out) :: surf_geo(:,:)

!  lon, lat = mean grid box longitude and latitude. (radians)
!  surf_geo = surface geopotential height. (meters/seconds)**2

integer :: i, j, n
real    :: pi, dx, dy, xxx, yyy

if(any(shape(surf_geo) /= (/size(lon),size(lat)/))) then
  call error_mesg ('simple_mountain in topography_mod', &
   'shape(surf_geo) is not equal to (/size(lon),size(lat)/)', FATAL)
endif

pi = 2*acos(0.)
call read_namelist

surf_geo(:,:) = 0.

do n = 1, maxmts
  if ( height(n) == 0. ) cycle
  do j=1,size(lat)
    dy = abs(lat(j) - olat(n))
    yyy = max(0., dy-rlat(n))/wlat(n)
    do i=1,size(lon)
      dx = abs(lon(i) - olon(n))
      dx = min(dx, abs(dx-2*pi)) ! To ensure that: -pi <= dx <= pi
      xxx = max(0., dx-rlon(n))/wlon(n)
      surf_geo(i,j) = surf_geo(i,j) + grav*height(n)*exp(-xxx**2 - yyy**2)
    enddo
  enddo
enddo

return
end subroutine simple_mountain

!#######################################################################

 function compute_real_topog (blon, blat, surf_geo)

   real, intent(in),  dimension(:)   :: blon, blat
   real, intent(out), dimension(:,:) :: surf_geo
   logical :: compute_real_topog

!  blon, blat = longitude and latitude in radians for grid box edges
!  surf_geo   = surface geopotential (meters/seconds)**2

   real :: zdata (ipts,jpts), hpi, del
!-----------------------------------------------------------------------

   if ( read_topog_file (zdata) ) then

      if ( any(shape(surf_geo) /= (/size(blon)-1,size(blat)-1/)) ) then
        call error_mesg('compute_real_topog','shape(surf_geo) is not&
          & equal to (/size(blon)-1,size(blat)-1/))', FATAL)
      endif

      hpi = acos(0.0)
      del = 2.*hpi/real(jpts)
      call horiz_interp ( zdata, 0.0, -hpi, del, del, blon, blat, surf_geo )
      surf_geo = surf_geo*grav

      compute_real_topog = .true.
   else
      compute_real_topog = .false.
   endif

!-----------------------------------------------------------------------

 end function compute_real_topog

!#######################################################################

function compute_ocean_mask(blon, blat, ocean_mask)

real   , intent(in),  dimension(:)   :: blon, blat
logical, intent(out), dimension(:,:) :: ocean_mask
logical :: compute_ocean_mask

real :: wdata(ipts,jpts), hpi, del
real, dimension(size(ocean_mask,1),size(ocean_mask,2)) :: wdata_interpolated
!-----------------------------------------------------------------------

if ( read_pctwater_file(wdata) ) then

  if ( any(shape(ocean_mask) /= (/size(blon)-1,size(blat)-1/)) ) then
    call error_mesg('compute_ocean_mask','shape(ocean_mask) is not&
      & equal to (/size(blon)-1,size(blat)-1/))', FATAL)
  endif

  hpi = acos(0.0)
  del = 2.*hpi/real(jpts)
  call horiz_interp ( wdata, 0.0, -hpi, del, del, blon, blat, wdata_interpolated )
  where (wdata_interpolated > 50.)
    ocean_mask = .true.
  else where
    ocean_mask = .false.
  end where

  compute_ocean_mask = .true.
else
  compute_ocean_mask = .false.
endif

return
end function compute_ocean_mask
!#######################################################################

 function compute_stdev (blon, blat, stdev)

   real, intent(in),  dimension(:)   :: blon, blat
   real, intent(out), dimension(:,:) :: stdev
   logical :: compute_stdev

!  blon, blat = longitude and latitude in radians for grid box edges
!  stdev      = standard deviation of surface geop height (meters)

   real :: zsurf (size(stdev,1),size(stdev,2))
   real :: zdata (ipts,jpts),  hpi, del
!-----------------------------------------------------------------------

   if ( read_topog_file (zdata) ) then

      hpi = acos(0.0)
      del = 2.*hpi/real(jpts)
      call horiz_interp ( zdata, 0.0, -hpi, del, del, blon, blat, zsurf )
      zdata = zdata * zdata
      call horiz_interp ( zdata, 0.0, -hpi, del, del, blon, blat, stdev )
      stdev = stdev - zsurf*zsurf
      where (stdev > 0.0)
         stdev = sqrt ( stdev )
      elsewhere
         stdev = 0.0
      endwhere
      compute_stdev = .true.
   else
      compute_stdev = .false.
   endif

!-----------------------------------------------------------------------

 end function compute_stdev

!#######################################################################

 function read_topog_file (zdata)

   real, intent(out) :: zdata (ipts,jpts)
   logical           :: read_topog_file
   integer           :: unit

   read_topog_file = .false.

!  ------- reads topography (hires navy) data set ------

   if ( file_exist('INPUT/hires_navy_topography.data') ) then
       unit = open_file ('INPUT/hires_navy_topography.data', &
                         action='read', form='ieee32')
!      --- read and convert topog from cm to m ---
       read (unit) zdata
       zdata = zdata*0.01
       call close_file (unit)
       read_topog_file = .true.
   endif

 end function read_topog_file

!#######################################################################

function read_pctwater_file (wdata)

real, intent(out) :: wdata(ipts,jpts)
logical           :: read_pctwater_file
integer           :: unit

read_pctwater_file = .false.

!  ------- reads percent water data ------

if ( file_exist('INPUT/pctwater.data') ) then
  unit = open_file ('INPUT/pctwater.data', action='read', form='ieee32')
  read (unit) wdata
  call close_file (unit)
  read_pctwater_file = .true.
endif

return
end function read_pctwater_file
!#######################################################################

subroutine read_namelist

   integer :: unit, ierr, io
   real    :: hpi, dtr, ght

!-------------- read namelist --------------

   if ( file_exist('input.nml')) then
      unit = open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=topography_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'topography_nml')
      enddo
 10   call close_file (unit)
   endif

!-------- write version and namelist --------

   unit = open_file ('logfile.out', action='append')
   if (get_my_pe() == 0) then
       write (unit, '(/,80("="),/(a))') trim(version),trim(tag)
       write (unit, nml=topography_nml)
   endif
   call close_file (unit)

!--------- convert degrees to radians -----------

   hpi = acos(0.0)
   dtr = hpi/90.

   olon = olon*dtr
   olat = olat*dtr
   wlon = wlon*dtr
   wlat = wlat*dtr
   rlon = rlon*dtr
   rlat = rlat*dtr

end subroutine read_namelist

!#######################################################################

end module topography_mod
