
module topography_mod

use gaussian_topog_mod, only: gaussian_topog_init, get_gaussian_topog
use   horiz_interp_mod, only: horiz_interp

use            fms_mod, only: file_exist, check_nml_error,               &
                              open_namelist_file, close_file, stdlog,    &
                              mpp_pe, mpp_root_pe, write_version_number, &
                              open_ieee32_file, error_mesg, FATAL

implicit none
private

public :: topography_init,                 &
          get_topog_mean, get_topog_stdev, &
          get_ocean_frac, get_ocean_mask,  &
          get_water_frac, get_water_mask,  &
          gaussian_topog_init, get_gaussian_topog

!-----------------------------------------------------------------------
!----- namelist ------

   character(len=128) :: topog_file = 'DATA/navy_topography.data', &
                         water_file = 'DATA/navy_pctwater.data'

   namelist /topography_nml/ topog_file, water_file

!-----------------------------------------------------------------------
! --- resolution of the topography data set ---

  integer :: unit
  integer :: ipts, jpts
  integer, parameter :: COMPUTE_STDEV = 123  ! use this flag to
                                             !   compute st dev

!-----------------------------------------------------------------------

 character(len=128) :: version = '$Id: topography.F90,v 1.5 2002/07/16 22:57:21 fms Exp $'
 character(len=128) :: tagname = '$Name: havana $'

 logical :: module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

 contains

!#######################################################################

   subroutine topography_init ()

     if ( module_is_initialized ) return

     call write_version_number (version,tagname)
     call read_namelist

     module_is_initialized = .TRUE.

   end subroutine topography_init

!#######################################################################

 function get_topog_mean (blon, blat, zmean)

   real, intent(in),  dimension(:)   :: blon, blat
   real, intent(out), dimension(:,:) :: zmean
   logical :: get_topog_mean

!  blon, blat = longitude and latitude in radians for grid box edges
!  zmean      = mean surface height (meters)

!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(zmean) /= (/size(blon)-1,size(blat)-1/)) ) &
        call error_mesg('get_topog_mean','shape(zmean) is not&
            & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   if ( open_topog_file(topog_file) ) then
       call interp_topog ( blon, blat, zmean )
       get_topog_mean = .true.
   else
       get_topog_mean = .false.
   endif

!-----------------------------------------------------------------------

 end function get_topog_mean

!#######################################################################
! returns standard deviation of "higher resolution" input data
! within the requested grid boxes

 function get_topog_stdev (blon, blat, stdev)

   real, intent(in),  dimension(:)   :: blon, blat
   real, intent(out), dimension(:,:) :: stdev
   logical :: get_topog_stdev

!  blon, blat = longitude and latitude in radians for grid box edges
!  stdev      = standard deviation of surface height (meters)

   real :: zsurf (size(stdev,1),size(stdev,2))
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(stdev) /= (/size(blon)-1,size(blat)-1/)) ) &
       call error_mesg('get_topog_stdev','shape(stdev) is not&
            & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   if ( open_topog_file(topog_file) ) then
       call interp_topog ( blon, blat, stdev, flag=COMPUTE_STDEV )
       get_topog_stdev = .true.
   else
       get_topog_stdev = .false.
   endif

!-----------------------------------------------------------------------

 end function get_topog_stdev

!#######################################################################
! returns fractional area covered by ocean

 function get_ocean_frac (blon, blat, ocean_frac)

 real, intent(in),  dimension(:)   :: blon, blat
 real, intent(out), dimension(:,:) :: ocean_frac
 logical :: get_ocean_frac

!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(ocean_frac) /= (/size(blon)-1,size(blat)-1/)) ) &
        call error_mesg('get_ocean_frac','shape(ocean_frac) is not&
                 & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   if ( open_topog_file(water_file) ) then
       call interp_water ( blon, blat, ocean_frac, do_ocean=.true. )
       get_ocean_frac = .true.
   else
       get_ocean_frac = .false.
   endif

!-----------------------------------------------------------------------

 end function get_ocean_frac

!#######################################################################

 function get_ocean_mask (blon, blat, ocean_mask)

 real   , intent(in),  dimension(:)   :: blon, blat
 logical, intent(out), dimension(:,:) :: ocean_mask
 logical :: get_ocean_mask

 real, dimension(size(ocean_mask,1),size(ocean_mask,2)) :: ocean_frac
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

 if ( get_ocean_frac(blon, blat, ocean_frac) ) then
   where (ocean_frac > 0.50)
     ocean_mask = .true.
   else where
     ocean_mask = .false.
   end where
   get_ocean_mask = .true.
 else
   get_ocean_mask = .false.
 endif

!-----------------------------------------------------------------------

 end function get_ocean_mask

!#######################################################################
! returns fractional area covered by water

 function get_water_frac (blon, blat, water_frac)

 real, intent(in),  dimension(:)   :: blon, blat
 real, intent(out), dimension(:,:) :: water_frac
 logical :: get_water_frac

!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(water_frac) /= (/size(blon)-1,size(blat)-1/)) ) &
        call error_mesg('get_water_frac','shape(water_frac) is not&
                 & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   if ( open_topog_file(water_file) ) then
       call interp_water ( blon, blat, water_frac )
       get_water_frac = .true.
   else
       get_water_frac = .false.
   endif

!-----------------------------------------------------------------------

 end function get_water_frac

!#######################################################################

 function get_water_mask (blon, blat, water_mask)

 real   , intent(in),  dimension(:)   :: blon, blat
 logical, intent(out), dimension(:,:) :: water_mask
 logical :: get_water_mask

 real, dimension(size(water_mask,1),size(water_mask,2)) :: water_frac
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

 if ( get_water_frac(blon, blat, water_frac) ) then
   where (water_frac > 0.50)
     water_mask = .true.
   else where
     water_mask = .false.
   end where
   get_water_mask = .true.
 else
   get_water_mask = .false.
 endif

!-----------------------------------------------------------------------

 end function get_water_mask

!#######################################################################
!##################   private interfaces below here   ##################
!#######################################################################

 function open_topog_file ( filename )
 character(len=*), intent(in) :: filename
 logical :: open_topog_file

  if ( file_exist(filename) ) then
       unit = open_ieee32_file (trim(filename), 'read')
       read (unit) ipts, jpts
       open_topog_file = .true.
  else
       open_topog_file = .false.
  endif

 end function open_topog_file

!#######################################################################

 subroutine interp_topog ( blon, blat, zout, flag )
 real   , intent(in)  :: blon(:), blat(:)
 real   , intent(out) :: zout(:,:)
 integer, intent(in), optional :: flag

 real :: xdat(ipts+1), ydat(jpts+1)
 real :: zdat(ipts,jpts)
 real :: zout2(size(zout,1),size(zout,2))

! note: ipts,jpts,unit are global

    read (unit) xdat, ydat    ! read lon/lat edges in radians
    read (unit) zdat          ! read land surface height in meters
    call close_file (unit)

    call horiz_interp ( zdat, xdat, ydat, blon, blat, zout )

! compute standard deviation if necessary
    if (present(flag)) then
       if (flag == COMPUTE_STDEV) then
           zdat = zdat*zdat
           call horiz_interp ( zdat, xdat, ydat, blon, blat, zout2 )
           zout = zout2 - zout*zout
           where (zout > 0.0)
             zout = sqrt ( zout )
           elsewhere
             zout = 0.0
           endwhere
       endif
    endif

 end subroutine interp_topog

!#######################################################################

 subroutine interp_water ( blon, blat, zout, do_ocean )
 real   , intent(in)  :: blon(:), blat(:)
 real   , intent(out) :: zout(:,:)
 logical, intent(in), optional :: do_ocean

 real :: xdat(ipts+1), ydat(jpts+1), zdat(ipts,jpts)

! note: ipts,jpts,unit are global

    read (unit) xdat, ydat    ! read lon/lat edges in radians
    read (unit) zdat          ! read fractional water
    call close_file (unit)

! only use designated ocean points
    if (present(do_ocean)) then
        if (do_ocean) call determine_ocean_points (zdat)
    endif

! interpolate onto output grid
    call horiz_interp ( zdat, xdat, ydat, blon, blat, zout )

 end subroutine interp_water

!#######################################################################

 subroutine determine_ocean_points ( pctwater )
 real, intent(inout) :: pctwater(:,:)
 logical :: ocean(size(pctwater,1),size(pctwater,2))
 integer :: i, j, m, n, im, ip, jm, jp, new 

 real :: ocean_pct_crit = .500

  ! resolution of the grid
    m = size(pctwater,1)
    n = size(pctwater,2)

  ! the 1/6 degree navy percent water data set
  ! designates ocean grid boxes as 100 percent water
  ! all other grid boxes have <= 99 percent water

  ! set a mask for ocean grid boxes
    ocean = (pctwater > .999) 
    new = count(ocean)

  ! set land grid boxes that have sufficient amount of water
  ! to ocean grid boxes when they are adjacent to ocean points
  ! iterate until there are no new ocean points
    do 
    if (new == 0) exit 
    new = 0 

       do j = 1, n
       do i = 1, m
          if (.not.ocean(i,j) .and. pctwater(i,j) > ocean_pct_crit) then
             im = i-1; ip = i+1; jm = j-1; jp = j+1
             if (im == 0)   im = m  
             if (ip == m+1) ip = 1
             if (jm == 0)   jm = 1
             if (jp == n+1) jp = n
           ! check the 8 grid boxes that surround this grid box
             if (ocean(im,j ) .or. ocean(ip,j ) .or. ocean(i ,jm) .or. ocean(i ,jp) .or. &
                 ocean(im,jm) .or. ocean(ip,jm) .or. ocean(ip,jp) .or. ocean(im,jp)) then
                 ocean(i,j) = .true.
                 new = new + 1
             endif
          endif
       enddo
       enddo
      !print *, 'new=',new

    enddo

  ! final step is to elimate water percentage if land
    where (.not.ocean) pctwater = 0.

 end subroutine determine_ocean_points

!#######################################################################
! reads the namelist file, write namelist to log file, 
! and initializes constants

subroutine read_namelist

   integer :: unit, ierr, io
   real    :: dtr, ght

!  read namelist

   if ( file_exist('input.nml')) then
      unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=topography_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'topography_nml')
      enddo
 10   call close_file (unit)
   endif

!  write version and namelist to log file

   if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=topography_nml)

end subroutine read_namelist

!#######################################################################

end module topography_mod
