module horiz_interp_mod

! <CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang </CONTACT>
! <CONTACT EMAIL="Bruce.Wyman@noaa.gov"> Bruce Wyman </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Performs spatial interpolation between grids.
! </OVERVIEW>

! <DESCRIPTION>
!     This module can interpolate data from rectangular/tripolar grid
!     to rectangular/tripolar grid. Three interpolation schems are used here.
!     when the source grid is tripolar grid, inverse of square distance weighted
!     scheme (spherical regrid) is used. When the source grid is rectangular 
!     longitude/latitude grid, any one of following three schemes can be applied: 
!     conservation scheme that conserves the area-weighed integral of the input field, 
!     bilinear interpolation that use the surround four source grid to interpolate onto
!     the destination grid, spherical regrid that use thes inverse of square distance
!     as weight. User can choose the interpolation method in the horiz_interp_init.
!     The default method is conservative scheme. When the source grid is tripolar grid,
!     only spherical regrid is allowed. When destination grid is tripolar, conservative 
!     scheme can not be used. 
!     There is an optional mask field for missing input data.
!     An optional output mask field may be used in conjunction with
!     the input mask to show where output data exists.
! </DESCRIPTION>

!-----------------------------------------------------------------------
!
!        Performs spatial interpolation between grids.
!
!-----------------------------------------------------------------------

use fms_mod,                    only: write_version_number
use mpp_mod,                    only: mpp_error, FATAL, stdout
use constants_mod,              only: pi
use horiz_interp_type_mod,      only: horiz_interp_type
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_init,  horiz_interp_conserve
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_end
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_init,  horiz_interp_bilinear
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_end
use horiz_interp_spherical_mod, only: horiz_interp_spherical_init,  horiz_interp_spherical
use horiz_interp_spherical_mod, only: horiz_interp_spherical_end

 implicit none
 private

!---- interfaces ----

 public   horiz_interp_type, horiz_interp, horiz_interp_init, &
          horiz_interp_end

! <INTERFACE NAME="horiz_interp_init">
!   <OVERVIEW>
!      Allocates space and initializes a derived-type variable
!      that contains pre-computed interpolation indices and weights.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Allocates space and initializes a derived-type variable
!      that contains pre-computed interpolation indices and weights
!      for improved performance of multiple interpolations between
!      the same grids. This routine does not need to be called if you
!      are doing a single grid-to-grid interpolation.
!   </DESCRIPTION>
!   <IN NAME="lon_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
!      Longitude (in radians) for source data grid. when lon_in is 1D, it is the longitude
!      edges and its value are for adjacent grid boxes and must increase 
!      in value. When lon_in is 2D, there are two cases: one is the longitude edges stored as
!      pairs for each grid box (when interp_method is "conservative"), the other is the longitude
!      of the center of each grid box. 
!   </IN>
!   <IN NAME="lat_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
!      Latitude (in radians) for source data grid. when lat_in is 1D, it is the latitude
!      edges and its value are for adjacent grid boxes and must increase 
!      in value. When lat_in is 2D, it is the longitude of the center of each grid box.
!      When interp_method is "conservative" or "bilinear", lon_in should be 1D.
!   </IN>
!   <IN NAME="lon_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
!      Longitude (in radians) for destination data grid. when lon_out is 1D, it is the longitude
!      edges and its value are for adjacent grid boxes and must increase 
!      in value. When lon_out is 2D, there are two cases: one is the longitude edges stored as
!      pairs for each grid box (when interp_method is "conservative"), the other is the longitude
!      of the center of each grid box (when interp_method is "bilinear"). 
!   </IN>
!   <IN NAME="lat_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
!      Latitude (in radians) for destination data grid. when lat_out is 1D, it is the latitude
!      edges and its value are for adjacent grid boxes and must increase 
!      in value. When lat_out is 2D, there are two cases: one is the latitude edges stored as
!      pairs for each grid box (when interp_method is "conservative"), the other is the latitude
!      of the center of each grid box (when interp_method is "bilinear").
!   </IN>
!   <IN NAME="verbose" TYPE="integer">
!      Integer flag that controls the amount of printed output.
!      verbose = 0, no output; = 1, min,max,means; = 2, still more
!   </IN>
!   <IN NAME="interp_method" TYPE="character(len=*)" > 
!      interpolation method, = "conservative", using conservation scheme,
!      = "bilinear", using bilinear interpolation, = "spherical",using spherical regrid.
!      when source grid is 1d, default value is "conservative"; when source grid is 2d,
!      default value is "spherical".
!   </IN>
!   <IN NAME = "src_modulo" >
!      Indicate the source data grid is cyclic or not.
!   </IN>
!   <OUT NAME="Interp" >
!      A derived-type variable containing indices and weights used for subsequent 
!      interpolations. To reinitialize this variable for a different grid-to-grid 
!      interpolation you must first use the "horiz_interp_end" interface.
!   </OUT>

 interface horiz_interp_init
    module procedure horiz_interp_init_1d     ! Source grid is 1d, destination grid is 1d
    module procedure horiz_interp_init_1d_src ! Source grid is 1d, destination grid is 2d
    module procedure horiz_interp_init_2d     ! Source grid is 2d, destination grid is 2d
    module procedure horiz_interp_init_1d_dst ! Source grid is 2d, destination grid is 1d
 end interface
! </INTERFACE>

! <INTERFACE NAME="horiz_interp">
!
!   <OVERVIEW>
!     Subroutine for performing the horizontal interpolation between two grids.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Subroutine for performing the horizontal interpolation between
!     two grids. There are two forms of this interface.
!     Form A requires first calling horiz_interp_init, while Form B
!     requires no initialization.
!   </DESCRIPTION>

!   <IN NAME="Interp" >
!     Derived-type variable containing interpolation indices and weights.
!     Returned by a previous call to horiz_interp_init.
!   </IN>
!   <IN NAME="data_in">
!      Input data on source grid.
!   </IN>
!   <IN NAME="verbose">
!      flag for the amount of print output.
!               verbose = 0, no output; = 1, min,max,means; = 2, still more
!   </IN>
!   <IN NAME="mask_in">
!      Input mask, must be the same size as the input data. The real value of
!      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points 
!      that should not be used or have missing data. It is Not needed for 
!      spherical regrid.
!   </IN>
!   <IN NAME="missing_value" >
!      Use the missing_value to indicate missing data.
!   </IN>
!   <IN NAME="missing_permit">
!      numbers of points allowed to miss for the bilinear interpolation. The value
!      should be between 0 and 3.
!   </IN>
!   <IN NAME="lon_in, lat_in" >
!      longitude and latitude (in radians) of source grid. More explanation can 
!      be found in the documentation of horiz_interp_init.
!   </IN>
!   <IN NAME="lon_out, lat_out" >
!      longitude and latitude (in radians) of destination grid. More explanation can 
!      be found in the documentation of horiz_interp_init.
!   </IN>
!   <OUT NAME="data_out">
!      Output data on destination grid.
!   </OUT>
!   <OUT NAME="mask_out">
!      Output mask that specifies whether data was computed.
!   </OUT>

!   <ERROR MSG="size of input array incorrect" STATUS="FATAL">
!      The input data array does not match the size of the input grid edges
!      specified. If you are using the initialization interface make sure you
!      have the correct grid size.
!   </ERROR>
!   <ERROR MSG="size of output array incorrect" STATUS="FATAL">
!      The output data array does not match the size of the input grid
!      edges specified. If you are using the initialization interface make
!      sure you have the correct grid size.
!   </ERROR>

 interface horiz_interp
    module procedure horiz_interp_base_2d
    module procedure horiz_interp_base_3d
    module procedure horiz_interp_solo_1d
    module procedure horiz_interp_solo_1d_src
    module procedure horiz_interp_solo_2d
    module procedure horiz_interp_solo_1d_dst
    module procedure horiz_interp_solo_old
 end interface
! </INTERFACE>

!---- derived-types ----


! parameter to determine interpolation method
 integer, parameter :: CONSERVE = 1
 integer, parameter :: BILINEAR = 2
 integer, parameter :: SPHERICA = 3

!-----------------------------------------------------------------------
 character(len=128) :: version = '$Id: horiz_interp.f90,v 12.0 2005/04/14 17:56:46 fms Exp $'
 character(len=128) :: tagname = '$Name: lima $'
 logical            :: do_vers = .true.
 logical            :: module_is_initialized = .FALSE.
!-----------------------------------------------------------------------

contains


!#######################################################################
!  <SUBROUTINE NAME="horiz_interp_init_1d" INTERFACE="horiz_interp_init">
!  <IN NAME="lon_in" TYPE="real" DIM="(:),(:,:)" UNITS="radians"></IN>
!  <IN NAME="lat_in" TYPE="real" DIM="(:),(:,:)"></IN>
!  <IN NAME="lon_out" TYPE="real" DIM="(:),(:,:)"></IN>
!  <IN NAME="lat_out" TYPE="real" DIM="(:),(:,:)"></IN>
!  <IN NAME="verbose" TYPE="integer, optional"></IN>
!  <IN NAME="interp_method" TYPE="character(len=*),optional"></IN>
!  <IN NAME="src_modulo" TYPE="logical, optional" > </IN>
!  <OUT NAME="Interp" TYPE="type(horiz_interp_type)"></OUT>

!<PUBLICROUTINE INTERFACE="horiz_interp_init">
  subroutine horiz_interp_init_1d (Interp, lon_in, lat_in, lon_out, lat_out,  &
       verbose, interp_method, num_nbrs, max_dist, src_modulo, &
       grid_at_center)
    !</PUBLICROUTINE>

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout)        :: Interp
    real, intent(in),  dimension(:)               :: lon_in , lat_in
    real, intent(in),  dimension(:)               :: lon_out, lat_out
    integer, intent(in),                 optional :: verbose
    character(len=*), intent(in),        optional :: interp_method
    integer, intent(in),                 optional :: num_nbrs
    real,    intent(in),                 optional :: max_dist
    logical, intent(in),                 optional :: src_modulo
    logical, intent(in),                 optional :: grid_at_center
    !-----------------------------------------------------------------------
    real, dimension(:,:), allocatable :: lon_src, lat_src, lon_dst, lat_dst
    real, dimension(:),   allocatable :: lon_src_1d, lat_src_1d
    integer                           :: i, j, nlon_in, nlat_in, nlon_out, nlat_out
    logical                           :: center
    character(len=40)                 :: method
    !-----------------------------------------------------------------------
    module_is_initialized = .true.
    !  write version number and tag name
    if (do_vers) then
       call write_version_number (version, tagname)
       do_vers = .false.
    endif

    method = 'conservative'
    if(present(interp_method)) method = interp_method

    select case (trim(method))
    case ("conservative")
       Interp%interp_method = CONSERVE
       nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
       allocate(lon_dst(nlon_out,2), lat_dst(nlat_out,2))
       do i=1,nlon_out
          lon_dst(i,1) = lon_out(i)
          lon_dst(i,2) = lon_out(i+1)
       enddo
       do j=1,nlat_out
          lat_dst(j,1) = lat_out(j)
          lat_dst(j,2) = lat_out(j+1)
       enddo
       call horiz_interp_conserve_init ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
            verbose)
       deallocate(lon_dst, lat_dst)
    case ("bilinear")
       Interp%interp_method = BILINEAR
       center = .false.
       if(present(grid_at_center) ) center = grid_at_center
       if(center) then
          nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i = 1, nlon_out
             lon_dst(i,:) = lon_out(i)
          enddo
          do j = 1, nlat_out
             lat_dst(:,j) = lat_out(j)
          enddo

          call horiz_interp_bilinear_init ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_dst, lat_dst)
       else
          nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i = 1, nlon_in
             lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
          enddo
          do j = 1, nlat_in
             lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
          enddo
          do i = 1, nlon_out
             lon_dst(i,:) = (lon_out(i) + lon_out(i+1)) * 0.5
          enddo
          do j = 1, nlat_out
             lat_dst(:,j) = (lat_out(j) + lat_out(j+1)) * 0.5
          enddo
          call horiz_interp_bilinear_init ( Interp, lon_src_1d, lat_src_1d, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst, lat_dst)
       endif
    case ("spherical")
       Interp%interp_method = SPHERICA
       nlon_in  = size(lon_in(:));   nlat_in  = size(lat_in(:))
       nlon_out  = size(lon_out(:)); nlat_out = size(lat_out(:))
       allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
       allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
       do i = 1, nlon_in
          lon_src(i,:) = lon_in(i)
       enddo
       do j = 1, nlat_in
          lat_src(:,j) = lat_in(j)
       enddo
       do i = 1, nlon_out
          lon_dst(i,:) = lon_out(i)
       enddo
       do j = 1, nlat_out
          lat_dst(:,j) = lat_out(j)
       enddo
       call horiz_interp_spherical_init ( Interp, lon_src, lat_src, lon_dst, lat_dst, &
            num_nbrs, max_dist, src_modulo)
       deallocate(lon_src, lat_src, lon_dst, lat_dst)
    case default
       call mpp_error(FATAL,'horiz_interp_mod: interp_method should be conservative, bilinear, spherical')
    end select

    !-----------------------------------------------------------------------

  end subroutine horiz_interp_init_1d
!  </SUBROUTINE>

!#######################################################################

 subroutine horiz_interp_init_1d_src (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                      verbose, interp_method, num_nbrs, max_dist, &
                                      src_modulo, grid_at_center )

   type(horiz_interp_type), intent(inout)        :: Interp
   real, intent(in),  dimension(:)               :: lon_in , lat_in
   real, intent(in),  dimension(:,:)             :: lon_out, lat_out
   integer, intent(in),                 optional :: verbose
   character(len=*), intent(in),        optional :: interp_method
   integer, intent(in),                 optional :: num_nbrs  ! minimum number of neighbors
   real,    intent(in),                 optional :: max_dist
   logical, intent(in),                 optional :: src_modulo
   logical, intent(in),                 optional :: grid_at_center

   real, dimension(:,:), allocatable :: lon_src, lat_src
   real, dimension(:),   allocatable :: lon_src_1d, lat_src_1d
   integer                           :: i, j, nlon_in, nlat_in
   character(len=40)                 :: method
   logical                           :: center
   !-----------------------------------------------------------------------
   module_is_initialized = .true.
   !  write version number and tag name
   if (do_vers) then
      call write_version_number (version, tagname)
      do_vers = .false.
   endif

   method = 'conservative'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      call horiz_interp_conserve_init ( Interp, lon_in, lat_in, lon_out, lat_out, &
           verbose )
   case ("bilinear")
      Interp%interp_method = BILINEAR
      center = .false.
      if(present(grid_at_center) ) center = grid_at_center
      if(center) then
         call horiz_interp_bilinear_init ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose, src_modulo )
      else
         nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
         allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
         do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
         enddo
         do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
         enddo
         call horiz_interp_bilinear_init ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
              verbose, src_modulo )
         deallocate(lon_src_1d,lat_src_1d)
      endif
   case ("spherical")
      Interp%interp_method = SPHERICA
      nlon_in  = size(lon_in(:));  nlat_in  = size(lat_in(:))
      allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
      do i = 1, nlon_in
         lon_src(i,:) = lon_in(i)
      enddo
      do j = 1, nlat_in
         lat_src(:,j) = lat_in(j)
      enddo
      call horiz_interp_spherical_init ( Interp, lon_src, lat_src, lon_out, lat_out, &
           num_nbrs, max_dist, src_modulo)
      deallocate(lon_src, lat_src)
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, spherical')
   end select

   !-----------------------------------------------------------------------

 end subroutine horiz_interp_init_1d_src

!#######################################################################

 subroutine horiz_interp_init_2d (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                  verbose, interp_method, num_nbrs, max_dist, src_modulo)
 type(horiz_interp_type), intent(inout)     :: Interp
 real, intent(in),  dimension(:,:)          :: lon_in , lat_in
 real, intent(in),  dimension(:,:)          :: lon_out, lat_out
 integer, intent(in),              optional :: verbose
 character(len=*), intent(in),     optional :: interp_method
 integer, intent(in),              optional :: num_nbrs
 real,    intent(in),              optional :: max_dist
 logical, intent(in),              optional :: src_modulo

 character(len=40) :: method
!-----------------------------------------------------------------------
   module_is_initialized = .true.   
!  write version number and tag name
   if (do_vers) then
      call write_version_number (version, tagname)
      do_vers = .false.
   endif

   method = 'bilinear'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("spherical") 
      Interp%interp_method = SPHERICA
      call horiz_interp_spherical_init ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                    num_nbrs, max_dist, src_modulo )
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_init ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                        verbose, src_modulo )
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select     

!-----------------------------------------------------------------------

 end subroutine horiz_interp_init_2d

!#######################################################################
 subroutine horiz_interp_init_1d_dst (Interp, lon_in, lat_in, lon_out, lat_out,   &
      verbose, interp_method, num_nbrs, max_dist, src_modulo)
   type(horiz_interp_type), intent(inout)     :: Interp
   real, intent(in),  dimension(:,:)          :: lon_in , lat_in
   real, intent(in),  dimension(:)            :: lon_out, lat_out
   integer, intent(in),              optional :: verbose
   character(len=*), intent(in),     optional :: interp_method
   integer, intent(in),              optional :: num_nbrs
   real,    intent(in),              optional :: max_dist
   logical, intent(in),              optional :: src_modulo

   character(len=40) :: method
   !-------------some local variables-----------------------------------------------
   integer                           :: i, j, nlon_out, nlat_out
   real, dimension(:,:), allocatable :: lon_dst, lat_dst
   !-----------------------------------------------------------------------
   module_is_initialized = .true.   
   !  write version number and tag name
   if (do_vers) then
      call write_version_number (version, tagname)
      do_vers = .false.
   endif

   method = 'bilinear'
   if(present(interp_method)) method = interp_method

   nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
   allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
   do i = 1, nlon_out
      lon_dst(i,:) = lon_out(i)
   enddo
   do j = 1, nlat_out
      lat_dst(:,j) = lat_out(j)
   enddo

   select case (trim(method))
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_init ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           verbose, src_modulo )
   case ("spherical")
      Interp%interp_method = SPHERICA
      call horiz_interp_spherical_init ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           num_nbrs, max_dist, src_modulo)
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select

   deallocate(lon_dst,lat_dst)

   !-----------------------------------------------------------------------

 end subroutine horiz_interp_init_1d_dst

!#######################################################################
! <SUBROUTINE NAME="horiz_interp_base_2d" INTERFACE="horiz_interp">
!   <IN NAME="Interp" TYPE="type(horiz_interp_type)"> </IN>
!   <IN NAME="data_in" TYPE="real" DIM="(:,:),(:,:,:)"> </IN>
!   <IN NAME="lon_in, lat_in" TYPE="real" DIM="(:),(:,:)"> </IN>
!   <IN NAME="lon_out, lat_out" TYPE="real" DIM="(:),(:,:)"> </IN>
!   <IN NAME="missing_value" TYPE="integer, optional" > </IN>
!   <IN NAME="missing_permit" TYPE="integer,optional" > </IN>
!   <IN NAME="verbose" TYPE="integer,optional"> </IN>
!   <IN NAME="mask_in" TYPE="real,optional" DIM="(:,:),(:,:,:)"> </IN>
!   <OUT NAME="data_out" TYPE="real" DIM="(:,:),(:,:,:)"> </OUT>
!   <OUT NAME="mask_out" TYPE="real,optional" DIM="(:,:),(:,:,:)"> </OUT>

!<PUBLICROUTINE INTERFACE="horiz_interp"> 
 subroutine horiz_interp_base_2d ( Interp, data_in, data_out, verbose, &
                                   mask_in, mask_out, missing_value, missing_permit )
!</PUBLICROUTINE>
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
      real, intent(in),                   optional :: missing_value
      integer, intent(in),                optional :: missing_permit
!-----------------------------------------------------------------------

   select case(Interp%interp_method)
   case(CONSERVE)
      call horiz_interp_conserve(Interp,data_in, data_out, verbose, mask_in, mask_out)
   case(BILINEAR)
      call horiz_interp_bilinear(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit )
   case(SPHERICA)
      call horiz_interp_spherical(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value )
   end select

   return

 end subroutine horiz_interp_base_2d
! </SUBROUTINE>

!#######################################################################

 subroutine horiz_interp_base_3d ( Interp, data_in, data_out, verbose, mask_in, mask_out, &
      missing_value, missing_permit  )
   !-----------------------------------------------------------------------
   !   overload of interface horiz_interp_base_2d
   !   uses 3d arrays for data and mask
   !   this allows for multiple interpolations with one call
   !-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in)           :: Interp
   real, intent(in),  dimension(:,:,:)            :: data_in
   real, intent(out), dimension(:,:,:)            :: data_out
   integer, intent(in),                  optional :: verbose
   real, intent(in),   dimension(:,:,:), optional :: mask_in
   real, intent(out),  dimension(:,:,:), optional :: mask_out
   real, intent(in),                     optional :: missing_value
   integer, intent(in),                  optional :: missing_permit
   !-----------------------------------------------------------------------
   integer :: n

   do n = 1, size(data_in,3)
      if (present(mask_in))then
         if(present(mask_out)) then
            call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), mask_out(:,:,n), &
                 missing_value, missing_permit )
         else
            call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
      else
         if(present(mask_out)) then
            call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_out=mask_out(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         else
            call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
     endif
   enddo
  
   return
!-----------------------------------------------------------------------
 end subroutine horiz_interp_base_3d

!#######################################################################
!<PUBLICROUTINE INTERFACE="horiz_interp"> 
 subroutine horiz_interp_solo_1d ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                   data_out, verbose, mask_in, mask_out,         &
                                   interp_method, missing_value, missing_permit, &
                                   num_nbrs, max_dist,src_modulo, grid_at_center  )              
!</PUBLICROUTINE>
!-----------------------------------------------------------------------
!   interpolates from a rectangular grid to rectangular grid.
!   interp_method can be the value conservative, bilinear or spherical.
!   horiz_interp_init don't need to be called before calling this routine.

!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:)   :: lon_in , lat_in
      real, intent(in),  dimension(:)   :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real, intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real, intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------

    call horiz_interp_init ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo, grid_at_center )

    call horiz_interp ( Interp, data_in, data_out, verbose,   &
                        mask_in, mask_out, missing_value, missing_permit )

    call horiz_interp_end ( Interp )
!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d

!#######################################################################

 subroutine horiz_interp_solo_1d_src ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                       data_out, verbose, mask_in, mask_out,         &
                                       interp_method, missing_value, missing_permit, &
                                       num_nbrs, max_dist, src_modulo, grid_at_center )
!-----------------------------------------------------------------------
!
!   interpolates from a uniformly spaced grid to any output grid.
!   interp_method can be the value "onservative","bilinear" or "spherical".
!   horiz_interp_init don't need to be called before calling this routine.
!
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:)   :: lon_in , lat_in
      real, intent(in),  dimension(:,:) :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real, intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real, intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------

    call horiz_interp_init ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo, grid_at_center )

    call horiz_interp ( Interp, data_in, data_out, verbose,   &
                        mask_in, mask_out, missing_value, missing_permit )

    call horiz_interp_end ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_src


!#######################################################################

 subroutine horiz_interp_solo_2d ( data_in, lon_in, lat_in, lon_out, lat_out, data_out, &
                                   verbose, mask_in, mask_out, interp_method, missing_value,&
                                   missing_permit, num_nbrs, max_dist, src_modulo  )
!-----------------------------------------------------------------------
!
!   interpolates from any grid to any grid. interp_method should be "spherical"
!   horiz_interp_init don't need to be called before calling this routine.
!
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:,:) :: lon_in , lat_in
      real, intent(in),  dimension(:,:) :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real, intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real, intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------

    call horiz_interp_init ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo )

    call horiz_interp ( Interp, data_in, data_out, verbose,   &
                        mask_in, mask_out, missing_value, missing_permit )

    call horiz_interp_end ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_2d

!#######################################################################

 subroutine horiz_interp_solo_1d_dst ( data_in, lon_in, lat_in, lon_out, lat_out, data_out,    &
                                       verbose, mask_in, mask_out,interp_method,missing_value, &
                                       missing_permit,  num_nbrs, max_dist, src_modulo)
!-----------------------------------------------------------------------
!
!   interpolates from any grid to rectangular longitude/latitude grid. 
!   interp_method should be "spherical".
!   horiz_interp_init don't need to be called before calling this routine.
!
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:,:) :: lon_in , lat_in
      real, intent(in),  dimension(:)   :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose  
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real, intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real, intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------

    call horiz_interp_init ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo )

    call horiz_interp ( Interp, data_in, data_out, verbose,   &
                        mask_in, mask_out, missing_value, missing_permit )

    call horiz_interp_end ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_dst

!#######################################################################

 subroutine horiz_interp_solo_old (data_in, wb, sb, dx, dy,  &
                                   lon_out, lat_out, data_out,  &
                                   verbose, mask_in, mask_out)

!-----------------------------------------------------------------------
!       Overloaded version of interface horiz_interp_solo_2
!
! input
!
!   data_in     Global input data stored from west to east (first dimension),
!               south to north (second dimension).  [real, dimension(:,:)]
!
!   wb          Longitude (in radians) that corresponds to western-most
!               boundary of grid box i=1 in array data_in.  [real]
!
!   sb          Latitude (in radians) that corresponds to southern-most
!               boundary of grid box j=1 in array data_in.  [real]
!
!   dx          Grid spacing (in radians) for the longitude axis (first
!               dimension) for the input data.  [real]
!
!   dy          Grid spacing (in radians) for the latitude axis (second
!               dimension) for the input data.  [real]
!
!   lon_out    The longitude edges (in radians) for output data grid boxes.
!               The values are for adjacent grid boxes and must increase in
!               value. If there are MLON grid boxes there must be MLON+1
!               edge values.  [real, dimension(:)]
!
!   lat_out    The latitude edges (in radians) for output data grid boxes.
!               The values are for adjacent grid boxes and may increase or
!               decrease in value. If there are NLAT grid boxes there must
!               be NLAT+1 edge values.  [real, dimension(:)]
!
! OUTPUT
!   data_out    Output data on the output grid defined by grid box
!               edges: blon_out and blat_out.  [real, dimension(:,:)]
!
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in)                  :: wb, sb, dx, dy
      real, intent(in),  dimension(:)   :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
!-----------------------------------------------------------------------
     real, dimension(size(data_in,1)+1)  :: blon_in
     real, dimension(size(data_in,2)+1)  :: blat_in
     integer :: i, j, nlon_in, nlat_in
     real    :: tpi
!-----------------------------------------------------------------------
 
   tpi = 2.*pi
   nlon_in = size(data_in,1)
   nlat_in = size(data_in,2)

   do i = 1, nlon_in+1
      blon_in(i) = wb + float(i-1)*dx
   enddo
      if (abs(blon_in(nlon_in+1)-blon_in(1)-tpi) < epsilon(blon_in)) &
              blon_in(nlon_in+1)=blon_in(1)+tpi

   do j = 2, nlat_in
      blat_in(j) = sb + float(j-1)*dy
   enddo
      blat_in(1)         = -0.5*pi
      blat_in(nlat_in+1) =  0.5*pi


   call horiz_interp_solo_1d (data_in, blon_in, blat_in,    &
                              lon_out, lat_out, data_out,   &
                              verbose, mask_in, mask_out    )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_old

!#######################################################################
! <SUBROUTINE NAME="horiz_interp_end">

!   <OVERVIEW>
!     Deallocates memory used by "horiz_interp_type" variables.
!       Must be called before reinitializing with horiz_interp_init.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Deallocates memory used by "horiz_interp_type" variables.
!     Must be called before reinitializing with horiz_interp_init.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call horiz_interp_end ( Interp )
!   </TEMPLATE>

!   <INOUT NAME="Interp" TYPE="horiz_interp_type">
!     A derived-type variable returned by previous call
!              to horiz_interp_init. The input variable must have
!              allocated arrays. The returned variable will contain
!              deallocated arrays.
!   </INOUT>

! </SUBROUTINE>

 subroutine horiz_interp_end ( Interp )

   type (horiz_interp_type), intent(inout) :: Interp

!-----------------------------------------------------------------------
!  releases space used by horiz_interp_type variables
!  must be called before re-initializing the same variable
!-----------------------------------------------------------------------
   select case(Interp % interp_method) 
   case (CONSERVE)
      call horiz_interp_conserve_end(Interp )
   case (BILINEAR)
      call horiz_interp_bilinear_end(Interp )
   case (SPHERICA)
      call horiz_interp_spherical_end(Interp )
   end select

!-----------------------------------------------------------------------

 end subroutine horiz_interp_end

 !#####################################################################

end module horiz_interp_mod

! <INFO>
!   <NOTE>             
!       Has not been checked with grids that do not cover the sphere.
!
!       Has not been checked with the optional mask arguments.
!
!       If a latitude or longitude index cannot be found the tolerance
!       used for making this determination may need to be increased.
!       This can be done by increasing the value of module variable
!       num_iters (default 4).
!   </NOTE>
!   <TESTPROGRAM>  
!     <PRE>
!       program test
!       use horiz_interp_mod
!       implicit none
!       integer, parameter :: nxi=177, nyi=91, nxo=133, nyo=77 ! resolution
!       real :: zi(nxi,nyi), zo(nxo,nyo)                       ! data
!       real :: xi(nxi+1), yi(nyi+1), xo(nxo+1), yo(nyo+1)     ! grid edges
!       real :: pi, tpi, hpi, dx, dy
!     
!       ! constants
!         hpi = acos(0.0)
!          pi = hpi*2.0
!         tpi = hpi*4.0
!     
!       ! grid setup: west to east, south to north
!         dx = tpi/real(nxi); call setaxis (0.,dx,xi);   xi(nxi+1) = xi(1)+tpi
!         dx = tpi/real(nxo); call setaxis (0.,dx,xo);   xo(nxo+1) = xo(1)+tpi
!         dy =  pi/real(nyi); call setaxis (-hpi,dy,yi); yi(nyi+1) = hpi
!         dy =  pi/real(nyo); call setaxis (-hpi,dy,yo); yo(nyo+1) = hpi
!     
!       ! random data on the input grid
!         call random_number (zi)
!     
!       ! interpolate (flipping y-axis)
!         call horiz_interp (zi(:,1:nyi:+1), xi, yi(1:nyi+1:+1), xo, yo(1:nyo+1:+1), zo, verbose=2)
!         call horiz_interp (zi(:,nyi:1:-1), xi, yi(nyi+1:1:-1), xo, yo(1:nyo+1:+1), zo, verbose=2)
!         call horiz_interp (zi(:,nyi:1:-1), xi, yi(nyi+1:1:-1), xo, yo(nyo+1:1:-1), zo, verbose=2)
!         call horiz_interp (zi(:,1:nyi:+1), xi, yi(1:nyi+1:+1), xo, yo(nyo+1:1:-1), zo, verbose=2)
!     
!       contains
!     ! set up a sequence of numbers
!         subroutine setaxis (xo,dx,x)
!         real, intent(in)  :: xo, dx
!         real, intent(out) :: x(:)
!         integer :: i
!           x(1) = xo
!           do i=2,size(x(:))
!             x(i) = x(i-1)+dx
!           enddo
!         end subroutine setaxis
!     
!       end program test
!     </PRE>
!   </TESTPROGRAM>
! </INFO>
