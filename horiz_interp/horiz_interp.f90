module horiz_interp_mod
! <CONTACT EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </CONTACT>

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

use fms_mod, only:       error_mesg, FATAL, mpp_pe, mpp_root_pe,  &
                         mpp_npes, write_version_number
use mpp_mod, only:       mpp_send, mpp_recv, mpp_sync_self
use constants_mod, only: pi

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
!      Longitude (in radians) for source data grid. when lon_in is 1D, it is the longitude
!      edges and its value are for adjacent grid boxes and must increase 
!      in value. When lon_in is 2D, there are two cases: one is the longitude edges stored as
!      pairs for each grid box (when interp_method is "conservative"), the other is the longitude
!      of the center of each grid box (when interp_method is "bilinear"). 
!   </IN>
!   <IN NAME="lat_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
!      Latitude (in radians) for source data grid. when lat_in is 1D, it is the latitude
!      edges and its value are for adjacent grid boxes and must increase 
!      in value. When lat_in is 2D, there are two cases: one is the latitude edges stored as
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

!<PUBLICTYPE >
 type horiz_interp_type
   private
   real,    dimension(:,:), pointer   :: faci, facj  !weights for conservative scheme
   integer, dimension(:,:), pointer   :: ilon, jlat  !indices for conservative scheme
   real,    dimension(:,:), pointer   :: area_src    !area of the source grid
   real,    dimension(:,:), pointer   :: area_dst    !area of the destination grid
   real,    dimension(:,:,:), pointer :: wti,wtj     !weights for bilinear interpolation 
   integer, dimension(:,:,:), pointer :: i_lon, j_lat!indices for bilinear interpolation 
                                                     !and spherical regrid
   real,    dimension(:,:,:), pointer :: src_dist    !distance between destination grid and 
                                                     !neighbor source grid.
   logical, dimension(:,:), pointer   :: found_neighbors    !indicate whether destination grid 
                                                            !has some source grid around it.
   real                               :: max_src_dist
   integer                            :: nlon_src, nlat_src !size of source grid
   integer                            :: nlon_dst, nlat_dst !size of destination grid
   character(len=64)                  :: interp_method      !interpolation method.
                                                            !="conservative", conservative scheme
                                                            !="bilinear", bilinear interpolation
                                                            !="spherical", spherical regrid
 end type
!</PUBLICTYPE>
! public data to determine interpolation method
 real, parameter :: max_dist_default = 0.17  ! radians
 integer, parameter :: num_nbrs_default = 4
 real, parameter :: epsln=1.e-10, large=1.e20
!-----------------------------------------------------------------------
 character(len=128) :: version = '$Id: horiz_interp.f90,v 1.6 2003/04/09 21:17:09 fms Exp $'
 character(len=128) :: tagname = '$Name: inchon $'
 logical :: do_vers = .true.
 logical :: module_is_initialized = .FALSE.
 integer :: num_iters = 4
 integer, parameter :: stdout = 6
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
                                  verbose, interp_method, num_nbrs, max_dist, src_modulo)
!</PUBLICROUTINE>

!-----------------------------------------------------------------------
 type(horiz_interp_type), intent(out)          :: Interp
 real, intent(in),  dimension(:)               :: lon_in , lat_in
 real, intent(in),  dimension(:)               :: lon_out, lat_out
 integer, intent(in),                 optional :: verbose
 character(len=*), intent(in),        optional :: interp_method
 integer, intent(in),                 optional :: num_nbrs
 real,    intent(in),                 optional :: max_dist
 logical, intent(in),                 optional :: src_modulo
!-----------------------------------------------------------------------
   real, dimension(:,:), allocatable :: lon_src, lat_src, lon_dst, lat_dst
   real, dimension(:),   allocatable :: lon_src_1d, lat_src_1d
   integer :: i, j, nlon_in, nlat_in, nlon_out, nlat_out
!-----------------------------------------------------------------------

   module_is_initialized = .true.
!  write version number and tag name
   if (do_vers) then
      call write_version_number (version, tagname)
      do_vers = .false.
   endif

   Interp%interp_method = "conservative"
   if(present(interp_method)) Interp%interp_method = interp_method

   select case (trim(Interp%interp_method))
   case ("conservative")
      nlon_out = size(lon_out)-1; nlat_out = size(lat_out)-1
      allocate(lon_dst(nlon_out,2), lat_dst(nlat_out,2))
      do i=1,nlon_out
         lon_dst(i,1) = lon_out(i)
         lon_dst(i,2) = lon_out(i+1)
      enddo
      do j=1,nlat_out
         lat_dst(j,1) = lat_out(j)
         lat_dst(j,2) = lat_out(j+1)
      enddo
      call horiz_interp_init_cons ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
                                    verbose)
      deallocate(lon_dst, lat_dst)
   case ("bilinear")
      nlon_in  = size(lon_in)-1;  nlat_in  = size(lat_in)-1
      nlon_out = size(lon_out)-1; nlat_out = size(lat_out)-1
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
      call horiz_interp_init_bili ( Interp, lon_src_1d, lat_src_1d, lon_dst, lat_dst, &
                                    verbose, src_modulo )
      deallocate(lon_src_1d, lat_src_1d, lon_dst, lat_dst)
   case ("spherical")
      nlon_in  = size(lon_in);   nlat_in  = size(lat_in)
      nlon_out  = size(lon_out); nlat_out = size(lat_out)
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
      call horiz_interp_init_sphe ( Interp, lon_src, lat_src, lon_dst, lat_dst, &
                                    num_nbrs, max_dist, src_modulo )
      deallocate(lon_src, lat_src, lon_dst, lat_dst)
   case default
      call error_handler("interp_method should be conservative, bilinear, spherical")
   end select

!-----------------------------------------------------------------------

 end subroutine horiz_interp_init_1d
!  </SUBROUTINE>

!#######################################################################

 subroutine horiz_interp_init_1d_src (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                      verbose, interp_method, num_nbrs, max_dist, src_modulo )

 type(horiz_interp_type), intent(out)          :: Interp
 real, intent(in),  dimension(:)               :: lon_in , lat_in
 real, intent(in),  dimension(:,:)             :: lon_out, lat_out
 integer, intent(in),                 optional :: verbose
 character(len=*), intent(in),        optional :: interp_method
 integer, intent(in),                 optional :: num_nbrs
 real,    intent(in),                 optional :: max_dist
 logical, intent(in),                 optional :: src_modulo

   real, dimension(:,:), allocatable :: lon_src, lat_src
   real, dimension(:),   allocatable :: lon_src_1d, lat_src_1d
   integer :: i, j, nlon_in, nlat_in
!-----------------------------------------------------------------------

   module_is_initialized = .true.
!  write version number and tag name
   if (do_vers) then
      call write_version_number (version, tagname)
      do_vers = .false.
   endif

   Interp%interp_method = "conservative"
   if(present(interp_method)) Interp%interp_method = interp_method

   select case (trim(Interp%interp_method))
   case ("conservative")
      call horiz_interp_init_cons ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                    verbose )
   case ("bilinear")
      nlon_in  = size(lon_in)-1;  nlat_in  = size(lat_in)-1
      allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
      do i = 1, nlon_in
         lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
      enddo
      do j = 1, nlat_in
         lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
      enddo
      call horiz_interp_init_bili ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
                                    verbose, src_modulo )
      deallocate(lon_src_1d,lat_src_1d)
   case ("spherical")
      nlon_in  = size(lon_in);  nlat_in  = size(lat_in)
      allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
      do i = 1, nlon_in
         lon_src(i,:) = lon_in(i)
      enddo
      do j = 1, nlat_in
         lat_src(:,j) = lat_in(j)
      enddo
      call horiz_interp_init_sphe ( Interp, lon_src, lat_src, lon_out, lat_out, &
                                    num_nbrs, max_dist, src_modulo )
      deallocate(lon_src, lat_src)
   case default
      call error_handler("interp_method should be conservative, bilinear, spherical")
   end select

!-----------------------------------------------------------------------

 end subroutine horiz_interp_init_1d_src

!#######################################################################

 subroutine horiz_interp_init_2d (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                  verbose, interp_method, num_nbrs, max_dist, src_modulo )
 type(horiz_interp_type), intent(out)       :: Interp
 real, intent(in),  dimension(:,:)          :: lon_in , lat_in
 real, intent(in),  dimension(:,:)          :: lon_out, lat_out
 integer, intent(in),              optional :: verbose
 character(len=*), intent(in),     optional :: interp_method
 integer, intent(in),              optional :: num_nbrs
 real,    intent(in),              optional :: max_dist
 logical, intent(in),              optional :: src_modulo

!-----------------------------------------------------------------------

   module_is_initialized = .true.   
!  write version number and tag name
   if (do_vers) then
      call write_version_number (version, tagname)
      do_vers = .false.
   endif

   Interp%interp_method = "spherical"
   if(present(interp_method)) Interp%interp_method = interp_method

   select case (trim(Interp%interp_method))
   case ("spherical") 
      call horiz_interp_init_sphe ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                    num_nbrs, max_dist, src_modulo )
   case default
      call error_handler("when source grid are 2d, interp_method should be spherical")
   end select     

!-----------------------------------------------------------------------

 end subroutine horiz_interp_init_2d

!#######################################################################
 subroutine horiz_interp_init_1d_dst (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                      verbose, interp_method, num_nbrs, max_dist, src_modulo )
 type(horiz_interp_type), intent(out)       :: Interp
 real, intent(in),  dimension(:,:)          :: lon_in , lat_in
 real, intent(in),  dimension(:)            :: lon_out, lat_out
 integer, intent(in),              optional :: verbose
 character(len=*), intent(in),     optional :: interp_method
 integer, intent(in),              optional :: num_nbrs
 real,    intent(in),              optional :: max_dist
 logical, intent(in),              optional :: src_modulo

!-------------some local variables-----------------------------------------------
   integer :: i, j, nlon_out, nlat_out
   real, dimension(:,:), allocatable :: lon_dst, lat_dst
!-----------------------------------------------------------------------
   module_is_initialized = .true.   
!  write version number and tag name
   if (do_vers) then
      call write_version_number (version, tagname)
      do_vers = .false.
   endif

   Interp%interp_method = "spherical"
   if(present(interp_method)) Interp%interp_method = interp_method


   select case (trim(Interp%interp_method))
   case ("spherical")
      nlon_out = size(lon_out); nlat_out = size(lat_out)
      allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
      do i = 1, nlon_out
         lon_dst(i,:) = lon_out(i)
      enddo
      do j = 1, nlat_out
         lat_dst(:,j) = lat_out(j)
      enddo
      call horiz_interp_init_sphe ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
                                    num_nbrs, max_dist, src_modulo )
      deallocate(lon_dst,lat_dst)
   case default
      call error_handler("when source grid are 2d, interp_method should be 3")
   end select     

!-----------------------------------------------------------------------

 end subroutine horiz_interp_init_1d_dst

!#######################################################################
 subroutine horiz_interp_init_cons ( Interp, lon_in, lat_in,   &
                                     lon_out, lat_out, verbose)

 type(horiz_interp_type), intent(out) :: Interp
 real, intent(in),  dimension(:)      :: lon_in , lat_in
 real, intent(in),  dimension(:,:)    :: lon_out, lat_out
 integer, intent(in),     optional    :: verbose

!-----------------------------------------------------------------------
   real, dimension(size(lat_out,1),size(lat_out,2)) :: sph
   real, dimension(size(lat_in))    :: slat_in
   real, dimension(size(lon_in,1)-1) :: dlon_in
   real, dimension(size(lon_in,1)-1) :: dsph_in
   real, dimension(size(lon_out,1))  :: dlon_out
   real, dimension(size(lat_out,1))  :: dsph_out
   real    :: blon, fac, hpi, tpi, eps
   integer :: i, j, m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
              np, iverbose, m2, n2, iter
   logical :: s2n
   character(len=64) :: mesg
!-----------------------------------------------------------------------
   iverbose = 0;  if (present(verbose)) iverbose = verbose
   
   allocate ( Interp % facj (size(lat_out,1),size(lat_out,2)),      &
              Interp % jlat (size(lat_out,1),size(lat_out,2)),      &
              Interp % faci (size(lon_out,1),size(lon_out,2)),      &
              Interp % ilon (size(lon_out,1),size(lon_out,2)),      &
              Interp % area_src (size(lon_in)-1, size(lat_in)-1),   &
              Interp % area_dst (size(lon_out,1),size(lat_out,1)) )
!-----------------------------------------------------------------------
   hpi = 0.5*pi
   tpi = 4.*hpi

   nlon_in = size(lon_in)-1;  nlat_in = size(lat_in)-1

! check size of input arguments

   if ( size(lon_out,2) /=2 .or. size(lat_out,2) /= 2 )  &
   call error_handler ('when using conservative scheme, dimension 2 of lon_out and/or lat_out must be 2')

!-----------------------------------------------------------------------
!  --- set-up for input grid boxes ---

   do j = 1, nlat_in+1
       slat_in(j) = sin(lat_in(j))
   enddo

   do j = 1, nlat_in
       dsph_in(j) = abs(slat_in(j+1)-slat_in(j))
   enddo

   do i = 1,nlon_in
       dlon_in(i) = abs(lon_in(i+1)-lon_in(i))
   enddo

!  set south to north flag
   s2n = .true.
   if (lat_in(1) > lat_in(nlat_in+1)) s2n = .false.

!-----------------------------------------------------------------------
!  --- set-up for output grid boxes ---

   nlon_out = size(lon_out,1);  nlat_out = size(lat_out,1)

   do n = 1, nlat_out
      dsph_out(n) = abs(sin(lat_out(n,2))-sin(lat_out(n,1)))
   enddo

   do m = 1,nlon_out
      dlon_out(m) = abs(lon_out(m,2)-lon_out(m,1))
   enddo

   Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
   Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out
!***********************************************************************

!------ set up latitudinal indexing ------
!------ make sure output grid goes south to north ------

 do n = 1, nlat_out
    if (lat_out(n,1) < lat_out(n,2)) then
       sph(n,1) = sin(lat_out(n,1))
       sph(n,2) = sin(lat_out(n,2))
    else
       sph(n,1) = sin(lat_out(n,2))
       sph(n,2) = sin(lat_out(n,1))
    endif
 enddo

 Interp%jlat = 0
 do n2 = 1, 2         ! looping on grid box edges
 do n = 1, nlat_out   ! looping on output latitudes
     eps = 0.0
 do iter=1,num_iters
! find indices from input latitudes
 do j = 1, nlat_in
    if ( (s2n .and. (slat_in(j)-sph(n,n2)) <= eps .and.   &
                    (sph(n,n2)-slat_in(j+1)) <= eps) .or. &
         (.not.s2n .and. (slat_in(j+1)-sph(n,n2)) <= eps .and.  &
                         (sph(n,n2)-slat_in(j)) <= eps) ) then
         Interp%jlat(n,n2) = j
       ! weight with sin(lat) to exactly conserve area-integral
         fac = (sph(n,n2)-slat_in(j))/(slat_in(j+1)-slat_in(j))
         if (s2n) then
           if (n2 == 1) Interp%facj(n,n2) = 1.0 - fac
           if (n2 == 2) Interp%facj(n,n2) = fac
         else
           if (n2 == 1) Interp%facj(n,n2) = fac
           if (n2 == 2) Interp%facj(n,n2) = 1.0 - fac
         endif
         exit
     endif
 enddo
     if ( Interp%jlat(n,n2) /= 0 ) exit
   ! did not find this output grid edge in the input grid
   ! increase tolerance for multiple passes
     eps  = epsilon(sph)*real(10**iter)
 enddo
   ! no match
     if ( Interp%jlat(n,n2) == 0 ) then
          write (mesg,710) n,sph(n,n2)
      710 format (': n,sph=',i3,f14.7,40x)
          call error_handler ('no latitude index found'//trim(mesg))
     endif
 enddo
 enddo

!------ set up longitudinal indexing ------

    Interp%ilon = 0
    do m2 = 1, 2         ! looping on grid box edges
    do m = 1, nlon_out   ! looping on output longitudes
        blon = lon_out(m,m2)
        if ( blon < lon_in(1)         ) blon = blon + tpi
        if ( blon > lon_in(nlon_in+1) ) blon = blon - tpi
        eps = 0.0
    do iter=1,num_iters
  ! find indices from input longitudes
    do i = 1, nlon_in
        if ( (lon_in(i)-blon) <= eps .and. &
             (blon-lon_in(i+1)) <= eps ) then
             Interp%ilon(m,m2) = i
             fac = (blon-lon_in(i))/(lon_in(i+1)-lon_in(i))
             if (m2 == 1) Interp%faci(m,m2) = 1.0 - fac
             if (m2 == 2) Interp%faci(m,m2) = fac
             exit
        endif
    enddo
       if ( Interp%ilon(m,m2) /= 0 ) exit
     ! did not find this output grid edge in the input grid
     ! increase tolerance for multiple passes
       eps  = epsilon(blon)*real(10**iter)
    enddo
     ! no match
       if ( Interp%ilon(m,m2) == 0 ) then
           print *, 'lon_out,blon,blon_in,eps=',  &
              lon_out(m,m2),blon,lon_in(1),lon_in(nlon_in+1),eps
           call error_handler ('no longitude index found')
       endif
    enddo
    enddo

!  --- area of input grid boxes ---

   do j = 1,nlat_in
   do i = 1,nlon_in
       Interp%area_src(i,j) = dlon_in(i) * dsph_in(j)
   enddo
   enddo

!  --- area of output grid boxes ---

   do n = 1, nlat_out
   do m = 1, nlon_out
       Interp%area_dst(m,n) = dlon_out(m) * dsph_out(n)
   enddo
   enddo

!-----------------------------------------------------------------------
! this output may be quite lengthy and is not recommended
! when using more than one processor
  if (iverbose > 2) then
      write (*,801) (i,Interp%ilon(i,1),Interp%ilon(i,2),  &
                       Interp%faci(i,1),Interp%faci(i,2),i=1,nlon_out)
      write (*,802) (j,Interp%jlat(j,1),Interp%jlat(j,2),  &
                       Interp%facj(j,1),Interp%facj(j,2),j=1,nlat_out)
 801  format (/,2x,'i',4x,'is',5x,'ie',4x,'facis',4x,'facie',  &
              /,(i4,2i7,2f10.5))
 802  format (/,2x,'j',4x,'js',5x,'je',4x,'facjs',4x,'facje',  &
              /,(i4,2i7,2f10.5))
  endif
!-----------------------------------------------------------------------

 end subroutine horiz_interp_init_cons

!#######################################################################
 subroutine horiz_interp_init_bili ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                     verbose, src_modulo )

!-----------------------------------------------------------------------
 type(horiz_interp_type), intent(out) :: Interp
      real, intent(in),  dimension(:)   :: lon_in , lat_in
      real, intent(in),  dimension(:,:) :: lon_out, lat_out
   integer, intent(in),                   optional :: verbose
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
  integer :: nlon_in, nlat_in, nlon_out, nlat_out, n, m,     &
             i, j, ie, is, je, js, istart, jstart, iend, jend, &
             ln_err, lt_err, warns
  real    :: wtw, wte, wts, wtn, lon, lat, tpi, hpi,           &
             glt_min, glt_max, gln_min, gln_max, min_lon, max_lon
  logical           :: src_is_modulo


  hpi = 0.5*pi
  tpi = 4.0*hpi
  glt_min = hpi
  glt_max = -hpi
  gln_min = tpi
  gln_max = -tpi
  min_lon = 0.0
  max_lon = tpi
  ln_err = 0
  lt_err = 0
!-----------------------------------------------------------------------

   allocate ( Interp % wti (size(lon_out,1),size(lon_out,2),2),   &
              Interp % wtj (size(lon_out,1),size(lon_out,2),2),   &
              Interp % i_lon (size(lon_out,1),size(lon_out,2),2), &
              Interp % j_lat (size(lon_out,1),size(lon_out,2),2))
!-----------------------------------------------------------------------
   warns = 0
   if(present(verbose)) warns = verbose
   src_is_modulo = .true. 
   if (present(src_modulo)) src_is_modulo = src_modulo

   if(size(lon_out,1) /= size(lat_out,1) .or. size(lon_out,2) /= size(lat_out,2) ) &
   call error_handler ('when using bilinear interplation, the output grids should be geographical grids')    
   nlon_in = size(lon_in)  ; nlat_in = size(lat_in)
   nlon_out = size(lon_out, 1); nlat_out = size(lon_out, 2)
   Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
   Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out

   if(lon_in(nlon_in) - lon_in(1) .gt. tpi + epsln) &
      call error_handler("The range of source grid longitude should be no larger than tpi")

   if(lon_in(1) .lt. 0.0) then
      min_lon = lon_in(1)
      max_lon = lon_in(nlon_in)
   endif

!     find longitude points of data within interval [-360., 360.]
   istart = 1
   do i=2,nlon_in
      if (lon_in(i-1) .lt. -tpi .and. lon_in(i) .ge. -tpi) istart = i
   enddo 
   iend = nlon_in
   do i=2,nlon_in
      if (lon_in(i-1) .lt. tpi  .and. lon_in(i) .ge. tpi) iend = i
   enddo

!     find latitude points of data within interval [-90., 90.]
   jstart = 1
   do j=2,nlat_in
      if (lat_in(j-1) .lt. -1.0 * hpi .and. lat_in(j) .ge. -1.0*hpi) jstart = j
   enddo 
   jend = nlat_in
   do j=2,nlat_in
      if (lat_in(j-1) .lt. hpi .and. lat_in(j) .ge.hpi ) jend = j
   enddo

   do n = 1, nlat_out
   do m = 1, nlon_out
      lon = lon_out(m,n)
      lat = lat_out(m,n)

      if(lon .lt. min_lon) then
         lon = lon + tpi
      else if(lon .gt. max_lon) then
         lon = lon - tpi
      endif
      ! when the input grid is in not cyclic, the output grid should located inside
      ! the input grid
      if(.not. src_is_modulo) then
         if((lon .lt. lon_in(1)) .or. (lon .gt. lon_in(nlon_in))) &
         call error_handler('when input grid is not modulo, output grid should locate inside input grid')
      endif

      glt_min = min(lat,glt_min);  glt_max = max(lat,glt_max)
      gln_min = min(lon,gln_min);  gln_max = max(lon,gln_max)

      is = indp(lon, lon_in(istart:iend) ) + istart - 1
      if( lon_in(is) .gt. lon ) is = max(is - 1,istart)
      if( lon_in(is) .eq. lon .and. is .eq. nlon_in) is = max(is - 1,istart)
      ie = min(is+1,iend)
      if(lon_in(is) .ne. lon_in(ie) .and. lon_in(is) .le. lon) then
         wtw = ( lon_in(ie) - lon) / (lon_in(ie) - lon_in(is) )
      else
!     east or west of the last data value. this could be because a
!     cyclic condition is needed or the dataset is too small. 
         ln_err = 1
         ie = istart
         is = iend
         if (lon_in(ie) .ge. lon ) then
            wtw = (lon_in(ie) -lon)/(lon_in(ie)-lon_in(is)+tpi+epsln)
         else
            wtw = (lon_in(ie) -lon+tpi+epsln)/(lon_in(ie)-lon_in(is)+tpi+epsln)
         endif   
      endif
      wte = 1. - wtw

      js = indp(lat, lat_in(jstart:jend) ) + jstart - 1

      if( lat_in(js) .gt. lat ) js = max(js - 1, jstart)
      if( lat_in(js) .eq. lat .and. js .eq. jend) js = max(js - 1, jstart)
      je = min(js + 1, jend)

      if ( lat_in(js) .ne. lat_in(je) .and. lat_in(js) .le. lat) then
         wts = ( lat_in(je) - lat )/(lat_in(je)-lat_in(js))
      else
!     north or south of the last data value. this could be because a
!     pole is not included in the data set or the dataset is too small.
!     in either case extrapolate north or south
         lt_err = 1
         wts = 1.
      endif

      wtn = 1. - wts
   
      Interp % i_lon (m,n,1) = is; Interp % i_lon (m,n,2) = ie
      Interp % j_lat (m,n,1) = js; Interp % j_lat (m,n,2) = je
      Interp % wti   (m,n,1) = wtw
      Interp % wti   (m,n,2) = wte
      Interp % wtj   (m,n,1) = wts
      Interp % wtj   (m,n,2) = wtn
      
   enddo
   enddo
  if (ln_err .eq. 1 .and. warns > 0) then
    write (stdout,'(/,(1x,a))')                                      &
      '==> Warning: the geographic data set does not extend far   ', &
      '             enough east or west - a cyclic boundary       ', &
      '             condition was applied. check if appropriate   '
    write (stdout,'(/,(1x,a,2f8.4))')                                &
      '    data required between longitudes:', gln_min, gln_max,     &
      '      data set is between longitudes:', lon_in(istart), lon_in(iend)
    warns = warns - 1
  endif
!
  if (lt_err .eq. 1 .and. warns > 0) then
    write (stdout,'(/,(1x,a))')                                     &
      '==> Warning: the geographic data set does not extend far   ',&
      '             enough north or south - extrapolation from    ',&
      '             the nearest data was applied. this may create ',&
      '             artificial gradients near a geographic pole   ' 
    write (stdout,'(/,(1x,a,2f8.4))')                             &
      '    data required between latitudes:', glt_min, glt_max,   &
      '      data set is between latitudes:', lat_in(jstart), lat_in(jend)
    warns = warns - 1
  endif

  return

 end subroutine horiz_interp_init_bili

!####################################################################

 subroutine horiz_interp_init_sphe(Interp, lon_in,lat_in,lon_out,lat_out, &
                                   num_nbrs, max_dist, src_modulo)
 type(horiz_interp_type), intent(out) :: Interp
 real, intent(in), dimension(:,:) :: lon_in, lat_in, lon_out, lat_out
 integer, intent(in), optional   :: num_nbrs
 real, optional :: max_dist
 logical, intent(in), optional :: src_modulo

 !------local variables ---------------------------------------
 integer :: i, j, n, m, k
 integer :: map_dst_xsize, map_dst_ysize, map_src_xsize, map_src_ysize, &
            map_src_size, num_neighbors
 real    :: sum, max_src_dist
 integer, dimension(:), allocatable     :: ilon, jlat
 integer, dimension(:,:,:), allocatable :: map_src_add
 logical, dimension(size(lon_out,1), size(lon_out,2)) :: map_found_neighbors
 real, dimension(:,:,:), allocatable    :: map_src_dist
 real, dimension(:,:,:), allocatable    :: map_src_wgt

 !--------------------------------------------------------------
 map_dst_xsize=size(lon_out,1);map_dst_ysize=size(lon_out,2)
 map_src_xsize=size(lon_in,1);map_src_ysize=size(lon_in,2)
 map_src_size = map_src_xsize*map_src_ysize

 num_neighbors = num_nbrs_default
 if(present(num_nbrs)) num_neighbors = num_nbrs

 max_src_dist = max_dist_default
 if (PRESENT(max_dist)) max_src_dist = max_dist
 Interp%max_src_dist = max_src_dist

 allocate(map_src_add(map_dst_xsize,map_dst_ysize,num_neighbors),  &
          map_src_dist(map_dst_xsize,map_dst_ysize,num_neighbors), &
          map_src_wgt(map_dst_xsize,map_dst_ysize,num_neighbors),  &
          ilon(num_neighbors),jlat(num_neighbors)  )

! allocate memory to data type
 allocate(Interp%i_lon(map_dst_xsize,map_dst_ysize,num_neighbors),    &
          Interp%j_lat(map_dst_xsize,map_dst_ysize,num_neighbors),    &
          Interp%src_dist(map_dst_xsize,map_dst_ysize,num_neighbors), &
          Interp%found_neighbors(map_dst_xsize,map_dst_ysize)       )

 map_src_wgt = 0.0

 !using radial_search to find the nearest points and corresponding distance.
 call radial_search(lon_in, lat_in, lon_out, lat_out,map_src_add, map_src_dist, &
                    map_found_neighbors, num_nbrs,max_dist,src_modulo )

 do j=1,map_dst_ysize
 do i=1,map_dst_xsize
    sum=0.0
    do n=1,num_neighbors
       if(map_src_add(i,j,n) == 0) then
          jlat(n) = 0; ilon(n) = 0
       else
          jlat(n) = map_src_add(i,j,n)/map_src_xsize + 1
          ilon(n) = map_src_add(i,j,n) - (jlat(n)-1)*map_src_xsize
          if(ilon(n) == 0) then
             jlat(n) = jlat(n) - 1
             ilon(n) = map_src_xsize
          endif
       endif
    enddo
    Interp%i_lon(i,j,:) = ilon(:)
    Interp%j_lat(i,j,:) = jlat(:)
 enddo
 enddo

 Interp%src_dist(:,:,:) = map_src_dist(:,:,:)
 Interp%found_neighbors(:,:) = map_found_neighbors(:,:)
 Interp%nlon_src = map_src_xsize; Interp%nlat_src = map_src_ysize
 Interp%nlon_dst = map_dst_xsize; Interp%nlat_dst = map_dst_ysize

 deallocate(map_src_add, map_src_dist, map_src_wgt, ilon, jlat)
 return


 end subroutine horiz_interp_init_sphe

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

   select case(trim(Interp%interp_method))
   case("conservative")
      call horiz_interp_cons(Interp,data_in, data_out, verbose, mask_in, mask_out)
   case("bilinear")
      call horiz_interp_bili(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit )
   case("spherical")
      call horiz_interp_sphe(Interp,data_in, data_out, verbose, mask_in, mask_out, &
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
        call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                                    verbose, missing_value = missing_value,  &
                                    missing_permit = missing_permit )
     endif
   enddo
  
   return
!-----------------------------------------------------------------------
 end subroutine horiz_interp_base_3d

!########################################################################
 subroutine horiz_interp_cons ( Interp, data_in, data_out, verbose, &
                                   mask_in, mask_out)
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
!----------local variables----------------------------------------------------
      integer :: i, j, m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
                 miss_in, miss_out, unit, is, ie, js, je,   &
                 np, npass, iverbose, m2, n2, pe, root_pe
      real    :: cph, dsum, wsum, avg_in, min_in, max_in,   &
                 avg_out, min_out, max_out, blon, eps, asum,   &
                 dwtsum, wtsum, arsum, hpi, tpi, dtr, dsph, fis, fie, fjs, fje
      character(len=64) :: mesg
!-----------------------------------------------------------------------
   iverbose = 0;  if (present(verbose)) iverbose = verbose
   pe = mpp_pe();   root_pe = mpp_pe()

   hpi = 0.5*pi; tpi = 4.*hpi
   dtr = hpi/90.
   eps = epsilon(wtsum)

   nlon_in  = Interp%nlon_src;  nlat_in  = Interp%nlat_src
   nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

 !  --- error checking ---
   if (size(data_in,1) /= nlon_in .or. size(data_in,2) /= nlat_in) &
   call error_handler ('size of input array incorrect')

   if (size(data_out,1) /= nlon_out .or. size(data_out,2) /= nlat_out) &
   call error_handler ('size of output array incorrect')

   if (present(mask_in)) then
      if ( count(mask_in < -.0001 .or. mask_in > 1.0001) > 0 ) &
      call error_handler ('input mask not between 0,1')
   endif

!-----------------------------------------------------------------------
!---- loop through output grid boxes ----

   data_out = 0.0
   do n = 1, nlat_out
    ! latitude window
    ! setup ascending latitude indices and weights
      if (Interp%jlat(n,1) <= Interp%jlat(n,2)) then
          js = Interp%jlat(n,1); je = Interp%jlat(n,2)
         fjs = Interp%facj(n,1); fje = Interp%facj(n,2)
      else
          js = Interp%jlat(n,2); je = Interp%jlat(n,1)
         fjs = Interp%facj(n,2); fje = Interp%facj(n,1)
      endif

   do m = 1, nlon_out
    ! longitude window
       is = Interp%ilon(m,1); ie = Interp%ilon(m,2)
      fis = Interp%faci(m,1); fie = Interp%faci(m,2)
      npass = 1
      dwtsum = 0.
       wtsum = 0.
       arsum = 0.

    ! wrap-around on input grid
    ! sum using 2 passes (pass 1: end of input grid)
      if ( ie < is ) then
           ie = nlon_in
          fie = 1.0
          npass = 2
      endif

      do np = 1, npass
       ! pass 2: beginning of input grid
         if ( np == 2 ) then
              is = 1
             fis = 1.0
              ie = Interp%ilon(m,2)
             fie = Interp%faci(m,2)
         endif

       ! summing data*weight and weight for single grid point
         if (present(mask_in)) then
            call data_sum_conserve ( data_in(is:ie,js:je), Interp%area_src(is:ie,js:je), &
                                     fis, fie, fjs,fje, dwtsum, wtsum, arsum, mask_in(is:ie,js:je)  )
         else
            call data_sum_conserve ( data_in(is:ie,js:je), Interp%area_src(is:ie,js:je), &
                                     fis, fie, fjs,fje,  dwtsum, wtsum, arsum    )
         endif
      enddo

      if (wtsum > eps) then
         data_out(m,n) = dwtsum/wtsum
         if (present(mask_out)) mask_out(m,n) = wtsum/arsum
      else
         data_out(m,n) = 0.
         if (present(mask_out)) mask_out(m,n) = 0.0
      endif

   enddo
   enddo

!***********************************************************************
! compute statistics: minimum, maximum, and mean
!-----------------------------------------------------------------------

 if (iverbose > 0) then

 ! compute statistics of input data

   call stats_type1(data_in, Interp%area_src, asum, dsum, wsum, min_in, max_in, miss_in, mask_in)
   ! diagnostic messages
   ! on the root_pe, we can calculate the global mean, minimum and maximum.
   if(pe == root_pe) then
      if (wsum > 0.0) then
         avg_in=dsum/wsum
      else
         print *, 'horiz_interp stats: input area equals zero '
         avg_in=0.0
      endif
     if (iverbose > 1) print '(2f16.11)', 'global sum area_in  = ',  asum, wsum
   endif
 
 ! compute statistics of output data
   call stats_type1(data_out, Interp%area_dst, asum, dsum, wsum, min_out, max_out, miss_out, mask_out)
   ! diagnostic messages
   if(pe == root_pe) then
      if (wsum > 0.0) then
         avg_out=dsum/wsum
      else
         print *, 'horiz_interp stats: output area equals zero '
         avg_out=0.0
      endif
      if (iverbose > 1) print '(2f16.11)', 'global sum area_out = ',  asum, wsum
   endif
    !---- output statistics ----
    ! the global mean, min and max are calculated on the root pe.
   if(pe == root_pe) then
      write (*,900)
      write (*,901)  min_in ,max_in ,avg_in
      if (present(mask_in))  write (*,903)  miss_in
      write (*,902)  min_out,max_out,avg_out
      if (present(mask_out)) write (*,903)  miss_out
   endif

 900  format (/,1x,10('-'),' output from horiz_interp ',10('-'))
 901  format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
 902  format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
 903  format ('          number of missing points = ',i6)

 endif

!-----------------------------------------------------------------------
 end subroutine horiz_interp_cons


!#######################################################################
 subroutine horiz_interp_bili ( Interp, data_in, data_out, verbose, mask_in,mask_out, &
                                missing_value, missing_permit)
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in)        :: Interp
   real, intent(in),  dimension(:,:)           :: data_in
   real, intent(out), dimension(:,:)           :: data_out
   integer, intent(in),               optional :: verbose
   real, intent(in), dimension(:,:),  optional :: mask_in
   real, intent(out), dimension(:,:), optional :: mask_out
   real, intent(in),                  optional :: missing_value
   integer, intent(in),               optional :: missing_permit
!-----------------------------------------------------------------------
      integer :: nlon_in, nlat_in, nlon_out, nlat_out, n, m,         &
                 is, ie, js, je, iverbose, max_missing, num_missing, &
                 miss_in, miss_out, pe, root_pe
      real    :: dwtsum, wtsum, min_in, max_in, avg_in, &
                 min_out, max_out, avg_out, wtw, wte, wts, wtn
      real,dimension(4) :: data, mask
      nlon_in  = Interp%nlon_src;  nlat_in  = Interp%nlat_src
      nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

   pe = mpp_pe(); root_pe = mpp_root_pe()
   iverbose = 0;  if (present(verbose)) iverbose = verbose

   max_missing = 0; if(present(missing_permit)) max_missing = missing_permit
   if(max_missing .gt. 3 .or. max_missing .lt. 0) call error_handler ('missing_permit should be between 0 and 3')

   if (size(data_in,1) /= nlon_in .or. size(data_in,2) /= nlat_in) &
   call error_handler ('size of input array incorrect')

   if (size(data_out,1) /= nlon_out .or. size(data_out,2) /= nlat_out) &
   call error_handler ('size of output array incorrect')

   data_out = 0.0

   do n = 1, nlat_out
   do m = 1, nlon_out

      is = Interp % i_lon (m,n,1); ie = Interp % i_lon (m,n,2)
      js = Interp % j_lat (m,n,1); je = Interp % j_lat (m,n,2)
      wtw = Interp % wti   (m,n,1)
      wte = Interp % wti   (m,n,2)
      wts = Interp % wtj   (m,n,1)
      wtn = Interp % wtj   (m,n,2)
      data(1) = data_in(is,js); data(2) = data_in(ie,js)
      data(3) = data_in(ie,je); data(4) = data_in(is,je)
      if(present(mask_in)) then
         mask(1) = mask_in(is,js); mask(2) = mask_in(ie,js)
         mask(3) = mask_in(ie,je); mask(4) = mask_in(is,je) 
      endif

      if(present(mask_in)) then  
         call data_sum_bilinear ( data,wtw,wte, wts,wtn, dwtsum, wtsum, num_missing, mask, missing_value)
      else
         call data_sum_bilinear ( data, wtw,wte, wts,wtn, dwtsum, wtsum, num_missing, missing_value = missing_value)
         !--- this will make sure it can reproduce old results.
         if(.not. present(missing_value)) wtsum = 1.0
      endif

     if(present(mask_out)) mask_out(m,n) = wtsum

     if(num_missing .gt. max_missing .or. wtsum .lt. epsln) then
        if(present(missing_value)) then
           data_out(m,n) = missing_value
        else
           data_out(m,n) = 0.0
        endif
        if(present(mask_out)) mask_out(m,n) = 0.0      
     else
        data_out(m,n) = dwtsum/wtsum
     endif
   enddo
   enddo

!***********************************************************************
! compute statistics: minimum, maximum, and mean
!-----------------------------------------------------------------------

 if (iverbose > 0) then

 ! compute statistics of input data

   call stats_type2 (data_in, min_in, max_in, avg_in, miss_in, missing_value, mask_in)

 ! compute statistics of output data
   call stats_type2 (data_out, min_out, max_out, avg_out, miss_out, missing_value, mask_out)

    !---- output statistics ----
    ! root_pe have the information of global mean, min and max
   if(pe == root_pe) then
      write (*,900)
      write (*,901)  min_in ,max_in, avg_in
      if (present(mask_in))  write (*,903)  miss_in
      write (*,902)  min_out,max_out,avg_out
      if (present(mask_out)) write (*,903)  miss_out
   endif
 900  format (/,1x,10('-'),' output from horiz_interp ',10('-'))
 901  format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
 902  format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
 903  format ('          number of missing points = ',i6)

 endif


   return

 end subroutine horiz_interp_bili

!#######################################################################

 subroutine horiz_interp_sphe ( Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value)
   type (horiz_interp_type), intent(in)        :: Interp
   real, intent(in),  dimension(:,:)           :: data_in
   real, intent(out), dimension(:,:)           :: data_out
   integer, intent(in),               optional :: verbose
   real, intent(in), dimension(:,:),  optional :: mask_in
   real, intent(out), dimension(:,:), optional :: mask_out
   real, intent(in),                  optional :: missing_value
   
   !--- some local variables ----------------------------------------
   real, dimension(:,:),   allocatable :: mask_src, mask_dst
   real, dimension(:,:,:), allocatable :: wt
   integer :: nlon_in, nlat_in, nlon_out, nlat_out, num_neighbors, &
              m, n, k, i, j, miss_in, miss_out, pe, root_pe,       &
              i1, i2, j1, j2, iverbose
   real    :: min_in, max_in, avg_in, min_out, max_out, avg_out, sum
   !-----------------------------------------------------------------

   pe = mpp_pe(); root_pe = mpp_root_pe()
   iverbose = 0;  if (present(verbose)) iverbose = verbose

   nlon_in  = Interp%nlon_src; nlat_in  = Interp%nlat_src
   nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst   
   num_neighbors = size(Interp%src_dist,3)

   if(size(data_in,1) .ne. nlon_in .or. size(data_in,2) .ne. nlat_in ) &
            call error_handler("size of input array incorrect")

   if(size(data_out,1) .ne. nlon_out .or. size(data_out,2) .ne. nlat_out ) &
            call error_handler("size of output array incorrect")

   allocate(mask_src(nlon_in, nlat_in), mask_dst(nlon_out, nlat_out), &
            wt(nlon_out, nlat_out,num_neighbors)  )
   mask_src = 1.0; mask_dst = 1.0
   if(present(mask_in)) mask_src = mask_in

   !first calculate destination mask and weight
   where(.not. Interp%found_neighbors ) mask_dst = 0.0

   do m=1,nlon_out
   do n=1,nlat_out
      ! neighbors are sorted nearest to farthest
      ! check nearest to see if it is a land point
      i1 = Interp%i_lon(m,n,1); j1 = Interp%j_lat(m,n,1)

      if( i1 .eq. 0 .and. j1 .eq. 0) then
         mask_dst(m,n) = 0.0
      else if (mask_src(i1,j1) .lt. 0.5) then
         mask_dst(m,n) = 0.0
      endif
      
      if(num_neighbors .gt. 1 ) then
         i2 = Interp%i_lon(m,n,2); j2 = Interp%j_lat(m,n,2)
      ! compare first 2 nearest neighbors -- if they are nearly
      ! equidistant then use this mask for robustness
         if(i2 .ne. 0 .and. j2 .ne. 0 .and. abs(Interp%src_dist(m,n,2)-Interp%src_dist(m,n,1)) .lt. epsln) then
            if((mask_src(i1,j1) .lt. 0.5))  mask_dst(m,n) = 0.0
         endif
      endif

      sum=0.0
      do k=1, num_neighbors
        if (Interp%src_dist(m,n,k) <= epsln) then
           wt(m,n,k) = large
           sum = sum + large
        else if(Interp%src_dist(m,n,k) <= Interp%max_src_dist .and.    &
               mask_src(Interp%i_lon(m,n,k),Interp%j_lat(m,n,k)) .ge. 0.5) then
           wt(m,n,k) = 1.0/Interp%src_dist(m,n,k)
           sum = sum+wt(m,n,k)
        else
           wt(m,n,k) = 0.0
        endif
     enddo
     if (sum > epsln) then
        do k = 1, num_neighbors
           wt(m,n,k) = wt(m,n,k)/sum
        enddo
     else
        mask_dst(m,n) = 0.0
     endif

   enddo
   enddo

   data_out = 0.0
   do m=1,nlon_out
   do n=1,nlat_out
      if(mask_dst(m,n) .gt. 0.5) then
         do k=1, num_neighbors
            i = Interp%i_lon(m,n,k)
            j = Interp%j_lat(m,n,k)
            data_out(m,n) = data_out(m,n)+data_in(i,j)*wt(m,n,k)
         enddo
      else
         if(present(missing_value)) then
            data_out(m,n) = missing_value
         else
            data_out(m,n) = 0.0
         endif
      endif
   enddo
   enddo

   if(present(mask_out)) mask_out = mask_dst

!***********************************************************************
! compute statistics: minimum, maximum, and mean
!-----------------------------------------------------------------------

 if (iverbose > 0) then

 ! compute statistics of input data

   call stats_type2 (data_in, min_in, max_in, avg_in, miss_in, missing_value, mask=mask_src)

 ! compute statistics of output data
   call stats_type2 (data_out, min_out, max_out, avg_out, miss_out, missing_value, mask=mask_dst)

    !---- output statistics ----
    ! root_pe have the information of global mean, min and max
   if(pe == root_pe) then
      write (*,900)
      write (*,901)  min_in ,max_in, avg_in
      if (present(mask_in))  write (*,903)  miss_in
      write (*,902)  min_out,max_out,avg_out
      if (present(mask_out)) write (*,903)  miss_out
   endif
 900  format (/,1x,10('-'),' output from horiz_interp ',10('-'))
 901  format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
 902  format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
 903  format ('          number of missing points = ',i6)

 endif

   !release the memory
   deallocate(mask_src, mask_dst, wt)

   return

 end subroutine horiz_interp_sphe


!#######################################################################
!<PUBLICROUTINE INTERFACE="horiz_interp"> 
 subroutine horiz_interp_solo_1d ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                   data_out, verbose, mask_in, mask_out,         &
                                   interp_method, missing_value, missing_permit, &
                                   num_nbrs, max_dist,src_modulo  )                                
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
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------

    call horiz_interp_init ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo )

    call horiz_interp ( Interp, data_in, data_out, verbose,   &
                        mask_in, mask_out, missing_value, missing_permit )

    call horiz_interp_end ( Interp )
!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d

!#######################################################################

 subroutine horiz_interp_solo_1d_src ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                       data_out, verbose, mask_in, mask_out,         &
                                       interp_method, missing_value, missing_permit, &
                                       num_nbrs, max_dist, src_modulo )
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
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------

    call horiz_interp_init ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo )

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
   select case(trim(Interp % interp_method)) 
   case ("conservative")
      deallocate ( Interp%area_src , Interp%area_dst , &
                Interp%facj, Interp%jlat,            &
                Interp%faci, Interp%ilon  )
   case ("bilinear")
      deallocate (Interp%wti,   Interp%wtj,              &
                  Interp%i_lon,    Interp%j_lat )
   case ("spherical")
      deallocate ( Interp%found_neighbors, Interp%src_dist, &
                   Interp%i_lon,    Interp%j_lat )
   end select

!-----------------------------------------------------------------------

 end subroutine horiz_interp_end

!#######################################################################
!---This statistics is for conservative scheme
 subroutine stats_type1 ( dat, area, asum, dsum, wsum, low, high, miss, mask )
 real,    intent(in)  :: dat(:,:), area(:,:)
 real,    intent(out) :: asum, dsum, wsum, low, high
 integer, intent(out) :: miss
 real,    intent(in), optional :: mask(:,:)

 integer :: pe, root_pe, npes, p, buffer_int(1)
 real    :: buffer_real(5)

   pe      = mpp_pe()
   root_pe = mpp_root_pe()
   npes    = mpp_npes()
  ! sum data, data*area; and find min,max on each pe.

   if (present(mask)) then
      asum = sum(area(:,:))
      dsum = sum(area(:,:)*dat(:,:)*mask(:,:))
      wsum = sum(area(:,:)*mask(:,:))
      miss = count(mask(:,:) <= 0.5)
      low  = minval(dat(:,:),mask=mask(:,:) > 0.5)
      high = maxval(dat(:,:),mask=mask(:,:) > 0.5)
   else
      asum = sum(area(:,:))
      dsum = sum(area(:,:)*dat(:,:))
      wsum = sum(area(:,:))
      miss = 0
      low  = minval(dat(:,:))
      high = maxval(dat(:,:))
   endif
  
   ! other pe send local min, max, avg to the root pe and 
   ! root pe receive these information
   call mpp_sync_self()   
   if(pe == root_pe) then
      do p = 1, npes - 1
         call mpp_recv(buffer_real,5,root_pe+p)
         asum = asum + buffer_real(1)
         dsum = dsum + buffer_real(2)
         wsum = wsum + buffer_real(3)
         low  = min(low, buffer_real(4))
         high = max(high, buffer_real(5))
         call mpp_recv(buffer_int,1,root_pe+p)
         miss = miss + buffer_int(1)
      enddo
   else
      buffer_real(1) = asum
      buffer_real(2) = dsum
      buffer_real(3) = wsum
      buffer_real(4) = low
      buffer_real(5) = high
      call mpp_send(buffer_real,5,root_pe)
      buffer_int(1) = miss
      call mpp_send(buffer_int,1,root_pe)
   endif
 end subroutine stats_type1
!#######################################################################
!---This statistics is for bilinear interpolation and spherical regrid.
 subroutine stats_type2 ( dat, low, high, avg, miss, missing_value, mask )
 real,    intent(in)  :: dat(:,:)
 real,    intent(out) :: low, high, avg
 integer, intent(out) :: miss
 real, intent(in), optional :: missing_value
 real,    intent(in), optional :: mask(:,:)

 real :: dsum, npts, buffer_real(3)
 integer :: pe, root_pe, npes, p, buffer_int(2)

   dsum = 0.0
   pe = mpp_pe()
   root_pe = mpp_root_pe()
   npes = mpp_npes()
   miss = 0

   if (present(missing_value)) then
      miss = count(dat(:,:) == missing_value)
      low  = minval(dat(:,:), dat(:,:) /= missing_value)
      high = maxval(dat(:,:), dat(:,:) /= missing_value)
      dsum = sum(dat(:,:), dat(:,:) /= missing_value)
   else if(present(mask)) then
      miss = count(mask(:,:) <= 0.5)
      low  = minval(dat(:,:),mask=mask(:,:) > 0.5)
      high = maxval(dat(:,:),mask=mask(:,:) > 0.5)
      dsum = sum(dat(:,:), mask=mask(:,:) > 0.5)
   else
      miss = 0
      low  = minval(dat(:,:))
      high = maxval(dat(:,:))
      dsum = sum(dat(:,:))
   endif
   avg = 0.0
   
   npts = size(dat) - miss
   if(pe == root_pe) then
      do p = 1, npes - 1  ! root_pe receive data from other pe
         call mpp_recv(buffer_real,3, p+root_pe)
         dsum = dsum + buffer_real(1)
         low  = min(low, buffer_real(2))
         high = max(high, buffer_real(3))
         call mpp_recv(buffer_int, 2, p+root_pe)
         miss = miss + buffer_int(1)
         npts = npts + buffer_int(2)
      enddo         
      if(npts == 0) then
         print*, 'Warning: no points is valid'
      else
         avg = dsum/real(npts)
      endif
    else   ! other pe send data to the root_pe.
      buffer_real(1) = dsum
      buffer_real(2) = low
      buffer_real(3) = high
      call mpp_send(buffer_real,3,root_pe)
      buffer_int(1) = miss
      buffer_int(2) = npts
      call mpp_send(buffer_int, 2, root_pe)
    endif

    return

 end subroutine stats_type2

!#######################################################################

 subroutine data_sum_conserve ( data, area, facis, facie, facjs, facje,  &
                       dwtsum, wtsum, arsum, mask )

!  sums up the data and weights for a single output grid box
!-----------------------------------------------------------------------
   real, intent(in), dimension(:,:) :: data, area
   real, intent(in)                 :: facis, facie, facjs, facje
   real, intent(inout)              :: dwtsum, wtsum, arsum
   real, intent(in), optional       :: mask(:,:)

!  fac__ = fractional portion of each boundary grid box included
!          in the integral
!  dwtsum = sum(data*area*mask)
!  wtsum  = sum(area*mask)
!  arsum  = sum(area)
!-----------------------------------------------------------------------
   real, dimension(size(area,1),size(area,2)) :: wt
   real    :: asum
   integer :: id, jd
!-----------------------------------------------------------------------

   id=size(area,1); jd=size(area,2)

   wt=area
   wt( 1,:)=wt( 1,:)*facis
   wt(id,:)=wt(id,:)*facie
   wt(:, 1)=wt(:, 1)*facjs
   wt(:,jd)=wt(:,jd)*facje

    asum = sum(wt)
   arsum = arsum + asum

   if (present(mask)) then
      wt = wt * mask
      dwtsum = dwtsum + sum(wt*data)
       wtsum =  wtsum + sum(wt)
   else
      dwtsum = dwtsum + sum(wt*data)
       wtsum =  wtsum + asum
   endif

!-----------------------------------------------------------------------

 end subroutine data_sum_conserve

!#######################################################################
 subroutine data_sum_bilinear ( data, wtw, wte, wts, wtn, dwtsum, wtsum, &
                                num_missing, mask, missing_value  )

!  sums up the data and weights for a single output grid box
!-----------------------------------------------------------------------
   real, intent(in), dimension(:) :: data
   real, intent(in)               :: wtw, wte, wts, wtn
   real, intent(inout)            :: dwtsum, wtsum
   integer, intent(out)           :: num_missing
   real, intent(in), dimension(:), optional     :: mask
   real, intent(in),               optional     :: missing_value

!  dwtsum = sum(data*wt*mask)
!  wtsum  = sum(wt*mask)

   integer :: i
   real :: wt1, wt2
!-----------------------------------------------------------------------

   if(size(data) .ne. 4 )      &
            call error_handler('size of data and weight is not correct')
   num_missing = 0
   dwtsum = 0.0
   wtsum  = 0.0

   do i = 1,size(data)
      if(present(missing_value) .and. data(i) .eq. missing_value) then
         num_missing = num_missing + 1
      else
         select case(i)
         case(1)
            wt1= wtw; wt2 = wts
         case(2)
            wt1= wte; wt2 = wts
         case(3)
            wt1= wte; wt2 = wtn
         case(4)
            wt1= wtw; wt2 = wtn
         end select

         if(present(mask)) then
            dwtsum = dwtsum + data(i)*mask(i)*wt1*wt2
            wtsum  = wtsum + mask(i)* wt1 * wt2
         else
            dwtsum = dwtsum + data(i) * wt1 * wt2
            wtsum  = wtsum + wt1 * wt2
         endif
      endif
   end do
!-----------------------------------------------------------------------

 end subroutine data_sum_bilinear

!#######################################################################

 subroutine error_handler ( message )
 character(len=*), intent(in) :: message

   call error_mesg ('horiz_interp_mod', message, FATAL)

 end subroutine error_handler

!#######################################################################


 function indp (value, array)
 integer                        :: indp
 real, dimension(:), intent(in) :: array
 real, intent(in)               :: value
!
!=======================================================================
!
!     indp = index of nearest data point within "array" corresponding to
!            "value".

!     inputs:
!     value  = arbitrary data...same units as elements in "array"
!     array  = array of data points  (must be monotonically increasing)

!     output:
!     indp =  index of nearest data point to "value"
!             if "value" is outside the domain of "array" then indp = 1
!             or "ia" depending on whether array(1) or array(ia) is
!             closest to "value"
!=======================================================================
!
 integer i, ii, ia
 logical keep_going
!
 ia = size(array)
  do i=2,ia
    if (array(i) .lt. array(i-1)) then
      write (stdout,*) &
     ' => Error: array must be monotonically increasing in "indp"' , &
     '           when searching for nearest element to value=',value
      write (stdout,*) '           array(i) < array(i-1) for i=',i 
      write (stdout,*) '           array(i) for i=1..ia follows:'
      call abort()
    endif
  enddo
  if (value .lt. array(1) .or. value .gt. array(ia)) then
    if (value .lt. array(1))  indp = 1
    if (value .gt. array(ia)) indp = ia
  else
    i=1
    keep_going = .true.
    do while (i .le. ia .and. keep_going)
      i = i+1
      if (value .le. array(i)) then
        indp = i
        if (array(i)-value .gt. value-array(i-1)) indp = i-1
        keep_going = .false.
      endif
    enddo
  endif
  return
 end function indp

!#######################################################################

 subroutine radial_search(x_src,y_src,x_dst,y_dst,map_src_add, map_src_dist, &
                          map_found_neighbors, num_nbrs,max_dist,src_modulo)

 ! x_dst,y_dst = destination grid lon,lat
 ! x_src,y_src = source grid lon,lat

 real, intent(in), dimension(:,:)       :: x_src, y_src, x_dst, y_dst
 integer, intent(out), dimension(:,:,:) :: map_src_add
 real,    intent(out), dimension(:,:,:) :: map_src_dist
 logical, intent(out), dimension(:,:)   :: map_found_neighbors
 integer, intent(in),          optional :: num_nbrs
 real, intent(in),             optional :: max_dist
 logical, intent(in),          optional :: src_modulo
 
 !-------local variables-----------------------------------------
 integer :: i, j, n, m, nx_dst, ny_dst, nx_src, ny_src, k, num
 integer :: bound, bound_start, bound_end, i0, j0, i_left, i_right
 real, dimension(:), allocatable :: theta_src,phi_src
 real, dimension(size(x_dst,1),size(x_dst,2)) :: theta_dst, phi_dst
 real :: min_theta_dst, max_theta_dst, min_phi_dst, max_phi_dst
 real :: min_theta_src, max_theta_src, min_phi_src, max_phi_src
 logical :: continue_search, found_neighbors, continue_radial_search, result, &
            src_is_modulo
 real :: d, sum, nearest,res, max_src_dist, tpi, hpi
 integer :: step, i_nearest, len, step_size, num_neighbors
 integer :: map_dst_xsize, map_dst_ysize, map_src_xsize, map_src_ysize, &
            map_dst_size, map_src_size
 !---------------------------------------------------------------

 tpi = 2.0*PI; hpi = 0.5*PI

 map_dst_xsize=size(x_dst,1);map_dst_ysize=size(x_dst,2)
 map_src_xsize=size(x_src,1);map_src_ysize=size(x_src,2)
 map_dst_size = map_dst_xsize*map_dst_ysize
 map_src_size = map_src_xsize*map_src_ysize

 num_neighbors = num_nbrs_default
 if(present(num_nbrs)) num_neighbors = num_nbrs
 if (num_neighbors <= 0) call error_handler('num_neighbors must be > 0') 

 max_src_dist = max_dist_default
 if (PRESENT(max_dist)) max_src_dist = max_dist

 src_is_modulo = .true.
 if (PRESENT(src_modulo)) src_is_modulo = src_modulo

 allocate(theta_src(map_src_size), phi_src(map_src_size))

 if (map_dst_size /= size(y_dst,1)*size(y_dst,2) .or. map_src_size /= size(y_src,1)*size(y_src,2)) &
      call error_handler( '=> grids not conformable')

 theta_src = reshape(x_src,(/map_src_size/))
 phi_src = reshape(y_src,(/map_src_size/))
 theta_dst(:,:) = x_dst(:,:)
 phi_dst(:,:) = y_dst(:,:)

 map_src_add = 0
 map_src_dist = large

 min_theta_dst=tpi;max_theta_dst=0.0;min_phi_dst=pi;max_phi_dst=-pi
 min_theta_src=tpi;max_theta_src=0.0;min_phi_src=pi;max_phi_src=-pi

 where(theta_dst<0.0)  theta_dst = theta_dst+tpi
 where(theta_dst>tpi)  theta_dst = theta_dst-tpi
 where(theta_src<0.0)  theta_src = theta_src+tpi
 where(theta_src>tpi)  theta_src = theta_src-tpi

 where(phi_dst < -hpi) phi_dst = -hpi
 where(phi_dst > hpi)  phi_dst =  hpi
 where(phi_src < -hpi) phi_src = -hpi
 where(phi_src > hpi)  phi_src =  hpi    

 do j=1,map_dst_ysize
 do i=1,map_dst_xsize
    min_theta_dst = min(min_theta_dst,theta_dst(i,j))
    max_theta_dst = max(max_theta_dst,theta_dst(i,j))
    min_phi_dst = min(min_phi_dst,phi_dst(i,j))
    max_phi_dst = max(max_phi_dst,phi_dst(i,j))
 enddo
 enddo

 do i=1,map_src_size
    min_theta_src = min(min_theta_src,theta_src(i))
    max_theta_src = max(max_theta_src,theta_src(i))
    min_phi_src = min(min_phi_src,phi_src(i))
    max_phi_src = max(max_phi_src,phi_src(i))
 enddo

 if (min_phi_dst < min_phi_src) print *, '=> WARNING:  latitute of dest grid exceeds src'
 if (max_phi_dst > max_phi_src) print *, '=> WARNING:  latitute of dest grid exceeds src'
 if (min_theta_dst < min_theta_src) print *, '=> WARNING : longitude of dest grid exceeds src'
 if (max_theta_dst > max_theta_src) print *, '=> WARNING : longitude of dest grid exceeds src'

 do j=1,map_dst_ysize
 do i=1,map_dst_xsize
    found_neighbors=.false.
    continue_search=.true.
    step = 1
    step_size = 1
    nearest = 1.e3
    do while (continue_search .and. step_size > 0)
       do while (step <= map_src_size .and. continue_search)
          ! count land points as nearest neighbors
          d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(step),phi_src(step))
          if (d <= max_src_dist) then
             found_neighbors = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                                     step,d,max_src_dist)
             if (found_neighbors) then
                n = 0
                i0 = mod(step,map_src_xsize)
                if (i0 == 0) i0 = map_src_xsize
                res = float(step)/float(map_src_xsize)
                j0 = ceiling(res)
                continue_radial_search = .true.
                do while (continue_radial_search)
                   n = n+1 ! radial counter
                   ! left boundary 
                   i_left = i0-n
                   if (i_left <= 0) then
                      if (src_is_modulo) then
                         i_left = map_src_xsize + i_left
                      else
                         i_left = 1
                      endif
                   endif
                   bound_start = max(j0-n-1,0)*map_src_xsize + i_left
                   bound_end = min(j0+n-1,map_src_ysize-1)*map_src_xsize + i_left
                   if (bound_end > map_src_size) call error_handler(" ")
                   bound = bound_start
                   continue_radial_search = .false.
                   do while (bound <= bound_end)
                      d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                      result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                                     bound,d,max_src_dist)
                      bound = bound + map_src_xsize
                      if (result) continue_radial_search = .true.
                      if (result) found_neighbors = .true.
                   enddo
                   ! right boundary 
                   i_right = i0+n
                   if (i_right > map_src_xsize) then
                      if (src_is_modulo) then
                         i_right = i_right - map_src_xsize
                      else
                         i_right = map_src_xsize
                      endif
                   endif
                   bound_start = max(j0-n-1,0)*map_src_xsize + i_right
                   bound_end = min(j0+n-1,map_src_ysize-1)*map_src_xsize + i_right
                   bound = bound_start
                   if (bound_end > map_src_size) call error_handler(' ')
                   do while (bound <= bound_end)
                      d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                      result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                                     bound,d,max_src_dist)
                      bound = bound + map_src_xsize
                      if (result) continue_radial_search = .true.
                      if (result) found_neighbors = .true.
                   enddo
                   ! bottom boundary 
                   bound_start = max(j0-n-1,0)*map_src_xsize + i_left 
                   bound_end =  max(j0-n-1,0)*map_src_xsize  + i_right 
                   if (bound_start > bound_end) then
                      bound_start = max(j0-n-1,0)*map_src_xsize + 1
                      bound_end = max(j0-n,1)*map_src_xsize
                   endif
                   bound = bound_start
                   if (bound_end > map_src_size) call error_handler(' ')
                   do while (bound <= bound_end)
                      d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                      result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                                     bound,d,max_src_dist)
                      bound = bound + 1
                      if (result) continue_radial_search = .true.
                      if (result) found_neighbors = .true.
                   enddo
                   ! top boundary 
                   bound_start = min(j0+n-1,map_src_ysize-1)*map_src_xsize + i_left
                   bound_end =  min(j0+n-1,map_src_ysize-1)*map_src_xsize + i_right
                   if (bound_start > bound_end) then
                      bound_start = min(j0+n-1,map_src_ysize-1)*map_src_xsize + 1
                      bound_end = min(j0+n,map_src_ysize-1)*map_src_xsize
                   endif
                   bound = bound_start
                   if (bound_end > map_src_size) call error_handler(' ')
                   do while (bound <= bound_end)
                      d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                      result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                                     bound,d,max_src_dist)
                      bound = bound + 1
                      if (result) continue_radial_search = .true.
                      if (result) found_neighbors = .true.
                   enddo
                enddo
                continue_search = .false. ! stop looking
             endif
          endif
          step=step+step_size
       enddo ! search loop
       step = 1
       step_size = step_size/2
    enddo
    map_found_neighbors(i,j) = found_neighbors
 enddo 
 enddo

 deallocate(theta_src,phi_src)

return

end subroutine radial_search

!########################################################################
function spherical_distance(theta1,phi1,theta2,phi2)

real, intent(in) :: theta1, phi1, theta2, phi2
real :: spherical_distance

real :: r1(3), r2(3), cross(3), s, dot, ang

! this is a simple, enough way to calculate distance on the sphere
! first, construct cartesian vectors r1 and r2
! then calculate the cross-product which is proportional to the area
! between the 2 vectors.  The angular distance is arcsin of the 
! distancealong the sphere
!
! theta is longitude and phi is latitude
!


r1(1) = cos(theta1)*cos(phi1);r1(2)=sin(theta1)*cos(phi1);r1(3)=sin(phi1)
r2(1) = cos(theta2)*cos(phi2);r2(2)=sin(theta2)*cos(phi2);r2(3)=sin(phi2)

cross(1) = r1(2)*r2(3)-r1(3)*r2(2)
cross(2) = r1(3)*r2(1)-r1(1)*r2(3)
cross(3) = r1(1)*r2(2)-r1(2)*r2(1)

s = sqrt(cross(1)**2.+cross(2)**2.+cross(3)**2.)

s = min(s,1.0-epsln)

dot = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)

if (dot > 0) then
    ang = asin(s)
else if (dot < 0) then
    ang = pi - asin(s)
else
    ang = pi/2.
endif

spherical_distance = abs(ang) ! in radians

return

end function spherical_distance

!#####################################################################
 function update_dest_neighbors(map_src_add, map_src_dist, src_add,d, max_src_dist)
 
 integer, intent(inout), dimension(:) :: map_src_add
 real, intent(inout),    dimension(:) :: map_src_dist
 integer, intent(in)                  :: src_add
 real, intent(in)                     :: d
 real, intent(in)                     :: max_src_dist

 logical :: update_dest_neighbors

 integer :: n,m, num_neighbors

 update_dest_neighbors = .false.
 num_neighbors = size(map_src_add) 

if (d .le. max_src_dist) then
   NLOOP : do n=1,num_neighbors
      DIST_CHK : if (d .lt. map_src_dist(n)) then
         if (n > 1 .and. src_add == map_src_add(n-1)) exit NLOOP
         do m=num_neighbors,n+1,-1
            map_src_add(m) = map_src_add(m-1)
            map_src_dist(m) = map_src_dist(m-1)
         enddo
         map_src_add(n) = src_add
         map_src_dist(n) = d
         update_dest_neighbors = .true.
         exit NLOOP ! n loop
      endif DIST_CHK
   end do NLOOP
endif

return

end function update_dest_neighbors

!#######################################################################

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
!           do i=2,size(x)
!             x(i) = x(i-1)+dx
!           enddo
!         end subroutine setaxis
!     
!       end program test
!     </PRE>
!   </TESTPROGRAM>
! </INFO>
