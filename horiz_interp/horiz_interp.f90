
module horiz_interp_mod

!-----------------------------------------------------------------------
!
!        Performs spatial interpolation between rectangular
!                latitude/longitude grids.
!
!-----------------------------------------------------------------------

use fms_mod, only:       error_mesg, FATAL,    &
                         mpp_pe, write_version_number
use constants_mod, only: pi

 implicit none
 private

!---- interfaces ----

 public   horiz_interp_type, horiz_interp, horiz_interp_init, &
          horiz_interp_end

 interface horiz_interp
    module procedure horiz_interp_base_2d
    module procedure horiz_interp_base_3d
    module procedure horiz_interp_solo_1
    module procedure horiz_interp_solo_2
    module procedure horiz_interp_solo_old
 end interface

 interface horiz_interp_init
    module procedure horiz_interp_init_1
    module procedure horiz_interp_init_2
 end interface

!---- derived-types ----

 type horiz_interp_type
   private
   real,    dimension(:), pointer :: dlon_in, dlon_out, & ! delta long
                                     dsph_in, dsph_out    ! delta sin(lat)
   real,    dimension(:,:), pointer :: faci, facj   ! weights
   integer, dimension(:,:), pointer :: ilon, jlat   ! indices
   real,    dimension(:,:,:), pointer :: wti,wtj    ! weights for global output grids
                                                    ! for bilinear interplation
   integer, dimension(:,:,:), pointer :: i_lon, j_lat! !indices for global output grids
                                                       ! for bilinear_interplation
   integer                           :: interp_method  ! determine the interpolation method,
                                                       ! =1 represents conservative scheme
                                                       ! =2 represents bilinear intrepolation
 end type

!-----------------------------------------------------------------------
 character(len=128) :: version = '$Id: horiz_interp.f90,v 1.5 2002/07/16 22:55:50 fms Exp $'
 character(len=128) :: tagname = '$Name: havana $'
 logical :: do_vers = .true.
 logical :: module_is_initialized = .FALSE.
 integer :: num_iters = 4
 integer, parameter :: stdout = 6
!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine horiz_interp_base_2d ( Interp, data_in, data_out, &
                                   verbose, mask_in, mask_out )

!-----------------------------------------------------------------------
!
!   interpolates from a uniformly spaced grid to any output grid
!   using an area weighing scheme to conserve the global integral
!   of the field being interpolated or bilinear interpolation.
!
!  input:
!  -----
!     Interp     Derived-type variable containing interpolation
!                indices and weights. Returned by a previous call
!                to horiz_interp_init.
!
!     data_in    Input data on input grid defined by grid box edges
!                initialized in variable Interp.
!                  [real, dimension(:,:)]
!
!  output:
!  ------
!     data_out   Output data on output grid defined by grid box edges
!                initialized in variable Interp.
!                  [real, dimension(:,:)]
!
!  optional
!  --------
!     verbose   flag for the amount of print output (integer scalar)
!                 =0 no output; =1 min,max,means; =2 still more
!
!     mask_in   input mask;  =0.0 for data points that should not
!               be used or have no data; has the same size as data_in
!     mask_out  output mask that specifies whether data was computed
!
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
!-----------------------------------------------------------------------
      real, dimension(size(Interp%dlon_in),  &
                      size(Interp%dsph_in))  :: area_in
      real, dimension(size(data_out,1), &
                      size(data_out,2)) :: area_out

      integer :: i, j, m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
                 miss_in, miss_out, unit, is, ie, js, je,   &
                 np, npass, iverbose, m2, n2, pe

      real :: cph, dsum, wsum, avg_in, min_in, max_in,   &
              avg_out, min_out, max_out, blon, eps,    &
              dwtsum, wtsum, arsum, hpi, tpi, dtr, dsph, fis, fie, fjs, fje

      character(len=64) :: mesg

!-----------------------------------------------------------------------
   iverbose = 0;  if (present(verbose)) iverbose = verbose

   if(Interp % interp_method == 2) then   ! =2 means using bilinear interpolation
      if(present(mask_in)) then
         if(present(mask_out)) then
            call new_horiz_interp_base_2d(Interp,data_in, data_out, mask_in, mask_out)
         else 
            call new_horiz_interp_base_2d(Interp,data_in, data_out, mask_in )
         endif
      else
         if(present(mask_out)) then
            call new_horiz_interp_base_2d(Interp,data_in, data_out, mask_out= mask_out)
         else
            call new_horiz_interp_base_2d(Interp,data_in, data_out)
         endif
      endif

      return
   endif

   pe = mpp_pe()

   hpi = 0.5*pi
   tpi = 4.*hpi
   dtr = hpi/90.
   eps = epsilon(wtsum)


!  --- area of input grid boxes ---

   nlon_in = size(Interp%dlon_in); nlat_in = size(Interp%dsph_in)

   do j = 1,nlat_in
   do i = 1,nlon_in
       area_in(i,j) = Interp % dlon_in(i) * Interp % dsph_in(j)
   enddo
   enddo

!  --- area of output grid boxes ---

   nlon_out = size(Interp%dlon_out); nlat_out = size(Interp%dsph_out)
 
   do n = 1, nlat_out
   do m = 1,nlon_out
       area_out(m,n) = Interp % dlon_out(m) * Interp % dsph_out(n)
   enddo
   enddo

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

   do n = 1, nlat_out

    ! latitude window
    ! setup ascending latitude indices and weights
      if (Interp%jlat(n,1) <= Interp%jlat(n,2)) then
          js = Interp%jlat(n,1)
          je = Interp%jlat(n,2)
         fjs = Interp%facj(n,1)
         fje = Interp%facj(n,2)
      else
          js = Interp%jlat(n,2)
          je = Interp%jlat(n,1)
         fjs = Interp%facj(n,2)
         fje = Interp%facj(n,1)
      endif

   do m = 1, nlon_out

    ! longitude window
       is = Interp%ilon(m,1)
       ie = Interp%ilon(m,2)
      fis = Interp%faci(m,1)
      fie = Interp%faci(m,2)
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
            call data_sum ( data_in(is:ie,js:je), area_in(is:ie,js:je), &
                            fis, fie, fjs,fje,                          &
                            dwtsum, wtsum, arsum, mask_in(is:ie,js:je)  )
         else
            call data_sum ( data_in(is:ie,js:je), area_in(is:ie,js:je), &
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
!
! compute statistics: minimum, maximum, and mean
!
! NOTE: will not correctly compute statistics on more than one processor
!
!-----------------------------------------------------------------------

 if (iverbose > 0) then

 ! compute statistics of input data

   call stats (data_in, area_in, dsum, wsum, min_in, max_in, miss_in, mask_in)
   ! diagnostic messages
   if (wsum > 0.0) then
      avg_in=dsum/wsum
   else
      print *, 'horiz_interp stats: input area equals zero on pe=', pe
      avg_in=0.0
   endif
   if (iverbose > 1) print '(a,i4,2f16.11)', 'pe, sum area_in  = ', pe, sum(area_in), wsum


 ! compute statistics of output data

   call stats (data_out, area_out, dsum, wsum, min_out, max_out, miss_out, mask_out)
   ! diagnostic messages
   if (wsum > 0.0) then
      avg_out=dsum/wsum
   else
      print *, 'horiz_interp stats: output area equals zero on pe=', pe
      avg_out=0.0
   endif
   if (iverbose > 1) print '(a,i4,2f16.11)', 'pe, sum area_out = ', pe, sum(area_out), wsum

    !---- output statistics ----

      write (*,900)
      write (*,901)  min_in ,max_in ,avg_in
      if (present(mask_in))  write (*,903)  miss_in
      write (*,902)  min_out,max_out,avg_out
      if (present(mask_out)) write (*,903)  miss_out

 900  format (/,1x,10('-'),' output from horiz_interp ',10('-'))
 901  format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
 902  format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
 903  format ('          number of missing points = ',i6)

 endif

!-----------------------------------------------------------------------

 end subroutine horiz_interp_base_2d


!#######################################################################

 subroutine new_horiz_interp_base_2d ( Interp, data_in, data_out, mask_in,mask_out)


!-----------------------------------------------------------------------
!
!   interpolates from a uniformly spaced grid to any output grid
!   using bilinear_interpolation.
!
!  input:
!  -----
!     Interp     Derived-type variable containing interpolation
!                indices and weights. Returned by a previous call
!                to horiz_interp_init.
!
!     data_in    Input data on input grid defined by grid box edges
!                initialized in variable Interp.
!                  [real, dimension(:,:)]
!
!  output:
!  ------
!     data_out   Output data on output grid defined by grid box edges
!                initialized in variable Interp.
!                  [real, dimension(:,:)]
!
!  optional
!  --------
!     verbose   flag for the amount of print output (integer scalar)
!                 =0 no output; =1 min,max,means; =2 still more
!
!     mask_in   input mask;  =0.0 for data points that should not
!               be used or have no data; has the same size as data_in
!     mask_out  output mask that specifies whether data was computed
!
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(out), dimension(:,:) :: data_out
      real, intent(in), dimension(:,:), optional :: mask_in
      real, intent(out), dimension(:,:), optional :: mask_out
!-----------------------------------------------------------------------
      integer :: nlon_in, nlat_in, nlon_out, nlat_out, n, m, &
                 is, ie, js, je
      real    :: wtw, wte, wts, wtn


      nlon_in = size(Interp%dlon_in)
      nlat_in = size(Interp%dsph_in) 

      nlon_out = size(Interp%wti,1)
      nlat_out = size(Interp%wti, 2)

   if (size(data_in,1) /= nlon_in .or. size(data_in,2) /= nlat_in) &
   call error_handler ('size of input array incorrect')

   if (size(data_out,1) /= nlon_out .or. size(data_out,2) /= nlat_out) &
   call error_handler ('size of output array incorrect')

   do n = 1, nlat_out
   do m = 1, nlon_out

   is = Interp % i_lon (m,n,1)
   ie = Interp % i_lon (m,n,2)
   js = Interp % j_lat (m,n,1)
   je = Interp % j_lat (m,n,2)

   wtw = Interp % wti   (m,n,1)
   wte = Interp % wti   (m,n,2)
   wts = Interp % wtj   (m,n,1)
   wtn = Interp % wtj   (m,n,2)

   if(present(mask_in)) then
      data_out(m,n) = data_in(is,js)*mask_in(is,js)*wtw*wts + &
                      data_in(ie,js)*mask_in(ie,js)*wte*wts + &
                      data_in(ie,je)*mask_in(ie,je)*wte*wtn + &
                      data_in(is,je)*mask_in(is,je)*wtw*wtn
   else
      data_out(m,n) = data_in(is,js)*wtw*wts + data_in(ie,js) * wte*wts +  &
                   data_in(ie,je)*wte*wtn + data_in(is,je) * wtw*wtn
   endif

   enddo
   enddo

   return

 end subroutine new_horiz_interp_base_2d


!#######################################################################

 subroutine stats ( dat, area, dsum, wsum, low, high, miss, mask )
 real,    intent(in)  :: dat(:,:), area(:,:)
 real,    intent(out) :: dsum, wsum, low, high
 integer, intent(out) :: miss
 real,    intent(in), optional :: mask(:,:)

 ! sum data, data*area; and find min,max

   if (present(mask)) then
      dsum = sum(area(:,:)*dat(:,:)*mask(:,:))
      wsum = sum(area(:,:)*mask(:,:))
      miss = count(mask(:,:) <= 0.5)
      low  = minval(dat(:,:),mask=mask(:,:) > 0.5)
      high = maxval(dat(:,:),mask=mask(:,:) > 0.5)
   else
      dsum = sum(area(:,:)*dat(:,:))
      wsum = sum(area(:,:))
      miss = 0
      low  = minval(dat(:,:))
      high = maxval(dat(:,:))
   endif

 end subroutine stats

!#######################################################################

 subroutine horiz_interp_init_2 ( Interp, blon_in, blat_in,   &
                                  lon_out, lat_out, verbose, bilinear_interp  )

!-----------------------------------------------------------------------
!
!   Allocates space and initializes a derived-type variable
!   that contains pre-computed interpolation indices and weights.
!
!  input:
!  -----
!
!   blon_in    The longitude edges (in radians) for input data grid boxes.
!              The values are for adjacent grid boxes and must increase in
!              value. If there are M longitude grid boxes there must be 
!              M+1 edge values.    [real, dimension(:)]
!
!   blat_in    The latitude edges (in radians) for input data grid boxes.
!              The values are for adjacent grid boxes and may increase or
!              decrease in value. If there are N latitude grid boxes there
!              must be N+1 edge values.    [real, dimension(:)]
!
!   lon_out  when using conservative scheme, 
!             The longitude edges (in radians) for output data grid boxes.
!              The edge values are stored as pairs for each grid box.  The
!              pairs do not have to be adjacent.  If there are M longitude
!              grid boxes in the output grid, then lon_out is dimensioned (M,2).
!                 [real, dimension(:,2)]
!            when using bilinear interpolation,
!              The geographical longitude (in radians) for output data grid boxes.
!   lat_out  when using conservative scheme, 
!              The latitude edges (in radians) for output data grid boxes.
!              The edge values are stored as pairs for each grid box.  The
!              pairs do not have to be adjacent.  If there are N latitude
!              grid boxes in the output grid, then lat_out is dimensioned (N,2).
!                 [real, dimension(:,2)]
!            when using bilinear interpolation,
!              The geographical latitude (in radians) for output data grid boxes.
!  input/output:
!  ------------
!     Interp     A derived-type variable containing interpolation
!                indices and weights.
!
!  optional
!  --------
!     verbose   flag for the amount of print output (integer scalar)
!                 =0 no output; =1 min,max,means; =2 still more
!
! bilinear_interp flag for the interplation method. Default value .false.
!                 indicates using conservative scheme, .true. indicates 
!                 using bilinear interpolation
!-----------------------------------------------------------------------
 type(horiz_interp_type), intent(inout) :: Interp
      real, intent(in),  dimension(:)   :: blon_in , blat_in
      real, intent(in),  dimension(:,:) :: lon_out, lat_out
   integer, intent(in),                   optional :: verbose
      logical, intent(in),                optional :: bilinear_interp
!-----------------------------------------------------------------------
      real, dimension(size(lat_out,1),size(lat_out,2)) :: sph
      real, dimension(size(blat_in)) :: slat_in

!-----------------------------------------------------------------------

   real    :: blon, fac, hpi, tpi, eps
   integer :: i, j, m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
              unit, np, npass, iverbose, m2, n2, iter
   logical :: s2n
   character(len=64) :: mesg
   logical           :: using_bilinear
!-----------------------------------------------------------------------
   iverbose = 0;  if (present(verbose)) iverbose = verbose
   using_bilinear = .false. 
   if (present(bilinear_interp)) using_bilinear = bilinear_interp

   Interp % interp_method = 1

   if(using_bilinear) then
      call horiz_interp_init_2_new(blon_in, blat_in,   &
                                          lon_out, lat_out, Interp, verbose=iverbose )
      Interp % interp_method = 2
      return
    endif

   allocate ( Interp % dlon_in  (size(blon_in)-1),  &
              Interp % dsph_in  (size(blat_in)-1),  &
              Interp % dlon_out (size(lon_out,1)), &
              Interp % dsph_out (size(lat_out,1)), &
              Interp % facj (size(lat_out,1),size(lat_out,2)), &
              Interp % jlat (size(lat_out,1),size(lat_out,2)), &
              Interp % faci (size(lon_out,1),size(lon_out,2)), &
              Interp % ilon (size(lon_out,1),size(lon_out,2))  )

!-----------------------------------------------------------------------


!  write version number and tag name
   if (do_vers) then
      call write_version_number (version, tagname)
      do_vers = .false.
   endif

   hpi = 0.5*pi
   tpi = 4.*hpi

   nlon_in = size(blon_in)-1;  nlat_in = size(blat_in)-1

! check size of input arguments

   if ( size(lon_out,2) /=2 .or. size(lat_out,2) /= 2 )  &
   call error_handler ('when using conservative scheme, dimension 2 of lon_out and/or lat_out must be 2')

!-----------------------------------------------------------------------
!  --- set-up for area of input grid boxes ---

   do j = 1, nlat_in+1
       slat_in(j) = sin(blat_in(j))
   enddo

   do j = 1, nlat_in
       Interp % dsph_in(j) = abs(slat_in(j+1)-slat_in(j))
   enddo

   do i = 1,nlon_in
       Interp % dlon_in(i) = abs(blon_in(i+1)-blon_in(i))
   enddo

!  set south to north flag
   s2n = .true.
   if (blat_in(1) > blat_in(nlat_in+1)) s2n = .false.

!-----------------------------------------------------------------------
!  --- set-up for area of output grid boxes ---

   nlon_out = size(lon_out,1);  nlat_out = size(lat_out,1)

   do n = 1, nlat_out
       Interp % dsph_out(n) = abs(sin(lat_out(n,2))-sin(lat_out(n,1)))
   enddo

   do m = 1,nlon_out
          Interp % dlon_out(m) = abs(lon_out(m,2)-lon_out(m,1))
   enddo

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
        if ( blon < blon_in(1)         ) blon = blon + tpi
        if ( blon > blon_in(nlon_in+1) ) blon = blon - tpi
        eps = 0.0
    do iter=1,num_iters
  ! find indices from input longitudes
    do i = 1, nlon_in
        if ( (blon_in(i)-blon) <= eps .and. &
             (blon-blon_in(i+1)) <= eps ) then
             Interp%ilon(m,m2) = i
             fac = (blon-blon_in(i))/(blon_in(i+1)-blon_in(i))
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
              lon_out(m,m2),blon,blon_in(1),blon_in(nlon_in+1),eps
           call error_handler ('no longitude index found')
       endif
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

 end subroutine horiz_interp_init_2

!#######################################################################

 subroutine horiz_interp_init_2_new ( blon_in, blat_in,   &
                                      lon_out, lat_out, Interp, verbose )

!-----------------------------------------------------------------------
!
!   Allocates space and initializes a derived-type variable
!   that contains pre-computed interpolation indices and weights.
!
!  input:
!  -----
!
!   blon_in    The longitude edges (in radians) for input data grid boxes.
!              The values are for adjacent grid boxes and must increase in
!              value. If there are M longitude grid boxes there must be 
!              M+1 edge values.    [real, dimension(:)]
!
!   blat_in    The latitude edges (in radians) for input data grid boxes.
!              The values are for adjacent grid boxes and may increase or
!              decrease in value. If there are N latitude grid boxes there
!              must be N+1 edge values.    [real, dimension(:)]
!
!   lon_out   The geographical longitude (in radians) for output data grid boxes.
!
!   lat_out   The geographical latitude  (in radians) for output data grid boxes.
!
!  input/output:
!  ------------
!     Interp     A derived-type variable containing interpolation
!                indices and weights.
!
!  optional
!  --------
!     verbose   flag for the amount of print output (integer scalar)
!                 =0 no output; =1 min,max,means; =2 still more
!
!-----------------------------------------------------------------------
 type(horiz_interp_type), intent(inout) :: Interp
      real, intent(in),  dimension(:)   :: blon_in , blat_in
      real, intent(in),  dimension(:,:) :: lon_out, lat_out
   integer, intent(in),                   optional :: verbose
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  integer :: nlon_in, nlat_in, nlon_out, nlat_out, n, m,     &
             i, j, ie, is, je, js, istart, jstart, iend, jend, &
             ln_err, lt_err, warns
  real    :: wtw, wte, wts, wtn, lon, lat, epsln, tpi, hpi,           &
             glt_min, glt_max, gln_min, gln_max
  real, dimension(size(blon_in)-1)  :: blon_c
  real, dimension(size(blat_in)-1)  :: blat_c

  hpi = 0.5*pi
  tpi = 4.0*hpi
  epsln = 1.0e-10
  glt_min = 90.
  glt_max = -90.
  gln_min = 360.
  gln_max = -360.
  ln_err = 0
  lt_err = 0
!-----------------------------------------------------------------------

   allocate ( Interp % dlon_in  (size(blon_in)-1),                &
              Interp % dsph_in  (size(blat_in)-1),                &
              Interp % wti (size(lon_out,1),size(lon_out,2),2),   &
              Interp % wtj (size(lon_out,1),size(lon_out,2),2),   &
              Interp % i_lon (size(lon_out,1),size(lon_out,2),2), &
              Interp % j_lat (size(lon_out,1),size(lon_out,2),2))
!-----------------------------------------------------------------------
   warns = 0
   if(present(verbose)) warns = verbose

!  write version number and tag name
   if (do_vers) then
      call write_version_number (version, tagname)
      do_vers = .false.
   endif

   if(size(lon_out,1) /= size(lat_out,1) .or. size(lon_out,2) /= size(lat_out,2) ) &
   call error_handler ('when using bilinear interplation, the output grids should be geographical grids')    
      nlon_in = size(blon_in)-1
      nlat_in = size(blat_in)-1

      do i =1 , nlon_in
         blon_c(i) = (blon_in(i) + blon_in(i+1))*0.5
      enddo
      do j =1 , nlat_in
         blat_c(j) = (blat_in(j) + blat_in(j+1))*0.5
      enddo 

!
!     find longitude points of data within interval [0., 360.]
  istart = 1
  do i=2,nlon_in
    if (blon_in(i-1) .lt. 0. .and. blon_in(i) .ge. 0.) istart = i
  enddo 
  iend = nlon_in
  do i=2,nlon_in
    if (blon_in(i-1) .lt.tpi  .and. blon_in(i) .ge. tpi) iend = i
  enddo
!
!     find latitude points of data within interval [-90., 90.]
  jstart = 1
  do j=2,nlat_in
    if (blat_in(j-1) .lt. -1.0 * hpi .and. blat_in(j) .ge. -1.0*hpi) jstart = j
  enddo 
  jend = nlat_in
  do j=2,nlat_in
    if (blat_in(j-1) .lt. hpi .and. blat_in(j) .ge.hpi ) jend = j
  enddo


   nlon_out = size(lon_out, 1)
   nlat_out = size(lon_out, 2)

   do n = 1, nlat_out
   do m = 1, nlon_out
      lon = lon_out(m,n)
      lat = lat_out(m,n)
      if(lon .lt. 0.) lon = lon + tpi
      if(lon .ge. tpi) lon = lon - tpi
      glt_min = min(lat,glt_min)
      glt_max = max(lat,glt_max)
      gln_min = min(lon,gln_min)
      gln_max = max(lon,gln_max)
      is = indp(lon, blon_c(istart),iend - istart + 1) + istart - 1

      if( blon_c(is) .gt. lon ) is = is - 1
      ie = is + 1
      if(is .ge. 1 .and. ie .le. nlon_in) then
         wtw = ( blon_c(ie) - lon) / (blon_c(ie) - blon_c(is) )
      else
!     east or west of the last data value. this could be because a
!     cyclic condition is needed or the dataset is too small. in either 
!     case apply a cyclic condition
         ln_err = 1
         ie = istart
         is = iend
         if (blon_c(ie) .ge. lon ) then
            wtw = (blon_c(ie) -lon)/(blon_c(ie)-blon_c(is)+tpi+epsln)
         else
            wtw = (blon_c(ie) -lon+tpi+epsln)/(blon_c(ie)-blon_c(is)+tpi+epsln)
         endif
       endif
    
       wte = 1. - wtw

      js = indp(lat, blat_c(jstart), jend - jstart + 1) + jstart - 1

      if( blat_c(js) .gt. lat ) js = max(js - 1, jstart)
      je = min(js + 1, jend)

      if ( blat_c(js) .ne. blat_c(je) .and. blat_c(js) .le. lat) then
         wts = ( blat_c(je) - lat )/(blat_c(je)-blat_c(js))
      else
!     north or south of the last data value. this could be because a
!     pole is not included in the data set or the dataset is too small.
!     in either case extrapolate north or south
         lt_err = 1
         wts = 1.
      endif

      wtn = 1. - wts
   
      Interp % i_lon (m,n,1) = is
      Interp % i_lon (m,n,2) = ie
      Interp % j_lat (m,n,1) = js
      Interp % j_lat (m,n,2) = je

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
    write (stdout,'(/,(1x,a,2f8.2))')                                &
      '    data required between longitudes:', gln_min, gln_max,     &
      '      data set is between longitudes:', blon_c(istart), blon_c(iend)
    warns = warns - 1
  endif
!
  if (lt_err .eq. 1 .and. warns > 0) then
    write (stdout,'(/,(1x,a))')                                     &
      '==> Warning: the geographic data set does not extend far   ',&
      '             enough north or south - extrapolation from    ',&
      '             the nearest data was applied. this may create ',&
      '             artificial gradients near a geographic pole   ' 
    write (stdout,'(/,(1x,a,2f8.2))')                             &
      '    data required between latitudes:', glt_min, glt_max,   &
      '      data set is between latitudes:', blat_c(jstart), blat_c(jend)
    warns = warns - 1
  endif

  return

 end subroutine horiz_interp_init_2_new

!#######################################################################

 subroutine data_sum ( data, area, facis, facie, facjs, facje,  &
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

 end subroutine data_sum

!#######################################################################

 subroutine horiz_interp_solo_1 ( data_in, blon_in, blat_in,    &
                                  blon_out, blat_out, data_out, &
                                  verbose, mask_in, mask_out    )

!-----------------------------------------------------------------------
!       Overloaded version of interface horiz_interp_solo_2
!
!   uses 1d arrays of adjacent values for blon_out and blat_out
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:)   :: blon_in , blat_in
      real, intent(in),  dimension(:)   :: blon_out, blat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
!-----------------------------------------------------------------------
     real, dimension(size(blon_out)-1,2) :: lonb
     real, dimension(size(blat_out)-1,2) :: latb
     integer :: i, j
!-----------------------------------------------------------------------
!-- create new grid box boundary format ---

   do i=1,size(blon_out)-1
     lonb(i,1) = blon_out(i)
     lonb(i,2) = blon_out(i+1)
   enddo
   do j=1,size(blat_out)-1
     latb(j,1) = blat_out(j)
     latb(j,2) = blat_out(j+1)
   enddo

   call horiz_interp_solo_2 ( data_in, blon_in, blat_in, &
                              lonb, latb, data_out,      &
                              verbose, mask_in, mask_out )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1

!#######################################################################

 subroutine horiz_interp_solo_2 ( data_in, blon_in, blat_in,    &
                                  blon_out, blat_out, data_out, &
                                  verbose, mask_in, mask_out    )

!-----------------------------------------------------------------------
!
!   interpolates from a uniformly spaced grid to any output grid
!   using an area weighing scheme to conserve the global integral
!   of the field being interpolated.
!
!  input:
!  -----
!     data_in   input data on long/lat grid
!                 [real, dimension(:,:)]
!
!     blon_in   contiguous longitude boundaries for input data
!               grid boxes, dimensioned by size(data_in,1)+1
!     blat_in   contiguous latitudes boundaries for input data
!               grid boxes, dimensioned by size(data_in,2)+1
!
!     blon_out  longitude boundaries for output data grid boxes,
!                 [real, dimension(size(data_out,1),2)]
!     blat_out  latitude boundaries for output data at grid boxes,
!                 [real, dimension(size(data_out,2),2)]
!
!  output:
!  ------
!     data_out   output data on long/lat grid
!                  [real, dimension(:,:)]
!
!  optional
!  --------
!     verbose   flag for the amount of print output (integer scalar)
!                 =0 no output; =1 min,max,means; =2 still more
!
!     mask_in   input mask;  =0.0 for data points that should not
!               be used or have no data; has the same size as data_in
!     mask_out  output mask that specifies whether data was computed
!
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:)   :: blon_in , blat_in
      real, intent(in),  dimension(:,:) :: blon_out, blat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out

    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------

    call horiz_interp_init_2 ( Interp, blon_in, blat_in,   &
                               blon_out, blat_out, verbose )

    call horiz_interp_base_2d ( Interp, data_in, data_out, &
                                verbose, mask_in, mask_out )

    call horiz_interp_end ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_2

!#######################################################################

 subroutine horiz_interp_init_1 ( Interp, blon_in, blat_in,   &
                                  blon_out, blat_out, verbose )

!-----------------------------------------------------------------------
!       Overloaded version of interface horiz_interp_init_2
!
!   uses 1d arrays of adjacent values for blon_out and blat_out
!-----------------------------------------------------------------------
 type(horiz_interp_type), intent(inout) :: Interp
      real, intent(in),  dimension(:)   :: blon_in , blat_in
      real, intent(in),  dimension(:)   :: blon_out, blat_out
   integer, intent(in),                   optional :: verbose
!-----------------------------------------------------------------------
     real, dimension(size(blon_out)-1,2) :: lonb
     real, dimension(size(blat_out)-1,2) :: latb
     integer :: i, j
!-----------------------------------------------------------------------
!-- create new grid box boundary format ---

   do i=1,size(blon_out)-1
     lonb(i,1) = blon_out(i)
     lonb(i,2) = blon_out(i+1)
   enddo
   do j=1,size(blat_out)-1
     latb(j,1) = blat_out(j)
     latb(j,2) = blat_out(j+1)
   enddo

   call horiz_interp_init_2 ( Interp, blon_in, blat_in, &
                              lonb, latb, verbose       )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_init_1

!#######################################################################

 subroutine horiz_interp_solo_old (data_in, wb, sb, dx, dy,  &
                                   blon_out, blat_out, data_out,  &
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
!   blon_out    The longitude edges (in radians) for output data grid boxes.
!               The values are for adjacent grid boxes and must increase in
!               value. If there are MLON grid boxes there must be MLON+1
!               edge values.  [real, dimension(:)]
!
!   blat_out    The latitude edges (in radians) for output data grid boxes.
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
      real, intent(in),  dimension(:)   :: blon_out, blat_out
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


   call horiz_interp_solo_1 (data_in, blon_in, blat_in,    &
                             blon_out, blat_out, data_out, &
                             verbose, mask_in, mask_out    )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_old

!#######################################################################

 subroutine horiz_interp_end ( Interp )

   type (horiz_interp_type), intent(inout) :: Interp

!-----------------------------------------------------------------------
!  releases space used by horiz_interp_type variables
!  must be called before re-initializing the same variable
!-----------------------------------------------------------------------
   if(Interp % interp_method ==1 ) then
      deallocate ( Interp % dlon_in , Interp % dsph_in , &
                Interp % dlon_out, Interp % dsph_out, &
                Interp % facj, Interp % jlat, &
                Interp % faci, Interp % ilon  )
   else if(Interp % interp_method ==2 ) then
      deallocate (Interp % dlon_in, Interp % dsph_in  ,  &
                  Interp % wti,         Interp % wtj ,   &
                  Interp % i_lon,    Interp % j_lat )
   endif

!-----------------------------------------------------------------------

 end subroutine horiz_interp_end

!#######################################################################

 subroutine horiz_interp_base_3d ( Interp, data_in, data_out, &
                                   verbose, mask_in, mask_out )

!-----------------------------------------------------------------------
!   overload of interface horiz_interp_base_2d
!
!   uses 3d arrays for data and mask
!   this allows for multiple interpolations with one call
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
      real, intent(in),  dimension(:,:,:) :: data_in
      real, intent(out), dimension(:,:,:) :: data_out
   integer, intent(in),                     optional :: verbose
      real, intent(in),   dimension(:,:,:), optional :: mask_in
      real, intent(out),  dimension(:,:,:), optional :: mask_out
!-----------------------------------------------------------------------
   integer :: n

   do n = 1, size(data_in,3)
     if (present(mask_in))then
       call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                                   verbose, mask_in(:,:,n), mask_out(:,:,n) )
     else
       call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                                   verbose )
     endif
   enddo
  

!-----------------------------------------------------------------------

 end subroutine horiz_interp_base_3d

!#######################################################################

 subroutine error_handler ( message )
 character(len=*), intent(in) :: message

   call error_mesg ('horiz_interp_mod', message, FATAL)

 end subroutine error_handler

!#######################################################################


function indp (value, array, ia)
integer             :: ia, indp
real, dimension(ia) :: array
real                :: value
!
!=======================================================================
!
!     indp = index of nearest data point within "array" corresponding to
!            "value".
!
!     inputs:
!
!     value  = arbitrary data...same units as elements in "array"
!     array  = array of data points  (must be monotonically increasing)
!     ia     = dimension of "array"
!
!     output:
!
!     indp =  index of nearest data point to "value"
!             if "value" is outside the domain of "array" then indp = 1
!             or "ia" depending on whether array(1) or array(ia) is
!             closest to "value"
!
!             note: if "array" is dimensioned array(0:ia) in the calling
!                   program, then the returned index should be reduced
!                   by one to account for the zero base.
!
!     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
!
!     example:
!
!     let model depths be defined by the following:
!     parameter (km=5)
!     dimension z(km)
!     data z /5.0, 10.0, 50.0, 100.0, 250.0/
!
!     k1 = indp (12.5, z, km)
!     k2 = indp (0.0, z, km)
!
!     k1 would be set to 2, & k2 would be set to 1 so that
!     z(k1) would be the nearest data point to 12.5 and z(k2) would
!     be the nearest data point to 0.0
!
!=======================================================================
!
!
integer i, ii
logical keep_going
!
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

end module horiz_interp_mod

