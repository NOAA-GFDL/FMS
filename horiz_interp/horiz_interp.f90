
module horiz_interp_mod

!-----------------------------------------------------------------------

 use utilities_mod, only:  open_file, print_version_number,  &
                           error_mesg, FATAL, get_my_pe

 implicit none
 private

!---- interfaces ----

 public   horiz_interp_type, horiz_interp, horiz_interp_init, &
          horiz_interp_end

 interface horiz_interp
    module procedure horiz_interp_0
    module procedure horiz_interp_1
    module procedure horiz_interp_2
    module procedure horiz_interp_old
 end interface

 interface horiz_interp_init
    module procedure horiz_interp_init_1
    module procedure horiz_interp_init_2
end interface

 type horiz_interp_type
   real,    dimension(:), pointer :: dlon_in, dlon_out, &
                                     dsph_in, dsph_out
   real,    dimension(:,:), pointer :: faci, facj
   integer, dimension(:,:), pointer :: ilon, jlat
 end type

!-----------------------------------------------------------------------
 character(len=4) :: vers_num = 'v2.0'
 logical :: do_vers = .true.
 integer :: num_iters = 4
!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine horiz_interp_2 ( data_in, blon_in, blat_in,    &
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
!     data_in     input data; dimensioned by nlon_in x nlat_in
!                      stored from south to north
!
!     blon_in   contiguous longitude boundaries for input data
!               grid boxes, dimensioned by size(data_in,1)+1
!     blat_in   contiguous latitudes boundaries for input data
!               grid boxes, dimensioned by size(data_in,2)+1
!
!     blon_out  longitude boundaries for output data grid boxes,
!                  dimensioned by size(data_out,1) x 2
!     blat_out  latitude boundaries for output data at grid boxes,
!                  dimensioned by size(data_out,2) x 2
!
!  output:
!  ------
!     data_out     output data
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

    call horiz_interp_0 ( Interp, data_in, data_out, &
                          verbose, mask_in, mask_out )

    call horiz_interp_end ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_2

!#######################################################################

 subroutine horiz_interp_0 ( Interp, data_in, data_out, &
                             verbose, mask_in, mask_out )

!-----------------------------------------------------------------------
!
!   interpolates from a uniformly spaced grid to any output grid
!   using an area weighing scheme to conserve the global integral
!   of the field being interpolated.
!
!  input:
!  -----
!     data_in     input data; dimensioned by nlon_in x nlat_in
!                      stored from south to north
!
!  output:
!  ------
!     data_out     output data
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
      real, dimension(size(Interp%dlon_out), &
                      size(Interp%dsph_out)) :: area_out

      integer :: i, j, m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
                 miss_in, miss_out, unit, is, ie, js, je,   &
                 np, npass, iverbose, m2, n2

      real :: cph, dsum, wsum, avg_in, min_in, max_in,   &
              avg_out, min_out, max_out, blon, eps,    &
              dwtsum, wtsum, hpie, tpie, dtr, dsph, fs, fe

      character(len=64) :: mesg

!-----------------------------------------------------------------------

   iverbose = 0;  if (present(verbose)) iverbose = verbose

   hpie = acos(0.0)
   tpie = 4.*hpie
   dtr  = hpie/90.
   eps  = epsilon(wtsum)


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

!-----------------------------------------------------------------------

   do n = 1, nlat_out
      js = Interp%jlat(n,1)
      je = Interp%jlat(n,2)

   do m = 1, nlon_out

      is = Interp%ilon(m,1)
      ie = Interp%ilon(m,2)
      fs = Interp%faci(m,1)
      fe = Interp%faci(m,2)
      npass = 1
      dwtsum = 0.
       wtsum = 0.

      if ( ie < is ) then
          ie = nlon_in
          fe = 1.0
          npass = 2
      endif

      do np = 1, npass

         if ( np == 2 ) then
             is = 1
             fs = 1.0
             ie = Interp%ilon(m,2)
             fe = Interp%faci(m,2)
         endif

         call data_sum ( data_in(is:ie,js:je), area_in(is:ie,js:je), &
                         fs, fe, Interp%facj(n,1),Interp%facj(n,2),  &
                         dwtsum, wtsum                               )

      enddo

      if (wtsum > eps) then
         data_out(m,n) = dwtsum/wtsum
         if (present(mask_out)) mask_out(m,n) = 1.0
      else
         data_out(m,n) = 0.
         if (present(mask_out)) mask_out(m,n) = 0.0
      endif

   enddo
   enddo

!***********************************************************************

!  compute statistics: minimum, maximum, and mean

!  WILL NOT WORK ON MORE THAN ONE PE
!-----------------------------------------------------------------------

 if (iverbose > 0) then

!--------- statistics of input data ----------

      dsum = sum(area_in(:,:)*data_in(:,:))
      wsum = sum(area_in(:,:))
      miss_in = 0
      min_in = minval(data_in(:,:))
      max_in = maxval(data_in(:,:))
      avg_in = dsum/wsum

   if (present(mask_in)) then
      miss_in = count(mask_in(:,:) <= 0.5)
      min_in=minval(data_in(:,:),mask=mask_in(:,:) > 0.5)
      max_in=maxval(data_in(:,:),mask=mask_in(:,:) > 0.5)
      avg_in=dsum/wsum
   else
      miss_in = 0
      min_in=minval(data_in(:,:))
      max_in=maxval(data_in(:,:))
      avg_in=dsum/wsum
   endif
      if (miss_in >= nlon_in*nlat_in) then
          call error_mesg ('horiz_interp', 'no input data', FATAL)
      endif
      if (iverbose > 1) print *, 'pe, sum area_in = ', get_my_pe(), wsum


!--------- statistics of output data ----------

      if (present(mask_out)) then
          avg_out=sum(area_out*mask_out*data_out)/sum(area_out*mask_out)
          min_out=minval(data_out,mask=mask_out > 0.5)
          max_out=maxval(data_out,mask=mask_out > 0.5)
          miss_out=count(mask_out <= 0.5)
      else
          avg_out=sum(area_out*data_out)/sum(area_out)
          min_out=minval(data_out)
          max_out=maxval(data_out)
          miss_out=0
      endif
      if (iverbose > 1) print *, 'pe, sum area_out = ',  &
                                         get_my_pe(), sum(area_out)

      write (*,900)
      write (*,901)  min_in ,max_in ,avg_in
      if (present(mask_in))  write (*,903)  miss_in
      write (*,902)  min_out,max_out,avg_out
      if (present(mask_out)) write (*,903)  miss_out

 900  format (/,1x,10('-'),' output from horiz_interp ',10('-'))
 901  format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f16.9)
 902  format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f16.9)
 903  format ('          number of missing points = ',i6)

 endif

!-----------------------------------------------------------------------

 end subroutine horiz_interp_0

!#######################################################################

 subroutine horiz_interp_init_2 ( Interp, blon_in, blat_in,   &
                                  blon_out, blat_out, verbose )

!-----------------------------------------------------------------------
!
!   interpolates from a uniformly spaced grid to any output grid
!   using an area weighing scheme to conserve the global integral
!   of the field being interpolated.
!
!  input:
!  -----
!
!     blon_in   contiguous longitude boundaries for input data
!               grid boxes, dimensioned by size(data_in,1)+1
!     blat_in   contiguous latitudes boundaries for input data 
!               grid boxes, dimensioned by size(data_in,2)+1
!
!     blon_out  longitude boundaries for output data grid boxes,
!                  dimensioned by size(data_out,1) x 2
!     blat_out  latitude boundaries for output data at grid boxes,
!                  dimensioned by size(data_out,2) x 2
!
!  input/output:
!  ------------
!     Interp     derived-type variable 
!
!  optional
!  --------
!     verbose   flag for the amount of print output (integer scalar)
!                 =0 no output; =1 min,max,means; =2 still more
!
!-----------------------------------------------------------------------
 type(horiz_interp_type), intent(inout) :: Interp
      real, intent(in),  dimension(:)   :: blon_in , blat_in
      real, intent(in),  dimension(:,:) :: blon_out, blat_out
   integer, intent(in),                   optional :: verbose
!-----------------------------------------------------------------------
      real, dimension(size(blat_out,1),size(blat_out,2)) :: ph

!-----------------------------------------------------------------------

   integer :: i, j, m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
              unit, np, npass, iverbose, m2, n2, iter

   real ::  blon, fac, hpie, tpie, eps

   character(len=64) :: mesg

!-----------------------------------------------------------------------

   allocate ( Interp % dlon_in  (size(blon_in)-1),  &
              Interp % dsph_in  (size(blat_in)-1),  &
              Interp % dlon_out (size(blon_out,1)), &
              Interp % dsph_out (size(blat_out,1)), &
              Interp % facj (size(blat_out,1),size(blat_out,2)), &
              Interp % jlat (size(blat_out,1),size(blat_out,2)), &
              Interp % faci (size(blon_out,1),size(blon_out,2)), &
              Interp % ilon (size(blon_out,1),size(blon_out,2))  )

!-----------------------------------------------------------------------

   if (do_vers) then
       unit = open_file ('logfile.out', action='append')
       call print_version_number (unit, 'horiz_interp', vers_num)
       close (unit)
       do_vers = .false.
   endif

   iverbose = 0;  if (present(verbose)) iverbose = verbose

   hpie = acos(0.0)
   tpie = 4.*hpie

   nlon_in = size(blon_in)-1;  nlat_in = size(blat_in)-1

!-----------------------------------------------------------------------
!  --- set-up for area of input grid boxes ---

   do j = 1, nlat_in
       Interp % dsph_in(j) = abs(sin(blat_in(j+1))-sin(blat_in(j)))
   enddo

   do i = 1,nlon_in
       Interp % dlon_in(i) = abs(blon_in(i+1)-blon_in(i))
   enddo

!-----------------------------------------------------------------------
!  --- set-up for area of output grid boxes ---

   nlon_out = size(blon_out,1);  nlat_out = size(blat_out,1)

   do n = 1, nlat_out
       Interp % dsph_out(n) = abs(sin(blat_out(n,2))-sin(blat_out(n,1)))
   enddo

   do m = 1,nlon_out
          Interp % dlon_out(m) = abs(blon_out(m,2)-blon_out(m,1))
   enddo

!***********************************************************************

!------ set up latitudinal indexing ------
!------ make sure output grid goes south to north ------

 do n = 1, nlat_out
    if (blat_out(n,1) < blat_out(n,2)) then
       ph(n,1) = blat_out(n,1)
       ph(n,2) = blat_out(n,2)
    else
       ph(n,1) = blat_out(n,2)
       ph(n,2) = blat_out(n,1)
    endif
 enddo

 Interp%jlat = 0
 do n2 = 1, 2
 do n = 1, nlat_out
     eps = 0.0
 do iter=1,num_iters
 do j = 1, nlat_in
    if ( (blat_in(j)-ph(n,n2)) <= eps .and.  &
         (ph(n,n2)-blat_in(j+1)) <= eps ) then
         Interp%jlat(n,n2) = j
         fac = (ph(n,n2)-blat_in(j))/(blat_in(j+1)-blat_in(j))
         if (n2 == 1) Interp%facj(n,n2) = 1.0 - fac
         if (n2 == 2) Interp%facj(n,n2) = fac
         exit
     endif
 enddo
     if ( Interp%jlat(n,n2) /= 0 ) exit
!    --- set tolerance for multiple passes ---
     eps  = epsilon(blon)*real(10**iter)
 enddo
     if ( Interp%jlat(n,n2) == 0 ) then
          write (mesg,888) n,ph(n,n2)
      888 format (': n,ph=',i3,f14.7,40x)
          call error_mesg ('horiz_interp_mod',  &
                   'no latitude index found'//trim(mesg), FATAL)
     endif
 enddo
 enddo

!------ set up longitudinal indexing ------

    Interp%ilon = 0
    do m2 = 1, 2
    do m = 1, nlon_out
        blon = blon_out(m,m2)
        if ( blon < blon_in(1)         ) blon = blon + tpie
        if ( blon > blon_in(nlon_in+1) ) blon = blon - tpie
        eps = 0.0
    do iter=1,num_iters
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
!      --- set tolerance for multiple passes ---
       eps  = epsilon(blon)*real(10**iter)
    enddo
       if ( Interp%ilon(m,m2) == 0 ) then
           print *, 'PE,blon_out,blon,blon_in,eps=', get_my_pe(), &
              blon_out(m,m2),blon,blon_in(1),blon_in(nlon_in+1),eps
           call error_mesg ('horiz_interp_mod', &
                            'no longitude index found', FATAL)
       endif
    enddo
    enddo

!-----------------------------------------------------------------------
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

 subroutine data_sum ( data, area, facis, facie, facjs, facje,  &
                       dwtsum, wtsum )

!-----------------------------------------------------------------------
   real, intent(in), dimension(:,:) :: data, area
   real, intent(in)                 :: facis, facie, facjs, facje
   real, intent(inout)              :: dwtsum, wtsum
!-----------------------------------------------------------------------
   real, dimension(size(area,1),size(area,2)) :: wt
   integer :: id, jd
!-----------------------------------------------------------------------

   id=size(area,1); jd=size(area,2)

   wt=area
   wt( 1,:)=wt( 1,:)*facis
   wt(id,:)=wt(id,:)*facie
   wt(:, 1)=wt(:, 1)*facjs
   wt(:,jd)=wt(:,jd)*facje

   dwtsum = dwtsum + sum(wt*data)
    wtsum =  wtsum + sum(wt)

!-----------------------------------------------------------------------

 end subroutine data_sum

!#######################################################################

 subroutine horiz_interp_1 (data_in, blon_in, blat_in, &
                            blon_out, blat_out, data_out,  &
                            verbose, mask_in, mask_out)

!-----------------------------------------------------------------------
!       Overloaded version of interface horiz_interp_2
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

   call horiz_interp_2 ( data_in, blon_in, blat_in, &
                         lonb, latb, data_out,      &
                         verbose, mask_in, mask_out )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_1

!#######################################################################

 subroutine horiz_interp_init_1 ( Interp, blon_in, blat_in,   &
                                  blon_out, blat_out, verbose )

!-----------------------------------------------------------------------
!       Overloaded version of interface horiz_interp_init_2
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

 subroutine horiz_interp_old (data_in, wb, sb, dx, dy,  &
                              blon_out, blat_out, data_out,  &
                              verbose, mask_in, mask_out)

!-----------------------------------------------------------------------
!       Overloaded version of interface horiz_interp_2
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
     real    :: pi, tpi
!-----------------------------------------------------------------------

   pi = 4.0*atan(1.0);  tpi = 2.*pi
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


   call horiz_interp_1 (data_in, blon_in, blat_in,    &
                        blon_out, blat_out, data_out, &
                        verbose, mask_in, mask_out    )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_old

!#######################################################################

 subroutine horiz_interp_end ( Interp )

   type (horiz_interp_type), intent(inout) :: Interp

!-----------------------------------------------------------------------

   deallocate ( Interp % dlon_in , Interp % dsph_in , &
                Interp % dlon_out, Interp % dsph_out, &
                Interp % facj, Interp % jlat, &
                Interp % faci, Interp % ilon  )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_end

!#######################################################################

end module horiz_interp_mod

