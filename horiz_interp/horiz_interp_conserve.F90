module horiz_interp_conserve_mod

  ! <CONTACT EMAIL="Bruce.Wyman@noaa.gov"> Bruce Wyman </CONTACT>
  ! <CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang </CONTACT>

  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

  ! <OVERVIEW>
  !   Performs spatial interpolation between grids using conservative interpolation
  ! </OVERVIEW>

  ! <DESCRIPTION>
  !     This module can interpolate data from regular rectangular grid
  !     to rectangular/tripolar grid. The interpolation scheme is area-averaging 
  !     conservative scheme. There is an optional mask field for missing input data.
  !     An optional output mask field may be used in conjunction with
  !     the input mask to show where output data exists.
  ! </DESCRIPTION>

  use mpp_mod,               only: mpp_send, mpp_recv, mpp_pe, mpp_root_pe
  use mpp_mod,               only: mpp_error, FATAL,  mpp_sync_self 
  use fms_mod,               only: write_version_number
  use constants_mod,         only: PI
  use horiz_interp_type_mod, only: horiz_interp_type


  implicit none
  private

  public :: horiz_interp_conserve_init, horiz_interp_conserve, horiz_interp_conserve_end

  integer :: pe, root_pe
  !-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: horiz_interp_conserve.F90,v 13.0 2006/03/28 21:39:36 fms Exp $'
  character(len=128) :: tagname = '$Name: memphis $'
  logical            :: do_vers = .true.

contains

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_conserve_init">

  !   <OVERVIEW>
  !      Initialization routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !      Allocates space and initializes a derived-type variable
  !      that contains pre-computed interpolation indices and weights.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_conserve_init ( Interp, lon_in, lat_in, lon_out, lat_out, verbose)

  !   </TEMPLATE>
  !   
  !   <IN NAME="lon_in" TYPE="real, dimension(:,:)" UNITS="radians">
  !      Longitude (in radians) for source data grid. 
  !   </IN>

  !   <IN NAME="lat_in" TYPE="real, dimension(:,:)" UNITS="radians">
  !      Latitude (in radians) for source data grid.
  !   </IN>

  !   <IN NAME="lon_out" TYPE="real, dimension(:,:)" UNITS="radians" >
  !      Longitude (in radians) for source data grid. 
  !   </IN>

  !   <IN NAME="lat_out" TYPE="real, dimension(:,:)" UNITS="radians" >
  !      Latitude (in radians) for source data grid. 
  !   </IN>

  !   <IN NAME="verbose" TYPE="integer, optional" >
  !      flag for the amount of print output.
  !   </IN>

  !   <INOUT NAME="Interp" TYPE="type(horiz_interp_type)" >
  !      A derived-type variable containing indices and weights used for subsequent 
  !      interpolations. To reinitialize this variable for a different grid-to-grid 
  !      interpolation you must first use the "horiz_interp_end" interface.
  !   </INOUT>
  subroutine horiz_interp_conserve_init ( Interp, lon_in, lat_in,   &
       lon_out, lat_out, verbose)

    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in),    dimension(:)      :: lon_in , lat_in
    real, intent(in),    dimension(:,:)    :: lon_out, lat_out
    integer, intent(in),       optional    :: verbose

    !-----------------------------------------------------------------------
    real, dimension(size(lat_out,1),size(lat_out,2)) :: sph
    real, dimension(size(lat_in(:)))    :: slat_in
    real, dimension(size(lon_in,1)-1) :: dlon_in
    real, dimension(size(lat_in,1)-1) :: dsph_in
    real, dimension(size(lon_out,1))  :: dlon_out
    real, dimension(size(lat_out,1))  :: dsph_out
    real    :: blon, fac, hpi, tpi, eps
    integer :: num_iters = 4
    integer :: i, j, m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
               iverbose, m2, n2, iter
    logical :: s2n
    character(len=64) :: mesg
    !-----------------------------------------------------------------------
    iverbose = 0;  if (present(verbose)) iverbose = verbose

    pe      = mpp_pe()
    root_pe = mpp_root_pe()

    if (do_vers) then
       call write_version_number (version, tagname)
       do_vers = .false.
    endif


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
         call mpp_error(FATAL, 'horiz_interp_conserve_mod: '// &
         'when using conservative scheme, dimension 2 of lon_out and/or lat_out must be 2')

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
710          format (': n,sph=',i3,f14.7,40x)
             call mpp_error(FATAL, 'horiz_interp_conserve_mod:no latitude index found'//trim(mesg))
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
             call mpp_error(FATAL, 'horiz_interp_conserve_mod: no longitude index found')
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
801    format (/,2x,'i',4x,'is',5x,'ie',4x,'facis',4x,'facie',  &
            /,(i4,2i7,2f10.5))
802    format (/,2x,'j',4x,'js',5x,'je',4x,'facjs',4x,'facje',  &
            /,(i4,2i7,2f10.5))
    endif
    !-----------------------------------------------------------------------

  end subroutine horiz_interp_conserve_init
  ! </SUBROUTINE>

  !########################################################################
  ! <SUBROUTINE NAME="horiz_interp_conserve">

  !   <OVERVIEW>
  !      Subroutine for performing the horizontal interpolation between two grids.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Subroutine for performing the horizontal interpolation between two grids. 
  !     horiz_interp_conserve_init must be called before calling this routine.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_conserve ( Interp, data_in, data_out, verbose, mask_in, mask_out)
  !   </TEMPLATE>
  !   
  !   <IN NAME="Interp" TYPE="type(horiz_interp_type)">
  !     Derived-type variable containing interpolation indices and weights.
  !     Returned by a previous call to horiz_interp_init.
  !   </IN>
  !   <IN NAME="data_in" TYPE="real, dimension(:,:)">
  !      Input data on source grid.
  !   </IN>

  !   <IN NAME="verbose" TYPE="integer, optional">
  !      flag for the amount of print output.
  !               verbose = 0, no output; = 1, min,max,means; = 2, still more
  !   </IN>
  !   <IN NAME="mask_in" TYPE="real, dimension(:,:),optional">
  !      Input mask, must be the same size as the input data. The real value of
  !      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points 
  !      that should not be used or have missing data.
  !   </IN>

  !   <OUT NAME="data_out" TYPE="real, dimension(:,:)">
  !      Output data on destination grid.
  !   </OUT>
  !   <OUT NAME="mask_out" TYPE="real, dimension(:,:),optional">
  !      Output mask that specifies whether data was computed.
  !   </OUT>

  subroutine horiz_interp_conserve ( Interp, data_in, data_out, verbose, &
       mask_in, mask_out)
    !-----------------------------------------------------------------------
    type (horiz_interp_type), intent(in) :: Interp
    real, intent(in),  dimension(:,:) :: data_in
    real, intent(out), dimension(:,:) :: data_out
    integer, intent(in),                   optional :: verbose
    real, intent(in),   dimension(:,:), optional :: mask_in
    real, intent(out),  dimension(:,:), optional :: mask_out
    !----------local variables----------------------------------------------------
    integer :: m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
         miss_in, miss_out, is, ie, js, je,   &
         np, npass, iverbose
    real    :: dsum, wsum, avg_in, min_in, max_in,   &
         avg_out, min_out, max_out, eps, asum,   &
         dwtsum, wtsum, arsum, fis, fie, fjs, fje
    !-----------------------------------------------------------------------
    iverbose = 0;  if (present(verbose)) iverbose = verbose

    eps = epsilon(wtsum)

    nlon_in  = Interp%nlon_src;  nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

    !  --- error checking ---
    if (size(data_in,1) /= nlon_in .or. size(data_in,2) /= nlat_in) &
         call mpp_error(FATAL, 'horiz_interp_conserve_mod: size of input array incorrect')

    if (size(data_out,1) /= nlon_out .or. size(data_out,2) /= nlat_out) &
         call mpp_error(FATAL, 'horiz_interp_conserve_mod: size of output array incorrect')

    if (present(mask_in)) then
       if ( count(mask_in < -.0001 .or. mask_in > 1.0001) > 0 ) &
            call mpp_error(FATAL, 'horiz_interp_conserve_mod: input mask not between 0,1')
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
                call data_sum ( data_in(is:ie,js:je), Interp%area_src(is:ie,js:je), &
                     fis, fie, fjs,fje, dwtsum, wtsum, arsum, mask_in(is:ie,js:je)  )
             else
                call data_sum ( data_in(is:ie,js:je), Interp%area_src(is:ie,js:je), &
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

       call stats(data_in, Interp%area_src, asum, dsum, wsum, min_in, max_in, miss_in, mask_in)
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
       call stats(data_out, Interp%area_dst, asum, dsum, wsum, min_out, max_out, miss_out, mask_out)
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

900    format (/,1x,10('-'),' output from horiz_interp ',10('-'))
901    format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
902    format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
903    format ('          number of missing points = ',i6)

    endif

    !-----------------------------------------------------------------------
  end subroutine horiz_interp_conserve
  ! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_conserve_end">

  !   <OVERVIEW>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_init.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_init.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_conserve_end ( Interp )
  !   </TEMPLATE>

  !   <INOUT NAME="Interp" TYPE="horiz_interp_type">
  !     A derived-type variable returned by previous call
  !     to horiz_interp_init. The input variable must have
  !     allocated arrays. The returned variable will contain
  !     deallocated arrays.
  !   </INOUT>

  subroutine horiz_interp_conserve_end ( Interp )

    type (horiz_interp_type), intent(inout) :: Interp

    deallocate ( Interp%area_src , Interp%area_dst , &
         Interp%facj, Interp%jlat, Interp%faci, Interp%ilon  )

  end subroutine horiz_interp_conserve_end
  ! </SUBROUTINE>

  !#######################################################################
  !---This statistics is for conservative scheme
  subroutine stats ( dat, area, asum, dsum, wsum, low, high, miss, mask )
    real,    intent(in)  :: dat(:,:), area(:,:)
    real,    intent(out) :: asum, dsum, wsum, low, high
    integer, intent(out) :: miss
    real,    intent(in), optional :: mask(:,:)

    integer :: pe, root_pe, npes, p, buffer_int(1)
    real    :: buffer_real(5)

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

    if(pe == root_pe) then
       do p = 1, npes - 1
          ! Force use of "scalar", integer pointer mpp interface
          call mpp_recv(buffer_real(1),glen=5,from_pe=root_pe+p)
          asum = asum + buffer_real(1)
          dsum = dsum + buffer_real(2)
          wsum = wsum + buffer_real(3)
          low  = min(low, buffer_real(4))
          high = max(high, buffer_real(5))
          call mpp_recv(buffer_int(1),glen=1,from_pe=root_pe+p)
          miss = miss + buffer_int(1)
       enddo
    else
       buffer_real(1) = asum
       buffer_real(2) = dsum
       buffer_real(3) = wsum
       buffer_real(4) = low
       buffer_real(5) = high
       ! Force use of "scalar", integer pointer mpp interface
       call mpp_send(buffer_real(1),plen=5,to_pe=root_pe)
       buffer_int(1) = miss
       call mpp_send(buffer_int(1),plen=1,to_pe=root_pe)
    endif

    call mpp_sync_self()   

  end subroutine stats

  !#######################################################################

  subroutine data_sum( data, area, facis, facie, facjs, facje,  &
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

end module horiz_interp_conserve_mod
