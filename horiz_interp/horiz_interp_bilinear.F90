module horiz_interp_bilinear_mod

  ! <CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang </CONTACT>

  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

  ! <OVERVIEW>
  !   Performs spatial interpolation between grids using bilinear interpolation
  ! </OVERVIEW>

  ! <DESCRIPTION>
  !     This module can interpolate data from regular rectangular grid
  !     to rectangular/tripolar grid. The interpolation scheme is bilinear interpolation.
  !     There is an optional mask field for missing input data.
  !     An optional output mask field may be used in conjunction with
  !     the input mask to show where output data exists.
  ! </DESCRIPTION>

  use mpp_mod,               only: mpp_error, FATAL, stdout, mpp_pe, mpp_root_pe
  use fms_mod,               only: write_version_number
  use constants_mod,         only: PI
  use horiz_interp_type_mod, only: horiz_interp_type, stats

  implicit none
  private


  public :: horiz_interp_bilinear_init, horiz_interp_bilinear, horiz_interp_bilinear_end

  real, parameter :: epsln=1.e-10
  integer         :: pe, root_pe
  !-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: horiz_interp_bilinear.F90,v 11.0 2004/09/28 19:59:42 fms Exp $'
  character(len=128) :: tagname = '$Name: khartoum $'
  logical            :: do_vers = .true.


contains

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_bilinear_init">

  !   <OVERVIEW>
  !      Initialization routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !      Allocates space and initializes a derived-type variable
  !      that contains pre-computed interpolation indices and weights.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_bilinear_init ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, src_modulo )

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

  !   <IN NAME="src_modulo" TYPE="logical, optional">
  !      logical variable to indicate if the boundary condition along zonal boundary
  !      is cyclic or not. When true, the zonal boundary condition is cyclic.
  !   </IN>

  !   <IN NAME="verbose" TYPE="integer, optional" >
  !      flag for the amount of print output.
  !   </IN>

  !   <INOUT NAME="Interp" TYPE="type(horiz_interp_type)" >
  !      A derived-type variable containing indices and weights used for subsequent 
  !      interpolations. To reinitialize this variable for a different grid-to-grid 
  !      interpolation you must first use the "horiz_interp_end" interface.
  !   </INOUT>
  subroutine horiz_interp_bilinear_init ( Interp, lon_in, lat_in, lon_out, lat_out, &
       verbose, src_modulo )

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout) :: Interp
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

    pe      = mpp_pe()
    root_pe = mpp_root_pe()

    if (do_vers) then
       call write_version_number (version, tagname)
       do_vers = .false.
    endif

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
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: when using bilinear interplation, '// &
         'the output grids should be geographical grids')    
    nlon_in = size(lon_in(:))  ; nlat_in = size(lat_in(:))
    nlon_out = size(lon_out, 1); nlat_out = size(lon_out, 2)
    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out

    if(lon_in(nlon_in) - lon_in(1) .gt. tpi + epsln) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: '// & 
         'The range of source grid longitude should be no larger than tpi')

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
                  call mpp_error(FATAL,'horiz_interp_bilinear_mod: ' //&
                  'when input grid is not modulo, output grid should locate inside input grid')
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
       write (stdout(),'(/,(1x,a))')                                      &
            '==> Warning: the geographic data set does not extend far   ', &
            '             enough east or west - a cyclic boundary       ', &
            '             condition was applied. check if appropriate   '
       write (stdout(),'(/,(1x,a,2f8.4))')                                &
            '    data required between longitudes:', gln_min, gln_max,     &
            '      data set is between longitudes:', lon_in(istart), lon_in(iend)
       warns = warns - 1
    endif
    !
    if (lt_err .eq. 1 .and. warns > 0) then
       write (stdout(),'(/,(1x,a))')                                     &
            '==> Warning: the geographic data set does not extend far   ',&
            '             enough north or south - extrapolation from    ',&
            '             the nearest data was applied. this may create ',&
            '             artificial gradients near a geographic pole   ' 
       write (stdout(),'(/,(1x,a,2f8.4))')                             &
            '    data required between latitudes:', glt_min, glt_max,   &
            '      data set is between latitudes:', lat_in(jstart), lat_in(jend)
       warns = warns - 1
    endif

    return

  end subroutine horiz_interp_bilinear_init
  ! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_bilinear">

  !   <OVERVIEW>
  !      Subroutine for performing the horizontal interpolation between two grids.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Subroutine for performing the horizontal interpolation between two grids. 
  !     horiz_interp_bilinear_init must be called before calling this routine.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_bilinear ( Interp, data_in, data_out, verbose, mask_in,mask_out, missing_value, missing_permit)
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
  !   <IN NAME="missing_value" TYPE="real, optional">
  !      Use the missing_value to indicate missing data.
  !   </IN>

  !   <IN NAME="missing_permit" TUPE="integer, optional">
  !      numbers of points allowed to miss for the bilinear interpolation. The value
  !      should be between 0 and 3.
  !   </IN>

  !   <OUT NAME="data_out" TYPE="real, dimension(:,:)">
  !      Output data on destination grid.
  !   </OUT>
  !   <OUT NAME="mask_out" TYPE="real, dimension(:,:),optional">
  !      Output mask that specifies whether data was computed.
  !   </OUT>

  subroutine horiz_interp_bilinear ( Interp, data_in, data_out, verbose, mask_in,mask_out, &
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
         miss_in, miss_out
    real    :: dwtsum, wtsum, min_in, max_in, avg_in, &
         min_out, max_out, avg_out, wtw, wte, wts, wtn
    real    :: mask(size(data_in,1), size(data_in,2) )

    num_missing = 0

    nlon_in  = Interp%nlon_src;  nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

    if(present(mask_in)) then
       mask = mask_in
    else
       mask = 1.0
    endif

    if (present(verbose)) then
       iverbose = verbose
    else
       iverbose = 0
    endif

    if(present(missing_permit)) then
       max_missing = missing_permit
    else
       max_missing = 0
    endif

    if(max_missing .gt. 3 .or. max_missing .lt. 0) call mpp_error(FATAL, &
         'horiz_interp_bilinear_mod: missing_permit should be between 0 and 3')

    if (size(data_in,1) /= nlon_in .or. size(data_in,2) /= nlat_in) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: size of input array incorrect')

    if (size(data_out,1) /= nlon_out .or. size(data_out,2) /= nlat_out) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: size of output array incorrect')

    do n = 1, nlat_out
       do m = 1, nlon_out
          is = Interp % i_lon (m,n,1); ie = Interp % i_lon (m,n,2)
          js = Interp % j_lat (m,n,1); je = Interp % j_lat (m,n,2)
          wtw = Interp % wti   (m,n,1)
          wte = Interp % wti   (m,n,2)
          wts = Interp % wtj   (m,n,1)
          wtn = Interp % wtj   (m,n,2)

          if(present(missing_value) ) then
             num_missing = 0
             if(data_in(is,js) == missing_value) then
                num_missing = num_missing+1
                mask(is,js) = 0.0
             endif
             if(data_in(ie,js) == missing_value) then
                num_missing = num_missing+1
                mask(ie,js) = 0.0
             endif
             if(data_in(ie,je) == missing_value) then
                num_missing = num_missing+1
                mask(ie,je) = 0.0
             endif
             if(data_in(is,je) == missing_value) then
                num_missing = num_missing+1
                mask(is,je) = 0.0
             endif
          endif

          dwtsum = data_in(is,js)*mask(is,js)*wtw*wts &
               + data_in(ie,js)*mask(ie,js)*wte*wts &
               + data_in(ie,je)*mask(ie,je)*wte*wtn &
               + data_in(is,je)*mask(is,je)*wtw*wtn 
          wtsum  = mask(is,js)*wtw*wts + mask(ie,js)*wte*wts  &
               + mask(ie,je)*wte*wtn + mask(is,je)*wtw*wtn

          !--- this will make sure the change will reproduce old results
          if(.not. present(mask_in) .and. .not. present(missing_value)) wtsum = 1.0

          if(num_missing .gt. max_missing ) then
             data_out(m,n) = missing_value
             if(present(mask_out)) mask_out(m,n) = 0.0
          else if(wtsum .lt. epsln) then 
             if(present(missing_value)) then
                data_out(m,n) = missing_value
             else
                data_out(m,n) = 0.0
             endif
             if(present(mask_out)) mask_out(m,n) = 0.0      
          else
             data_out(m,n) = dwtsum/wtsum
             if(present(mask_out)) mask_out(m,n) = wtsum
          endif
       enddo
    enddo
    !***********************************************************************
    ! compute statistics: minimum, maximum, and mean
    !-----------------------------------------------------------------------
    if (iverbose > 0) then

       ! compute statistics of input data

       call stats (data_in, min_in, max_in, avg_in, miss_in, missing_value, mask_in)

       ! compute statistics of output data
       call stats (data_out, min_out, max_out, avg_out, miss_out, missing_value, mask_out)

       !---- output statistics ----
       ! root_pe have the information of global mean, min and max
       if(pe == root_pe) then
          write (*,900)
          write (*,901)  min_in ,max_in, avg_in
          if (present(mask_in))  write (*,903)  miss_in
          write (*,902)  min_out,max_out,avg_out
          if (present(mask_out)) write (*,903)  miss_out
       endif
900    format (/,1x,10('-'),' output from horiz_interp ',10('-'))
901    format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
902    format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
903    format ('          number of missing points = ',i6)

    endif

    return

  end subroutine horiz_interp_bilinear
  ! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_bilinear_end">

  !   <OVERVIEW>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_init.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_init.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_bilinear_end ( Interp )
  !   </TEMPLATE>

  !   <INOUT NAME="Interp" TYPE="horiz_interp_type">
  !     A derived-type variable returned by previous call
  !     to horiz_interp_init. The input variable must have
  !     allocated arrays. The returned variable will contain
  !     deallocated arrays.
  !   </INOUT>

  subroutine horiz_interp_bilinear_end( Interp )

    type (horiz_interp_type), intent(inout) :: Interp

    deallocate (Interp%wti,   Interp%wtj, Interp%i_lon,    Interp%j_lat )


  end subroutine horiz_interp_bilinear_end
  ! </SUBROUTINE>

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
    ia = size(array(:))
    do i=2,ia
       if (array(i) .lt. array(i-1)) then
          write (stdout(),*) &
               ' => Error: array must be monotonically increasing in "indp"' , &
               '           when searching for nearest element to value=',value
          write (stdout(),*) '           array(i) < array(i-1) for i=',i 
          write (stdout(),*) '           array(i) for i=1..ia follows:'
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

  !######################################################################

end module horiz_interp_bilinear_mod
