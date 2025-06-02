!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @addtogroup horiz_interp_bilinear_mod
!> @{
  subroutine HORIZ_INTERP_BILINEAR_NEW_1D_ ( Interp, lon_in, lat_in, lon_out, lat_out, &
       verbose, src_modulo )

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout) :: Interp
    real(FMS_HI_KIND_), intent(in),  dimension(:)        :: lon_in , lat_in
    real(FMS_HI_KIND_), intent(in),  dimension(:,:)      :: lon_out, lat_out
    integer, intent(in),          optional :: verbose
    logical, intent(in),          optional :: src_modulo

    logical :: src_is_modulo
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, n, m
    integer :: ie, is, je, js, ln_err, lt_err, warns, iunit
    real(FMS_HI_KIND_)    :: wtw, wte, wts, wtn, lon, lat, tpi, hpi
    real(FMS_HI_KIND_)    :: glt_min, glt_max, gln_min, gln_max, min_lon, max_lon
    integer,parameter     :: kindl = FMS_HI_KIND_
    logical               :: decreasing_lat !< .True. if latitude is monotically decreasing
    logical               :: decreasing_lon !< .True. if longitude is monotically decreasing

    warns = 0
    if(present(verbose)) warns = verbose
    src_is_modulo = .true.
    if (present(src_modulo)) src_is_modulo = src_modulo

    decreasing_lat = .false.
    if (lat_in(1) > lat_in(2)) decreasing_lat = .true.

    decreasing_lon = .false.
    if (lon_in(1) > lon_in(2)) decreasing_lon = .true.

    hpi = 0.5_kindl * real(pi, FMS_HI_KIND_)
    tpi = 4.0_kindl * hpi
    glt_min = hpi
    glt_max = -hpi
    gln_min = tpi
    gln_max = -tpi
    min_lon = 0.0_kindl
    max_lon = tpi
    ln_err = 0
    lt_err = 0
    !-----------------------------------------------------------------------

    allocate ( Interp % HI_KIND_TYPE_ % wti (size(lon_out,1),size(lon_out,2),2),   &
               Interp % HI_KIND_TYPE_ % wtj (size(lon_out,1),size(lon_out,2),2),   &
               Interp % i_lon (size(lon_out,1),size(lon_out,2),2), &
               Interp % j_lat (size(lon_out,1),size(lon_out,2),2))
    !-----------------------------------------------------------------------

    nlon_in = size(lon_in(:))  ; nlat_in = size(lat_in(:))
    nlon_out = size(lon_out, 1); nlat_out = size(lon_out, 2)
    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out

    if(src_is_modulo) then
       if(lon_in(nlon_in) - lon_in(1) .gt. tpi + real(epsln, FMS_HI_KIND_)) &
            call mpp_error(FATAL,'horiz_interp_bilinear_mod: '// &
            'The range of source grid longitude should be no larger than tpi')

       if(lon_in(1) .lt. 0.0_kindl .OR. lon_in(nlon_in) > tpi ) then
          min_lon = lon_in(1)
          max_lon = lon_in(nlon_in)
       endif
    endif

    do n = 1, nlat_out
       do m = 1, nlon_out
          lon = lon_out(m,n)
          lat = lat_out(m,n)

          if(src_is_modulo) then
             if(lon .lt. min_lon) then
                lon = lon + tpi
             else if(lon .gt. max_lon) then
                lon = lon - tpi
             endif
          else  ! when the input grid is in not cyclic, the output grid should located inside
             ! the input grid
             if((lon .lt. lon_in(1)) .or. (lon .gt. lon_in(nlon_in))) &
                  call mpp_error(FATAL,'horiz_interp_bilinear_mod: ' //&
                  'when input grid is not modulo, output grid should locate inside input grid')
          endif

          glt_min = min(lat,glt_min);  glt_max = max(lat,glt_max)
          gln_min = min(lon,gln_min);  gln_max = max(lon,gln_max)

          is = nearest_index(lon, lon_in )
          if (decreasing_lon) then
            ! Lon_in is increasing
            ! This is so that is is the lower bound.
            ! For example, if the array is [50 40 30 20 10] and lon is 11, `is` is going to be 5 from `nearest_index`
            ! but it needs to be 4 so that it can use the data at lon_in(4) and lon_in(5)
            if( lon_in(is) .lt. lon ) is = max(is-1,1)
          else
            ! Lon_in is increasing
            ! This is so that is is the lower bound.
            ! For example, if the array is [10 20 30 40 50] and lon is 49, `is` is going to be 5 from `nearest_index`
            ! but it needs to be 4 so that it can use the data at lon_in(4) and lon_in(5)
            if( lon_in(is) .gt. lon ) is = max(is-1,1)
          endif
          if( lon_in(is) .eq. lon .and. is .eq. nlon_in) is = max(is - 1,1)
          ie = min(is+1,nlon_in)
          if(lon_in(is) .ne. lon_in(ie) .and. (decreasing_lon .or. lon_in(is) .le. lon)) then
             wtw = ( lon_in(ie) - lon) / (lon_in(ie) - lon_in(is) )
          else
             !     east or west of the last data value. this could be because a
             !     cyclic condition is needed or the dataset is too small.
             ln_err = 1
             ie = 1
             is = nlon_in
             if (lon_in(ie) .ge. lon ) then
                wtw = (lon_in(ie) -lon)/(lon_in(ie)-lon_in(is)+tpi+real(epsln,FMS_HI_KIND_))
             else
                wtw = (lon_in(ie)-lon+tpi+real(epsln,FMS_HI_KIND_))/(lon_in(ie)-lon_in(is)+tpi+real(epsln,FMS_HI_KIND_))
             endif
          endif
          wte = 1.0_kindl - wtw

          js = nearest_index(lat, lat_in )
          if (decreasing_lat) then
            ! Lat_in is decreasing
            if( lat_in(js) .lt. lat ) js = max(js - 1, 1)
          else
            ! Lat_in is increasing
            if( lat_in(js) .gt. lat ) js = max(js - 1, 1)
          endif
          if( lat_in(js) .eq. lat .and. js .eq. nlat_in) js = max(js - 1, 1)
          je = min(js + 1, nlat_in)

          if ( lat_in(js) .ne. lat_in(je) .and. (decreasing_lat .or. lat_in(js) .le. lat)) then
             wts = ( lat_in(je) - lat )/(lat_in(je)-lat_in(js))
          else
             !     north or south of the last data value. this could be because a
             !     pole is not included in the data set or the dataset is too small.
             !     in either case extrapolate north or south
             lt_err = 1
             wts = 1.0_kindl
          endif

          wtn = 1.0_kindl - wts

          Interp % i_lon (m,n,1) = is; Interp % i_lon (m,n,2) = ie
          Interp % j_lat (m,n,1) = js; Interp % j_lat (m,n,2) = je
          Interp % HI_KIND_TYPE_ % wti   (m,n,1) = wtw
          Interp % HI_KIND_TYPE_ % wti   (m,n,2) = wte
          Interp % HI_KIND_TYPE_ % wtj   (m,n,1) = wts
          Interp % HI_KIND_TYPE_ % wtj   (m,n,2) = wtn

       enddo
    enddo

    iunit = stdout()

    if (ln_err .eq. 1 .and. warns > 0) then
       write (iunit,'(/,(1x,a))')                                      &
            '==> Warning: the geographic data set does not extend far   ', &
            '             enough east or west - a cyclic boundary       ', &
            '             condition was applied. check if appropriate   '
       write (iunit,'(/,(1x,a,2f8.4))')                                &
            '    data required between longitudes:', gln_min, gln_max,     &
            '      data set is between longitudes:', lon_in(1), lon_in(nlon_in)
       warns = warns - 1
    endif

    if (lt_err .eq. 1 .and. warns > 0) then
       write (iunit,'(/,(1x,a))')                                     &
            '==> Warning: the geographic data set does not extend far   ',&
            '             enough north or south - extrapolation from    ',&
            '             the nearest data was applied. this may create ',&
            '             artificial gradients near a geographic pole   '
       write (iunit,'(/,(1x,a,2f8.4))')                             &
            '    data required between latitudes:', glt_min, glt_max,   &
            '      data set is between latitudes:', lat_in(1), lat_in(nlat_in)
    endif
    Interp% HI_KIND_TYPE_ % is_allocated = .true.
    Interp% interp_method = BILINEAR

    return

  end subroutine HORIZ_INTERP_BILINEAR_NEW_1D_

  !#######################################################################

  !> Initialization routine.
  !!
  !> Allocates space and initializes a derived-type variable
  !! that contains pre-computed interpolation indices and weights.
  subroutine HORIZ_INTERP_BILINEAR_NEW_2D_ ( Interp, lon_in, lat_in, lon_out, lat_out, &
       verbose, src_modulo, new_search, no_crash_when_not_found )

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout) :: Interp !< A derived type variable containing indices
                                          !! and weights for subsequent interpolations. To
                                          !! reinitialize for different grid-to-grid interpolation
                                          !! @ref horiz_interp_del must be used first.
    real(FMS_HI_KIND_), intent(in),  dimension(:,:)      :: lon_in !< Latitude (radians) for source data grid
    real(FMS_HI_KIND_), intent(in),  dimension(:,:)      :: lat_in !< Longitude (radians) for source data grid
    real(FMS_HI_KIND_), intent(in),  dimension(:,:)      :: lon_out !< Longitude (radians) for output data grid
    real(FMS_HI_KIND_), intent(in),  dimension(:,:)      :: lat_out !< Latitude (radians) for output data grid
    integer, intent(in),          optional :: verbose !< flag for amount of print output
    logical, intent(in),          optional :: src_modulo !< indicates if the boundary condition
                                          !! along zonal boundary is cyclic or not. Cyclic when true
    logical, intent(in),          optional :: new_search
    logical, intent(in),          optional :: no_crash_when_not_found
    integer                                :: warns
    logical                                :: src_is_modulo
    integer                                :: nlon_in, nlat_in, nlon_out, nlat_out
    integer                                :: m, n, is, ie, js, je, num_solution
    real(FMS_HI_KIND_)                                   :: lon, lat, quadra, x, y, y1, y2
    real(FMS_HI_KIND_)                                   :: a1, b1, c1, d1, a2, b2, c2, d2, a, b, c
    real(FMS_HI_KIND_)                                   :: lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4
    real(FMS_HI_KIND_)                                   :: tpi, lon_min, lon_max
    real(FMS_HI_KIND_)                                   :: epsln2
    logical                                :: use_new_search, no_crash

    integer, parameter :: kindl=FMS_HI_KIND_

    tpi = 2.0_kindl * real(pi, FMS_HI_KIND_)

    warns = 0
    if(present(verbose)) warns = verbose
    src_is_modulo = .true.
    if (present(src_modulo)) src_is_modulo = src_modulo
    use_new_search = .false.
    if (present(new_search)) use_new_search = new_search
    no_crash = .false.
    if(present(no_crash_when_not_found)) no_crash = no_crash_when_not_found

    ! make sure lon and lat has the same dimension
    if(size(lon_out,1) /= size(lat_out,1) .or. size(lon_out,2) /= size(lat_out,2) ) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: when using bilinear ' // &
         'interplation, the output grids should be geographical grids')

    if(size(lon_in,1) /= size(lat_in,1) .or. size(lon_in,2) /= size(lat_in,2) ) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: when using bilinear '// &
         'interplation, the input grids should be geographical grids')

    !--- get the grid size
    nlon_in  = size(lon_in,1) ; nlat_in  = size(lat_in,2)
    nlon_out = size(lon_out,1); nlat_out = size(lon_out,2)
    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out

    allocate ( Interp % HI_KIND_TYPE_ % wti (size(lon_out,1),size(lon_out,2),2),   &
               Interp % HI_KIND_TYPE_ % wtj (size(lon_out,1),size(lon_out,2),2),   &
               Interp % i_lon (size(lon_out,1),size(lon_out,2),2), &
               Interp % j_lat (size(lon_out,1),size(lon_out,2),2))

    !--- first fine the neighbor points for the destination points.
    if(use_new_search) then
       epsln2 = real(epsln,FMS_HI_KIND_)* 1.0e5_kindl
       call FIND_NEIGHBOR_NEW_(Interp, lon_in, lat_in, lon_out, lat_out, src_is_modulo, no_crash)
    else
       epsln2 = real(epsln,FMS_HI_KIND_)
       call FIND_NEIGHBOR_(Interp, lon_in, lat_in, lon_out, lat_out, src_is_modulo)
    endif

    !***************************************************************************
    !         Algorithm explanation (from disscussion with Steve Garner )      *
    !                                                                          *
    !    lon(x,y) = a1*x + b1*y + c1*x*y + d1         (1)                      *
    !    lat(x,y) = a2*x + b2*y + c2*x*y + d2         (2)                      *
    !    f (x,y) = a3*x + b3*y + c3*x*y + d3          (3)                      *
    !    with x and y is between 0 and 1.                                      *
    !    lon1 = lon(0,0) = d1,          lat1 = lat(0,0) = d2                   *
    !    lon2 = lon(1,0) = a1+d1,       lat2 = lat(1,0) = a2+d2                *
    !    lon3 = lon(1,1) = a1+b1+c1+d1, lat3 = lat(1,1) = a2+b2+c2+d2          *
    !    lon4 = lon(0,1) = b1+d1,       lat4 = lat(0,1) = b2+d2                *
    !    where (lon1,lat1),(lon2,lat2),(lon3,lat3),(lon4,lat4) represents      *
    !    the four corners starting from the left lower corner of grid box      *
    !    that encloses a destination grid ( the rotation direction is          *
    !    counterclockwise ). With these conditions, we get                     *
    !    a1 = lon2-lon1,           a2 = lat2-lat1                              *
    !    b1 = lon4-lon1,           b2 = lat4-lat1                              *
    !    c1 = lon3-lon2-lon4+lon1, c2 = lat3-lat2-lat4+lat1                    *
    !    d1 = lon1                 d2 = lat1                                   *
    !    So given any point (lon,lat), from equation (1) and (2) we can        *
    !    solve (x,y).                                                          *
    !    From equation (3)                                                     *
    !    f1 = f(0,0) = d3,          f2 = f(1,0) = a3+d3                        *
    !    f3 = f(1,1) = a3+b3+c3+d3, f4 = f(0,1) = b3+d3                        *
    !    we obtain                                                             *
    !    a3 = f2-f1,       b3 = f4-f1                                          *
    !    c3 = f3-f2-f4+f1, d3 = f1                                             *
    !    at point (lon,lat) ---> (x,y)                                         *
    !    f(x,y) = (f2-f1)x + (f4-f1)y + (f3-f2-f4+f1)xy + f1                   *
    !           = f1*(1-x)*(1-y) + f2*x*(1-y) + f3*x*y + f4*y*(1-x)            *
    !    wtw=1-x; wte=x; wts=1-y; xtn=y                                        *
    !                                                                          *
    !***************************************************************************

    lon_min = minval(lon_in);
    lon_max = maxval(lon_in);
    !--- calculate the weight
    do n = 1, nlat_out
       do m = 1, nlon_out
          lon = lon_out(m,n)
          lat = lat_out(m,n)
          if(lon .lt. lon_min) then
             lon = lon + tpi
          else if(lon .gt. lon_max) then
             lon = lon - tpi
          endif
          is = Interp%i_lon(m,n,1); ie = Interp%i_lon(m,n,2)
          js = Interp%j_lat(m,n,1); je = Interp%j_lat(m,n,2)
          if( is == DUMMY) cycle
          lon1 = lon_in(is,js); lat1 = lat_in(is,js);
          lon2 = lon_in(ie,js); lat2 = lat_in(ie,js);
          lon3 = lon_in(ie,je); lat3 = lat_in(ie,je);
          lon4 = lon_in(is,je); lat4 = lat_in(is,je);
          if(lon .lt. lon_min) then
             lon1 = lon1 -tpi; lon4 = lon4 - tpi
          else if(lon .gt. lon_max) then
             lon2 = lon2 +tpi; lon3 = lon3 + tpi
          endif
          a1 = lon2-lon1
          b1 = lon4-lon1
          c1 = lon1+lon3-lon4-lon2
          d1 = lon1
          a2 = lat2-lat1
          b2 = lat4-lat1
          c2 = lat1+lat3-lat4-lat2
          d2 = lat1
          !--- the coefficient of the quadratic equation
          a  = b2*c1-b1*c2
          b  = a1*b2-a2*b1+c1*d2-c2*d1+c2*lon-c1*lat
          c  = a2*lon-a1*lat+a1*d2-a2*d1
          quadra = b*b-4._kindl*a*c
          if(abs(quadra) < real(epsln, FMS_HI_KIND_)) quadra = 0.0_kindl
          if(quadra < 0.0_kindl) call mpp_error(FATAL, &
               "horiz_interp_bilinear_mod: No solution existed for this quadratic equation")
          if ( abs(a) .lt. epsln2) then  ! a = 0 is a linear equation
             if( abs(b) .lt. real(epsln,FMS_HI_KIND_)) call mpp_error(FATAL, &
                  "horiz_interp_bilinear_mod: no unique solution existed for this linear equation")
             y = -c/b
          else
             y1 = 0.5_kindl*(-b+sqrt(quadra))/a
             y2 = 0.5_kindl*(-b-sqrt(quadra))/a
             if(abs(y1) < epsln2) y1 = 0.0_kindl
             if(abs(y2) < epsln2) y2 = 0.0_kindl
             if(abs(1.0_kindl-y1) < epsln2) y1 = 1.0_kindl
             if(abs(1.0_kindl-y2) < epsln2) y2 = 1.0_kindl
             num_solution = 0
             if(y1 >= 0.0_kindl .and. y1 <= 1.0_kindl) then
                y = y1
                num_solution = num_solution +1
             endif
             if(y2 >= 0.0_kindl .and. y2 <= 1.0_kindl) then
                y = y2
                num_solution = num_solution + 1
             endif
             if(num_solution == 0) then
                call mpp_error(FATAL, "horiz_interp_bilinear_mod: No solution found")
             else if(num_solution == 2) then
                call mpp_error(FATAL, "horiz_interp_bilinear_mod: Two solutions found")
             endif
           endif
           if(abs(a1+c1*y) < real(epsln,FMS_HI_KIND_)) call mpp_error(FATAL, &
               "horiz_interp_bilinear_mod: the denomenator is 0")
           if(abs(y) < epsln2) y = 0.0_kindl
           if(abs(1.0_kindl-y) < epsln2) y = 1.0_kindl
           x = (lon-b1*y-d1)/(a1+c1*y)
           if(abs(x) < epsln2) x = 0.0_kindl
           if(abs(1.0_kindl-x) < epsln2) x = 1.0_kindl
           ! x and y should be between 0 and 1.
           !! Added for ECDA
           if(use_new_search) then
             if (x < 0.0_kindl) x = 0.0_kindl ! snz
             if (y < 0.0_kindl) y = 0.0_kindl ! snz
             if (x > 1.0_kindl) x = 1.0_kindl
             if (y > 1.0_kindl) y = 1.0_kindl
           endif
           if( x>1.0_kindl .or. x<0.0_kindl .or. y>1.0_kindl .or. y < 0.0_kindl) &
                  call mpp_error(FATAL, "horiz_interp_bilinear_mod: weight should be between 0 and 1")
           Interp % HI_KIND_TYPE_ % wti(m,n,1)=1.0_kindl-x
           Interp % HI_KIND_TYPE_ % wti(m,n,2)=x
           Interp % HI_KIND_TYPE_ % wtj(m,n,1)=1.0_kindl-y
           Interp % HI_KIND_TYPE_ % wtj(m,n,2)=y
       enddo
    enddo

    Interp% HI_KIND_TYPE_ % is_allocated = .true.
    Interp% interp_method = BILINEAR
  end subroutine

  !#######################################################################
  !> this routine will search the source grid to fine the grid box that encloses
  !! each destination grid.
  subroutine FIND_NEIGHBOR_ ( Interp, lon_in, lat_in, lon_out, lat_out, src_modulo )
    type(horiz_interp_type), intent(inout) :: Interp
    real(FMS_HI_KIND_), intent(in),       dimension(:,:) :: lon_in , lat_in
    real(FMS_HI_KIND_), intent(in),       dimension(:,:) :: lon_out, lat_out
    logical,                 intent(in)    :: src_modulo
    integer                                :: nlon_in, nlat_in, nlon_out, nlat_out
    integer                                :: max_step, n, m, l, i, j, ip1, jp1, step
    integer                                :: is, js, jstart, jend, istart, iend, npts
    integer, allocatable, dimension(:)     :: ilon, jlat
    real(FMS_HI_KIND_)                     :: lon_min, lon_max, lon, lat, tpi
    logical                                :: found
    real(FMS_HI_KIND_)                     :: lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4

    integer, parameter :: kindl=FMS_HI_KIND_

    tpi = 2.0_kindl*real(pi, FMS_HI_KIND_)
    nlon_in  = size(lon_in,1) ; nlat_in  = size(lat_in,2)
    nlon_out = size(lon_out,1); nlat_out = size(lon_out,2)

    lon_min = minval(lon_in);
    lon_max = maxval(lon_in);

    max_step = max(nlon_in,nlat_in) ! can be adjusted if needed
    allocate(ilon(5*max_step), jlat(5*max_step) )

    do n = 1, nlat_out
       do m = 1, nlon_out
          found = .false.
          lon = lon_out(m,n)
          lat = lat_out(m,n)

          if(src_modulo) then
             if(lon .lt. lon_min) then
                lon = lon + tpi
             else if(lon .gt. lon_max) then
                lon = lon - tpi
             endif
          else
             if(lon .lt. lon_min .or. lon .gt. lon_max ) &
             call mpp_error(FATAL,'horiz_interp_bilinear_mod: ' //&
                  'when input grid is not modulo, output grid should locate inside input grid')
          endif
          !--- search for the surrounding four points locatioon.
          if(m==1 .and. n==1) then
             J_LOOP: do j = 1, nlat_in-1
                do i = 1, nlon_in
                   ip1 = i+1
                   jp1 = j+1
                   if(i==nlon_in) then
                      if(src_modulo)then
                         ip1 = 1
                      else
                         cycle
                      endif
                   endif
                   lon1 = lon_in(i,  j);   lat1 = lat_in(i,j)
                   lon2 = lon_in(ip1,j);   lat2 = lat_in(ip1,j)
                   lon3 = lon_in(ip1,jp1); lat3 = lat_in(ip1,jp1)
                   lon4 = lon_in(i,  jp1); lat4 = lat_in(i,  jp1)

                   if(lon .lt. lon_min .or. lon .gt. lon_max) then
                      if(i .ne. nlon_in) then
                         cycle
                      else
                         if(lon .lt. lon_min) then
                             lon1 = lon1 -tpi; lon4 = lon4 - tpi
                         else if(lon .gt. lon_max) then
                             lon2 = lon2 +tpi; lon3 = lon3 + tpi
                         endif
                      endif
                   endif

                   if(lat .ge. intersect(lon1,lat1,lon2,lat2,lon))then ! south
                      if(lon .le. intersect(lat2,lon2,lat3,lon3,lat))then ! east
                         if(lat .le. intersect(lon3,lat3,lon4,lat4,lon))then ! north
                            if(lon .ge. intersect(lat4,lon4,lat1,lon1,lat))then  ! west
                               found = .true.
                               Interp % i_lon (m,n,1) = i; Interp % i_lon (m,n,2) = ip1
                               Interp % j_lat (m,n,1) = j; Interp % j_lat (m,n,2) = jp1
                               exit J_LOOP
                            endif
                         endif
                      endif
                   endif
                enddo
             enddo J_LOOP
          else
             step = 0
             do while ( .not. found .and. step .lt. max_step )
                !--- take the adajcent point as the starting point
                if(m == 1) then
                   is = Interp % i_lon (m,n-1,1)
                   js = Interp % j_lat (m,n-1,1)
                else
                   is = Interp % i_lon (m-1,n,1)
                   js = Interp % j_lat (m-1,n,1)
                endif
                if(step==0) then
                   npts = 1
                   ilon(1) = is
                   jlat(1) = js
                else
                   npts = 0
                   !--- bottom boundary
                   jstart = max(js-step,1)
                   jend   = min(js+step,nlat_in)

                   do l = -step, step
                      i = is+l
                      if(src_modulo)then
                         if( i < 1) then
                            i = i + nlon_in
                         else if (i > nlon_in) then
                            i = i - nlon_in
                         endif
                         if( i < 1 .or. i > nlon_in) call mpp_error(FATAL, &
                              'horiz_interp_bilinear_mod: max_step is too big, decrease max_step' )
                      else
                         if( i < 1 .or. i > nlon_in) cycle
                      endif

                      npts       = npts + 1
                      ilon(npts) = i
                      jlat(npts) = jstart
                   enddo

                   !--- right and left boundary -----------------------------------------------
                   istart = is - step
                   iend   = is + step
                   if(src_modulo) then
                      if( istart < 1)       istart = istart + nlon_in
                      if( iend   > nlon_in) iend   = iend   - nlon_in
                   else
                      istart = max(istart,1)
                      iend   = min(iend, nlon_in)
                   endif
                   do l = -step, step
                      j = js+l
                         if( j < 1 .or. j > nlat_in .or. j==jstart .or. j==jend) cycle
                         npts = npts+1
                         ilon(npts) = istart
                         jlat(npts) = j
                         npts = npts+1
                         ilon(npts) = iend
                         jlat(npts) = j
                  end do

                   !--- top boundary

                   do l = -step, step
                      i = is+l
                      if(src_modulo)then
                         if( i < 1) then
                            i = i + nlon_in
                         else if (i > nlon_in) then
                            i = i - nlon_in
                         endif
                         if( i < 1 .or. i > nlon_in) call mpp_error(FATAL, &
                              'horiz_interp_bilinear_mod: max_step is too big, decrease max_step' )
                      else
                         if( i < 1 .or. i > nlon_in) cycle
                      endif

                      npts       = npts + 1
                      ilon(npts) = i
                      jlat(npts) = jend
                   enddo


                end if

                !--- find the surrouding points
                do l = 1, npts
                   i = ilon(l)
                   j = jlat(l)
                   ip1 = i+1
                   if(ip1>nlon_in) then
                      if(src_modulo) then
                         ip1 = 1
                      else
                         cycle
                      endif
                   endif
                   jp1 = j+1
                   if(jp1>nlat_in) cycle
                   lon1 = lon_in(i,  j);   lat1 = lat_in(i,j)
                   lon2 = lon_in(ip1,j);   lat2 = lat_in(ip1,j)
                   lon3 = lon_in(ip1,jp1); lat3 = lat_in(ip1,jp1)
                   lon4 = lon_in(i,  jp1); lat4 = lat_in(i,  jp1)

                   if(lon .lt. lon_min .or. lon .gt. lon_max) then
                      if(i .ne. nlon_in) then
                         cycle
                      else
                         if(lon .lt. lon_min) then
                             lon1 = lon1 -tpi; lon4 = lon4 - tpi
                         else if(lon .gt. lon_max) then
                             lon2 = lon2 +tpi; lon3 = lon3 + tpi
                         endif
                      endif
                   endif

                   if(lat .ge. intersect(lon1,lat1,lon2,lat2,lon))then ! south
                      if(lon .le. intersect(lat2,lon2,lat3,lon3,lat))then ! east
                         if(lat .le. intersect(lon3,lat3,lon4,lat4,lon))then !north
                            if(lon .ge. intersect(lat4,lon4,lat1,lon1,lat))then ! west
                               found = .true.
                               is=i; js=j
                               Interp % i_lon (m,n,1) = i; Interp % i_lon (m,n,2) = ip1
                               Interp % j_lat (m,n,1) = j; Interp % j_lat (m,n,2) = jp1
                               exit
                            endif
                         endif
                      endif
                   endif
                enddo
                step = step + 1
             enddo
          endif
          if(.not.found) then
              print *,'lon,lat=',lon*180.0_kindl/real(PI,FMS_HI_KIND_),lat*180.0_kindl/real(PI,FMS_HI_KIND_)
              print *,'npts=',npts
              print *,'is,ie= ',istart,iend
              print *,'js,je= ',jstart,jend
              print *,'lon_in(is,js)=',lon_in(istart,jstart)*180.0_kindl/real(PI,FMS_HI_KIND_)
              print *,'lon_in(ie,js)=',lon_in(iend,jstart)*180.0_kindl/real(PI,FMS_HI_KIND_)
              print *,'lat_in(is,js)=',lat_in(istart,jstart)*180.0_kindl/real(PI,FMS_HI_KIND_)
              print *,'lat_in(ie,js)=',lat_in(iend,jstart)*180.0_kindl/real(PI,FMS_HI_KIND_)
              print *,'lon_in(is,je)=',lon_in(istart,jend)*180.0_kindl/real(PI,FMS_HI_KIND_)
              print *,'lon_in(ie,je)=',lon_in(iend,jend)*180.0_kindl/real(PI,FMS_HI_KIND_)
              print *,'lat_in(is,je)=',lat_in(istart,jend)*180.0_kindl/real(PI,FMS_HI_KIND_)
              print *,'lat_in(ie,je)=',lat_in(iend,jend)*180.0_kindl/real(PI,FMS_HI_KIND_)

             call mpp_error(FATAL, &
                  'FIND_NEIGHBOR_: the destination point is not inside the source grid' )
          endif
       enddo
    enddo

  end subroutine

  !#######################################################################

  !> The function will return true if the point x,y is inside a polygon, or
  !! false if it is not.  If the point is exactly on the edge of a polygon,
  !! the function will return .true.
  function INSIDE_POLYGON_(polyx, polyy, x, y)
     real(FMS_HI_KIND_), dimension(:), intent(in) :: polyx !< longitude coordinates of corners
     real(FMS_HI_KIND_), dimension(:), intent(in) :: polyy !< latitude coordinates of corners
     real(FMS_HI_KIND_),               intent(in) :: x !< x coordinate of point to be tested
     real(FMS_HI_KIND_),               intent(in) :: y !< y coordinate of point to be tested
     logical                        :: INSIDE_POLYGON_
     integer                        :: i, j, nedges
     real(FMS_HI_KIND_)                           :: xx

     INSIDE_POLYGON_ = .false.
     nedges = size(polyx(:))
     j = nedges
     do i = 1, nedges
        if( (polyy(i) < y .AND. polyy(j) >= y) .OR. (polyy(j) < y .AND. polyy(i) >= y) ) then
           xx = polyx(i)+(y-polyy(i))/(polyy(j)-polyy(i))*(polyx(j)-polyx(i))
           if( xx == x ) then
             INSIDE_POLYGON_ = .true.
             return
           else if( xx < x ) then
             INSIDE_POLYGON_ = .not. INSIDE_POLYGON_
           endif
        endif
        j = i
     enddo

     return

  end function

  !#######################################################################
  !> this routine will search the source grid to fine the grid box that encloses
  !! each destination grid.
  subroutine FIND_NEIGHBOR_NEW_( Interp, lon_in, lat_in, lon_out, lat_out, src_modulo, no_crash )
    type(horiz_interp_type), intent(inout) :: Interp
    real(FMS_HI_KIND_), intent(in),       dimension(:,:) :: lon_in , lat_in
    real(FMS_HI_KIND_), intent(in),       dimension(:,:) :: lon_out, lat_out
    logical,                 intent(in)    :: src_modulo, no_crash
    integer                                :: nlon_in, nlat_in, nlon_out, nlat_out
    integer                                :: max_step, n, m, l, i, j, ip1, jp1, step
    integer                                :: is, js, jstart, jend, istart, iend, npts
    integer, allocatable, dimension(:)     :: ilon, jlat
    real(FMS_HI_KIND_)                                   :: lon_min, lon_max, lon, lat, tpi
    logical                                :: found
    real(FMS_HI_KIND_)                                   :: polyx(4), polyy(4)
    real(FMS_HI_KIND_)                                   :: min_lon, min_lat, max_lon, max_lat

    integer, parameter :: step_div=8, kindl = FMS_HI_KIND_

    tpi = 2.0_kindl * real(pi, FMS_HI_KIND_)
    nlon_in  = size(lon_in,1) ; nlat_in  = size(lat_in,2)
    nlon_out = size(lon_out,1); nlat_out = size(lon_out,2)

    lon_min = minval(lon_in);
    lon_max = maxval(lon_in);

    max_step = min(nlon_in,nlat_in)/step_div ! can be adjusted if needed
    allocate(ilon(step_div*max_step), jlat(step_div*max_step) )

    do n = 1, nlat_out
       do m = 1, nlon_out
          found = .false.
          lon = lon_out(m,n)
          lat = lat_out(m,n)

          if(src_modulo) then
             if(lon .lt. lon_min) then
                lon = lon + tpi
             else if(lon .gt. lon_max) then
                lon = lon - tpi
             endif
          else
             if(lon .lt. lon_min .or. lon .gt. lon_max ) &
             call mpp_error(FATAL,'horiz_interp_bilinear_mod: ' //&
                  'when input grid is not modulo, output grid should locate inside input grid')
          endif
          !--- search for the surrounding four points locatioon.
          if(m==1 .and. n==1) then
             J_LOOP: do j = 1, nlat_in-1
                do i = 1, nlon_in
                   ip1 = i+1
                   jp1 = j+1
                   if(i==nlon_in) then
                      if(src_modulo)then
                         ip1 = 1
                      else
                         cycle
                      endif
                   endif

                   polyx(1) = lon_in(i,  j);   polyy(1) = lat_in(i,j)
                   polyx(2) = lon_in(ip1,j);   polyy(2) = lat_in(ip1,j)
                   polyx(3) = lon_in(ip1,jp1); polyy(3) = lat_in(ip1,jp1)
                   polyx(4) = lon_in(i,  jp1); polyy(4) = lat_in(i,  jp1)
                   if(lon .lt. lon_min .or. lon .gt. lon_max) then
                      if(i .ne. nlon_in) then
                         cycle
                      else
                         if(lon .lt. lon_min) then
                             polyx(1) = polyx(1) -tpi; polyx(4) = polyx(4) - tpi
                         else if(lon .gt. lon_max) then
                             polyx(2) = polyx(2) +tpi; polyx(3) = polyx(3) + tpi
                         endif
                      endif
                   endif

                   min_lon = minval(polyx)
                   max_lon = maxval(polyx)
                   min_lat = minval(polyy)
                   max_lat = maxval(polyy)
!                   if( lon .GE. min_lon .AND. lon .LE. max_lon .AND. &
!                       lat .GE. min_lat .AND. lat .LE. max_lat ) then
!                      print*, 'i =', i, 'j = ', j
!                      print '(5f15.11)', lon, polyx
!                      print '(5f15.11)', lat, polyy
!                   endif

                   if(INSIDE_POLYGON_(polyx, polyy, lon, lat)) then
                      found = .true.
!                      print*, " found ", i, j
                      Interp % i_lon (m,n,1) = i; Interp % i_lon (m,n,2) = ip1
                      Interp % j_lat (m,n,1) = j; Interp % j_lat (m,n,2) = jp1
                      exit J_LOOP
                   endif
                enddo
             enddo J_LOOP
          else
             step = 0
             do while ( .not. found .and. step .lt. max_step )
                !--- take the adajcent point as the starting point
                if(m == 1) then
                   is = Interp % i_lon (m,n-1,1)
                   js = Interp % j_lat (m,n-1,1)
                else
                   is = Interp % i_lon (m-1,n,1)
                   js = Interp % j_lat (m-1,n,1)
                endif
                if(step==0) then
                   npts = 1
                   ilon(1) = is
                   jlat(1) = js
                else
                   npts = 0
                   !--- bottom and top boundary
                   jstart = max(js-step,1)
                   jend   = min(js+step,nlat_in)

                   do l = -step, step
                      i = is+l
                      if(src_modulo)then
                         if( i < 1) then
                            i = i + nlon_in
                         else if (i > nlon_in) then
                            i = i - nlon_in
                         endif
                         if( i < 1 .or. i > nlon_in) call mpp_error(FATAL, &
                              'horiz_interp_bilinear_mod: max_step is too big, decrease max_step' )
                      else
                         if( i < 1 .or. i > nlon_in) cycle
                      endif

                      npts       = npts + 1
                      ilon(npts) = i
                      jlat(npts) = jstart
                      npts       = npts + 1
                      ilon(npts) = i
                      jlat(npts) = jend
                   enddo

                   !--- right and left boundary -----------------------------------------------
                   istart = is - step
                   iend   = is + step
                   if(src_modulo) then
                      if( istart < 1)       istart = istart + nlon_in
                      if( iend   > nlon_in) iend   = iend   - nlon_in
                   else
                      istart = max(istart,1)
                      iend   = min(iend, nlon_in)
                   endif
                   do l = -step, step
                      j = js+l
                         if( j < 1 .or. j > nlat_in) cycle
                         npts = npts+1
                         ilon(npts) = istart
                         jlat(npts) = j
                         npts = npts+1
                         ilon(npts) = iend
                         jlat(npts) = j
                  end do
                end if

                !--- find the surrouding points
                do l = 1, npts
                   i = ilon(l)
                   j = jlat(l)
                   ip1 = i+1
                   if(ip1>nlon_in) then
                      if(src_modulo) then
                         ip1 = 1
                      else
                         cycle
                      endif
                   endif
                   jp1 = j+1
                   if(jp1>nlat_in) cycle
                   polyx(1) = lon_in(i,  j);   polyy(1) = lat_in(i,j)
                   polyx(2) = lon_in(ip1,j);   polyy(2) = lat_in(ip1,j)
                   polyx(3) = lon_in(ip1,jp1); polyy(3) = lat_in(ip1,jp1)
                   polyx(4) = lon_in(i,  jp1); polyy(4) = lat_in(i,  jp1)
                   if(INSIDE_POLYGON_(polyx, polyy, lon, lat)) then
                      found = .true.
                      Interp % i_lon (m,n,1) = i; Interp % i_lon (m,n,2) = ip1
                      Interp % j_lat (m,n,1) = j; Interp % j_lat (m,n,2) = jp1
                      exit
                   endif
                enddo
                step = step + 1
             enddo
          endif
          if(.not.found) then
             if(no_crash) then
                Interp % i_lon (m,n,1:2) = DUMMY
                Interp % j_lat (m,n,1:2) = DUMMY
                print*,'lon,lat=',lon,lat ! snz
             else
                call mpp_error(FATAL, &
                    'horiz_interp_bilinear_mod: the destination point is not inside the source grid' )
             endif
          endif
       enddo
    enddo

  end subroutine

  !#######################################################################
  function INTERSECT_(x1, y1, x2, y2, x)
     real(FMS_HI_KIND_), intent(in) :: x1, y1, x2, y2, x
     real(FMS_HI_KIND_)             :: INTERSECT_

     INTERSECT_ = (y2-y1)*(x-x1)/(x2-x1) + y1

  return

  end function INTERSECT_

  !#######################################################################

  !> Subroutine for performing the horizontal interpolation between two grids
  !!
  !! @ref horiz_interp_bilinear_new must be called before calling this routine.
  subroutine HORIZ_INTERP_BILINEAR_ ( Interp, data_in, data_out, verbose, mask_in,mask_out, &
       missing_value, missing_permit, new_handle_missing )
    !-----------------------------------------------------------------------
    type (horiz_interp_type), intent(in)        :: Interp !< Derived type variable containing
                                               !! interpolation indices and weights. Returned by a
                                               !! previous call to horiz_interp_bilinear_new
    real(FMS_HI_KIND_), intent(in),  dimension(:,:)           :: data_in !< input data on source grid
    real(FMS_HI_KIND_), intent(out), dimension(:,:)           :: data_out !< output data on source grid
    integer, intent(in),               optional :: verbose !< 0 = no output; 1 = min,max,means; 2 =
                                                           !! all output
    real(FMS_HI_KIND_), intent(in), dimension(:,:),  optional :: mask_in !< Input mask, must be the same size as
                                    !! the input data. The real(FMS_HI_KIND_) value of mask_in must be in the
                                    !! range (0.,1.). Set mask_in=0.0 for data points
                                    !! that should not be used or have missing data
    real(FMS_HI_KIND_), intent(out), dimension(:,:), optional :: mask_out !< output mask that specifies whether
                                                            !! data was computed
    real(FMS_HI_KIND_), intent(in),                  optional :: missing_value
    integer, intent(in),               optional :: missing_permit
    logical, intent(in),               optional :: new_handle_missing
    !-----------------------------------------------------------------------
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, n, m,         &
         is, ie, js, je, iverbose, max_missing, num_missing, &
         miss_in, miss_out, iunit
    real(FMS_HI_KIND_)    :: dwtsum, wtsum, min_in, max_in, avg_in, &
         min_out, max_out, avg_out, wtw, wte, wts, wtn
    real(FMS_HI_KIND_)    :: mask(size(data_in,1), size(data_in,2) )
    logical :: set_to_missing, is_missing(4), new_handler
    real(FMS_HI_KIND_)    :: f1, f2, f3, f4, middle, w, s
    integer, parameter    :: kindl = FMS_HI_KIND_

    num_missing = 0

    nlon_in  = Interp%nlon_src;  nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

    if(present(mask_in)) then
       mask = mask_in
    else
       mask = 1.0_kindl
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

    if(present(new_handle_missing)) then
       new_handler = new_handle_missing
    else
       new_handler = .false.
    endif

    if(max_missing .gt. 3 .or. max_missing .lt. 0) call mpp_error(FATAL, &
         'horiz_interp_bilinear_mod: missing_permit should be between 0 and 3')

    if (size(data_in,1) /= nlon_in .or. size(data_in,2) /= nlat_in) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: size of input array incorrect')

    if (size(data_out,1) /= nlon_out .or. size(data_out,2) /= nlat_out) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: size of output array incorrect')

    if(new_handler) then
       if( .not. present(missing_value) )  call mpp_error(FATAL, &
            "horiz_interp_bilinear_mod: misisng_value must be present when new_handle_missing is .true.")
       if( present(mask_in) ) call mpp_error(FATAL, &
            "horiz_interp_bilinear_mod: mask_in should not be present when new_handle_missing is .true.")
       do n = 1, nlat_out
          do m = 1, nlon_out
             is = Interp % i_lon (m,n,1); ie = Interp % i_lon (m,n,2)
             js = Interp % j_lat (m,n,1); je = Interp % j_lat (m,n,2)
             wtw = Interp % HI_KIND_TYPE_ % wti   (m,n,1)
             wte = Interp % HI_KIND_TYPE_ % wti   (m,n,2)
             wts = Interp % HI_KIND_TYPE_ % wtj   (m,n,1)
             wtn = Interp % HI_KIND_TYPE_ % wtj   (m,n,2)

             is_missing = .false.
             num_missing = 0
             set_to_missing = .false.
             if(data_in(is,js) == missing_value) then
                num_missing = num_missing+1
                is_missing(1) = .true.
                if(wtw .GE. 0.5_kindl .AND. wts .GE. 0.5_kindl) set_to_missing = .true.
             endif

             if(data_in(ie,js) == missing_value) then
                num_missing = num_missing+1
                is_missing(2) = .true.
                if(wte .GE. 0.5_kindl .AND. wts .GE. 0.5_kindl ) set_to_missing = .true.
             endif
             if(data_in(ie,je) == missing_value) then
                num_missing = num_missing+1
                is_missing(3) = .true.
                if(wte .GE. 0.5_kindl .AND. wtn .GE. 0.5_kindl ) set_to_missing = .true.
             endif
             if(data_in(is,je) == missing_value) then
                num_missing = num_missing+1
                is_missing(4) = .true.
                if(wtw .GE. 0.5_kindl .AND. wtn .GE. 0.5_kindl) set_to_missing = .true.
             endif

             if( num_missing == 4 .OR. set_to_missing ) then
                data_out(m,n) = missing_value
                if(present(mask_out)) mask_out(m,n) = 0.0_kindl
                 cycle
             else if(num_missing == 0) then
                f1 = data_in(is,js)
                f2 = data_in(ie,js)
                f3 = data_in(ie,je)
                f4 = data_in(is,je)
                w = wtw
                s = wts
             else if(num_missing == 3) then  !--- three missing value
                if(.not. is_missing(1) ) then
                   data_out(m,n) = data_in(is,js)
                else if(.not. is_missing(2) ) then
                   data_out(m,n) = data_in(ie,js)
                else if(.not. is_missing(3) ) then
                   data_out(m,n) = data_in(ie,je)
                else if(.not. is_missing(4) ) then
                   data_out(m,n) = data_in(is,je)
                endif
                if(present(mask_out) ) mask_out(m,n) = 1.0_kindl
                cycle
             else   !--- one or two missing value
                if( num_missing == 1) then
                   if( is_missing(1) .OR. is_missing(3) ) then
                      middle = 0.5_kindl *(data_in(ie,js)+data_in(is,je))
                   else
                      middle = 0.5_kindl *(data_in(is,js)+data_in(ie,je))
                   endif
                else ! num_missing = 2
                   if( is_missing(1) .AND. is_missing(2) ) then
                      middle = 0.5_kindl *(data_in(ie,je)+data_in(is,je))
                   else if( is_missing(1) .AND. is_missing(3) ) then
                      middle = 0.5_kindl *(data_in(ie,js)+data_in(is,je))
                   else if( is_missing(1) .AND. is_missing(4) ) then
                      middle = 0.5_kindl *(data_in(ie,js)+data_in(ie,je))
                   else if( is_missing(2) .AND. is_missing(3) ) then
                      middle = 0.5_kindl *(data_in(is,js)+data_in(is,je))
                   else if( is_missing(2) .AND. is_missing(4) ) then
                      middle = 0.5_kindl*(data_in(is,js)+data_in(ie,je))
                   else if( is_missing(3) .AND. is_missing(4) ) then
                      middle = 0.5_kindl*(data_in(is,js)+data_in(ie,js))
                   endif
                endif

                if( wtw .GE. 0.5_kindl .AND. wts .GE. 0.5_kindl ) then  ! zone 1
                   w = 2.0_kindl*(wtw-0.5_kindl)
                   s = 2.0_kindl*(wts-0.5_kindl)
                   f1 = data_in(is,js)
                   if(is_missing(2)) then
                      f2 = f1
                   else
                      f2 = 0.5_kindl*(data_in(is,js)+data_in(ie,js))
                   endif
                   f3 = middle
                   if(is_missing(4)) then
                      f4 = f1
                   else
                      f4 = 0.5_kindl*(data_in(is,js)+data_in(is,je))
                   endif
                else if( wte .GE. 0.5_kindl .AND. wts .GE. 0.5_kindl ) then  ! zone 2
                   w = 2.0_kindl*(1.0_kindl-wte)
                   s = 2.0_kindl*(wts-0.5_kindl)
                   f2 = data_in(ie,js)
                   if(is_missing(1)) then
                      f1 = f2
                   else
                      f1 = 0.5_kindl*(data_in(is,js)+data_in(ie,js))
                   endif
                   f4 = middle
                   if(is_missing(3)) then
                      f3 = f2
                   else
                      f3 = 0.5_kindl*(data_in(ie,js)+data_in(ie,je))
                   endif
                else if( wte .GE. 0.5_kindl .AND. wtn .GE. 0.5_kindl ) then  ! zone 3
                   w = 2.0_kindl*(1.0_kindl-wte)
                   s = 2.0_kindl*(1.0_kindl-wtn)
                   f3 = data_in(ie,je)
                   if(is_missing(2)) then
                      f2 = f3
                   else
                      f2 = 0.5_kindl*(data_in(ie,js)+data_in(ie,je))
                   endif
                   f1 = middle
                   if(is_missing(4)) then
                      f4 = f3
                   else
                      f4 = 0.5_kindl*(data_in(ie,je)+data_in(is,je))
                   endif
                else if( wtw .GE. 0.5_kindl .AND. wtn .GE. 0.5_kindl ) then  ! zone 4
                   w = 2.0_kindl*(wtw-0.5_kindl)
                   s = 2.0_kindl*(1.0_kindl-wtn)
                   f4 = data_in(is,je)
                   if(is_missing(1)) then
                      f1 = f4
                   else
                      f1 = 0.5_kindl*(data_in(is,js)+data_in(is,je))
                   endif
                   f2 = middle
                   if(is_missing(3)) then
                      f3 = f4
                   else
                      f3 = 0.5_kindl*(data_in(ie,je)+data_in(is,je))
                   endif
                else
                   call mpp_error(FATAL, &
                      "horiz_interp_bilinear_mod: the point should be in one of the four zone")
                endif
             endif

             data_out(m,n) = f3 + (f4-f3)*w + (f2-f3)*s + ((f1-f2)+(f3-f4))*w*s
            if(present(mask_out)) mask_out(m,n) = 1.0_kindl
          enddo
       enddo
    else
       do n = 1, nlat_out
          do m = 1, nlon_out
             is = Interp % i_lon (m,n,1); ie = Interp % i_lon (m,n,2)
             js = Interp % j_lat (m,n,1); je = Interp % j_lat (m,n,2)
             wtw = Interp % HI_KIND_TYPE_ % wti   (m,n,1)
             wte = Interp % HI_KIND_TYPE_ % wti   (m,n,2)
             wts = Interp % HI_KIND_TYPE_ % wtj   (m,n,1)
             wtn = Interp % HI_KIND_TYPE_ % wtj   (m,n,2)

             if(present(missing_value) ) then
                num_missing = 0
                if(data_in(is,js) == missing_value) then
                   num_missing = num_missing+1
                   mask(is,js) = 0.0_kindl
                endif
                if(data_in(ie,js) == missing_value) then
                   num_missing = num_missing+1
                   mask(ie,js) = 0.0_kindl
                endif
                if(data_in(ie,je) == missing_value) then
                   num_missing = num_missing+1
                   mask(ie,je) = 0.0_kindl
                endif
                if(data_in(is,je) == missing_value) then
                   num_missing = num_missing+1
                   mask(is,je) = 0.0_kindl
                endif
             endif

             dwtsum = data_in(is,js)*mask(is,js)*wtw*wts &
                  + data_in(ie,js)*mask(ie,js)*wte*wts &
                  + data_in(ie,je)*mask(ie,je)*wte*wtn &
                  + data_in(is,je)*mask(is,je)*wtw*wtn
             wtsum  = mask(is,js)*wtw*wts + mask(ie,js)*wte*wts  &
                  + mask(ie,je)*wte*wtn + mask(is,je)*wtw*wtn

             if(.not. present(mask_in) .and. .not. present(missing_value)) wtsum = 1.0_kindl

             if(num_missing .gt. max_missing ) then
                data_out(m,n) = missing_value
                if(present(mask_out)) mask_out(m,n) = 0.0_kindl
             else if(wtsum .lt. real(epsln, FMS_HI_KIND_)) then
                if(present(missing_value)) then
                   data_out(m,n) = missing_value
                else
                   data_out(m,n) = 0.0_kindl
                endif
                if(present(mask_out)) mask_out(m,n) = 0.0_kindl
             else
                data_out(m,n) = dwtsum/wtsum
                if(present(mask_out)) mask_out(m,n) = wtsum
             endif
          enddo
       enddo
    endif
    !***********************************************************************
    ! compute statistics: minimum, maximum, and mean
    !-----------------------------------------------------------------------
    if (iverbose > 0) then

       ! compute statistics of input data

       call stats (data_in, min_in, max_in, avg_in, miss_in, missing_value, mask_in)

       ! compute statistics of output data
       call stats (data_out, min_out, max_out, avg_out, miss_out, missing_value, mask_out)

       !---- output statistics ----
       iunit = stdout()
       write (iunit,900)
       write (iunit,901)  min_in ,max_in, avg_in
       if (present(mask_in))  write (iunit,903)  miss_in
       write (iunit,902)  min_out,max_out,avg_out
       if (present(mask_out)) write (iunit,903)  miss_out

900    format (/,1x,10('-'),' output from horiz_interp ',10('-'))
901    format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
902    format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
903    format ('          number of missing points = ',i6)

    endif

    return

  end subroutine

  !> Subroutine for reading a weight file and use it to fill in the horiz interp type
  !! for the bilinear interpolation method.
  subroutine HORIZ_INTERP_READ_WEIGHTS_BILINEAR_(Interp, weight_filename, lon_out, lat_out, lon_in, lat_in, &
                                                 weight_file_source, interp_method, isw, iew, jsw, jew, nglon, nglat)
    type(horiz_interp_type), intent(inout) :: Interp             !< Horiz interp time to fill
    character(len=*),        intent(in)    :: weight_filename    !< Name of the weight file
    real(FMS_HI_KIND_), target,     intent(in)    :: lat_out(:,:)       !< Output (model) latitude
    real(FMS_HI_KIND_), target,      intent(in)    :: lon_out(:,:)       !< Output (model) longitude
    real(FMS_HI_KIND_),      intent(in)    :: lat_in(:)          !< Input (data) latitude
    real(FMS_HI_KIND_),      intent(in)    :: lon_in(:)          !< Input (data) longitude
    character(len=*),        intent(in)    :: weight_file_source !< Source of the weight file
    character(len=*),        intent(in)    :: interp_method      !< The interp method to use
    integer,                 intent(in)    :: isw, iew, jsw, jew !< Starting and ending indices of the compute domain
    integer,                 intent(in)    :: nglon              !< Number of longitudes in the global domain
    integer,                 intent(in)    :: nglat              !< Number of latitudes in the globl domain


    real(FMS_HI_KIND_), allocatable :: var(:,:,:)     !< Dummy variable to read the indices and weight into
    type(FmsNetcdfFile_t)           :: weight_fileobj !< FMS2io fileob for the weight file
    integer                         :: nlon           !< Number of longitudes in the model grid as read
                                                      !! from the weight file
    integer                         :: nlat           !< Number of latitude in the model grid as read
                                                      !! from the weight file

    if (.not. open_file(weight_fileobj, weight_filename, "read" )) &
      call mpp_error(FATAL, "Error opening the weight file:"//&
        &trim(weight_filename))

    !< Check that weight file has the correct dimensions
    select case (trim(weight_file_source))
    case ("fregrid")
      call get_dimension_size(weight_fileobj, "nlon", nlon)
      if (nlon .ne. nglon) &
        call mpp_error(FATAL, "The nlon from the weight file is not the same as in the input grid."//&
                             &" From weight file:"//string(nlon)//" from input grid:"//string(size(lon_out,1)))
      call get_dimension_size(weight_fileobj, "nlat", nlat)
      if (nlat .ne. nglat) &
        call mpp_error(FATAL, "The nlat from the weight file is not the same as in the input grid."//&
                   &" From weight file:"//string(nlat)//" from input grid:"//string(size(lon_out,2)))
    case default
      call mpp_error(FATAL, trim(weight_file_source)//&
        &" is not a supported weight file source. fregrid is the only supported weight file source." )
    end select

    Interp%nlon_src = size(lon_in(:)) ;  Interp%nlat_src = size(lat_in(:))
    Interp%nlon_dst = size(lon_out,1); Interp%nlat_dst = size(lon_out,2)

    allocate ( Interp % HI_KIND_TYPE_ % wti (Interp%nlon_dst,Interp%nlat_dst,2),   &
              Interp % HI_KIND_TYPE_ % wtj (Interp%nlon_dst,Interp%nlat_dst,2),   &
              Interp % i_lon (Interp%nlon_dst,Interp%nlat_dst,2), &
              Interp % j_lat (Interp%nlon_dst,Interp%nlat_dst,2))


    !! Three is for lon, lat, tile
    !! Currently, interpolation is only supported from lat,lon input data
    allocate(var(Interp%nlon_dst,Interp%nlat_dst, 3))
    call read_data(weight_fileobj, "index", var, corner=(/isw, jsw, 1/), edge_lengths=(/iew-isw+1, jew-jsw+1, 3/))

    !! Each point has a lon (i), and lat(j) index
    !! From there the four corners are (i,j), (i,j+1) (i+1) (i+1,j+1)
    Interp % i_lon (:,:,1) = var(:,:,1)
    Interp % i_lon (:,:,2) = Interp % i_lon (:,:,1) + 1
    where (Interp % i_lon (:,:,2) > size(lon_in(:))) Interp % i_lon (:,:,2) = 1

    Interp % j_lat (:,:,1) = var(:,:,2)
    Interp % j_lat (:,:,2) = Interp % j_lat (:,:,1) + 1
    where (Interp % j_lat (:,:,2) > size(lat_in(:))) Interp % j_lat (:,:,2) = 1

    deallocate(var)

    allocate(var(Interp%nlon_dst,Interp%nlat_dst, 4))
    call read_data(weight_fileobj, "weight", var, corner=(/isw, jsw, 1/), edge_lengths=(/iew-isw+1, jew-jsw+1, 4/))

    !! The weights for the four corners
    !! var(:,:,1) -> (i,j)
    !! var(:,:,2) -> (i,j+1)
    !! var(:,:,3) -> (i+1,j)
    !! var(:,:,4) -> (i+1,j+1)
    Interp % HI_KIND_TYPE_ % wti = var(:,:,1:2)
    Interp % HI_KIND_TYPE_ % wtj = var(:,:,3:4)
    deallocate(var)

    Interp% HI_KIND_TYPE_ % is_allocated = .true.
    Interp% interp_method = BILINEAR
    Interp% I_am_initialized = .True.
    call close_file(weight_fileobj)
  end subroutine HORIZ_INTERP_READ_WEIGHTS_BILINEAR_

!> @}
