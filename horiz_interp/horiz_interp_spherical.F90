module horiz_interp_spherical_mod

  ! <CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matthew Harrison </CONTACT>
  ! <CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang </CONTACT>

  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

  ! <OVERVIEW>
  !   Performs spatial interpolation between grids using inverse-distance-weighted scheme.
  ! </OVERVIEW>

  ! <DESCRIPTION>
  !     This module can interpolate data from rectangular/tripolar grid
  !     to rectangular/tripolar grid. The interpolation scheme is inverse-distance-weighted 
  !     scheme.    There is an optional mask field for missing input data.
  !     An optional output mask field may be used in conjunction with
  !     the input mask to show where output data exists.
  ! </DESCRIPTION>

  use mpp_mod,               only : mpp_error, FATAL, WARNING, stdout
  use mpp_mod,               only : mpp_root_pe, mpp_pe
  use fms_mod,               only : write_version_number
  use constants_mod,         only : pi
  use horiz_interp_type_mod, only : horiz_interp_type, stats

  implicit none
  private


  public :: horiz_interp_spherical_init, horiz_interp_spherical, horiz_interp_spherical_end

  real,    parameter :: max_dist_default = 0.17  ! radians
  integer, parameter :: num_nbrs_default = 4
  real,    parameter :: large=1.e20
  real,    parameter :: epsln=1.e-10

  integer            :: pe, root_pe
  logical            :: search_all = .false.

  !-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: horiz_interp_spherical.F90,v 11.0 2004/09/28 19:59:57 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'
  logical            :: do_vers = .true.
  logical            :: module_is_initialized = .FALSE.

contains

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_spherical_init">

  !   <OVERVIEW>
  !      Initialization routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !      Allocates space and initializes a derived-type variable
  !      that contains pre-computed interpolation indices and weights.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_spherical_init(Interp, lon_in,lat_in,lon_out,lat_out, num_nbrs, max_dist, src_modulo)
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

  !   <IN NAME="num_nbrs" TYPE="integer, optional">
  !     Number of nearest neighbors for regridding. When number of neighbors within
  !     the radius max_dist ( namelist variable) is less than num_nbrs, All the neighbors 
  !     will be used to interpolate onto destination grid. when number of neighbors within
  !     the radius max_dist ( namelist variable) is greater than num_nbrs, at least "num_nbrs"
  !     neighbors will be used to remap onto destination grid. 
  !   </IN>

  !   <IN NAME="max_dist" TYPE="real, optional" UNITS="radians">
  !      Maximum region of influence around destination grid points.
  !   </IN>

  !   <IN NAME="src_modulo" TYPE="logical, optional">
  !      logical variable to indicate if the boundary condition along zonal boundary
  !      is cyclic or not. When true, the zonal boundary condition is cyclic.
  !   </IN>

  !   <INOUT NAME="Interp" TYPE="type(horiz_interp_type)">
  !      A derived-type variable containing indices and weights used for subsequent 
  !      interpolations. To reinitialize this variable for a different grid-to-grid 
  !      interpolation you must first use the "horiz_interp_end" interface.
  !   </INOUT>

  subroutine horiz_interp_spherical_init(Interp, lon_in,lat_in,lon_out,lat_out, &
       num_nbrs, max_dist, src_modulo)
    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in),       dimension(:,:) :: lon_in, lat_in, lon_out, lat_out
    integer, intent(in),        optional   :: num_nbrs
    real, optional,             intent(in) :: max_dist
    logical,          intent(in), optional :: src_modulo

    !------local variables ---------------------------------------
    integer, parameter :: max_neighbors = 400 
    integer :: i, j, n, m, k
    integer :: map_dst_xsize, map_dst_ysize, map_src_xsize, map_src_ysize
    integer :: map_src_size, num_neighbors
    real    :: sum, max_src_dist
    integer, dimension(:), allocatable     :: ilon, jlat
    integer, dimension(:,:,:), allocatable :: map_src_add
    integer, dimension(:,:),   allocatable :: num_found
    real, dimension(:,:,:), allocatable    :: map_src_dist

    !--------------------------------------------------------------

    pe      = mpp_pe()
    root_pe = mpp_root_pe()

    if (do_vers) then
       call write_version_number (version, tagname)
       do_vers = .false.
    endif

    map_dst_xsize=size(lon_out,1);map_dst_ysize=size(lon_out,2)
    map_src_xsize=size(lon_in,1);map_src_ysize=size(lon_in,2)
    map_src_size = map_src_xsize*map_src_ysize

    num_neighbors = num_nbrs_default
    if(present(num_nbrs)) num_neighbors = num_nbrs

    max_src_dist = max_dist_default
    if (PRESENT(max_dist)) max_src_dist = max_dist
    Interp%max_src_dist = max_src_dist

    allocate(map_src_add(map_dst_xsize,map_dst_ysize,max_neighbors),    &
         map_src_dist(map_dst_xsize,map_dst_ysize,max_neighbors),   &
         num_found(map_dst_xsize,map_dst_ysize),                    &
         ilon(max_neighbors),jlat(max_neighbors)  )

    ! allocate memory to data type
    allocate(Interp%i_lon(map_dst_xsize,map_dst_ysize,max_neighbors), &
         Interp%j_lat(map_dst_xsize,map_dst_ysize,max_neighbors),     &
         Interp%src_dist(map_dst_xsize,map_dst_ysize,max_neighbors),  &
         Interp%num_found(map_dst_xsize,map_dst_ysize)       )

    map_src_add         = 0
    map_src_dist        = large
    num_found           = 0

    !using radial_search to find the nearest points and corresponding distance.

    if(search_all) then   ! normally full_search is false, set true to test the radial_search
       call full_search(lon_in, lat_in, lon_out, lat_out,map_src_add, map_src_dist, &
            num_found, num_nbrs, max_dist)
    else
       call radial_search(lon_in, lat_in, lon_out, lat_out,map_src_add, map_src_dist, &
            num_found, num_nbrs,max_dist,src_modulo)    
    endif

    do j=1,map_dst_ysize
       do i=1,map_dst_xsize
          do n=1,num_found(i,j)
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
          Interp%i_lon(i,j,:)         = ilon(:)
          Interp%j_lat(i,j,:)         = jlat(:)
          Interp%num_found(i,j)       = num_found(i,j)
          Interp%src_dist(i,j,:)      = map_src_dist(i,j,:)
       enddo
    enddo

    Interp%nlon_src = map_src_xsize; Interp%nlat_src = map_src_ysize
    Interp%nlon_dst = map_dst_xsize; Interp%nlat_dst = map_dst_ysize

    deallocate(map_src_add, map_src_dist, ilon, jlat)
    return

  end subroutine horiz_interp_spherical_init
  ! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_spherical">

  !   <OVERVIEW>
  !      Subroutine for performing the horizontal interpolation between two grids.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Subroutine for performing the horizontal interpolation between two grids. 
  !     horiz_interp_spherical_init must be called before calling this routine.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_spherical( Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value)
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
  !   <OUT NAME="data_out" TYPE="real, dimension(:,:)">
  !      Output data on destination grid.
  !   </OUT>
  !   <OUT NAME="mask_out" TYPE="real, dimension(:,:),optional">
  !      Output mask that specifies whether data was computed.
  !   </OUT>

  subroutine horiz_interp_spherical( Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value)
    type (horiz_interp_type), intent(in)        :: Interp
    real, intent(in),  dimension(:,:)           :: data_in
    real, intent(out), dimension(:,:)           :: data_out
    integer, intent(in),               optional :: verbose
    real, intent(in), dimension(:,:),  optional :: mask_in
    real, intent(out), dimension(:,:), optional :: mask_out
    real, intent(in),                  optional :: missing_value

    !--- some local variables ----------------------------------------
    real, dimension(Interp%nlon_dst, Interp%nlat_dst,size(Interp%src_dist,3)) :: wt
    real, dimension(Interp%nlon_src, Interp%nlat_src) :: mask_src
    real, dimension(Interp%nlon_dst, Interp%nlat_dst) :: mask_dst
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, num_found
    integer :: m, n, i, j, k, miss_in, miss_out, i1, i2, j1, j2, iverbose
    real    :: min_in, max_in, avg_in, min_out, max_out, avg_out, sum
    !-----------------------------------------------------------------

    iverbose = 0;  if (present(verbose)) iverbose = verbose

    nlon_in  = Interp%nlon_src; nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst   

    if(size(data_in,1) .ne. nlon_in .or. size(data_in,2) .ne. nlat_in ) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: size of input array incorrect')

    if(size(data_out,1) .ne. nlon_out .or. size(data_out,2) .ne. nlat_out ) & 
         call mpp_error(FATAL,'horiz_interp_spherical_mod: size of output array incorrect')

    mask_src = 1.0; mask_dst = 1.0
    if(present(mask_in)) mask_src = mask_in

    do n=1,nlat_out
       do m=1,nlon_out
          ! neighbors are sorted nearest to farthest
          ! check nearest to see if it is a land point
          num_found = Interp%num_found(m,n) 
          if(num_found == 0 ) then
             mask_dst(m,n) = 0.0
          else          
             i1 = Interp%i_lon(m,n,1); j1 = Interp%j_lat(m,n,1) 
             if (mask_src(i1,j1) .lt. 0.5) then
                mask_dst(m,n) = 0.0
             endif

             if(num_found .gt. 1 ) then
                i2 = Interp%i_lon(m,n,2); j2 = Interp%j_lat(m,n,2)
                ! compare first 2 nearest neighbors -- if they are nearly
                ! equidistant then use this mask for robustness
                if(abs(Interp%src_dist(m,n,2)-Interp%src_dist(m,n,1)) .lt. epsln) then
                   if((mask_src(i1,j1) .lt. 0.5))  mask_dst(m,n) = 0.0
                endif
             endif

             sum=0.0
             do k=1, num_found
                if(mask_src(Interp%i_lon(m,n,k),Interp%j_lat(m,n,k)) .lt. 0.5 ) then
                   wt(m,n,k) = 0.0
                else
                   if (Interp%src_dist(m,n,k) <= epsln) then
                      wt(m,n,k) = large
                      sum = sum + large
                   else if(Interp%src_dist(m,n,k) <= Interp%max_src_dist ) then
                      wt(m,n,k) = 1.0/Interp%src_dist(m,n,k)
                      sum = sum+wt(m,n,k)
                   else
                      wt(m,n,k) = 0.0
                   endif
                endif
             enddo
             if (sum > epsln) then
                do k = 1, num_found
                   wt(m,n,k) = wt(m,n,k)/sum
                enddo
             else
                mask_dst(m,n) = 0.0
             endif
          endif
       enddo
    enddo

    data_out = 0.0
    do n=1,nlat_out
       do m=1,nlon_out
          if(mask_dst(m,n) .gt. 0.5) then
             do k=1, Interp%num_found(m,n)
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

       call stats (data_in, min_in, max_in, avg_in, miss_in, missing_value, mask=mask_src)

       ! compute statistics of output data
       call stats (data_out, min_out, max_out, avg_out, miss_out, missing_value, mask=mask_dst)

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
  end subroutine horiz_interp_spherical
  ! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_spherical_end">

  !   <OVERVIEW>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_init.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_init.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_spherical_end ( Interp )
  !   </TEMPLATE>

  !   <INOUT NAME="Interp" TYPE="horiz_interp_type">
  !     A derived-type variable returned by previous call
  !     to horiz_interp_init. The input variable must have
  !     allocated arrays. The returned variable will contain
  !     deallocated arrays.
  !   </INOUT>


  subroutine horiz_interp_spherical_end( Interp )

    type (horiz_interp_type), intent(inout) :: Interp

    deallocate ( Interp%src_dist, Interp%num_found, Interp%i_lon, Interp%j_lat )

  end subroutine horiz_interp_spherical_end
  ! </SUBROUTINE>

  !#######################################################################


  subroutine radial_search(x_src,y_src,x_dst,y_dst, map_src_add, map_src_dist, &
       num_found, num_nbrs,max_dist,src_modulo)
    real,    intent(in),    dimension(:,:) :: x_src, y_src
    real,    intent(in),    dimension(:,:) :: x_dst, y_dst
    integer, intent(out), dimension(:,:,:) :: map_src_add
    real,    intent(out), dimension(:,:,:) :: map_src_dist
    integer, intent(inout), dimension(:,:) :: num_found
    integer, intent(in),          optional :: num_nbrs
    real,    intent(in),          optional :: max_dist
    logical, intent(in),          optional :: src_modulo

    !---------- local variables ----------------------------------------
    integer, parameter :: max_nbrs = 50
    real, dimension(size(x_src,1)*size(x_src,2)) :: theta_src,phi_src
    real, dimension(size(x_dst,1),size(x_dst,2)) :: theta_dst, phi_dst
    integer :: i, j, jj, i0, j0, n, l,i_left, j_left, i_right, j_right
    integer :: map_dst_xsize, map_dst_ysize, map_src_xsize, map_src_ysize
    integer :: i_left1, i_left2, i_right1, i_right2, num_neighbors
    integer :: map_src_size, map_dst_size, step, step_size, bound, bound_start, bound_end
    logical :: continue_search, result, continue_radial_search, src_is_modulo
    real    :: d, res, max_src_dist, tpi, hpi 
    real    :: min_theta_dst, max_theta_dst, min_phi_dst, max_phi_dst
    real    :: min_theta_src, max_theta_src, min_phi_src, max_phi_src 
    !------------------------------------------------------------------

    num_neighbors = num_nbrs_default
    if(present(num_nbrs)) num_neighbors = num_nbrs
    if (num_neighbors <= 0) call mpp_error(FATAL,'horiz_interp_spherical_mod: num_neighbors must be > 0') 

    max_src_dist = max_dist_default
    if (PRESENT(max_dist)) max_src_dist = max_dist

    src_is_modulo = .true.
    if (PRESENT(src_modulo)) src_is_modulo = src_modulo

    tpi = 2.0*PI; hpi = 0.5*PI

    map_dst_xsize=size(x_dst,1);map_dst_ysize=size(x_dst,2)
    map_src_xsize=size(x_src,1);map_src_ysize=size(x_src,2)
    map_dst_size = map_dst_xsize*map_dst_ysize
    map_src_size = map_src_xsize*map_src_ysize

    if (map_dst_size /= size(y_dst,1)*size(y_dst,2) .or. map_src_size /= size(y_src,1)*size(y_src,2)) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: grids not conformable')

    theta_src      = reshape(x_src,(/map_src_size/))
    phi_src        = reshape(y_src,(/map_src_size/))
    theta_dst(:,:) = x_dst(:,:)
    phi_dst(:,:)   = y_dst(:,:)

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
          continue_search=.true.
          step = 1
          step_size = sqrt(real(map_src_size) )
          do while (continue_search .and. step_size > 0)
             do while (step <= map_src_size .and. continue_search)
                ! count land points as nearest neighbors
                d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(step),phi_src(step))
                if (d <= max_src_dist) then
                   result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                        step,d,max_src_dist, num_found(i,j), num_neighbors )
                   if (result) then
                      n = 0
                      i0 = mod(step,map_src_xsize)

                      if (i0 == 0) i0 = map_src_xsize
                      res = float(step)/float(map_src_xsize)
                      j0 = ceiling(res)
                      continue_radial_search = .true.
                      do while (continue_radial_search)
                         continue_radial_search = .false.
                         n = n+1 ! radial counter
                         if(n > max_nbrs) exit
                         ! ************** left boundary *******************************
                         i_left = i0-n
                         if (i_left <= 0) then
                            if (src_is_modulo) then
                               i_left = map_src_xsize + i_left
                            else
                               i_left = 1
                            endif
                         endif

                         do l = 0, 2*n
                            jj = j0 - n - 1 + l
                            if( jj < 0) then
                               bound = ( 1 - jj )*map_src_xsize - i_left
                            else if ( jj >= map_src_ysize ) then
                               bound = ( 2*map_src_ysize - jj ) * map_src_xsize - i_left
                            else
                               bound = jj * map_src_xsize + i_left
                            endif

                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                            result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                 bound,d,max_src_dist, num_found(i,j), num_neighbors)

                            if (result) continue_radial_search = .true.
                         enddo

                         ! ***************************right boundary ******************************* 
                         i_right = i0+n
                         if (i_right > map_src_xsize) then
                            if (src_is_modulo) then
                               i_right = i_right - map_src_xsize
                            else
                               i_right = map_src_xsize
                            endif
                         endif

                         do l = 0, 2*n
                            jj = j0 - n - 1 + l
                            if( jj < 0) then
                               bound = ( 1 - jj )*map_src_xsize - i_right
                            else if ( jj >= map_src_ysize ) then
                               bound = ( 2*map_src_ysize - jj) * map_src_xsize - i_right

                            else
                               bound = jj * map_src_xsize + i_right
                            endif

                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                            result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                 bound,d,max_src_dist, num_found(i,j), num_neighbors)

                            if (result) continue_radial_search = .true.
                         enddo

                         ! ************************* bottom boundary **********************************
                         i_left2 = 0
                         if( i_left > i_right) then
                            i_left1 = 1
                            i_right1 = i_right
                            i_left2 = i_left
                            i_right2 = map_src_xsize
                         else
                            i_left1 = i_left
                            i_right1 = i_right                            
                         endif

                         jj = j0 - n - 1
                         if( jj < 0 ) then
                            bound_start = ( 1 - jj)*map_src_xsize - i_right1
                            bound_end   = ( 1 - jj)*map_src_xsize - i_left1
                         else
                            bound_start = jj * map_src_xsize + i_left1
                            bound_end = jj * map_src_xsize + i_right1
                         endif

                         bound = bound_start
                         do while (bound <= bound_end)
                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                            result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                 bound,d,max_src_dist, num_found(i,j), num_neighbors)
                            bound = bound + 1
                            if (result) continue_radial_search = .true.
                         enddo

                         if(i_left2 > 0 ) then
                            if( jj < 0 ) then
                               bound_start = ( 1 - jj)*map_src_xsize - i_right2
                               bound_end   = ( 1 - jj)*map_src_xsize - i_left2
                            else
                               bound_start = jj * map_src_xsize + i_left2
                               bound_end = jj * map_src_xsize + i_right2
                            endif

                            bound = bound_start
                            do while (bound <= bound_end)
                               d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                               result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d,max_src_dist, num_found(i,j), num_neighbors)
                               bound = bound + 1

                               if (result) continue_radial_search = .true.
                            enddo
                         endif

                         ! ************************** top boundary ************************************
                         jj = j0 + n - 1
                         if( jj >= map_src_ysize) then
                            bound_start = ( 2*map_src_ysize - jj ) * map_src_xsize - i_right1
                            bound_end   = ( 2*map_src_ysize - jj ) * map_src_xsize - i_left1
                         else
                            bound_start = jj * map_src_xsize + i_left1
                            bound_end = jj * map_src_xsize + i_right1
                         endif

                         bound = bound_start
                         do while (bound <= bound_end)
                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                            result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                 bound,d,max_src_dist, num_found(i,j), num_neighbors)
                            bound = bound + 1
                            if (result) continue_radial_search = .true.
                         enddo

                         if(i_left2 > 0) then
                            if( jj >= map_src_ysize) then
                               bound_start = ( 2*map_src_ysize - jj ) * map_src_xsize - i_right2
                               bound_end   = ( 2*map_src_ysize - jj ) * map_src_xsize - i_left2
                            else
                               bound_start = jj * map_src_xsize + i_left2
                               bound_end = jj * map_src_xsize + i_right2
                            endif

                            bound = bound_start
                            do while (bound <= bound_end)
                               d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound),phi_src(bound))
                               result = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d,max_src_dist, num_found(i,j), num_neighbors)
                               bound = bound + 1
                               if (result) continue_radial_search = .true.
                            enddo
                         endif

                      enddo
                      continue_search = .false. ! stop looking
                   endif
                endif
                step=step+step_size
             enddo ! search loop
             step = 1
             step_size = step_size/2
          enddo
       enddo
    enddo

    return

  end subroutine radial_search


  !#####################################################################

  function update_dest_neighbors(map_src_add, map_src_dist, src_add,d, max_src_dist, num_found, min_nbrs)

    integer, intent(inout), dimension(:) :: map_src_add
    real, intent(inout),    dimension(:) :: map_src_dist
    integer, intent(in)                  :: src_add
    real, intent(in)                     :: d
    real, intent(in)                     :: max_src_dist
    integer, intent(inout)               :: num_found
    integer, intent(in)                  :: min_nbrs

    logical :: update_dest_neighbors, already_exist = .false.

    integer :: n,m, max_neighbors

    update_dest_neighbors = .false.
    max_neighbors = size(map_src_add) 

    if (d .le. max_src_dist) then
       n = 0
       NLOOP : do while ( n .le. num_found )
          n = n + 1
          DIST_CHK : if (d .le. map_src_dist(n)) then
             do m=n,num_found
                if (src_add == map_src_add(m)) then
                   already_exist = .true.
                   exit NLOOP
                endif
             enddo
             if(num_found < max_neighbors) then
                num_found = num_found + 1
             else
                call mpp_error(FATAL,'update_dest_neighbors: '// &
                     'number of neighbor points found is greated than maxium neighbor points' )
             endif
             do m=num_found,n+1,-1
                map_src_add(m) = map_src_add(m-1)
                map_src_dist(m) = map_src_dist(m-1)
             enddo
             map_src_add(n) = src_add
             map_src_dist(n) = d
             update_dest_neighbors = .true.
             if( num_found > min_nbrs ) then
                if( map_src_dist(num_found) > map_src_dist(num_found-1) ) then
                   num_found = num_found - 1
                endif
                if( map_src_dist(min_nbrs+1) > map_src_dist(min_nbrs) ) then
                   num_found = min_nbrs
                endif
             endif
             exit NLOOP ! n loop
          endif DIST_CHK
       end do NLOOP
       if(already_exist) return

       if( .not. update_dest_neighbors ) then
          if( num_found < min_nbrs ) then
             num_found               = num_found + 1
             update_dest_neighbors   = .true.
             map_src_add(num_found)  = src_add
             map_src_dist(num_found) = d
          endif
       endif

    endif

    return

  end function update_dest_neighbors

  !########################################################################
!  function spherical_distance(theta1,phi1,theta2,phi2)

!    real, intent(in) :: theta1, phi1, theta2, phi2
!    real :: spherical_distance

!    real :: r1(3), r2(3), cross(3), s, dot, ang

    ! this is a simple, enough way to calculate distance on the sphere
    ! first, construct cartesian vectors r1 and r2
    ! then calculate the cross-product which is proportional to the area
    ! between the 2 vectors.  The angular distance is arcsin of the 
    ! distancealong the sphere
    !
    ! theta is longitude and phi is latitude
    !


!    r1(1) = cos(theta1)*cos(phi1);r1(2)=sin(theta1)*cos(phi1);r1(3)=sin(phi1)
!    r2(1) = cos(theta2)*cos(phi2);r2(2)=sin(theta2)*cos(phi2);r2(3)=sin(phi2)

!    cross(1) = r1(2)*r2(3)-r1(3)*r2(2)
!    cross(2) = r1(3)*r2(1)-r1(1)*r2(3)
!    cross(3) = r1(1)*r2(2)-r1(2)*r2(1)

!    s = sqrt(cross(1)**2.+cross(2)**2.+cross(3)**2.)

!    s = min(s,1.0-epsln)

!    dot = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)

!    if (dot > 0) then
!       ang = asin(s)
!    else if (dot < 0) then
!       ang = pi + asin(s)  !?  original is pi - asin(s)
!    else
!       ang = pi/2.
!    endif

!    spherical_distance = abs(ang) ! in radians

!    return

!  end function spherical_distance
  ! The great cycle distance
  function spherical_distance(theta1,phi1,theta2,phi2)

    real, intent(in) :: theta1, phi1, theta2, phi2
    real :: spherical_distance, dot

    if(theta1 == theta2 .and. phi1 == phi2) then
        spherical_distance = 0.0
        return
    endif
  
    dot = cos(phi1)*cos(phi2)*cos(theta1-theta2) + sin(phi1)*sin(phi2)
    if(dot > 1 ) dot = 1.
    if(dot < -1) dot = -1
    spherical_distance = acos(dot)

    return

  end function spherical_distance


  !#######################################################################

  subroutine full_search(x_src,y_src,x_dst,y_dst,map_src_add, map_src_dist,num_found, num_nbrs,max_dist)
    real, intent(in), dimension(:,:)       :: x_src, y_src, x_dst, y_dst
    integer, intent(out), dimension(:,:,:) :: map_src_add
    real,    intent(out), dimension(:,:,:) :: map_src_dist
    integer, intent(out), dimension(:,:)   :: num_found
    integer, intent(in),          optional :: num_nbrs
    real, intent(in),             optional :: max_dist

    integer :: i,j,num_neighbors, map_src_size, step
    integer :: map_dst_xsize,map_dst_ysize, map_src_xsize,map_src_ysize 
    real    :: max_src_dist, tpi, hpi, d
    logical :: found

    real, dimension(size(x_dst,1),size(x_dst,2)) :: theta_dst, phi_dst
    real, dimension(size(x_src,1)*size(x_src,2)) :: theta_src, phi_src


    tpi = 2.0*PI; hpi = 0.5*PI

    map_dst_xsize=size(x_dst,1);map_dst_ysize=size(x_dst,2)
    map_src_xsize=size(x_src,1);map_src_ysize=size(x_src,2)
    map_src_size = map_src_xsize*map_src_ysize

    num_neighbors = num_nbrs_default
    if(present(num_nbrs)) num_neighbors = num_nbrs
    if (num_neighbors <= 0) call mpp_error(FATAL,'horiz_interp_spherical: num_neighbors must be > 0') 

    max_src_dist = max_dist_default
    if (PRESENT(max_dist)) max_src_dist = max_dist

    theta_dst = x_dst; phi_dst = y_dst
    theta_src = reshape(x_src,(/map_src_size/))
    phi_src   = reshape(y_src,(/map_src_size/))

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
          do step = 1, map_src_size
             d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(step),phi_src(step))
             found = update_dest_neighbors(map_src_add(i,j,:),map_src_dist(i,j,:), &
                  step,d,max_src_dist, num_found(i,j), num_neighbors )
          enddo
       enddo
    enddo

  end subroutine full_search

  !#######################################################################


end module horiz_interp_spherical_mod
