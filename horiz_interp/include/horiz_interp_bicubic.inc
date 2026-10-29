!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
  !> @parblock
  !! Computes weights and mapping indices for bicubic interpation from a rectilinear input grid
  !! specified as 1D arrays to a general output grid specified as 2D arrays.
  !! @endparblock
  subroutine HORIZ_INTERP_BICUBIC_NEW_1D_S_ ( Interp, lon_in, lat_in, lon_out, lat_out, &
       verbose, src_modulo )

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout) :: Interp
      !< will be populated with interpolation weights and mapping indices
    real(FMS_HI_KIND_), intent(in), dimension(:) :: lon_in
      !< are the longitude coordinates for the rectilinear input grid in radians
    real(FMS_HI_KIND_), intent(in), dimension(:) :: lat_in
      !< are the latitude coordinates for the rectilinear input grid in radians
    real(FMS_HI_KIND_), intent(in), dimension(:,:) :: lon_out
      !< are the longitude coordinates for the output grid in radians
    real(FMS_HI_KIND_), intent(in), dimension(:,:) :: lat_out
      !< are the latitude coordinates for the output grid in radians
    integer, intent(in), optional :: verbose
      !< sets the verbosity level for output
    logical, intent(in), optional :: src_modulo
      !< is true if boundary condition along zonal boundary is cyclic
    integer :: i, j, ip1, im1, jp1, jm1
    integer :: nlon_in, nlat_in, nlon_out, nlat_out
    integer :: jcl, jcu, icl, icu, jj
    real(FMS_HI_KIND_) :: xz, yz
    integer :: iunit
    integer, parameter :: kindl = FMS_HI_KIND_
      !< real size at compile time

    if(present(verbose)) verbose_bicubic = verbose
    ! Retain the historical periodic longitude behavior even when src_modulo is false.

    if(size(lon_out,1) /= size(lat_out,1) .or. size(lon_out,2) /= size(lat_out,2) ) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: when using bilinear ' // &
         'interpolation, the output grids should be geographical grids')

    !--- get the grid size
    nlon_in  = size(lon_in)   ; nlat_in  = size(lat_in)
    nlon_out = size(lon_out,1); nlat_out = size(lat_out,2)
    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out
!   use wti(:,:,1) for x-derivative, wti(:,:,2) for y-derivative, wti(:,:,3) for xy-derivative
    allocate ( Interp%HI_KIND_TYPE_%wti    (nlon_in, nlat_in, 3) )
    allocate ( Interp%HI_KIND_TYPE_%lon_in (nlon_in) )
    allocate ( Interp%HI_KIND_TYPE_%lat_in (nlat_in) )
    allocate ( Interp%HI_KIND_TYPE_%rat_x  (nlon_out, nlat_out) )
    allocate ( Interp%HI_KIND_TYPE_%rat_y  (nlon_out, nlat_out) )
    allocate ( Interp%i_lon  (nlon_out, nlat_out, 2) )
    allocate ( Interp%j_lat  (nlon_out, nlat_out, 2) )

    Interp%HI_KIND_TYPE_%lon_in = lon_in
    Interp%HI_KIND_TYPE_%lat_in = lat_in

    if ( verbose_bicubic > 0 ) then
       iunit = stdout()
       write (iunit,'(/,"Initialising bicubic interpolation, interface horiz_interp_bicubic_new_1d_s")')
       write (iunit,'(/," Longitude of coarse grid points (radian): xc(i) i=1, ",i4)') Interp%nlon_src
       write (iunit,'(1x,10f10.4)') (Interp%HI_KIND_TYPE_%lon_in(jj),jj=1,Interp%nlon_src)
       write (iunit,'(/," Latitude of coarse grid points (radian):  yc(j) j=1, ",i4)') Interp%nlat_src
       write (iunit,'(1x,10f10.4)') (Interp%HI_KIND_TYPE_%lat_in(jj),jj=1,Interp%nlat_src)
       do i=1, Interp%nlat_dst
         write (iunit,*)
         write (iunit,'(/," Longitude of fine grid points (radian): xf(i) i=1, ",i4)') Interp%nlat_dst
         write (iunit,'(1x,10f10.4)') (lon_out(jj,i),jj=1,Interp%nlon_dst)
       enddo
       do i=1, Interp%nlon_dst
         write (iunit,*)
         write (iunit,'(/," Latitude of fine grid points (radian):  yf(j) j=1, ",i4)') Interp%nlon_dst
         write (iunit,'(1x,10f10.4)') (lat_out(i,jj),jj=1,Interp%nlat_dst)
       enddo
    endif


!---------------------------------------------------------------------------
!     Find the x-derivative. Use central differences and forward or
!     backward steps at the boundaries

    do j=1,nlat_in
      do i=1,nlon_in
        ip1=min(i+1,nlon_in)
        im1=max(i-1,1)
        Interp%HI_KIND_TYPE_%wti(i,j,1) = 1.0_kindl/(Interp%HI_KIND_TYPE_%lon_in(ip1)-Interp%HI_KIND_TYPE_%lon_in(im1))
      enddo
    enddo


!---------------------------------------------------------------------------

!     Find the y-derivative. Use central differences and forward or
!     backward steps at the boundaries
      do j=1,nlat_in
        jp1=min(j+1,nlat_in)
        jm1=max(j-1,1)
        do i=1,nlon_in
          Interp%HI_KIND_TYPE_%wti(i,j,2) =1.0_kindl/(Interp%HI_KIND_TYPE_%lat_in(jp1)-Interp%HI_KIND_TYPE_%lat_in(jm1))
        enddo
      enddo

!---------------------------------------------------------------------------

!     Find the xy-derivative. Use central differences and forward or
!     backward steps at the boundaries
      do j=1,nlat_in
        jp1=min(j+1,nlat_in)
        jm1=max(j-1,1)
        do i=1,nlon_in
          ip1=min(i+1,nlon_in)
          im1=max(i-1,1)
          Interp%HI_KIND_TYPE_%wti(i,j,3) = 1.0_kindl / &
                                              ((Interp%HI_KIND_TYPE_%lon_in(ip1)-Interp%HI_KIND_TYPE_%lon_in(im1)) * &
                                               (Interp%HI_KIND_TYPE_%lat_in(jp1)-Interp%HI_KIND_TYPE_%lat_in(jm1)))
        enddo
      enddo
!---------------------------------------------------------------------------
!     Now for each point at the dest-grid find the boundary points of
!     the source grid
      do j=1, nlat_out
        do i=1,nlon_out
          yz  = lat_out(i,j)
          xz  = lon_out(i,j)

          call bicubic_interval(lat_in, yz, .false., jcl, jcu, Interp%HI_KIND_TYPE_%rat_y(i,j))
          call bicubic_interval(lon_in, xz, .true., icl, icu, Interp%HI_KIND_TYPE_%rat_x(i,j))
          Interp%j_lat(i,j,1) = jcl
          Interp%j_lat(i,j,2) = jcu
          Interp%i_lon(i,j,1) = icl
          Interp%i_lon(i,j,2) = icu
        enddo
      enddo
    Interp% HI_KIND_TYPE_ % is_allocated = .true.
    Interp%interp_method = BICUBIC
  end subroutine HORIZ_INTERP_BICUBIC_NEW_1D_S_

  !> @parblock
  !! Subroutine to compute weignts and mapping indices for bicubic interpolation
  !! from a rectilinear input grid to a rectilinear output grid, both supplied
  !! as 1D arrays of latitude and longitude gridpoints.
  !! @endparblock
  subroutine HORIZ_INTERP_BICUBIC_NEW_1D_ ( Interp, lon_in, lat_in, lon_out, lat_out, &
       verbose, src_modulo )

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout) :: Interp
      !< will be populated with interpolation weights and mapping indices
    real(FMS_HI_KIND_), intent(in), dimension(:) :: lon_in
      !< are the longitude coordinates for the rectilinear input grid in radians
    real(FMS_HI_KIND_), intent(in), dimension(:) :: lat_in
      !< are the latitude coordinates for the rectilinear input grid in radians
    real(FMS_HI_KIND_), intent(in), dimension(:) :: lon_out
      !< are the longitude coordinates for the rectilinear output grid in radians
    real(FMS_HI_KIND_), intent(in), dimension(:) :: lat_out
      !< are the latitude coordinates for the rectilinear output grid in radians
    integer, intent(in), optional :: verbose
      !< sets the verbosity level for output
    logical, intent(in), optional :: src_modulo
      !< is true if boundary condition along zonal boundary is cyclic
    integer :: i, j, ip1, im1, jp1, jm1
    integer :: nlon_in, nlat_in, nlon_out, nlat_out
    integer :: jcl, jcu, icl, icu, jj
    real(FMS_HI_KIND_) :: xz, yz, latitude_ratio
    integer :: iunit
    integer, parameter :: kindl = FMS_HI_KIND_ ! real size at compile time

    if(present(verbose)) verbose_bicubic = verbose
    ! Retain the historical periodic longitude behavior even when src_modulo is false.

    !--- get the grid size
    nlon_in  = size(lon_in) ; nlat_in  = size(lat_in)
    nlon_out = size(lon_out); nlat_out = size(lat_out)
    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out
    allocate ( Interp%HI_KIND_TYPE_%wti     (nlon_in, nlat_in, 3) )
    allocate ( Interp%HI_KIND_TYPE_%lon_in  (nlon_in) )
    allocate ( Interp%HI_KIND_TYPE_%lat_in  (nlat_in) )
    allocate ( Interp%HI_KIND_TYPE_%rat_x   (nlon_out, nlat_out) )
    allocate ( Interp%HI_KIND_TYPE_%rat_y   (nlon_out, nlat_out) )
    allocate ( Interp%i_lon   (nlon_out, nlat_out, 2) )
    allocate ( Interp%j_lat   (nlon_out, nlat_out, 2) )

    Interp%HI_KIND_TYPE_%lon_in = lon_in
    Interp%HI_KIND_TYPE_%lat_in = lat_in

    if ( verbose_bicubic > 0 ) then
       iunit = stdout()
       write (iunit,'(/,"Initialising bicubic interpolation, interface HORIZ_INTERP_BICUBIC_NEW_1D_")')
       write (iunit,'(/," Longitude of coarse grid points (radian): xc(i) i=1, ",i4)') Interp%nlon_src
       write (iunit,'(1x,10f10.4)') (Interp%HI_KIND_TYPE_%lon_in(jj),jj=1,Interp%nlon_src)
       write (iunit,'(/," Latitude of coarse grid points (radian):  yc(j) j=1, ",i4)') Interp%nlat_src
       write (iunit,'(1x,10f10.4)') (Interp%HI_KIND_TYPE_%lat_in(jj),jj=1,Interp%nlat_src)
       write (iunit,*)
       write (iunit,'(/," Longitude of fine grid points (radian): xf(i) i=1, ",i4)') Interp%nlat_dst
       write (iunit,'(1x,10f10.4)') (lon_out(jj),jj=1,Interp%nlon_dst)
       write (iunit,'(/," Latitude of fine grid points (radian):  yf(j) j=1, ",i4)') Interp%nlon_dst
       write (iunit,'(1x,10f10.4)') (lat_out(jj),jj=1,Interp%nlat_dst)
    endif


!---------------------------------------------------------------------------
!     Find the x-derivative. Use central differences and forward or
!     backward steps at the boundaries

    do j=1,nlat_in
      do i=1,nlon_in
        ip1=min(i+1,nlon_in)
        im1=max(i-1,1)
        Interp%HI_KIND_TYPE_%wti(i,j,1) = 1.0_kindl /(lon_in(ip1)-lon_in(im1))
      enddo
    enddo


!---------------------------------------------------------------------------

!     Find the y-derivative. Use central differences and forward or
!     backward steps at the boundaries
      do j=1,nlat_in
        jp1=min(j+1,nlat_in)
        jm1=max(j-1,1)
        do i=1,nlon_in
          Interp%HI_KIND_TYPE_%wti(i,j,2) = 1.0_kindl /(lat_in(jp1)-lat_in(jm1))
        enddo
      enddo

!---------------------------------------------------------------------------

!     Find the xy-derivative. Use central differences and forward or
!     backward steps at the boundaries
      do j=1,nlat_in
        jp1=min(j+1,nlat_in)
        jm1=max(j-1,1)
        do i=1,nlon_in
          ip1=min(i+1,nlon_in)
          im1=max(i-1,1)
          Interp%HI_KIND_TYPE_%wti(i,j,3) = 1.0_kindl /((lon_in(ip1)-lon_in(im1))*(lat_in(jp1)-lat_in(jm1)))
        enddo
      enddo
!---------------------------------------------------------------------------
!     Now for each point at the dest-grid find the boundary points of
!     the source grid
      do j=1, nlat_out
        yz  = lat_out(j)
        call bicubic_interval(lat_in, yz, .false., jcl, jcu, latitude_ratio)
        do i=1,nlon_out
          xz = lon_out(i)
          Interp%HI_KIND_TYPE_%rat_y(i,j) = latitude_ratio
          call bicubic_interval(lon_in, xz, .true., icl, icu, Interp%HI_KIND_TYPE_%rat_x(i,j))
          Interp%j_lat(i,j,1) = jcl
          Interp%j_lat(i,j,2) = jcu
          Interp%i_lon(i,j,1) = icl
          Interp%i_lon(i,j,2) = icu
        enddo
      enddo
    Interp% HI_KIND_TYPE_ % is_allocated = .true.
    Interp%interp_method = BICUBIC

  end subroutine HORIZ_INTERP_BICUBIC_NEW_1D_

  !> @parblock
  !! Interpolates 2D data using the weights and mapping indices stored in Interp
  !! @endparblock
  subroutine HORIZ_INTERP_BICUBIC_( Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value, &
                                 &  missing_permit)
    type (horiz_interp_type), intent(in) :: Interp
      !< contains the interpolation weights and mapping indices
    real(FMS_HI_KIND_), intent(in), dimension(:,:) :: data_in
      !< is the data on the source grid to interpolate
    real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
      !< will contain the interpolated data on the output grid
    integer, intent(in), optional :: verbose
      !< sets the verbosity level for output
    real(FMS_HI_KIND_), intent(in), dimension(:,:), optional :: mask_in
      !< is not used
    real(FMS_HI_KIND_), intent(out), dimension(:,:), optional :: mask_out
      !< is populated to be 1 for all points
    real(FMS_HI_KIND_), intent(in), optional :: missing_value
      !< is not used
    integer, intent(in), optional :: missing_permit
      !< is not used
    real(FMS_HI_KIND_) :: yz, ycu, ycl
    real(FMS_HI_KIND_) :: xz, xcu, xcl
    real(FMS_HI_KIND_) :: val, val1, val2
    real(FMS_HI_KIND_), dimension(4) :: y, y1, y2, y12
    integer :: icl, icu, jcl, jcu
    integer :: iclp1, icup1, jclp1, jcup1
    integer :: iclm1, icum1, jclm1, jcum1
    integer :: i,j
    integer, parameter :: kindl = FMS_HI_KIND_
      !< set kind size at compile time

    if ( present(verbose) ) verbose_bicubic = verbose

    do j=1, Interp%nlat_dst
      do i=1, Interp%nlon_dst
        yz  = Interp%HI_KIND_TYPE_%rat_y(i,j)
        xz  = Interp%HI_KIND_TYPE_%rat_x(i,j)
        jcl = Interp%j_lat(i,j,1)
        jcu = Interp%j_lat(i,j,2)
        icl = Interp%i_lon(i,j,1)
        icu = Interp%i_lon(i,j,2)
        ! Derivatives use storage neighbors and signed coordinate differences.
        ! A physical west/east pair has decreasing coordinate values only at the seam.
        xcl = Interp%HI_KIND_TYPE_%lon_in(icl)
        xcu = Interp%HI_KIND_TYPE_%lon_in(icu)
        if (xcu < xcl) xcu = xcu + real(tpi, FMS_HI_KIND_)
        iclp1 = min(icl+1,Interp%nlon_src)
        icum1 = max(icu-1,1)
        iclm1 = max(icl-1,1)
        icup1 = min(icu+1,Interp%nlon_src)
        jclp1 = min(jcl+1,Interp%nlat_src)
        jclm1 = max(jcl-1,1)
        jcup1 = min(jcu+1,Interp%nlat_src)
        jcum1 = max(jcu-1,1)
        ycl = Interp%HI_KIND_TYPE_%lat_in(jcl)
        ycu = Interp%HI_KIND_TYPE_%lat_in(jcu)
!        xcl = Interp%HI_KIND_TYPE_%lon_in(icl)
!        xcu = Interp%HI_KIND_TYPE_%lon_in(icu)
        y(1)  =  data_in(icl,jcl)
        y(2)  =  data_in(icu,jcl)
        y(3)  =  data_in(icu,jcu)
        y(4)  =  data_in(icl,jcu)
        y1(1) = ( data_in(iclp1,jcl) - data_in(iclm1,jcl) ) * Interp%HI_KIND_TYPE_%wti(icl,jcl,1)
        y1(2) = ( data_in(icup1,jcl) - data_in(icum1,jcl) ) * Interp%HI_KIND_TYPE_%wti(icu,jcl,1)
        y1(3) = ( data_in(icup1,jcu) - data_in(icum1,jcu) ) * Interp%HI_KIND_TYPE_%wti(icu,jcu,1)
        y1(4) = ( data_in(iclp1,jcu) - data_in(iclm1,jcu) ) * Interp%HI_KIND_TYPE_%wti(icl,jcu,1)
        y2(1) = ( data_in(icl,jclp1) - data_in(icl,jclm1) ) * Interp%HI_KIND_TYPE_%wti(icl,jcl,2)
        y2(2) = ( data_in(icu,jclp1) - data_in(icu,jclm1) ) * Interp%HI_KIND_TYPE_%wti(icu,jcl,2)
        y2(3) = ( data_in(icu,jcup1) - data_in(icu,jcum1) ) * Interp%HI_KIND_TYPE_%wti(icu,jcu,2)
        y2(4) = ( data_in(icl,jcup1) - data_in(icl,jcum1) ) * Interp%HI_KIND_TYPE_%wti(icl,jcu,2)
        y12(1)= ( data_in(iclp1,jclp1) + data_in(iclm1,jclm1) - data_in(iclm1,jclp1) &
                - data_in(iclp1,jclm1) ) * Interp%HI_KIND_TYPE_%wti(icl,jcl,3)
        y12(2)= ( data_in(icup1,jclp1) + data_in(icum1,jclm1) - data_in(icum1,jclp1) &
                - data_in(icup1,jclm1) ) * Interp%HI_KIND_TYPE_%wti(icu,jcl,3)
        y12(3)= ( data_in(icup1,jcup1) + data_in(icum1,jcum1) - data_in(icum1,jcup1) &
                - data_in(icup1,jcum1) ) * Interp%HI_KIND_TYPE_%wti(icu,jcu,3)
        y12(4)= ( data_in(iclp1,jcup1) + data_in(iclm1,jcum1) - data_in(iclm1,jcup1) &
                - data_in(iclp1,jcum1) ) * Interp%HI_KIND_TYPE_%wti(icl,jcu,3)

        call bcuint(y,y1,y2,y12,xcl,xcu,ycl,ycu,xz,yz,val,val1,val2)
        data_out   (i,j) = val
        if(present(mask_out)) mask_out(i,j) = 1.0_kindl
!!        dff_x(i,j) = val1 ! x-derivative of fc at the fine grid
!!        dff_y(i,j) = val2 ! y-derivative of fc at the fine grid
      enddo
    enddo
  return
  end subroutine HORIZ_INTERP_BICUBIC_

!---------------------------------------------------------------------------

   !> @parblock
   !! Internally used subroutine to interpolate data.  Adopted from
   !! Numerical Recipes for Fortran (1986-1992)
   !! @endparblock
   subroutine BCUINT_(y,y1,y2,y12,x1l,x1u,x2l,x2u,t,u,ansy,ansy1,ansy2)
      real(FMS_HI_KIND_) ansy,ansy1,ansy2,x1l,x1u,x2l,x2u,y(4),y1(4),y12(4),y2(4)
!     uses BCUCOF_
      integer i
      integer, parameter :: kindl = FMS_HI_KIND_
        !< compiled kind size
      real(FMS_HI_KIND_) t,u,c(4,4)
      call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
      ansy=0.0_kindl
      ansy2=0.0_kindl
      ansy1=0.0_kindl
      do i=4,1,-1
        ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
!        ansy2=t*ansy2+(3.*c(i,4)*u+2.*c(i,3))*u+c(i,2)
!        ansy1=u*ansy1+(3.*c(4,i)*t+2.*c(3,i))*t+c(2,i)
      enddo
!      ansy1=ansy1/(x1u-x1l) ! could be used for accuracy checks
!      ansy2=ansy2/(x2u-x2l) ! could be used for accuracy checks
      return
!  (c) copr. 1986-92 numerical recipes software -3#(-)f.
   end subroutine BCUINT_
!---------------------------------------------------------------------------

   !> @parblock
   !! Internally used subroutine called by bcuint_r4/r8 to compute bicubic coefficients.
   !! Adopted from Numerical Recipes for Fortran (1986-1992)
   !! @endparblock
   subroutine BCUCOF_(y,y1,y2,y12,d1,d2,c)
      real(FMS_HI_KIND_) d1,d2,c(4,4),y(4),y1(4),y12(4),y2(4)
      integer i,j,k,l
      real(FMS_HI_KIND_) d1d2,xx,cl(16),x(16)
      integer, parameter :: kindl = FMS_HI_KIND_
        !< compiled kind type
      !! n*0.0 represents n consecutive 0.0's
      real(FMS_HI_KIND_), save, dimension(16,16) :: wt
        !< weights use
      data wt/1.0_kindl, 0.0_kindl, -3.0_kindl, 2.0_kindl, 4*0.0_kindl, -3.0_kindl, 0.0_kindl, 9.0_kindl, -6.0_kindl, &
              2.0_kindl, 0.0_kindl, -6.0_kindl, 4.0_kindl, 8*0.0_kindl, 3.0_kindl, 0.0_kindl, -9.0_kindl, 6.0_kindl,  &
              -2.0_kindl, 0.0_kindl, 6.0_kindl, -4.0_kindl, 10*0.0_kindl, 9.0_kindl, -6.0_kindl, 2*0.0_kindl,         &
              -6.0_kindl, 4.0_kindl, 2*0.0_kindl, 3.0_kindl, -2.0_kindl, 6*0.0_kindl, -9.0_kindl, 6.0_kindl,          &
              2*0.0_kindl, 6.0_kindl, -4.0_kindl, 4*0.0_kindl, 1.0_kindl, 0.0_kindl, -3.0_kindl, 2.0_kindl,-2.0_kindl,&
              0.0_kindl, 6.0_kindl, -4.0_kindl, 1.0_kindl, 0.0_kindl, -3.0_kindl, 2.0_kindl, 8*0.0_kindl, -1.0_kindl, &
              0.0_kindl, 3.0_kindl, -2.0_kindl, 1.0_kindl, 0.0_kindl, -3.0_kindl, 2.0_kindl, 10*0.0_kindl, -3.0_kindl,&
              2.0_kindl, 2*0.0_kindl, 3.0_kindl, -2.0_kindl, 6*0.0_kindl, 3.0_kindl, -2.0_kindl, 2*0.0_kindl,         &
              -6.0_kindl, 4.0_kindl, 2*0.0_kindl, 3.0_kindl, -2.0_kindl, 0.0_kindl, 1.0_kindl, -2.0_kindl, 1.0_kindl, &
              5*0.0_kindl, -3.0_kindl, 6.0_kindl, -3.0_kindl, 0.0_kindl, 2.0_kindl, -4.0_kindl, 2.0_kindl,9*0.0_kindl,&
              3.0_kindl, -6.0_kindl, 3.0_kindl, 0.0_kindl, -2.0_kindl, 4.0_kindl, -2.0_kindl, 10*0.0_kindl,-3.0_kindl,&
              3.0_kindl, 2*0.0_kindl, 2.0_kindl, -2.0_kindl, 2*0.0_kindl, -1.0_kindl, 1.0_kindl,6*0.0_kindl,3.0_kindl,&
              -3.0_kindl, 2*0.0_kindl, -2.0_kindl, 2.0_kindl, 5*0.0_kindl, 1.0_kindl, -2.0_kindl, 1.0_kindl,0.0_kindl,&
              -2.0_kindl, 4.0_kindl, -2.0_kindl, 0.0_kindl, 1.0_kindl, -2.0_kindl, 1.0_kindl, 9*0.0_kindl, -1.0_kindl,&
              2.0_kindl, -1.0_kindl, 0.0_kindl, 1.0_kindl, -2.0_kindl, 1.0_kindl, 10*0.0_kindl, 1.0_kindl, -1.0_kindl,&
              2*0.0_kindl, -1.0_kindl, 1.0_kindl, 6*0.0_kindl, -1.0_kindl, 1.0_kindl, 2*0.0_kindl, 2.0_kindl,         &
              -2.0_kindl, 2*0.0_kindl, -1.0_kindl, 1.0_kindl/



      d1d2=d1*d2
      do i=1,4
        x(i)=y(i)
        x(i+4)=y1(i)*d1
        x(i+8)=y2(i)*d2
        x(i+12)=y12(i)*d1d2
      enddo
      do i=1,16
        xx=0.0_kindl
        do k=1,16
          xx=xx+wt(i,k)*x(k)
        enddo
        cl(i)=xx
      enddo
      l=0
      do i=1,4
        do j=1,4
          l=l+1
          c(i,j)=cl(l)
        enddo
      enddo
      return
!  (c) copr. 1986-92 numerical recipes software -3#(-)f.
   end subroutine BCUCOF_

!-----------------------------------------------------------------------

   !> Select physical lower/upper neighbors, independently of storage direction.
   !! Latitude is clamped; periodic longitude is reduced to [minimum,minimum+2*pi).
   subroutine BICUBIC_INTERVAL_(xc, xf, cyclic, lower, upper, ratio)
     real(FMS_HI_KIND_), intent(in) :: xc(:), xf
     logical, intent(in) :: cyclic
     integer, intent(out) :: lower, upper
     real(FMS_HI_KIND_), intent(out) :: ratio
     integer :: first, last, step, ii, next
     real(FMS_HI_KIND_) :: x
     first = 1
     last = size(xc)
     step = 1
     if (xc(first) > xc(last)) then
       first = size(xc)
       last = 1
       step = -1
     endif
     x = xf
     if (cyclic) then
       x = xc(first) + modulo(x-xc(first), real(tpi,FMS_HI_KIND_))
       if (x > xc(last)) then
         lower = last
         upper = first
         ratio = (x-xc(last))/(xc(first)+real(tpi,FMS_HI_KIND_)-xc(last))
         return
       endif
     endif
     lower = first
     upper = first
     ratio = real(0,FMS_HI_KIND_)
     if (x <= xc(first)) return
     lower = last
     upper = last
     if (x >= xc(last)) return
     do ii = first, last-step, step
       next = ii+step
       if (x < xc(next)) then
         lower = ii
         upper = next
         ratio = (x-xc(ii))/(xc(next)-xc(ii))
         return
       endif
     enddo
   end subroutine BICUBIC_INTERVAL_

!-----------------------------------------------------------------------

    !> @parblock
    !! Unused subroutine.
    !! @endparblock
    subroutine FILL_XY_(fi, ics, ice, jcs, jce, mask, maxpass)
      integer, intent(in) :: ics,ice,jcs,jce
      real(FMS_HI_KIND_), intent(inout) :: fi(ics:ice,jcs:jce)
      real(FMS_HI_KIND_), intent(in), optional :: mask(ics:ice,jcs:jce)
      integer, intent(in) :: maxpass
      real(FMS_HI_KIND_) :: work_old(ics:ice,jcs:jce)
      real(FMS_HI_KIND_) :: work_new(ics:ice,jcs:jce)
      logical :: ready
      integer, parameter :: kindl = FMS_HI_KIND_
      real(FMS_HI_KIND_), parameter :: blank = real(-1.e30, FMS_HI_KIND_)
      real(FMS_HI_KIND_) :: tavr
      integer :: ipass
      integer :: inl, inr, jnl, jnu, i, j, is, js, iavr


      ready = .false.

      work_new(:,:) = fi(:,:)
      work_old(:,:) = work_new(:,:)
      ipass = 0
      if ( present(mask) ) then
         do while (.not.ready)
           ipass = ipass+1
           ready = .true.
           do j=jcs, jce
             do i=ics, ice
               if (work_old(i,j).le.blank) then
                 tavr=0.0_kindl
                 iavr=0
                 inl = max(i-1,ics)
                 inr = min(i+1,ice)
                 jnl = max(j-1,jcs)
                 jnu = min(j+1,jce)
                 do js=jnl,jnu
                   do is=inl,inr
                     if (work_old(is,js) .ne. blank .and. mask(is,js).ne.0.0_kindl) then
                       tavr = tavr + work_old(is,js)
                       iavr = iavr+1
                     endif
                   enddo
                 enddo
                 if (iavr.gt.0) then
                   if (iavr.eq.1) then
! spreading is not allowed if the only valid neighbor is a corner point
! otherwise an ill posed cellular automaton is established leading to
! a spreading of constant values in diagonal direction
! if all corner points are blanked the valid neighbor must be a direct one
! and spreading is allowed
                     if (work_old(inl,jnu).eq.blank.and.&
                         work_old(inr,jnu).eq.blank.and.&
                         work_old(inr,jnl).eq.blank.and.&
                         work_old(inl,jnl).eq.blank) then
                           work_new(i,j)=tavr/real(iavr,FMS_HI_KIND_)
                           ready = .false.
                     endif
                  else
                    work_new(i,j)=tavr/real(iavr,FMS_HI_KIND_)
                    ready = .false.
                  endif
                endif
              endif
            enddo ! j
          enddo   ! i
! save changes made during this pass to work_old
          work_old(:,:)=work_new(:,:)
          if(ipass.eq.maxpass) ready=.true.
        enddo !while (.not.ready)
        fi(:,:) = work_new(:,:)
      else
         do while (.not.ready)
           ipass = ipass+1
           ready = .true.
           do j=jcs, jce
             do i=ics, ice
               if (work_old(i,j).le.blank) then
                 tavr=0.0_kindl
                 iavr=0
                 inl = max(i-1,ics)
                 inr = min(i+1,ice)
                 jnl = max(j-1,jcs)
                 jnu = min(j+1,jce)
                 do is=inl,inr
                   do js=jnl,jnu
                     if (work_old(is,js).gt.blank) then
                       tavr = tavr + work_old(is,js)
                       iavr = iavr+1
                     endif
                   enddo
                 enddo
                 if (iavr.gt.0) then
                   if (iavr.eq.1) then
! spreading is not allowed if the only valid neighbor is a corner point
! otherwise an ill posed cellular automaton is established leading to
! a spreading of constant values in diagonal direction
! if all corner points are blanked the valid neighbor must be a direct one
! and spreading is allowed
                     if (work_old(inl,jnu).le.blank.and. &
                         work_old(inr,jnu).le.blank.and. &
                         work_old(inr,jnl).le.blank.and. &
                         work_old(inl,jnl).le.blank) then
                           work_new(i,j)=tavr/real(iavr,FMS_HI_KIND_)
                           ready = .false.
                     endif
                  else
                    work_new(i,j)=tavr/real(iavr,FMS_HI_KIND_)
                    ready = .false.
                  endif
                endif
              endif
            enddo ! j
          enddo   ! i
! save changes made during this pass to work_old
          work_old(:,:)=work_new(:,:)
          if(ipass.eq.maxpass) ready=.true.
        enddo !while (.not.ready)
        fi(:,:) = work_new(:,:)
      endif
      return
    end subroutine FILL_XY_
