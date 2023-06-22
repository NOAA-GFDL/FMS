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
!> @addtogroup horiz_interp_bicubic_mod
!> @{

  !> @brief Creates a new @ref horiz_interp_type
  !!
  !> Allocates space and initializes a derived-type variable
  !! that contains pre-computed interpolation indices and weights.
  subroutine HORIZ_INTERP_BICUBIC_NEW_1D_S_ ( Interp, lon_in, lat_in, lon_out, lat_out, &
       verbose, src_modulo )

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout) :: Interp !< A derived-type variable containing indices
                                       !! and weights used for subsequent interpolations. To
                                       !! reinitialize this variable for a different grid-to-grid
                                       !! interpolation you must first use the
                                       !! @ref HORIZ_INTERP_BICUBIC_NEW__del interface.
    real(FMS_HI_KIND_), intent(in),  dimension(:)        :: lon_in !< Longitude (radians) for source data grid
    real(FMS_HI_KIND_), intent(in),  dimension(:)        :: lat_in !< Latitude (radians) for source data grid
    real(FMS_HI_KIND_), intent(in),  dimension(:,:)      :: lon_out !< Longitude (radians) for output data grid
    real(FMS_HI_KIND_), intent(in),  dimension(:,:)      :: lat_out !< Latitude (radians) for output data grid
    integer, intent(in),          optional :: verbose !< flag for print output amount
    logical, intent(in),          optional :: src_modulo !< indicates if the boundary condition along
                                       !! zonal boundary is cyclic or not. Zonal boundary condition
                                       !!is cyclic when true
    integer                                :: i, j, ip1, im1, jp1, jm1
    logical                                :: src_is_modulo
    integer                                :: nlon_in, nlat_in, nlon_out, nlat_out
    integer                                :: jcl, jcu, icl, icu, jj
    real(FMS_HI_KIND_)                     :: xz, yz
    integer                                :: unit
    integer, parameter   :: kindl = FMS_HI_KIND_ !< real size at compile time

    if(present(verbose)) verbose_bicubic = verbose
    src_is_modulo = .false.
    if (present(src_modulo)) src_is_modulo = src_modulo

    if(size(lon_out,1) /= size(lat_out,1) .or. size(lon_out,2) /= size(lat_out,2) ) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: when using bilinear ' // &
         'interplation, the output grids should be geographical grids')

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
       unit = stdout()
       write (unit,'(/,"Initialising bicubic interpolation, interface horiz_interp_bicubic_new_1d_s")')
       write (unit,'(/," Longitude of coarse grid points (radian): xc(i) i=1, ",i4)') Interp%nlon_src
       write (unit,'(1x,10f10.4)') (Interp%HI_KIND_TYPE_%lon_in(jj),jj=1,Interp%nlon_src)
       write (unit,'(/," Latitude of coarse grid points (radian):  yc(j) j=1, ",i4)') Interp%nlat_src
       write (unit,'(1x,10f10.4)') (Interp%HI_KIND_TYPE_%lat_in(jj),jj=1,Interp%nlat_src)
       do i=1, Interp%nlat_dst
         write (unit,*)
         write (unit,'(/," Longitude of fine grid points (radian): xf(i) i=1, ",i4)') Interp%nlat_dst
         write (unit,'(1x,10f10.4)') (lon_out(jj,i),jj=1,Interp%nlon_dst)
       enddo
       do i=1, Interp%nlon_dst
         write (unit,*)
         write (unit,'(/," Latitude of fine grid points (radian):  yf(j) j=1, ",i4)') Interp%nlon_dst
         write (unit,'(1x,10f10.4)') (lat_out(i,jj),jj=1,Interp%nlat_dst)
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

          jcl = 0
          jcu = 0
          if( yz .le. Interp%HI_KIND_TYPE_%lat_in(1) ) then
             jcl = 1
             jcu = 1
          else if( yz .ge. Interp%HI_KIND_TYPE_%lat_in(nlat_in) ) then
             jcl = nlat_in
             jcu = nlat_in
          else
             jcl = indl(Interp%HI_KIND_TYPE_%lat_in, yz)
             jcu = indu(Interp%HI_KIND_TYPE_%lat_in, yz)
          endif

          icl = 0
          icu = 0
          !--- cyclic condition, do we need to use do while
          if( xz .gt. Interp%HI_KIND_TYPE_%lon_in(nlon_in) ) xz = xz - real(tpi,FMS_HI_KIND_)
          if( xz .le. Interp%HI_KIND_TYPE_%lon_in(1) ) xz = xz + real(tpi,FMS_HI_KIND_)
          if( xz .ge. Interp%HI_KIND_TYPE_%lon_in(nlon_in) ) then
            icl = nlon_in
            icu = 1
            Interp%HI_KIND_TYPE_%rat_x(i,j) = (xz - Interp%HI_KIND_TYPE_%lon_in(icl))/(Interp%HI_KIND_TYPE_%lon_in(icu)&
                                             & - Interp%HI_KIND_TYPE_%lon_in(icl) + real(tpi,FMS_HI_KIND_))
          else
            icl = indl(Interp%HI_KIND_TYPE_%lon_in, xz)
            icu = indu(Interp%HI_KIND_TYPE_%lon_in, xz)
            Interp%HI_KIND_TYPE_%rat_x(i,j) = (xz - Interp%HI_KIND_TYPE_%lon_in(icl))/(Interp%HI_KIND_TYPE_%lon_in(icu)&
                                             & - Interp%HI_KIND_TYPE_%lon_in(icl))
          endif
          Interp%j_lat(i,j,1) = jcl
          Interp%j_lat(i,j,2) = jcu
          Interp%i_lon(i,j,1) = icl
          Interp%i_lon(i,j,2) = icu
          if(jcl == jcu) then
             Interp%HI_KIND_TYPE_%rat_y(i,j) = 0.0_kindl
          else
             Interp%HI_KIND_TYPE_%rat_y(i,j) = (yz-Interp%HI_KIND_TYPE_%lat_in(jcl))/(Interp%HI_KIND_TYPE_%lat_in(jcu)&
                                              & - Interp%HI_KIND_TYPE_%lat_in(jcl))
          endif
!          if(yz.gt.Interp%HI_KIND_TYPE_%lat_in(jcu)) call mpp_error(FATAL, ' HORIZ_INTERP_BICUBIC_NEW_1D_S_:
!          yf < ycl, no valid boundary point')
!          if(yz.lt.Interp%HI_KIND_TYPE_%lat_in(jcl)) call mpp_error(FATAL, ' HORIZ_INTERP_BICUBIC_NEW_1D_S_:
!          yf > ycu, no valid boundary point')
!          if(xz.gt.Interp%HI_KIND_TYPE_%lon_in(icu)) call mpp_error(FATAL, ' HORIZ_INTERP_BICUBIC_NEW_1D_S_:
!          xf < xcl, no valid boundary point')
!          if(xz.lt.Interp%HI_KIND_TYPE_%lon_in(icl)) call mpp_error(FATAL, ' HORIZ_INTERP_BICUBIC_NEW_1D_S_:
!          xf > xcu, no valid boundary point')
        enddo
      enddo
  end subroutine HORIZ_INTERP_BICUBIC_NEW_1D_S_

  !> @brief Creates a new @ref horiz_interp_type
  !!
  !> Allocates space and initializes a derived-type variable
  !! that contains pre-computed interpolation indices and weights.
  subroutine HORIZ_INTERP_BICUBIC_NEW_1D_ ( Interp, lon_in, lat_in, lon_out, lat_out, &
       verbose, src_modulo )

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout) :: Interp
    real(FMS_HI_KIND_), intent(in),  dimension(:)        :: lon_in , lat_in
    real(FMS_HI_KIND_), intent(in),  dimension(:)        :: lon_out, lat_out
    integer, intent(in),          optional :: verbose
    logical, intent(in),          optional :: src_modulo
    integer                                :: i, j, ip1, im1, jp1, jm1
    logical                                :: src_is_modulo
    integer                                :: nlon_in, nlat_in, nlon_out, nlat_out
    integer                                :: jcl, jcu, icl, icu, jj
    real(FMS_HI_KIND_)                                   :: xz, yz
    integer                                :: unit
    integer, parameter                     :: kindl = FMS_HI_KIND_ !< real size at compile time

    if(present(verbose)) verbose_bicubic = verbose
    src_is_modulo = .false.
    if (present(src_modulo)) src_is_modulo = src_modulo

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
       unit = stdout()
       write (unit,'(/,"Initialising bicubic interpolation, interface HORIZ_INTERP_BICUBIC_NEW_1D_")')
       write (unit,'(/," Longitude of coarse grid points (radian): xc(i) i=1, ",i4)') Interp%nlon_src
       write (unit,'(1x,10f10.4)') (Interp%HI_KIND_TYPE_%lon_in(jj),jj=1,Interp%nlon_src)
       write (unit,'(/," Latitude of coarse grid points (radian):  yc(j) j=1, ",i4)') Interp%nlat_src
       write (unit,'(1x,10f10.4)') (Interp%HI_KIND_TYPE_%lat_in(jj),jj=1,Interp%nlat_src)
       write (unit,*)
       write (unit,'(/," Longitude of fine grid points (radian): xf(i) i=1, ",i4)') Interp%nlat_dst
       write (unit,'(1x,10f10.4)') (lon_out(jj),jj=1,Interp%nlon_dst)
       write (unit,'(/," Latitude of fine grid points (radian):  yf(j) j=1, ",i4)') Interp%nlon_dst
       write (unit,'(1x,10f10.4)') (lat_out(jj),jj=1,Interp%nlat_dst)
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
        jcl = 0
        jcu = 0
        if( yz .le. lat_in(1) ) then
           jcl = 1
           jcu = 1
        else if( yz .ge. lat_in(nlat_in) ) then
           jcl = nlat_in
           jcu = nlat_in
        else
           jcl = indl(lat_in, yz)
           jcu = indu(lat_in, yz)
        endif
        do i=1,nlon_out
          xz = lon_out(i)
          icl = 0
          icu = 0
         !--- cyclic condition, do we need to use do while
          if( xz .gt. lon_in(nlon_in) ) xz = xz - real(tpi,FMS_HI_KIND_)
          if( xz .le. lon_in(1) ) xz = xz + real(tpi, FMS_HI_KIND_)
          if( xz .ge. lon_in(nlon_in) ) then
            icl = nlon_in
            icu = 1
            Interp%HI_KIND_TYPE_%rat_x(i,j) = (xz - Interp%HI_KIND_TYPE_%lon_in(icl))/(Interp%HI_KIND_TYPE_%lon_in(icu)&
                                            & - Interp%HI_KIND_TYPE_%lon_in(icl) + real(tpi,FMS_HI_KIND_))
          else
            icl = indl(lon_in, xz)
            icu = indu(lon_in, xz)
            Interp%HI_KIND_TYPE_%rat_x(i,j) = (xz - Interp%HI_KIND_TYPE_%lon_in(icl))/(Interp%HI_KIND_TYPE_%lon_in(icu)&
                                            & - Interp%HI_KIND_TYPE_%lon_in(icl))
          endif
          icl = indl(lon_in, xz)
          icu = indu(lon_in, xz)
          Interp%j_lat(i,j,1) = jcl
          Interp%j_lat(i,j,2) = jcu
          Interp%i_lon(i,j,1) = icl
          Interp%i_lon(i,j,2) = icu
          if(jcl == jcu) then
             Interp%HI_KIND_TYPE_%rat_y(i,j) = 0.0_kindl
          else
             Interp%HI_KIND_TYPE_%rat_y(i,j) = (yz- Interp%HI_KIND_TYPE_%lat_in(jcl))/(Interp%HI_KIND_TYPE_%lat_in(jcu)&
                                             & - Interp%HI_KIND_TYPE_%lat_in(jcl))
          endif
!          if(yz.gt.lat_in(jcu)) call mpp_error(FATAL, ' HORIZ_INTERP_BICUBIC_NEW_1D_: yf <
!          ycl, no valid boundary point')
!          if(yz.lt.lat_in(jcl)) call mpp_error(FATAL, ' HORIZ_INTERP_BICUBIC_NEW_1D_: yf >
!          ycu, no valid boundary point')
!          if(xz.gt.lon_in(icu)) call mpp_error(FATAL, ' HORIZ_INTERP_BICUBIC_NEW_1D_: xf <
!          xcl, no valid boundary point')
!          if(xz.lt.lon_in(icl)) call mpp_error(FATAL, ' HORIZ_INTERP_BICUBIC_NEW_1D_: xf >
!          xcu, no valid boundary point')
        enddo
      enddo

  end subroutine HORIZ_INTERP_BICUBIC_NEW_1D_

  !> @brief Perform bicubic horizontal interpolation
  subroutine HORIZ_INTERP_BICUBIC_NEW_( Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value, &
                                 &  missing_permit)
    type (horiz_interp_type), intent(in)        :: Interp
    real(FMS_HI_KIND_), intent(in),  dimension(:,:)           :: data_in
    real(FMS_HI_KIND_), intent(out), dimension(:,:)           :: data_out
    integer, intent(in),               optional :: verbose
    real(FMS_HI_KIND_), intent(in),  dimension(:,:), optional :: mask_in
    real(FMS_HI_KIND_), intent(out), dimension(:,:), optional :: mask_out
    real(FMS_HI_KIND_), intent(in),                  optional :: missing_value
    integer, intent(in),               optional :: missing_permit
    real(FMS_HI_KIND_) :: yz, ycu, ycl
    real(FMS_HI_KIND_) :: xz, xcu, xcl
    real(FMS_HI_KIND_) :: val, val1, val2
    real(FMS_HI_KIND_), dimension(4) :: y, y1, y2, y12
    integer :: icl, icu, jcl, jcu
    integer :: iclp1, icup1, jclp1, jcup1
    integer :: iclm1, icum1, jclm1, jcum1
    integer :: i,j
    integer, parameter :: kindl = FMS_HI_KIND_ !< set kind size at compile time

    if ( present(verbose) ) verbose_bicubic = verbose

    do j=1, Interp%nlat_dst
      do i=1, Interp%nlon_dst
        yz  = Interp%HI_KIND_TYPE_%rat_y(i,j)
        xz  = Interp%HI_KIND_TYPE_%rat_x(i,j)
        jcl = Interp%j_lat(i,j,1)
        jcu = Interp%j_lat(i,j,2)
        icl = Interp%i_lon(i,j,1)
        icu = Interp%i_lon(i,j,2)
        if( icl > icu ) then
          iclp1 = icu
          icum1 = icl
          xcl = Interp%HI_KIND_TYPE_%lon_in(icl)
          xcu = Interp%HI_KIND_TYPE_%lon_in(icu)+real(tpi, FMS_HI_KIND_)
        else
          iclp1 = min(icl+1,Interp%nlon_src)
          icum1 = max(icu-1,1)
          xcl = Interp%HI_KIND_TYPE_%lon_in(icl)
          xcu = Interp%HI_KIND_TYPE_%lon_in(icu)
        endif
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
!!        dff_x(i,j) = val1
!!        dff_y(i,j) = val2
      enddo
    enddo
  return
  end subroutine HORIZ_INTERP_BICUBIC_NEW_

!---------------------------------------------------------------------------

   subroutine BCUINT_(y,y1,y2,y12,x1l,x1u,x2l,x2u,t,u,ansy,ansy1,ansy2)
      real(FMS_HI_KIND_) ansy,ansy1,ansy2,x1l,x1u,x2l,x2u,y(4),y1(4),y12(4),y2(4)
!     uses BCUCOF_
      integer i
      integer, parameter :: kindl = FMS_HI_KIND_ !< compiled kind size
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

   subroutine BCUCOF_(y,y1,y2,y12,d1,d2,c)
      real(FMS_HI_KIND_) d1,d2,c(4,4),y(4),y1(4),y12(4),y2(4)
      integer i,j,k,l
      real(FMS_HI_KIND_) d1d2,xx,cl(16),x(16)
      integer, parameter :: kindl = FMS_HI_KIND_!< compiled kind type
      !! n*0.0 represents n consecutive 0.0's
      real(FMS_HI_KIND_), save, dimension(16,16) :: wt !< weights use
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

!! TODO These routines are redundant, we can find the lower neighbor and add 1
!> find the lower neighbour of xf in field xc, return is the index
    function INDL_(xc, xf)
    real(FMS_HI_KIND_), intent(in) :: xc(1:)
    real(FMS_HI_KIND_), intent(in) :: xf
    integer             :: INDL_
    integer             :: ii
       INDL_ = 1
       do ii=1, size(xc)
         if(xc(ii).gt.xf) return
         INDL_ = ii
       enddo
       call mpp_error(FATAL,'Error in INDL_')
    return
    end function INDL_

!-----------------------------------------------------------------------

!> find the upper neighbour of xf in field xc, return is the index
    function INDU_(xc, xf)
    real(FMS_HI_KIND_), intent(in) :: xc(1:)
    real(FMS_HI_KIND_), intent(in) :: xf
    integer             :: INDU_
    integer             :: ii
       do ii=1, size(xc)
         INDU_ = ii
         if(xc(ii).gt.xf) return
       enddo
       call mpp_error(FATAL,'Error in INDU_')
    return
    end function INDU_

!-----------------------------------------------------------------------

    subroutine FILL_XY_(fi, ics, ice, jcs, jce, mask, maxpass)
      integer, intent(in)        :: ics,ice,jcs,jce
      real(FMS_HI_KIND_), intent(inout)        :: fi(ics:ice,jcs:jce)
      real(FMS_HI_KIND_), intent(in), optional :: mask(ics:ice,jcs:jce)
      integer, intent(in)        :: maxpass
      real(FMS_HI_KIND_)                       :: work_old(ics:ice,jcs:jce)
      real(FMS_HI_KIND_)                       :: work_new(ics:ice,jcs:jce)
      logical :: ready
      integer, parameter :: kindl = FMS_HI_KIND_
      real(FMS_HI_KIND_), parameter    :: blank = real(-1.e30, FMS_HI_KIND_)
      real(FMS_HI_KIND_)    :: tavr
      integer :: ipass
      integer :: inl, inr, jnl, jnu, i, j, is, js,  iavr


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
!> @}
