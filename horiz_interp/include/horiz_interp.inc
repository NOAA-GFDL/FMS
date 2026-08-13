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
  !! Wrapper to compute weights for conservative, bilinear, bicubic, or spherical
  !! interpolation for rectangular input and output grids that are specified as
  !! 1D arrays.  The 1D array of gridpoints are used to construct 2D grids before
  !! calling either horiz_interp_conserve_new, horiz_interp_bilinear_new, horiz_interp_bicubic_new,
  !! or horiz_interp_spherical_new.
  !! @endparblock
  subroutine HORIZ_INTERP_NEW_1D_ (Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                                  interp_method, num_nbrs, max_dist, src_modulo,     &
                                  grid_at_center, mask_in, mask_out)

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout)  :: Interp
      !< is a derived type that will be populated with interpolation weights and indices
    real(FMS_HI_KIND_), intent(in),  dimension(:) :: lon_in
      !< are the longitude coordinates [radians] for the rectangular input grid
    real(FMS_HI_KIND_), intent(in),  dimension(:) :: lat_in
      !< are the latitude coordinates [radians] for the rectangular input grid
    real(FMS_HI_KIND_), intent(in),  dimension(:) :: lon_out
      !< are the longitude coordinates [radians] for the rectangular output grid
    real(FMS_HI_KIND_), intent(in), dimension(:) :: lat_out
      !< are the latitude coordinates [radians] for the rectangular output grid
    integer, intent(in), optional :: verbose
      !< sets the verbosity level
    character(len=*), intent(in), optional :: interp_method
      !< is the interpolation method; can be "conservative" (default), "bilinear", "bicubic", or "spherical"
    integer, intent(in), optional :: num_nbrs
      !< is the number of nearest neighbor for spherical interp; defaults to 4
    real(FMS_HI_KIND_),  intent(in), optional :: max_dist
      !< is the maximum region of influence around output grid cell for spherical interp; defaults to 0.1 radians
    logical, intent(in), optional :: src_modulo
      !< is true if boundary condition along zonal boundary is cyclic
    logical, intent(in), optional :: grid_at_center
      !< is true if the grid points are the grid center points for bilinear interp; if false, computes center points
    real(FMS_HI_KIND_), intent(in), dimension(:,:), optional :: mask_in
      !< is an unused variable
    real(FMS_HI_KIND_), intent(inout),dimension(:,:),  optional :: mask_out
      !< is an unused variable

    real(FMS_HI_KIND_), dimension(:,:), allocatable :: lon_src, lat_src, lon_dst, lat_dst
    real(FMS_HI_KIND_), dimension(:),   allocatable :: lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d
    integer  :: i, j, nlon_in, nlat_in, nlon_out, nlat_out
    logical  :: center
    character(len=40)  :: method
    integer, parameter  :: kindl = FMS_HI_KIND_ !> real kind size currently compiling
    !-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method

    select case (trim(method))
    case ("conservative")
       Interp%interp_method = CONSERVE
       call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose)
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

          call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_dst, lat_dst)
       else
          nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i = 1, nlon_in
             lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_in
             lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
          enddo
          do i = 1, nlon_out
             lon_dst(i,:) = (lon_out(i) + lon_out(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_out
             lat_dst(:,j) = (lat_out(j) + lat_out(j+1)) * 0.5_kindl
          enddo
          call horiz_interp_bilinear_new ( Interp, lon_src_1d, lat_src_1d, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst, lat_dst)
       endif
    case ("bicubic")
       Interp%interp_method = BICUBIC
       center = .false.
       if(present(grid_at_center) ) center = grid_at_center
       !No need to expand to 2d, horiz_interp_bicubic_new does 1d-1d
       if(center) then
          call horiz_interp_bicubic_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
            verbose, src_modulo)
       else
          nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst_1d(nlon_out), lat_dst_1d(nlat_out))
          do i = 1, nlon_in
             lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_in
             lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
          enddo
          do i = 1, nlon_out
             lon_dst_1d(i) = (lon_out(i) + lon_out(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_out
             lat_dst_1d(j) = (lat_out(j) + lat_out(j+1)) * 0.5_kindl
          enddo
          call horiz_interp_bicubic_new ( Interp, lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d, &
               verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d)
       endif
    case ("spherical")
       Interp%interp_method = SPHERICAL
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
       call horiz_interp_spherical_new ( Interp, lon_src, lat_src, lon_dst, lat_dst, &
            num_nbrs, max_dist, src_modulo)
       deallocate(lon_src, lat_src, lon_dst, lat_dst)
    case default
       call mpp_error(FATAL,'horiz_interp_mod: interp_method should be conservative, bilinear, bicubic, spherical')
    end select

    !-----------------------------------------------------------------------
    Interp% HI_KIND_TYPE_ % is_allocated = .true.
    Interp%I_am_initialized = .true.

  end subroutine HORIZ_INTERP_NEW_1D_

!#######################################################################
  !> @parblock
  !! Wrapper to compute weights for conservative, bilinear, bicubic, or spherical
  !! interpolation for rectangular input grids specified as 1D arrays and 2D output grids
  !! of any kind.  The 1D array of input gridpoints will be used to construct
  !! 2D grids before calling either horiz_interp_conserve_new, horiz_interp_bilinear_new,
  !! horiz_interp_bicubic_new, or horiz_interp_spherical_new.
  !! @endparblock
 subroutine HORIZ_INTERP_NEW_1D_SRC_ (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                     verbose, interp_method, num_nbrs, max_dist, &
                                     src_modulo, grid_at_center, mask_in, mask_out, is_latlon_out )

   type(horiz_interp_type), intent(inout) :: Interp
     !< is a derived type that will be populated with interpolation weights and indices
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lon_in
     !< are the longitude coordinates [radians] for the rectangular input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lat_in
     !< are the latitude coordinates [radians] for the rectangular input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_out
     !< are the 2D lontitude coordinates [radians] for the output grid
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lat_out
     !< are the 2D latitude coordinates [radians] for the output grid
   integer, intent(in), optional :: verbose
     !< sets the verbosity level
   character(len=*), intent(in), optional :: interp_method
     !< is the interpolation method; can be "conservative" (default), "bilinear", "bicubic", or "spherical"
   integer, intent(in), optional :: num_nbrs
     !< is the number of nearest neighbor for spherical interp; defaults to 4
   real(FMS_HI_KIND_), intent(in), optional :: max_dist
     !< is the maximum region of influence around output grid cell for spherical interp; defaults to 0.1 radians
   logical, intent(in), optional :: src_modulo
     !< is true if boundary condition along zonal boundary is cyclic
   logical, intent(in), optional :: grid_at_center
     !< is true if the grid points are the grid center points for bilinear interp; if false, computes center points
   real(FMS_HI_KIND_), intent(in), dimension(:,:), optional :: mask_in
     !< used in conservative interp to mask input grid cells (skip interp for masked cells)
   real(FMS_HI_KIND_), intent(out),dimension(:,:), optional :: mask_out
     !< will be populated with output grid cell areas computed from exchange grid for conserve interp
   logical, intent(in), optional :: is_latlon_out
     !< is true if output grid is a lat-lon grid

   real(FMS_HI_KIND_), dimension(:,:), allocatable :: lon_src, lat_src
   real(FMS_HI_KIND_), dimension(:),   allocatable :: lon_src_1d, lat_src_1d
   integer  :: i, j, nlon_in, nlat_in
   character(len=40)  :: method
   logical  :: center
   logical  :: dst_is_latlon
   integer, parameter  :: kindl = FMS_HI_KIND_ !< real kind size currently compiling
   !-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'conservative'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      !--- check to see if the source grid is regular lat-lon grid or not.
      if(PRESENT(is_latlon_out)) then
         dst_is_latlon = is_latlon_out
      else
         dst_is_latlon = is_lat_lon(lon_out, lat_out)
      end if
      if(dst_is_latlon ) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -.0001_kindl) .or. ANY(mask_in > 1.0001_kindl)) &
               call mpp_error(FATAL, &
                  'horiz_interp_conserve_new_1d_src(horiz_interp_conserve_mod): input mask not between 0,1')
            allocate(Interp%HI_KIND_TYPE_%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%HI_KIND_TYPE_%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out(:,1), lat_out(1,:), &
              verbose=verbose )
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if
   case ("bilinear")
      Interp%interp_method = BILINEAR
      center = .false.
      if(present(grid_at_center) ) center = grid_at_center
      if(center) then
         call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose, src_modulo )
      else
         nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
         allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
         do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
         enddo
         do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
         enddo
         call horiz_interp_bilinear_new ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
              verbose, src_modulo )
         deallocate(lon_src_1d,lat_src_1d)
      endif
   case ("bicubic")
      Interp%interp_method = BICUBIC
      center = .false.
      if(present(grid_at_center) ) center = grid_at_center
      if(center) then
        call horiz_interp_bicubic_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose, src_modulo )
      else
         nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
         allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
         do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
         enddo
         do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
         enddo
           call horiz_interp_bicubic_new ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
              verbose, src_modulo )
         deallocate(lon_src_1d,lat_src_1d)
      endif
   case ("spherical")
      Interp%interp_method = SPHERICAL
      nlon_in  = size(lon_in(:));  nlat_in  = size(lat_in(:))
      allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
      do i = 1, nlon_in
         lon_src(i,:) = lon_in(i)
      enddo
      do j = 1, nlat_in
         lat_src(:,j) = lat_in(j)
      enddo
      call horiz_interp_spherical_new ( Interp, lon_src, lat_src, lon_out, lat_out, &
           num_nbrs, max_dist, src_modulo)
      deallocate(lon_src, lat_src)
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, bicubic, spherical')
   end select

   !-----------------------------------------------------------------------
   Interp% HI_KIND_TYPE_ % is_allocated = .true.
   Interp%I_am_initialized = .true.

 end subroutine HORIZ_INTERP_NEW_1D_SRC_

!#######################################################################
!> @parblock
!! Wrapper to compute weights for conservative, bilinear, or spherical
!! interpolation.  Both input and output grids are represented as 2D arrays and
!! do not need to be rectangular.  Bicubic interpolation is not supported in this subroutine.
!! @endparblock
 subroutine HORIZ_INTERP_NEW_2D_ (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                 verbose, interp_method, num_nbrs, max_dist, &
                                 src_modulo, mask_in, mask_out, is_latlon_in, is_latlon_out  )
 type(horiz_interp_type), intent(inout)  :: Interp
   !< is a derived type that will be populated with interpolation weights and indices
 real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_in
   !< are the longitude coordinates of the 2D input grid
 real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lat_in
   !< are the latitude coordinates of the 2D input grid
 real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_out
   !< are the longitude coordinates of the 2D output grid
 real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lat_out
   !< are the latitude coordinates of the 2D output grid
 integer, intent(in), optional :: verbose
   !< sets the verbosity level
 character(len=*), intent(in), optional :: interp_method
   !< is the interpolation method; can be "conservative", "bilinear", "bicubic", or "spherical"; defaults to "bilinear"
 integer, intent(in), optional :: num_nbrs
   !< is the number of nearest neighbor for spherical interp; defaults to 4
 real(FMS_HI_KIND_), intent(in), optional :: max_dist
   !< is the maximum region of influence around output grid cell for spherical interp; defaults to 0.1 radians
 logical, intent(in), optional :: src_modulo
   !< is true if boundary condition along zonal boundary is cyclic
 real(FMS_HI_KIND_), intent(in), dimension(:,:), optional :: mask_in
   !< used in conservative interp to mask input grid cells (skip interp for masked cells)
 real(FMS_HI_KIND_), intent(out),dimension(:,:), optional :: mask_out
   !< will be populated with output grid cell areas computed from exchange grid for conserve interp
 logical, intent(in), optional :: is_latlon_in
   !< indicates if the input grid is a lat-lon grid; used for conservative interpolation
 logical, intent(in), optional :: is_latlon_out
   !< indicates if the output grid is a lat-lon grid; used for conservative interpolation

 logical  :: src_is_latlon, dst_is_latlon
 character(len=40) :: method
 integer, parameter  :: kindl = FMS_HI_KIND_ !< real kind size currently compiling
!-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'bilinear'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      if(PRESENT(is_latlon_in)) then
         src_is_latlon = is_latlon_in
      else
         src_is_latlon = is_lat_lon(lon_in, lat_in)
      end if
      if(PRESENT(is_latlon_out)) then
         dst_is_latlon = is_latlon_out
      else
         dst_is_latlon = is_lat_lon(lon_out, lat_out)
      end if
      if(src_is_latlon .AND. dst_is_latlon) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -0.0001_kindl) .or. ANY(mask_in > 1.0001_kindl)) then
                call mpp_error(FATAL, 'horiz_interp_conserve_new_2d(horiz_interp_conserve_mod):' // &
                                       ' input mask not between 0,1')
            endif
            allocate(Interp%HI_KIND_TYPE_%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%HI_KIND_TYPE_%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out(:,1), lat_out(1,:), &
              verbose=verbose )
      else if(src_is_latlon) then
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      else if(dst_is_latlon) then
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out(:,1), lat_out(1,:), &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if

   case ("spherical")
      Interp%interp_method = SPHERICAL
      call horiz_interp_spherical_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                    num_nbrs, max_dist, src_modulo )
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                        verbose, src_modulo )
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select

   Interp% HI_KIND_TYPE_ % is_allocated = .true.
   Interp%I_am_initialized = .true.

 end subroutine HORIZ_INTERP_NEW_2D_

!#######################################################################
 !> @parblock
 !! Wrapper to compute weights for conservative, bilinear, or spherical
 !! interpolation for any 2D input grids and rectangular output grids specified
 !! as 1D arrays.  The 1D array of input gridpoints will be used to construct
 !! 2D grids before calling either horiz_interp_conserve_new, horiz_interp_bilinear_new,
 !! or horiz_interp_spherical_new.
 !! @endparblock
 subroutine HORIZ_INTERP_NEW_1D_DST_ (Interp, lon_in, lat_in, lon_out, lat_out,   &
      verbose, interp_method, num_nbrs, max_dist, src_modulo, mask_in, mask_out, is_latlon_in )

   type(horiz_interp_type), intent(inout) :: Interp
     !< is a derived type that will be populated with interpolation weights and indices
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_in
     !< are the longitude coordinates of the 2D input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lat_in
     !< are the latitude coordinates of the 2D input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lon_out
     !< are the longitude coordinates of a rectangular output grid
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lat_out
     !< are the latitude coordinates of a rectangular output grid
   integer, intent(in), optional :: verbose
     !< sets the verbosity level
   character(len=*), intent(in), optional :: interp_method
     !< is the interpolation method; can be "conservative", "bilinear", or "spherical"; defaults to "bilinear"
   integer, intent(in), optional :: num_nbrs
     !< is the number of nearest neighbor for spherical interp; defaults to 4
   real(FMS_HI_KIND_), intent(in), optional :: max_dist
     !< is the maximum region of influence around output grid cell for spherical interp; defaults to 0.1 radians
   logical, intent(in), optional :: src_modulo
     !< is true if boundary condition along zonal boundary is cyclic
   real(FMS_HI_KIND_), intent(in), dimension(:,:), optional :: mask_in
     !< used in conservative interp to mask input grid cells (skip interp for masked cells)
   real(FMS_HI_KIND_), intent(out),dimension(:,:), optional :: mask_out
     !< will be populated with output grid cell areas computed from exchange grid for conserve interp
   logical, intent(in), optional :: is_latlon_in
     !< indicates if the input grid is a lat-lon grid; used for conservative interpolation

   character(len=40) :: method
   integer, parameter  :: kindl = FMS_HI_KIND_ !< real kind size currently compiling
   !-------------some local variables-----------------------------------------------
   integer  :: i, j, nlon_out, nlat_out
   real(FMS_HI_KIND_), dimension(:,:), allocatable :: lon_dst, lat_dst
   logical  :: src_is_latlon
   !-----------------------------------------------------------------------
   call horiz_interp_init

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
   case ("conservative")
      Interp%interp_method = CONSERVE
      if(PRESENT(is_latlon_in)) then
         src_is_latlon = is_latlon_in
      else
         src_is_latlon = is_lat_lon(lon_in, lat_in)
      end if

      if(src_is_latlon) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -0.0001_kindl) .or. ANY(mask_in > 1.0001_kindl)) &
                 call mpp_error(FATAL, &
                     'horiz_interp_conserve_new_1d_dst(horiz_interp_conserve_mod): input mask not between 0,1')
            allocate(Interp%HI_KIND_TYPE_%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%HI_KIND_TYPE_%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out, lat_out, &
              verbose=verbose)
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           verbose, src_modulo )
   case ("spherical")
      Interp%interp_method = SPHERICAL
      call horiz_interp_spherical_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           num_nbrs, max_dist, src_modulo)
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select

   deallocate(lon_dst,lat_dst)

   Interp% HI_KIND_TYPE_ % is_allocated = .true.
   Interp%I_am_initialized = .true.

 end subroutine HORIZ_INTERP_NEW_1D_DST_

!#######################################################################

 !> @parblock
 !! Interpolates 2D input data using the weights stored in Interp.
 !! @endparblock
 subroutine HORIZ_INTERP_BASE_2D_ ( Interp, data_in, data_out, verbose, &
                                   mask_in, mask_out, missing_value, missing_permit, &
                                   err_msg, new_missing_handle )
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
     !< contains the interpolation weights and indices
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
     !< is the data to interpolate
   real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
     !< will contain the interpolated data
   integer, intent(in), optional :: verbose
     !< sets the verbosity level
   real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
     !< is not use unless Interp%version=1 for conservative interp
   real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
     !< is not used unless Interp%version=1 for conservative interp
   real(FMS_HI_KIND_), intent(in), optional :: missing_value
     !< is the missing_value for bilinear or spherical interp
   integer, intent(in), optional :: missing_permit
     !< is the maximum number of missing values allowed for bilinear
   character(len=*), intent(out), optional :: err_msg
     !< will contain the error message
   logical, intent(in), optional :: new_missing_handle
     !< turns on alternative missing value handling for bilinear interp
!-----------------------------------------------------------------------
   if(present(err_msg)) err_msg = ''
   if(.not.Interp%I_am_initialized) then
     if(fms_error_handler('horiz_interp','The horiz_interp_type variable is not initialized',err_msg)) return
   endif

   select case(Interp%interp_method)
   case(CONSERVE)
      call horiz_interp_conserve(Interp,data_in, data_out, verbose, mask_in, mask_out)
   case(BILINEAR)
      call horiz_interp_bilinear(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit, new_missing_handle )
   case(BICUBIC)
      call horiz_interp_bicubic(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit )
   case(SPHERICAL)
      call horiz_interp_spherical(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value )
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, bicubic, spherical')
   end select

   return

 end subroutine HORIZ_INTERP_BASE_2D_

!#######################################################################

 !> @parblock
 !! Interpolates data for each 2D slice in data_in(:,:,n) for each n using the weights
 !! stored in Interp.
 !! @endparblock
 subroutine HORIZ_INTERP_BASE_3D_ ( Interp, data_in, data_out, verbose, mask_in, mask_out, &
      missing_value, missing_permit, err_msg  )
   !-----------------------------------------------------------------------
   !   overload of interface HORIZ_INTERP_BASE_2D_
   !   uses 3d arrays for data and mask
   !   this allows for multiple interpolations with one call
   !-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in)  :: Interp
     !< contains the interpolation weights and indices
   real(FMS_HI_KIND_), intent(in),  dimension(:,:,:)  :: data_in
     !< is the data to interpolate
   real(FMS_HI_KIND_), intent(out), dimension(:,:,:)  :: data_out
     !< will contain the interpolated data
   integer, intent(in), optional :: verbose
     !< sets the verbosity level
   real(FMS_HI_KIND_), intent(in),   dimension(:,:,:), optional :: mask_in
     !< is not use unless Interp%version=1 for conservative interp
   real(FMS_HI_KIND_), intent(out),  dimension(:,:,:), optional :: mask_out
     !< is not used unless Interp%version=1 for conservative interp
   real(FMS_HI_KIND_), intent(in), optional :: missing_value
     !< is the missing_value for bilinear or spherical interp
   integer, intent(in), optional :: missing_permit
     !< is the maximum number of missing values allowed for bilinear interp
   character(len=*), intent(out), optional :: err_msg
     !< will contain the error message
   !-----------------------------------------------------------------------
   integer :: n

   if(present(err_msg)) err_msg = ''
   if(.not.Interp%I_am_initialized) then
     if(fms_error_handler('horiz_interp','The horiz_interp_type variable is not initialized',err_msg)) return
   endif

   do n = 1, size(data_in,3)
      if (present(mask_in))then
         if(present(mask_out)) then
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), mask_out(:,:,n), &
                 missing_value, missing_permit )
         else
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
      else
         if(present(mask_out)) then
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_out=mask_out(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         else
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
     endif
   enddo

   return
!-----------------------------------------------------------------------
 end subroutine HORIZ_INTERP_BASE_3D_

!#######################################################################

 !> @parblock
 !! Compute weights and interpolate 2D data from a rectangular grid to a rectangular grid that
 !! are provided as 1D arrays.  The 1D arrays of gridpoints will be used to construct the 2D grids
 !! in horiz_interp_new that will be called to generate weights for interp_methods "conservative" (default),
 !! "bilinear", "bicubic", or "spherical" before calling horiz_interp to interpolate data.
 !! @endparblock
 subroutine HORIZ_INTERP_SOLO_1D_ ( data_in, lon_in, lat_in, lon_out, lat_out,   &
                                   data_out, verbose, mask_in, mask_out,         &
                                   interp_method, missing_value, missing_permit, &
                                   num_nbrs, max_dist,src_modulo, grid_at_center  )
!-----------------------------------------------------------------------
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
     !< is the data to be interpolated
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lon_in
     !< are the longitude coordinates of rectangular input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lat_in
     !< are the latitude coordinates of rectangular input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lon_out
     !< are the longitude coordinates of rectangular output grid
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lat_out
     !< are the latitude coordinates of rectangular output grid
   real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
     !< will contain the interpolated data
   integer, intent(in), optional :: verbose
     !< is the verbosity level
   real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
     !< used in conservative interp to mask input grid cells (skip interp for masked cells)
   real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
     !< will be populated with output grid cell areas computed from exchange grid for conserve interp
   character(len=*), intent(in), optional :: interp_method
     !< is the interp method; can be "conservative" (default), "bilinear", "bicubic", or "spherical"
   real(FMS_HI_KIND_), intent(in), optional :: missing_value
     !< is the missing_value for bilinear or spherical interp
   integer, intent(in), optional :: missing_permit
     !< is the maximum number of missing values allowed for bilinear interp
   integer, intent(in), optional :: num_nbrs
     !< is the number of nearest neighbor for spherical interp; defaults to 4
   real(FMS_HI_KIND_), intent(in), optional :: max_dist
     !< is the maximum region of influence around output grid cell for spherical interp; defaults to 0.1 radians
   logical, intent(in), optional :: src_modulo
     !< is true if boundary condition along zonal boundary is cyclic
   logical, intent(in), optional :: grid_at_center
     !< is true if the grid points are the grid center points for bilinear interp; if false, computes center points
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------
    call horiz_interp_init

    call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo, grid_at_center )

    call horiz_interp ( Interp, data_in, data_out, verbose,   &
                        mask_in, mask_out, missing_value, missing_permit )

    call horiz_interp_del ( Interp )
!-----------------------------------------------------------------------

 end subroutine HORIZ_INTERP_SOLO_1D_

!#######################################################################

 !> @parblock
 !! Subroutine to compute weights and interpolate 2D data where the input grid is
 !! a rectangular grid specified as 1D arrays and the output grid as 2D arrays.
 !! The 1D array of gridpoints will be used to construct the 2D grid in horiz_interp_new
 !! that will be called to generate weights for interp_methods "conservative" (default),
 !! "bilinear", "bicubic", or "spherical" before calling horiz_interp to interpolate data.
 !! @endparblock
 subroutine HORIZ_INTERP_SOLO_1D_SRC_ ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                       data_out, verbose, mask_in, mask_out,         &
                                       interp_method, missing_value, missing_permit, &
                                       num_nbrs, max_dist, src_modulo, grid_at_center )
!-----------------------------------------------------------------------
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
     !< is the data to be interpolated
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lon_in
     !< are the longitude coordinates of a rectangular input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:) :: lat_in
     !< are the latitude coordinates of a rectangular input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_out
     !< are the longitude coordinates of the 2D output grid
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lat_out
     !< are the latitude coordinates of the 2D output grid
   real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
     !< will contain the interpolated data
   integer, intent(in), optional :: verbose
     !< sets the verbosity level
   real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
     !< used in conservative interp to mask input grid cells (skip interp for masked cells)
   real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
     !< will be populated with output grid cell areas computed from exchange grid in conserve interp
   character(len=*), intent(in), optional :: interp_method
     !< is the interp method; can be "conservative" (default), "bilinear", "bicubic", or "spherical"
   real(FMS_HI_KIND_), intent(in), optional :: missing_value
     !< is the missing_value for bilinear or spherical interp
   integer, intent(in), optional :: missing_permit
     !< is the maximum number of missing values allowed for bilinear or spherical interp
   integer, intent(in), optional :: num_nbrs
     !< is the number of nearest neighbor for spherical interp; defaults to 4
   real(FMS_HI_KIND_), intent(in), optional :: max_dist
     !< is the maximum region of influence around output grid cell for spherical interp; defaults to 0.1 radians
   logical, intent(in), optional :: src_modulo
     !< is true if boundary condition along zonal boundary is cyclic
   logical, intent(in), optional :: grid_at_center
     !< is true if the grid points are the grid center points for bilinear interp; if false, computes center points

!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical  :: dst_is_latlon
   character(len=128)  :: method
!-----------------------------------------------------------------------
    call horiz_interp_init
    method = 'conservative'
    if(present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    if(trim(method) == 'conservative') dst_is_latlon = is_lat_lon(lon_out, lat_out)

    if(dst_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               grid_at_center, is_latlon_out = dst_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               grid_at_center, mask_in, mask_out, is_latlon_out = dst_is_latlon)

       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine HORIZ_INTERP_SOLO_1D_SRC_


!#######################################################################

 !> @parblock
 !! Compute weights and interpolate 2D by first calling horiz_interp_new to generate weights for
 !! interp_method "conservative", "bilinear" (default), or "spherical" before calling horiz_interp
 !! to interpolate data.
 !! @endparblock
 subroutine HORIZ_INTERP_SOLO_2D_ ( data_in, lon_in, lat_in, lon_out, lat_out, data_out, &
                                   verbose, mask_in, mask_out, interp_method, missing_value,&
                                   missing_permit, num_nbrs, max_dist, src_modulo  )
!-----------------------------------------------------------------------
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
     !< is the data to be interpolated
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_in
     !< are the longitude coordinates of 2D input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lat_in
     !< are the latitude coordinates of 2D input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_out
     !< are the longitude coordinates of the 2D output grid
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lat_out
     !< are the latitude coordinates of the 2D output grid
   real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
     !< will contain the interpolated data
   integer, intent(in), optional :: verbose
     !< sets the verbosity level
   real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
     !< used in conservative interp to mask input grid cells (skip interp for masked cells)
   real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
     !< will be populated with output grid cell areas computed from exchange grid for conserve interp
   character(len=*), intent(in), optional :: interp_method
     !< is the interp method; can be "conservative", "bilinear" (default), "bicubic", or "spherical"
   real(FMS_HI_KIND_), intent(in), optional :: missing_value
     !< is the missing_value for bilinear or spherical interp
   integer, intent(in), optional :: missing_permit
     !< is the maximum number of missing values allowed for bilinear interp
   integer, intent(in), optional :: num_nbrs
     !< is the number of nearest neighbor for spherical interp; defaults to 4
   real(FMS_HI_KIND_), intent(in), optional :: max_dist
     !< is the maximum region of influence around output grid cell for spherical interp; defaults to 0.1 radians
   logical, intent(in), optional :: src_modulo
     !< is true if boundary condition along zonal boundary is cyclic
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical  :: dst_is_latlon, src_is_latlon
   character(len=128)  :: method
!-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    src_is_latlon = .true.
    if(trim(method) == 'conservative') then
       dst_is_latlon = is_lat_lon(lon_out, lat_out)
       src_is_latlon = is_lat_lon(lon_in, lat_in)
    end if

    if(dst_is_latlon .and. src_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               is_latlon_in=dst_is_latlon, is_latlon_out = dst_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               mask_in, mask_out, &
                               is_latlon_in=dst_is_latlon, is_latlon_out = dst_is_latlon)
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine HORIZ_INTERP_SOLO_2D_

!#######################################################################

 !> @parblock
  !! Subroutine to compute weights and interpolate 2D data where the output grid is
 !! a rectangular grid specified as 1D arrays and the output grid as 2D arrays.
 !! The 1D array of gridpoints will be used to construct the 2D grid in horiz_interp_new
 !! that will be called to generate weights for interp_methods "conservative",
 !! "bilinear" (default), "bicubic", or "spherical" before calling horiz_interp to interpolate data.
 !! @endparblock
 subroutine HORIZ_INTERP_SOLO_1D_DST_ ( data_in, lon_in, lat_in, lon_out, lat_out, data_out,    &
                                       verbose, mask_in, mask_out,interp_method,missing_value, &
                                       missing_permit,  num_nbrs, max_dist, src_modulo)
!-----------------------------------------------------------------------
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
     !< is the data to be interpolated
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_in
     !< are the longitude coordinates of the 2D input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lat_in
     !< are the latitude coordinates of the 2D input grid
   real(FMS_HI_KIND_), intent(in),  dimension(:)  :: lon_out
   !< are the longitude coordinates of a lat-lon output grid
   real(FMS_HI_KIND_), intent(in),  dimension(:)  :: lat_out
     !< are the latitude coordinates of a lat-lon output grid
   real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
     !< will contain the interpolated data
   integer, intent(in), optional :: verbose
     !< sets the verbosity level
   real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
     !< used in conservative interp to mask input grid cells (skip interp for masked cells)
   real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
     !< will be populated with output grid cell areas computed from exchange grid for conserve interp
   character(len=*), intent(in), optional :: interp_method
     !< is the interpolation method; can be "conservative", "bilinear" (default), "bicubic", or "spherical"
   real(FMS_HI_KIND_), intent(in), optional :: missing_value
     !< is the missing_value for bilinear and spherical_interp
   integer, intent(in), optional :: missing_permit
     !< is the maximum number of missing values allowed for bilinear interp
   integer, intent(in), optional :: num_nbrs
     !< is the number of nearest neighbor for spherical interp; defaults to 4
   real(FMS_HI_KIND_), intent(in), optional :: max_dist
     !< is the maximum region of influence around output grid cell for spherical interp; defaults to 0.1 radians
   logical, intent(in), optional :: src_modulo
     !< is true if boundary condition along zonal boundary is cyclic
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical  :: src_is_latlon
   character(len=128)  :: method
!-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method
    src_is_latlon = .true.
    if(trim(method) == 'conservative') src_is_latlon = is_lat_lon(lon_in, lat_in)

    if(src_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               is_latlon_in = src_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               mask_in, mask_out, is_latlon_in = src_is_latlon)

       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine HORIZ_INTERP_SOLO_1D_DST_

!#######################################################################

!> Overloaded version of interface horiz_interp_solo_2
 subroutine HORIZ_INTERP_SOLO_OLD_ (data_in, wb, sb, dx, dy,  &
                                   lon_out, lat_out, data_out,  &
                                   verbose, mask_in, mask_out)

!-----------------------------------------------------------------------
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in !< Global input data stored from west to east
                                        !! (1st dimension), south to north (2nd dimension)
      real(FMS_HI_KIND_), intent(in)  :: wb !< Longitude (radians) that correspond to western-most
                                              !! boundary of grid box j=1 in array data_in
      real(FMS_HI_KIND_), intent(in)  :: sb !< Latitude (radians) that correspond to western-most
                                              !! boundary of grid box j=1 in array data_in
      real(FMS_HI_KIND_), intent(in)  :: dx !< Grid spacing (in radians) for the longitude axis
                                              !! (first dimension) for the input data
      real(FMS_HI_KIND_), intent(in)  :: dy !< Grid spacing (in radians) for the latitude axis
                                              !! (first dimension) for the input data
      real(FMS_HI_KIND_), intent(in),  dimension(:)  :: lon_out !< The longitude edges (in radians) for output
                                        !! data grid boxes. The values are for adjacent grid boxes
                                        !! and must increase in value. If there are MLON grid boxes
                                        !! there must be MLON+1 edge values
      real(FMS_HI_KIND_), intent(in),  dimension(:)  :: lat_out !< The latitude edges (in radians) for output
                                        !! data grid boxes. The values are for adjacent grid boxes
                                        !! and may increase or decrease in value. If there are NLAT
                                        !! grid boxes there must be NLAT+1 edge values
      real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out !< Output data on the output grid defined by grid box
   integer, intent(in), optional :: verbose
      real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
        !< used in conservative interp to mask input grid cells (skip interp for masked cells)
      real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
        !< will be populated with output grid cell areas computed from exchange grid for conserve interp
!-----------------------------------------------------------------------
     real(FMS_HI_KIND_), dimension(size(data_in,1)+1)  :: blon_in
     real(FMS_HI_KIND_), dimension(size(data_in,2)+1)  :: blat_in
     integer :: i, j, nlon_in, nlat_in
     real(FMS_HI_KIND_)  :: tpi
     integer, parameter  :: kindl = FMS_HI_KIND_ !< real size at compile time
!-----------------------------------------------------------------------
   call horiz_interp_init

   tpi = 2.0_kindl * real(pi, FMS_HI_KIND_)
   nlon_in = size(data_in,1)
   nlat_in = size(data_in,2)

   do i = 1, nlon_in+1
      blon_in(i) = wb + real(i-1, FMS_HI_KIND_)*dx
   enddo
      if (abs(blon_in(nlon_in+1)-blon_in(1)-tpi) < epsilon(blon_in)) &
              blon_in(nlon_in+1)=blon_in(1)+tpi

   do j = 2, nlat_in
      blat_in(j) = sb + real(j-1, FMS_HI_KIND_)*dy
   enddo
      blat_in(1)         = -0.5_kindl * real(pi, FMS_HI_KIND_)
      blat_in(nlat_in+1) =  0.5_kindl * real(pi, FMS_HI_KIND_)


   call horiz_interp_solo_1d (data_in, blon_in, blat_in,    &
                              lon_out, lat_out, data_out,   &
                              verbose, mask_in, mask_out    )

!-----------------------------------------------------------------------

 end subroutine HORIZ_INTERP_SOLO_OLD_

!#######################################################################

 !> @parblock
 !! Function to determine if a grid is a lat-lon grid.  A lat-lon grid
 !! is defined where the latitude values are constant for a longitudinal line
 !! and longitude values are constant along each latitudinal line.
 !! @endparblock
 function IS_LAT_LON_(lon, lat)
    real(FMS_HI_KIND_), dimension(:,:), intent(in) :: lon
      !< are the longitude points of a grid
    real(FMS_HI_KIND_), dimension(:,:), intent(in) :: lat
      !< are the latitude points of a grid
    logical  :: IS_LAT_LON_
    integer  :: i, j, nlon, nlat, num

    IS_LAT_LON_ = .true.
    nlon = size(lon,1)
    nlat = size(lon,2)
    LOOP_LAT: do j = 1, nlat
       do i = 2, nlon
          if(lat(i,j) .NE. lat(1,j)) then
             IS_LAT_LON_ = .false.
             exit LOOP_LAT
          end if
       end do
    end do LOOP_LAT

    if(IS_LAT_LON_) then
       LOOP_LON: do i = 1, nlon
          do j = 2, nlat
             if(lon(i,j) .NE. lon(i,1)) then
                IS_LAT_LON_ = .false.
                exit LOOP_LON
             end if
          end do
       end do LOOP_LON
    end if

    num = 0
    if(IS_LAT_LON_) num = 1
    call mpp_min(num)
    if(num == 1) then
       IS_LAT_LON_ = .true.
    else
       IS_LAT_LON_ = .false.
    end if

    return
 end function IS_LAT_LON_

 !> @parblock
 !! Subroutine to populate Interp with weights and indices read in from the weight file.
 !! Only supports weight files generated by fregrid for bilinear interpolation.
 !! @endparblock
 subroutine HORIZ_INTERP_READ_WEIGHTS_(Interp, weight_filename, lon_out, lat_out, lon_in, lat_in, &
                                       weight_file_source, interp_method, isw, iew, jsw, jew, nglon, nglat)
   type(horiz_interp_type), intent(inout) :: Interp
     !< will be populated with interpolation weights and indices read in from the weight file
   character(len=*), intent(in)  :: weight_filename
     !< is the weight file containing interpolation weights
   real(FMS_HI_KIND_), intent(in)  :: lat_out(:,:)
     !< are the latitude coordiantes of output grid
   real(FMS_HI_KIND_), intent(in)  :: lon_out(:,:)
     !< are the longitude coordinates of output grid
   real(FMS_HI_KIND_), intent(in)  :: lat_in(:)
     !< are the latitude coordinates of the input grid
   real(FMS_HI_KIND_), intent(in)  :: lon_in(:)
     !< are the longitude coordinates of the input grid
   character(len=*), intent(in)  :: weight_file_source
     !< Source of the weight file, only accepts "fregrid" at this time
   character(len=*), intent(in)  :: interp_method
     !< is the interpolation method, only accepts bilinear at this time
   integer, intent(in) :: isw
     !< starting x index in compute domain
   integer, intent(in) :: iew
     !< ending x index in compute domain
   integer, intent(in) :: jsw
     !< starting y index in compute domain
   integer, intent(in) :: jew
     !< ending y index in compute domain
   integer, intent(in) :: nglon
     !< is the number of longitudinal cells in the global domain
   integer, intent(in) :: nglat
     !< is the number of latitudinal cells in the global domain

   integer  :: i, j                 ! For do loops
   integer  :: nlon_in              ! Number of longitude in the data
   integer  :: nlat_in              ! Number of latitude in the data grid
   real(FMS_HI_KIND_), allocatable  :: lon_src_1d(:)        ! Center points of the longitude data grid
   real(FMS_HI_KIND_), allocatable  :: lat_src_1d(:)        ! Center points of the latitude data grid
   integer, parameter  :: kindl = FMS_HI_KIND_ ! real kind size currently compiling

   select case (trim(interp_method))
   case ("bilinear")
     !! This is to reproduce the behavior in horiz_interp_new
     !! The subroutine assumes that the data grid (lon_in, lat_in) are
     !! the edges and not the centers.
     !! Data_override passes in the edges, which are calculated using the axis_edges subroutine
     nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
     allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
     do i = 1, nlon_in
       lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
     enddo
     do j = 1, nlat_in
       lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
     enddo

     call horiz_interp_read_weights_bilinear(Interp, weight_filename, lon_out, lat_out, &
      lon_src_1d, lat_src_1d, weight_file_source, interp_method, &
      isw, iew, jsw, jew, nglon, nglat)
     deallocate(lon_src_1d,lat_src_1d)
   case default
     call mpp_error(FATAL, "Reading weight from file is not supported for the "//&
      trim(interp_method)//" method. It is currently only supported for bilinear")
   end select
 end subroutine HORIZ_INTERP_READ_WEIGHTS_
