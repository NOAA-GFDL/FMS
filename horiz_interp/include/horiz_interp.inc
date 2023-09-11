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
!> @addtogroup horiz_interp_mod
!> @{
  !> @brief Creates a 1D @ref horiz_interp_type with the given parameters
  subroutine HORIZ_INTERP_NEW_1D_ (Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                                  interp_method, num_nbrs, max_dist, src_modulo,     &
                                  grid_at_center, mask_in, mask_out)

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout)        :: Interp
    real(FMS_HI_KIND_), intent(in),  dimension(:)               :: lon_in , lat_in
    real(FMS_HI_KIND_), intent(in),  dimension(:)               :: lon_out, lat_out
    integer, intent(in),                 optional :: verbose
    character(len=*), intent(in),        optional :: interp_method
    integer, intent(in),                 optional :: num_nbrs
    real(FMS_HI_KIND_),    intent(in),                 optional :: max_dist
    logical, intent(in),                 optional :: src_modulo
    logical, intent(in),                 optional :: grid_at_center
    real(FMS_HI_KIND_), intent(in), dimension(:,:),    optional :: mask_in  !< dummy variable
    real(FMS_HI_KIND_), intent(inout),dimension(:,:),  optional :: mask_out !< dummy variable
    !-----------------------------------------------------------------------
    real(FMS_HI_KIND_), dimension(:,:), allocatable :: lon_src, lat_src, lon_dst, lat_dst
    real(FMS_HI_KIND_), dimension(:),   allocatable :: lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d
    integer                           :: i, j, nlon_in, nlat_in, nlon_out, nlat_out
    logical                           :: center
    character(len=40)                 :: method
    integer, parameter                :: kindl = FMS_HI_KIND_ !> real kind size currently compiling
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
       Interp%interp_method = SPHERICA
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

 subroutine HORIZ_INTERP_NEW_1D_SRC_ (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                     verbose, interp_method, num_nbrs, max_dist, &
                                     src_modulo, grid_at_center, mask_in, mask_out, is_latlon_out )

   type(horiz_interp_type), intent(inout)        :: Interp
   real(FMS_HI_KIND_), intent(in),  dimension(:)               :: lon_in , lat_in
   real(FMS_HI_KIND_), intent(in),  dimension(:,:)             :: lon_out, lat_out
   integer, intent(in),                 optional :: verbose
   character(len=*), intent(in),        optional :: interp_method
   integer, intent(in),                 optional :: num_nbrs  !< minimum number of neighbors
   real(FMS_HI_KIND_),    intent(in),                 optional :: max_dist
   logical, intent(in),                 optional :: src_modulo
   logical, intent(in),                 optional :: grid_at_center
   real(FMS_HI_KIND_), intent(in), dimension(:,:),    optional :: mask_in
   real(FMS_HI_KIND_), intent(out),dimension(:,:),    optional :: mask_out
   logical, intent(in),                 optional :: is_latlon_out

   real(FMS_HI_KIND_), dimension(:,:), allocatable :: lon_src, lat_src
   real(FMS_HI_KIND_), dimension(:),   allocatable :: lon_src_1d, lat_src_1d
   integer                           :: i, j, nlon_in, nlat_in
   character(len=40)                 :: method
   logical                           :: center
   logical                           :: dst_is_latlon
   integer, parameter                :: kindl = FMS_HI_KIND_ !< real kind size currently compiling
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
      Interp%interp_method = SPHERICA
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

 subroutine HORIZ_INTERP_NEW_2D_ (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                 verbose, interp_method, num_nbrs, max_dist, &
                                 src_modulo, mask_in, mask_out, is_latlon_in, is_latlon_out  )
 type(horiz_interp_type), intent(inout)     :: Interp
 real(FMS_HI_KIND_), intent(in),  dimension(:,:)          :: lon_in , lat_in
 real(FMS_HI_KIND_), intent(in),  dimension(:,:)          :: lon_out, lat_out
 integer, intent(in),              optional :: verbose
 character(len=*), intent(in),     optional :: interp_method
 integer, intent(in),              optional :: num_nbrs
 real(FMS_HI_KIND_),    intent(in),              optional :: max_dist
 logical, intent(in),              optional :: src_modulo
 real(FMS_HI_KIND_), intent(in), dimension(:,:), optional :: mask_in
 real(FMS_HI_KIND_), intent(out),dimension(:,:), optional :: mask_out
 logical, intent(in),              optional :: is_latlon_in, is_latlon_out
 logical           :: src_is_latlon, dst_is_latlon
 character(len=40) :: method
 integer, parameter                :: kindl = FMS_HI_KIND_ !< real kind size currently compiling
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
      Interp%interp_method = SPHERICA
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
 subroutine HORIZ_INTERP_NEW_1D_DST_ (Interp, lon_in, lat_in, lon_out, lat_out,   &
      verbose, interp_method, num_nbrs, max_dist, src_modulo, mask_in, mask_out, is_latlon_in )
   type(horiz_interp_type), intent(inout)     :: Interp
   real(FMS_HI_KIND_), intent(in),  dimension(:,:)          :: lon_in , lat_in
   real(FMS_HI_KIND_), intent(in),  dimension(:)            :: lon_out, lat_out
   integer, intent(in),              optional :: verbose
   character(len=*), intent(in),     optional :: interp_method
   integer, intent(in),              optional :: num_nbrs
   real(FMS_HI_KIND_),    intent(in),              optional :: max_dist
   logical, intent(in),              optional :: src_modulo
   real(FMS_HI_KIND_), intent(in), dimension(:,:), optional :: mask_in
   real(FMS_HI_KIND_), intent(out),dimension(:,:), optional :: mask_out
   logical, intent(in),              optional :: is_latlon_in

   character(len=40) :: method
   integer, parameter                :: kindl = FMS_HI_KIND_ !< real kind size currently compiling
   !-------------some local variables-----------------------------------------------
   integer                           :: i, j, nlon_out, nlat_out
   real(FMS_HI_KIND_), dimension(:,:), allocatable :: lon_dst, lat_dst
   logical                           :: src_is_latlon
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
      Interp%interp_method = SPHERICA
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

 subroutine HORIZ_INTERP_BASE_2D_ ( Interp, data_in, data_out, verbose, &
                                   mask_in, mask_out, missing_value, missing_permit, &
                                   err_msg, new_missing_handle )
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
      real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
      real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
      real(FMS_HI_KIND_), intent(in),                   optional :: missing_value
      integer, intent(in),                optional :: missing_permit
   character(len=*), intent(out),         optional :: err_msg
      logical, intent(in),                optional :: new_missing_handle
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
   case(SPHERICA)
      call horiz_interp_spherical(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value )
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, bicubic, spherical')
   end select

   return

 end subroutine HORIZ_INTERP_BASE_2D_

!#######################################################################

 !> Overload of interface HORIZ_INTERP_BASE_2D_
 !! uses 3d arrays for data and mask
 !! this allows for multiple interpolations with one call
 subroutine HORIZ_INTERP_BASE_3D_ ( Interp, data_in, data_out, verbose, mask_in, mask_out, &
      missing_value, missing_permit, err_msg  )
   !-----------------------------------------------------------------------
   !   overload of interface HORIZ_INTERP_BASE_2D_
   !   uses 3d arrays for data and mask
   !   this allows for multiple interpolations with one call
   !-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in)           :: Interp
   real(FMS_HI_KIND_), intent(in),  dimension(:,:,:)            :: data_in
   real(FMS_HI_KIND_), intent(out), dimension(:,:,:)            :: data_out
   integer, intent(in),                  optional :: verbose
   real(FMS_HI_KIND_), intent(in),   dimension(:,:,:), optional :: mask_in
   real(FMS_HI_KIND_), intent(out),  dimension(:,:,:), optional :: mask_out
   real(FMS_HI_KIND_), intent(in),                     optional :: missing_value
   integer, intent(in),                  optional :: missing_permit
   character(len=*), intent(out),        optional :: err_msg
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

!> Interpolates from a rectangular grid to rectangular grid.
!! interp_method can be the value conservative, bilinear or spherical.
!! horiz_interp_new don't need to be called before calling this routine.
 subroutine HORIZ_INTERP_SOLO_1D_ ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                   data_out, verbose, mask_in, mask_out,         &
                                   interp_method, missing_value, missing_permit, &
                                   num_nbrs, max_dist,src_modulo, grid_at_center  )
!-----------------------------------------------------------------------
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
      real(FMS_HI_KIND_), intent(in),  dimension(:)   :: lon_in , lat_in
      real(FMS_HI_KIND_), intent(in),  dimension(:)   :: lon_out, lat_out
      real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
      real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(FMS_HI_KIND_), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(FMS_HI_KIND_), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center
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

!> Interpolates from a uniformly spaced grid to any output grid.
!! interp_method can be the value "onservative","bilinear" or "spherical".
!! horiz_interp_new don't need to be called before calling this routine.
 subroutine HORIZ_INTERP_SOLO_1D_SRC_ ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                       data_out, verbose, mask_in, mask_out,         &
                                       interp_method, missing_value, missing_permit, &
                                       num_nbrs, max_dist, src_modulo, grid_at_center )
!-----------------------------------------------------------------------
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
      real(FMS_HI_KIND_), intent(in),  dimension(:)   :: lon_in , lat_in
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_out, lat_out
      real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
      real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(FMS_HI_KIND_), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(FMS_HI_KIND_), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center

!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: dst_is_latlon
   character(len=128)       :: method
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

!> Interpolates from any grid to any grid. interp_method should be "spherical"
!! horiz_interp_new don't need to be called before calling this routine.
 subroutine HORIZ_INTERP_SOLO_2D_ ( data_in, lon_in, lat_in, lon_out, lat_out, data_out, &
                                   verbose, mask_in, mask_out, interp_method, missing_value,&
                                   missing_permit, num_nbrs, max_dist, src_modulo  )
!-----------------------------------------------------------------------
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_in , lat_in
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_out, lat_out
      real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
      real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(FMS_HI_KIND_), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(FMS_HI_KIND_), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: dst_is_latlon, src_is_latlon
   character(len=128)       :: method
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

!>   interpolates from any grid to rectangular longitude/latitude grid.
!!   interp_method should be "spherical".
!!   horiz_interp_new don't need to be called before calling this routine.
 subroutine HORIZ_INTERP_SOLO_1D_DST_ ( data_in, lon_in, lat_in, lon_out, lat_out, data_out,    &
                                       verbose, mask_in, mask_out,interp_method,missing_value, &
                                       missing_permit,  num_nbrs, max_dist, src_modulo)
!-----------------------------------------------------------------------
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: data_in
      real(FMS_HI_KIND_), intent(in),  dimension(:,:) :: lon_in , lat_in
      real(FMS_HI_KIND_), intent(in),  dimension(:)   :: lon_out, lat_out
      real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
      real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(FMS_HI_KIND_), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(FMS_HI_KIND_), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: src_is_latlon
   character(len=128)       :: method
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
      real(FMS_HI_KIND_), intent(in)                  :: wb !< Longitude (radians) that correspond to western-most
                                              !! boundary of grid box j=1 in array data_in
      real(FMS_HI_KIND_), intent(in)                  :: sb !< Latitude (radians) that correspond to western-most
                                              !! boundary of grid box j=1 in array data_in
      real(FMS_HI_KIND_), intent(in)                  :: dx !< Grid spacing (in radians) for the longitude axis
                                              !! (first dimension) for the input data
      real(FMS_HI_KIND_), intent(in)                  :: dy !< Grid spacing (in radians) for the latitude axis
                                              !! (first dimension) for the input data
      real(FMS_HI_KIND_), intent(in),  dimension(:)   :: lon_out !< The longitude edges (in radians) for output
                                        !! data grid boxes. The values are for adjacent grid boxes
                                        !! and must increase in value. If there are MLON grid boxes
                                        !! there must be MLON+1 edge values
      real(FMS_HI_KIND_), intent(in),  dimension(:)   :: lat_out !< The latitude edges (in radians) for output
                                        !! data grid boxes. The values are for adjacent grid boxes
                                        !! and may increase or decrease in value. If there are NLAT
                                        !! grid boxes there must be NLAT+1 edge values
      real(FMS_HI_KIND_), intent(out), dimension(:,:) :: data_out !< Output data on the output grid defined by grid box
   integer, intent(in),                   optional :: verbose
      real(FMS_HI_KIND_), intent(in),   dimension(:,:), optional :: mask_in
      real(FMS_HI_KIND_), intent(out),  dimension(:,:), optional :: mask_out
!-----------------------------------------------------------------------
     real(FMS_HI_KIND_), dimension(size(data_in,1)+1)  :: blon_in
     real(FMS_HI_KIND_), dimension(size(data_in,2)+1)  :: blat_in
     integer :: i, j, nlon_in, nlat_in
     real(FMS_HI_KIND_)    :: tpi
     integer, parameter    :: kindl = FMS_HI_KIND_ !< real size at compile time
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


 !####################################################################
 function IS_LAT_LON_(lon, lat)
    real(FMS_HI_KIND_), dimension(:,:), intent(in) :: lon, lat
    logical                          :: IS_LAT_LON_
    integer                          :: i, j, nlon, nlat, num

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
!> @}
