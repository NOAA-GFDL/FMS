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
!> @defgroup topography_mod topography_mod
!> @ingroup topography
!> @brief Routines for creating land surface topography fields and land-water masks
!! for latitude-longitude grids.
!> @author Bruce Wyman
!!
!! This module generates realistic mountains and land-water masks
!! on a specified latitude-longitude grid by interpolating from the
!! 1/6 degree Navy mean topography and percent water data sets.
!! The fields that can be generated are mean and standard deviation
!! of topography within the specified grid boxes; and land-ocean (or
!! water) mask and land-ocean (or water) fractional area.
!!
!! The interpolation scheme conserves the area-weighted average
!! of the input data by using module horiz_interp.
!!
!! The interfaces get_gaussian_topog and gaussian_topog_init are documented
!! in @ref gaussian_topog_mod

!#######################################################################

 !> @brief Returns a "realistic" mean surface height field.
 !!
 !> Returns realistic mountains on a latitude-longtude grid.
 !! The returned field is the mean topography for the given grid boxes.
 !! Computed using a conserving area-weighted interpolation.
 !! The current input data set is the 1/6 degree Navy mean topography.
 !!
 !! @returns A logical value of true is returned if the surface height field was successfully
 !! created. A value of false may be returned if the input topography data set was not readable.
 !!
 !! @throws FATAL, shape(zmean) is not equal to (/size(blon)-1,size(blat)-1/)
 !! Check the input grid size and output field size.

 function GET_TOPOG_MEAN_1D_(blon, blat, zmean)

   real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blon !< Longitude (radians) at grid box boundaries
   real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blat !< Latitude (radians) at grid box boundaries
   real(kind=FMS_TOP_KIND_), intent(out), dimension(:,:) :: zmean !< Mean surface height(meters). Size must be
                                              !! size(blon)-1 by size(blat)-1
   logical :: GET_TOPOG_MEAN_1D_
   integer, parameter :: lkind = FMS_TOP_KIND_

!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(zmean(:,:)) /= (/size(blon(:))-1,size(blat(:))-1/)) ) &
        call error_mesg('get_topog_mean_1d','shape(zmean) is not&
            & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   GET_TOPOG_MEAN_1D_ = open_topog_file()

   if (GET_TOPOG_MEAN_1D_) call interp_topog(blon, blat, zmean)

!-----------------------------------------------------------------------

 end function GET_TOPOG_MEAN_1D_

!############################################################

 function GET_TOPOG_MEAN_2D_(blon, blat, zmean)

   real(kind=FMS_TOP_KIND_), intent(in),  dimension(:,:) :: blon, blat
   real(kind=FMS_TOP_KIND_), intent(out), dimension(:,:) :: zmean
   logical :: GET_TOPOG_MEAN_2D_
   integer, parameter :: lkind = FMS_TOP_KIND_
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(zmean(:,:)) /= (/size(blon,1)-1,size(blon,2)-1/)) .or. &
        any(shape(zmean(:,:)) /= (/size(blat,1)-1,size(blat,2)-1/)) ) &
        call error_mesg('get_topog_mean_2d','shape(zmean) is not&
            & equal to (/size(blon,1)-1,size(blon,2)-1/))', FATAL)

   GET_TOPOG_MEAN_2D_ = open_topog_file()

   if (GET_TOPOG_MEAN_2D_) call interp_topog(blon, blat, zmean)
!-----------------------------------------------------------------------

 end function GET_TOPOG_MEAN_2D_

!#######################################################################

 !> @brief Returns a standard deviation of higher resolution topography with
 !! the given model grid boxes.
 !!
 !> Returns the standard deviation of the "finer" input topography data set,
 !! currently the Navy 1/6 degree mean topography data, within the
 !! boundaries of the given input grid.
 !!
 !! @returns A logical value of true if the output field was successfully created and false
 !! if the input topography data set was not readable.
 function GET_TOPOG_STDEV_1D_(blon, blat, stdev)

   real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blon !< Longitude (radians) at grid box boundaries
   real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blat !< Latitude (radians) at grid box boundaries
   real(kind=FMS_TOP_KIND_), intent(out), dimension(:,:) :: stdev !< The standard deviation of surface height (in
                                     !! meters) within given input model grid boxes. Size must be
                                     !! size(blon)-1 by size(blat)-1
   logical :: GET_TOPOG_STDEV_1D_
   integer, parameter :: lkind = FMS_TOP_KIND_
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(stdev(:,:)) /= (/size(blon(:))-1.0_lkind,size(blat(:))-1.0_lkind/)) ) &
       call error_mesg('get_topog_stdev','shape(stdev) is not&
            & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   GET_TOPOG_STDEV_1D_ = open_topog_file()

   if (GET_TOPOG_STDEV_1D_) call interp_topog(blon, blat, &
              stdev, flag=COMPUTE_STDEV)

!-----------------------------------------------------------------------

 end function GET_TOPOG_STDEV_1D_

!#######################################################################

 function GET_TOPOG_STDEV_2D_(blon, blat, stdev)

   real(FMS_TOP_KIND_), intent(in),  dimension(:,:) :: blon, blat
   real(FMS_TOP_KIND_), intent(out), dimension(:,:) :: stdev
   logical :: GET_TOPOG_STDEV_2D_
   integer, parameter :: lkind = FMS_TOP_KIND_
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(stdev(:,:)) /= (/size(blon,1)-1.0_lkind,size(blon,2)-1.0_lkind/)) .or. &
        any(shape(stdev(:,:)) /= (/size(blat,1)-1.0_lkind,size(blat,2)-1.0_lkind/)) ) &
        call error_mesg('get_topog_stdev_2d','shape(stdev) is not&
            & equal to (/size(blon,1)-1,size(blon,2)-1/))', FATAL)

   GET_TOPOG_STDEV_2D_ = open_topog_file()

   if (GET_TOPOG_STDEV_2D_) call interp_topog(blon, blat, &
              stdev, flag=COMPUTE_STDEV)
!-----------------------------------------------------------------------

 end function GET_TOPOG_STDEV_2D_

!#######################################################################

 !> @brief Returns fractional area covered by ocean in a grid box.
 !> @returns A logical value of true if the output field was successfully created. A value of false
 !! may be returned if the Navy 1/6 degree percent water data set was not readable.
 function GET_OCEAN_FRAC_1D_(blon, blat, ocean_frac)

 real(FMS_TOP_KIND_), intent(in),  dimension(:)   :: blon !< Longitude (radians) at grid box boundaries
 real(FMS_TOP_KIND_), intent(in),  dimension(:)   :: blat !< Latitude (radians) at grid box boundaries
 real(FMS_TOP_KIND_), intent(out), dimension(:,:) :: ocean_frac !< The fractional amount (0-1) of ocean in a grid
                                     !! box. The size must be size(blon)-1 by size(blat)-1
 logical :: GET_OCEAN_FRAC_1D_
 integer, parameter :: lkind = FMS_TOP_KIND_
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(ocean_frac(:,:)) /= (/size(blon(:))-1.0_lkind,size(blat(:))-1.0_lkind/)) ) &
        call error_mesg('get_ocean_frac','shape(ocean_frac) is not&
                 & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   GET_OCEAN_FRAC_1D_ = open_water_file()
   if(GET_OCEAN_FRAC_1D_)  call interp_water( blon, blat, &
                ocean_frac, do_ocean=.true. )

!-----------------------------------------------------------------------

 end function GET_OCEAN_FRAC_1D_

!#######################################################################

 function GET_OCEAN_FRAC_2D_(blon, blat, ocean_frac)

 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:,:) :: blon, blat
 real(kind=FMS_TOP_KIND_), intent(out), dimension(:,:) :: ocean_frac
 logical :: GET_OCEAN_FRAC_2D_
 integer, parameter :: lkind = FMS_TOP_KIND_
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(ocean_frac(:,:)) /= (/size(blon,1)-1.0_lkind,size(blon,2)-1.0_lkind/)) .or. &
        any(shape(ocean_frac(:,:)) /= (/size(blat,1)-1.0_lkind,size(blat,2)-1.0_lkind/)) ) &
        call error_mesg('get_ocean_frac_2d','shape(ocean_frac) is not&
            & equal to (/size(blon,1)-1,size(blon,2)-1/))', FATAL)

   GET_OCEAN_FRAC_2D_ = open_water_file()
   if(GET_OCEAN_FRAC_2D_)  call interp_water( blon, blat, &
                ocean_frac, do_ocean=.true. )

!-----------------------------------------------------------------------

 end function GET_OCEAN_FRAC_2D_

!#######################################################################

 !> @brief Returns a land-ocean mask in a grid box.
 !> @returns A logical value of true if the output field was successfully created. A value of false
 !! may be returned if the Navy 1/6 degree percent water data set was not readable.
 function GET_OCEAN_MASK_1D_(blon, blat, ocean_mask)

 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blon !< Longitude (radians) at grid box boundaries
 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blat !< Latitude (radians) at grid box boundaries
 logical, intent(out), dimension(:,:) :: ocean_mask !< Mask for ocean in a grid box.
                                                 !! The size must be size(blon)-1 by size(blat)-1
 logical :: GET_OCEAN_MASK_1D_
 real(kind=FMS_TOP_KIND_), dimension(size(ocean_mask,1),size(ocean_mask,2)) :: ocean_frac
 integer, parameter :: lkind = FMS_TOP_KIND_
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if (get_ocean_frac(blon, blat, ocean_frac) ) then
     where (ocean_frac > 0.50_lkind)
       ocean_mask = .true.
     elsewhere
       ocean_mask = .false.
     end where
     GET_OCEAN_MASK_1D_ = .true.
   else
    GET_OCEAN_MASK_1D_ = .false.
   endif
!-----------------------------------------------------------------------

 end function GET_OCEAN_MASK_1D_

!#######################################################################

 function GET_OCEAN_MASK_2D_(blon, blat, ocean_mask)

 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:,:) :: blon, blat
 logical, intent(out), dimension(:,:) :: ocean_mask
 logical :: GET_OCEAN_MASK_2D_
 real(kind=FMS_TOP_KIND_), dimension(size(ocean_mask,1),size(ocean_mask,2)) :: ocean_frac
 integer, parameter :: lkind = FMS_TOP_KIND_
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if (get_ocean_frac(blon, blat, ocean_frac) ) then
     where (ocean_frac > 0.50_lkind)
       ocean_mask = .true.
     elsewhere
       ocean_mask = .false.
     end where
     GET_OCEAN_MASK_2D_ = .true.
   else
     GET_OCEAN_MASK_2D_ = .false.
   endif

!-----------------------------------------------------------------------

 end function GET_OCEAN_MASK_2D_

 !> @brief Returns the percent of water in a grid box.
 !> @returns A logical value of true if the output field was successfully created. A value of false
 !! may be returned if the Navy 1/6 degree percent water data set was not readable.
 !!
 !! @throws FATAL, shape(water_frac) is not equal to (/size(blon)-1,size(blat)-1/)
 !! Check the input grid size and output field size.
 function GET_WATER_FRAC_1D_(blon, blat, water_frac)
 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blon !< The longitude (in radians) at grid box boundaries.
 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blat !< The latitude (in radians) at grid box boundaries.
 real(kind=FMS_TOP_KIND_), intent(out), dimension(:,:) :: water_frac !< The fractional amount (0 to 1) of water in a
                          !! grid box. The size of this field must be size(blon)-1 by size(blat)-1.
 logical :: GET_WATER_FRAC_1D_
 integer, parameter :: lkind = FMS_TOP_KIND_

!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(water_frac(:,:)) /= (/size(blon(:))-1.0_lkind,size(blat(:))-1.0_lkind/)) ) &
        call error_mesg('get_water_frac_1d','shape(water_frac) is not&
                 & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   GET_WATER_FRAC_1D_ = open_water_file()
   if(GET_WATER_FRAC_1D_) call interp_water( blon, blat, water_frac )

!-----------------------------------------------------------------------

 end function GET_WATER_FRAC_1D_

!#######################################################################

 function GET_WATER_FRAC_2D_(blon, blat, water_frac)

 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:,:)   :: blon !< The longitude (in radians) at grid box boundaries.
 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:,:)   :: blat !< The latitude (in radians) at grid box boundaries.
 real(kind=FMS_TOP_KIND_), intent(out), dimension(:,:) :: water_frac !< The fractional amount (0 to 1) of water in a
                          !! grid box. The size of this field must be size(blon)-1 by size(blat)-1.
 logical :: GET_WATER_FRAC_2D_
 integer, parameter :: lkind = FMS_TOP_KIND_

!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(water_frac(:,:)) /= (/size(blon,1)-1.0_lkind,size(blon,2)-1.0_lkind/)) .or. &
        any(shape(water_frac(:,:)) /= (/size(blat,1)-1.0_lkind,size(blat,2)-1.0_lkind/)) ) &
        call error_mesg('get_water_frac_2d','shape(water_frac) is not&
            & equal to (/size(blon,1)-1,size(blon,2)-1/))', FATAL)

   GET_WATER_FRAC_2D_ = open_water_file()
   if(GET_WATER_FRAC_2D_) call interp_water( blon, blat, water_frac )

!-----------------------------------------------------------------------

 end function GET_WATER_FRAC_2D_

!#######################################################################

 !> @brief Returns a land-water mask in the given model grid boxes.
 !> @return A logical value of TRUE is returned if the output field
 !! was successfully created. A value of FALSE may be returned
 !! if the Navy 1/6 degree percent water data set was not readable.
 function GET_WATER_MASK_1D_(blon, blat, water_mask)

 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blon !< The longitude (in radians) at grid box boundaries.
 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:)   :: blat !< The latitude (in radians) at grid box boundaries.
 logical, intent(out), dimension(:,:) :: water_mask !< A binary mask for water (true) or land (false).
                                     !! The size of this field must be size(blon)-1 by size(blat)-1.
 logical :: GET_WATER_MASK_1D_

 real(kind=FMS_TOP_KIND_), dimension(size(water_mask,1),size(water_mask,2)) :: water_frac
 integer, parameter :: lkind = FMS_TOP_KIND_
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if (get_water_frac(blon, blat, water_frac) ) then
     where (water_frac > 0.50_lkind)
       water_mask = .true.
     elsewhere
       water_mask = .false.
     end where
     GET_WATER_MASK_1D_ = .true.
   else
    GET_WATER_MASK_1D_ = .false.
   endif
!-----------------------------------------------------------------------

 end function GET_WATER_MASK_1D_

!#######################################################################

 function GET_WATER_MASK_2D_(blon, blat, water_mask)

 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:,:)   :: blon !< The longitude (in radians) at grid box boundaries.
 real(kind=FMS_TOP_KIND_), intent(in),  dimension(:,:)   :: blat !< The latitude (in radians) at grid box boundaries.
 logical, intent(out), dimension(:,:) :: water_mask !< A binary mask for water (true) or land (false).
                                     !! The size of this field must be size(blon)-1 by size(blat)-1.
 logical :: GET_WATER_MASK_2D_
 real(kind=FMS_TOP_KIND_), dimension(size(water_mask,1),size(water_mask,2)) :: water_frac
 integer, parameter :: lkind = FMS_TOP_KIND_
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if (get_water_frac(blon, blat, water_frac) ) then
     where (water_frac > 0.50_lkind)
       water_mask = .true.
     elsewhere
       water_mask = .false.
     end where
     GET_WATER_MASK_2D_ = .true.
   else
     GET_WATER_MASK_2D_ = .false.
   endif

!-----------------------------------------------------------------------

 end function GET_WATER_MASK_2D_

!#######################################################################

subroutine INTERP_TOPOG_1D_( blon, blat, zout, flag)
real(kind=FMS_TOP_KIND_), intent(in)  :: blon(:), blat(:)
real(kind=FMS_TOP_KIND_), intent(out) :: zout(:,:)
integer, intent(in), optional :: flag

real(kind=FMS_TOP_KIND_) :: xdat(ipts+1), ydat(jpts+1)
real(kind=FMS_TOP_KIND_) :: zdat(ipts,jpts)
real(kind=FMS_TOP_KIND_) :: zout2(size(zout,1),size(zout,2))
integer, parameter       :: lkind = FMS_TOP_KIND_

   call input_data(TOPOG_INDEX, xdat, ydat, zdat)

   call horiz_interp( zdat, xdat, ydat, blon, blat, zout )

! compute standard deviation if necessary.
   if (present(flag)) then
      if (flag == COMPUTE_STDEV) then
          zdat = zdat*zdat
          call horiz_interp ( zdat, xdat, ydat, blon, blat, zout2 )
          zout = zout2 - zout*zout
           where (zout > 0.0_lkind)
             zout = sqrt ( zout )
           elsewhere
             zout = 0.0_lkind
           endwhere
       endif
    endif

end subroutine INTERP_TOPOG_1D_

!#######################################################################

subroutine INTERP_TOPOG_2D_( blon, blat, zout, flag )
real(kind=FMS_TOP_KIND_), intent(in)  :: blon(:,:), blat(:,:)
real(kind=FMS_TOP_KIND_), intent(out) :: zout(:,:)
integer, intent(in), optional :: flag

real(kind=FMS_TOP_KIND_) :: xdat(ipts+1), ydat(jpts+1)
real(kind=FMS_TOP_KIND_) :: zdat(ipts,jpts)
real(kind=FMS_TOP_KIND_) :: zout2(size(zout,1),size(zout,2))
integer, parameter       :: lkind = FMS_TOP_KIND_
integer :: js, je
type (horiz_interp_type) :: Interp

    call input_data(TOPOG_INDEX, xdat, ydat, zdat)

    call find_indices( minval(blat), maxval(blat), ydat, js, je )

    call horiz_interp_new ( Interp, xdat, ydat(js:je+1), blon, blat )
    call horiz_interp     ( Interp, zdat(:,js:je), zout )

! compute standard deviation if necessary
    if (present(flag)) then
       if (flag == COMPUTE_STDEV) then
           zdat = zdat*zdat
           call horiz_interp ( Interp, zdat(:,js:je), zout2 )
           zout = zout2 - zout*zout
           where (zout > 0.0_lkind)
             zout = sqrt ( zout )
           elsewhere
             zout = 0.0_lkind
           endwhere
       endif
    endif

    call horiz_interp_del( Interp )

end subroutine INTERP_TOPOG_2D_

!#######################################################################

subroutine FIND_INDICES_( ybeg, yend, ydat, js, je )
real(kind=FMS_TOP_KIND_),    intent(in)  :: ybeg, yend, ydat(:)
integer, intent(out) :: js, je
integer :: j
integer, parameter :: lkind = FMS_TOP_KIND_

   js = 1
   do j = 1, size(ydat(:))-1.0_lkind
      if (ybeg >= ydat(j) .and. ybeg <= ydat(j+1.0_lkind)) then
         js = j
         exit
      endif
   enddo

   je = size(ydat(:))-1.0_lkind
   do j = js, size(ydat(:))-1.0_lkind
      if (yend >= ydat(j) .and. yend <= ydat(j+1.0_lkind)) then
         je = j
         exit
      endif
   enddo

!print '(a,i2,2(a,f10.5),2(a,i4))', "PE=",mpp_pe(),"  phs=",ybeg,"  phn=",yend,"  js=",js,"  je=",je

end subroutine FIND_INDICES_

!#######################################################################
subroutine INPUT_DATA_( indx, xdat, ydat, zdat )
integer, intent(in) :: indx
real(kind=FMS_TOP_KIND_), intent(out) :: xdat(ipts+1), ydat(jpts+1), zdat(ipts,jpts)

if( file_is_opened(indx) ) then
    call read_data(fileobj(indx), 'xdat', xdat)
    call read_data(fileobj(indx), 'ydat', ydat)
    call read_data(fileobj(indx), 'zdat', zdat)
endif

end subroutine INPUT_DATA_

!#######################################################################

subroutine INTERP_WATER_1D_( blon, blat, zout, do_ocean )
real(kind=FMS_TOP_KIND_), intent(in)  :: blon(:), blat(:)
real(kind=FMS_TOP_KIND_), intent(out) :: zout(:,:)
logical, intent(in), optional :: do_ocean
real(kind=FMS_TOP_KIND_) :: xdat(ipts+1), ydat(jpts+1), zdat(ipts,jpts)
    call input_data ( WATER_INDEX, xdat, ydat, zdat )

! only use designated ocean points
if (present(do_ocean)) then
    if (do_ocean) call determine_ocean_points (zdat)
endif

! interpolate onto output grid
call horiz_interp ( zdat, xdat, ydat, blon, blat, zout )

end subroutine INTERP_WATER_1D_

!#######################################################################

subroutine INTERP_WATER_2D_( blon, blat, zout, do_ocean )
real(kind=FMS_TOP_KIND_), intent(in)  :: blon(:,:), blat(:,:)
real(kind=FMS_TOP_KIND_), intent(out) :: zout(:,:)
logical, intent(in), optional :: do_ocean
real(kind=FMS_TOP_KIND_) :: xdat(ipts+1), ydat(jpts+1), zdat(ipts,jpts)

call input_data ( WATER_INDEX, xdat, ydat, zdat )

! only use designated ocean points
if (present(do_ocean)) then
    if (do_ocean) call determine_ocean_points (zdat)
endif

! interpolate onto output grid
call horiz_interp ( zdat, xdat, ydat, blon, blat, zout )

end subroutine INTERP_WATER_2D_

!#######################################################################

subroutine DETERMINE_OCEAN_POINTS_(pctwater)
real(kind=FMS_TOP_KIND_), intent(inout) :: pctwater(:,:)
logical :: ocean(size(pctwater,1),size(pctwater,2))
integer :: i, j, m, n, im, ip, jm, jp, new
integer, parameter :: lkind = FMS_TOP_KIND_
real(kind=FMS_TOP_KIND_) :: ocean_pct_crit = real(0.500_lkind)


! resolution of the grid
m = size(pctwater,1)
n = size(pctwater,2)

! the 1/6 degree navy percent water data set
! designates ocean grid boxes as 100 percent water
! all other grid boxes have <= 99 percent water

! set a mask for ocean grid boxes
ocean = (pctwater > 0.999_lkind)
new = count(ocean)

! set land grid boxes that have sufficient amount of water
! to ocean grid boxes when they are adjacent to ocean points
! iterate until there are no new ocean points
do
    if (new == 0) exit
        new = 0

       do j = 1, n
         do i = 1, m
            if (.not.ocean(i,j) .and. pctwater(i,j) > ocean_pct_crit) then
                 im = i-1; ip = i+1; jm = j-1; jp = j+1
                if (im == 0)   im = m
                if (ip == m+1) ip = 1
                if (jm == 0)   jm = 1
                if (jp == n+1) jp = n
    ! check the 8 grid boxes that surround this grid box
                if (ocean(im,j ) .or. ocean(ip,j ) .or. ocean(i ,jm) .or. ocean(i ,jp) .or. &
                    ocean(im,jm) .or. ocean(ip,jm) .or. ocean(ip,jp) .or. ocean(im,jp)) then
                    ocean(i,j) = .true.
                    new = new + 1
                endif
            endif
          enddo
        enddo
      !print *, 'new=',new

enddo

! final step is to elimate water percentage if land
where (.not.ocean) pctwater = 0.0_lkind

end subroutine DETERMINE_OCEAN_POINTS_

!> @}
! close documentation grouping
