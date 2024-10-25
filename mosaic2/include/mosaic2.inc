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

!> @file

!> @brief Get exchange grid information from mosaic xgrid file.
!> Example usage:
!!
!!           call get_mosaic_xgrid(fileobj, nxgrid, i1, j1, i2, j2, area)
!!
  subroutine GET_MOSAIC_XGRID_(fileobj, i1, j1, i2, j2, area, ibegin, iend)
    type(FmsNetcdfFile_t), intent(in) :: fileobj !> The file that contains exchange grid information.
    integer,       intent(inout) :: i1(:), j1(:), i2(:), j2(:) !> i and j indices for grids 1 and 2
    real(kind=FMS_MOS_KIND_),  intent(inout) :: area(:) !> area of the exchange grid. The area is scaled to
                                                        !! represent unit earth area
    integer, optional, intent(in) :: ibegin, iend

    integer                 :: start(4), nread(4), istart
    real(kind=FMS_MOS_KIND_),  dimension(2, size(i1(:))) :: tile1_cell, tile2_cell
    integer                 :: nxgrid, n
    real(kind=r8_kind) :: garea
    real(kind=r8_kind) :: get_global_area

    garea = get_global_area() !< get_global_area returns a r8_kind

    ! When start and nread present, make sure nread(1) is the same as the size of the data
    if(present(ibegin) .and. present(iend)) then
       istart = ibegin
       nxgrid = iend - ibegin + 1
       if(nxgrid .NE. size(i1(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(i1(:))")
       if(nxgrid .NE. size(j1(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(j1(:))")
       if(nxgrid .NE. size(i2(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(i2(:))")
       if(nxgrid .NE. size(j2(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(j2(:))")
       if(nxgrid .NE. size(area(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(area(:))")
    else
       istart = 1
       nxgrid = size(i1(:))
    endif

    start  = 1; nread = 1
    start(1) = istart; nread(1) = nxgrid

    call read_data(fileobj, 'xgrid_area', area, corner=start, edge_lengths=nread)

    start = 1; nread = 1
    nread(1) = 2
    start(2) = istart; nread(2) = nxgrid

    call read_data(fileobj, 'tile1_cell', tile1_cell, corner=start, edge_lengths=nread)
    call read_data(fileobj, 'tile2_cell', tile2_cell, corner=start, edge_lengths=nread)

     do n = 1, nxgrid
       i1(n) = int(tile1_cell(1,n))
       j1(n) = int(tile1_cell(2,n))
       i2(n) = int(tile2_cell(1,n))
       j2(n) = int(tile2_cell(2,n))
       area(n) = real( real(area(n),r8_kind)/garea, FMS_MOS_KIND_ )
    end do

    return

  end subroutine GET_MOSAIC_XGRID_
  !###############################################################################
  !> @brief Calculate grid cell area.
  !> Calculate the grid cell area. The purpose of this routine is to make
  !! sure the consistency between model grid area and exchange grid area.
  !> @param lon geographical longitude of grid cell vertices.
  !> @param lat geographical latitude of grid cell vertices.
  !> @param[inout] area grid cell area.
  !> <br>Example usage:
  !!            call calc_mosaic_grid_area(lon, lat, area)
  subroutine CALC_MOSAIC_GRID_AREA_(lon, lat, area)
     real(kind=FMS_MOS_KIND_), dimension(:,:), intent(in)  :: lon
     real(kind=FMS_MOS_KIND_), dimension(:,:), intent(in)  :: lat
     real(kind=FMS_MOS_KIND_), dimension(:,:), intent(out) :: area
     integer :: nlon, nlat

     real(r8_kind) :: area_r8(size(area,1),size(area,2))

     nlon = size(area,1)
     nlat = size(area,2)
     ! make sure size of lon, lat and area are consitency
     if( size(lon,1) .NE. nlon+1 .OR. size(lat,1) .NE. nlon+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,1) and size(lat,1) should equal to size(area,1)+1")
     if( size(lon,2) .NE. nlat+1 .OR. size(lat,2) .NE. nlat+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,2) and size(lat,2) should equal to size(area,2)+1")

     ! get_grid_area only accepts double precision data
     call get_grid_area( nlon, nlat, real(lon,r8_kind), real(lat,r8_kind), area_r8)

     area=real(area_r8,FMS_MOS_KIND_)

  end subroutine CALC_MOSAIC_GRID_AREA_
  !###############################################################################
  !> @brief Calculate grid cell area using great cirlce algorithm
  !> Calculate the grid cell area. The purpose of this routine is to make
  !! sure the consistency between model grid area and exchange grid area.
  !> @param lon geographical longitude of grid cell vertices.
  !> @param lat geographical latitude of grid cell vertices.
  !> @param[inout] area grid cell area.
  !> <br>Example usage:
  !!            call calc_mosaic_grid_great_circle_area(lon, lat, area)
  subroutine CALC_MOSAIC_GRID_GREAT_CIRCLE_AREA_(lon, lat, area)
     real(kind=FMS_MOS_KIND_), dimension(:,:), intent(in)    :: lon
     real(kind=FMS_MOS_KIND_), dimension(:,:), intent(in)    :: lat
     real(kind=FMS_MOS_KIND_), dimension(:,:), intent(inout) :: area
     integer :: nlon, nlat

     real(r8_kind) :: area_r8(size(area,1),size(area,2))

     nlon = size(area,1)
     nlat = size(area,2)
     ! make sure size of lon, lat and area are consitency
     if( size(lon,1) .NE. nlon+1 .OR. size(lat,1) .NE. nlon+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,1) and size(lat,1) should equal to size(area,1)+1")
     if( size(lon,2) .NE. nlat+1 .OR. size(lat,2) .NE. nlat+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,2) and size(lat,2) should equal to size(area,2)+1")

     ! get_grid_great_circle_area only accepts r8_kind arguments
     call get_grid_great_circle_area( nlon, nlat, real(lon,r8_kind), real(lat,r8_kind), area_r8)

     area=real(area_r8, FMS_MOS_KIND_)

   end subroutine CALC_MOSAIC_GRID_GREAT_CIRCLE_AREA_
  !#####################################################################
  !> This function check if a point (lon1,lat1) is inside a polygon (lon2(:), lat2(:))
  !! lon1, lat1, lon2, lat2 are in radians.
  function IS_INSIDE_POLYGON_(lon1, lat1, lon2, lat2 )
     real(kind=FMS_MOS_KIND_), intent(in) :: lon1, lat1
     real(kind=FMS_MOS_KIND_), intent(in) :: lon2(:), lat2(:)
     logical :: IS_INSIDE_POLYGON_
     integer :: npts, isinside
     integer :: inside_a_polygon

     npts = size(lon2(:))

     !> inside_a_polygon function only accepts r8_kind real variables

     isinside = inside_a_polygon(real(lon1,r8_kind), real(lat1,r8_kind), npts, real(lon2,r8_kind), real(lat2,r8_kind))
     if(isinside == 1) then
        IS_INSIDE_POLYGON_ = .TRUE.
     else
        IS_INSIDE_POLYGON_ = .FALSE.
     endif

     return

  end function IS_INSIDE_POLYGON_
!> @}
