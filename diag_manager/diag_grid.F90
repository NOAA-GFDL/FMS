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

!> @defgroup diag_grid_mod diag_grid_mod
!> @ingroup diag_manager
!> @brief diag_grid_mod is a set of procedures to work with the
!!   model's global grid to allow regional output.
!!
!> @author Seth Underwood seth.underwood@noaa.gov
!!
!! <TT>diag_grid_mod</TT> contains useful utilities for dealing
!!   with, mostly, regional output for grids other than the standard
!!   lat/lon grid.  This module contains three public procedures <TT>
!!   diag_grid_init</TT>, which is shared globably in the <TT>
!!   diag_manager_mod</TT>, <TT>diag_grid_end</TT> which will free
!!   up memory used during the register field calls, and
!!   <TT>get_local_indexes</TT>.  The <TT>send_global_grid</TT>
!!   procedure is called by the model that creates the global grid.
!!   <TT>send_global_grid</TT> needs to be called before any fields
!!   are registered that will output only regions.  <TT>get_local_indexes</TT>
!!   is to be called by the <TT>diag_manager_mod</TT> to discover the
!!   global indexes defining a subregion on the tile.

!> @file
!> @brief File for @ref diag_grid_mod

MODULE diag_grid_mod
use platform_mod

  ! <INFO>
  !   <FUTURE>
  !     Multi-tile regional output in the cubed sphere.
  !   </FUTURE>
  !   <FUTURE>
  !     Single grid in the tri-polar grid.
  !   </FUTURE>
  !   <FUTURE>
  !     Multi-tile regional output in the tri-polar grid.
  !   </FUTURE>
  !   <FUTURE>
  !     Regional output using array masking.  This should allow
  !     regional output to work on any current or future grid.
  !   </FUTURE>
  ! </INFO>

  USE constants_mod, ONLY: DEG_TO_RAD, RAD_TO_DEG, RADIUS
  USE fms_mod, ONLY: write_version_number, error_mesg, WARNING, FATAL,&
       & mpp_pe
  USE mpp_mod, ONLY: mpp_root_pe, mpp_npes, mpp_max, mpp_min
  USE mpp_domains_mod, ONLY: domain2d, mpp_get_tile_id,&
       & mpp_get_ntile_count, mpp_get_compute_domains

  IMPLICIT NONE

  ! Include variable "version" to be written to log file.
#include<file_version.h>

  !> @brief Private type to hold the model's global grid data, and other grid information for use
  !! in this module.
  !> @ingroup diag_grid_mod
  type, private :: diag_global_grid_type
     REAL, allocatable, DIMENSION(:,:) :: glo_lat !< The latitude values on the global grid.
     REAL, allocatable, DIMENSION(:,:) :: glo_lon !< The longitude values on the global grid.
     REAL, allocatable, DIMENSION(:,:) :: aglo_lat !< The latitude values on the global a-grid.  Here we expect isc-1:iec+1 and
                                                   !! jsc=1:jec+1 to be passed in.
     REAL, allocatable, DIMENSION(:,:) :: aglo_lon !< The longitude values on the global a-grid.  Here we expec isc-1:iec+j and
                                                   !! jsc-1:jec+1 to be passed in.
     INTEGER :: myXbegin !< The starting index of the compute domain on the current PE.
     INTEGER :: myYbegin !< The starting index of the compute domain on the current PE.
     INTEGER :: dimI !< The dimension of the global grid in the 'i' / longitudal direction.
     INTEGER :: dimJ !< The dimension of the global grid in the 'j' / latitudal direction.
     INTEGER :: adimI !< The dimension of the global a-grid in the 'i' / longitudal direction.  Again,
                      !! the expected dimension for diag_grid_mod is isc-1:iec+1.
     INTEGER :: adimJ !< The dimension of the global a-grid in the 'j' / latitudal direction.  Again,
                      !! the expected dimension for diag_grid_mod is jsc-1:jec+1.
     INTEGER :: tile_number !< The tile the <TT>glo_lat</TT> and <TT>glo_lon</TT> define.
     INTEGER :: ntiles !< The number of tiles.
     INTEGER :: peStart !< The starting PE number for the current tile.
     INTEGER :: peEnd !< The ending PE number for the current tile.
     CHARACTER(len=128) :: grid_type !< The global grid type.
  END TYPE diag_global_grid_type

  !> @brief Private type to hold the corresponding (x,y,z) location for a (lat,lon)
  !! location.
  !> @ingroup diag_grid_mod
  type, private :: point
     REAL :: x !< The x value of the (x,y,z) coordinates.
     REAL :: y !< The y value of the (x,y,z) coordinates.
     REAL :: z !< The z value of the (x,y,z) coordinates.
  END TYPE point

!> @addtogroup diag_grid_mod
!> @{

  TYPE(diag_global_grid_type) :: diag_global_grid !< Variable to hold the global grid data

  LOGICAL :: diag_grid_initialized = .FALSE. !< Indicates if the diag_grid_mod has been initialized.

  PRIVATE
  PUBLIC :: diag_grid_init, diag_grid_end, get_local_indexes,  &
            get_local_indexes2

CONTAINS

  !> @brief Send the global grid to the <TT>diag_manager_mod</TT> for
  !!   regional output.
  !!
  !> In order for the diag_manager to do regional output for grids
  !! other than the standard lat/lon grid, the <TT>
  !! diag_manager_mod</TT> needs to know the the latitude and
  !! longitude values for the entire global grid.  This procedure
  !! is the mechanism the models will use to share their grid with
  !! the diagnostic manager.
  !!
  !! This procedure needs to be called after the grid is created,
  !! and before the first call to register the fields.
  SUBROUTINE diag_grid_init(domain, glo_lat, glo_lon, aglo_lat, aglo_lon)
    TYPE(domain2d), INTENT(in) :: domain !< The domain to which the grid data corresponds.
    REAL, INTENT(in), DIMENSION(:,:) :: glo_lat !< The latitude information for the grid tile.
    REAL, INTENT(in), DIMENSION(:,:) :: glo_lon !< The longitude information for the grid tile.
    REAL, INTENT(in), DIMENSION(:,:) :: aglo_lat !< The latitude information for the a-grid tile.
    REAL, INTENT(in), DIMENSION(:,:) :: aglo_lon !< The longitude information for the a-grid tile.

    INTEGER, DIMENSION(1) :: tile
    INTEGER :: ntiles
    INTEGER :: stat
    INTEGER :: i_dim, j_dim
    INTEGER :: ai_dim, aj_dim
    INTEGER, DIMENSION(2) :: latDim, lonDim
    INTEGER, DIMENSION(2) :: alatDim, alonDim
    INTEGER :: myPe, npes, npesPerTile
    INTEGER, ALLOCATABLE, DIMENSION(:) :: xbegin, xend, ybegin, yend

    ! Write the file version to the logfile
    CALL write_version_number("DIAG_GRID_MOD", version)

    ! Verify all allocatable / pointers for diag_global_grid hare not
    ! allocated / associated.
    IF ( ALLOCATED(xbegin) ) DEALLOCATE(xbegin)
    IF ( ALLOCATED(ybegin) ) DEALLOCATE(ybegin)
    IF ( ALLOCATED(xend) ) DEALLOCATE(xend)
    IF ( ALLOCATED(yend) ) DEALLOCATE(yend)

    ! What is my PE
    myPe = mpp_pe() -mpp_root_pe() + 1

    ! Get the domain/pe layout, and allocate the [xy]begin|end arrays/pointers
    npes = mpp_npes()
    ALLOCATE(xbegin(npes), &
         &   ybegin(npes), &
         &   xend(npes), &
         &   yend(npes), STAT=stat)
    IF ( stat .NE. 0 ) THEN
       CALL error_mesg('diag_grid_mod::diag_grid_init',&
            &'Could not allocate memory for the compute grid indices&
            &.', FATAL)
    END IF

    ! Get tile information
    ntiles = mpp_get_ntile_count(domain)
    tile = mpp_get_tile_id(domain)

    ! Number of PEs per tile
    npesPerTile = npes / ntiles
    diag_global_grid%peEnd = npesPerTile * tile(1)
    diag_global_grid%peStart = diag_global_grid%peEnd - npesPerTile + 1

    ! Get the compute domains
    CALL mpp_get_compute_domains(domain,&
         & XBEGIN=xbegin, XEND=xend,&
         & YBEGIN=ybegin, YEND=yend)

    ! Module initialized
    diag_grid_initialized = .TRUE.

    ! Get the size of the grids
    latDim = SHAPE(glo_lat)
    lonDim = SHAPE(glo_lon)
    IF (  (latDim(1) == lonDim(1)) .AND.&
         &(latDim(2) == lonDim(2)) ) THEN
       i_dim = latDim(1)
       j_dim = latDim(2)
    ELSE
       CALL error_mesg('diag_grid_mod::diag_grid_init',&
            &'glo_lat and glo_lon must be the same shape.', FATAL)
    END IF

    ! Same thing for the a-grid
    alatDim = SHAPE(aglo_lat)
    alonDim = SHAPE(aglo_lon)
    IF (  (alatDim(1) == alonDim(1)) .AND. &
         &(alatDim(2) == alonDim(2)) ) THEN
       IF ( tile(1) == 4 .OR. tile(1) == 5 ) THEN
          ! These tiles need to be transposed.
          ai_dim = alatDim(2)
          aj_dim = alatDim(1)
       ELSE
          ai_dim = alatDim(1)
          aj_dim = alatDim(2)
       END IF
    ELSE
       CALL error_mesg('diag_grid_mod::diag_grid_init',&
            & "a-grid's glo_lat and glo_lon must be the same shape.", FATAL)
    END IF

    ! Allocate the grid arrays
    IF (   allocated(diag_global_grid%glo_lat) .OR.&
         & allocated(diag_global_grid%glo_lon) ) THEN
       IF ( mpp_pe() == mpp_root_pe() ) &
            & CALL error_mesg('diag_grid_mod::diag_grid_init',&
            &'The global grid has already been initialized', WARNING)
    ELSE
       ALLOCATE(diag_global_grid%glo_lat(i_dim,j_dim),&
            &   diag_global_grid%glo_lon(i_dim,j_dim), STAT=stat)
       IF ( stat .NE. 0 ) THEN
          CALL error_mesg('diag_grid_mod::diag_grid_init',&
               &'Could not allocate memory for the global grid.', FATAL)
       END IF
    END IF

    ! Same thing for the a-grid
    IF (   allocated(diag_global_grid%aglo_lat) .OR.&
         & allocated(diag_global_grid%aglo_lon) ) THEN
       IF ( mpp_pe() == mpp_root_pe() ) &
            & CALL error_mesg('diag_grid_mod::diag_grid_init',&
            &'The global a-grid has already been initialized', WARNING)
    ELSE
       ALLOCATE(diag_global_grid%aglo_lat(0:ai_dim-1,0:aj_dim-1),&
            &   diag_global_grid%aglo_lon(0:ai_dim-1,0:aj_dim-1), STAT=stat)
       IF ( stat .NE. 0 ) THEN
          CALL error_mesg('diag_global_mod::diag_grid_init',&
               &'Could not allocate memory for the global a-grid', FATAL)
       END IF
    END IF

    ! Set the values for diag_global_grid

    ! If we are on tile 4 or 5, we need to transpose the grid to get
    ! this to work.
    IF ( tile(1) == 4 .OR. tile(1) == 5 ) THEN
       diag_global_grid%aglo_lat = TRANSPOSE(aglo_lat)
       diag_global_grid%aglo_lon = TRANSPOSE(aglo_lon)
    ELSE
       diag_global_grid%aglo_lat = aglo_lat
       diag_global_grid%aglo_lon = aglo_lon
    END IF
    diag_global_grid%glo_lat = glo_lat
    diag_global_grid%glo_lon = glo_lon
    diag_global_grid%dimI = i_dim
    diag_global_grid%dimJ = j_dim
    diag_global_grid%adimI = ai_dim
    diag_global_grid%adimJ = aj_dim
    !--- For the nested model, the nested region only has 1 tile ( ntiles = 1) but
    !--- the tile_id is 7 for the nested region. In the routine get_local_indexes,
    !--- local variables ijMin and ijMax have dimesnion (ntiles) and will access
    !--- ijMin(diag_global_grid%tile_number,:). For the nested region, ntiles = 1 and
    !--- diag_global_grid%tile_number = 7 will cause out of bounds. So need to
    !--- set diag_global_grid%tile_number = 1 when ntiles = 1 for the nested model.
    if(ntiles == 1) then
       diag_global_grid%tile_number = 1
    else
       diag_global_grid%tile_number = tile(1)
    endif
    diag_global_grid%ntiles = ntiles
    diag_global_grid%myXbegin = xbegin(myPe)
    diag_global_grid%myYbegin = ybegin(myPe)

    ! Unallocate arrays used here
    DEALLOCATE(xbegin)
    DEALLOCATE(ybegin)
    DEALLOCATE(xend)
    DEALLOCATE(yend)
  END SUBROUTINE diag_grid_init

  !> @brief Unallocate the diag_global_grid variable.
  !!
  !> The <TT>diag_global_grid</TT> variable is only needed during
  !! the register field calls, and then only if there are fields
  !! requestion regional output.  Once all the register fields
  !! calls are complete (before the first <TT>send_data</TT> call
  !! this procedure can be called to free up memory.
  SUBROUTINE diag_grid_end()

    IF ( diag_grid_initialized ) THEN
       ! De-allocate grid
       IF ( allocated(diag_global_grid%glo_lat) ) THEN
          DEALLOCATE(diag_global_grid%glo_lat)
       ELSE IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_grid_mod::diag_grid_end',&
               &'diag_global_grid%glo_lat was not allocated.', WARNING)
       END IF

       IF ( allocated(diag_global_grid%glo_lon) ) THEN
          DEALLOCATE(diag_global_grid%glo_lon)
       ELSE IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_grid_mod::diag_grid_end',&
               &'diag_global_grid%glo_lon was not allocated.', WARNING)
       END IF
       ! De-allocate a-grid
       IF ( allocated(diag_global_grid%aglo_lat) ) THEN
          DEALLOCATE(diag_global_grid%aglo_lat)
       ELSE IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_grid_mod::diag_grid_end',&
               &'diag_global_grid%aglo_lat was not allocated.', WARNING)
       END IF

       IF ( allocated(diag_global_grid%aglo_lon) ) THEN
          DEALLOCATE(diag_global_grid%aglo_lon)
       ELSE IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_grid_mod::diag_grid_end',&
               &'diag_global_grid%aglo_lon was not allocated.', WARNING)
       END IF

       diag_grid_initialized = .FALSE.
    END IF
  END SUBROUTINE diag_grid_end

  !> @brief Find the local start and local end indexes on the local PE
  !!   for regional output.
  !!
  !> Given a defined region, find the local indexes on the local
  !!   PE surrounding the region.
  SUBROUTINE get_local_indexes(latStart, latEnd, lonStart, lonEnd,&
       & istart, iend, jstart, jend)
    REAL, INTENT(in) :: latStart !< lat start angles
    REAL, INTENT(in) :: lonStart !< lon start angles
    REAL, INTENT(in) :: latEnd !< lat end angles
    REAL, INTENT(in) :: lonEnd !< lon end angles
    INTEGER, INTENT(out) :: istart !< i start indexes
    INTEGER, INTENT(out) :: jstart !< j start indexes
    INTEGER, INTENT(out) :: iend !< i end indexes
    INTEGER, INTENT(out) :: jend !< j end indexes

    REAL, ALLOCATABLE, DIMENSION(:,:) :: delta_lat, delta_lon, grid_lon

    REAL, DIMENSION(4) :: dists_lon, dists_lat
    REAL :: lonEndAdj, my_lonStart, my_lonEnd

    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ijMin, ijMax
    INTEGER :: myTile, ntiles, i, j, k, dimI, dimJ, istat
    INTEGER :: count

    LOGICAL :: onMyPe

    !For cfsite potential fix.
    INTEGER :: minI
    INTEGER :: minJ
    REAL :: minimum_distance
    REAL :: global_min_distance
    INTEGER :: rank_buf

    IF ( .NOT. diag_grid_initialized )&
         & CALL error_mesg('diag_grid_mod::get_local_indexes',&
         &'Module not initialized, first initialze module with a call &
         &to diag_grid_init', FATAL)

    ! Make adjustment for negative longitude values
    if ( lonStart < 0. ) then
       my_lonStart = lonStart + 360.
    else
       my_lonStart = lonStart
    end if
    if ( lonEnd < 0. ) then
       my_lonEnd = lonEnd + 360.
    else
       my_lonEnd = lonEnd
    end if

    IF (latStart .EQ. latEnd .AND. my_lonStart .EQ. my_lonEnd) THEN

        !For a single point, use the a-grid longitude and latitude
        !values.

        myTile = diag_global_grid%tile_number
        ntiles = diag_global_grid%ntiles

        allocate(ijMin(ntiles,2))
        ijMin = 0

        !Find the i,j indices of the a-grid point nearest to the
        !my_lonStart,latStart point.
        CALL find_nearest_agrid_index(latStart, &
                                      my_lonStart, &
                                      minI, &
                                      minJ, &
                                      minimum_distance)

        !Find the minimum distance across all ranks.
        global_min_distance = minimum_distance
        CALL mpp_min(global_min_distance)

        !In the case of a tie (i.e. two ranks with exactly the same
        !minimum distance), use the i,j values from the larger rank id.
        IF (global_min_distance .EQ. minimum_distance) THEN
            rank_buf = mpp_pe()
        ELSE
            rank_buf = -1
        ENDIF
        CALL mpp_max(rank_buf)

        !Sanity check.
        IF (rank_buf .EQ. -1) THEN
            CALL error_mesg("get_local_indexes", &
                            "No rank has minimum distance.", &
                            FATAL)
        ENDIF

        IF (rank_buf .EQ. mpp_pe()) THEN
            ijMin(mytile,1) = minI + diag_global_grid%myXbegin - 1
            ijMin(mytile,2) = minJ + diag_global_grid%myYbegin - 1
        ENDIF

        DO i = 1,ntiles
            CALL mpp_max(ijMin(i,1))
            CALL mpp_max(ijMin(i,2))
        ENDDO

        istart = ijMin(mytile,1)
        jstart = ijMin(mytile,2)
        iend = istart
        jend = jstart

        DEALLOCATE(ijMin)
    ELSE

        myTile = diag_global_grid%tile_number
        ntiles = diag_global_grid%ntiles

        ! Arrays to home min/max for each tile
        ALLOCATE(ijMin(ntiles,2), STAT=istat)
        IF ( istat .NE. 0 )&
             & CALL error_mesg('diag_grid_mod::get_local_indexes',&
             &'Cannot allocate ijMin index array', FATAL)
        ALLOCATE(ijMax(ntiles,2), STAT=istat)
        IF ( istat .NE. 0 )&
             & CALL error_mesg('diag_grid_mod::get_local_indexes',&
             &'Cannot allocate ijMax index array', FATAL)
        ijMin = 0
        ijMax = 0

        ! There will be four points to define a region, find all four.
        ! Need to call the correct function depending on if the tile is a
        ! pole tile or not.
       dimI = diag_global_grid%dimI
       dimJ = diag_global_grid%dimJ

       ! Build the delta array
       ALLOCATE(delta_lat(dimI,dimJ), STAT=istat)
       IF ( istat .NE. 0 )&
            & CALL error_mesg('diag_grid_mod::get_local_indexes',&
            &'Cannot allocate latitude delta array', FATAL)
       ALLOCATE(delta_lon(dimI,dimJ), STAT=istat)
       IF ( istat .NE. 0 )&
            & CALL error_mesg('diag_grid_mod::get_local_indexes',&
            &'Cannot allocate longitude delta array', FATAL)
       DO j=1, dimJ
          DO i=1, dimI
             count = 0
             dists_lon = 0.
             dists_lat = 0.
             IF ( i < dimI ) THEN
                dists_lon(1) = ABS(diag_global_grid%glo_lon(i+1,j) - diag_global_grid%glo_lon(i,j))
                dists_lat(1) = ABS(diag_global_grid%glo_lat(i+1,j) - diag_global_grid%glo_lat(i,j))
                count = count+1
             END IF
             IF ( j < dimJ ) THEN
                dists_lon(2) = ABS(diag_global_grid%glo_lon(i,j+1) - diag_global_grid%glo_lon(i,j))
                dists_lat(2) = ABS(diag_global_grid%glo_lat(i,j+1) - diag_global_grid%glo_lat(i,j))
                count = count+1
             END IF
             IF ( i > 1 ) THEN
                dists_lon(3) = ABS(diag_global_grid%glo_lon(i,j) - diag_global_grid%glo_lon(i-1,j))
                dists_lat(3) = ABS(diag_global_grid%glo_lat(i,j) - diag_global_grid%glo_lat(i-1,j))
                count = count+1
             END IF
             IF ( j > 1 ) THEN
                dists_lon(4) = ABS(diag_global_grid%glo_lon(i,j) - diag_global_grid%glo_lon(i,j-1))
                dists_lat(4) = ABS(diag_global_grid%glo_lat(i,j) - diag_global_grid%glo_lat(i,j-1))
                count = count+1
             END IF

             ! Fix wrap around problem
             DO k=1, 4
                IF ( dists_lon(k) > 180.0 ) THEN
                   dists_lon(k) = 360.0 - dists_lon(k)
                END IF
             END DO
             delta_lon(i,j) = SUM(dists_lon)/real(count)
             delta_lat(i,j) = SUM(dists_lat)/real(count)
          END DO
       END DO

       ijMin = HUGE(1)
       ijMax = -HUGE(1)

       ! Adjusted longitude array
       ALLOCATE(grid_lon(dimI,dimJ), STAT=istat)
       IF ( istat .NE. 0 )&
            & CALL error_mesg('diag_grid_mod::get_local_indexes',&
            &'Cannot allocate temporary longitude array', FATAL)
       grid_lon = diag_global_grid%glo_lon

       ! Make adjustments where required
       IF ( my_lonStart > my_lonEnd ) THEN
          WHERE ( grid_lon < my_lonStart )
             grid_lon = grid_lon + 360.0
          END WHERE
          lonEndAdj = my_lonEnd + 360.0
       ELSE
          lonEndAdj = my_lonEnd
       END IF

       DO j=1, dimJ-1
          DO i=1, dimI-1
             onMyPe = .false.
             IF ( latStart-delta_lat(i,j) <= diag_global_grid%glo_lat(i,j) .AND.&
                  & diag_global_grid%glo_lat(i,j) < latEnd+delta_lat(i,j) ) THEN
                ! Short-cut for the poles
                IF ( (ABS(latStart)-delta_lat(i,j) <= 90.0 .AND.&
                     & 90.0 <= ABS(latEnd)+delta_lat(i,j)) .AND.&
                     & ABS(diag_global_grid%glo_lat(i,j)) == 90.0 ) THEN
                   onMyPe = .TRUE.
                ELSE IF ( (my_lonStart-delta_lon(i,j) <= grid_lon(i,j) .AND.&
                     & grid_lon(i,j) < lonEndAdj+delta_lon(i,j)) ) THEN
                   onMyPe = .TRUE.
                ELSE
                   onMyPe = .FALSE.
                END IF
                IF ( onMyPe ) THEN
                   ijMin(myTile,1) = MIN(ijMin(myTile,1),i + diag_global_grid%myXbegin - 1)
                   ijMax(myTile,1) = MAX(ijMax(myTile,1),i + diag_global_grid%myXbegin - 1)
                   ijMin(myTile,2) = MIN(ijMin(myTile,2),j + diag_global_grid%myYbegin - 1)
                   ijMax(myTile,2) = MAX(ijMax(myTile,2),j + diag_global_grid%myYbegin - 1)
                END IF
             END IF
          END DO
       END DO
       DEALLOCATE(delta_lon)
       DEALLOCATE(delta_lat)
       DEALLOCATE(grid_lon)

       ! Global min/max reduce
       DO i=1, ntiles
          CALL mpp_min(ijMin(i,1))
          CALL mpp_max(ijMax(i,1))
          CALL mpp_min(ijMin(i,2))
          CALL mpp_max(ijMax(i,2))
       END DO

       IF ( ijMin(myTile,1) == HUGE(1) .OR. ijMax(myTile,1) == -HUGE(1) ) THEN
          ijMin(myTile,1) = 0
          ijMax(myTile,1) = 0
       END IF
       IF ( ijMin(myTile,2) == HUGE(1) .OR. ijMax(myTile,2) == -HUGE(1) ) THEN
          ijMin(myTile,2) = 0
          ijMax(myTile,2) = 0
       END IF

       istart = ijMin(myTile,1)
       jstart = ijMin(myTile,2)
       iend = ijMax(myTile,1)
       jend = ijMax(myTile,2)

       DEALLOCATE(ijMin)
       DEALLOCATE(ijMax)
    END IF

  END SUBROUTINE get_local_indexes

  !> @brief Find the indices of the nearest grid point of the a-grid to the
  !!   specified (lon,lat) location on the local PE. if desired point not
  !!   within domain of local PE, return (0,0) as the indices.
  SUBROUTINE get_local_indexes2(lat, lon, iindex, jindex)
    REAL, INTENT(in) :: lat !< lat location
    REAL, INTENT(in) :: lon !< lon location
    INTEGER, INTENT(out) :: iindex !< i indexes
    INTEGER, INTENT(out) :: jindex !< j indexes

    INTEGER  :: indexes(2)

    IF ( .NOT. diag_grid_initialized )&
         & CALL error_mesg('diag_grid_mod::get_local_indexes2',&
         &'Module not initialized, first initialze module with a call &
         &to diag_grid_init', FATAL)

    indexes = 0

    IF ( MOD(diag_global_grid%tile_number,3) == 0 ) THEN
       IF ( lat > 30.0 .AND. diag_global_grid%tile_number == 3 ) THEN
          indexes(:) = find_pole_index_agrid(lat,lon)
       ELSE IF ( lat < -30.0 .AND. diag_global_grid%tile_number == 6 ) THEN
          indexes(:) = find_pole_index_agrid(lat,lon)
       ENDIF
    ELSE
       indexes(:) = find_equator_index_agrid(lat,lon)
    END IF

    iindex = indexes(1)
    jindex = indexes(2)
    IF (iindex ==  diag_global_grid%adimI -1 .OR.&
        jindex ==  diag_global_grid%adimJ -1 ) THEN
      iindex = 0
      jindex = 0
    ENDIF

  END SUBROUTINE get_local_indexes2

  !> @fn pure elemental real rad2deg(real angle)
  !> @brief Convert an angle in radian to degrees.
  !!
  !> Given a scalar, or an array of angles in radians this
  !! function will return a scalar or array (of the same
  !! dimension) of angles in degrees.
  !! @return Scalar or array (depending on the size of angle) of angles in
  !! degrees.
  PURE ELEMENTAL REAL FUNCTION rad2deg(angle)
    REAL, INTENT(in) :: angle !< Scalar or array of angles in radians.

    rad2deg = RAD_TO_DEG * angle
  END FUNCTION rad2deg

  !> @brief Convert an angle in degrees to radians.
  !!
  !> Given a scalar, or an array of angles in degrees this
  !!   function will return a scalar or array (of the same
  !!   dimension) of angles in radians.
  !! @return Scalar or array (depending on the size of angle) of angles in
  !!   radians.
  PURE ELEMENTAL REAL FUNCTION deg2rad(angle)
    REAL, INTENT(in) :: angle !< Scalar or array of angles in degrees.

    deg2rad = DEG_TO_RAD * angle
  END FUNCTION deg2rad

  !> @brief Return the closest index (i,j) to the given (lat,lon) point.
  !!
  !> This function searches a pole a-grid tile looking for the grid point
  !!   closest to the give (lat, lon) location, and returns the i
  !!   and j indexes of the point.
  !! @return The (i, j) location of the closest grid to the given (lat,
  !!   lon) location.
  PURE FUNCTION find_pole_index_agrid(lat, lon)
    INTEGER, DIMENSION(2) :: find_pole_index_agrid !< The (i, j) location of the closest grid to the given (lat,
                                                   !! lon) location.
    REAL, INTENT(in) :: lat !< Latitude location
    REAL, INTENT(in) :: lon !< Longitude location

    INTEGER :: indxI !< Indexes to be returned.
    INTEGER :: indxJ !< Indexes to be returned.
    INTEGER :: dimI !< Size of the grid dimensions
    INTEGER :: dimJ !< Size of the grid dimensions
    INTEGER :: i !< Count indexes
    INTEGER :: j !< Count indexes
    INTEGER :: nearestCorner !< index of the nearest corner.
    INTEGER, DIMENSION(4,2) :: ijArray !< indexes of the cornerPts and pntDistances arrays
    REAL :: llat, llon
    REAL :: maxCtrDist !< maximum distance to center of grid
    REAL, DIMENSION(4) :: pntDistances !< distance from origPt to corner
    TYPE(point) :: origPt !< Original point
    REAL, DIMENSION(4,2) :: cornerPts !< Corner points using (lat,lon)
    TYPE(point), DIMENSION(9) :: points !< xyz of 8 nearest neighbors
    REAL, DIMENSION(9) :: distSqrd !< distance between origPt and points(:)

    ! Set the inital fail values for indxI and indxJ
    indxI = 0
    indxJ = 0

    dimI = diag_global_grid%adimI
    dimJ = diag_global_grid%adimJ

    ! Since the poles have an non-unique longitude value, make a small correction if looking for one of the poles.
    IF ( lat == 90.0 ) THEN
       llat = lat - .1
    ELSE IF ( lat == -90.0 ) THEN
       llat = lat + .1
    ELSE
       llat = lat
    END IF
    llon = lon

    origPt = latlon2xyz(llat,llon)

    iLoop: DO i=0, dimI-2
       jLoop: DO j = 0, dimJ-2
          cornerPts = RESHAPE( (/ diag_global_grid%aglo_lat(i,  j),  diag_global_grid%aglo_lon(i,  j),&
               &                  diag_global_grid%aglo_lat(i+1,j+1),diag_global_grid%aglo_lon(i+1,j+1),&
               &                  diag_global_grid%aglo_lat(i+1,j),  diag_global_grid%aglo_lon(i+1,j),&
               &                  diag_global_grid%aglo_lat(i,  j+1),diag_global_grid%aglo_lon(i,  j+1) /),&
               &               (/ 4, 2 /), ORDER=(/2,1/) )
          ! Find the maximum half distance of the corner points
          maxCtrDist = MAX(gCirDistance(cornerPts(1,1),cornerPts(1,2), cornerPts(2,1),cornerPts(2,2)),&
               &           gCirDistance(cornerPts(3,1),cornerPts(3,2), cornerPts(4,1),cornerPts(4,2)))

          ! Find the distance of the four corner points to the point of interest.
          pntDistances = gCirDistance(cornerPts(:,1),cornerPts(:,2), llat,llon)

          IF ( (MINVAL(pntDistances) <= maxCtrDist) .AND. (i*j.NE.0) ) THEN
             ! Set up the i,j index array
             ijArray = RESHAPE( (/ i, j, i+1, j+1, i+1, j, i, j+1 /), (/ 4, 2 /), ORDER=(/2,1/) )

             ! the nearest point index
             nearestCorner = MINLOC(pntDistances,1)

             indxI = ijArray(nearestCorner,1)
             indxJ = ijArray(nearestCorner,2)

             EXIT iLoop
          END IF
       END DO jLoop
    END DO iLoop


    ! Make sure we have indexes in the correct range
    valid: IF (  (indxI <= 0 .OR. dimI-1 <= indxI) .OR. &
         &       (indxJ <= 0 .OR. dimJ-1 <= indxJ) ) THEN
       indxI = 0
       indxJ = 0
    ELSE ! indxI and indxJ are valid.
       ! Since we are looking for the closest grid point to the
       ! (lat,lon) point, we need to check the surrounding
       ! points.  The indexes for the variable points are as follows
       !
       ! 1---4---7
       ! |   |   |
       ! 2---5---8
       ! |   |   |
       ! 3---6---9

       ! Set the 'default' values for points(:) x,y,z to some large
       ! value.
       DO i=1, 9
          points(i)%x = 1.0e20
          points(i)%y = 1.0e20
          points(i)%z = 1.0e20
       END DO

       ! All the points around the i,j indexes
       points(1) = latlon2xyz(diag_global_grid%aglo_lat(indxI-1,indxJ+1),&
            &                 diag_global_grid%aglo_lon(indxI-1,indxJ+1))
       points(2) = latlon2xyz(diag_global_grid%aglo_lat(indxI-1,indxJ),&
            &                 diag_global_grid%aglo_lon(indxI-1,indxJ))
       points(3) = latlon2xyz(diag_global_grid%aglo_lat(indxI-1,indxJ-1),&
            &                 diag_global_grid%aglo_lon(indxI-1,indxJ-1))
       points(4) = latlon2xyz(diag_global_grid%aglo_lat(indxI,  indxJ+1),&
            &                 diag_global_grid%aglo_lon(indxI,  indxJ+1))
       points(5) = latlon2xyz(diag_global_grid%aglo_lat(indxI,  indxJ),&
            &                 diag_global_grid%aglo_lon(indxI,  indxJ))
       points(6) = latlon2xyz(diag_global_grid%aglo_lat(indxI,  indxJ-1),&
            &                 diag_global_grid%aglo_lon(indxI,  indxJ-1))
       points(7) = latlon2xyz(diag_global_grid%aglo_lat(indxI+1,indxJ+1),&
            &                 diag_global_grid%aglo_lon(indxI+1,indxJ+1))
       points(8) = latlon2xyz(diag_global_grid%aglo_lat(indxI+1,indxJ),&
            &                 diag_global_grid%aglo_lon(indxI+1,indxJ))
       points(9) = latlon2xyz(diag_global_grid%aglo_lat(indxI+1,indxJ-1),&
            &                 diag_global_grid%aglo_lon(indxI+1,indxJ-1))


       ! Calculate the distance squared between the points(:) and the origPt
       distSqrd = distanceSqrd(origPt, points)

       SELECT CASE (MINLOC(distSqrd,1))
       CASE ( 1 )
          indxI = indxI-1
          indxJ = indxJ+1
       CASE ( 2 )
          indxI = indxI-1
          indxJ = indxJ
       CASE ( 3 )
          indxI = indxI-1
          indxJ = indxJ-1
       CASE ( 4 )
          indxI = indxI
          indxJ = indxJ+1
       CASE ( 5 )
          indxI = indxI
          indxJ = indxJ
       CASE ( 6 )
          indxI = indxI
          indxJ = indxJ-1
       CASE ( 7 )
          indxI = indxI+1
          indxJ = indxJ+1
       CASE ( 8 )
          indxI = indxI+1
          indxJ = indxJ
       CASE ( 9 )
          indxI = indxI+1
          indxJ = indxJ-1
       CASE DEFAULT
          indxI = 0
          indxJ = 0
       END SELECT
    END IF valid

    ! Set the return value for the funtion
    find_pole_index_agrid = (/indxI, indxJ/)
  END FUNCTION find_pole_index_agrid

  !> @brief Return the closest index (i,j) to the given (lat,lon) point.
  !!
  !> This function searches a equator grid tile looking for the grid point
  !!   closest to the give (lat, lon) location, and returns the i
  !!   and j indexes of the point.
  !! @return The (i, j) location of the closest grid to the given (lat,
  !!   lon) location.
  PURE FUNCTION find_equator_index_agrid(lat, lon)
    INTEGER, DIMENSION(2) :: find_equator_index_agrid !< The (i, j) location of the closest grid to the given (lat,
                                                      !! lon) location.
    REAL, INTENT(in) :: lat !< Latitude location
    REAL, INTENT(in) :: lon !< Longitude location

    INTEGER :: indxI, indxJ !< Indexes to be returned.
    INTEGER :: indxI_tmp !< Hold the indxI value if on tile 3 or 4
    INTEGER :: dimI, dimJ !< Size of the grid dimensions
    INTEGER :: i,j !< Count indexes
    INTEGER :: jstart, jend, nextj !< j counting variables
    TYPE(point) :: origPt !< Original point
    TYPE(point), DIMENSION(4) :: points !< xyz of 8 nearest neighbors
    REAL, DIMENSION(4) :: distSqrd !< distance between origPt and points(:)

    ! Set the inital fail values for indxI and indxJ
    indxI = 0
    indxJ = 0

    dimI = diag_global_grid%adimI
    dimJ = diag_global_grid%adimJ

    ! check to see if the 'fix' for the latitude index is needed
    IF ( diag_global_grid%aglo_lat(1,1) > &
         &diag_global_grid%aglo_lat(1,2) ) THEN
       ! reverse the j search
       jstart = dimJ-1
       jend = 1
       nextj = -1
    ELSE
       jstart = 0
       jend = dimJ-2
       nextJ = 1
    END IF

    ! find the I index
    iLoop: DO i=0, dimI-2
       IF (   diag_global_grid%aglo_lon(i,0) >&
            & diag_global_grid%aglo_lon(i+1,0) ) THEN
          ! We are at the 0 longitudal line
          IF (   (diag_global_grid%aglo_lon(i,0) <= lon .AND. lon <= 360.) .OR.&
               & (0. <= lon .AND. lon < diag_global_grid%aglo_lon(i+1, 0)) ) THEN
             indxI = i
             EXIT iLoop
          END IF
       ELSEIF ( diag_global_grid%aglo_lon(i,0) <= lon .AND.&
            &   lon <= diag_global_grid%aglo_lon(i+1,0) ) THEN
          indxI = i
          EXIT iLoop
       END IF
    END DO iLoop

    ! Find the J index
    IF ( indxI > 0 ) THEN
       jLoop: DO j=jstart, jend, nextj
          IF (   diag_global_grid%aglo_lat(indxI,j) <= lat .AND.&
               & lat <= diag_global_grid%aglo_lat(indxI,j+nextj) ) THEN
             indxJ = j
             EXIT jLoop
          END IF
       END DO jLoop
    END IF

    ! Make sure we have indexes in the correct range
    valid: IF ( (indxI <= 0 .OR. dimI-1 < indxI) .OR. &
         &      (indxJ <= 0 .OR. dimJ-1 < indxJ) ) THEN
       indxI = 0
       indxJ = 0
    ELSE ! indxI and indxJ are valid.
       ! Since we are looking for the closest grid point to the
       ! (lat,lon) point, we need to check the surrounding
       ! points.  The indexes for the variable points are as follows
       !
       ! 1---3
       ! |   |
       ! 2---4

       ! The original point
       origPt = latlon2xyz(lat,lon)

       ! Set the 'default' values for points(:) x,y,z to some large
       ! value.
       DO i=1, 4
          points(i)%x = 1.0e20
          points(i)%y = 1.0e20
          points(i)%z = 1.0e20
       END DO

       ! The original point
       origPt = latlon2xyz(lat,lon)

       points(1) = latlon2xyz(diag_global_grid%aglo_lat(indxI,indxJ),&
            &                 diag_global_grid%aglo_lon(indxI,indxJ))
       points(2) = latlon2xyz(diag_global_grid%aglo_lat(indxI,indxJ+nextj),&
            &                 diag_global_grid%aglo_lon(indxI,indxJ+nextj))
       points(3) = latlon2xyz(diag_global_grid%aglo_lat(indxI+1,indxJ+nextj),&
            &                 diag_global_grid%aglo_lon(indxI+1,indxJ+nextj))
       points(4) = latlon2xyz(diag_global_grid%aglo_lat(indxI+1,indxJ),&
            &                 diag_global_grid%aglo_lon(indxI+1,indxJ))

       ! Find the distance between the original point and the four
       ! grid points
       distSqrd = distanceSqrd(origPt, points)

       SELECT CASE (MINLOC(distSqrd,1))
       CASE ( 1 )
          indxI = indxI;
          indxJ = indxJ;
       CASE ( 2 )
          indxI = indxI;
          indxJ = indxJ+nextj;
       CASE ( 3 )
          indxI = indxI+1;
          indxJ = indxJ+nextj;
       CASE ( 4 )
          indxI = indxI+1;
          indxJ = indxJ;
       CASE DEFAULT
          indxI = 0;
          indxJ = 0;
       END SELECT

       ! If we are on tile 3 or 4, then the indxI and indxJ are
       ! reversed due to the transposed grids.
       IF (   diag_global_grid%tile_number == 4 .OR.&
            & diag_global_grid%tile_number == 5 ) THEN
          indxI_tmp = indxI
          indxI = indxJ
          indxJ = indxI_tmp
       END IF
    END IF valid

    ! Set the return value for the function
    find_equator_index_agrid = (/indxI, indxJ/)
  END FUNCTION find_equator_index_agrid

  !> @brief Return the (x,y,z) position of a given (lat,lon) point.
  !!
  !> Given a specific (lat, lon) point on the Earth, return the
  !!   corresponding (x,y,z) location.  The return of latlon2xyz
  !!   will be either a scalar or an array of the same size as lat
  !!   and lon.
  !! @return The return of latlon2xyz
  !!   will be either a scalar or an array of the same size as lat
  !!   and lon.
  PURE ELEMENTAL TYPE(point) FUNCTION latlon2xyz(lat, lon)
    REAL, INTENT(in) :: lat !< The latitude of the (x,y,z) location to find.  <TT>lat</TT>
                            !! can be either a scalar or array.  <TT>lat</TT> must be of the
                            !! same rank / size as <TT>lon</TT>.  This function assumes
                            !! <TT>lat</TT> is in the range [-90,90].
    REAL, INTENT(in) :: lon !< The longitude of the (x,y,z) location to find.  <TT>lon</TT>
                            !! can be either a scalar or array.  <TT>lon</TT> must be of the
                            !! same rank / size as <TT>lat</TT>.  This function assumes
                            !! <TT>lon</TT> is in the range [0,360].

    ! lat/lon angles in radians
    REAL :: theta !< lat angles in radians
    REAL :: phi   !< lon angles in radians

    ! Convert the lat lon values to radians The lat values passed in
    ! are in the range [-90,90], but we need to have a radian range
    ! [0,pi], where 0 is at the north pole.  This is the reason for
    ! the subtraction from 90
    theta = deg2rad(90.-lat)
    phi = deg2rad(lon)

    ! Calculate the x,y,z point
    latlon2xyz%x = RADIUS * SIN(theta) * COS(phi)
    latlon2xyz%y = RADIUS * SIN(theta) * SIN(phi)
    latlon2xyz%z = RADIUS * COS(theta)
  END FUNCTION latlon2xyz

  !> @brief Find the distance between two points in the Cartesian
  !!   coordinate space.
  !!
  !> <TT>distanceSqrd</TT> will find the distance squared between
  !!   two points in the xyz coordinate space.  <TT>pt1</TT> and <TT>
  !!   pt2</TT> can either be both scalars, both arrays of the same
  !!   size, or one a scalar and one an array.  The return value
  !!   will be a scalar or array of the same size as the input array.
  !! @return The return value will be a scalar or array of the same size as the input array.
  PURE ELEMENTAL REAL FUNCTION distanceSqrd(pt1, pt2)
    TYPE(point), INTENT(in) :: pt1, pt2

    distanceSqrd = (pt1%x-pt2%x)**2 +&
         &         (pt1%y-pt2%y)**2 +&
         &         (pt1%z-pt2%z)**2
  END FUNCTION distanceSqrd

  !> @brief Find the distance, along the geodesic, between two points.
  !!
  !> <TT>gCirDistance</TT> will find the distance, along the geodesic, between two points defined
  !!   by the (lat,lon) position of each point.
  !! @return real
  PURE ELEMENTAL REAL FUNCTION gCirDistance(lat1, lon1, lat2, lon2)
    REAL, INTENT(in) :: lat1, lat2, lon1, lon2

    REAL :: theta1, theta2
    REAL :: deltaLambda !< Difference in longitude angles, in radians.
    REAL :: deltaTheta !< Difference in latitude angels, in radians.

    theta1 = deg2rad(lat1)
    theta2 = deg2rad(lat2)
    deltaLambda = deg2rad(lon2-lon1)
    deltaTheta = deg2rad(lat2-lat1)

    gCirDistance = RADIUS * 2. * ASIN(SQRT((SIN(deltaTheta/2.))**2 + COS(theta1)*COS(theta2)*(SIN(deltaLambda/2.))**2))
  END FUNCTION gCirDistance

  !> @brief Find the i,j indices and distance of the a-grid point nearest to
  !! the inputted lat,lon point.
  SUBROUTINE find_nearest_agrid_index(lat, &
                                      lon, &
                                      minI, &
                                      minJ, &
                                      minimum_distance)

    !Inputs/outputs
    REAL,INTENT(IN) :: lat
    REAL,INTENT(IN) :: lon
    INTEGER,INTENT(OUT) :: minI
    INTEGER,INTENT(OUT) :: minJ
    REAL,INTENT(OUT) :: minimum_distance

    !Local variables
    REAL :: llat
    REAL :: llon
    INTEGER :: j
    INTEGER :: i
    REAL :: dist

    !Since the poles have an non-unique longitude value, make a small
    !correction if looking for one of the poles.
    IF (lat .EQ. 90.0) THEN
       llat = lat - .1
    ELSEIF (lat .EQ. -90.0) THEN
       llat = lat + .1
    ELSE
       llat = lat
    END IF
    llon = lon

    !Loop through non-halo points.  Calculate the distance
    !between each a-grid point and the point that we
    !are seeking.  Store the minimum distance and its
    !corresponding i,j indices.
    minI = 0
    minJ = 0
    minimum_distance = 2.0*RADIUS*3.141592653
    DO j = 1,diag_global_grid%adimJ-2
        DO i = 1,diag_global_grid%adimI-2
            dist = gCirDistance(llat, &
                                llon, &
                                diag_global_grid%aglo_lat(i,j), &
                                diag_global_grid%aglo_lon(i,j))
            IF (dist .LT. minimum_distance) THEN

                !These number shouldn't be hardcoded, but they have to
                !match the ones in diag_grid_init.
                if (diag_global_grid%tile_number .eq. 4 .or. &
                        diag_global_grid%tile_number .eq. 5) then

                    !Because of transpose in diag_grid_init.
                    minI = j
                    minJ = i

                else
                    minI = i
                    minJ = j
                endif
                minimum_distance = dist
            ENDIF
        ENDDO
    ENDDO

    !Check that valid i,j indices have been found.
    IF (minI .EQ. 0 .OR. minJ .EQ. 0) THEN
        call error_mesg("find_nearest_agrid_index", &
                        "A minimum distance was not found.", &
                        FATAL)
    ENDIF

  END SUBROUTINE find_nearest_agrid_index

END MODULE diag_grid_mod
!> @}
! close documentation grouping
