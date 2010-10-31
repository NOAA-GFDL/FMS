#include <fms_platform.h>

MODULE diag_grid_mod
  ! <CONTACT EMAIL="seth.underwood@noaa.gov">
  !   Seth Underwood
  ! </CONTACT>
  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/" />
  ! <OVERVIEW>
  !   <TT>diag_grid_mod</TT> is a set of procedures to work with the
  !   model's global grid to allow regional output.
  ! </OVERVIEW>
  ! <DESCRIPTION>
  !   <TT>diag_grid_mod</TT> contains useful utilities for dealing
  !   with, mostly, regional output for grids other than the standard
  !   lat/lon grid.  This module contains three public procedures <TT>
  !   diag_grid_init</TT>, which is shared globably in the <TT>
  !   diag_manager_mod</TT>, <TT>diag_grid_end</TT> which will free
  !   up memory used during the register field calls, and
  !   <TT>get_local_indexes</TT>.  The <TT>send_global_grid</TT>
  !   procedure is called by the model that creates the global grid.
  !   <TT>send_global_grid</TT> needs to be called before any fields
  !   are registered that will output only regions.  <TT>get_local_indexes</TT>
  !   is to be called by the <TT>diag_manager_mod</TT> to discover the
  !   global indexes defining a subregion on the tile.
  !
  !   <B>Change Log</B>
  !   <DL>
  !     <DT>September 2009</DT>
  !     <DD>
  !       <UL>
  !         <LI>Single point region in Cubed Sphere</LI>
  !         <LI>Single tile regions in the cubed sphere</LI>
  !       </UL>
  !     </DD> 
  !   </DL>
  ! </DESCRIPTION>

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
  USE mpp_mod, ONLY: mpp_root_pe, mpp_npes, mpp_max
  Use mpp_domains_mod, ONLY: domain2d, mpp_get_tile_id,&
       & mpp_get_ntile_count, mpp_get_compute_domains

  IMPLICIT NONE

  ! Parameters
  CHARACTER(len=128), PARAMETER :: version =&
       & '$Id: diag_grid.F90,v 18.0.2.6 2010/04/12 21:25:34 sdu Exp $'
  CHARACTER(len=128), PARAMETER :: tagname =&
       & '$Name: riga_201012 $'

  ! Derived data types
  ! <PRIVATE>
  ! <TYPE NAME="diag_global_grid_type">
  !   <DESCRIPTION>
  !     Contains the model's global grid data, and other grid information.
  !   </DESCRIPTION>
  !   <DATA NAME="glo_lat" TYPE="REAL, _ALLOCATABLE, DIMENSION(:,:)">
  !     The latitude values on the global grid.
  !   </DATA>
  !   <DATA NAME="glo_lon" TYPE="REAL, _ALLOCATABLE, DIMENSION(:,:)">
  !     The longitude values on the global grid.
  !   </DATA>
  !   <DATA NAME="aglo_lat" TYPE="REAL, _ALLOCATABLE, DIMENSION(:,:)">
  !     The latitude values on the global a-grid.  Here we expect isc-1:iec+1 and
  !     jsc=1:jec+1 to be passed in.
  !   </DATA>
  !   <DATA NAME="aglo_lon" TYPE="REAL, _ALLOCATABLE, DIMENSION(:,:)">
  !     The longitude values on the global a-grid.  Here we expec isc-1:iec+j and
  !     jsc-1:jec+1 to be passed in.
  !   </DATA>
  !   <DATA NAME="myXbegin" TYPE="INTEGER">
  !     The starting index of the compute domain on the current PE.
  !   </DATA>
  !   <DATA NAME="myYbegin" TYPE="INTEGER">
  !     The starting index of the compute domain on the cureent PE.
  !   </DATA>
  !   <DATA NAME="dimI" TYPE="INTEGER">
  !     The dimension of the global grid in the 'i' / longitudal direction.
  !   </DATA>
  !   <DATA NAME="dimJ" TYPE="INTEGER">
  !     The dimension of the global grid in the 'j' / latitudal direction.
  !   </DATA>
  !   <DATA NAME="adimI" TYPE="INTEGER">
  !     The dimension of the global a-grid in the 'i' / longitudal direction.  Again,
  !     the expected dimension for diag_grid_mod is isc-1:iec+1.
  !   </DATA>
  !   <DATA NAME="adimJ" TYPE="INTEGER">
  !     The dimension of the global a-grid in the 'j' / latitudal direction.  Again,
  !     the expected dimension for diag_grid_mod is jsc-1:jec+1.
  !   </DATA>
  !   <DATA NAME="tile_number" TYPE="INTEGER">
  !     The tile the <TT>glo_lat</TT> and <TT>glo_lon</TT> define.
  !   </DATA>
  !   <DATA NAME="ntimes" TYPE="INTEGER">
  !     The number of tiles.
  !   </DATA>
  !   <DATA NAME="peStart" TYPE="INTEGER">
  !     The starting PE number for the current tile.
  !   </DATA>
  !   <DATA NAME="peEnd" TYPE="INTEGER">
  !     The ending PE number for the current tile.
  !   </DATA>
  !   <DATA NAME="grid_type" TYPE="CHARACTER(len=128)">
  !     The global grid type.
  !   </DATA>
  TYPE :: diag_global_grid_type
     REAL, _ALLOCATABLE, DIMENSION(:,:) :: glo_lat, glo_lon
     REAL, _ALLOCATABLE, DIMENSION(:,:) :: aglo_lat, aglo_lon
     INTEGER :: myXbegin, myYbegin
     INTEGER :: dimI, dimJ
     INTEGER :: adimI, adimJ
     INTEGER :: tile_number
     INTEGER :: ntiles
     INTEGER :: peStart, peEnd
     CHARACTER(len=128) :: grid_type
  END TYPE diag_global_grid_type
  ! </TYPE>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <TYPE NAME="point">
  !   <DESCRIPTION>
  !      Private point type to hold the (x,y,z) location for a (lat,lon)
  !      location.
  !   </DESCRIPTION>
  !   <DATA NAME="x" TYPE="REAL">
  !     The x value of the (x,y,z) coordinates.
  !   </DATA>
  !   <DATA NAME="y" TYPE="REAL">
  !     The y value of the (x,y,z) coordinates.
  !   </DATA>
  !   <DATA NAME="z" TYPE="REAL">
  !     The z value of the (x,y,z) coordinates.
  !   </DATA>
  TYPE :: point
     REAL :: x,y,z
  END TYPE point
  ! </TYPE>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <DATA NAME="diag_global_grid" TYPE="TYPE(diag_global_grid_type)">
  !   Variable to hold the global grid data
  ! </DATA>
  ! </PRIVATE>
  TYPE(diag_global_grid_type) :: diag_global_grid

  ! <PRIVATE>
  ! <DATA NAME="diag_grid_initialized" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   Indicates if the diag_grid_mod has been initialized.
  ! </DATA>
  ! </PRIVATE>
  LOGICAL :: diag_grid_initialized = .FALSE.

  PRIVATE
  PUBLIC :: diag_grid_init, diag_grid_end, get_local_indexes,  &
            get_local_indexes2

CONTAINS

  ! <SUBROUTINE NAME="diag_grid_init">
  !   <OVERVIEW>
  !     Send the global grid to the <TT>diag_manager_mod</TT> for
  !     regional output. 
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE diag_grid_init(domain, glo_lat, glo_lon, aglo_lat, aglo_lon)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     In order for the diag_manager to do regional output for grids
  !     other than the standard lat/lon grid, the <TT>
  !     diag_manager_mod</TT> needs to know the the latitude and
  !     longitude values for the entire global grid.  This procedure
  !     is the mechanism the models will use to share their grid with
  !     the diagnostic manager.
  !     
  !     This procedure needs to be called after the grid is created,
  !     and before the first call to register the fields.
  !   </DESCRIPTION>
  !   <IN NAME="domain" TYPE="INTEGER">
  !     The domain to which the grid data corresponds.
  !   </IN>
  !   <IN NAME="glo_lat" TYPE="REAL, DIMENSION(:,:)">
  !     The latitude information for the grid tile.
  !   </IN>
  !   <IN NAME="glo_lon" TYPE="REAL, DIMENSION(:,:)">
  !     The longitude information for the grid tile.
  !   </IN>
  !   <IN NAME="aglo_lat" TYPE="REAL, DIMENSION(:,:)">
  !     The latitude information for the a-grid tile.
  !   </IN>
  !   <IN NAME="aglo_lon" TYPE="REAL, DIMENSION(:,:)">
  !     The longitude information for the a-grid tile.
  !   </IN>
  SUBROUTINE diag_grid_init(domain, glo_lat, glo_lon, aglo_lat, aglo_lon)
    TYPE(domain2d), INTENT(in) :: domain
    REAL, INTENT(in), DIMENSION(:,:) :: glo_lat, glo_lon
    REAL, INTENT(in), DIMENSION(:,:) :: aglo_lat, aglo_lon

    INTEGER, DIMENSION(1) :: tile
    INTEGER :: ntiles
    INTEGER :: stat
    INTEGER :: i_dim, j_dim
    INTEGER :: ai_dim, aj_dim
    INTEGER, DIMENSION(2) :: latDim, lonDim
    INTEGER, DIMENSION(2) :: alatDim, alonDim
    INTEGER :: myPe, npes, npesPerTile
    INTEGER, ALLOCATABLE, DIMENSION(:) :: xbegin, xend, ybegin, yend

    ! Write the version and tagname to the logfile
    CALL write_version_number(version, tagname)

    ! Verify all allocatable / pointers for diag_global_grid hare not
    ! allocated / associated.
    IF ( ALLOCATED(xbegin) ) DEALLOCATE(xbegin)
    IF ( ALLOCATED(ybegin) ) DEALLOCATE(ybegin)
    IF ( ALLOCATED(xend) ) DEALLOCATE(xend)
    IF ( ALLOCATED(yend) ) DEALLOCATE(yend)

    ! What is my PE
    myPe = mpp_pe() + 1
    
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
       IF ( tile(1) == 4 .OR. tile(1) == 5 ) THEN
          ! These tiles need to be transposed.
          i_dim = latDim(2)
          j_dim = latDim(1)
       ELSE 
          i_dim = latDim(1)
          j_dim = latDim(2)
       END IF
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
    IF (   _ALLOCATED(diag_global_grid%glo_lat) .OR.&
         & _ALLOCATED(diag_global_grid%glo_lon) ) THEN
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
    IF (   _ALLOCATED(diag_global_grid%aglo_lat) .OR.&
         & _ALLOCATED(diag_global_grid%aglo_lon) ) THEN
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
       diag_global_grid%glo_lat = TRANSPOSE(glo_lat)
       diag_global_grid%glo_lon = TRANSPOSE(glo_lon)
       diag_global_grid%aglo_lat = TRANSPOSE(aglo_lat)
       diag_global_grid%aglo_lon = TRANSPOSE(aglo_lon)
    ELSE
       diag_global_grid%glo_lat = glo_lat
       diag_global_grid%glo_lon = glo_lon
       diag_global_grid%aglo_lat = aglo_lat
       diag_global_grid%aglo_lon = aglo_lon
    END IF
    diag_global_grid%dimI = i_dim
    diag_global_grid%dimJ = j_dim
    diag_global_grid%adimI = ai_dim
    diag_global_grid%adimJ = aj_dim
    diag_global_grid%tile_number = tile(1)
    diag_global_grid%ntiles = ntiles
    diag_global_grid%myXbegin = xbegin(myPe)
    diag_global_grid%myYbegin = ybegin(myPe)

    ! Unallocate arrays used here
    DEALLOCATE(xbegin)
    DEALLOCATE(ybegin)
    DEALLOCATE(xend)
    DEALLOCATE(yend)
  END SUBROUTINE diag_grid_init
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="diag_grid_end">
  !   <OVERVIEW>
  !     Unallocate the diag_global_grid variable.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE diag_grid_end()
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     The <TT>diag_global_grid</TT> variable is only needed during
  !     the register field calls, and then only if there are fields
  !     requestion regional output.  Once all the register fields
  !     calls are complete (before the first <TT>send_data</TT> call
  !     this procedure can be called to free up memory.
  !   </DESCRIPTION>
  SUBROUTINE diag_grid_end()
    
    IF ( diag_grid_initialized ) THEN
       ! De-allocate grid
       IF ( _ALLOCATED(diag_global_grid%glo_lat) ) THEN
          DEALLOCATE(diag_global_grid%glo_lat)
       ELSE IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_grid_mod::diag_grid_end',&
               &'diag_global_grid%glo_lat was not allocated.', WARNING)
       END IF
       
       IF ( _ALLOCATED(diag_global_grid%glo_lon) ) THEN
          DEALLOCATE(diag_global_grid%glo_lon)
       ELSE IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_grid_mod::diag_grid_end',&
               &'diag_global_grid%glo_lon was not allocated.', WARNING)
       END IF
       ! De-allocate a-grid
       IF ( _ALLOCATED(diag_global_grid%aglo_lat) ) THEN
          DEALLOCATE(diag_global_grid%aglo_lat)
       ELSE IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_grid_mod::diag_grid_end',&
               &'diag_global_grid%aglo_lat was not allocated.', WARNING)
       END IF
       
       IF ( _ALLOCATED(diag_global_grid%aglo_lon) ) THEN
          DEALLOCATE(diag_global_grid%aglo_lon)
       ELSE IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_grid_mod::diag_grid_end',&
               &'diag_global_grid%aglo_lon was not allocated.', WARNING)
       END IF
       
       diag_grid_initialized = .FALSE.
    END IF
  END SUBROUTINE diag_grid_end
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="get_local_indexes">
  !   <OVERVIEW>
  !     Find the local start and local end indexes on the local PE
  !     for regional output.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_local_indexes(latStart, latEnd, lonStart,
  !     lonEnd, istart, iend, jstart, jend)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Given a defined region, find the local indexes on the local
  !     PE surrounding the region.
  !   </DESCRIPTION>
  !   <IN NAME="latStart" TYPE="REAL">
  !     The minimum latitude value defining the region.  This value
  !     must be less than latEnd, and be in the range [-90,90]
  !   </IN>
  !   <IN NAME="latEnd" TYPE="REAL">
  !     The maximum latitude value defining the region.  This value
  !     must be greater than latStart, and be in the range [-90,90]
  !   </IN>
  !   <IN NAME="lonStart" TYPE="REAL">
  !     The western most longitude value defining the region.
  !     Possible ranges are either [-180,180] or [0,360].
  !   </IN>
  !   <IN NAME="lonEnd" TYPE="REAL">
  !     The eastern most longitude value defining the region.
  !     Possible ranges are either [-180,180] or [0,360].
  !   </IN>
  !   <OUT NAME="istart" TYPE="INTEGER">
  !     The local start index on the local PE in the 'i' direction.
  !   </OUT>
  !   <OUT NAME="iend" TYPE="INTEGER">
  !     The local end index on the local PE in the 'i' direction.
  !   </OUT>
  !   <OUT NAME="jstart" TYPE="INTEGER">
  !     The local start index on the local PE in the 'j' direction.
  !   </OUT>
  !   <OUT NAME="jend" TYPE="INTEGER">
  !     The local end index on the local PE in the 'j' direction.
  !   </OUT>
  SUBROUTINE get_local_indexes(latStart, latEnd, lonStart, lonEnd,&
       & istart, iend, jstart, jend)
    REAL, INTENT(in) :: latStart, lonStart !< lat/lon start angles
    REAL, INTENT(in) :: latEnd, lonEnd !< lat/lon end angles
    INTEGER, INTENT(out) :: istart, jstart !< i/j start indexes
    INTEGER, INTENT(out) :: iend, jend !< i/j end indexes

    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: indexes
    INTEGER :: myTile, ntiles, i, j

    IF ( .NOT. diag_grid_initialized )&
         & CALL error_mesg('diag_grid_mod::get_local_indexes',&
         &'Module not initialized, first initialze module with a call &
         &to diag_grid_init', FATAL)
    
    myTile = diag_global_grid%tile_number
    ntiles = diag_global_grid%ntiles

    ! Allocate the indexes array, and initialize to zero.  (Useful for
    ! reduction later on.)
    IF ( ALLOCATED(indexes) ) DEALLOCATE(indexes)
    ALLOCATE(indexes(ntiles,4,2))
    indexes = 0 

    ! There will be four points to define a region, find all four.
    ! Need to call the correct function depending on if the tile is a
    ! pole tile or not.
    !
    ! Also, if looking for a single point, then use the a-grid
    IF ( latStart == latEnd .AND. lonStart == lonEnd ) THEN
       ! single point
       IF ( MOD(diag_global_grid%tile_number,3) == 0 ) THEN
          indexes(myTile,1,:) = find_pole_index_agrid(latStart,lonStart)
          indexes(myTile,2,:) = indexes(myTile,1,:)
          indexes(myTile,3,:) = indexes(myTile,1,:)
          indexes(myTile,4,:) = indexes(myTile,1,:)
       ELSE
          indexes(myTile,1,:) = find_equator_index_agrid(latStart,lonStart)
          indexes(myTile,2,:) = indexes(myTile,1,:)
          indexes(myTile,3,:) = indexes(myTile,1,:)
          indexes(myTile,4,:) = indexes(myTile,1,:)
       END IF
    ELSE
       ! multi-point
       IF ( MOD(diag_global_grid%tile_number,3) == 0 ) THEN
          ! Pole tile
          indexes(myTile,1,:) = find_pole_index(latStart, lonStart)
          indexes(myTile,2,:) = find_pole_index(latStart, lonEnd)
          indexes(myTile,3,:) = find_pole_index(latEnd, lonStart)
          indexes(myTile,4,:) = find_pole_index(latEnd, lonEnd)
       ELSE
          indexes(myTile,1,:) = find_equator_index(latStart, lonStart)
          indexes(myTile,2,:) = find_equator_index(latStart, lonEnd)
          indexes(myTile,3,:) = find_equator_index(latEnd, lonStart)
          indexes(myTile,4,:) = find_equator_index(latEnd, lonEnd)
       END IF
    END IF
    
    WHERE ( indexes(:,:,1) .NE. 0 )
       indexes(:,:,1) = indexes(:,:,1) + diag_global_grid%myXbegin - 1
    END WHERE
    WHERE ( indexes(:,:,2) .NE. 0 )
       indexes(:,:,2) = indexes(:,:,2) + diag_global_grid%myYbegin - 1
    END WHERE
    
    DO j = 1, 6 ! Each tile.
       DO i = 1, 4
          CALL mpp_max(indexes(j,i,1))
          CALL mpp_max(indexes(j,i,2))
       END DO
    END DO
    
    ! Are there any indexes found on this tile?
    ! Check if all points are on this tile
    ! Works since the find index functions return 0 if not found.
    IF (   PRODUCT(indexes(myTile,:,1)) /= 0 .OR.&
         & PRODUCT(indexes(myTile,:,2)) /= 0  ) THEN
       istart = MINVAL(indexes(myTile,:,1))
       jstart = MINVAL(indexes(myTile,:,2))
       iend = MAXVAL(indexes(myTile,:,1))
       jend = MAXVAL(indexes(myTile,:,2))
    ELSE
       istart = 0
       jstart = 0
       iend = 0
       jend = 0
    END IF

    DEALLOCATE(indexes)
  END SUBROUTINE get_local_indexes
  ! </SUBROUTINE>
  
  ! <SUBROUTINE NAME="get_local_indexes2">
  !   <OVERVIEW>
  !     Find the indices of the nearest grid point of the a-grid to the 
  !     specified (lon,lat) location on the local PE. if desired point not 
  !     within domain of local PE, return (0,0) as the indices. 
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_local_indexes2 (lat, lon, iindex, jindex)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Given a specified location, find the nearest a-grid indices on 
  !     the local PE.
  !   </DESCRIPTION>
  !   <IN NAME="lat" TYPE="REAL">
  !     The requested latitude.  This value must be in the range [-90,90]
  !   </IN>
  !   <IN NAME="lon" TYPE="REAL">
  !     The requested longitude.
  !     Possible ranges are either [-180,180] or [0,360].
  !   </IN>
  !   <OUT NAME="iindex" TYPE="INTEGER">
  !     The local index on the local PE in the 'i' direction.
  !   </OUT>
  !   <OUT NAME="jindex" TYPE="INTEGER">
  !     The local index on the local PE in the 'j' direction.
  !   </OUT>
  SUBROUTINE get_local_indexes2(lat, lon, iindex, jindex)
    REAL, INTENT(in) :: lat, lon !< lat/lon location    
    INTEGER, INTENT(out) :: iindex, jindex !< i/j indexes

    INTEGER  :: indexes(2)
    INTEGER :: i, j

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
    if (iindex ==  diag_global_grid%adimI -1 .or.&
        jindex ==  diag_global_grid%adimJ -1 ) then
      iindex = 0
      jindex = 0
    endif
           
  END SUBROUTINE get_local_indexes2
  ! </SUBROUTINE>

  ! <PRIVATE>
  ! <FUNCTION NAME="rad2deg">
  !   <OVERVIEW>
  !     Convert and angle in radian to degrees.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE ELEMENTAL REAL FUNCTION rad2deg(angle)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Given a scalar, or an array of angles in radians this
  !     function will return a scalar or array (of the same
  !     dimension) of angles in degrees.
  !   </DESCRIPTION>
  !   <IN NAME="angle" TYPE="REAL">
  !     Scalar or array of angles in radians.
  !   </IN>
  !   <OUT NAME="rad2deg" TYPE="REAL">
  !     Scalar or array (depending on the size of angle) of angles in
  !     degrees.
  !   </OUT>
  PURE ELEMENTAL REAL FUNCTION rad2deg(angle)
    REAL, INTENT(in) :: angle

    rad2deg = RAD_TO_DEG * angle
  END FUNCTION rad2deg
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="deg2rad">
  !   <OVERVIEW>
  !     Convert an angle in degrees to radians.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE ELEMENTAL REAL FUNCTION deg2rad(angle)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Given a scalar, or an array of angles in degrees this
  !     function will return a scalar or array (of the same
  !     dimension) of angles in radians.
  !   </DESCRIPTION>
  !   <IN NAME="angle" TYPE="REAL">
  !     Scalar or array of angles in degrees.
  !   </IN>
  PURE ELEMENTAL REAL FUNCTION deg2rad(angle)
    REAL, INTENT(in) :: angle

    deg2rad = DEG_TO_RAD * angle
  END FUNCTION deg2rad
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="find_pole_index">
  !   <OVERVIEW>
  !     Return the closest index (i,j) to the given (lat,lon) point.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE FUNCTION find_pole_index(lat, lon)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     This function searches a pole grid tile looking for the grid point
  !     closest to the give (lat, lon) location, and returns the i
  !     and j indexes of the point.
  !   </DESCRIPTION>
  !   <IN NAME="lat" TYPE="REAL">
  !     Latitude location
  !   </IN>
  !   <IN NAME="lon" TYPE="REAL">
  !     Longitude location
  !   </IN>
  !   <OUT NAME="find_pole_index" TYPE="INTEGER, DIMENSION(2)">
  !     The (i, j) location of the closest grid to the given (lat,
  !     lon) location.
  !   </OUT>
  PURE FUNCTION find_pole_index(lat, lon)
    INTEGER, DIMENSION(2) :: find_pole_index
    REAL, INTENT(in) :: lat, lon
    
    INTEGER :: indxI, indxJ !< Indexes to be returned.
    INTEGER :: dimI, dimJ !< Size of the grid dimensions
    INTEGER :: i,j !< Count indexes
    INTEGER :: nearestCorner !< index of the nearest corner
    INTEGER , DIMENSION(4,2) :: ijArray !< indexes of the cornerPts and pntDistances arrays
    REAL :: llat, llon !< Corrected lat and lon location (if looking for pole point.)
    REAL :: maxCtrDist !< maximum distance to the origPt to corner
    REAL, DIMENSION(4) :: pntDistances !< distance from origPt to corner
    REAL, DIMENSION(4,2) :: cornerPts !< Corner points using (lat,lon)
    TYPE(point) :: origPt !< Original point
    TYPE(point), DIMENSION(9) :: points !< xyz of 8 nearest neighbors
    REAL, DIMENSION(9) :: distSqrd !< distance between origPt and points(:)

    ! Set the inital fail values for indxI and indxJ
    indxI = 0
    indxJ = 0
    
    dimI = diag_global_grid%dimI
    dimJ = diag_global_grid%dimJ

    ! Since the poles have an non-unique longitude value, make a small correction if looking for one of the poles.
    IF ( lat == 90.0 ) THEN
       llat = lat - .1
    ELSE IF ( lat == -90.0 ) THEN
       llat = lat + .1
    ELSE
       llat = lat
    END IF
    llon = lon

    iLoop: DO i=1, dimI-1
       jLoop: DO j = 1, dimJ-1
          ! Get the lat,lon for the four corner points.
          cornerPts = RESHAPE( (/ diag_global_grid%glo_lat(i,  j),  diag_global_grid%glo_lon(i,  j),&
               &                  diag_global_grid%glo_lat(i+1,j+1),diag_global_grid%glo_lon(i+1,j+1),&
               &                  diag_global_grid%glo_lat(i+1,j),  diag_global_grid%glo_lon(i+1,j),&
               &                  diag_global_grid%glo_lat(i,  j+1),diag_global_grid%glo_lon(i,  j+1) /),&
               &               (/ 4, 2 /), ORDER=(/2,1/) )

          ! Find the maximum half distance of the corner points
          maxCtrDist = MAX(gCirDistance(cornerPts(1,1),cornerPts(1,2), cornerPts(2,1),cornerPts(2,2)),&
               &           gCirDistance(cornerPts(3,1),cornerPts(3,2), cornerPts(4,1),cornerPts(4,2)))/2
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
    valid: IF (  (indxI <= 0 .OR. dimI < indxI) .OR. &
         &       (indxJ <= 0 .OR. dimJ < indxJ) ) THEN
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

       ! The original point
       origPt = latlon2xyz(lat,lon)

       ! Set the 'default' values for points(:) x,y,z to some large
       ! value.
       DO i=1, 9
          points(i)%x = 1.0e20
          points(i)%y = 1.0e20
          points(i)%z = 1.0e20
       END DO

       ! All the points around the i,j indexes
       IF ( indxI > 1 ) THEN
          points(1) = latlon2xyz(diag_global_grid%glo_lat(indxI-1,indxJ+1),&
               &                 diag_global_grid%glo_lon(indxI-1,indxJ+1))
          points(2) = latlon2xyz(diag_global_grid%glo_lat(indxI-1,indxJ),&
               &                 diag_global_grid%glo_lon(indxI-1,indxJ))
          IF ( indxJ > 1 ) THEN
             points(3) = latlon2xyz(diag_global_grid%glo_lat(indxI-1,indxJ-1),&
                  &                 diag_global_grid%glo_lon(indxI-1,indxJ-1))
          END IF
       END IF
       points(4) = latlon2xyz(diag_global_grid%glo_lat(indxI,  indxJ+1),&
            &                 diag_global_grid%glo_lon(indxI,  indxJ+1))
       points(5) = latlon2xyz(diag_global_grid%glo_lat(indxI,  indxJ),&
            &                 diag_global_grid%glo_lon(indxI,  indxJ))
       IF ( indxJ > 1 ) THEN
          points(6) = latlon2xyz(diag_global_grid%glo_lat(indxI,  indxJ-1),&
               &                 diag_global_grid%glo_lon(indxI,  indxJ-1))
       END IF
       points(7) = latlon2xyz(diag_global_grid%glo_lat(indxI+1,indxJ+1),&
            &                 diag_global_grid%glo_lon(indxI+1,indxJ+1))
       points(8) = latlon2xyz(diag_global_grid%glo_lat(indxI+1,indxJ),&
            &                 diag_global_grid%glo_lon(indxI+1,indxJ))
       IF ( indxJ > 1 ) THEN
          points(9) = latlon2xyz(diag_global_grid%glo_lat(indxI+1,indxJ-1),&
               &                 diag_global_grid%glo_lon(indxI+1,indxJ-1))
       END IF
          
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
    find_pole_index = (/indxI, indxJ/)
  END FUNCTION find_pole_index
  ! </FUNCTION>
  ! </PRIVATE>
  
  ! <PRIVATE>
  ! <FUNCTION NAME="find_equator_index">
  !   <OVERVIEW>
  !     Return the closest index (i,j) to the given (lat,lon) point.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE FUNCTION find_equator_index(lat, lon)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     This function searches a equator grid tile looking for the grid point
  !     closest to the give (lat, lon) location, and returns the i
  !     and j indexes of the point.
  !   </DESCRIPTION>
  !   <IN NAME="lat" TYPE="REAL">
  !     Latitude location
  !   </IN>
  !   <IN NAME="lon" TYPE="REAL">
  !     Longitude location
  !   </IN>
  !   <OUT NAME="find_equator_index" TYPE="INTEGER, DIMENSION(2)">
  !     The (i, j) location of the closest grid to the given (lat,
  !     lon) location.
  !   </OUT>
  PURE FUNCTION find_equator_index(lat, lon)
    INTEGER, DIMENSION(2) :: find_equator_index
    REAL, INTENT(in) :: lat, lon
    
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
    
    dimI = diag_global_grid%dimI
    dimJ = diag_global_grid%dimJ

    ! check to see if the 'fix' for the latitude index is needed
    IF ( diag_global_grid%glo_lat(1,1) > &
         &diag_global_grid%glo_lat(1,2) ) THEN
       ! reverse the j search
       jstart = dimJ
       jend = 2
       nextj = -1
    ELSE
       jstart = 1
       jend = dimJ-1
       nextJ = 1
    END IF

    ! find the I index
    iLoop: DO i=1, dimI-1
       IF (   diag_global_grid%glo_lon(i,1) >&
            & diag_global_grid%glo_lon(i+1,1) ) THEN
          ! We are at the 0 longitudal line
          IF (   (diag_global_grid%glo_lon(i,1) <= lon .AND. lon <= 360) .OR.&
               & (0 <= lon .AND. lon < diag_global_grid%glo_lon(i+1, 1)) ) THEN
             indxI = i
             EXIT iLoop
          END IF
       ELSEIF ( diag_global_grid%glo_lon(i,1) <= lon .AND.&
            &   lon <= diag_global_grid%glo_lon(i+1,1) ) THEN
          indxI = i
          EXIT iLoop
       END IF
    END DO iLoop
    
    ! Find the J index
    IF ( indxI > 0 ) THEN
       jLoop: DO j=jstart, jend, nextj
          IF (   diag_global_grid%glo_lat(indxI,j) <= lat .AND.&
               & lat <= diag_global_grid%glo_lat(indxI,j+nextj) ) THEN
             indxJ = j
             EXIT jLoop
          END IF
       END DO jLoop
    END IF

    ! Make sure we have indexes in the correct range
    valid: IF ( (indxI <= 0 .OR. dimI < indxI) .OR. &
         &      (indxJ <= 0 .OR. dimJ < indxJ) ) THEN
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

       points(1) = latlon2xyz(diag_global_grid%glo_lat(indxI,indxJ),&
            &                 diag_global_grid%glo_lon(indxI,indxJ))
       points(2) = latlon2xyz(diag_global_grid%glo_lat(indxI,indxJ+nextj),&
            &                 diag_global_grid%glo_lon(indxI,indxJ+nextj))
       points(3) = latlon2xyz(diag_global_grid%glo_lat(indxI+1,indxJ+nextj),&
            &                 diag_global_grid%glo_lon(indxI+1,indxJ+nextj))
       points(4) = latlon2xyz(diag_global_grid%glo_lat(indxI+1,indxJ),&
            &                 diag_global_grid%glo_lon(indxI+1,indxJ))
  
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
    find_equator_index = (/indxI, indxJ/)
  END FUNCTION find_equator_index
  ! </FUNCTION>
  ! </PRIVATE>
  
  ! <PRIVATE>
  ! <FUNCTION NAME="find_pole_index_agrid">
  !   <OVERVIEW>
  !     Return the closest index (i,j) to the given (lat,lon) point.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE FUNCTION find_pole_index_agrid(lat, lon)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     This function searches a pole a-grid tile looking for the grid point
  !     closest to the give (lat, lon) location, and returns the i
  !     and j indexes of the point.
  !   </DESCRIPTION>
  !   <IN NAME="lat" TYPE="REAL">
  !     Latitude location
  !   </IN>
  !   <IN NAME="lon" TYPE="REAL">
  !     Longitude location
  !   </IN>
  !   <OUT NAME="find_pole_index" TYPE="INTEGER, DIMENSION(2)">
  !     The (i, j) location of the closest grid to the given (lat,
  !     lon) location.
  !   </OUT>
  PURE FUNCTION find_pole_index_agrid(lat, lon)
    INTEGER, DIMENSION(2) :: find_pole_index_agrid
    REAL, INTENT(in) :: lat, lon
    
    INTEGER :: indxI, indxJ !< Indexes to be returned.
    INTEGER :: dimI, dimJ !< Size of the grid dimensions
    INTEGER :: i,j !< Count indexes
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
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="find_equator_index_agrid">
  !   <OVERVIEW>
  !     Return the closest index (i,j) to the given (lat,lon) point.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE FUNCTION find_equator_index_agrid(lat, lon)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     This function searches a equator grid tile looking for the grid point
  !     closest to the give (lat, lon) location, and returns the i
  !     and j indexes of the point.
  !   </DESCRIPTION>
  !   <IN NAME="lat" TYPE="REAL">
  !     Latitude location
  !   </IN>
  !   <IN NAME="lon" TYPE="REAL">
  !     Longitude location
  !   </IN>
  !   <OUT NAME="find_equator_index" TYPE="INTEGER, DIMENSION(2)">
  !     The (i, j) location of the closest grid to the given (lat,
  !     lon) location.
  !   </OUT>
  PURE FUNCTION find_equator_index_agrid(lat, lon)
    INTEGER, DIMENSION(2) :: find_equator_index_agrid
    REAL, INTENT(in) :: lat, lon
    
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
          IF (   (diag_global_grid%aglo_lon(i,0) <= lon .AND. lon <= 360) .OR.&
               & (0 <= lon .AND. lon < diag_global_grid%aglo_lon(i+1, 0)) ) THEN
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
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="latlon2xyz">
  !   <OVERVIEW>
  !     Return the (x,y,z) position of a given (lat,lon) point.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE ELEMENTAL TYPE(point) FUNCTION latlon2xyz(lat, lon)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Given a specific (lat, lon) point on the Earth, return the
  !     corresponding (x,y,z) location.  The return of latlon2xyz
  !     will be either a scalar or an array of the same size as lat
  !     and lon.
  !   </DESCRIPTION>
  !   <IN NAME="lat" TYPE="REAL">
  !     The latitude of the (x,y,z) location to find.  <TT>lat</TT>
  !     can be either a scalar or array.  <TT>lat</TT> must be of the
  !     same rank / size as <TT>lon</TT>.  This function assumes
  !     <TT>lat</TT> is in the range [-90,90].
  !   </IN>
  !   <IN NAME="lon" TYPE="REAL">
  !     The longitude of the (x,y,z) location to find.  <TT>lon</TT>
  !     can be either a scalar or array.  <TT>lon</TT> must be of the
  !     same rank / size as <TT>lat</TT>.  This function assumes
  !     <TT>lon</TT> is in the range [0,360].
  !   </IN>
  PURE ELEMENTAL TYPE(point) FUNCTION latlon2xyz(lat, lon)
    REAL, INTENT(in) :: lat, lon

    ! lat/lon angles in radians
    REAL :: theta, phi

    ! Convert the lat lon values to radians The lat values passed in
    ! are in the range [-90,90], but we need to have a radian range
    ! [0,pi], where 0 is at the north pole.  This is the reason for
    ! the subtraction from 90
    theta = deg2rad(90-lat)
    phi = deg2rad(lon)

    ! Calculate the x,y,z point
    latlon2xyz%x = RADIUS * SIN(theta) * COS(phi)
    latlon2xyz%y = RADIUS * SIN(theta) * SIN(phi)
    latlon2xyz%z = RADIUS * COS(theta)
  END FUNCTION latlon2xyz
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="distanceSqrd">
  !   <OVERVIEW>
  !     Find the distance between two points in the Cartesian
  !     coordinate space.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE ELEMENTAL REAL FUNCTION distanceSqrd(pt1, pt2)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>distanceSqrd</TT> will find the distance squared between
  !     two points in the xyz coordinate space.  <TT>pt1</TT> and <TT>
  !     pt2</TT> can either be both scalars, both arrays of the same
  !     size, or one a scalar and one an array.  The return value
  !     will be a scalar or array of the same size as the input array.
  !   </DESCRIPTION>
  !   <IN NAME="pt1" TYPE="TYPE(POINT)" />
  !   <IN NAME="pt2" TYPE="TYPE(POINT)" />
  PURE ELEMENTAL REAL FUNCTION distanceSqrd(pt1, pt2)
    TYPE(point), INTENT(in) :: pt1, pt2

    distanceSqrd = (pt1%x-pt2%x)**2 +&
         &         (pt1%y-pt2%y)**2 +&
         &         (pt1%z-pt2%z)**2
  END FUNCTION distanceSqrd
  ! </FUNCTION>
  ! </PRIVATE>

  ! <FUNCTION NAME="gCirDistance">
  !   <OVERVIEW>
  !     Find the distance, along the geodesic, between two points.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE ELEMENTAL REAL FUNCTION gCirDistance(lat1, lon1, lat2, lon2)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>aCirDistance</TT> will find the distance, along the geodesic, between two points defined by the (lat,lon) position of
  !     each point.
  !   </DESCRIPTION>
  !   <IN NAME="lat1" TYPE="REAL">Latitude of the first point</IN>
  !   <IN NAME="lon1" TYPE="REAL">Longitude of the first point</IN>
  !   <IN NAME="lat2" TYPE="REAL">Latitude of the second point</IN>
  !   <IN NAME="lon2" TYPE="REAL">Longitude of the second point</IN>
  PURE ELEMENTAL REAL FUNCTION gCirDistance(lat1, lon1, lat2, lon2)
    REAL, INTENT(in) :: lat1, lat2, lon1, lon2

    REAL :: theta1, theta2
    REAL :: deltaLambda !< Difference in longitude angles, in radians.
    REAL :: deltaTheta !< Difference in latitude angels, in radians.

    theta1 = deg2rad(lat1)
    theta2 = deg2rad(lat2)
    deltaLambda = deg2rad(lon2-lon1)
    deltaTheta = deg2rad(lat2-lat1)

    gCirDistance = RADIUS * 2 * ASIN(SQRT((SIN(deltaTheta/2))**2 + COS(theta1)*COS(theta2)*(SIN(deltaLambda/2))**2))
  END FUNCTION gCirDistance
  ! </FUNCTION>
END MODULE diag_grid_mod
