module exchange_mod

use utilities_mod, only: file_exist, open_file, check_nml_error, &
                         close_file
use mpp_mod

implicit none

private

!=========================================================================

! Module defines two accessible types; one that stores information about
! an exchange grid and one that provides a map from a component model grid
! to the cells of an exchange grid.
public exchange_map_type, boundary_map_type

! Subroutines and functions used for data exchange
public put_exchange_grid, get_exchange_grid

! Subroutines for getting information about exchange grid
public get_exchange_grid_size

! Subroutine that allows for dynamic partitioning of model subgrids
public set_frac_area

! Additional public interfaces available for initializing and setting up grids
public init_boundary_map, lon_lat_size, lon_lat_map
public complete_side1_boundary_map, complete_side2_boundary_map

! NO WRITING FACILITIES AVAILABLE FOR MPP AT THIS POINT, so no restarts
! Public interfaces for writing out and reading in maps
!public write_exchange_grid_map, write_boundary_map, read_exchange_grid_map, &
!   read_boundary_map

!==========================================================================

!---- namelist - with default values

logical :: make_exchange_reproduce = .true.
 
!      true => the summation of the gather operation done in get_exchange
!               is done in a magnitude sorted order so that results are
!               reproducible between runs and on any number of PEs. This can
!               be quite expensive, especially for small numbers of PEs and
!               can be turned off using the namelist for economy. For safety
!               the default is true.
!      false => the summation is not done in a reproducing fashion.

namelist /exchange_nml/ make_exchange_reproduce

!--- version number and tag name ----

character(len=128) :: version = '$Id: exchange.F90,v 1.3 2000/08/04 20:01:01 fms Exp $'
character(len=128) :: tag = '$Name: eugene $'


!==========================================================================

! Type that maps from model cells to exchange cells
type boundary_map_type
   private
   type(exchange_map_type), pointer :: ex
   integer :: side, id, total_target_cells
   logical :: allocated
   integer, pointer :: ex_cells_per_pe(:)
! New stuff in boundary_map for revision
   integer, pointer :: x(:), y(:), part(:), ex_cell(:)
   integer, pointer :: start(:), len(:), ex_start(:), ex_len(:)
   integer, pointer :: ex_ind(:)
end type boundary_map_type

! Back map from exchange grid to model cells
type exchange_map_type
   private
   integer :: max_size, size
! All second components are 2, one for each side
   real, pointer :: total_area(:), area(:), frac_area(:, :)
   integer, pointer :: bd_map_id(:, :), x(:, :), y(:, :), part(:, :), pe(:, :)
end type exchange_map_type
   
! Type to buffer data for transmits until all is sent
type transmit_buffer_type
   real, pointer :: data(:)
end type transmit_buffer_type

! VB: num_pes will come from above everywhere? How to set up buffer?
type (transmit_buffer_type), allocatable :: buffer(:)
logical :: first_lon_lat_size = .true.
logical :: read_namelist = .true.

!============================================================================

! All interfaces for exchanges are available for component model grids that
! are either 2d (2 space dimensions) or 3d (2 space dimensions plus model
! partions; i.e. ice and ice free in ice model).

! NOTE: May want to make this more efficient for two-d by adding extra code
interface put_exchange_grid
   module procedure put_exchange_grid_2d
   module procedure put_exchange_grid_3d
end interface

interface get_exchange_grid
   module procedure get_exchange_grid_2d
   module procedure get_exchange_grid_3d
end interface

! Grid overlap routines have been overloaded to allow efficient use by spectral
! models with grid decompositions that may not be latitudinally contiguous.
interface lon_lat_size
   module procedure lon_lat1_size
   module procedure lon_lat2_size
end interface

interface lon_lat_map
   module procedure lon_lat1_map
   module procedure lon_lat2_map
end interface


!============================================================================

! next_id is used as a registry for unique id's for boundary maps
integer :: next_id = 1

! Pi needed
real, parameter :: pi = 3.1415927

! Tolerance for exchange grid overlap
real, parameter :: tol = 0.0001

! Tolerance for sum of fractional areas being unity 
real, parameter :: area_tol = 0.000001

contains

!===========================================================================

subroutine init_boundary_map(bd_map, side, map)

!----------------------------------------------------------------------------
! Each boundary map must be initialized before use to set up storage.
! This routine also links them to the appropriate exchange map.
! Side specifies which side of the exchange grid the current boundary map
! is on. Side 1 maps are set up so that the exchange grid storage and the
! component model storage are guaranteed to be on the same pe, thus 
! eliminating the need for communication on puts and gets. Side 2 maps
! make no such guarantee and may require communication.
!----------------------------------------------------------------------------

implicit none

type (boundary_map_type), intent(inout) :: bd_map
type (exchange_map_type), intent(inout), target :: map
integer, intent(in) :: side
integer :: num_pes, unit, ierr, io

! Get unique boundary map id
bd_map%id = next_id
next_id = next_id + 1
! Set number of cells to 0 and current number of cells computed to 0
bd_map%total_target_cells = 0

! This map has not yet been set up
bd_map%allocated = .false.

! Set side for this component model on this exchange grid
if(side /= 1 .and. side /= 2) call fatal_error('init_boundary_map: side must be 1 or 2')
bd_map%side = side

! Initialize the size accumulation for the exchange grid
! Doesn't matter that this may be repeatedly set to 0
map%max_size = 0
! map%size -1 is a flag for lon_lat_map to call init_exchange_map
map%size = -1

! Link the boundary map to this exchange grid
bd_map%ex => map


!----------- Added March 2, 2000 -------------------
! Allocate storage for start and length of data for each remote pe
num_pes = mpp_npes()
allocate(bd_map%start(0:num_pes - 1), bd_map%len(0:num_pes - 1),  &
   bd_map%ex_start(0:num_pes - 1), bd_map%ex_len(0:num_pes - 1))
!----------------------------------------------------------------

!----------- Added 11 April, 2000 -------------------------
! Read namelist for run time control on first call here
if(read_namelist) then
   if ( file_exist( 'input.nml' ) ) then
      unit = open_file ( file = 'input.nml', action = 'read'  )
      ierr = 1
      do while ( ierr /= 0 )
         read ( unit,  nml = exchange_nml, iostat = io, end = 10 )
         ierr = check_nml_error ( io, 'land_model_nml' )
      enddo
10    continue
      call close_file( unit )
   endif

!  write the namelist to a log file
   unit = open_file ('logfile.out', action='append')
   if ( mpp_pe() == 0 ) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit, nml=exchange_nml)
   endif
   call close_file (unit)

! Name list has been read
   read_namelist = .false.
endif

end subroutine init_boundary_map

!===========================================================================

subroutine lon_lat1_size(blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map)

!---------------------------------------------------------------------------
! Computes the number of overlapping exchange grid cells formed by two lon
! lat grids. The lon lat grids are specified by their longitude and latitude
! boundaries which are passed in blon1, blat1 and blon2, blat2. Mask arrays,
! mask1 and mask2, indicate cells of the component model that are not used
! for exchange purposes (for instance cells in an ocean model that are over
! land). This routine works for contiguous lon/lat boundary specifications.
! This routine repacks the lon lat boundary arrays into two dimensional
! arrays of boundary pairs that allow non-continguous lon lat grids to be
! represented and calls lon_lat2_size to compute a variety of information
! about the grid overlap that is stored in the boundary map structure.
! Each component model grid has a number of partitions associated with it
! specified by num_part1, num_part2. 
!---------------------------------------------------------------------------

implicit none

real, intent(in) :: blon1(:), blat1(:), blon2(:), blat2(:)
logical, intent(in), dimension(:, :) :: mask1, mask2
integer, intent(in) :: num_part1, num_part2
type (boundary_map_type), intent(inout) :: bd_map

! Need paired lons/latitudes for call to lon_lat2_size
real :: lat_pair1(size(blat1) - 1, 2), lat_pair2(size(blat2) - 1, 2)
real :: lon_pair1(size(blon1) - 1, 2), lon_pair2(size(blon2) - 1, 2)
integer :: i

do i = 1, size(blon1) - 1
   lon_pair1(i, 1) = blon1(i)
   lon_pair1(i, 2) = blon1(i + 1)
end do

do i = 1, size(blon2) - 1
   lon_pair2(i, 1) = blon2(i)
   lon_pair2(i, 2) = blon2(i + 1)
end do

do i = 1, size(blat1) - 1
   lat_pair1(i, 1) = blat1(i)
   lat_pair1(i, 2) = blat1(i + 1)
end do

do i = 1, size(blat2) - 1
   lat_pair2(i, 1) = blat2(i)
   lat_pair2(i, 2) = blat2(i + 1)
end do

call lon_lat2_size(lon_pair1, lat_pair1, lon_pair2, lat_pair2, mask1, mask2, &
   num_part1, num_part2, bd_map)

end subroutine lon_lat1_size

!===========================================================================

subroutine lon_lat2_size(blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map)

! Does mpp version of lon_lat_size;

implicit none

real, intent(in) :: blon1(:, :), blon2(:, :)    ! Second dimension is 2
real, intent(in) :: blat1(:, :), blat2(:, :)    ! Second dimension is 2
logical, intent(in), dimension(:, :) :: mask1, mask2
integer, intent(in) :: num_part1, num_part2
type (boundary_map_type), intent(inout) :: bd_map

integer :: this_pe, num_pes
logical :: call_size = .true.

! Get number of pes and identity of this pe
this_pe = mpp_pe()
num_pes = mpp_npes()

! If this is first lon_lat_size call need to set up mpp buffers
if(first_lon_lat_size) then
   first_lon_lat_size = .false.
   allocate(buffer(0:num_pes-1))
endif


! Communication for setup is repeated rather than bufferd
call size_map(call_size, blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map, bd_map, this_pe, num_pes)

end subroutine lon_lat2_size

!===========================================================================
subroutine one_lon_lat_size(blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map)

!----------------------------------------------------------------------------
!
! SHOULD ONLY BE CALLED FOR SIDE 1 MODELS AS LEADING ARGUMENTS
!
! Given arrays of the longitude and latitude boundaries for two rectangular
! grids, and the number of partitions for the component models corresponding
! to each grid. The masks, mask1 and mask2, are
! set to false for component grid cells that should NOT be used in the
! computation.  This duplicates much of the computation in lon_lat_map
! but can be called first to allow for appropriate storage allocation in
! the map structure.  This subroutine can be called repeatedly to build
! up an exchange grid that connects several different component models 
! associated with the same exchange grid corresponding to the 
! exchange_map_type map.  
!  FOR NOW, ASSUME THAT THIS IS CALLED BY THE MODEL COMPONENT ON THE CURRENT
!  PE IN TURN FOR THE OTHER MODEL ON ALL PROCESSORS SEQUENTIALLY.
!---------------------------------------------------------------------------

implicit none

real, intent(in) :: blon1(:, :), blat1(:, :), blon2(:, :), blat2(:, :)
logical, intent(in), dimension(:, :) :: mask1, mask2
integer, intent(in) :: num_part1, num_part2 
type (boundary_map_type), intent(inout) :: bd_map

integer :: i1, j1, i2, j2, np
real :: lon1(2), lat1(2), lon2(2), lat2(2)

type (exchange_map_type), pointer :: map

integer :: lon_ov(size(blon1, 1) + 1 + size(blon2, 1) + 1, 2), lon_next, i

! ONLY CALLED FOR SIDE 1 check for error
if(bd_map%side /= 1) &
   call fatal_error('one_lon_lat_size; only call for side 1 bd_maps')

! Get map pointer
map => bd_map%ex

! Number of exchange grid cells needed to handle partitions for each overlap
np = num_part1 * num_part2

! Find all longitude overlap pairs
lon_next = 0
do i1 = 1, size(blon1, 1)
   lon1 = blon1(i1, 1:2)
   do i2 = 1, size(blon2, 1)
      lon2 = blon2(i2, 1:2)
      if(lon_overlap(lon1, lon2)) then
         lon_next = lon_next + 1
         if(lon_next > size(lon_ov)) &
            call fatal_error('one_lon_lat_size: lon_next too big')
         lon_ov(lon_next, 1) = i1
         lon_ov(lon_next, 2) = i2
       end if
   end do
end do

! Find all latitude overlap pairs
do j1 = 1, size(blat1, 1)
   lat1 = blat1(j1, 1:2)
   do j2 = 1, size(blat2, 1)
      lat2 = blat2(j2, 1:2)
      if(lat_overlap(lat1, lat2)) then
! Compute the latitudinal area of the exchange grid cell
         do i = 1, lon_next
            i1 = lon_ov(i, 1)
            i2 = lon_ov(i, 2)
            if(mask1(i1, j1) .and. mask2(i2, j2)) then 
! Increment space in exchange map and boundary map for side 1 boundary map calls
               map%max_size = map%max_size + np
               bd_map%total_target_cells = bd_map%total_target_cells + np
            endif
         end do
      end if
   end do
end do

end subroutine one_lon_lat_size

!===========================================================================

subroutine lon_lat1_map(blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map, bd_map2)

! Called to allocate storage for boundary map pointers
! lon_lat_size must be called for all bd_map pairs before repeating these
! calls for lon_lat_map

implicit none

real, intent(in) :: blon1(:), blat1(:), blon2(:), blat2(:)
logical, intent(in), dimension(:, :) :: mask1, mask2
integer, intent(in) :: num_part1, num_part2
type (boundary_map_type), intent(inout) :: bd_map
type (boundary_map_type), intent(in) :: bd_map2

! Need paired lons/latitudes for call to lon_lat2_map
real :: lat_pair1(size(blat1) - 1, 2), lat_pair2(size(blat2) - 1, 2)
real :: lon_pair1(size(blon1) - 1, 2), lon_pair2(size(blon2) - 1, 2)
integer :: i

do i = 1, size(blon1) - 1
   lon_pair1(i, 1) = blon1(i)
   lon_pair1(i, 2) = blon1(i + 1)
end do

do i = 1, size(blon2) - 1
   lon_pair2(i, 1) = blon2(i)
   lon_pair2(i, 2) = blon2(i + 1)
end do

do i = 1, size(blat1) - 1
   lat_pair1(i, 1) = blat1(i)
   lat_pair1(i, 2) = blat1(i + 1)
end do

do i = 1, size(blat2) - 1
   lat_pair2(i, 1) = blat2(i)
   lat_pair2(i, 2) = blat2(i + 1)
end do

call lon_lat2_map(lon_pair1, lat_pair1, lon_pair2, lat_pair2, &
   mask1, mask2, num_part1, num_part2, bd_map, bd_map2)

end subroutine lon_lat1_map

!===========================================================================

subroutine lon_lat2_map(blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map, bd_map2)

! Called to allocate storage for boundary map pointers
! lon_lat_size must be called for all bd_map pairs before repeating these
! calls for lon_lat_map

implicit none

! Second dimension on lon lats is 2
real, intent(in) :: blon1(:, :), blat1(:, :), blon2(:, :), blat2(:, :)
logical, intent(in), dimension(:, :) :: mask1, mask2
integer, intent(in) :: num_part1, num_part2
type (boundary_map_type), intent(inout) :: bd_map
type (boundary_map_type), intent(in) :: bd_map2

integer :: this_pe, num_pes
logical :: call_size = .false.

! Get number of pes and identity of this pe
this_pe = mpp_pe()
num_pes = mpp_npes()

call size_map(call_size, blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map, bd_map2, this_pe, num_pes)

end subroutine lon_lat2_map

!===========================================================================

subroutine one_lon_lat_map(blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map, bd_map2, this_pe, other_pe)

!---------------------------------------------------------------------------
! Given arrays of the longitude and latitude boundaries for two rectangular
! grids, and the number of partitions for the component models corresponding
! to each grid, generates the appropriate entries in an exchange_map
! and in each of the boundary maps for the two components models.
! The areas computed are percentage surface area of a unit sphere.  
! The masks, mask1 and mask2, are set to false for component grid cells 
! that should NOT be used in the computation.  Compare to lon_lat_size
! which just figures out storage needs.
!  FOR NOW, ASSUME THAT THIS IS CALLED BY THE MODEL COMPONENT ON THE CURRENT
!  PE IN TURN FOR THE OTHER MODEL ON ALL PROCESSORS SEQUENTIALLY.
!---------------------------------------------------------------------------

implicit none

real, intent(in) :: blon1(:, :), blat1(:, :), blon2(:, :), blat2(:, :)
logical, intent(in), dimension(:, :) :: mask1, mask2
integer, intent(in) :: num_part1, num_part2, this_pe, other_pe
type (boundary_map_type), intent(inout) :: bd_map
type (boundary_map_type), intent(in) :: bd_map2

real :: lon1(2), lat1(2), lon2(2), lat2(2), lon_area, lat_area, area, x
integer :: i1, j1, i2, j2, k, m, side, other_side
type (exchange_map_type), pointer :: map
real :: lon_ov_area(size(blon1, 1) + 1 + size(blon2, 1) + 1), tlon1(2), tlon2(2)
integer :: lon_ov(size(blon1, 1) + 1 + size(blon2, 1) + 1, 2), lon_next, i

! ONLY CALLED FOR SIDE 1 check for error
if(bd_map%side /= 1) &
   call fatal_error('one_lon_lat_map; only call for side 1 bd_maps')

! Get pointer to exchange map
map => bd_map%ex

side = 1
other_side = 2

! Make sure that the exchange grid map has been initialized
if(map%size == -1) call init_exchange_map(map)

! Make sure that the side 1 boundary map has storage initialized
if(.not. bd_map%allocated) call allocate_boundary_map(bd_map)

! Find all longitude overlap pairs
lon_next = 0
do i1 = 1, size(blon1, 1)
   lon1 = blon1(i1, 1:2)
   do i2 = 1, size(blon2, 1)
      lon2 = blon2(i2, 1:2)
      if(lon_overlap(lon1, lon2)) then
         lon_next = lon_next + 1
         if(lon_next > size(lon_ov)) &
            call fatal_error('one_lon_lat_map: lon_next too big')
         lon_ov(lon_next, 1) = i1
         lon_ov(lon_next, 2) = i2
! Compute the longitudinal area of the exchange grid cell
         x = min(lon1(1),lon2(1))
         tlon1 = lon1 -x
         tlon2 = lon2 -x
         if(tlon1(2) >= 2.0*pi) tlon1 = tlon1 - 2.0*pi
         if(tlon2(2) >= 2.0*pi) tlon2 = tlon2 - 2.0*pi
         lon_ov_area(lon_next) = min(tlon1(2),tlon2(2)) - max(tlon1(1),tlon2(1))
       end if
   end do
end do

! Find all latitude overlap pairs
do j1 = 1, size(blat1, 1)
   lat1 = blat1(j1, 1:2)
   do i = 1, lon_next
      i1 = lon_ov(i, 1)
      i2 = lon_ov(i, 2)
      do j2 = 1, size(blat2, 1)
         lat2 = blat2(j2, 1:2)
         if(lat_overlap(lat1, lat2)) then
! Compute the latitudinal area of the exchange grid cell
            lat_area =  sin(min(lat1(2), lat2(2))) -sin(max(lat1(1), lat2(1)))
            if(mask1(i1, j1) .and. mask2(i2, j2)) then
! Compute the area of the exchange grid cell
               area = lon_ov_area(i) * lat_area
! Loop through total number of partitions and load map
               do k = 1, num_part1
                  do m = 1, num_part2
! Increment number of exchange cells set up in this map
                     map%size = map%size + 1
! Violating the following error would indicate serious problems in setup
                     if(map%size > map%max_size) &
                        call fatal_error('one_lon_lat_map: Max map_size exceeded')
! Initialize area and total area for this exchange cell
                     map%total_area(map%size) = area
                     map%area(map%size) = area
! Initialized value for fractional areas is 1, okay for no partitions
                     map%frac_area(map%size, :) = 1.0
! Set model ids for the two bounding models
                     map%bd_map_id(map%size, side) = bd_map%id
                     map%bd_map_id(map%size, other_side) = bd_map2%id
! Initialize exchange grid back pointers
                     map%x(map%size, side) = i1
                     map%x(map%size, other_side) = i2
                     map%y(map%size, side) = j1
                     map%y(map%size, other_side) = j2
                     map%part(map%size, side) = k
                     map%part(map%size, other_side) = m
                     map%pe(map%size, side) = this_pe
                     map%pe(map%size, other_side) = other_pe
                  end do
               end do
            end if
         end if
      end do
   end do
end do

end subroutine one_lon_lat_map

!============================================================================

subroutine size_map(call_size, blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map, bd_map2, this_pe, num_pes)

! Computes overlaps between pairs of component model grids. Called by
! lon_lat_size with call_size true. In that case, the size of storage needed
! is computed but allocations for pointers are not done. Also called by
! lon_lat_map (call_size false) in which case pointer storage is 
! allocated.

implicit none

logical, intent(in) :: call_size
real, intent(in) :: blon1(:, :), blon2(:, :)    ! Second dimension is 2
real, intent(in) :: blat1(:, :), blat2(:, :)    ! Second dimension is 2
logical, intent(in), dimension(:, :) :: mask1, mask2
integer, intent(in) :: num_part1, num_part2, this_pe, num_pes
type (boundary_map_type), intent(inout) :: bd_map
type (boundary_map_type), intent(in) :: bd_map2

integer :: numq(2, 0:num_pes - 1), temp(2)
integer :: i, j, k, index, send_length, rec_length, from_pe, to_pe
real, pointer :: send_data(:), rec_data(:)
real, pointer :: rlon2(:, :), rlat2(:, :)
logical, pointer :: rmask2(:, :)

! Get number of lats/lons for each other processor's grid 2 information
numq(1, this_pe) = size(blon2, 1)
numq(2, this_pe) = size(blat2, 1)

if(size(mask2, 1)/=numq(1, this_pe) .or. size(mask2, 2)/=numq(2, this_pe)) &
   call fatal_error('size_map; mask2 size inconsistent with blon2 or blat2')

! Broadcast the lon and lat sizes
do from_pe = 0, num_pes - 1
   call mpp_transmit( numq(1:2,from_pe), 2, ALL_PES, numq(1:2,from_pe), 2, from_pe )
end do
      
! Load the lons and lats and the masks for transmission
! Send a pair of lons, a pair of lats, and a lon by lat mask
send_length = numq(1, this_pe) * 2 + numq(2, this_pe) * 2 + &
   numq(1, this_pe) * numq(2, this_pe)
! Allocate temporary space for these quantities
allocate(send_data(send_length))

! Load up the lower lon bounds, then the upper lon bounds
send_data(1:numq(1, this_pe)) = blon2(:, 1)
send_data(numq(1, this_pe) + 1 : 2 * numq(1, this_pe)) = blon2(:, 2)

! Load up the lower lat bounds, then the upper lat bounds
send_data(2*numq(1, this_pe) + 1 : 2*numq(1, this_pe) + numq(2, this_pe)) = &
   blat2(:, 1)
send_data(2*numq(1, this_pe) + numq(2, this_pe) + 1 : &
   2*numq(1, this_pe) + 2*numq(2, this_pe)) = blat2(:, 2)

! Loading the logical mask for sending as real; load as false, toggle if true
index = 2*numq(1, this_pe) + 2*numq(2, this_pe) + 1
send_data(index:) = 0.0
do i = 1, numq(1, this_pe)
   do j = 1, numq(2, this_pe)
      if(mask2(i, j)) send_data(index) = 1.0
      index = index + 1 
   end do
end do
! Development check
if(index /= send_length + 1) call fatal_error('size_map:inconsistent index')

! Broadcast the lons and lats and the masks
do from_pe = 0, num_pes - 1
! Compute length of data to be recieved from from_pe
   rec_length = 2*numq(1, from_pe) + 2*numq(2, from_pe) + &
        numq(1, from_pe) * numq(2, from_pe)
   if( this_pe.EQ.from_pe .AND. rec_length.NE.send_length )call fatal_error( 'size_map: rec_length.NE.send_length on this_pe.' )
! WARNING: THERE IS REAL CHANCE OF HEAP FRAGMENTATION HERE>>>
   allocate(rec_data(rec_length))
   allocate(rlon2(numq(1, from_pe),2), rlat2(numq(2, from_pe),2))
   allocate( rmask2(numq(1, from_pe), numq(2, from_pe)))
   call mpp_sync()
   call mpp_transmit( send_data, rec_length, ALL_PES, rec_data, rec_length, from_pe )
! Load the remote side2 pe's lon, lat and mask data
   rlon2(:, 1) = rec_data(1:numq(1, from_pe))
   rlon2(:, 2) = rec_data(numq(1, from_pe) + 1 : 2*numq(1, from_pe))
   rlat2(:, 1) = rec_data(2*numq(1, from_pe) + 1 : 2*numq(1, from_pe) + &
        numq(2, from_pe))
   rlat2(:, 2) = rec_data(2*numq(1, from_pe) + numq(2, from_pe) + 1 : &
        2*numq(1, from_pe) + 2*numq(2, from_pe))
! Load the logical mask array; set to .false. toggle if true
   rmask2 = .false.
   index = 2*numq(1, from_pe) + 2*numq(2, from_pe) + 1
   do j = 1, numq(1, from_pe)
      do k = 1, numq(2, from_pe)
         if(rec_data(index) == 1.0) rmask2(j, k) = .true.
         index = index + 1
      end do
   end do
   if(index/=rec_length + 1) call fatal_error('size_map:inconsistnt index')

! Compute the overlaps
   if(call_size) then
       call one_lon_lat_size(blon1, blat1, rlon2, rlat2, mask1, rmask2, &
            num_part1, num_part2, bd_map)
   else
       call one_lon_lat_map(blon1, blat1, rlon2, rlat2, mask1, rmask2, &
            num_part1, num_part2, bd_map, bd_map2, this_pe, from_pe)
   endif
   deallocate(rec_data, rlon2, rlat2, rmask2)
end do

! Free up the broadcast buffers
deallocate(send_data)

end subroutine size_map

!============================================================================

function lat_overlap(lat1, lat2)

!----------------------------------------------------------------------------
! Determines if two latitude bands overlap.  A tolerance value is included 
! to take care of round-off errors in boundary computations.  
!----------------------------------------------------------------------------

implicit none

real, intent(in), dimension(2) :: lat1, lat2
logical :: lat_overlap

if(lat1(1) >= lat2(2) - tol .or. lat1(2) <= lat2(1) + tol) then
   lat_overlap = .false.
else
   lat_overlap = .true.
end if

end function lat_overlap

!============================================================================

function lon_overlap(ilon1, ilon2)

!----------------------------------------------------------------------------
! Determines if two longitude bands overlap.  A tolerance value is included 
! to take care of round-off errors in boundary computations.  The incoming 
! longitude values are assumed to be in the range 0 to 4Pi
!----------------------------------------------------------------------------

implicit none

real, intent(in), dimension(2) :: ilon1, ilon2
logical :: lon_overlap

real, dimension(2) :: lon1, lon2
real :: x

! Copy inputs to allow modification for wraparound
lon1 = ilon1;     lon2 = ilon2

! Shift the intervals so that one starts at 0
x = min(lon1(1),lon2(1))
lon1 = lon1 -x
lon2 = lon2 -x
! If the other one is still exceeding 2pi, move it down
if(lon1(2) >= 2.0*pi) lon1 = lon1 - 2.0*pi
if(lon2(2) >= 2.0*pi) lon2 = lon2 - 2.0*pi

!Check for overlap with rounding tolerance
if(lon1(1) >= lon2(2) - tol .or. lon1(2) <= lon2(1) + tol) then
   lon_overlap = .false.
else
   lon_overlap = .true.
end if

end function lon_overlap

!===========================================================================

subroutine init_exchange_map(map)

!---------------------------------------------------------------------------
! Initializes an exchange grid map.
!---------------------------------------------------------------------------

implicit none

type(exchange_map_type), intent(inout), target :: map
integer :: num

! Number of exchange grid cells has been accumulated
num = map%max_size

! There are num cells
allocate(map%total_area(num), map%area(num), map%frac_area(num, 2), &
   map%bd_map_id(num, 2), map%x(num, 2), map%y(num, 2), &
   map%part(num, 2), map%pe(num, 2))

! max_size stores total number of cells, size has how many currently filled
map%size = 0

end subroutine init_exchange_map

!===========================================================================

function get_exchange_grid_size(bd_map)

!---------------------------------------------------------------------------
! Gets size of total exchange grid for all boundary maps
! A single boundary map is passed as argument but only to allow a pointer
! to the exchange map to be accessed.
!---------------------------------------------------------------------------

implicit none

type (boundary_map_type), intent(in) :: bd_map
integer :: get_exchange_grid_size

get_exchange_grid_size = bd_map%ex%max_size

end function get_exchange_grid_size

!===========================================================================

subroutine allocate_boundary_map(bd_map)

!----------------------------------------------------------------------------
! Initializes a boundary map, bd_map, for a component model grid
!----------------------------------------------------------------------------

implicit none

type(boundary_map_type), intent(inout) :: bd_map

! Allocations for revised algorithm added 29 Feb., 2000 --------------------
allocate(bd_map%x(bd_map%total_target_cells), &
   bd_map%y(bd_map%total_target_cells),bd_map%part(bd_map%total_target_cells), &
   bd_map%ex_cell(bd_map%total_target_cells))
!---------------------------------------------------------------------------

! Set allocation flag
bd_map%allocated = .true.

end subroutine allocate_boundary_map

!===========================================================================

subroutine complete_side1_boundary_map(bd_map)

! Called for side 1 boundary maps after all side 1 calls to lon_lat_map
! have been completed. 

implicit none

type(boundary_map_type), intent(inout) :: bd_map
integer :: count, i, this_pe, side, num_pes
type (exchange_map_type), pointer :: map

if(bd_map%side /= 1) &
   call fatal_error('complete_side1_boundary_map: called for side2')

! Get map pointer
map => bd_map%ex

! CODE SECTION ADDED FOR NEW IMPLEMENTATION, 29 FEB., 2000 ----------------
do i = 1, bd_map%total_target_cells
   bd_map%x(i) = map%x(i, 1)
   bd_map%y(i) = map%y(i, 1)
   bd_map%part(i) = map%part(i, 1)
   bd_map%ex_cell(i) = i
end do

! On side 1 pe, target exchange cells are ALL on this pe
this_pe = mpp_pe()
bd_map%start = 0
bd_map%len = 0
bd_map%start(this_pe) = 1
bd_map%len(this_pe) = bd_map%total_target_cells

! Loop through each exchange grid cell; find out what pe has the
! target component model cell for this side of the exchange grid,
! and increment a counter of how many cells on each pe.
side = 1

num_pes = mpp_npes()
allocate(bd_map%ex_cells_per_pe(0:num_pes - 1))
bd_map%ex_cells_per_pe = 0
do i = 1, map%size
   if(map%bd_map_id(i, side) == bd_map%id) &
   bd_map%ex_cells_per_pe(map%pe(i, side)) = bd_map%ex_cells_per_pe(map%pe(i, side)) + 1
end do

end subroutine complete_side1_boundary_map

!===========================================================================

subroutine complete_side2_boundary_map(bd_map, num_x, num_y, num_parts)

implicit none

! Called for side 2 boundary maps after all side 1 calls to lon_lat_map
! have been completed. One cannot guarantee that side 2 component model
! cells and the corresponding exchange grid cells are on the same pe.
! Need to send around the information to all needed pes.

integer, intent(in) :: num_x, num_y, num_parts
type (boundary_map_type), intent(inout) :: bd_map

type (exchange_map_type), pointer :: ex_map

integer, allocatable :: send_count(:), receive_count(:)
integer, allocatable :: send_x(:), send_y(:), send_part(:), &
   send_ex_cell(:)
integer, allocatable :: x(:), y(:), part(:), ex_cell(:), ex_pe(:)
integer :: this_pe, num_pes, pe_dist
integer :: total_target_cells, pe, index, start, end, i, j, k, from_pe, to_pe
integer :: ipe, total_ex_cells, side
integer, allocatable :: pe_index(:) 
real :: temp

! Get number of pes and identity of this pe
this_pe = mpp_pe()
num_pes = mpp_npes()
allocate(send_count(0:num_pes - 1), receive_count(0:num_pes - 1), pe_index(0:num_pes - 1))

! Sets up the boundary map pointers for side 2 boundary maps
if(bd_map%side /= 2) &
   call fatal_error('complete_side2_boundary_map: called for side1')

! Get map pointer
ex_map => bd_map%ex

! Allocate storage space for this side 2 boundary map
call allocate_boundary_map(bd_map)


! Initialize the number of exchange grid cells per pe
! Get map pointer
side = 2

! Loop through each exchange grid cell; find out what pe has the
! target component model cell for this side of the exchange grid,
! and increment a counter of how many cells on each pe.
side = 2

allocate(bd_map%ex_cells_per_pe(0:num_pes - 1))
bd_map%ex_cells_per_pe = 0
do i = 1, ex_map%size
   if(ex_map%bd_map_id(i, side) == bd_map%id) &
   bd_map%ex_cells_per_pe(ex_map%pe(i, side)) = bd_map%ex_cells_per_pe(ex_map%pe(i, side)) + 1
end do

! Get number of cells this pe's exchange grid sends to all other pes
send_count = ex_cells_per_pe(bd_map, num_pes)

! Just copy for info on the same pe
receive_count(this_pe) = send_count(this_pe)

! Loop through to send info from each pe to all other pes in turn
do pe_dist = 1,num_pes-1
   to_pe = mod( this_pe+pe_dist, num_pes )
   from_pe = mod( this_pe+num_pes-pe_dist, num_pes )
   call mpp_transmit( send_count(to_pe), 1, to_pe, receive_count(from_pe), 1, from_pe )
end do

! Compute the total number of target cells for this boundary map 
total_target_cells = sum(receive_count)
bd_map%total_target_cells = total_target_cells
! Allocate buffers plus temporary storage for sort
allocate(ex_pe(total_target_cells), x(total_target_cells), &
   y(total_target_cells), ex_cell(total_target_cells), part(total_target_cells))

! Storage allocation
allocate(bd_map%x(bd_map%total_target_cells), &
   bd_map%y(bd_map%total_target_cells),bd_map%part(bd_map%total_target_cells), &
   bd_map%ex_cell(bd_map%total_target_cells))

! Load up information on where the target exchange cells hang out
do i = 0, num_pes - 1
   bd_map%len(i) = receive_count(i)
end do
bd_map%start(0) = 1
do i = 1, num_pes - 1
   bd_map%start(i) = bd_map%start(i - 1) + bd_map%len(i - 1)
end do

! Initial index values for gathering the transmissions from other pes
start = 1
end = 0

! Send to each pe the number a quadruplet of boundary grid x, y, part and
! corresponding exchange grid cell number
! Loop through to send data to each processor in turn
! Currently sending to self also , may want to optimize this
! pe loop states that this pe is transmitting to all others
do from_pe = 0, num_pes - 1

! Send to all pes if it's my_pe's turn
   if(from_pe == this_pe) then
      do to_pe = 0, num_pes - 1
         allocate(send_x(send_count(to_pe)), send_y(send_count(to_pe)), &
            send_part(send_count(to_pe)), send_ex_cell(send_count(to_pe)))
         index = 1
! Load up the tranmsit buffers for the current target pe
         do i = 1, ex_map%size
            if(ex_map%pe(i, bd_map%side) == to_pe .and. &
               ex_map%bd_map_id(i, bd_map%side) == bd_map%id) then
               send_x(index) = ex_map%x(i, bd_map%side)
               send_y(index) = ex_map%y(i, bd_map%side)
               send_part(index) = ex_map%part(i, bd_map%side)
               send_ex_cell(index) = i
               index = index + 1
            end if
         end do
! Development error check on index here
         if(index /= send_count(to_pe) + 1) call fatal_error &
            ('complete_side2_boundary_map: index and send_count inconsistent')

! THIS MUST AGGREGATE AND BUFFER
         call complete_bd_transmit(to_pe, send_x, send_y, send_part, &
            send_ex_cell)
         deallocate(send_x, send_y, send_part, send_ex_cell)
      end do
   end if

! Every pe is receiving a message from this one (self-transmit for now)
   start = end + 1
   end = end + receive_count(from_pe)
   call complete_bd_receive(from_pe, x(start:end), y(start:end), &
      part(start:end), ex_cell(start:end))
   call mpp_sync
! Take care of pe number
   ex_pe(start:end) = from_pe
end do

! Development error check
if(end /= total_target_cells) &
   call fatal_error('complete_boundary_map: end /= to total target cells')

! Free up the communication buffering
call end_exchange(num_pes)

! Need to grab stuff in order by pe; NO NEED TO SORT BY OTHER KEYS???
pe_index = 0
pe_index(0 : num_pes - 1) = bd_map%start(0 : num_pes - 1)

do i = 1, total_target_cells
   ipe = ex_pe(i)
   index = pe_index(ipe)
   pe_index(ipe) = pe_index(ipe) + 1
   bd_map%x(index) = x(i)
   bd_map%y(index) = y(i)
   bd_map%part(index) = part(i)
   bd_map%ex_cell(index) = ex_cell(i)
end do

! Development check, make sure indices aren't wrapping
do i = 0, num_pes - 1
   if(pe_index(i) .ne. bd_map%start(i) + bd_map%len(i)) then
      write(*, *) 'FAILED in pe_index check'
      write(*, *) 'FAIL:', i, pe_index(i), bd_map%start(i), bd_map%len(i)
      stop
   endif
end do

! Block for enhancing get efficiency, index to exchange
! First, figure out how many exchange cells on this pe go each remote pe
send_count = ex_cells_per_pe(bd_map, num_pes)
bd_map%ex_start(0) = 1
do i = 1, num_pes - 1
   bd_map%ex_start(i) = bd_map%ex_start(i - 1) + send_count(i - 1)
end do   
do i = 0, num_pes - 1
   bd_map%ex_len(i) = send_count(i)
end do
total_ex_cells = sum(send_count(0:num_pes - 1))
allocate(bd_map%ex_ind(total_ex_cells))

! Loop through the exchange grid to find all the cells that need to go off PE
do to_pe = 0, num_pes - 1
   index = bd_map%ex_start(to_pe)
   do i = 1, ex_map%size
      if(ex_map%pe(i, bd_map%side) == to_pe .and. &
         ex_map%bd_map_id(i, bd_map%side) == bd_map%id) then
         bd_map%ex_ind(index) = i
         index = index + 1
      endif
   end do
end do

deallocate(send_count, receive_count)

end subroutine complete_side2_boundary_map

!===========================================================================

subroutine get_exchange_grid_2d(ex_data, model_data, bd_map)

! Gets data from exchange grid for a component model with no partitions

implicit none

real, intent(inout) :: model_data(:, :)
type (boundary_map_type), intent(inout) :: bd_map
real, intent(in) :: ex_data(:)

real :: model_data3(size(model_data, 1), size(model_data, 2), 1)

call get_exchange_grid_3d(ex_data, model_data3, bd_map)

model_data(:, :) = model_data3(:, :, 1)

end subroutine get_exchange_grid_2d

!===========================================================================

subroutine get_exchange_grid_3d(ex_data, model_data, bd_map)

! Gets data from exchange grid for a given component model whose boundary
! map is specified.

implicit none

real, intent(inout) :: model_data(:, :, :)
type (boundary_map_type), intent(inout) :: bd_map
real, intent(in) :: ex_data(:)

integer :: this_pe, num_pes
type (exchange_map_type), pointer :: ex_map

integer, allocatable :: send_count(:), receive_count(:)
real, allocatable :: send_data(:), send_area(:), rec_data(:), rec_area(:)
integer, allocatable :: send_to_x(:), send_to_y(:), send_to_part(:)
integer, allocatable :: x(:), y(:), part(:)
integer :: i, j, k, index, cell, m, ind_num, side, from_pe, to_pe
integer :: rec_total, start, end
integer :: send_count_max, recv_count_max, send_length, recv_length, pe_dist
real, allocatable, dimension(:) :: send_buffer, recv_buffer

! Get map pointer
ex_map => bd_map%ex

! Do get for  side1 pes
if(bd_map%side == 1) then
   call get_sum(bd_map%x,bd_map%y,bd_map%part, ex_data, ex_map%area, model_data)
   return
end if

! Rest of code is for side 2 pes which may require communication
! Get number of pes and identity of this pe
this_pe = mpp_pe()
num_pes = mpp_npes()
allocate(send_count(0:num_pes - 1), receive_count(0:num_pes - 1))

! First determine how much data is going to each processor from exchange
! grid cells that correspond to this bd_map. I.E. how many exchange grid
! values will my pe receive from each other pe for this get.
receive_count = bd_map%len
rec_total = sum(receive_count(0:num_pes - 1))
recv_count_max = maxval(receive_count)
! Allocate the buffers for receiving data
allocate(rec_data(rec_total), rec_area(rec_total), x(rec_total), y(rec_total), &
   part(rec_total))
allocate( recv_buffer(recv_count_max*5) )

! Which side of exchange grid is this component model on?
side = bd_map%side

! Next determine how much data is going to each processor?  I.E. how many
! exchange grid values will my processor send to every other processor
! for this get, (where are the model cells corresponding to my pes
! exchange grid cells).
send_count = ex_cells_per_pe(bd_map, num_pes)
send_count_max = maxval(send_count)
allocate( send_buffer(send_count_max*5) )

! Can I do this, or should old values be kept???
model_data = 0.0

! Set up indices for reading into the receive buffers
start = 1
end = 0
! ringwise communication
do pe_dist = 0, num_pes-1
   to_pe = mod( this_pe+pe_dist, num_pes )
   from_pe = mod( this_pe+num_pes-pe_dist, num_pes )
   if(send_count(to_pe) /= 0) then
       call mpp_sync_self()
       do i = 1, bd_map%ex_len(to_pe)
          index = bd_map%ex_ind(bd_map%ex_start(to_pe) + i - 1)
          send_buffer(5*i - 4) = ex_data(index)
          send_buffer(5*i - 3) = ex_map%area(index)
          send_buffer(5*i - 2) = ex_map%x(index, side)
          send_buffer(5*i - 1) = ex_map%y(index, side)
          send_buffer(5*i)     = ex_map%part(index, side)
       end do 
       send_length = send_count(to_pe)*5
   else
       to_pe = NULL_PE
       send_length = 1
   end if
   if(receive_count(from_pe) /= 0) then
       recv_length = receive_count(from_pe)*5
   else
       recv_length = 1
       from_pe = NULL_PE
   end if
   call mpp_transmit( send_buffer, send_length, to_pe, recv_buffer, recv_length, from_pe )
   if( from_pe.NE.NULL_PE )then
! Update pointers into receive buffer
       start = bd_map%start(from_pe)
       end = start + bd_map%len(from_pe) - 1
       rec_data(start:end) = recv_buffer(1:recv_length:5)
       rec_area(start:end) = recv_buffer(2:recv_length:5)
       x(start:end) = recv_buffer(3:recv_length:5)
       y(start:end) = recv_buffer(4:recv_length:5)
       part(start:end) = recv_buffer(5:recv_length:5)
   end if
end do
call mpp_sync_self()
deallocate( send_buffer, recv_buffer )

! Do weighted sum of the incoming data into the appropriate target cells
call get_sum(x, y, part, rec_data, rec_area, model_data)

! Free up the local buffering
deallocate(rec_data, rec_area, x, y, part, send_count, receive_count)

end subroutine get_exchange_grid_3d

!===========================================================================

subroutine get_sum(x, y, part, dat, area, sum)

implicit none

integer, intent(in) :: x(:), y(:), part(:)
real, intent(in) :: dat(:), area(:)
real, intent(out) :: sum(:, :, :)

integer :: nx, ny, npart, i, j, k
real :: total_area(size(sum, 1), size(sum, 2), size(sum, 3))

integer :: cx(size(x)), cy(size(y)), cpart(size(part))
real :: cdat(size(dat)), carea(size(area))
logical :: problem( size(sum, 1), size(sum, 2), size(sum, 3) )

problem=.FALSE.
nx = size(sum, 1)
ny = size(sum, 2)
npart = size(sum, 3)

sum = 0.0
total_area = 0.0

! Want quick efficient sum if no sort
if(.not. make_exchange_reproduce) then
   do i = 1, size(x)
      sum(x(i), y(i), part(i)) = sum(x(i), y(i), part(i)) + dat(i) * area(i)   
      total_area(x(i), y(i), part(i)) = total_area(x(i), y(i), part(i)) +area(i)
   end do
   
!   do i = 1, nx
!      do j = 1, ny
!         do k = 1, npart
!            if(total_area(i, j, k) == 0.0) then
!               if(sum(i, j, k) /= 0.0) &
!                  call fatal_error('get_sum: non-zero area')
!            else
!               sum(i, j, k) = sum(i, j, k) / total_area(i, j, k)
!            end if
!         end do
!      end do
!   end do
   problem = .FALSE.
   where( total_area.NE.0. )
       sum = sum/total_area
   elsewhere
       problem = sum.NE.0.
   end where
   if( any(problem) )call fatal_error('get_sum: non-zero area')
   return
end if

! NOTICE : ORDER ON THESE MUST BE FIXED TO REPRODUCE; SORT REQUIRED
! THIS IS POTENTIALLY EXTREMELY COSTLY: MAY WANT TO TURN IT OFF
cdat = dat * area
carea = area
cx = x
cy = y
cpart = part
call get_exchange_sort(cdat, carea, cx, cy, cpart, size(sum, 1), &
   size(sum, 2), size(sum, 3))

do i = 1, size(x)
   sum(cx(i), cy(i), cpart(i)) = sum(cx(i), cy(i), cpart(i)) + cdat(i)
   total_area(cx(i), cy(i), cpart(i)) = total_area(cx(i), cy(i), cpart(i)) + carea(i)
end do

!do i = 1, nx
!   do j = 1, ny
!      do k = 1, npart
!         if(total_area(i, j, k) == 0.0) then
!            if(sum(i, j, k) /= 0.0) &
!               call fatal_error('get_sum: non-zero area')
!         else
!            sum(i, j, k) = sum(i, j, k) / total_area(i, j, k)
!         end if
!      end do
!   end do
!end do
problem = .FALSE.
where( total_area.NE.0. )
    sum = sum/total_area
elsewhere
    problem = sum.NE.0.
end where
if( any(problem) )call fatal_error('get_sum: non-zero area')

end subroutine get_sum

!===========================================================================

subroutine get_exchange_sort(dat, area, x, y, part, num_x, num_y, num_part)

! Sorts all 5 fields by dat value from low to high
! Sort is needed to guarantee reproducible get_exchange_grid sums

implicit none

real, intent(inout) :: dat(:), area(:)
integer, intent(in) :: num_x, num_y, num_part
integer, intent(inout) :: x(:), y(:), part(:)

integer :: num(num_x, num_y, num_part), start(num_x, num_y, num_part)
real :: tdat(size(dat)), tarea(size(dat))
integer :: tx(size(dat)), ty(size(dat)), tpart(size(part))
integer :: i, j, k, next, ind, st, en

! Get count of how many values to sort for each cell
num(:, :, :) = 0
do i = 1, size(dat)
   num(x(i), y(i), part(i)) = num(x(i), y(i), part(i)) + 1
end do

! Next flat storage location to be used for this model_cell
next = 1
do i = 1, num_x
   do j = 1, num_y
      do k = 1, num_part
         start(i, j, k) = next
         next = next + num(i, j, k)
      end do
   end do
end do
! Development check, can be removed
if(next /= size(dat) + 1) then
   write(*, *) 'ERROR, next doesnt work'
   stop
end if

! Could do this in a more storage efficient fashion if needed
do i = 1, size(dat)
   ind = start(x(i), y(i), part(i))
   tdat(ind) = dat(i)
   tarea(ind) = area(i)
   tx(ind) = x(i)
   ty(ind) = y(i)
   tpart(ind) = part(i)
   start(x(i), y(i), part(i)) = start(x(i), y(i), part(i)) + 1
end do

! Now sort only within each cell
do i = 1, num_x
   do j = 1, num_y
      do k = 1, num_part
! Compute start and end for this cell; remember start has been incremented
         en = start(i, j, k) - 1
         st = start(i, j, k) - num(i, j, k)
         if(en - st > 1) call single_cell_sort(tdat(st:en), tarea(st:en), &
            tx(st:en), ty(st:en), tpart(st:en))
      end do
   end do
end do

! Copy sorted arrays back to inputs
dat = tdat
area = tarea
x = tx
y = ty
part = tpart

end subroutine get_exchange_sort

!===========================================================================

subroutine single_cell_sort(dat, area, x, y, part)

! Sorts all 5 fields by dat value from low to high
! Sort is needed to guarantee reproducible get_exchange_grid sums

implicit none

real, intent(inout) :: dat(:), area(:)
integer, intent(inout) :: x(:), y(:), part(:)

integer :: i, j, ti
real :: tr

! For now do silly bubble sort, make more efficient as needed

do i = 1, size(dat) - 1
   do j = i + 1, size(dat)
      if(dat(i) > dat(j)) then
! Switch dat, area, x, y, part
         tr = dat(i)
         dat(i) = dat(j)
         dat(j) = tr

         tr = area(i)
         area(i) = area(j)
         area(j) = tr

         ti = x(i)
         x(i) = x(j)
         x(j) = ti

         ti = y(i)
         y(i) = y(j)
         y(j) = ti

         ti = part(i)
         part(i) = part(j)
         part(j) = ti
      end if
   end do
end do

end subroutine single_cell_sort

!===========================================================================

subroutine put_exchange_grid_2d(model_data, ex_data, bd_map, avail)

! Puts data from a component model without partitions onto the exchange grid

implicit none

real, intent(in) :: model_data(:, :)
type (boundary_map_type), intent(inout) :: bd_map
real, intent(inout) :: ex_data(:)
logical, optional, intent(inout) :: avail(:)

real :: model_data3(size(model_data, 1), size(model_data, 2), 1)

model_data3(:, :, 1) = model_data(:, :)

call put_exchange_grid_3d(model_data3, ex_data, bd_map, avail)

end subroutine put_exchange_grid_2d

!===========================================================================

subroutine put_exchange_grid_3d(model_data, ex_data, bd_map, avail)

! Puts data from a component model grid onto the exchange grid

implicit none

real, intent(in) :: model_data(:, :, :)
type (boundary_map_type), intent(inout) :: bd_map
real, intent(inout) :: ex_data(:)
logical, optional, intent(inout) :: avail(:)

type (exchange_map_type), pointer :: ex_map

integer, allocatable :: send_count(:), receive_count(:)
integer :: num_pes, this_pe
real, allocatable :: send_data(:), rec_data(:)
integer, allocatable :: send_to_cell(:), cell(:)
integer :: i, j, k, index, m, ind_num, side, from_pe, to_pe
integer :: send_count_max, recv_count_max, send_length, recv_length, pe_dist
real, allocatable, dimension(:) :: send_buffer, recv_buffer

! Get map pointer
ex_map => bd_map%ex

! Get number of pes and identity of this pe
this_pe = mpp_pe()
num_pes = mpp_npes()

! Block for side1 
if(bd_map%side == 1) then
   do i = 1, bd_map%total_target_cells
      ex_data(i) = model_data(bd_map%x(i), bd_map%y(i), bd_map%part(i))
   end do
! If avail is present, need to set it to true if area is not 0
   if(present(avail)) then
      do i = 1, bd_map%total_target_cells
         if(ex_map%area(i) /= 0.0) avail(i) = .true.
      end do
   endif
   return
end if

! Rest of subroutine is for side 2, may require communication
! Allocate buffer storage
allocate(send_count(0:num_pes - 1), receive_count(0:num_pes - 1))

! First determine how much data is going to each processor from this one?
send_count = bd_map%len
send_count_max = maxval(send_count)
! Allocate storage to do the gather of data to be sent to this pe
allocate( send_buffer(send_count_max*2) )

! Next determine how much data is coming from each processor
side = bd_map%side
receive_count = ex_cells_per_pe(bd_map, num_pes)
recv_count_max = maxval(receive_count)
allocate(rec_data(recv_count_max), cell(recv_count_max), recv_buffer(recv_count_max*2))

! ringwise communication
do pe_dist = 0, num_pes-1
   to_pe = mod( this_pe+pe_dist, num_pes )
   from_pe = mod( this_pe+num_pes-pe_dist, num_pes )
   if(send_count(to_pe) /= 0) then
      call mpp_sync_self()
      index = bd_map%start(to_pe) - 1
      do i = 1, send_count(to_pe)
         send_buffer(2*i-1) = model_data(bd_map%x(index + i), &
            bd_map%y(index + i), bd_map%part(index + i))
         send_buffer(2*i) = bd_map%ex_cell(index + i)
      end do
      send_length = send_count(to_pe) * 2

   else
       to_pe = NULL_PE
       send_length = 1
   end if
   if(receive_count(from_pe) /= 0) then
       recv_length = receive_count(from_pe)*2
   else
       recv_length = 1
       from_pe = NULL_PE
   end if
   call mpp_transmit( send_buffer, send_length, to_pe, recv_buffer, recv_length, from_pe )
   if( from_pe.NE.NULL_PE )then
       rec_data(1:recv_length/2) = recv_buffer(1:recv_length:2)
       cell(1:recv_length/2) = recv_buffer(2:recv_length:2)
! Put the data in the appropriate exchange cell and set avail to true
      do i = 1, recv_length/2
         ex_data(cell(i)) = rec_data(i)
      end do
! Set avail if there is data and exchange cell has non-zero area
      if(present(avail)) then
         do i = 1, recv_length/2
            if(ex_map%area(cell(i)) /= 0.0) avail(cell(i)) = .true.
         end do
      end if
   end if
end do
call mpp_sync_self()
deallocate( send_buffer, recv_buffer, rec_data, cell, send_count,receive_count )

end subroutine put_exchange_grid_3d

!===========================================================================

subroutine set_frac_area(frac_area, bd_map)

! WARNING: MAY NOT HAVE BEEN ROBUSTLY TESTED

! Sets the fractional area of the different partitions associated with a
! particular boundary map. The exchange grid data in the corresponding
! exchange map must be updated to reflect this area change. The logic
! here is almost identical to put_exchange_grid.

implicit none

real, intent(in) :: frac_area(:, :, :)
type (boundary_map_type), intent(inout) :: bd_map

integer :: num_pes, this_pe
type (exchange_map_type), pointer :: ex_map
integer :: i, j, k, m, ind_num, side, index, from_pe, to_pe
integer, allocatable :: send_count(:), receive_count(:)

real, allocatable :: send_frac(:)
integer, allocatable :: send_to_cell(:), cell(:)
real, allocatable :: rec_frac(:)
integer :: send_count_max, recv_count_max, send_length, recv_length, pe_dist
real, allocatable, dimension(:) :: send_buffer, recv_buffer

! Get pointer to corresponding exchange grid map
ex_map => bd_map%ex

! Block for side1
if(bd_map%side == 1) then
   do i = 1, bd_map%total_target_cells
      ex_map%frac_area(i, 1) = frac_area(bd_map%x(i), bd_map%y(i), bd_map%part(i))
      ex_map%area(i) = ex_map%total_area(i) * ex_map%frac_area(i, 1) * &
         ex_map%frac_area(i, 2)
   end do
   return
end if

! Get number of pes and identity of this pe
this_pe = mpp_pe()
num_pes = mpp_npes()
allocate(send_count(0:num_pes - 1), receive_count(0:num_pes - 1))

! Error check for fractional area of cell not adding up to 1.0
if(any(abs(sum(frac_area, 3) - 1.0) > area_tol)) &
   call fatal_error('set_frac_area; sum of fractions is not 1')

! First determine how much data is going to each processor from this one
send_count = bd_map%len
send_count_max = maxval(send_count)
! Allocate storage to do the gather of data to be sent to this pe
allocate( send_buffer(send_count_max*2) )

! Next determine how much data is coming from each processor to this one.
side = bd_map%side
receive_count = ex_cells_per_pe(bd_map, num_pes)
recv_count_max = maxval(receive_count)
allocate( recv_buffer(recv_count_max*2), rec_frac(recv_count_max), cell(recv_count_max))

do pe_dist = 0, num_pes-1
   to_pe = mod( this_pe+pe_dist, num_pes )
   from_pe = mod( this_pe+num_pes-pe_dist, num_pes )
   if(send_count(to_pe) /= 0) then
      call mpp_sync_self()
      index = bd_map%start(to_pe) - 1
      do i = 1, send_count(to_pe)
         send_buffer(2*i-1) = frac_area(bd_map%x(index + i), &
            bd_map%y(index + i), bd_map%part(index + i))
         send_buffer(2*i) = bd_map%ex_cell(index + i)
      end do
      send_length = send_count(to_pe) * 2
   else
       to_pe = NULL_PE
       send_length = 1
   end if
   if(receive_count(from_pe) /= 0) then
       recv_length = receive_count(from_pe)*2
   else
       recv_length = 1
       from_pe = NULL_PE
   end if
   call mpp_transmit( send_buffer, send_length, to_pe, recv_buffer, recv_length, from_pe )
   if( from_pe.NE.NULL_PE )then
       rec_frac(1:recv_length/2) = recv_buffer(1:recv_length:2)
       cell(1:recv_length/2) = recv_buffer(2:recv_length:2)
! Put the data in the appropriate exchange cell and set avail to true
      do i = 1, recv_length/2
! Update the fractional area for this side and total area of the exchange cell
         ex_map%frac_area(cell(i), side) = rec_frac(i)
         ex_map%area(cell(i)) = ex_map%total_area(cell(i)) * &
            ex_map%frac_area(cell(i), 1) * ex_map%frac_area(cell(i), 2)
      end do
   end if
end do
call mpp_sync_self()
deallocate( send_buffer, recv_buffer, rec_frac, cell, send_count, receive_count)

end subroutine set_frac_area

!===========================================================================

function ex_cells_per_pe(bd_map, num_pes)

! Given a boundary map and the total number of pes determines
! the number of target component model cells that come from this
! exchange map for each processor and returns this count in an array.
! Finds the exchange grid map corresponding to this boundary map.
! Loop through the exchange grid to find out where its data goes
! for this boundary maps model (this side of exchange grid). 
! To do this, loop through each cell in the exchange grid and
! see what pe has the corresponding model cell. Increment a count
! for that pe (the returned) value. 

implicit none

type(boundary_map_type), intent(inout) :: bd_map
integer, intent(in) :: num_pes
integer :: ex_cells_per_pe(0:num_pes - 1)

integer :: i, side
type (exchange_map_type), pointer :: ex_map

! If not first call, just return data from bd_map
   ex_cells_per_pe = bd_map%ex_cells_per_pe
   return

end function ex_cells_per_pe

!===========================================================================

subroutine complete_bd_transmit(pe, x, y, part, ex_cell)

implicit none

integer, intent(in) :: pe, x(:), y(:), part(:), ex_cell(:)

integer :: length, end, start, i

length = size(x) + size(y) + size(part) + size(ex_cell)
allocate(buffer(pe)%data(length))

! Don't send zero length messages
if(length == 0) return

! Load all the data into a single real buffer for transmission to min latency
end = size(x)
buffer(pe)%data(1:end) = real(x)
start = end + 1
end = end + size(y)
buffer(pe)%data(start:end) = real(y)
start = end + 1
end = end + size(part)
buffer(pe)%data(start:end) = real(part)
start = end + 1
end = end + size(ex_cell)
buffer(pe)%data(start:end) = real(ex_cell)

! Call Balaji's send call
call mpp_transmit(buffer(pe)%data, length, pe, buffer(pe)%data, length, NULL_PE)

end subroutine complete_bd_transmit

!===========================================================================

subroutine complete_bd_receive(from_pe, x, y, part, ex_cell)

implicit none

integer, intent(in) :: from_pe
integer, intent(out) ::  x(:), y(:), part(:), ex_cell(:)

integer :: length, end, start
real, allocatable :: temp(:)

length = size(x) + size(y) + size(part) + size(ex_cell)
! Don't receive zero length messages
if(length == 0) return

allocate(temp(length))

! Call Balaji's receive call
call mpp_transmit(temp, length, NULL_PE, temp, length, from_pe)

! Unpack the data
end = size(x)
x = int(temp(1:end))
start = end + 1
end = end + size(y)
y = int(temp(start:end))
start = end + 1
end = end + size(part)
part = int(temp(start:end))
start = end + 1
end = end + size(ex_cell)
ex_cell = int(temp(start:end))

deallocate(temp)

end subroutine complete_bd_receive

!===========================================================================

subroutine end_exchange(num_pes)

implicit none

integer, intent(in) :: num_pes

integer :: i

! Need to synch at this point to avoid releasing buffers before they
! are transmitted; should be able to do this smarter at some point.
call mpp_sync_self()

! Release buffers from get_exchange or put_exchange calls
do i = 0, num_pes - 1
   if(associated(buffer(i)%data)) deallocate(buffer(i)%data)
end do

end subroutine end_exchange

!===========================================================================

subroutine fatal_error(message)

character (len = *), intent(in) :: message

! Call the mpp error handling facility with status FATAL
call mpp_error(FATAL, message)

end subroutine fatal_error

!===========================================================================

end module exchange_mod
