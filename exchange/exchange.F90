! NOTE: lon_lat routines need to have additional interface to help spectral
! modelers.

! NOTE: Need to consider supporting partial grids.

module exchange_mod

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

! NO WRITING FACILITIES AVAILABLE FOR MPP AT THIS POINT
! Public interfaces for writing out and reading in maps
!public write_exchange_grid_map, write_boundary_map, read_exchange_grid_map, &
!   read_boundary_map

!==========================================================================

! Type that maps from model cells to exchange cells
type boundary_map_type
   type(exchange_map_type), pointer :: ex
   integer :: side
   logical :: allocated
   integer :: model_id, num_x, num_y, num_parts, total_target_cells
   integer, pointer :: ind(:, :, :), num_targets(:, :, :)
   integer :: count
   integer, pointer :: pe(:), num(:), next(:)
   logical :: ex_cells_init = .false., model_cells_init = .false.
   integer, pointer :: ex_cells_per_pe(:), model_cells_per_pe(:)
end type boundary_map_type

! Back map from exchange grid to model cells
type exchange_map_type
   integer :: max_size, size
! All second components are 2, one for each side
   real, pointer :: total_area(:), area(:), frac_area(:, :)
   integer, pointer :: model_id(:, :), x(:, :), y(:, :), part(:, :), pe(:, :)
end type exchange_map_type
   
! Type to buffer data for transmits until all is sent
type transmit_buffer_type
   real, pointer :: data(:)
end type transmit_buffer_type

! VB: num_pes will come from above everywhere? How to set up buffer?
type (transmit_buffer_type), allocatable :: buffer(:)
logical :: first_lon_lat_size = .true.

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
! Boundary maps must be initialized before use; F90 sure is silly sometimes.
! Also link them to the appropriate exchange_map_type
!----------------------------------------------------------------------------

implicit none

type (boundary_map_type), intent(inout) :: bd_map
type (exchange_map_type), intent(inout), target :: map
integer, intent(in) :: side


! Get unique boundary map id
bd_map%model_id = next_id
next_id = next_id + 1
! Set number of cells to 0 and current number of cells computed to 0
bd_map%count = 0
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

end subroutine init_boundary_map

!===========================================================================

subroutine lon_lat1_size(blon1, blat1, blon2, blat2, mask1, mask2, &
   num_part1, num_part2, bd_map)

! Does mpp version of lon_lat_size; for continguous lon/lat boundary pairs
! passed in 1d array

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
if(.not. bd_map%allocated) then 
   call allocate_boundary_map(bd_map, size(blon1, 1), size(blat1, 1), num_part1)
   bd_map%allocated = .true.
endif

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
                     map%model_id(map%size, side) = bd_map%model_id
                     map%model_id(map%size, other_side) = bd_map2%model_id
! Initialize exchange grid back pointers
                     map%x(map%size, side) = i1
                     map%x(map%size, other_side) = i2
                     map%y(map%size, side) = j1
                     map%y(map%size, other_side) = j2
                     map%part(map%size, side) = k
                     map%part(map%size, other_side) = m
                     map%pe(map%size, side) = this_pe
                     map%pe(map%size, other_side) = other_pe
! Link in the boundary cells
                     call bd_map_add(bd_map, i1, j1, k, m, map%size, this_pe)
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
! lon_lat_size with call_size true. In that case, the size of storeage needed
! is computed but allocations for pointers are not done. Also called by
! lon_lat_map (call_size false) in which case pointer storeage is 
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
   map%model_id(num, 2), map%x(num, 2), map%y(num, 2), &
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

subroutine allocate_boundary_map(bd_map, num_x, num_y, num_parts)

!----------------------------------------------------------------------------
! Initializes a boundary map, bd_map, for a component model grid with
! space dimensions x by y and num_parts partitions.
!----------------------------------------------------------------------------

implicit none

type(boundary_map_type), intent(inout) :: bd_map
integer, intent(in) :: num_x, num_y, num_parts


bd_map%num_x = num_x
bd_map%num_y = num_y
bd_map%num_parts = num_parts

! Allocate the storage for the pointers to the exchange grid and
! Allocate the storage for the indirection arrays; next needed for initial setup
allocate(bd_map%ind(num_x, num_y, num_parts), &
   bd_map%num_targets(num_x, num_y, num_parts), &
   bd_map%pe(bd_map%total_target_cells), bd_map%num(bd_map%total_target_cells), &
   bd_map%next(bd_map%total_target_cells))

! Set allocation flag
bd_map%allocated = .true.

! Nullify all (integer implemented) pointers
bd_map%ind = -1
bd_map%num_targets = 0

! Illegal values for pe and and cell num and next pointer
bd_map%pe = -1
bd_map%num = -1
bd_map%next = -1

end subroutine allocate_boundary_map

!===========================================================================

subroutine bd_map_add(bd_map, i, j, k, m, ex_num, pe)

! GIven a boundary map cell and the number of the corresponding exchange cell
! as well as the pe on which it resides, set up the storage and index info
! in the boundary map structure. This is called only for side 1 boundary
! maps and only from a single place inside lon_lat_map. 

implicit none

type(boundary_map_type), intent(inout) :: bd_map
integer, intent(in) :: i, j, k, m, ex_num, pe

integer :: n, index

! Increment number of cells added to this map
bd_map%count = bd_map%count + 1

! If this is first exchange cell for this model cell start pointer chain
if(bd_map%ind(i, j, k) == -1) then
   bd_map%ind(i, j, k) = bd_map%count
else
! Follow the pointer chain to link in newest cell
   index = bd_map%ind(i, j, k)
   do n = 1, bd_map%num_targets(i, j, k) - 1
      index = bd_map%next(index)
   end do
! Now link in the next piece of the chain
bd_map%next(index) = bd_map%count
endif

! Put pointers into indirect chains
bd_map%num_targets(i, j, k) = bd_map%num_targets(i, j, k) + 1
bd_map%pe(bd_map%count) = pe
bd_map%num(bd_map%count) = ex_num

end subroutine bd_map_add

!===========================================================================

subroutine complete_side1_boundary_map(bd_map)

! Called for side 1 boundary maps after all side 1 calls to lon_lat_map
! have been completed. 

implicit none

! Copy the pe and num data to temporary storage and then reorder

type(boundary_map_type), intent(inout) :: bd_map

integer, dimension(bd_map%total_target_cells) :: temp_pe, temp_num
integer :: temp_ind(size(bd_map%ind, 1), size(bd_map%ind,2), size(bd_map%ind,3))
integer :: count, i, j, k, index, n


if(bd_map%side /= 1) &
   call fatal_error('complete_side1_boundary_map: called for side2')

temp_pe = bd_map%pe
temp_num = bd_map%num
temp_ind = bd_map%ind

count = 0

! Loop through each boundary map cell and follow pointers to all its data
do i = 1, bd_map%num_x
   do j = 1, bd_map%num_y
      do k = 1, bd_map%num_parts
         bd_map%ind(i, j, k) = count + 1
         index = temp_ind(i, j, k)
         do n = 1, bd_map%num_targets(i, j, k)
            count = count + 1
            bd_map%pe(count) = temp_pe(index)
            bd_map%num(count) = temp_num(index)
            index = bd_map%next(index)
         end do
      end do
   end do
end do

if(count /= bd_map%total_target_cells) &
   call fatal_error('complete_side1_boundary_map: cells copied not equal to total_target_cells')

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
integer, allocatable :: send_x(:), send_y(:), send_part(:), send_num(:)
integer, allocatable :: x(:), y(:), part(:)
integer :: this_pe, num_pes, pe_dist
integer :: total_target_cells, pe, index, start, end, i, j, k, from_pe, to_pe
real :: temp


! Get number of pes and identity of this pe
this_pe = mpp_pe()
num_pes = mpp_npes()
allocate(send_count(0:num_pes - 1), receive_count(0:num_pes - 1))

! Sets up the boundary map pointers for side 2 boundary maps
if(bd_map%side /= 2) &
   call fatal_error('complete_side2_boundary_map: called for side1')

! Get map pointer
ex_map => bd_map%ex

! Allocate storage space for this side 2 boundary map
call allocate_boundary_map(bd_map, num_x, num_y, num_parts)

! Get number of cells this pe's exchange grid sends to all other pes
send_count = ex_cells_per_pe(bd_map, num_pes)

! Just copy for info on the same pe
receive_count(this_pe) = send_count(this_pe)

! Loop through to send info from each pe to all other pes in turn
!vb: this cannot easily be converted to a broadcast since the message
!is different for each (from_pe,to_pe) combination.
!do from_pe = 0, num_pes - 1
!   if(this_pe == from_pe) then
!      do to_pe = 0, num_pes - 1
!         if(to_pe /= from_pe) call mpp_transmit(send_count(to_pe), 1, to_pe, &
!            send_count(to_pe), 1, NULL_PE)
!      end do
!      call mpp_sync_self()
!   else
!      call mpp_transmit(receive_count(from_pe), 1, NULL_PE, &
!         receive_count(from_pe), 1, from_pe)
!   end if
!end do
do pe_dist = 1,num_pes-1
   to_pe = mod( this_pe+pe_dist, num_pes )
   from_pe = mod( this_pe+num_pes-pe_dist, num_pes )
   call mpp_transmit( send_count(to_pe), 1, to_pe, receive_count(from_pe), 1, from_pe )
end do

!?????????? What about problem of zero length messages??? Just don't send them.

! Compute the total number of target cells for this boundary map 
total_target_cells = sum(receive_count)
bd_map%total_target_cells = total_target_cells
! Allocate buffers plus temporary storeage for sort
allocate(bd_map%pe(total_target_cells), bd_map%num(total_target_cells), &
   x(total_target_cells), y(total_target_cells), part(total_target_cells))

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
            send_part(send_count(to_pe)), send_num(send_count(to_pe)))
         index = 1
! Load up the tranmsit buffers for the current target pe
         do i = 1, ex_map%size
            if(ex_map%pe(i, bd_map%side) == to_pe .and. &
               ex_map%model_id(i, bd_map%side) == bd_map%model_id) then
               send_x(index) = ex_map%x(i, bd_map%side)
               send_y(index) = ex_map%y(i, bd_map%side)
               send_part(index) = ex_map%part(i, bd_map%side)
               send_num(index) = i
               index = index + 1
            end if
         end do
! Development error check on index here
         if(index /= send_count(to_pe) + 1) call fatal_error &
            ('complete_side2_boundary_map: index and send_count inconsistent')

! THIS MUST AGGREGATE AND BUFFER
         call complete_bd_transmit(to_pe, send_x, send_y, send_part, send_num)
         deallocate(send_x, send_y, send_part, send_num)
      end do
   end if

! Every pe is receiving a message from this one (self-transmit for now)
   start = end + 1
   end = end + receive_count(from_pe)
   call complete_bd_receive(from_pe, x(start:end), y(start:end), &
      part(start:end), bd_map%num(start:end))
   call mpp_sync
! Take care of pe number
   bd_map%pe(start:end) = from_pe
end do

! Development error check
if(end /= total_target_cells) &
   call fatal_error('complete_boundary_map: end /= to total target cells')

! Free up the communication buffering
call end_exchange(num_pes)

! When this is all done, need to sort so that cells are ordered by x, y, part
call cell_sort(x, y, part, bd_map%pe, bd_map%num, num_x, num_y, num_parts)

! Put this sorted data into the tables in boundary map, pe and num are done
do i = 1, total_target_cells
   bd_map%num_targets(x(i), y(i), part(i)) = &
      bd_map%num_targets(x(i), y(i), part(i)) + 1
   if(bd_map%num_targets(x(i), y(i), part(i)) == 1) &
      bd_map%ind(x(i), y(i), part(i)) = i
end do

deallocate(send_count, receive_count)

end subroutine complete_side2_boundary_map

!===========================================================================

subroutine cell_sort(x, y, part, pe, num, num_x, num_y, num_parts)

! Sorts by x, y, part and moves along the corresponding pe and num fields 
! for setting up ordered side 2 boundary maps.

implicit none

integer, intent(inout) :: x(:), y(:), part(:), pe(:), num(:)
integer, intent(in) :: num_x, num_y, num_parts

integer :: sx(size(x)), sy(size(y)), spart(size(part)), spe(size(pe)), &
   snum(size(num))
integer :: first(num_x, num_y, num_parts), num_cells(num_x, num_y, num_parts)

integer :: i, j, k, temp, sloc, next


! Get count of cells per each index set
num_cells = 0
do i = 1, size(x)
   num_cells(x(i), y(i), part(i)) = num_cells(x(i), y(i), part(i)) + 1
end do

! Next, convert to first available slot in sorted array
next = 1
do j = 1, num_y
   do i = 1, num_x
      do k = 1, num_parts
         first(i, j, k) = next
         next = next + num_cells(i, j, k)
      end do
   end do
end do

! Loop through the cells and put them in sorted order
do i = 1, size(x)
   sloc = first(x(i), y(i), part(i))
   sx(sloc) = x(i)
   sy(sloc) = y(i)
   spart(sloc) = part(i)
   spe(sloc) = pe(i)
   snum(sloc) = num(i)
   first(x(i), y(i), part(i)) = first(x(i), y(i), part(i)) + 1
end do








! Activate results of new sort routine
x = sx
y = sy
part = spart
pe = spe
num = snum










! Skip the old sort routine

!if(1 == 1) goto 10


! For now do silly bubble sort, make more efficient as needed

! 20 do i = 1, size(x) - 1
!   do j = i + 1, size(x)
!      if(y(i) > y(j) .or. (y(i) == y(j) .and. x(i) > x(j)) .or. &
!         (y(i) == y(j) .and. x(i) == x(j) .and. part(i) > part(j))) then
! Switch x, y, part, pe and num
!            temp = x(i)
!            x(i) = x(j)
!            x(j) = temp
!            temp = y(i)
!            y(i) = y(j)
!            y(j) = temp
!            temp = part(i)
!            part(i) = part(j)
!            part(j) = temp
!            temp = pe(i)
!            pe(i) = pe(j)
!            pe(j) = temp
!            temp = num(i)
!            num(i) = num(j)
!            num(j) = temp
!      end if
!   end do
!end do

! 10 continue

!do i = 1, size(x)
!   write(*, *) 'sort ', x(i), y(i), part(i), pe(i), num(i)
!end do

end subroutine cell_sort

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
real, allocatable :: send_data(:), send_area(:)
real, allocatable :: rec_data(:), rec_area(:)
integer, allocatable :: send_to_x(:), send_to_y(:), send_to_part(:)
integer, allocatable :: x(:), y(:), part(:)
integer :: i, j, k, index, cell, m, ind_num, side, from_pe, to_pe
integer :: rec_total, start, end
real ::total_area(size(model_data, 1), size(model_data, 2), size(model_data, 3))
integer :: send_count_max, recv_count_max, send_length, recv_length, pe_dist
real, allocatable, dimension(:) :: send_buffer, recv_buffer


! Get map pointer
ex_map => bd_map%ex

! Get number of pes and identity of this pe
this_pe = mpp_pe()
num_pes = mpp_npes()
allocate(send_count(0:num_pes - 1), receive_count(0:num_pes - 1))

! First determine how much data is going to each processor from exchange
! grid cells that correspond to this bd_map. I.E. how many exchange grid
! values will my pe receive from each other pe for this get.
receive_count = model_cells_per_pe(bd_map, num_pes)
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
total_area = 0.0

! Set up indices for reading into the receive buffers
start = 1
end = 0
! ringwise communication
do pe_dist = 0, num_pes-1
   to_pe = mod( this_pe+pe_dist, num_pes )
   from_pe = mod( this_pe+num_pes-pe_dist, num_pes )
   if(send_count(to_pe) /= 0) then
       call mpp_sync_self()
       index = 1
! Load up the transfer data (this is incredibly inefficient, fix it?)
       do i = 1, ex_map%size
! If this exchange map cell goes to pe, add its info to send buffers
! Need target pe number correct and that data is for this component model
          if(ex_map%pe(i, side) == to_pe .and. &
               ex_map%model_id(i, side) == bd_map%model_id) then
              send_buffer(index) = ex_data(i)
              send_buffer(index+1) = ex_map%area(i)
              send_buffer(index+2) = ex_map%x(i, side)
              send_buffer(index+3) = ex_map%y(i, side)
              send_buffer(index+4) = ex_map%part(i, side)
              index = index + 5
          end if
       end do
       send_length = send_count(to_pe)*5
! Development error check on index here
       if(index /= send_count(to_pe)*5 + 1) &
            call fatal_error('get_exchange_grid_3d: index doesnt fit send_count')
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
       start = end + 1
       end = end + recv_length/5
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
! NOTICE : ORDER ON THESE MUST BE FIXED TO REPRODUCE;
! THIS IS POTENTIALLY COSTLY: MAY WANT TO TURN IT OFF
call get_exchange_sort(rec_data, rec_area, x, y, part, size(model_data, 1), &
   size(model_data, 2), size(model_data, 3))

do i = 1, rec_total
   model_data(x(i), y(i), part(i)) = model_data(x(i), y(i), part(i)) + &
      rec_data(i) * rec_area(i)
   total_area(x(i), y(i), part(i)) = &
      total_area(x(i), y(i), part(i)) + rec_area(i)
end do

! Free up the local buffering
deallocate(rec_data, rec_area, x, y, part)

! Now do the weight normalization
do i = 1, size(model_data, 1)
   do j = 1, size(model_data, 2)
      do k = 1, size(model_data, 3)
         if(total_area(i, j, k) == 0.0) then
            if(model_data(i, j, k) /= 0.0) &
               call fatal_error('get_exchange_grid_3d: non-zero area')
         else
            model_data(i, j, k) = model_data(i, j, k) / total_area(i, j, k)
         end if
      end do
   end do
end do

deallocate(send_count, receive_count)

end subroutine get_exchange_grid_3d

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
      if(dat(i) > dat(j) .or. (dat(i) == dat(j) .and. area(i) > area(j))) then
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
allocate(send_count(0:num_pes - 1), receive_count(0:num_pes - 1))

! First determine how much data is going to each processor from this one?
send_count = model_cells_per_pe(bd_map, num_pes)
send_count_max = maxval(send_count)
! Allocate storage to do the gather of data to be sent to this pe
allocate( send_buffer(send_count_max*2) )

! Next determine how much data is coming from each processor
side = bd_map%side
receive_count = ex_cells_per_pe(bd_map, num_pes)
recv_count_max = maxval(receive_count)
allocate(rec_data(recv_count_max), cell(recv_count_max))
allocate( recv_buffer(recv_count_max*2) )

! ringwise communication
do pe_dist = 0, num_pes-1
   to_pe = mod( this_pe+pe_dist, num_pes )
   from_pe = mod( this_pe+num_pes-pe_dist, num_pes )
   if(send_count(to_pe) /= 0) then
       call mpp_sync_self()
       index = 1
! Loop through each boundary map cell (x, y, part)
       do i = 1, bd_map%num_x
          do j = 1, bd_map%num_y
             do k = 1, bd_map%num_parts
! Loop through each exchange cell that is a target from this boundary cell
                do m = 1, bd_map%num_targets(i, j, k)
! Find section of indirection array for these exchange cells
                   ind_num = bd_map%ind(i, j, k) + m - 1
! If the pe number of these exchange cells is current, put them in the message
                   if(bd_map%pe(ind_num) == to_pe) then
                       send_buffer(index) = model_data(i, j, k)
                       send_buffer(index+1) = bd_map%num(ind_num)
                       index = index + 2
                   end if
                end do
             end do
          end do
       end do
       send_length = send_count(to_pe)*2
! Development error check on index here
       if(index /= send_count(to_pe)*2 + 1) &
            call fatal_error('put_exchange_grid_3d: index doesnt fit send_count')
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
         if(present(avail)) avail(cell(i)) = .true.
      end do
   end if
end do
call mpp_sync_self()
deallocate( send_buffer, recv_buffer )
deallocate(rec_data, cell)

deallocate(send_count, receive_count)

end subroutine put_exchange_grid_3d

!===========================================================================

subroutine set_frac_area(frac_area, bd_map)

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


! Get number of pes and identity of this pe
this_pe = mpp_pe()
num_pes = mpp_npes()
allocate(send_count(0:num_pes - 1), receive_count(0:num_pes - 1))

! Error check for fractional area of cell not adding up to 1.0
if(any(abs(sum(frac_area, 3) - 1.0) > area_tol)) &
   call fatal_error('set_frac_area; sum of fractions is not 1')

! Get pointer to corresponding exchange grid map
ex_map => bd_map%ex

! First determine how much data is going to each processor from this one
send_count = model_cells_per_pe(bd_map, num_pes)
send_count_max = maxval(send_count)
! Allocate storage to do the gather of data to be sent to this pe
allocate( send_buffer(send_count_max*2) )

! Next determine how much data is coming from each processor to this one.
side = bd_map%side
receive_count = ex_cells_per_pe(bd_map, num_pes)
recv_count_max = maxval(receive_count)
allocate( recv_buffer(recv_count_max*2) )
allocate(rec_frac(recv_count_max), cell(recv_count_max))

do pe_dist = 0, num_pes-1
   to_pe = mod( this_pe+pe_dist, num_pes )
   from_pe = mod( this_pe+num_pes-pe_dist, num_pes )
   if(send_count(to_pe) /= 0) then
       call mpp_sync_self()
       index = 1
! Loop through each boundary map cell (x, y, part)
       do i = 1, bd_map%num_x
          do j = 1, bd_map%num_y
             do k = 1, bd_map%num_parts
! Loop through each exchange cell that is a target from this boundary cell
                do m = 1, bd_map%num_targets(i, j, k)
! Find section of indirection array for these exchange cells
                   ind_num = bd_map%ind(i, j, k) + m - 1
! If the pe number of these exchange cells is current, put them in the message
                   if(bd_map%pe(ind_num) == to_pe) then
                       send_buffer(index) = frac_area(i, j, k)
                       send_buffer(index+1) = bd_map%num(ind_num)
                       index = index + 2
                   end if
                end do
             end do
          end do
       end do
       send_length = send_count(to_pe)*2
! Development error check on index here
       if(index /= send_count(to_pe)*2 + 1) &
            call fatal_error('put_exchange_grid_3d: index doesnt fit send_count')
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
deallocate( send_buffer, recv_buffer )
deallocate(rec_frac, cell)

deallocate(send_count, receive_count)

end subroutine set_frac_area

!===========================================================================

function model_cells_per_pe(bd_map, num_pes)

! Given a boundary map and the total number of pes, determine how many
! exchange grid cells corresponding to model cells in this bd_map
! reside on each pe. 

! Could enhance efficiency by putting this in table on first call and
! using that thereafter?

implicit none

type(boundary_map_type), intent(inout) :: bd_map
integer, intent(in) :: num_pes
integer :: model_cells_per_pe(0:num_pes - 1)

integer i, j, k, m, ind_num

! If not first call, just return data from bd_map
if(bd_map%model_cells_init) then
   model_cells_per_pe = bd_map%model_cells_per_pe
   return
endif

model_cells_per_pe = 0
! Loop through each bd_map cell (x, y, num_parts)
do i = 1, bd_map%num_x
   do j = 1, bd_map%num_y
      do k = 1, bd_map%num_parts
! Each bd_map cell points to 1 or more exchange cells, loop through these
         do m = 1, bd_map%num_targets(i, j, k)
! The ind array has an index into the pe array for each cell
            ind_num = bd_map%ind(i, j, k) + m - 1
! Increment a count for the pe on which this exchange grid cell resides
            model_cells_per_pe(bd_map%pe(ind_num)) = &
               model_cells_per_pe(bd_map%pe(ind_num)) + 1
         end do
      end do
   end do
end do

! Load values into the bd_map storage and set initialized to true
allocate(bd_map%model_cells_per_pe(0:num_pes - 1))
do i = 0, num_pes - 1
   bd_map%model_cells_per_pe = model_cells_per_pe
end do
bd_map%model_cells_init = .true.

end function model_cells_per_pe

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
if(bd_map%ex_cells_init) then
   ex_cells_per_pe = bd_map%ex_cells_per_pe
   return
endif

! Get map pointer
ex_map => bd_map%ex
side = bd_map%side

ex_cells_per_pe = 0

! Loop through each exchange grid cell; find out what pe has the 
! target component model cell for this side of the exchange grid,
! and increment a counter of how many cells on each pe.
do i = 1, ex_map%size
   if(ex_map%model_id(i, side) == bd_map%model_id) &
   ex_cells_per_pe(ex_map%pe(i, side)) = ex_cells_per_pe(ex_map%pe(i, side)) + 1
end do

! Load values into the bd_map storage and set initialized to true
allocate(bd_map%ex_cells_per_pe(0:num_pes - 1))
do i = 0, num_pes - 1
   bd_map%ex_cells_per_pe = ex_cells_per_pe
end do
bd_map%ex_cells_init = .true.

end function ex_cells_per_pe

!===========================================================================

subroutine complete_bd_transmit(pe, x, y, part, num)

implicit none

integer, intent(in) :: pe, x(:), y(:), part(:), num(:)

integer :: length, end, start, i

length = size(x) + size(y) + size(part) + size(num)
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
end = end + size(num)
buffer(pe)%data(start:end) = real(num)

! Call Balaji's send call
call mpp_transmit(buffer(pe)%data, length, pe, buffer(pe)%data, length, NULL_PE)

end subroutine complete_bd_transmit

!===========================================================================

subroutine complete_bd_receive(from_pe, x, y, part, num)

implicit none

integer, intent(in) :: from_pe
integer, intent(out) ::  x(:), y(:), part(:), num(:)

integer :: length, end, start
real, allocatable :: temp(:)

length = size(x) + size(y) + size(part) + size(num)
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
end = end + size(num)
num = int(temp(start:end))

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
