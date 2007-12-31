!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! xgrid_mod - implements exchange grids.  An exchange grid is the grid whose
!             boundary set is the union of the boundaries of the participating
!             grids.  The exchange grid is the coarsest grid that is a 
!             refinement of each of the participating grids.  Every exchange
!             grid cell is a subarea of one and only one cell in each of the
!             participating grids.  The exchange grid has two purposes:
!
!               (1) The exchange cell areas are used as weights for
!                   conservative interpolation between model grids.
!
!               (2) Computation of surface fluxes takes place on it,
!                   thereby using the finest scale data obtainable.
!
!             The exchange cells are the 2D intersections between cells of the
!             participating grids.  They are computed elsewhere and are
!             read here from a NetCDF grid file as a sequence of quintuples
!             (i and j on each of two grids and the cell area).
!
!             Each processing element (PE) computes a subdomain of each of the
!             participating grids as well as a subset of the exchange cells.
!             The geographic regions corresponding to these subdomains will,
!             in general, not be the same so communication must occur between
!             the PEs.  The scheme for doing this is as follows.  A distinction
!             is drawn between the participating grids.  There is a single
!             "side 1" grid and it does not have partitions (sub-grid surface
!             types).  There are one or more "side 2" grids and they may have
!             more than 1 partition.  In standard usage, the atmosphere grid is
!             on side 1 and the land and sea ice grids are on side 2.  The set
!             of exchange cells computed on a PE corresponds to its side 2
!             geographic region(s).  Communication between the PEs takes place
!             on the side 1 grid.  Note:  this scheme does not generally allow
!             reproduction of answers across varying PE counts.  This is
!             because, in the side 1 "get", exchange cells are first summed
!             locally onto a side 1 grid, then these side 1 contributions are
!             further summed after they have been communicated to their target
!             PE.  For the make_exchange_reproduce option, a special side 1 get
!             is used.  This get communicates individual exchange cells.  The
!             cells are summed in the order they appear in the grid spec. file.
!                                    Michael Winton (Michael.Winton@noaa.gov) Oct 2001
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module xgrid_mod

! <CONTACT EMAIL="Michael.Winton@noaa.gov">
!   Michael Winton
! </CONTACT>
! <CONTACT EMAIL="Zhi.Liang@noaa.gov">
!   Zhi Liang
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!    <TT>xgrid_mod</TT> implements exchange grids for coupled models running on
!     multiple processors.  An exchange grid is formed from the union of
!     the bounding lines of the two (logically rectangular) participating
!     grids.  The exchange grid is therefore the coarsest grid that is a
!     refinement of both participating grids.  Exchange grids are used for
!     two purposes by coupled models:  (1) conservative interpolation of fields
!     between models uses the exchange grid cell areas as weights and
!     (2) the surface flux calculation takes place on the exchange grid thereby
!     using the finest scale data available.  <TT>xgrid_mod</TT> uses a NetCDF grid
!     specification file containing the grid cell overlaps in combination with
!     the <LINK SRC="ftp://ftp.gfdl.gov/pub/vb/mpp/mpp_domains.F90">
!     <TT>mpp_domains</TT></LINK> domain decomposition information to determine 
!     the grid and processor connectivities.
! </OVERVIEW>

! <DESCRIPTION>
!     <TT>xgrid_mod</TT> is initialized with a list of model identifiers (three characters
!     each), a list of <TT>mpp_domains</TT> domain data structures, and a grid specification
!     file name.  The first element in the lists refers to the "side one" grid.
!     The remaining elements are on "side two".  Thus, there may only be a single
!     side one grid and it is further restricted to have no partitions (sub-grid
!     areal divisions).  In standard usage, the atmosphere model is on side one
!     and the land and sea ice models are on side two.  <TT>xgrid_mod</TT> performs
!     interprocessor communication on the side one grid.  Exchange grid variables
!     contain no data for zero sized partitions.  The size and format of exchange
!     grid variables change every time the partition sizes or number of partitions
!     are modified with a <TT>set_frac_area</TT> call on a participating side two grid.
!     Existing exchange grid variables cannot be properly interpreted after
!     that time; new ones must be allocated and assigned with the <TT>put_to_xgrid</TT>
!     call.
! </DESCRIPTION>

! <DATA NAME="xmap_type"  TYPE=""  >
!   The fields of xmap_type are all private.
! </DATA>

! <DATASET NAME="">
!     <TT>xgrid_mod</TT> reads a NetCDF grid specification file to determine the
!     grid and processor connectivities.  The exchange grids are defined
!     by a sequence of quintuples:  the <TT>i/j</TT> indices of the intersecting
!     cells of the two participating grids and their areal overlap.
!     The names of the five fields are generated automatically from the
!     three character ids of the participating grids.  For example, if
!     the side one grid id is "ATM" and the side two grid id is "OCN",
!     <TT>xgrid_mod</TT> expects to find the following five fields in the grid
!     specification file:  <TT>I_ATM_ATMxOCN, J_ATM_ATMxOCN, I_OCN_ATMxOCN,
!     J_OCN_ATMxOCN, and AREA_ATMxOCN</TT>.  These fields may be generated
!     by the <TT>make_xgrids</TT> utility.
! </DATASET>
use       fms_mod,   only: file_exist, open_namelist_file, check_nml_error,  &
                           error_mesg, close_file, FATAL, NOTE, stdlog,      &
                           write_version_number, read_data, field_exist,     &
                           field_size, lowercase, string,                    &
                           get_mosaic_tile_grid
use mpp_mod,         only: mpp_npes, mpp_pe, mpp_root_pe, mpp_send, mpp_recv, &
                           mpp_sync_self, stdout
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                           Domain2d, mpp_global_sum, mpp_update_domains,    &
                           mpp_modify_domain, mpp_get_data_domain, XUPDATE, &
                           YUPDATE, mpp_get_current_ntile, mpp_get_tile_id, &
                           mpp_get_ntile_count, mpp_get_tile_list,          &
                           mpp_get_global_domains
use mpp_io_mod,      only: mpp_open, MPP_MULTI, MPP_SINGLE, MPP_OVERWR
use constants_mod,   only: PI
use mosaic_mod,          only: get_mosaic_xgrid, get_mosaic_xgrid_size
use stock_constants_mod, only: ISTOCK_TOP, ISTOCK_BOTTOM, ISTOCK_SIDE, STOCK_NAMES, STOCK_UNITS, NELEMS, stocks_file

implicit none
private

public xmap_type, setup_xmap, set_frac_area, put_to_xgrid, get_from_xgrid, &
       xgrid_count, some, conservation_check, xgrid_init, &
       AREA_ATM_SPHERE, AREA_LND_SPHERE, AREA_OCN_SPHERE, &
       AREA_ATM_MODEL, AREA_LND_MODEL, AREA_OCN_MODEL, &
       get_ocean_model_area_elements
!--- paramters that determine the remapping method
integer, parameter :: FIRST_ORDER        = 1
integer, parameter :: SECOND_ORDER       = 2
integer, parameter :: SECOND_ORDER_MERID = 3
integer, parameter :: SECOND_ORDER_ZONAL = 4 
integer, parameter :: VERSION1           = 1 ! grid spec file
integer, parameter :: VERSION2           = 2 ! mosaic grid file

! <NAMELIST NAME="xgrid_nml">
!   <DATA NAME="make_exchange_reproduce" TYPE="logical"  DEFAULT=".false.">
!     Set to .true. to make <TT>xgrid_mod</TT> reproduce answers on different
!     numbers of PEs.  This option has a considerable performance impact.
!   </DATA>
!   <DATA NAME="interp_method" TYPE="character(len=64)"  DEFAULT=" 'first_order' ">
!     exchange grid interpolation method. It has four options: 
!     "first_order", "second_order", "second_order_merid", "second_order_zonal".
!   </DATA>
logical :: make_exchange_reproduce = .false. ! exactly same on different # PEs
character(len=64) :: interp_method = 'first_order'
logical :: debug_stocks = .false. 
namelist /xgrid_nml/ make_exchange_reproduce, interp_method, debug_stocks
! </NAMELIST>
logical :: init = .true.
integer :: remapping_method

! Area elements used inside each model
real, allocatable, dimension(:,:) :: AREA_ATM_MODEL, AREA_LND_MODEL, AREA_OCN_MODEL
! Area elements based on a the spherical model used by the ICE layer
real, allocatable, dimension(:,:) :: AREA_ATM_SPHERE, AREA_LND_SPHERE, AREA_OCN_SPHERE

! <INTERFACE NAME="put_to_xgrid">

!   <OVERVIEW>
!     Scatters data from model grid onto exchange grid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Scatters data from model grid onto exchange grid.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call put_to_xgrid(d, grid_id, x, xmap, remap_order)
!   </TEMPLATE>
!   <IN NAME="d"  TYPE="real"  > </IN>
!   <IN NAME="grid_id"  TYPE=" character(len=3)"  > </IN>
!   <INOUT NAME="x"  TYPE="real"  > </INOUT>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>
!   <IN NAME="remap_method" TYPE="integer,optional">
!     exchange grid interpolation method. It has four possible values: 
!     FIRST_ORDER (=1), SECOND_ORDER(=2), SECOND_ORDER_MERID(=3) and
!     SECOND_ORDER_ZONAL(=4). Default value is FIRST_ORDER.
!   </IN>
interface put_to_xgrid
  module procedure put_side1_to_xgrid
  module procedure put_side2_to_xgrid
end interface
! </INTERFACE>

! <INTERFACE NAME="get_from_xgrid">

!   <OVERVIEW>
!     Sums data from exchange grid to model grid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Sums data from exchange grid to model grid.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_from_xgrid(d, grid_id, x, xmap)
!   </TEMPLATE>
!   <IN NAME="x"  TYPE="real"  > </IN>
!   <IN NAME="grid_id"  TYPE=" character(len=3)"  > </IN>
!   <OUT NAME="d"  TYPE="real"  > </OUT>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>
interface get_from_xgrid
  module procedure get_side1_from_xgrid
  module procedure get_side2_from_xgrid
end interface
! </INTERFACE>

! <INTERFACE NAME="conservation_check">

!   <OVERVIEW>
!     Returns three numbers which are the global sum of a variable.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns three numbers which are the global sum of a
!     variable (1) on its home model grid, (2) after interpolation to the other
!     side grid(s), and (3) after re_interpolation back onto its home side grid(s).
!     Conservation_check must be called by all PEs to work properly.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call conservation_check(d, grid_id, xmap,remap_order)
!   </TEMPLATE>
!   <IN NAME="d"  TYPE="real" DIM="(:,:)" > </IN>
!   <IN NAME="grid_id"  TYPE="character(len=3)"  > </IN>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>
!   <OUT NAME="" TYPE="real" DIM="3">The global sum of a variable.</OUT>
!   <IN NAME="remap_method" TYPE="integer,optional">
!   </IN>
interface conservation_check
  module procedure conservation_check_side1
  module procedure conservation_check_side2
end interface
! </INTERFACE>

!--- private interface
interface grad_zonal
  module procedure grad_zonal_1d
  module procedure grad_zonal_2d
end interface

interface grad_merid
  module procedure grad_merid_1d
  module procedure grad_merid_2d
end interface

type xcell_type
  integer :: i1, j1, i2, j2 ! indices of cell in model arrays on both sides
  integer :: pe             ! other side pe that has this cell
  integer :: tile           ! tile index of side 1 mosaic.
  real    :: area           ! geographic area of exchange cell
!  real    :: area1_ratio     !(= x_area/grid1_area), will be added in the future to improve efficiency
!  real    :: area2_ratio     !(= x_area/grid2_area), will be added in the future to improve efficiency
  real    :: di, dj         ! Weight for the gradient of flux
end type xcell_type

type grid_type
  character(len=3)                :: id                               ! grid identifier
  integer                         :: ntile                            ! number of tiles in mosaic
  integer                         :: ni, nj                           ! max of global size of all the tiles
  integer, pointer, dimension(:)  :: tile =>NULL()                    ! tile id ( pe index )
  integer, pointer, dimension(:)  :: is =>NULL(), ie =>NULL()         ! domain - i-range (pe index)
  integer, pointer, dimension(:)  :: js =>NULL(), je =>NULL()         ! domain - j-range (pe index)
  integer, pointer                :: is_me =>NULL(),  ie_me =>NULL()  ! my domain - i-range
  integer, pointer                :: js_me =>NULL(),  je_me =>NULL()  ! my domain - j-range
  integer                         :: isd_me, ied_me                   ! my data domain - i-range
  integer                         :: jsd_me, jed_me                   ! my data domain - j-range
  integer, pointer                :: tile_me                          ! my tile id
  integer                         :: im , jm , km                     ! global domain range
  real, pointer, dimension(:)     :: lon =>NULL(), lat =>NULL()       ! center of global grids
  real, pointer, dimension(:,:)   :: geolon=>NULL(), geolat=>NULL()   ! geographical grid center
  real, pointer, dimension(:,:,:) :: frac_area =>NULL()               ! partition fractions
  real, pointer, dimension(:,:)   :: area =>NULL()                    ! cell area
  real, pointer, dimension(:,:)   :: area_inv =>NULL()                ! 1 / area for normalization
  integer                         :: first, last                      ! xgrid index range
  integer                         :: size                             ! # xcell patterns
  type(xcell_type), pointer, dimension(:) :: x =>NULL()               ! xcell patterns
  integer                         :: size_repro                       ! # side 1 patterns for repro
  type(xcell_type), pointer, dimension(:) :: x_repro =>NULL()         ! side 1 patterns for repro
  type(Domain2d) :: domain                                            ! used for conservation checks
  type(Domain2d) :: domain_with_halo                                  ! used for second order remapping
  logical        :: is_latlon                                         ! indicate if the grid is lat-lon or not.
end type grid_type

type x1_type
  integer :: i, j
  real    :: area   ! (= geographic area * frac_area)
!  real    :: area_ratio !(= x1_area/grid1_area) ! will be added in the future to improve efficiency
  real    :: di, dj ! weight for the gradient of flux
  integer :: tile           ! tile index of side 1 mosaic.
end type x1_type

type x2_type
  integer :: i, j, k
  real    :: area   ! geographic area of exchange cell
!  real    :: area_ratio !(=x2_area/grid2_area )  ! will be added in the future to improve efficiency
end type x2_type

type xmap_type
  private
  integer :: size            ! # of exchange grid cells with area > 0 on this pe

  integer :: me, npes, root_pe
  logical, pointer, dimension(:) :: your1my2  =>NULL()! true if side 1 domain on
                                                      ! indexed pe overlaps side 2
                                                      ! domain on this pe
  logical, pointer, dimension(:) :: your2my1 =>NULL() ! true if a side 2 domain on
                                                      ! indexed pe overlaps side 1
                                                      ! domain on this pe

  type (grid_type), pointer, dimension(:) :: grids =>NULL() ! 1st grid is side 1;
                                                            ! rest on side 2
  !
  ! Description of the individual exchange grid cells (index is cell #)
  !
  type(x1_type), pointer, dimension(:) :: x1 =>NULL() ! side 1 info
  type(x2_type), pointer, dimension(:) :: x2 =>NULL() ! side 2 info

  real, pointer,    dimension(:) :: send_buffer =>NULL() ! for non-blocking sends
  real, pointer,    dimension(:) :: recv_buffer =>NULL() ! for non-blocking recv
  integer, pointer, dimension(:) :: send_count_repro =>NULL()
  integer, pointer, dimension(:) :: recv_count_repro  =>NULL()
  integer :: version                                  ! version of xgrids. version=VERSION! is for grid_spec file 
                                                      ! and version=VERSION2 is for mosaic grid.
end type xmap_type

!-----------------------------------------------------------------------
 character(len=128) :: version = '$Id: xgrid.F90,v 15.0 2007/08/14 04:13:45 fms Exp $'
 character(len=128) :: tagname = '$Name: omsk_2007_12 $'

 real, parameter                              :: EPS = 1.0e-10
 logical :: module_is_initialized = .FALSE.

 ! The following is required to compute stocks of water, heat, ...

  integer, parameter :: NSIDES  = 3         ! top, bottom, side
 
  type stock_type
     real                   :: q_start    = 0.0    ! stock at start
     ! delta_t * surf integr of flux
     ! these are the stock increments at the present time, one for
     ! each side (ISTOCK_TOP, ISTOCK_BOTTOM, ISTOCK_SIDE)
     real                   :: dq(NSIDES) = 0.0    ! stock increments at present time      
  end type stock_type

  interface stock_move
     module procedure stock_move_3d, stock_move_2d
  end interface

  public stock_move, stock_type, stock_print, get_index_range, stock_integrate_2d

contains

!#######################################################################

logical function in_box(i, j, is, ie, js, je)
integer :: i, j, is, ie, js, je

  in_box = (i>=is) .and. (i<=ie) .and. (j>=js) .and. (j<=je)
end function in_box

!#######################################################################

! <SUBROUTINE NAME="xgrid_init">

!   <OVERVIEW>
!     Initialize the xgrid_mod. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Initialization routine for the xgrid module. It reads the xgrid_nml,  
!     writes the version information and xgrid_nml to the log file.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call xgrid_init ( )
!   </TEMPLATE>
!   <OUT NAME="remap_method" TYPE="integer">
!     exchange grid interpolation method. It has four possible values: 
!     FIRST_ORDER (=1), SECOND_ORDER(=2), SECOND_ORDER_MERID(=3) and
!     SECOND_ORDER_ZONAL(=4)
!   </OUT>
subroutine xgrid_init(remap_method) 
  integer, intent(out) :: remap_method

  integer :: unit, ierr, io

  if (module_is_initialized) return
  module_is_initialized = .TRUE.

  if ( file_exist( 'input.nml' ) ) then
      unit = open_namelist_file ( )
      ierr = 1
      do while ( ierr /= 0 )
        read ( unit,  nml = xgrid_nml, iostat = io, end = 10 )
        ierr = check_nml_error ( io, 'xgrid_nml' )
      enddo
  10 continue
      call close_file ( unit )
  endif

!--------- write version number and namelist ------------------
  call write_version_number (version, tagname)

  unit = stdlog ( )
  if ( mpp_pe() == mpp_root_pe() ) write (unit,nml=xgrid_nml)
  call close_file (unit)

!--------- check interp_method has suitable value

  select case(trim(interp_method))
  case('first_order')
     remap_method = FIRST_ORDER
  case('second_order')
     remap_method = SECOND_ORDER
  case('second_order_merid')
     remap_method = SECOND_ORDER_MERID
  case('second_order_zonal')
     remap_method = SECOND_ORDER_ZONAL
  case default
     call error_mesg('xgrid_mod', ' nml interp_method = ' //trim(interp_method)// &
      ' is not a valid namelist option', FATAL)
  end select
  
  remapping_method = remap_method

end subroutine xgrid_init
! </SUBROUTINE>

!#######################################################################

subroutine load_xgrid (xmap, grid, grid_file, grid1_id, grid_id, tile1, tile2, use_higher_order, complete)
type(xmap_type), intent(inout)         :: xmap
type(grid_type), intent(inout)         :: grid
character(len=*), intent(in)           :: grid_file
character(len=3), intent(in)           :: grid1_id, grid_id
integer,          intent(in)           :: tile1, tile2
logical,        intent(in)             :: use_higher_order
logical,        intent(in)             :: complete   

  integer, allocatable, dimension(:) :: i1, j1, i2, j2            ! xgrid quintuples
  real,    allocatable, dimension(:) :: area, di, dj              ! from grid file
  type (grid_type), pointer, save    :: grid1 =>NULL()
  integer                            :: l, ll, ll_repro, p, siz(4), nxgrid, size_prev
  type(xcell_type), allocatable      :: x_local(:)
  integer                            :: size_repro

  grid1 => xmap%grids(1)
  select case(xmap%version)
  case(VERSION1)
     call field_size(grid_file, 'AREA_'//grid1_id//'x'//grid_id, siz)
     nxgrid = siz(1);
     if(nxgrid==0 .or. siz(4) == 0 ) return
     allocate(i1(nxgrid), j1(nxgrid), i2(nxgrid), j2(nxgrid), area(nxgrid))
     call read_data(grid_file, 'I_'//grid1_id//'_'//grid1_id//'x'//grid_id, i1, no_domain=.true.)
     call read_data(grid_file, 'J_'//grid1_id//'_'//grid1_id//'x'//grid_id, j1, no_domain=.true.)
     call read_data(grid_file, 'I_'//grid_id//'_'//grid1_id//'x'//grid_id, i2, no_domain=.true.)
     call read_data(grid_file, 'J_'//grid_id//'_'//grid1_id//'x'//grid_id, j2, no_domain=.true.)
     call read_data(grid_file, 'AREA_'//grid1_id//'x'//grid_id, area, no_domain=.true.)
     if(use_higher_order) then
        allocate(di(nxgrid), dj(nxgrid))
        call read_data(grid_file, 'DI_'//grid1_id//'x'//grid_id, di, no_domain=.true.)
        call read_data(grid_file, 'DJ_'//grid1_id//'x'//grid_id, dj, no_domain=.true.)
     end if
  case(VERSION2)
     !--- max_size is the exchange grid size between super grid.
     nxgrid = get_mosaic_xgrid_size(grid_file)
     allocate(i1(nxgrid), j1(nxgrid), i2(nxgrid), j2(nxgrid), area(nxgrid) )
     if(use_higher_order) then
        allocate(di(nxgrid), dj(nxgrid))
        call get_mosaic_xgrid(grid_file, i1, j1, i2, j2, area, di, dj)
     else
        call get_mosaic_xgrid(grid_file, i1, j1, i2, j2, area)
     end if
  end select

  size_repro = 0
  if(grid1%tile_me == tile1) then
     do l=1,nxgrid
        if (in_box(i1(l), j1(l), grid1%is_me, grid1%ie_me, grid1%js_me, grid1%je_me) ) then
           grid1%area(i1(l),j1(l)) = grid1%area(i1(l),j1(l))+area(l)
           do p=0,xmap%npes-1
              if (grid%tile(p) == tile2) then
                 if (in_box(i2(l), j2(l), grid%is(p), grid%ie(p), grid%js(p), grid%je(p)))  then
                    xmap%your2my1(p) = .true.
                 end if
              end if
           end do
           size_repro = size_repro + 1
        end if
     end do
  end if

  size_prev = grid%size

  if(grid%tile_me == tile2) then
     do l=1,nxgrid
        if (in_box(i2(l), j2(l), grid%is_me, grid%ie_me, grid%js_me, grid%je_me) ) then
           grid%size = grid%size + 1
           grid%area(i2(l),j2(l)) = grid%area(i2(l),j2(l))+area(l)
           do p=0,xmap%npes-1
              if(grid1%tile(p) == tile1) then
                 if (in_box(i1(l), j1(l), grid1%is(p), grid1%ie(p), &
                      grid1%js(p), grid1%je(p))) then
                    xmap%your1my2(p) = .true.
                 end if
              end if
           end do
        end if
     end do
  end if

  if(grid%size > size_prev) then
     if(size_prev > 0) then ! need to extend data
        allocate(x_local(size_prev))
        x_local = grid%x
        if(ASSOCIATED(grid%x)) deallocate(grid%x)
        allocate( grid%x( grid%size ) )
        grid%x(1:size_prev) = x_local
        deallocate(x_local)
     else
        allocate( grid%x( grid%size ) )
        grid%x%di = 0; grid%x%dj = 0
     end if
  end if

  ll = size_prev
  if( grid%tile_me == tile2 ) then ! me is tile2
     do l=1,nxgrid
        if (in_box(i2(l), j2(l), grid%is_me, grid%ie_me, grid%js_me, grid%je_me)) then
           ! insert in this grids cell pattern list and add area to side 2 area
           ll = ll + 1
           grid%x(ll)%i1   = i1(l); grid%x(ll)%i2   = i2(l)
           grid%x(ll)%j1   = j1(l); grid%x(ll)%j2   = j2(l)
           grid%x(ll)%tile = tile1
           grid%x(ll)%area = area(l)
           if(use_higher_order) then
              grid%x(ll)%di  = di(l)
              grid%x(ll)%dj  = dj(l)
           end if

           if (make_exchange_reproduce) then
              do p=0,xmap%npes-1
                 if(grid1%tile(p) == tile1) then
                    if (in_box(i1(l), j1(l), grid1%is(p), grid1%ie(p), &
                         grid1%js(p), grid1%je(p))) then
                       grid%x(ll)%pe = p + xmap%root_pe
                    end if
                 end if
              end do
           end if ! make_exchange reproduce
        end if
     end do
  end if

  if (make_exchange_reproduce .and. grid1%tile_me == tile1 .and. size_repro > 0) then
     ll_repro = grid%size_repro
     grid%size_repro = ll_repro + size_repro
     if(ll_repro > 0) then  ! extend data
        allocate(x_local(ll_repro))
        x_local = grid%x_repro
        if(ASSOCIATED(grid%x_repro)) deallocate(grid%x_repro)
        allocate( grid%x_repro(grid%size_repro ) )
        grid%x_repro(1:ll_repro) = x_local
        deallocate(x_local)
     else
        allocate( grid%x_repro( grid%size_repro ) )
        grid%x_repro%di = 0; grid%x_repro%dj = 0
     end if
     do l=1,nxgrid
        if (in_box(i1(l),j1(l), grid1%is_me,grid1%ie_me, grid1%js_me,grid1%je_me) ) then
           ll_repro = ll_repro + 1
           grid%x_repro(ll_repro)%i1   = i1(l); grid%x_repro(ll_repro)%i2   = i2(l)
           grid%x_repro(ll_repro)%j1   = j1(l); grid%x_repro(ll_repro)%j2   = j2(l)
           grid%x_repro(ll_repro)%tile = tile1
           grid%x_repro(ll_repro)%area = area(l)
           if(use_higher_order) then
              grid%x_repro(ll_repro)%di  = di(l)
              grid%x_repro(ll_repro)%dj  = dj(l)
           end if
           do p=0,xmap%npes-1
              if(grid%tile(p) == tile2) then
                 if (in_box(i2(l), j2(l), grid%is(p), grid%ie(p), &
                      grid%js(p), grid%je(p))) then
                    grid%x_repro(ll_repro)%pe = p + xmap%root_pe
                 end if
              end if
           end do
        end if ! make_exchange_reproduce
     end do
  end if

  if(complete ) then  
     grid%area_inv = 0.0;
     where (grid%area>0.0) grid%area_inv = 1.0/grid%area
  end if

  deallocate(i1, j1, i2, j2, area)
  if(use_higher_order) deallocate(di, dj)

end subroutine load_xgrid

!#######################################################################
!
! get_grid - read the center point of the grid from grid_spec.nc.
!          - only the grid at the side 1 is needed, so we only read 
!          - atm and land grid
!
!

subroutine get_grid(grid, grid_id, grid_file, grid_version)
  type(grid_type), intent(inout) :: grid
  character(len=3), intent(in)   :: grid_id
  character(len=*), intent(in)   :: grid_file
  integer,          intent(in)   :: grid_version

  real, dimension(grid%im) :: lonb
  real, dimension(grid%jm) :: latb
  real, allocatable        :: tmpx(:,:), tmpy(:,:)
  real                     :: d2r
  integer                  :: is, ie, js, je, nlon, nlat, siz(4), i, j

  d2r = PI/180.0

  call mpp_get_compute_domain(grid%domain, is, ie, js, je)

  select case(grid_version)
  case(VERSION1)
     allocate(grid%lon(grid%im), grid%lat(grid%jm))
     if(grid_id == 'ATM') then
        call read_data(grid_file, 'xta', lonb, no_domain=.true.)
        call read_data(grid_file, 'yta', latb, no_domain=.true.)

        if(.not. allocated(AREA_ATM_MODEL)) then
           allocate(AREA_ATM_MODEL(is:ie, js:je))
           call get_area_elements(grid_file, 'AREA_ATM_MODEL', grid%domain, AREA_ATM_MODEL)
        endif
        if(.not. allocated(AREA_ATM_SPHERE)) then
           allocate(AREA_ATM_SPHERE(is:ie, js:je))
           call get_area_elements(grid_file, 'AREA_ATM', grid%domain, AREA_ATM_SPHERE)
        endif
     else if(grid_id == 'LND') then
        call read_data(grid_file, 'xtl', lonb, no_domain=.true.)
        call read_data(grid_file, 'ytl', latb, no_domain=.true.)
        if(.not. allocated(AREA_LND_MODEL)) then
           allocate(AREA_LND_MODEL(is:ie, js:je))
           call get_area_elements(grid_file, 'AREA_LND_MODEL', grid%domain, AREA_LND_MODEL)
        endif
        if(.not. allocated(AREA_LND_SPHERE)) then
           allocate(AREA_LND_SPHERE(is:ie, js:je))
           call get_area_elements(grid_file, 'AREA_LND', grid%domain, AREA_LND_SPHERE)
        endif
     else if(grid_id == 'OCN' ) then
        if(.not. allocated(AREA_OCN_SPHERE)) then
           allocate(AREA_OCN_SPHERE(is:ie, js:je))
           call get_area_elements(grid_file, 'AREA_OCN', grid%domain, AREA_OCN_SPHERE)
        endif
     endif
     !--- second order remapping suppose second order
     if(grid_id == 'LND' .or. grid_id == 'ATM') then
        grid%lon   = lonb * d2r
        grid%lat   = latb * d2r
     endif
     grid%is_latlon = .true.
  case(VERSION2)
     call field_size(grid_file, 'area', siz)
     nlon = siz(1); nlat = siz(2)
     if( mod(nlon,2) .NE. 0) call error_mesg('xgrid_mod',  &
          'flux_exchange_mod: atmos supergrid longitude size can not be divided by 2', FATAL)
     if( mod(nlat,2) .NE. 0) call error_mesg('xgrid_mod',  &
          'flux_exchange_mod: atmos supergrid latitude size can not be divided by 2', FATAL)
     nlon = nlon/2
     nlat = nlat/2
     if(nlon .NE. grid%im .OR. nlat .NE. grid%jm) call error_mesg('xgrid_mod', &
         'grid size in tile_file does not match the global grid size', FATAL)
     allocate(tmpx(nlon*2+1, nlat*2+1), tmpy(nlon*2+1, nlat*2+1))
     call read_data( grid_file, 'x', tmpx, no_domain=.true.)
     call read_data( grid_file, 'y', tmpy, no_domain=.true.) 
     if( grid_id == 'LND' .or. grid_id == 'ATM') then
        if(is_lat_lon(tmpx, tmpy) ) then
           allocate(grid%lon(grid%im), grid%lat(grid%jm))
           do i = 1, grid%im
              grid%lon(i) = tmpx(2*i,2) * d2r
           end do
           do j = 1, grid%jm
              grid%lat(j) = tmpy(2, 2*j) * d2r
           end do
           grid%is_latlon = .true.
        else
           allocate(grid%geolon(grid%isd_me:grid%ied_me, grid%jsd_me:grid%jed_me))
           allocate(grid%geolat(grid%isd_me:grid%ied_me, grid%jsd_me:grid%jed_me))
           grid%geolon = 1e10
           grid%geolat = 1e10
           !--- area_ocn_sphere, area_lnd_sphere, area_atm_sphere is not been defined.
           do j = grid%js_me,grid%je_me
              do i = grid%is_me,grid%ie_me
                 grid%geolon(i, j) = tmpx(i*2,j*2)*d2r
                 grid%geolat(i, j) = tmpy(i*2,j*2)*d2r
              end do
           end do
           call mpp_update_domains(grid%geolon, grid%domain)
           call mpp_update_domains(grid%geolat, grid%domain)
           grid%is_latlon = .false.
        end if
     end if
     deallocate(tmpx, tmpy)
  end select

  return

end subroutine get_grid
  
!#######################################################################
! Read the area elements from NetCDF file
subroutine get_area_elements(file, name, domain, data)
  character(len=*), intent(in) :: file
  character(len=*), intent(in) :: name
  type(domain2d),   intent(in) :: domain
  real, intent(out)            :: data(:,:)

  if(field_exist(file, name)) then
     call read_data(file, name, data, domain)
  else
     call error_mesg('xgrid_mod', 'no field named '//trim(name)//' in grid file '//trim(file)// &
                     ' Will set data to negative values...', NOTE)
     ! area elements no present in grid_spec file, set to negative values....
     data = -1
  endif    

end subroutine get_area_elements

!#######################################################################
! Read the OCN model area elements from NetCDF file
! <SUBROUTINE NAME="get_ocean_model_area_elements">

!   <OVERVIEW>
!      Read Ocean area element data.
!   </OVERVIEW>
!   <DESCRIPTION>
!      If available in the NetCDF file, this routine will read the 
!      AREA_OCN_MODEL field and load the data into global AREA_OCN_MODEL.
!      If not available, then the array AREA_OCN_MODEL will be left
!      unallocated. Must be called by all PEs.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_ocean_model_area_elements(ocean_domain, grid_file)
!   </TEMPLATE>

!   <IN NAME="ocean_domain" TYPE="type(Domain2d)"> </IN>
!   <IN NAME="grid_file" TYPE="character(len=*)" > </IN>
subroutine get_ocean_model_area_elements(domain, grid_file)

  type(Domain2d), intent(in) :: domain
  character(len=*), intent(in) :: grid_file
  integer :: is, ie, js, je

  if(allocated(AREA_OCN_MODEL)) return

  call mpp_get_compute_domain(domain, is, ie, js, je)
  ! allocate even if ie<is, ... in which case the array will have zero size
  ! but will still return .T. for allocated(...)
  allocate(AREA_OCN_MODEL(is:ie, js:je))
  if(ie < is .or. je < js ) return


  if(field_exist(grid_file, 'AREA_OCN_MODEL') )then
     call read_data(grid_file, 'AREA_OCN_MODEL', AREA_OCN_MODEL, domain)
  else
     deallocate(AREA_OCN_MODEL)
  endif


end subroutine get_ocean_model_area_elements
! </SUBROUTINE>
!#######################################################################

! <SUBROUTINE NAME="setup_xmap">

!   <OVERVIEW>
!      Sets up exchange grid connectivity using grid specification file and
!      processor domain decomposition. 
!   </OVERVIEW>
!   <DESCRIPTION>
!      Sets up exchange grid connectivity using grid specification file and
!      processor domain decomposition. Initializes xmap.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call setup_xmap(xmap, grid_ids, grid_domains, grid_file)
!   </TEMPLATE>

!   <IN NAME="grid_ids" TYPE="character(len=3)" DIM="(:)"> </IN>
!   <IN NAME="grid_domains" TYPE="type(Domain2d)" DIM="(:)"> </IN>
!   <IN NAME="grid_file" TYPE="character(len=*)" > </IN>
!   <OUT NAME="xmap" TYPE="xmap_type"  > </OUT>

subroutine setup_xmap(xmap, grid_ids, grid_domains, grid_file )
  type (xmap_type),                          intent(inout) :: xmap
  character(len=3), dimension(:),            intent(in ) :: grid_ids
!!$  type(Domain2d), dimension(size(grid_ids(:))), intent(in ) :: grid_domains
  type(Domain2d), dimension(:), intent(in ) :: grid_domains
  character(len=*)                         , intent(in ) :: grid_file

  integer :: g,     p, send_size, recv_size, i, siz(4)
  integer :: unit, nxgrid_file, i1, i2, i3, tile1, tile2, j
  type (grid_type), pointer, save :: grid =>NULL(), grid1 =>NULL()
  real, dimension(3) :: xxx
  real, dimension(:,:), allocatable   :: check_data
  real, dimension(:,:,:), allocatable :: check_data_3D
  integer, dimension(:),  allocatable :: nilist, njlist
  character(len=256)                  :: xgrid_file, xgrid_name
  character(len=256)                  :: tile_file, mosaic_file
  character(len=256)                  :: mosaic1, mosaic2, contact
  character(len=256)                  :: tile1_name, tile2_name
  character(len=256),     allocatable :: tile1_list(:), tile2_list(:)
  logical :: use_higher_order = .false.

  if(interp_method .ne. 'first_order')  use_higher_order = .true.

  xmap%me   = mpp_pe  ()
  xmap%npes = mpp_npes()
  xmap%root_pe = mpp_root_pe()

  allocate( xmap%grids(1:size(grid_ids(:))) )

  allocate ( xmap%your1my2(0:xmap%npes-1), xmap%your2my1(0:xmap%npes-1) )

  xmap%your1my2 = .false.; xmap%your2my1 = .false.;

!  check the exchange grid file version to be used by checking the field in the file
  if(field_exist(grid_file, "AREA_ATMxOCN" ) ) then
     xmap%version = VERSION1
  else if(field_exist(grid_file, "ocn_mosaic_file" ) ) then
     xmap%version = VERSION2
  else
     call error_mesg('xgrid_mod', 'both AREA_ATMxOCN and ocn_mosaic_file does not exist in '//trim(grid_file), FATAL)
  end if

  if(xmap%version==VERSION1) then
     call error_mesg('xgrid_mod', 'reading exchange grid information from grid spec file', NOTE)
  else
     call error_mesg('xgrid_mod', 'reading exchange grid information from mosaic grid file', NOTE)
  end if

  allocate ( nilist(0:xmap%npes-1), njlist(0:xmap%npes-1) )
  do g=1,size(grid_ids(:))
     grid => xmap%grids(g)
     if (g==1) grid1 => xmap%grids(g)
     grid%id     = grid_ids    (g)
     grid%domain = grid_domains(g)

     allocate ( grid%is(0:xmap%npes-1), grid%ie(0:xmap%npes-1) )
     allocate ( grid%js(0:xmap%npes-1), grid%je(0:xmap%npes-1) )
     allocate ( grid%tile(0:xmap%npes-1) )
     call mpp_get_compute_domains(grid%domain, xbegin=grid%is, xend=grid%ie, &
          ybegin=grid%js, yend=grid%je  )
     call mpp_get_global_domains(grid%domain, xsize=nilist, ysize=njlist)
     grid%ni = maxval(nilist)
     grid%nj = maxval(njlist)
     
     call mpp_get_tile_list(grid%domain, grid%tile)
     grid%ntile = mpp_get_ntile_count(grid%domain)
     ! make sure the grid%tile are between 1 and ntile 
     do p = 0, xmap%npes-1
        if(grid%tile(p) > grid%ntile .or. grid%tile(p) < 1) call error_mesg('xgrid_mod', &
                 'tile id should between 1 and ntile', FATAL)
     end do 

     grid%is_me => grid%is(xmap%me-xmap%root_pe); grid%ie_me => grid%ie(xmap%me-xmap%root_pe)
     grid%js_me => grid%js(xmap%me-xmap%root_pe); grid%je_me => grid%je(xmap%me-xmap%root_pe)
     grid%tile_me => grid%tile(xmap%me-xmap%root_pe)
     call mpp_get_data_domain(grid%domain, grid%isd_me, grid%ied_me, grid%jsd_me, grid%jed_me)
     if( use_higher_order ) then
        call mpp_modify_domain(grid%domain, grid%domain_with_halo, whalo=1, ehalo=1, shalo=1, nhalo=1)
        call mpp_get_data_domain(grid%domain_with_halo, grid%isd_me, grid%ied_me, grid%jsd_me, grid%jed_me) 
     end if

     !--- The starting index of compute domain may not start at 1.
     grid%im = maxval(grid%ie) - minval(grid%is) + 1
     grid%jm = maxval(grid%je) - minval(grid%js) + 1
     grid%km = 1

     allocate( grid%area    (grid%is_me:grid%ie_me, grid%js_me:grid%je_me) )
     allocate( grid%area_inv(grid%is_me:grid%ie_me, grid%js_me:grid%je_me) )
     grid%area       = 0.0
     grid%size       = 0
     grid%size_repro = 0

     ! get the center point of the grid box
     select case(xmap%version)
     case(VERSION1)
        call get_grid(grid, grid_ids(g), grid_file, xmap%version)
     case(VERSION2)
        call read_data(grid_file, lowercase(grid_ids(g))//'_mosaic_file', mosaic_file)      
        call get_mosaic_tile_grid(tile_file, 'INPUT/'//trim(mosaic_file), grid%domain)
        call get_grid(grid, grid_ids(g), tile_file, xmap%version)
        if(use_higher_order .AND. grid%id == 'ATM' .AND. (.NOT.grid%is_latlon) ) call error_mesg( &
            'xgrid_mod', 'higher order flux exchange is not implemented for non-latlon grid', FATAL)
     end select

     if (g>1) then
        allocate( grid%frac_area(grid%is_me:grid%ie_me, grid%js_me:grid%je_me, grid%km) )
        grid%frac_area = 1.0

        ! load exchange cells, sum grid cell areas, set your1my2/your2my1
        select case(xmap%version)
        case(VERSION1)
           call load_xgrid (xmap, grid, grid_file, grid_ids(1), grid_ids(g), 1, 1, use_higher_order, .true. )
        case(VERSION2)
           select case(grid_ids(1))
           case( 'ATM' )
              xgrid_name = 'a'
           case( 'LND' )
              xgrid_name = 'l'
           case default 
              call error_mesg('xgrid_mod', 'grid_ids(1) should be ATM or LND', FATAL)
           end select
           select case(grid_ids(g))
           case( 'LND' )
              xgrid_name = trim(xgrid_name)//'Xl_file'
           case( 'OCN' )
              xgrid_name = trim(xgrid_name)//'Xo_file'
           case default 
              call error_mesg('xgrid_mod', 'grid_ids(g) should be LND or OCN', FATAL)
           end select       
           ! get the tile list for each mosaic
           call read_data(grid_file, lowercase(grid_ids(1))//'_mosaic_file', mosaic1) 
           call read_data(grid_file, lowercase(grid_ids(g))//'_mosaic_file', mosaic2) 
           mosaic1 = 'INPUT/'//trim(mosaic1)
           mosaic2 = 'INPUT/'//trim(mosaic2)
           allocate(tile1_list(grid1%ntile), tile2_list(grid%ntile) )
           do j = 1, grid1%ntile
              call read_data(mosaic1, 'gridtiles', tile1_list(j), level=j)
           end do
           do j = 1, grid%ntile
              call read_data(mosaic2, 'gridtiles', tile2_list(j), level=j)
           end do
           call field_size(grid_file, xgrid_name, siz)
           nxgrid_file = siz(2)
           ! loop through all the exchange grid file
           do i = 1, nxgrid_file
              call read_data(grid_file, xgrid_name, xgrid_file, level = i)
              xgrid_file = 'INPUT/'//trim(xgrid_file) 
              if( .NOT. file_exist(xgrid_file) )call error_mesg('xgrid_mod', &
                   'file '//trim(xgrid_file)//' does not exist, check your xgrid file.', FATAL)
              
              ! find the tile number of side 1 and side 2 mosaic, which is contained in field contact
              call read_data(xgrid_file, "contact", contact)
              i1 = index(contact, ":")
              i2 = index(contact, "::")
              i3 = index(contact, ":", back=.true. )
              if(i1 == 0 .OR. i2 == 0) call error_mesg('xgrid_mod', &
                   'field contact in file '//trim(xgrid_file)//' should contains ":" and "::" ', FATAL)
              if(i1 == i3) call error_mesg('xgrid_mod', &
                   'field contact in file '//trim(xgrid_file)//' should contains two ":"', FATAL)
              tile1_name = contact(i1+1:i2-1)
              tile2_name = contact(i3+1:len_trim(contact))
              tile1 = 0; tile2 = 0
              do j = 1, grid1%ntile
                 if( tile1_name == tile1_list(j) ) then
                    tile1 = j
                    exit
                 end if
              end do
              do j = 1, grid%ntile
                 if( tile2_name == tile2_list(j) ) then
                    tile2 = j
                    exit
                 end if
              end do
              if(tile1 == 0) call error_mesg('xgrid_mod', &
                   trim(tile1_name)//' is not a tile of mosaic '//trim(mosaic1), FATAL)
              if(tile2 == 0) call error_mesg('xgrid_mod', &
                   trim(tile2_name)//' is not a tile of mosaic '//trim(mosaic2), FATAL)
              call load_xgrid (xmap, grid, xgrid_file, grid_ids(1), grid_ids(g), tile1, tile2, use_higher_order, i==nxgrid_file)
           end do
           deallocate(tile1_list, tile2_list)
        end select
     end if
  end do

  deallocate(nilist, njlist)
  grid1%area_inv = 0.0;
  where (grid1%area>0.0)
     grid1%area_inv = 1.0/grid1%area
  end where

  xmap%your1my2(xmap%me-xmap%root_pe) = .false. ! this is not necessarily true but keeps
  xmap%your2my1(xmap%me-xmap%root_pe) = .false. ! a PE from communicating with itself

  send_size = sum((grid1%ie-grid1%is+1)*(grid1%je-grid1%js+1))
  send_size = max(send_size, grid1%im*grid1%jm)
  recv_size = maxval((grid1%ie-grid1%is+1)*(grid1%je-grid1%js+1) )
  if (make_exchange_reproduce) then
     allocate( xmap%send_count_repro(0:xmap%npes-1) )
     allocate( xmap%recv_count_repro(0:xmap%npes-1) )
     xmap%send_count_repro = 0
     xmap%recv_count_repro = 0
     do g=2,size(xmap%grids(:))
        do p=0,xmap%npes-1
           if(xmap%grids(g)%size >0) &
                xmap%send_count_repro(p) = xmap%send_count_repro(p) &
                +count(xmap%grids(g)%x      (:)%pe==p+xmap%root_pe)
           if(xmap%grids(g)%size_repro >0) &
                xmap%recv_count_repro(p) = xmap%recv_count_repro(p) &
                +count(xmap%grids(g)%x_repro(:)%pe==p+xmap%root_pe)
        end do
     end do
     send_size = max(send_size, sum(xmap%send_count_repro))
  end if
  allocate (xmap%send_buffer(send_size))
  allocate (xmap%recv_buffer(recv_size))

  call mpp_open( unit, 'xgrid.out', action=MPP_OVERWR, threading=MPP_MULTI, &
       fileset=MPP_SINGLE, nohdrs=.TRUE. )  

  write( unit,* )xmap%grids(:)%id, ' GRID: PE ', xmap%me, ' #XCELLS=', &
       xmap%grids(2:size(xmap%grids(:)))%size, ' #COMM. PARTNERS=', &
       count(xmap%your1my2), '/', count(xmap%your2my1), &
       pack((/(p+xmap%root_pe,p=0,xmap%npes-1)/), xmap%your1my2),  &
       '/', pack((/(p+xmap%root_pe,p=0,xmap%npes-1)/), xmap%your2my1)

  allocate( xmap%x1(1:sum(xmap%grids(2:size(xmap%grids(:)))%size)) )
  allocate( xmap%x2(1:sum(xmap%grids(2:size(xmap%grids(:)))%size)) )

  call regen(xmap)

  xxx = conservation_check(grid1%area*0+1.0, grid1%id, xmap)
  write(stdout(),* )"Checked data is array of constant 1"
  write(stdout(),* )grid1%id,'(',xmap%grids(:)%id,')=', xxx 

  do g=2,size(xmap%grids(:))
     xxx = conservation_check(xmap%grids(g)%frac_area*0+1.0, xmap%grids(g)%id, xmap )
     write( stdout(),* )xmap%grids(g)%id,'(',xmap%grids(:)%id,')=', xxx 
  enddo
  ! create an random number 2d array
  if(grid1%id == "ATM") then
     allocate(check_data(size(grid1%area,1), size(grid1%area,2)))
     call random_number(check_data)

     !--- second order along both zonal and meridinal direction
     xxx = conservation_check(check_data, grid1%id, xmap,  remap_method = remapping_method )
     write( stdout(),* ) &
          "Checked data is array of random number between 0 and 1 using "//trim(interp_method)
     write( stdout(),* )grid1%id,'(',xmap%grids(:)%id,')=', xxx 

     deallocate(check_data)
     do g=2,size(xmap%grids(:))
        allocate(check_data_3d(size(xmap%grids(g)%frac_area,1),size(xmap%grids(g)%frac_area,2), &
             size(xmap%grids(g)%frac_area,3) )) 
        call random_number(check_data_3d)
        xxx = conservation_check(check_data_3d, xmap%grids(g)%id, xmap,  remap_method = remapping_method )
        write( stdout(),* )xmap%grids(g)%id,'(',xmap%grids(:)%id,')=', xxx
        deallocate( check_data_3d)
     end do
  endif

  call close_file (unit)

end subroutine setup_xmap
! </SUBROUTINE>

!#######################################################################


subroutine regen(xmap)
type (xmap_type), intent(inout) :: xmap

  integer :: g, l, i, j, k, max_size

  max_size = 0;
  do g=2,size(xmap%grids(:))
    max_size = max_size + xmap%grids(g)%size * xmap%grids(g)%km
  end do
  if (max_size>size(xmap%x1(:))) then
    deallocate(xmap%x1)
    deallocate(xmap%x2)
    allocate( xmap%x1(1:max_size) )
    allocate( xmap%x2(1:max_size) )
  end if

  xmap%size = 0
  do g=2,size(xmap%grids(:))
    xmap%grids(g)%first = xmap%size + 1;
    do l=1,xmap%grids(g)%size
      i = xmap%grids(g)%x(l)%i2
      j = xmap%grids(g)%x(l)%j2
      do k=1,xmap%grids(g)%km
        if (xmap%grids(g)%frac_area(i,j,k)/=0.0) then
          xmap%size = xmap%size+1
          xmap%x1(xmap%size)%i    = xmap%grids(g)%x(l)%i1
          xmap%x1(xmap%size)%j    = xmap%grids(g)%x(l)%j1
          xmap%x1(xmap%size)%tile = xmap%grids(g)%x(l)%tile
          xmap%x1(xmap%size)%area = xmap%grids(g)%x(l)%area &
                                   *xmap%grids(g)%frac_area(i,j,k)
          xmap%x1(xmap%size)%di   = xmap%grids(g)%x(l)%di 
          xmap%x1(xmap%size)%dj   = xmap%grids(g)%x(l)%dj 
          xmap%x2(xmap%size)%i    = xmap%grids(g)%x(l)%i2
          xmap%x2(xmap%size)%j    = xmap%grids(g)%x(l)%j2
          xmap%x2(xmap%size)%k    = k
          xmap%x2(xmap%size)%area = xmap%grids(g)%x(l)%area
        end if
      end do
    end do
    xmap%grids(g)%last = xmap%size
  end do
end subroutine regen

!#######################################################################

! <SUBROUTINE NAME="set_frac_area">

!   <OVERVIEW>
!     Changes sub-grid portion areas and/or number.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Changes sub-grid portion areas and/or number.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call set_frac_area(f, grid_id, xmap)
!   </TEMPLATE>

!   <IN NAME="f" TYPE="real" DIM="(:,:,:)"> </IN>
!   <IN NAME="grid_id" TYPE="character(len=3)" > </IN>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>

subroutine set_frac_area(f, grid_id, xmap)
real, dimension(:,:,:), intent(in   ) :: f
character(len=3),       intent(in   ) :: grid_id
type (xmap_type),       intent(inout) :: xmap

  integer :: g
  type(grid_type), pointer, save :: grid =>NULL()

  if (grid_id==xmap%grids(1)%id) call error_mesg ('xgrid_mod',  &
                                   'set_frac_area called on side 1 grid', FATAL)
  do g=2,size(xmap%grids(:))
    grid => xmap%grids(g)
    if (grid_id==grid%id) then
      if (size(f,3)/=size(grid%frac_area,3)) then
        deallocate (grid%frac_area)
        grid%km = size(f,3);
        allocate( grid%frac_area(grid%is_me:grid%ie_me, grid%js_me:grid%je_me, &
                                                                      grid%km) )
      end if
      grid%frac_area = f;
      call regen(xmap)
      return;
    end if
  end do

  call error_mesg ('xgrid_mod', 'set_frac_area: could not find grid id', FATAL)

end subroutine  set_frac_area
! </SUBROUTINE>



!#######################################################################

! <FUNCTION NAME="xgrid_count">

!   <OVERVIEW>
!     Returns current size of exchange grid variables.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns current size of exchange grid variables.
!   </DESCRIPTION>
!   <TEMPLATE>
!     xgrid_count(xmap)
!   </TEMPLATE>

!   <IN NAME="xmap" TYPE="xmap_type" > </IN>
!   <OUT NAME="xgrid_count"  TYPE="integer"  > </OUT>

integer function xgrid_count(xmap)
type (xmap_type), intent(inout) :: xmap

  xgrid_count = xmap%size
end function xgrid_count
! </FUNCTION>

!#######################################################################

! <SUBROUTINE NAME="put_side1_to_xgrid" INTERFACE="put_to_xgrid">
!   <IN NAME="d"  TYPE="real" DIM="(:,:)" > </IN>
!   <IN NAME="grid_id"  TYPE=" character(len=3)"  > </IN>
!   <INOUT NAME="x"  TYPE="real" DIM="(:)" > </INOUT>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>
!   <IN NAME="remap_method" TYPE="integer,optional"></IN>

subroutine put_side1_to_xgrid(d, grid_id, x, xmap, remap_method)
real, dimension(:,:), intent(in   )    :: d
character(len=3),     intent(in   )    :: grid_id
real, dimension(:),   intent(inout)    :: x
type (xmap_type),     intent(inout)    :: xmap
integer, intent(in), optional          :: remap_method

  integer :: g, method

  method = FIRST_ORDER      ! default
  if(present(remap_method)) method = remap_method

  if (grid_id==xmap%grids(1)%id) then
       if(method == FIRST_ORDER) then
          call put_1_to_xgrid_order_1(d, x, xmap)
       else 
          if(grid_id .NE. 'ATM') call error_mesg ('xgrid_mod',  &
                       "second order put_to_xgrid should only be applied to 'ATM' model, "//&
                       "contact developer", FATAL)
          call put_1_to_xgrid_order_2(d, x, xmap, method )
       endif
    return;
  end if

  do g=2,size(xmap%grids(:))
    if (grid_id==xmap%grids(g)%id)    &
      call error_mesg ('xgrid_mod',  &
                       'put_to_xgrid expects a 3D side 2 grid', FATAL)
  end do

  call error_mesg ('xgrid_mod', 'put_to_xgrid: could not find grid id', FATAL)

end subroutine put_side1_to_xgrid
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="put_side2_to_xgrid" INTERFACE="put_to_xgrid">
!   <IN NAME="d"  TYPE="real" DIM="(:,:,:)" > </IN>
!   <IN NAME="grid_id"  TYPE=" character(len=3)"  > </IN>
!   <INOUT NAME="x"  TYPE="real" DIM="(:)" > </INOUT>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>

subroutine put_side2_to_xgrid(d, grid_id, x, xmap)
real, dimension(:,:,:), intent(in   ) :: d
character(len=3),       intent(in   ) :: grid_id
real, dimension(:),     intent(inout) :: x
type (xmap_type),       intent(inout) :: xmap

  integer :: g

  if (grid_id==xmap%grids(1)%id) &
    call error_mesg ('xgrid_mod',  &
                     'put_to_xgrid expects a 2D side 1 grid', FATAL)

  do g=2,size(xmap%grids(:))
    if (grid_id==xmap%grids(g)%id) then
         call put_2_to_xgrid(d, xmap%grids(g), x, xmap)
      return;
    end if
  end do

  call error_mesg ('xgrid_mod', 'put_to_xgrid: could not find grid id', FATAL)

end subroutine put_side2_to_xgrid
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="get_side1_from_xgrid" INTERFACE="get_from_xgrid">
!   <IN NAME="x"  TYPE="real" DIM="(:)" > </IN>
!   <IN NAME="grid_id"  TYPE=" character(len=3)"  > </IN>
!   <OUT NAME="d"  TYPE="real" DIM="(:,:)" > </OUT>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>

subroutine get_side1_from_xgrid(d, grid_id, x, xmap)
real, dimension(:,:), intent(  out) :: d
character(len=3),     intent(in   ) :: grid_id
real, dimension(:),   intent(in   ) :: x
type (xmap_type),     intent(inout) :: xmap

  integer :: g

  if (grid_id==xmap%grids(1)%id) then
    if (make_exchange_reproduce) then
      call get_1_from_xgrid_repro(d, x, xmap)
    else
      call get_1_from_xgrid(d, x, xmap)
    end if
    return;
  end if
  
  do g=2,size(xmap%grids(:))
    if (grid_id==xmap%grids(g)%id) &
      call error_mesg ('xgrid_mod',  & 
                       'get_from_xgrid expects a 3D side 2 grid', FATAL)
  end do
  
  call error_mesg ('xgrid_mod', 'get_from_xgrid: could not find grid id', FATAL)

end subroutine get_side1_from_xgrid
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="get_side2_from_xgrid" INTERFACE="get_from_xgrid">
!   <IN NAME="x"  TYPE="real" DIM="(:)" > </IN>
!   <IN NAME="grid_id"  TYPE=" character(len=3)"  > </IN>
!   <OUT NAME="d"  TYPE="real" DIM="(:,:,:)" > </OUT>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>

subroutine get_side2_from_xgrid(d, grid_id, x, xmap)
real, dimension(:,:,:), intent(  out) :: d
character(len=3),       intent(in   ) :: grid_id
real, dimension(:),     intent(in   ) :: x
type (xmap_type),       intent(in   ) :: xmap

  integer :: g

  if (grid_id==xmap%grids(1)%id) &
    call error_mesg ('xgrid_mod',  &
                     'get_from_xgrid expects a 2D side 1 grid', FATAL)
  
  do g=2,size(xmap%grids(:))
    if (grid_id==xmap%grids(g)%id) then
      call get_2_from_xgrid(d, xmap%grids(g), x, xmap)
      return;
    end if
  end do
  
  call error_mesg ('xgrid_mod', 'get_from_xgrid: could not find grid id', FATAL)

end subroutine get_side2_from_xgrid
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="some">

!   <OVERVIEW>
!     Returns logical associating exchange grid cells with given side two grid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns logical associating exchange grid cells with given side two grid.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call some(xmap, some_arr, grid_id)
!   </TEMPLATE>

!   <IN NAME="xmap"  TYPE="xmap_type"  ></IN>
!   <IN NAME="grid_id"  TYPE="character(len=3)"  ></IN>
!   <OUT NAME="some_arr"  TYPE="logical" DIM="(xmap%size)" >
!     logical associating exchange grid cells with given side 2 grid.
!   </OUT>

subroutine some(xmap, some_arr, grid_id)
type (xmap_type),           intent(in) :: xmap
character(len=3), optional, intent(in) :: grid_id
logical, dimension(xmap%size), intent(out) :: some_arr

  integer :: g

  if (.not.present(grid_id)) then
    some_arr = .true.
    return;
  end if

  if (grid_id==xmap%grids(1)%id) &
    call error_mesg ('xgrid_mod', 'some expects a side 2 grid id', FATAL)
  
  do g=2,size(xmap%grids(:))
    if (grid_id==xmap%grids(g)%id) then
      some_arr = .false.
      some_arr(xmap%grids(g)%first:xmap%grids(g)%last) = .true.;
      return;
    end if
  end do
  
  call error_mesg ('xgrid_mod', 'some could not find grid id', FATAL)

end subroutine some
! </SUBROUTINE>

!#######################################################################

subroutine put_2_to_xgrid(d, grid, x, xmap)
type (grid_type),                                intent(in) :: grid
real, dimension(grid%is_me:grid%ie_me, &
                grid%js_me:grid%je_me, grid%km), intent(in) :: d
real, dimension(:    ), intent(inout) :: x
type (xmap_type),       intent(in   ) :: xmap

  integer                 ::   l

  do l=grid%first,grid%last
    x(l) = d(xmap%x2(l)%i,xmap%x2(l)%j,xmap%x2(l)%k)
  end do
end subroutine put_2_to_xgrid

!#######################################################################

subroutine get_2_from_xgrid(d, grid, x, xmap)
type (grid_type),                                intent(in ) :: grid
real, dimension(grid%is_me:grid%ie_me, &
                grid%js_me:grid%je_me, grid%km), intent(out) :: d
real, dimension(:),     intent(in   ) :: x
type (xmap_type),       intent(in   ) :: xmap

  integer                 :: l, k

  d = 0.0
  do l=grid%first,grid%last
    d(xmap%x2(l)%i,xmap%x2(l)%j,xmap%x2(l)%k) = &
            d(xmap%x2(l)%i,xmap%x2(l)%j,xmap%x2(l)%k) + xmap%x2(l)%area*x(l)
  end do
  !
  !  normalize with side 2 grid cell areas
  !
  do k=1,size(d,3)
    d(:,:,k) = d(:,:,k) * grid%area_inv
  end do
end subroutine get_2_from_xgrid

!#######################################################################

function get_side_1(pe, im, jm)
integer, intent(in)    :: pe, im, jm
real, dimension(im,jm) :: get_side_1




!  call mpp_recv(buf, im*jm, pe)
!  l = 0
!  do j=1,jm; do i=1,im;
!    l = l + 1
!    get_side_1(i,j) = buf(l)
!  end do; end do
  ! Force use of "scalar", integer pointer mpp interface.
  call mpp_recv( get_side_1(1,1), glen=im*jm, from_pe=pe )
end function get_side_1

!#######################################################################

subroutine put_1_to_xgrid_order_1(d, x, xmap)
real, dimension(:,:), intent(in   ) :: d
real, dimension(:  ), intent(inout) :: x
type (xmap_type),     intent(inout) :: xmap

  integer :: i, is, ie, im, j, js, je, jm, p, l, tile
  real, dimension(xmap%grids(1)%ni,xmap%grids(1)%nj,xmap%grids(1)%ntile) :: dg
  type (grid_type), pointer, save :: grid1 =>NULL()

  grid1 => xmap%grids(1)
  is = grid1%is_me; ie = grid1%ie_me;
  js = grid1%js_me; je = grid1%je_me;
  tile = grid1%tile_me
  dg(is:ie,js:je,tile) = d;

  im = ie-is+1; jm = je-js+1;
  l = 0
  call mpp_sync_self()          !Balaji
  do j=1,jm; do i=1,im;
    l = l + 1;
    xmap%send_buffer(l) =  d(i,j)
  end do; end do;
  do p=0,xmap%npes-1
    if (xmap%your2my1(p)) then
      ! Force use of "scalar", integer pointer mpp interface.
      call mpp_send(xmap%send_buffer(1), plen=im*jm, to_pe=p+xmap%root_pe);
    end if
  end do
  do p=0,xmap%npes-1
    if (xmap%your1my2(p)) then
      is = grid1%is(p); ie = grid1%ie(p);
      js = grid1%js(p); je = grid1%je(p);
      tile = grid1%tile(p)
      dg(is:ie,js:je,tile) = get_side_1(p+xmap%root_pe,ie-is+1,je-js+1);
    end if
  end do
  do l=1,xmap%size
    x(l) =  dg(xmap%x1(l)%i,xmap%x1(l)%j,xmap%x1(l)%tile)
  end do

!  call mpp_sync_self
end subroutine put_1_to_xgrid_order_1

!#######################################################################


subroutine put_1_to_xgrid_order_2(d, x, xmap, remap_method)
  real, dimension(:,:), intent(in   ) :: d
  real, dimension(:  ), intent(inout) :: x
  type (xmap_type),     intent(inout) :: xmap
  integer,              intent(in)    :: remap_method

  integer :: i, is, ie, im, j, js, je, jm, p, l, isd, jsd, tile
  real, dimension(xmap%grids(1)%im,xmap%grids(1)%jm,xmap%grids(1)%ntile) :: dg
  real, dimension(xmap%grids(1)%im,xmap%grids(1)%jm,xmap%grids(1)%ntile) :: grad_x, grad_y
  real, dimension(xmap%grids(1)%isd_me:xmap%grids(1)%ied_me,xmap%grids(1)%jsd_me:xmap%grids(1)%jed_me) :: tmp
  type (grid_type), pointer, save :: grid1 =>NULL()
  integer        :: num_block, send_size, recv_size

  grid1 => xmap%grids(1)
  is = grid1%is_me;   ie = grid1%ie_me
  js = grid1%js_me;   je = grid1%je_me
  isd = grid1%isd_me
  jsd = grid1%jsd_me
  im = ie-is+1;       jm = je-js+1
  tile = grid1%tile_me
  dg(is:ie,js:je,tile) = d

  ! first get the halo of data
  tmp(is:ie,js:je) = d(:,:)

  if (remap_method == SECOND_ORDER_ZONAL) then
     call mpp_update_domains(tmp,grid1%domain_with_halo, flags=XUPDATE)
  else if ( remap_method == SECOND_ORDER_MERID ) then
     call mpp_update_domains(tmp,grid1%domain_with_halo, flags=YUPDATE)
  else
     call mpp_update_domains(tmp,grid1%domain_with_halo)
  endif

  num_block = 1
  if( remap_method .ne. SECOND_ORDER_ZONAL ) then   ! second_order_merid or second_order
     if(grid1%is_latlon) then
       grad_y(is:ie,js:je,tile) = grad_merid(tmp, grid1%lat, is, ie, js, je, isd, jsd)
     else
       grad_y(is:ie,js:je,tile) = grad_merid(tmp, grid1%geolat, is, ie, js, je, isd, jsd)
     end if
     num_block = num_block + 1
  endif

  if (remap_method .ne. SECOND_ORDER_MERID) then ! second_order_zonal or second_order
     if(grid1%is_latlon) then
        grad_x(is:ie,js:je,tile) = grad_zonal(tmp, grid1%lon, grid1%lat, is, ie, js, je, isd, jsd)
     else
        grad_x(is:ie,js:je,tile) = grad_zonal(tmp, grid1%geolon, grid1%geolat, is, ie, js, je, isd, jsd)
     end if
     num_block = num_block + 1
  endif

  call mpp_sync_self()          !Balaji

  send_size = num_block*im*jm
  ! if size of send_buffer is not enough, need to reallocate send_buffer
  if(size(xmap%send_buffer(:)) .lt. send_size) then
     deallocate(xmap%send_buffer)
     allocate(xmap%send_buffer(send_size))
  endif

  l = 0
  do j=js,je; do i=is,ie
     l = l + 1
     xmap%send_buffer(l) =  tmp(i,j)
  end do; end do

  if(remap_method .ne. SECOND_ORDER_ZONAL) then
     do j=js,je; do i=is,ie
        l = l + 1
        xmap%send_buffer(l) = grad_y(i,j,tile)
     end do; end do
  endif

  if (remap_method .ne. SECOND_ORDER_MERID) then
     do j=js,je; do i=is,ie
        l = l + 1
        xmap%send_buffer(l) = grad_x(i,j,tile)
     end do; end do
  endif

  do p=0,xmap%npes-1
     if (xmap%your2my1(p)) then
        ! Force use of "scalar", integer pointer mpp interface.
        call mpp_send(xmap%send_buffer(1), plen=send_size, to_pe=p+xmap%root_pe);
     end if
  end do

  do p=0,xmap%npes-1
     if (xmap%your1my2(p)) then
        is = grid1%is(p);  ie = grid1%ie(p)
        js = grid1%js(p);  je = grid1%je(p)
        tile = grid1%tile(p)
        recv_size = num_block*(ie-is+1)*(je-js+1)
        if(size(xmap%recv_buffer(:)) .lt. recv_size) then
           deallocate(xmap%recv_buffer)
           allocate(xmap%recv_buffer(recv_size))
        endif
        call mpp_recv(xmap%recv_buffer(1), glen = recv_size, from_pe = p+xmap%root_pe)
        l = 0
        do j = js,je; do i=is,ie
           l = l + 1
           dg(i,j,tile) = xmap%recv_buffer(l)
        enddo; enddo
        if(remap_method .ne. SECOND_ORDER_ZONAL) then
           do j = js,je; do i=is,ie
              l = l + 1
              grad_y(i,j,tile) = xmap%recv_buffer(l)
           enddo; enddo
        endif

        if (remap_method .ne. SECOND_ORDER_MERID) then
           do j = js,je; do i=is,ie
              l = l + 1
              grad_x(i,j,tile) = xmap%recv_buffer(l)
           enddo; enddo
        endif
     end if
  end do

  do l=1,xmap%size
     tile = xmap%x1(l)%tile
     x(l) =  dg(xmap%x1(l)%i,xmap%x1(l)%j,tile)
     if(remap_method .ne. SECOND_ORDER_ZONAL) then
        x(l) = x(l) + grad_y(xmap%x1(l)%i,xmap%x1(l)%j,tile ) *xmap%x1(l)%dj
     endif
     if (remap_method .ne. SECOND_ORDER_MERID) then
        x(l) = x(l) + grad_x(xmap%x1(l)%i,xmap%x1(l)%j,tile ) *xmap%x1(l)%di
     endif
  end do

end subroutine put_1_to_xgrid_order_2

!#######################################################################

subroutine get_1_from_xgrid(d, x, xmap)
real, dimension(:,:), intent(out)   :: d
real, dimension(:  ), intent(in )   :: x
type (xmap_type),     intent(inout) :: xmap

  real, dimension(xmap%grids(1)%ni,xmap%grids(1)%nj,xmap%grids(1)%ntile), target :: dg
  integer :: i, is, ie, im, j, js, je, jm, l, le, p, tile
  real             , pointer, save :: dgp =>NULL()
  type (grid_type) , pointer, save :: grid1 =>NULL()

  grid1 => xmap%grids(1)

  dg = 0.0;
  do l=1,xmap%size
    dgp => dg(xmap%x1(l)%i,xmap%x1(l)%j,xmap%x1(l)%tile)
    dgp =  dgp + xmap%x1(l)%area*x(l)
  end do

  le = 0;
  call mpp_sync_self()          !Balaji
  do p=0,xmap%npes-1
    if (xmap%your1my2(p)) then
      l = le + 1;
      is = grid1%is(p); ie = grid1%ie(p);
      js = grid1%js(p); je = grid1%je(p);
      tile = grid1%tile(p)
      do j=js,je; do i=is,ie;
        le = le + 1
        xmap%send_buffer(le) = dg(i,j,tile)
      end do; end do;
      ! Force use of "scalar", integer pointer mpp interface.
      call mpp_send(xmap%send_buffer(l), plen=le-l+1, to_pe=p+xmap%root_pe);
    end if
  end do
  d = dg(grid1%is_me:grid1%ie_me,grid1%js_me:grid1%je_me,grid1%tile_me);
  im = grid1%ie_me-grid1%is_me+1;
  jm = grid1%je_me-grid1%js_me+1;
  do p=0,xmap%npes-1
    if (xmap%your2my1(p)) d = d + get_side_1(p+xmap%root_pe,im,jm)
  end do
  !
  ! normalize with side 1 grid cell areas
  !
  d = d * grid1%area_inv

!  call mpp_sync_self
end subroutine get_1_from_xgrid

!#######################################################################

subroutine get_1_from_xgrid_repro(d, x, xmap)
type (xmap_type), intent(inout)                  :: xmap
real, dimension(xmap%grids(1)%is_me:xmap%grids(1)%ie_me, &
                xmap%grids(1)%js_me:xmap%grids(1)%je_me), intent(out) :: d
real, dimension(:  ), intent(in ) :: x

  real,    dimension(:), allocatable :: x_psum
  integer, dimension(:), allocatable :: pe_psum
  integer :: l1, l2, l3, g, i, j, k, p
  integer, dimension(0:xmap%npes-1) :: pl
  type (grid_type), pointer, save :: grid =>NULL()

  allocate ( x_psum  (sum(xmap%send_count_repro)) )
  allocate ( pe_psum (sum(xmap%send_count_repro)) )
  x_psum = 0.0
  l1 = 0 ! index into partition summed exchange grid variable
  l2 = 0 ! index into exchange grid variable
  do g=2,size(xmap%grids(:))
    do l3=1,xmap%grids(g)%size ! index into this side 2 grid's patterns
      l1 = l1 + 1
      do k=1,xmap%grids(g)%km
        i = xmap%grids(g)%x(l3)%i2
        j = xmap%grids(g)%x(l3)%j2
        if (xmap%grids(g)%frac_area(i,j,k)/=0.0) then
          l2 = l2 + 1
          x_psum (l1) = x_psum(l1) + xmap%x1(l2)%area * x(l2)
          pe_psum(l1) = xmap%grids(g)%x(l3)%pe
        end if
      end do
    end do
  end do
  l2 = 0;
  call mpp_sync_self()          !Balaji
  do p=0,xmap%npes-1
    l1 = l2 + 1
    l2 = l2 + xmap%send_count_repro(p)
    if (xmap%send_count_repro(p)>0) then ! can send to myself
      xmap%send_buffer(l1:l2) = pack(x_psum, pe_psum==p+xmap%root_pe)
      ! Force use of "scalar", integer pointer mpp interface.
      call mpp_send(xmap%send_buffer(l1), plen=l2-l1+1, to_pe=p+xmap%root_pe);
    end if
  end do
  deallocate ( x_psum, pe_psum)
  allocate ( x_psum (sum(xmap%recv_count_repro)) )
  l2 = 0;
  do p=0,xmap%npes-1
    l1 = l2 + 1
    l2 = l2 + xmap%recv_count_repro(p)
    if (xmap%recv_count_repro(p)>0) then ! can receive from myself
      ! Force use of "scalar", integer pointer mpp interface.
      call mpp_recv(x_psum(l1), glen=l2-l1+1, from_pe=p+xmap%root_pe);
      pl(p) = l1
    end if
  end do
  d = 0.0
  do g=2,size(xmap%grids(:))
    grid => xmap%grids(g)
    do l3=1,grid%size_repro ! index into side1 grid's patterns
      i = grid%x_repro(l3)%i1
      j = grid%x_repro(l3)%j1
      d(i,j) = d(i,j) + x_psum(pl(grid%x_repro(l3)%pe-xmap%root_pe))
      pl(grid%x_repro(l3)%pe-xmap%root_pe) = pl(grid%x_repro(l3)%pe-xmap%root_pe) + 1
    end do
  end do
  deallocate ( x_psum )
  !
  ! normalize with side 1 grid cell areas
  !
  d = d * xmap%grids(1)%area_inv

!  call mpp_sync_self
end subroutine get_1_from_xgrid_repro

!#######################################################################

! <FUNCTION NAME="conservation_check_side1" INTERFACE="conservation_check">
!   <IN NAME="d"  TYPE="real" DIM="(:,:)" > </IN>
!   <IN NAME="grid_id"  TYPE="character(len=3)"  > </IN>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>
!   <OUT NAME="conservation_check_side1" TYPE="real" DIM="dimension(3)" > </OUT>
!   <IN NAME="remap_method" TYPE="integer,optional"></IN>
! conservation_check - returns three numbers which are the global sum of a
! variable (1) on its home model grid, (2) after interpolation to the other
! side grid(s), and (3) after re_interpolation back onto its home side grid(s).
!
function conservation_check_side1(d, grid_id, xmap,remap_method) ! this one for 1->2->1
real, dimension(:,:),    intent(in   ) :: d
character(len=3),        intent(in   ) :: grid_id
type (xmap_type),        intent(inout) :: xmap
real, dimension(3)                     :: conservation_check_side1
integer, intent(in), optional :: remap_method


  real, dimension(xmap%size) :: x_over, x_back
  real, dimension(size(d,1),size(d,2)) :: d1
  real, dimension(:,:,:), allocatable  :: d2
  integer                              :: g
  type (grid_type), pointer, save      :: grid1 =>NULL(), grid2 =>NULL()

  grid1 => xmap%grids(1)
  conservation_check_side1(1) = mpp_global_sum(grid1%domain, grid1%area*d)
  conservation_check_side1(2) = 0.0
  call put_to_xgrid (d, grid1%id, x_over, xmap, remap_method)    ! put from side 1
  do g=2,size(xmap%grids(:))
    grid2 => xmap%grids(g)
    allocate (d2 (grid2%is_me:grid2%ie_me, grid2%js_me:grid2%je_me,  grid2%km) )
    call get_from_xgrid (d2, grid2%id, x_over, xmap) ! get onto side 2's
    conservation_check_side1(2) = conservation_check_side1(2) + &
      mpp_global_sum( grid2%domain, grid2%area * sum(grid2%frac_area*d2,DIM=3) )
    call put_to_xgrid (d2, grid2%id, x_back, xmap) ! put from side 2's
    deallocate (d2)
  end do
  call get_from_xgrid(d1, grid1%id, x_back, xmap)  ! get onto side 1
  conservation_check_side1(3) = mpp_global_sum(grid1%domain, grid1%area*d1)
end function conservation_check_side1
! </FUNCTION>

!#######################################################################
!
! conservation_check - returns three numbers which are the global sum of a
! variable (1) on its home model grid, (2) after interpolation to the other
! side grid(s), and (3) after re_interpolation back onto its home side grid(s).
!
! <FUNCTION NAME="conservation_check_side2" INTERFACE="conservation_check">
!   <IN NAME="d"  TYPE="real" DIM="(:,:,:)" > </IN>
!   <IN NAME="grid_id"  TYPE="character(len=3)"  > </IN>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>
!   <OUT NAME="conservation_check_side2" TYPE="real" DIM="dimension(3)" > </OUT>

function conservation_check_side2(d, grid_id, xmap,remap_method) ! this one for 2->1->2
real, dimension(:,:,:), intent(in   )  :: d
character(len=3),       intent(in   )  :: grid_id
type (xmap_type),       intent(inout)  :: xmap
real, dimension(3)                     :: conservation_check_side2
integer, intent(in), optional :: remap_method


  real, dimension(xmap%size) :: x_over, x_back
  real, dimension(:,:  ), allocatable :: d1
  real, dimension(:,:,:), allocatable :: d2
  integer                             :: g
  type (grid_type), pointer, save     :: grid1 =>NULL(), grid2 =>NULL()

  grid1 => xmap%grids(1)
  do g = 2,size(xmap%grids(:))
    grid2 => xmap%grids(g)
    if (grid_id==grid2%id) then
      conservation_check_side2(1) = mpp_global_sum( grid2%domain, &
                                     grid2%area * sum(grid2%frac_area*d,DIM=3) )
      call put_to_xgrid(d, grid_id, x_over, xmap)  ! put from this side 2
    else
      call put_to_xgrid(0 * grid2%frac_area, grid2%id, x_over, xmap) ! zero rest
    end if
  end do

  allocate ( d1(size(grid1%area,1),size(grid1%area,2)) )
  call get_from_xgrid(d1, grid1%id, x_over, xmap)  ! get onto side 1
  conservation_check_side2(2) = mpp_global_sum(grid1%domain, grid1%area*d1)
  call put_to_xgrid(d1,  grid1%id, x_back, xmap,remap_method)   ! put from side 1
  deallocate ( d1 )

  conservation_check_side2(3) = 0.0;
  do g = 2,size(xmap%grids(:))
    grid2 => xmap%grids(g)
    allocate ( d2 ( size(grid2%frac_area, 1), size(grid2%frac_area, 2),  &
                                              size(grid2%frac_area, 3) ) )
    call get_from_xgrid(d2,  grid2%id, x_back, xmap) ! get onto side 2's
    conservation_check_side2(3) = conservation_check_side2(3)                  &
                                 +mpp_global_sum( grid2%domain,                &
                                    grid2%area * sum(grid2%frac_area*d2,DIM=3) )
    deallocate ( d2 )
  end do
  
end function conservation_check_side2
! </FUNCTION>

!#######################################################################

! This function is used to calculate the gradient along zonal direction.
! Maybe need to setup a limit for the gradient. The grid is assumeed 
! to be regular lat-lon grid

function grad_zonal_1d(d, lon, lat, is, ie, js, je, isd, jsd) 

  integer,                    intent(in) :: isd, jsd
  real, dimension(isd:,jsd:), intent(in) :: d
  real, dimension(:),         intent(in) :: lon
  real, dimension(:),         intent(in) :: lat
  integer,                    intent(in) :: is, ie, js, je 
  real, dimension(is:ie,js:je)           :: grad_zonal_1d
  real                                   :: dx, costheta
  integer                                :: i, j, ip1, im1

  !  calculate the gradient of the data on each grid
  do i = is, ie
     if(i == 1) then
        ip1 = i+1; im1 = i
     else if(i==size(lon(:)) ) then
        ip1 = i; im1 = i-1
     else
        ip1 = i+1; im1 = i-1
     endif
     dx = lon(ip1) - lon(im1)
     if(abs(dx).lt.EPS )  call error_mesg('xgrids_mod(grad_zonal_1d)', 'Improper grid size in lontitude', FATAL)
     if(dx .gt. PI)  dx = dx - 2.0* PI
     if(dx .lt. -PI) dx = dx + 2.0* PI
     do j = js, je
        costheta = cos(lat(j))
        if(abs(costheta) .lt. EPS) call error_mesg('xgrids_mod(grad_zonal_1d)', 'Improper latitude grid', FATAL)
        grad_zonal_1d(i,j) = (d(ip1,j)-d(im1,j))/(dx*costheta)
     enddo
  enddo

  return

end function grad_zonal_1d

!#######################################################################

! This function is used to calculate the gradient along zonal direction.
! Maybe need to setup a limit for the gradient. The grid can be any
! logically rectangular grid.

function grad_zonal_2d(d, lon, lat, is, ie, js, je, isd, jsd) 

  integer,                    intent(in) :: isd, jsd
  real, dimension(isd:,jsd:), intent(in) :: d
  real, dimension(isd:,jsd:), intent(in) :: lon
  real, dimension(isd:,jsd:), intent(in) :: lat
  integer,                    intent(in) :: is, ie, js, je 
  real, dimension(is:ie,js:je)           :: grad_zonal_2d
  real                                   :: dx, costheta
  integer                                :: i, j, ip1, im1

  !  calculate the gradient of the data on each grid
  do j = js, je
     do i = is, ie
        ip1 = i+1; im1 = i-1
        dx = lon(ip1,j) - lon(im1,j)
        if(abs(dx).lt.EPS )  call error_mesg('xgrids_mod(grad_zonal_2d)', 'Improper grid size in lontitude', FATAL)
        if(dx .gt. PI)  dx = dx - 2.0* PI
        if(dx .lt. -PI) dx = dx + 2.0* PI
        costheta = cos(lat(i,j))
        if(abs(costheta) .lt. EPS) call error_mesg('xgrids_mod(grad_zonal_2d)', 'Improper latitude grid', FATAL)
        grad_zonal_2d(i,j) = (d(ip1,j)-d(im1,j))/(dx*costheta)
     enddo
  enddo

  return

end function grad_zonal_2d


!#######################################################################

! This function is used to calculate the gradient along meridinal direction.
! Maybe need to setup a limit for the gradient. regular lat-lon grid are assumed

function grad_merid_1d(d, lat, is, ie, js, je, isd, jsd) 
  integer,                    intent(in) :: isd, jsd
  real, dimension(isd:,jsd:), intent(in) :: d
  real, dimension(:),         intent(in) :: lat
  integer,                    intent(in) :: is, ie, js, je 
  real, dimension(is:ie,js:je)           :: grad_merid_1d
  real                                   :: dy
  integer                                :: i, j, jp1, jm1

  !  calculate the gradient of the data on each grid
  do j = js, je
     if(j == 1) then
        jp1 = j+1; jm1 = j
     else if(j == size(lat(:)) ) then
        jp1 = j;   jm1 = j-1
     else
        jp1 = j+1; jm1 = j-1
     endif
     dy = lat(jp1) - lat(jm1)
     if(abs(dy).lt.EPS) call error_mesg('xgrids_mod(grad_merid_1d)', 'Improper grid size in latitude', FATAL)

     do i = is, ie
        grad_merid_1d(i,j) = (d(i,jp1) - d(i,jm1))/dy
     enddo
  enddo

  return
end function grad_merid_1d

!#######################################################################

! This function is used to calculate the gradient along meridinal direction.
! Maybe need to setup a limit for the gradient. Grid can be any logically
! rectangular grid.

function grad_merid_2d(d, lat, is, ie, js, je, isd, jsd) 
  integer,                    intent(in) :: isd, jsd
  real, dimension(isd:,jsd:), intent(in) :: d
  real, dimension(isd:,jsd:), intent(in) :: lat
  integer,                    intent(in) :: is, ie, js, je 
  real, dimension(is:ie,js:je)           :: grad_merid_2d
  real                                   :: dy
  integer                                :: i, j, jp1, jm1

  !  calculate the gradient of the data on each grid
  do j = js, je
     jp1 = j+1; jm1 = j-1
     do i = is, ie
        dy = lat(i, jp1) - lat(i, jm1)
        if(abs(dy).lt.EPS) call error_mesg('xgrids_mod(grad_merid_2d)', 'Improper grid size in latitude', FATAL)
        grad_merid_2d(i,j) = (d(i,jp1) - d(i,jm1))/dy
     enddo
  enddo

  return
end function grad_merid_2d

!#######################################################################
subroutine get_index_range(xmap, grid_index, is, ie, js, je, km)

  type(xmap_type), intent(in)     :: xmap
  integer, intent(in)             :: grid_index
  integer, intent(out)            :: is, ie, js, je, km

  is = xmap % grids(grid_index) % is_me
  ie = xmap % grids(grid_index) % ie_me
  js = xmap % grids(grid_index) % js_me
  je = xmap % grids(grid_index) % je_me
  km = xmap % grids(grid_index) % km
  
end subroutine get_index_range
!#######################################################################

subroutine stock_move_3d(from, to, grid_index, data, xmap, &
     & delta_t, from_side, to_side, radius, verbose, ier)

  ! this version takes rank 3 data, it can be used to compute the flux on anything but the 
  ! first grid, which typically is on the atmos side.
  ! note that "from" and "to" are optional, the stocks will be subtracted, resp. added, only
  ! if these are present.

  use mpp_mod, only : mpp_sum
  use mpp_domains_mod, only : domain2D, mpp_redistribute, mpp_get_compute_domain

  type(stock_type), intent(inout), optional :: from, to
  integer, intent(in)             :: grid_index        ! grid index
  real, intent(in)                :: data(:,:,:)  ! data array is 3d
  type(xmap_type), intent(in)     :: xmap
  real, intent(in)                :: delta_t
  integer, intent(in)             :: from_side, to_side ! ISTOCK_TOP, ISTOCK_BOTTOM, or ISTOCK_SIDE
  real, intent(in)                :: radius       ! earth radius
  character(len=*), intent(in), optional      :: verbose
  integer, intent(out)            :: ier

  real    :: from_dq, to_dq

  ier = 0
  if(grid_index == 1) then
     ! data has rank 3 so grid index must be > 1
     ier = 1
     return
  endif

  if(.not. associated(xmap%grids) ) then
     ier = 2
     return
  endif

     from_dq = delta_t * 4*PI*radius**2 * sum( sum(xmap%grids(grid_index)%area * &
          & sum(xmap%grids(grid_index)%frac_area * data, DIM=3), DIM=1))
     to_dq = from_dq

  ! update only if argument is present.
  if(present(to  )) to   % dq(  to_side) = to   % dq(  to_side) + to_dq
  if(present(from)) from % dq(from_side) = from % dq(from_side) - from_dq

  if(present(verbose).and.debug_stocks) then
     call mpp_sum(from_dq)
     call mpp_sum(to_dq)
     from_dq = from_dq/(4*PI*radius**2)
     to_dq   = to_dq  /(4*PI*radius**2)
     if(mpp_pe()==mpp_root_pe()) then
        write(stocks_file,'(a,es19.12,a,es19.12,a)') verbose, from_dq,' [*/m^2]'
     endif
  endif

end subroutine stock_move_3d

!...................................................................

subroutine stock_move_2d(from, to, grid_index, data, xmap, &
     & delta_t, from_side, to_side, radius, verbose, ier)

  ! this version takes rank 2 data, it can be used to compute the flux on the atmos side
  ! note that "from" and "to" are optional, the stocks will be subtracted, resp. added, only
  ! if these are present.

  use mpp_mod, only : mpp_sum
  use mpp_domains_mod, only : domain2D, mpp_redistribute, mpp_get_compute_domain

  type(stock_type), intent(inout), optional :: from, to
  integer, optional, intent(in)   :: grid_index
  real, intent(in)                :: data(:,:)    ! data array is 2d
  type(xmap_type), intent(in)     :: xmap
  real, intent(in)                :: delta_t
  integer, intent(in)             :: from_side, to_side ! ISTOCK_TOP, ISTOCK_BOTTOM, or ISTOCK_SIDE
  real, intent(in)                :: radius       ! earth radius
  character(len=*), intent(in)    :: verbose
  integer, intent(out)            :: ier

  real    :: to_dq, from_dq

  ier = 0

  if(.not. associated(xmap%grids) ) then
     ier = 3
     return
  endif

  if( .not. present(grid_index) .or. grid_index==1 ) then

     ! only makes sense if grid_index == 1
     from_dq = delta_t * 4*PI*radius**2 * sum(sum(xmap%grids(1)%area * data, DIM=1))
     to_dq = from_dq

  else

     ier = 4
     return

  endif

  ! update only if argument is present.
  if(present(to  )) to   % dq(  to_side) = to   % dq(  to_side) + to_dq
  if(present(from)) from % dq(from_side) = from % dq(from_side) - from_dq

  if(debug_stocks) then
     call mpp_sum(from_dq)
     call mpp_sum(to_dq)
     from_dq = from_dq/(4*PI*radius**2)
     to_dq   = to_dq  /(4*PI*radius**2)
     if(mpp_pe()==mpp_root_pe()) then
        write(stocks_file,'(a,es19.12,a,es19.12,a)') verbose, from_dq,' [*/m^2]'
     endif
  endif

end subroutine stock_move_2d

!#######################################################################
subroutine stock_integrate_2d(data, xmap, delta_t, radius, res, ier)

  ! surface/time integral of a 2d array

  use mpp_mod, only : mpp_sum

  real, intent(in)                :: data(:,:)    ! data array is 2d
  type(xmap_type), intent(in)     :: xmap
  real, intent(in)                :: delta_t
  real, intent(in)                :: radius       ! earth radius
  real, intent(out)               :: res
  integer, intent(out)            :: ier

  real :: tmp

  ier = 0
  res = 0

  if(.not. associated(xmap%grids) ) then
     ier = 6
     return
  endif

  res = delta_t * 4*PI*radius**2 * sum(sum(xmap%grids(1)%area * data, DIM=1))

end subroutine stock_integrate_2d
!#######################################################################

!#######################################################################



subroutine stock_print(stck, Time, comp_name, index, ref_value, radius, pelist)

  use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_sum
  use time_manager_mod, only : time_type, get_time
  use diag_manager_mod, only : register_diag_field,send_data

  type(stock_type), intent(in)  :: stck
  type(time_type), intent(in)   :: Time
  character(len=*)              :: comp_name
  integer, intent(in)           :: index     ! to map stock element (water, heat, ..) to a name
  real, intent(in)              :: ref_value ! the stock value returned by the component per PE
  real, intent(in)              :: radius
  integer, intent(in), optional :: pelist(:)

  real :: f_value, c_value, planet_area
  character(len=80) :: formatString
  integer :: iday, isec, hours
  integer :: diagID, compInd
  integer, dimension(NELEMS,4), save :: f_valueDiagID = -1
  integer, dimension(NELEMS,4), save :: c_valueDiagID = -1
  integer, dimension(NELEMS,4), save :: fmc_valueDiagID = -1

  real :: diagField
  logical :: used
  character(len=30) :: field_name, units

  f_value = sum(stck % dq)
  c_value = ref_value      - stck % q_start
  if(present(pelist)) then
     call mpp_sum(f_value, pelist=pelist)
     call mpp_sum(c_value, pelist=pelist)
  else
     call mpp_sum(f_value)
     call mpp_sum(c_value)
  endif

  if(mpp_pe() == mpp_root_pe()) then
     ! normalize to 1 earth m^2
     planet_area = 4*PI*radius**2
     f_value       = f_value     / planet_area
     c_value       = c_value     / planet_area

     if(comp_name == 'ATM') compInd = 1
     if(comp_name == 'LND') compInd = 2
     if(comp_name == 'ICE') compInd = 3
     if(comp_name == 'OCN') compInd = 4


     if(f_valueDiagID(index,compInd) == -1) then
        field_name = trim(comp_name) // trim(STOCK_NAMES(index))
        field_name  = trim(field_name) // 'StocksChange_Flux'
        units = trim(STOCK_UNITS(index))
        f_valueDiagID(index,compInd) = register_diag_field('stock_print', field_name, Time, &
             units=units)
     endif

     if(c_valueDiagID(index,compInd) == -1) then
        field_name = trim(comp_name) // trim(STOCK_NAMES(index))
        field_name = trim(field_name) // 'StocksChange_Comp'
        units = trim(STOCK_UNITS(index))
        c_valueDiagID(index,compInd) = register_diag_field('stock_print', field_name, Time, &
             units=units)
     endif

     if(fmc_valueDiagID(index,compInd) == -1) then
        field_name = trim(comp_name) // trim(STOCK_NAMES(index))
        field_name = trim(field_name) // 'StocksChange_Diff'
        units = trim(STOCK_UNITS(index))
        fmc_valueDiagID(index,compInd) = register_diag_field('stock_print', field_name, Time, &
             units=units)
     endif

     DiagID=f_valueDiagID(index,compInd)
     diagField = f_value
     if (DiagID > 0)  used = send_data(DiagID, diagField, Time)
     DiagID=c_valueDiagID(index,compInd)
     diagField = c_value
     if (DiagID > 0)  used = send_data(DiagID, diagField, Time)
     DiagID=fmc_valueDiagID(index,compInd)
     diagField = f_value-c_value
     if (DiagID > 0)  used = send_data(DiagID, diagField, Time)


     call get_time(Time, isec, iday)
     hours = iday*24 + isec/3600
     formatString = '(a,a,a,i16,2x,es22.15,2x,es22.15,2x,es22.15)'
     write(stocks_file,formatString) trim(comp_name),STOCK_NAMES(index),STOCK_UNITS(index) &
          ,hours,f_value,c_value,f_value-c_value

  endif


end subroutine stock_print


!###############################################################################
function is_lat_lon(lon, lat)
  real, dimension(:,:), intent(in) :: lon, lat
  logical                          :: is_lat_lon
  integer                          :: i, j, nlon, nlat

  is_lat_lon = .true.
  nlon = size(lon,1)
  nlat = size(lon,2)
  do j = 1, nlat
     do i = 2, nlon
        if(lat(i,j) .NE. lat(1,j)) then
           is_lat_lon = .false.
           return
        end if
     end do
  end do

  do i = 1, nlon
     do j = 2, nlat
        if(lon(i,j) .NE. lon(i,1)) then
           is_lat_lon = .false.
           return
        end if
     end do
  end do

  return
end function is_lat_lon

end module xgrid_mod


! <INFO>

!   <REFERENCE>   
!      A <LINK SRC="http://www.gfdl.noaa.gov/~mw/docs/grid_coupling.html"> guide </LINK>to grid coupling in FMS.
!   </REFERENCE>
!   <REFERENCE>
!      A simple xgrid <LINK SRC="http://www.gfdl.gov/~mw/docs/xgrid_example.f90.txt"> example. </LINK>
!   </REFERENCE>

! </INFO>


!======================================================================================
! standalone unit test

#ifdef TEST_XGRID
! Now only test some simple test, will test cubic grid mosaic in the future.

program xgrid_test

  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_error, FATAL
  use mpp_domains_mod, only : mpp_define_domains, mpp_define_layout, mpp_domains_exit
  use mpp_domains_mod, only : mpp_get_compute_domain, domain2d, mpp_domains_init
  use mpp_domains_mod, only : mpp_define_mosaic_pelist, mpp_define_mosaic
  use fms_mod,         only : fms_init, file_exist, field_size, open_namelist_file
  use fms_mod,         only : check_nml_error, close_file, read_data, stdout, fms_end
  use xgrid_mod,       only : xgrid_init, setup_xmap, put_to_xgrid, get_from_xgrid
  use xgrid_mod,       only : xmap_type, xgrid_count
  use mosaic_mod,      only : get_mosaic_ntiles, get_mosaic_grid_sizes


implicit none

  real, parameter :: EPSLN = 1.0e-10
  integer :: grid_version = 2         ! = 1, read 'INPUT/grid_spec.nc'
                                      ! = 2, read 'INPUT/mosaic.nc'
  namelist /xgrid_test_nml/ grid_version

  integer              :: remap_method
  integer              :: pe, npes, ierr, nml_unit, io, n
  integer              :: siz(4), ntile_lnd, ntile_atm, ntile_ocn, ncontact
  integer, allocatable :: layout(:,:), global_indices(:,:)
  integer, allocatable :: atm_nx(:), atm_ny(:), ocn_nx(:), ocn_ny(:), lnd_nx(:), lnd_ny(:)
  integer, allocatable :: pe_start(:), pe_end(:)
  integer, allocatable :: dummy(:)
  character(len=256)   :: grid_file = "INPUT/grid_spec.nc"
  character(len=256)   :: atm_mosaic, ocn_mosaic, lnd_mosaic
  character(len=256)   :: atm_mosaic_file, ocn_mosaic_file, lnd_mosaic_file
  type(domain2d)       :: Atm_domain, Ocn_domain, Lnd_domain
  type(xmap_type)      :: Xmap


  call fms_init
  call mpp_domains_init
  call xgrid_init(remap_method)

  npes   = mpp_npes()
  pe     = mpp_pe()

 if (file_exist('input.nml')) then
   ierr=1
   nml_unit = open_namelist_file()
   do while (ierr /= 0)
     read(nml_unit, nml=xgrid_test_nml, iostat=io, end=10)
          ierr = check_nml_error(io, 'xgrid_test_nml')
   enddo
10 call close_file(nml_unit)
 endif

  select case(grid_version)
  case ( 1 )
     allocate(atm_nx(1), atm_ny(1))
     allocate(lnd_nx(1), lnd_ny(1))
     allocate(ocn_nx(1), ocn_ny(1))
     allocate(layout(1,2))
     call field_size(grid_file, "AREA_ATM", siz )
     atm_nx = siz(1); atm_ny = siz(2)
     call field_size(grid_file, "AREA_OCN", siz )
     ocn_nx = siz(1); ocn_ny = siz(2)
     call field_size(grid_file, "AREA_LND", siz )
     lnd_nx = siz(1); lnd_ny = siz(2)
     call mpp_define_layout( (/1,atm_nx,1,atm_ny/), npes, layout(1,:)) 
     call mpp_define_domains( (/1,atm_nx,1,atm_ny/), layout(1,:), Atm_domain)
     call mpp_define_layout( (/1,lnd_nx,1,lnd_ny/), npes, layout(1,:)) 
     call mpp_define_domains( (/1,lnd_nx,1,lnd_ny/), layout(1,:), Lnd_domain) 
     call mpp_define_layout( (/1,ocn_nx,1,ocn_ny/), npes, layout(1,:))
     call mpp_define_domains( (/1,ocn_nx,1,ocn_ny/), layout(1,:), Ocn_domain)
     deallocate(layout)
  case ( 2 )
     !--- Get the mosaic data of each component model 
     call read_data(grid_file, 'atm_mosaic', atm_mosaic)
     call read_data(grid_file, 'lnd_mosaic', lnd_mosaic)
     call read_data(grid_file, 'ocn_mosaic', ocn_mosaic)
     atm_mosaic_file = 'INPUT/'//trim(atm_mosaic)//'.nc'
     lnd_mosaic_file = 'INPUT/'//trim(lnd_mosaic)//'.nc'
     ocn_mosaic_file = 'INPUT/'//trim(ocn_mosaic)//'.nc'

     ntile_lnd = get_mosaic_ntiles(lnd_mosaic_file);
     ntile_ocn = get_mosaic_ntiles(ocn_mosaic_file);
     ntile_atm = get_mosaic_ntiles(atm_mosaic_file);
     if(ntile_lnd > 1) call mpp_error(FATAL,  &
         'xgrid_test: there is more than one tile in lnd_mosaic, which is not implemented yet')
     if(ntile_ocn > 1) call mpp_error(FATAL,  &
           'xgrid_test: there is more than one tile in ocn_mosaic, which is not implemented yet')

     write(stdout(),*)" There is ", ntile_atm, " tiles in atmos mosaic"
     write(stdout(),*)" There is ", ntile_lnd, " tiles in land  mosaic"
     write(stdout(),*)" There is ", ntile_ocn, " tiles in ocean mosaic"
     allocate(atm_nx(ntile_atm), atm_ny(ntile_atm))
     allocate(lnd_nx(ntile_ocn), lnd_ny(ntile_lnd))
     allocate(ocn_nx(ntile_ocn), ocn_ny(ntile_ocn))

     call get_mosaic_grid_sizes(atm_mosaic_file, atm_nx, atm_ny)
     call get_mosaic_grid_sizes(lnd_mosaic_file, lnd_nx, lnd_ny)
     call get_mosaic_grid_sizes(ocn_mosaic_file, ocn_nx, ocn_ny)

     ! since no update is needed, no need to get the contacts 
     ncontact = 0

     if(mod(npes, ntile_atm) .NE. 0 ) call mpp_error(FATAL,"npes should be divided by ntile_atm")

     allocate(pe_start(ntile_atm), pe_end(ntile_atm) )
     allocate(global_indices(4, ntile_atm), layout(2,ntile_atm))
     call mpp_define_mosaic_pelist( atm_nx*atm_ny, pe_start, pe_end)
     do n = 1, ntile_atm
        global_indices(:,n) = (/1, atm_nx(n), 1, atm_ny(n)/)
        call mpp_define_layout( global_indices(:,n), pe_end(n)-pe_start(n)+1, layout(:,n))
     end do
 
     call mpp_define_mosaic(global_indices, layout, Atm_domain, ntile_atm, ncontact, dummy, dummy, &
                            dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, pe_start, pe_end)
     deallocate( pe_start, pe_end, global_indices, layout )

     allocate(pe_start(ntile_lnd), pe_end(ntile_lnd) )
     allocate(global_indices(4,ntile_lnd), layout(2,ntile_lnd))
     call mpp_define_mosaic_pelist( lnd_nx*lnd_ny, pe_start, pe_end)
     do n = 1, ntile_lnd
        global_indices(:,n) = (/1, lnd_nx(n), 1, lnd_ny(n)/)
        call mpp_define_layout( global_indices(:,n), pe_end(n)-pe_start(n)+1, layout(:,n))
     end do
 
     call mpp_define_mosaic(global_indices, layout, Lnd_domain, ntile_lnd, ncontact, dummy, dummy, &
                            dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, pe_start, pe_end)
     deallocate( pe_start, pe_end, global_indices, layout )

     allocate(pe_start(ntile_ocn), pe_end(ntile_ocn) )
     allocate(global_indices(4, ntile_ocn), layout(2, ntile_ocn))
     call mpp_define_mosaic_pelist( ocn_nx*ocn_ny, pe_start, pe_end)
     do n = 1, ntile_ocn
        global_indices(:,n) = (/1, ocn_nx(n), 1, ocn_ny(n)/)
        call mpp_define_layout( global_indices(:,n), pe_end(n)-pe_start(n)+1, layout(:,n))
     end do
 
     call mpp_define_mosaic(global_indices, layout, Ocn_domain, ntile_ocn, ncontact, dummy, dummy, &
                            dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, pe_start, pe_end)
     deallocate( pe_start, pe_end, global_indices, layout )

  case default 
     call mpp_error(FATAL, "xgrid_test: nml grid_version should be 1 or 2")
  end select

  deallocate(atm_nx, atm_ny, lnd_nx, lnd_ny, ocn_nx, ocn_ny)

  !--- conservation check is done in setup_xmap. 
  call setup_xmap(Xmap, (/ 'ATM', 'OCN', 'LND' /), (/ Atm_domain, Ocn_domain, Lnd_domain /), grid_file)

  write(stdout(),*) "************************************************************************"
  write(stdout(),*) "***********      Finish running program test_xgrid         *************"
  write(stdout(),*) "************************************************************************"

  call mpp_domains_exit
  call fms_end

end program xgrid_test

! end of TEST_XGRID
#endif


#ifdef _XGRID_MAIN
! to compile on Altix:
! setenv FMS /net2/ap/regression/ia64/10-Aug-2006/CM2.1U_Control-1990_E1.k_dyn30pe/exec
! ifort -fpp -r8 -i4 -g -check all -D_XGRID_MAIN -I $FMS xgrid.f90 $FMS/stock_constants.o $FMS/fms*.o $FMS/mpp*.o $FMS/constants.o $FMS/time_manager.o $FMS/memutils.o $FMS/threadloc.o -L/usr/local/lib -lnetcdf -L/usr/lib -lmpi -lsma
! mpirun -np 30 a.out
program main
  use mpp_mod
  use fms_mod
  use mpp_domains_mod
  use xgrid_mod, only : xmap_type, setup_xmap, stock_move, stock_type, get_index_range
  use stock_constants_mod, only : ISTOCK_TOP, ISTOCK_BOTTOM, ISTOCK_SIDE, ISTOCK_WATER, ISTOCK_HEAT, NELEMS
  use constants_mod, only       : PI
  implicit none

  type(xmap_type)  :: xmap_sfc, xmap_runoff
  integer          :: npes, pe, root, i, nx, ny, ier
  integer          :: patm_beg, patm_end, pocn_beg, pocn_end
  integer          :: is, ie, js, je, km, index_ice, index_lnd
  integer          :: layout(2)
  type(stock_type), save :: Atm_stock(NELEMS), Ice_stock(NELEMS), &
       &                    Lnd_stock(NELEMS), Ocn_stock(NELEMS)
  type(domain2D)   :: Atm_domain, Ice_domain, Lnd_domain, Ocn_domain
  logical, pointer :: maskmap(:,:)
  real, allocatable :: data2d(:,:), data3d(:,:,:)
  real              :: dt, dq_tot_atm, dq_tot_ice, dq_tot_lnd, dq_tot_ocn



  call fms_init

  npes   = mpp_npes()
  pe     = mpp_pe()
  root   = mpp_root_pe()
  patm_beg = 0
  patm_end = npes/2 - 1
  pocn_beg = patm_end + 1
  pocn_end = npes - 1

  if(npes /= 30) call mpp_error(FATAL,'must run unit test on 30 pes')

  call mpp_domains_init ! (MPP_DEBUG)

  call mpp_declare_pelist( (/ (i, i=patm_beg, patm_end) /), 'atm_lnd_ice pes' ) 
  call mpp_declare_pelist( (/ (i, i=pocn_beg, pocn_end) /), 'ocn pes' ) 

  index_ice = 2 ! 2nd exchange grid
  index_lnd = 3 ! 3rd exchange grid

  dt = 1.0

  if(pe < 15) then

     call mpp_set_current_pelist( (/ (i, i=patm_beg, patm_end) /) )

     ! Lnd
     nx = 144
     ny = 90
     layout = (/ 5, 3 /)
     call mpp_define_domains( (/1,nx, 1,ny/), layout, Lnd_domain, &
          & xflags = CYCLIC_GLOBAL_DOMAIN, name = 'LAND MODEL' )

     ! Atm
     nx = 144
     ny = 90
     layout = (/1, 15/)
     call mpp_define_domains( (/1,nx, 1,ny/), layout, Atm_domain)

     ! Ice
     nx = 360
     ny = 200
     layout = (/15, 1/)
     call mpp_define_domains( (/1,nx, 1,ny/), layout, Ice_domain, name='ice_nohalo' )

     ! Build exchange grid
     call setup_xmap(xmap_sfc, (/ 'ATM', 'OCN', 'LND' /), &
          &                    (/ Atm_domain, Ice_domain, Lnd_domain /), &
          &                    "INPUT/grid_spec.nc")

!  call setup_xmap(xmap_sfc, (/ 'LND', 'OCN' /), &
!       &                    (/  Lnd_domain, Ice_domain /), &
!       &                    "INPUT/grid_spec.nc")


     ! Atm -> Ice

     i = index_ice

     call get_index_range(xmap=xmap_sfc, grid_index=i, is=is, ie=ie, js=js, je=je, km=km)

     allocate(data3d(is:ie, js:je, km))
     data3d(:,:,1   ) = 1.0/(4.0*PI)
     data3d(:,:,2:km) = 0.0
     call stock_move(from=Atm_stock(ISTOCK_WATER), to=Ice_stock(ISTOCK_WATER), &
          & grid_index=i, data=data3d, xmap=xmap_sfc, &
          & delta_t=dt, from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, radius=1.0, ier=ier)
     deallocate(data3d)

     ! Atm -> Lnd

     i = index_lnd

     call get_index_range(xmap=xmap_sfc, grid_index=i, is=is, ie=ie, js=js, je=je, km=km)
     
     allocate(data3d(is:ie, js:je, km))
     data3d(:,:,1   ) = 1.0/(4.0*PI)
     data3d(:,:,2:km) = 0.0
     call stock_move(from=Atm_stock(ISTOCK_WATER), to=Lnd_stock(ISTOCK_WATER), &
          & grid_index=i, data=data3d, xmap=xmap_sfc, &
          & delta_t=dt, from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, radius=1.0, ier=ier)
     deallocate(data3d)

  else ! pes: 15...29

     call mpp_set_current_pelist( (/ (i, i=pocn_beg, pocn_end) /) )

     ! Ocn
     nx = 360
     ny = 200
     layout = (/ 5, 3 /)
     call mpp_define_domains( (/1,nx,1,ny/), layout, Ocn_domain, name='ocean model')

  endif

  ! Ice -> Ocn (same grid different layout)

  i = index_ice

  if( pe < pocn_beg ) then 

     call get_index_range(xmap=xmap_sfc, grid_index=i, is=is, ie=ie, js=js, je=je, km=km)

     allocate(data3d(is:ie, js:je, km))
     data3d(:,:,1   ) = 1.0/(4.0*PI)
     data3d(:,:,2:km) = 0.0
  else
     is = 0
     ie = 0
     js = 0
     je = 0
     km = 0
     allocate(data3d(is:ie, js:je, km))
  endif

  call stock_move(from=Ice_stock(ISTOCK_WATER), to=Ocn_stock(ISTOCK_WATER), &
       & grid_index=i, data=data3d(:,:,1), xmap=xmap_sfc, &
       & delta_t=dt, from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, radius=1.0, ier=ier)
  deallocate(data3d)

  ! Sum across sides and PEs

  dq_tot_atm = sum(Atm_stock(ISTOCK_WATER)%dq)
  call mpp_sum(dq_tot_atm)

  dq_tot_lnd = sum(Lnd_stock(ISTOCK_WATER)%dq)
  call mpp_sum(dq_tot_lnd)

  dq_tot_ice = sum(Ice_stock(ISTOCK_WATER)%dq)
  call mpp_sum(dq_tot_ice)

  dq_tot_ocn = sum(Ocn_stock(ISTOCK_WATER)%dq)
  call mpp_sum(dq_tot_ocn)

  if(pe==root) then
     write(*,'(a,4f10.7,a,e10.2)') ' Total delta_q(water) Atm/Lnd/Ice/Ocn: ', &
          & dq_tot_atm, dq_tot_lnd, dq_tot_ice, dq_tot_ocn, &
          & ' residue: ', dq_tot_atm + dq_tot_lnd + dq_tot_ice + dq_tot_ocn
  endif

 
end program main
! end of _XGRID_MAIN
#endif
