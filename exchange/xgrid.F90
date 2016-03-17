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

#include <fms_platform.h>

use       fms_mod,   only: file_exist, open_namelist_file, check_nml_error,  &
                           error_mesg, close_file, FATAL, NOTE, stdlog,      &
                           write_version_number, read_data, field_exist,     &
                           field_size, lowercase, string,                    &
                           get_mosaic_tile_grid
use     fms_io_mod,  only: get_var_att_value
use mpp_mod,         only: mpp_npes, mpp_pe, mpp_root_pe, mpp_send, mpp_recv, &
                           mpp_sync_self, stdout, mpp_max, EVENT_RECV,        &
                           mpp_get_current_pelist, mpp_clock_id, mpp_min,     &
                           mpp_alltoall,                                      &
                           mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC,    &
                           COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4,    &
                           COMM_TAG_5, COMM_TAG_6, COMM_TAG_7, COMM_TAG_8,    &
                           COMM_TAG_9, COMM_TAG_10
use mpp_mod,         only: input_nml_file, mpp_set_current_pelist, mpp_sum, mpp_sync
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                           Domain2d, mpp_global_sum, mpp_update_domains,    &
                           mpp_modify_domain, mpp_get_data_domain, XUPDATE, &
                           YUPDATE, mpp_get_current_ntile, mpp_get_tile_id, &
                           mpp_get_ntile_count, mpp_get_tile_list,          &
                           mpp_get_global_domain, Domain1d,                 &
                           mpp_deallocate_domain, mpp_define_domains,       &
                           mpp_get_domain_npes, mpp_get_domain_root_pe,     &
                           mpp_domain_is_initialized, mpp_broadcast_domain, &
                           mpp_get_domain_pelist, mpp_compute_extent
use mpp_io_mod,      only: mpp_open, MPP_MULTI, MPP_SINGLE, MPP_OVERWR
use constants_mod,   only: PI, RADIUS
use mosaic_mod,          only: get_mosaic_xgrid, get_mosaic_xgrid_size, &
                               get_mosaic_ntiles, get_mosaic_ncontacts, &
                               get_mosaic_contact, get_mosaic_grid_sizes

use stock_constants_mod, only: ISTOCK_TOP, ISTOCK_BOTTOM, ISTOCK_SIDE, STOCK_NAMES, &
                               STOCK_UNITS, NELEMS, stocks_file, stock_type
use gradient_mod,        only: gradient_cubic

implicit none
private

public xmap_type, setup_xmap, set_frac_area, put_to_xgrid, get_from_xgrid, &
       xgrid_count, some, conservation_check, xgrid_init, &
       AREA_ATM_SPHERE, AREA_LND_SPHERE, AREA_OCN_SPHERE, &
       AREA_ATM_MODEL, AREA_LND_MODEL, AREA_OCN_MODEL, &
       get_ocean_model_area_elements, grid_box_type,   &
       get_xmap_grid_area

!--- paramters that determine the remapping method
integer, parameter :: FIRST_ORDER        = 1
integer, parameter :: SECOND_ORDER       = 2
integer, parameter :: VERSION1           = 1 ! grid spec file
integer, parameter :: VERSION2           = 2 ! mosaic grid file
integer, parameter :: MAX_FIELDS         = 80

! <NAMELIST NAME="xgrid_nml">
!   <DATA NAME="make_exchange_reproduce" TYPE="logical"  DEFAULT=".false.">
!     Set to .true. to make <TT>xgrid_mod</TT> reproduce answers on different
!     numbers of PEs.  This option has a considerable performance impact.
!   </DATA>
!   <DATA NAME="interp_method" TYPE="character(len=64)"  DEFAULT=" 'first_order' ">
!     exchange grid interpolation method. It has two options: 
!     "first_order", "second_order".
!   </DATA>
!   <DATA NAME="xgrid_log" TYPE="logical"  DEFAULT=" .false. ">
!     Outputs exchange grid information to xgrid.out.<pe> for debug/diag purposes.
!   </DATA>
!   <DATA NAME="nsubset" TYPE="integer" DEFAULT="0">
!     number of processors to read exchange grid information. Those processors that read
!     the exchange grid information will send data to other processors to prepare for flux exchange.
!     Default value is 0. When nsubset is 0, each processor will read part of the exchange grid 
!     information. The purpose of this namelist is to improve performance of setup_xmap when running
!     on highr processor count and solve receiving size mismatch issue on high processor count.
!     Try to set nsubset = mpp_npes/MPI_rank_per_node.
!   </DATA>
logical :: make_exchange_reproduce = .false. ! exactly same on different # PEs
logical :: xgrid_log = .false. 
character(len=64) :: interp_method = 'first_order'
logical :: debug_stocks = .false. 
logical :: xgrid_clocks_on = .false.
logical :: monotonic_exchange = .false.
integer :: nsubset = 0 ! 0 means mpp_npes()
logical :: do_alltoall = .true.
logical :: do_alltoallv = .false.
namelist /xgrid_nml/ make_exchange_reproduce, interp_method, debug_stocks, xgrid_log, xgrid_clocks_on, &
    monotonic_exchange, nsubset, do_alltoall, do_alltoallv
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
!     FIRST_ORDER (=1), SECOND_ORDER(=2). Default value is FIRST_ORDER.
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

type xcell_type
  integer :: i1, j1, i2, j2 ! indices of cell in model arrays on both sides
  integer :: recv_pos       ! position in the receive buffer.
  integer :: pe             ! other side pe that has this cell
  integer :: tile           ! tile index of side 1 mosaic.
  real    :: area           ! geographic area of exchange cell
!  real    :: area1_ratio     !(= x_area/grid1_area), will be added in the future to improve efficiency
!  real    :: area2_ratio     !(= x_area/grid2_area), will be added in the future to improve efficiency
  real    :: di, dj         ! Weight for the gradient of flux
  real    :: scale
end type xcell_type

type grid_box_type
   real, dimension(:,:),   pointer :: dx     => NULL()
   real, dimension(:,:),   pointer :: dy     => NULL()
   real, dimension(:,:),   pointer :: area   => NULL()
   real, dimension(:),     pointer :: edge_w => NULL()
   real, dimension(:),     pointer :: edge_e => NULL()
   real, dimension(:),     pointer :: edge_s => NULL()
   real, dimension(:),     pointer :: edge_n => NULL()
   real, dimension(:,:,:), pointer :: en1    => NULL()
   real, dimension(:,:,:), pointer :: en2    => NULL()
   real, dimension(:,:,:), pointer :: vlon   => NULL()
   real, dimension(:,:,:), pointer :: vlat   => NULL()
end type grid_box_type

type grid_type
  character(len=3)                :: id                               ! grid identifier
  integer                         :: npes                             ! number of processor on this grid.
  logical                         :: on_this_pe                       ! indicate the domain is defined on this pe
  integer                         :: root_pe                          ! indicate the root pe of the domain
  integer, pointer, dimension(:)  :: pelist                           ! pelist of the domain
  integer                         :: ntile                            ! number of tiles in mosaic
  integer                         :: ni, nj                           ! max of global size of all the tiles
  integer, pointer, dimension(:)  :: tile =>NULL()                    ! tile id ( pe index )
  integer, pointer, dimension(:)  :: is =>NULL(), ie =>NULL()         ! domain - i-range (pe index)
  integer, pointer, dimension(:)  :: js =>NULL(), je =>NULL()         ! domain - j-range (pe index)
  integer, pointer                :: is_me =>NULL(),  ie_me =>NULL()  ! my domain - i-range
  integer, pointer                :: js_me =>NULL(),  je_me =>NULL()  ! my domain - j-range
  integer                         :: isd_me, ied_me                   ! my data domain - i-range
  integer                         :: jsd_me, jed_me                   ! my data domain - j-range
  integer                         :: nxd_me, nyd_me                   ! data domain size
  integer                         :: nxc_me, nyc_me                   ! compute domain size
  integer, pointer                :: tile_me                          ! my tile id
  integer                         :: im , jm , km                     ! global domain range
  real, pointer, dimension(:)     :: lon =>NULL(), lat =>NULL()       ! center of global grids
  real, pointer, dimension(:,:)   :: geolon=>NULL(), geolat=>NULL()   ! geographical grid center
  real, pointer, dimension(:,:,:) :: frac_area =>NULL()               ! partition fractions
  real, pointer, dimension(:,:)   :: area =>NULL()                    ! cell area
  real, pointer, dimension(:,:)   :: area_inv =>NULL()                ! 1 / area for normalization
  integer                         :: first, last                      ! xgrid index range
  integer                         :: first_get, last_get              ! xgrid index range for get_2_from_xgrid
  integer                         :: size                             ! # xcell patterns
  type(xcell_type), pointer       :: x(:) =>NULL()                    ! xcell patterns
  integer                         :: size_repro                       ! # side 1 patterns for repro
  type(xcell_type), pointer       :: x_repro(:) =>NULL()              ! side 1 patterns for repro
  type(Domain2d)                  :: domain                           ! used for conservation checks
  type(Domain2d)                  :: domain_with_halo                 ! used for second order remapping
  logical                         :: is_latlon                        ! indicate if the grid is lat-lon grid or not.
  type(grid_box_type)             :: box                              ! used for second order remapping.
end type grid_type

type x1_type
  integer :: i, j
  real    :: area   ! (= geographic area * frac_area)
!  real    :: area_ratio !(= x1_area/grid1_area) ! will be added in the future to improve efficiency
  real    :: di, dj ! weight for the gradient of flux
  integer :: tile           ! tile index of side 1 mosaic.
  integer :: pos
end type x1_type

type x2_type
  integer :: i, j, k, pos
  real    :: area   ! geographic area of exchange cell
!  real    :: area_ratio !(=x2_area/grid2_area )  ! will be added in the future to improve efficiency
end type x2_type

type overlap_type
   integer          :: count
   integer          :: pe
   integer          :: buffer_pos
   integer, _ALLOCATABLE :: i(:) _NULL
   integer, _ALLOCATABLE :: j(:) _NULL
   integer, _ALLOCATABLE :: g(:) _NULL
   integer, _ALLOCATABLE :: xLoc(:) _NULL
   integer, _ALLOCATABLE :: tile(:) _NULL
   real,    _ALLOCATABLE :: di(:) _NULL
   real,    _ALLOCATABLE :: dj(:) _NULL
end type overlap_type

type comm_type
  integer                         :: nsend, nrecv
  integer                         :: sendsize, recvsize
  integer,            pointer, dimension(:) :: unpack_ind=>NULL()
  type(overlap_type), pointer, dimension(:) :: send=>NULL()
  type(overlap_type), pointer, dimension(:) :: recv=>NULL()  
end type comm_type

type xmap_type
  private
  integer :: size            ! # of exchange grid cells with area > 0 on this pe
  integer :: size_put1       ! # of exchange grid cells for put_1_to_xgrid
  integer :: size_get2       ! # of exchange grid cells for get_2_to_xgrid
  integer :: me, npes, root_pe
  logical, pointer, dimension(:) :: your1my2  =>NULL()! true if side 1 domain on
                                                      ! indexed pe overlaps side 2
                                                      ! domain on this pe
  logical, pointer, dimension(:) :: your2my1 =>NULL() ! true if a side 2 domain on
                                                      ! indexed pe overlaps side 1
                                                      ! domain on this pe
  integer, pointer, dimension(:) :: your2my1_size=>NULL() ! number of exchange grid of 
                                                          ! a side 2 domain on
                                                          ! indexed pe overlaps side 1
                                                          ! domain on this pe

  type (grid_type), pointer, dimension(:) :: grids =>NULL() ! 1st grid is side 1;
                                                            ! rest on side 2
  !
  ! Description of the individual exchange grid cells (index is cell #)
  !
  type(x1_type), pointer, dimension(:) :: x1 =>NULL() ! side 1 info
  type(x1_type), pointer, dimension(:) :: x1_put =>NULL() ! side 1 info
  type(x2_type), pointer, dimension(:) :: x2 =>NULL() ! side 2 info
  type(x2_type), pointer, dimension(:) :: x2_get =>NULL() ! side 2 info

  integer, pointer, dimension(:) :: send_count_repro =>NULL()
  integer, pointer, dimension(:) :: recv_count_repro  =>NULL()
  integer                        :: send_count_repro_tot ! sum(send_count_repro)
  integer                        :: recv_count_repro_tot ! sum(recv_count_repro)
  integer :: version                                  ! version of xgrids. version=VERSION! is for grid_spec file 
                                                      ! and version=VERSION2 is for mosaic grid.
  integer, pointer, dimension(:) :: ind_get1 =>NULL() ! indx for side1 get and side2 put.
  integer, pointer, dimension(:) :: ind_put1 =>NULL() ! indx for side1 put and side 2get.
  type(comm_type), pointer       :: put1 =>NULL()      ! for put_1_to_xgrid
  type(comm_type), pointer       :: get1 =>NULL()      ! for get_1_from_xgrid
  type(comm_type), pointer       :: get1_repro =>NULL()! for get_1_from_xgrid_repro
end type xmap_type

!-----------------------------------------------------------------------
! Include variable "version" to be written to log file.
#include<file_version.h>

 real, parameter                              :: EPS = 1.0e-10
 real, parameter                              :: LARGE_NUMBER = 1.e20
 logical :: module_is_initialized = .FALSE.
 integer :: id_put_1_to_xgrid_order_1 = 0
 integer :: id_put_1_to_xgrid_order_2 = 0
 integer :: id_get_1_from_xgrid = 0
 integer :: id_get_1_from_xgrid_repro = 0
 integer :: id_get_2_from_xgrid = 0
 integer :: id_put_2_to_xgrid = 0
 integer :: id_setup_xmap = 0
 integer :: id_load_xgrid1, id_load_xgrid2, id_load_xgrid3
 integer :: id_load_xgrid4, id_load_xgrid5
 integer :: id_load_xgrid, id_set_comm, id_regen, id_conservation_check


 ! The following is for nested model
 integer :: nnest=0, tile_nest, tile_parent
 integer :: is_nest=0, ie_nest=0, js_nest=0, je_nest=0
 integer :: is_parent=0, ie_parent=0, js_parent=0, je_parent=0 

 ! The following is required to compute stocks of water, heat, ...

  interface stock_move
     module procedure stock_move_3d, stock_move_2d
  end interface

  public stock_move, stock_type, stock_print, get_index_range, stock_integrate_2d
  public FIRST_ORDER, SECOND_ORDER

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
!     FIRST_ORDER (=1), SECOND_ORDER(=2).
!   </OUT>
subroutine xgrid_init(remap_method) 
  integer, intent(out) :: remap_method

  integer :: unit, ierr, io, out_unit

  if (module_is_initialized) return
  module_is_initialized = .TRUE.


#ifdef INTERNAL_FILE_NML
      read (input_nml_file, xgrid_nml, iostat=io)
      ierr = check_nml_error ( io, 'xgrid_nml' )
#else
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
#endif

!--------- write version number and namelist ------------------
  call write_version_number("XGRID_MOD", version)

  unit = stdlog ( )
  out_unit = stdout()
  if ( mpp_pe() == mpp_root_pe() ) write (unit,nml=xgrid_nml)
  call close_file (unit)

!--------- check interp_method has suitable value
!--- when monotonic_exchange is true, interp_method must be second order.

  select case(trim(interp_method))
  case('first_order')
     remap_method = FIRST_ORDER
     if( monotonic_exchange ) call error_mesg('xgrid_mod', &
         'xgrid_nml monotonic_exchange must be .false. when interp_method = first_order', FATAL)  
     write(out_unit,*)"NOTE from xgrid_mod: use first_order conservative exchange"
  case('second_order')
     if(monotonic_exchange) then
        write(out_unit,*)"NOTE from xgrid_mod: use monotonic second_order conservative exchange"
     else
        write(out_unit,*)"NOTE from xgrid_mod: use second_order conservative exchange"
     endif
     remap_method = SECOND_ORDER
  case default
     call error_mesg('xgrid_mod', ' nml interp_method = ' //trim(interp_method)// &
      ' is not a valid namelist option', FATAL)
  end select

  if(xgrid_clocks_on) then
     id_put_1_to_xgrid_order_1 = mpp_clock_id("put_1_to_xgrid_order_1", flags=MPP_CLOCK_SYNC)
     id_put_1_to_xgrid_order_2 = mpp_clock_id("put_1_to_xgrid_order_2", flags=MPP_CLOCK_SYNC)
     id_get_1_from_xgrid       = mpp_clock_id("get_1_from_xgrid", flags=MPP_CLOCK_SYNC) 
     id_get_1_from_xgrid_repro = mpp_clock_id("get_1_from_xgrid_repro", flags=MPP_CLOCK_SYNC)
     id_get_2_from_xgrid       = mpp_clock_id("get_2_from_xgrid", flags=MPP_CLOCK_SYNC) 
     id_put_2_to_xgrid         = mpp_clock_id("put_2_to_xgrid", flags=MPP_CLOCK_SYNC) 
     id_setup_xmap             = mpp_clock_id("setup_xmap", flags=MPP_CLOCK_SYNC) 
     id_set_comm               = mpp_clock_id("set_comm") 
     id_regen                  = mpp_clock_id("regen") 
     id_conservation_check     = mpp_clock_id("conservation_check")
     id_load_xgrid             = mpp_clock_id("load_xgrid") 
     id_load_xgrid1            = mpp_clock_id("load_xgrid1")
     id_load_xgrid2            = mpp_clock_id("load_xgrid2")
     id_load_xgrid3            = mpp_clock_id("load_xgrid3")
     id_load_xgrid4            = mpp_clock_id("load_xgrid4")
     id_load_xgrid5            = mpp_clock_id("load_xgrid5")
  endif

  remapping_method = remap_method

end subroutine xgrid_init
! </SUBROUTINE>

!#######################################################################

subroutine load_xgrid (xmap, grid, grid_file, grid1_id, grid_id, tile1, tile2, use_higher_order)
type(xmap_type), intent(inout)         :: xmap
type(grid_type), intent(inout)         :: grid
character(len=*), intent(in)           :: grid_file
character(len=3), intent(in)           :: grid1_id, grid_id
integer,          intent(in)           :: tile1, tile2
logical,        intent(in)             :: use_higher_order

  integer, pointer,       dimension(:)   :: i1=>NULL(), j1=>NULL()
  integer, pointer,       dimension(:)   :: i2=>NULL(), j2=>NULL()            
  real,    pointer,       dimension(:)   :: di=>NULL(), dj=>NULL()
  real,    pointer,       dimension(:)   :: area =>NULL()     
  integer, pointer,       dimension(:)   :: i1_tmp=>NULL(), j1_tmp=>NULL()
  integer, pointer,       dimension(:)   :: i2_tmp=>NULL(), j2_tmp=>NULL()
  real,    pointer,       dimension(:)   :: di_tmp=>NULL(), dj_tmp=>NULL()
  real,    pointer,       dimension(:)   :: area_tmp =>NULL()
  integer, pointer,       dimension(:)   :: i1_side1=>NULL(), j1_side1=>NULL()
  integer, pointer,       dimension(:)   :: i2_side1=>NULL(), j2_side1=>NULL()            
  real,    pointer,       dimension(:)   :: di_side1=>NULL(), dj_side1=>NULL()
  real,    pointer,       dimension(:)   :: area_side1 =>NULL()    

  real,    allocatable, dimension(:,:) :: tmp
  real,    allocatable, dimension(:)   :: send_buffer, recv_buffer
  type (grid_type),   pointer, save    :: grid1 =>NULL()
  integer                              :: l, ll, ll_repro, p, siz(4), nxgrid, size_prev
  type(xcell_type),   allocatable      :: x_local(:)
  integer                              :: size_repro, out_unit
  logical                              :: scale_exist = .false.
  logical                              :: is_distribute = .false.
  real,    allocatable,   dimension(:) :: scale
  real                                 :: garea
  integer                              :: npes, isc, iec, nxgrid_local, pe, nxgrid_local_orig
  integer                              :: nxgrid1, nxgrid2, nset1, nset2, ndivs, cur_ind
  integer                              :: pos, nsend, nrecv, l1, l2, n, mypos, m
  integer                              :: start(4), nread(4)
  logical                              :: found
  character(len=128)                   :: attvalue
  integer, dimension(0:xmap%npes-1)    :: pelist
  logical, dimension(0:xmap%npes-1)    :: subset_rootpe
  integer, dimension(0:xmap%npes-1)    :: nsend1, nsend2, nrecv1, nrecv2
  integer, dimension(0:xmap%npes-1)    :: send_cnt, recv_cnt
  integer, dimension(0:xmap%npes-1)    :: send_buffer_pos, recv_buffer_pos
  integer, dimension(0:xmap%npes-1)    :: ibegin, iend, pebegin, peend
  integer, dimension(2*xmap%npes)      :: ibuf1, ibuf2
  integer, dimension(0:xmap%npes-1)    :: pos_x, y2m1_size
  integer, allocatable,   dimension(:) :: y2m1_pe
  integer, pointer, save               :: iarray(:), jarray(:)
  integer, allocatable, save           :: pos_s(:)
  integer, pointer,       dimension(:) :: iarray2(:)=>NULL(), jarray2(:)=>NULL()
  logical                              :: last_grid
  integer                              :: nxgrid1_old

  scale_exist = .false.
  grid1 => xmap%grids(1)
  out_unit = stdout()
  npes     = xmap%npes
  pe       = mpp_pe()
  mypos = mpp_pe()-mpp_root_pe()

  call mpp_get_current_pelist(pelist)
  !--- make sure npes = pelist(npes-1) - pelist(0) + 1
  if( npes .NE. pelist(npes-1) - pelist(0) + 1 ) then
     print*, "npes =", npes, ", pelist(npes-1)=", pelist(npes-1), ", pelist(0)=", pelist(0)
     call error_mesg('xgrid_mod', 'npes .NE. pelist(npes-1) - pelist(0)', FATAL)
  endif

  select case(xmap%version)
  case(VERSION1)
     call field_size(grid_file, 'AREA_'//grid1_id//'x'//grid_id, siz)
     nxgrid = siz(1);
     if(nxgrid .LE. 0) return
  case(VERSION2)
     !--- max_size is the exchange grid size between super grid.
     nxgrid = get_mosaic_xgrid_size(grid_file)
     if(nxgrid .LE. 0) return
  end select

  !--- define a domain to read exchange grid.
  if(nxgrid > npes) then
     ndivs = npes
     if(nsubset >0 .AND. nsubset < npes) ndivs = nsubset     
     call mpp_compute_extent( 1, nxgrid, ndivs, ibegin, iend)
     if(npes == ndivs) then
        p = mpp_pe()-mpp_root_pe()
        isc = ibegin(p)
        iec = iend(p)
        subset_rootpe(:) = .true.
     else
        isc = 0; iec = -1
        call mpp_compute_extent(pelist(0), pelist(npes-1), ndivs, pebegin, peend)
        do n = 0, ndivs-1
           if(pe == pebegin(n)) then
              isc = ibegin(n)
              iec = iend(n)
              exit
           endif
        enddo
        cur_ind = 0
        subset_rootpe(:) = .false.

        do n = 0, npes-1
           if(pelist(n) == pebegin(cur_ind)) then
              subset_rootpe(n) = .true.
              cur_ind = cur_ind+1
              if(cur_ind == ndivs) exit
           endif
        enddo
     endif
     is_distribute = .true.
  else
     is_distribute = .false.
     isc = 1; iec = nxgrid
  endif

  nset1 = 5
  nset2 = 5
  if(use_higher_order) then
     nset1 = nset1 + 2
     nset2 = nset2 + 2
  end if
  if(scale_exist) nset2 = nset1 + 1

  call mpp_clock_begin(id_load_xgrid1)
  if(iec .GE. isc) then
     nxgrid_local = iec - isc + 1
     allocate(i1_tmp(isc:iec), j1_tmp(isc:iec), i2_tmp(isc:iec), j2_tmp(isc:iec), area_tmp(isc:iec) )
     if(use_higher_order) allocate(di_tmp(isc:iec), dj_tmp(isc:iec))

     start = 1; nread = 1

     select case(xmap%version)
     case(VERSION1)
        start(1) = isc; nread(1) = nxgrid_local
        allocate(tmp(nxgrid_local,1))
        call read_data(grid_file, 'I_'//grid1_id//'_'//grid1_id//'x'//grid_id, tmp, start, nread, no_domain=.TRUE.)
        i1_tmp = tmp(:,1)
        call read_data(grid_file, 'J_'//grid1_id//'_'//grid1_id//'x'//grid_id, tmp, start, nread, no_domain=.TRUE.)
        j1_tmp = tmp(:,1)
        call read_data(grid_file, 'I_'//grid_id//'_'//grid1_id//'x'//grid_id, tmp, start, nread, no_domain=.TRUE.)
        i2_tmp = tmp(:,1)
        call read_data(grid_file, 'J_'//grid_id//'_'//grid1_id//'x'//grid_id, tmp, start, nread, no_domain=.TRUE.)
        j2_tmp = tmp(:,1)
        call read_data(grid_file, 'AREA_'//grid1_id//'x'//grid_id, tmp, start, nread, no_domain=.TRUE.)
        area_tmp = tmp(:,1)
        if(use_higher_order) then
           call read_data(grid_file, 'DI_'//grid1_id//'x'//grid_id, tmp, start, nread, no_domain=.TRUE.)
           di_tmp = tmp(:,1)
           call read_data(grid_file, 'DJ_'//grid1_id//'x'//grid_id, tmp, start, nread, no_domain=.TRUE.)
           dj_tmp = tmp(:,1)
        end if
        deallocate(tmp)
     case(VERSION2)
        nread(1) = 2; start(2) = isc; nread(2) = nxgrid_local
        allocate(tmp(2, isc:iec))
        call read_data(grid_file, "tile1_cell", tmp, start, nread, no_domain=.TRUE.)
        i1_tmp(isc:iec) = tmp(1, isc:iec)
        j1_tmp(isc:iec) = tmp(2, isc:iec)
        call read_data(grid_file, "tile2_cell", tmp, start, nread, no_domain=.TRUE.)
        i2_tmp(isc:iec) = tmp(1, isc:iec)
        j2_tmp(isc:iec) = tmp(2, isc:iec)
        if(use_higher_order) then
           call read_data(grid_file, "tile1_distance", tmp, start, nread, no_domain=.TRUE.)
           di_tmp(isc:iec) = tmp(1, isc:iec)
           dj_tmp(isc:iec) = tmp(2, isc:iec)
        end if
        start = 1; nread = 1
        start(1) = isc; nread(1) = nxgrid_local
        deallocate(tmp)
        allocate(tmp(isc:iec,1) )
        call read_data(grid_file, "xgrid_area", tmp(:,1:1), start, nread, no_domain=.TRUE.)
        ! check the units of "xgrid_area 
        call get_var_att_value(grid_file, "xgrid_area", "units", attvalue)
        if( trim(attvalue) == 'm2' ) then
           garea = 4.0*PI*RADIUS*RADIUS;
           area_tmp = tmp(:,1)/garea
        else if( trim(attvalue) == 'none' ) then
           area_tmp = tmp(:,1)
        else 
           call error_mesg('xgrid_mod', 'In file '//trim(grid_file)//', xgrid_area units = '// &
                trim(attvalue)//' should be "m2" or "none"', FATAL)
        endif

        !--- if field "scale" exist, read this field. Normally this 
        !--- field only exist in landXocean exchange grid cell.
        if(grid1_id == 'LND' .AND. grid_id == 'OCN') then
           if(field_exist(grid_file, "scale")) then
              allocate(scale(isc:iec))
              write(out_unit, *)"NOTE from load_xgrid(xgrid_mod): field 'scale' exist in the file "// &
                   trim(grid_file)//", this field will be read and the exchange grid cell area will be multiplied by scale"
              call read_data(grid_file, "scale", tmp, start, nread, no_domain=.TRUE.)
              scale = tmp(:,1)
              scale_exist = .true.
           endif
        endif
        deallocate(tmp)
     end select

     !---z1l: The following change is for the situation that some processor is masked out.
     !---loop through all the pe to see if side 1 and side of each exchange grid is on some processor
     nxgrid_local_orig = nxgrid_local
     allocate(i1(isc:iec), j1(isc:iec), i2(isc:iec), j2(isc:iec), area(isc:iec) )
     if(use_higher_order) allocate(di(isc:iec), dj(isc:iec))
     pos = isc-1

     do l = isc, iec
        found = .false.
        !--- first check if the exchange grid is on one of side 1 processor
        do p = 0, npes - 1
           if(grid1%tile(p) == tile1) then
              if (in_box(i1_tmp(l), j1_tmp(l), grid1%is(p), grid1%ie(p), grid1%js(p), grid1%je(p)))  then
                 found = .true.
                 exit
              endif
           endif
        enddo
        !--- Then check if the exchange grid is on one of side 2 processor
        if( found ) then
           do p = 0, npes - 1
              if(grid%tile(p) == tile2) then
                 if (in_box(i2_tmp(l), j2_tmp(l), grid%is(p), grid%ie(p), grid%js(p), grid%je(p)))  then
                    pos = pos+1
                    i1(pos) = i1_tmp(l)
                    j1(pos) = j1_tmp(l)
                    i2(pos) = i2_tmp(l)
                    j2(pos) = j2_tmp(l)
                    area(pos) = area_tmp(l)
                    if(use_higher_order) then
                       di(pos) = di_tmp(l)
                       dj(pos) = dj_tmp(l)
                    endif
                    exit
                 endif
              endif
           enddo
        endif
     enddo

     deallocate(i1_tmp, i2_tmp, j1_tmp, j2_tmp, area_tmp)
     if(use_higher_order) deallocate( di_tmp, dj_tmp)
     iec = pos
     if(iec .GE. isc) then
        nxgrid_local = iec - isc + 1
     else
        nxgrid_local = 0
     endif


  else
     nxgrid_local = 0
     nxgrid_local_orig = 0
  endif

  call mpp_clock_end(id_load_xgrid1)

  if(is_distribute) then
     !--- Since the xgrid is distributed according to side 2 grid. Send all the xgrid to its own side 2. 
     !--- Also need to send the xgrid to its own side 1 for the reproducing ability between processor count.
     !--- first find out number of points need to send to other pe and fill the send buffer.
     nsend1(:) = 0; nrecv1(:) = 0
     nsend2(:) = 0; nrecv2(:) = 0
     ibuf1(:)= 0; ibuf2(:)= 0

     call mpp_clock_begin(id_load_xgrid2)

     if(nxgrid_local>0) then
        allocate( send_buffer(nxgrid_local * (nset1+nset2)) )
        pos = 0
        do p = 0, npes - 1
           send_buffer_pos(p) = pos
           if(grid%tile(p) == tile2) then
              do l = isc, iec
                 if (in_box(i2(l), j2(l), grid%is(p), grid%ie(p), grid%js(p), grid%je(p)))  then
                    nsend2(p) = nsend2(p) + 1
                    send_buffer(pos+1) = i1(l)
                    send_buffer(pos+2) = j1(l)
                    send_buffer(pos+3) = i2(l)
                    send_buffer(pos+4) = j2(l)
                    send_buffer(pos+5) = area(l)
                    if(use_higher_order) then
                       send_buffer(pos+6) = di(l)
                       send_buffer(pos+7) = dj(l)
                    endif
                    if(scale_exist) send_buffer(pos+nset2) = scale(l)
                    pos = pos + nset2
                 endif
              enddo
           endif
           if(grid1%tile(p) == tile1) then
              do l = isc, iec
                 if (in_box(i1(l), j1(l), grid1%is(p), grid1%ie(p), grid1%js(p), grid1%je(p)))  then
                    nsend1(p) = nsend1(p) + 1
                    send_buffer(pos+1) = i1(l)
                    send_buffer(pos+2) = j1(l)
                    send_buffer(pos+3) = i2(l)
                    send_buffer(pos+4) = j2(l)
                    send_buffer(pos+5) = area(l)
                    if(use_higher_order) then
                       send_buffer(pos+6) = di(l)
                       send_buffer(pos+7) = dj(l)
                    endif
                    pos = pos + nset1
                 endif
              enddo
           endif
        enddo
     endif
     call mpp_clock_end(id_load_xgrid2)

     !--- send the size of the data on side 1 to be sent over.
     call mpp_clock_begin(id_load_xgrid3)

     if (do_alltoall) then
        do p = 0, npes-1
                ibuf1(2*p+1) = nsend1(p)
                ibuf1(2*p+2) = nsend2(p)
        enddo
        call mpp_alltoall(ibuf1, 2, ibuf2, 2)
     else
        do n = 0, npes-1
           p = mod(mypos+npes-n, npes)
           if(.not. subset_rootpe(p)) cycle
           call mpp_recv( ibuf2(2*p+1), glen=2, from_pe=pelist(p), block=.FALSE., tag=COMM_TAG_1)
        enddo

        if(nxgrid_local_orig>0) then
           do n = 0, npes-1
              p = mod(mypos+n, npes)
              ibuf1(2*p+1) = nsend1(p)
              ibuf1(2*p+2) = nsend2(p)
              call mpp_send( ibuf1(2*p+1), plen=2, to_pe=pelist(p), tag=COMM_TAG_1)
           enddo
        endif
        call mpp_sync_self(check=EVENT_RECV)
     endif
     do p = 0, npes-1
        nrecv1(p) = ibuf2(2*p+1)
        nrecv2(p) = ibuf2(2*p+2)
     enddo
    
     if(.not. do_alltoall) call mpp_sync_self()
     call mpp_clock_end(id_load_xgrid3)
     call mpp_clock_begin(id_load_xgrid4)
     pos = 0
     do p = 0, npes - 1
        recv_buffer_pos(p) = pos
        pos = pos + nrecv1(p) * nset1 + nrecv2(p) * nset2
     end do

     !--- now get the data
     nxgrid1 = sum(nrecv1)
     nxgrid2 = sum(nrecv2)
     if(nxgrid1>0 .OR. nxgrid2>0) allocate(recv_buffer(nxgrid1*nset1+nxgrid2*nset2))

     if (do_alltoallv) then
        ! Construct the send and receive counters
        send_cnt(:) = nset1 * nsend1(:) + nset2 * nsend2(:)
        recv_cnt(:) = nset1 * nrecv1(:) + nset2 * nrecv2(:)

        call mpp_alltoall(send_buffer, send_cnt, send_buffer_pos, &
                          recv_buffer, recv_cnt, recv_buffer_pos)
     else
        do n = 0, npes-1
           p = mod(mypos+npes-n, npes)
           nrecv = nrecv1(p)*nset1+nrecv2(p)*nset2
           if(nrecv==0) cycle
           pos = recv_buffer_pos(p)
           call mpp_recv(recv_buffer(pos+1), glen=nrecv, from_pe=pelist(p), &
                         block=.FALSE., tag=COMM_TAG_2)
        end do

        do n = 0, npes-1
           p = mod(mypos+n, npes)
           nsend = nsend1(p)*nset1 + nsend2(p)*nset2
           if(nsend==0) cycle
           pos = send_buffer_pos(p)
           call mpp_send(send_buffer(pos+1), plen=nsend, to_pe=pelist(p), &
                         tag=COMM_TAG_2)
        end do
        call mpp_sync_self(check=EVENT_RECV)
     end if
     call mpp_clock_end(id_load_xgrid4)
     !--- unpack buffer.
     if( nxgrid_local>0) then
        deallocate(i1,j1,i2,j2,area)
     endif 

     allocate(i1(nxgrid2), j1(nxgrid2))
     allocate(i2(nxgrid2), j2(nxgrid2))
     allocate(area(nxgrid2))
     allocate(i1_side1(nxgrid1), j1_side1(nxgrid1))
     allocate(i2_side1(nxgrid1), j2_side1(nxgrid1))
     allocate(area_side1(nxgrid1))
     if(use_higher_order) then
        if(nxgrid_local>0) deallocate(di,dj)
        allocate(di      (nxgrid2), dj      (nxgrid2))
        allocate(di_side1(nxgrid1), dj_side1(nxgrid1))
     endif
     if(scale_exist) then
        if(nxgrid_local>0)deallocate(scale)
        allocate(scale(nxgrid2))
     endif
     pos = 0
     l1 = 0; l2 = 0
     do p = 0,npes-1
        do n = 1, nrecv2(p)
           l2 = l2+1
           i1(l2) = recv_buffer(pos+1)
           j1(l2) = recv_buffer(pos+2)
           i2(l2) = recv_buffer(pos+3)
           j2(l2) = recv_buffer(pos+4)
           area(l2) = recv_buffer(pos+5)
           if(use_higher_order) then
              di(l2) = recv_buffer(pos+6)
              dj(l2) = recv_buffer(pos+7)
           endif
           if(scale_exist)scale(l2) = recv_buffer(pos+nset2)
           pos = pos + nset2
        enddo
        do n = 1, nrecv1(p)
           l1 = l1+1
           i1_side1(l1) = recv_buffer(pos+1)
           j1_side1(l1) = recv_buffer(pos+2)
           i2_side1(l1) = recv_buffer(pos+3)
           j2_side1(l1) = recv_buffer(pos+4)
           area_side1(l1) = recv_buffer(pos+5)
           if(use_higher_order) then
              di_side1(l1) = recv_buffer(pos+6)
              dj_side1(l1) = recv_buffer(pos+7)
           endif
           pos = pos + nset1
        enddo
     enddo
     call mpp_sync_self()
     if(allocated(send_buffer)) deallocate(send_buffer)
     if(allocated(recv_buffer)) deallocate(recv_buffer)

  else
     nxgrid1 = nxgrid
     nxgrid2 = nxgrid
     i1_side1 => i1; j1_side1 => j1
     i2_side1 => i2; j2_side1 => j2
     area_side1 => area
     if(use_higher_order) then
        di_side1 => di
        dj_side1 => dj
     endif
  endif

  call mpp_clock_begin(id_load_xgrid5)


  size_prev = grid%size

  if(grid%tile_me == tile2) then
     do l=1,nxgrid2
        if (in_box(i2(l), j2(l), grid%is_me, grid%ie_me, grid%js_me, grid%je_me) ) then
           grid%size = grid%size + 1
           ! exclude the area overlapped with parent grid 
           if( grid1_id .NE. "ATM" .OR. tile1 .NE. tile_parent .OR.  &
                   .NOT. in_box(i1(l), j1(l), is_parent, ie_parent, js_parent, je_parent) ) then
              grid%area(i2(l),j2(l)) = grid%area(i2(l),j2(l))+area(l)
           endif
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
        grid%x%di = 0.0; grid%x%dj = 0.0
     end if
  end if

  ll = size_prev
  if( grid%tile_me == tile2 ) then ! me is tile2
     do l=1,nxgrid2
        if (in_box(i2(l), j2(l), grid%is_me, grid%ie_me, grid%js_me, grid%je_me)) then
           ! insert in this grids cell pattern list and add area to side 2 area
           ll = ll + 1
           grid%x(ll)%i1   = i1(l); grid%x(ll)%i2   = i2(l)
           grid%x(ll)%j1   = j1(l); grid%x(ll)%j2   = j2(l)
           grid%x(ll)%tile = tile1
           grid%x(ll)%area = area(l)
           if(scale_exist) then
              grid%x(ll)%scale = scale(l)
           else
              grid%x(ll)%scale = 1.0
           endif
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

  if(grid%id == xmap%grids(size(xmap%grids(:)))%id) then
     last_grid = .true.
  else 
     last_grid = .false.
  endif

  size_repro = 0
  if(grid1%tile_me == tile1) then
     if(associated(iarray)) then
        nxgrid1_old = size(iarray(:))
     else
        nxgrid1_old = 0
     endif

     allocate(y2m1_pe(nxgrid1))
     if(.not. last_grid ) allocate(pos_s(0:xmap%npes-1))
     y2m1_pe = -1
     if(nxgrid1_old > 0) then
        do p=0,xmap%npes-1
           y2m1_size(p) = xmap%your2my1_size(p)
        enddo
     else
        y2m1_size = 0
     endif

     do l=1,nxgrid1
        if (in_box(i1_side1(l), j1_side1(l), grid1%is_me, grid1%ie_me, grid1%js_me, grid1%je_me) ) then
           grid1%area(i1_side1(l),j1_side1(l)) = grid1%area(i1_side1(l),j1_side1(l))+area_side1(l)
           do p=0,xmap%npes-1
              if (grid%tile(p) == tile2) then
                 if (in_box(i2_side1(l), j2_side1(l), grid%is(p), grid%ie(p), grid%js(p), grid%je(p)))  then
                    xmap%your2my1(p) = .true.
                    y2m1_pe(l) = p
                    y2m1_size(p) = y2m1_size(p) + 1
                 endif
              endif
           enddo
           size_repro = size_repro + 1
        endif
     enddo
     pos_x = 0
     do p = 1, npes-1
        pos_x(p) = pos_x(p-1) + y2m1_size(p-1)
     enddo

     if(.not. last_grid) pos_s(:) = pos_x(:)

     if(nxgrid1_old > 0) then
        y2m1_size(:) = xmap%your2my1_size(:)
        iarray2 => iarray
        jarray2 => jarray
        allocate(iarray(nxgrid1+nxgrid1_old), jarray(nxgrid1+nxgrid1_old))
        ! copy the i-j index
        do p=0,xmap%npes-1
           do n = 1, xmap%your2my1_size(p)
              iarray(pos_x(p)+n) = iarray2(pos_s(p)+n)
              jarray(pos_x(p)+n) = jarray2(pos_s(p)+n)
           enddo
        enddo    
        deallocate(iarray2, jarray2)        
     else
        allocate(iarray(nxgrid1), jarray(nxgrid1))
        iarray(:) = 0
        jarray(:) = 0
        y2m1_size(:) = 0
     endif

     do l=1,nxgrid1
        p = y2m1_pe(l)
        if(p<0) cycle
        found = .false.
        if(y2m1_size(p) > 0) then
           pos = pos_x(p)+y2m1_size(p)
           if( i1_side1(l) == iarray(pos) .AND. j1_side1(l) == jarray(pos) ) then
              found = .true.
           else
              !---may need to replace with a fast search algorithm
              do n = 1, y2m1_size(p)
                 pos = pos_x(p)+n
                 if(i1_side1(l) == iarray(pos) .AND. j1_side1(l) == jarray(pos)) then
                    found = .true.
                    exit
                 endif
              enddo
           endif
        endif
        if( (.NOT. found) .OR. monotonic_exchange ) then
           y2m1_size(p) = y2m1_size(p)+1
           pos = pos_x(p)+y2m1_size(p)
           iarray(pos) = i1_side1(l)
           jarray(pos) = j1_side1(l)
        endif
     end do
     xmap%your2my1_size(:) =  y2m1_size(:)
     deallocate(y2m1_pe)
     if(last_grid) then
        deallocate(iarray, jarray)
        if(allocated(pos_s)) deallocate(pos_s)
     end if
  end if

  if (grid1%tile_me == tile1 .and. size_repro > 0) then
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
        grid%x_repro%di = 0.0; grid%x_repro%dj = 0.0
     end if
     do l=1,nxgrid1
        if (in_box(i1_side1(l),j1_side1(l), grid1%is_me,grid1%ie_me, grid1%js_me,grid1%je_me) ) then
           ll_repro = ll_repro + 1
           grid%x_repro(ll_repro)%i1   = i1_side1(l); grid%x_repro(ll_repro)%i2   = i2_side1(l)
           grid%x_repro(ll_repro)%j1   = j1_side1(l); grid%x_repro(ll_repro)%j2   = j2_side1(l)
           grid%x_repro(ll_repro)%tile = tile1
           grid%x_repro(ll_repro)%area = area_side1(l)
           if(use_higher_order) then
              grid%x_repro(ll_repro)%di  = di_side1(l)
              grid%x_repro(ll_repro)%dj  = dj_side1(l)
           end if

           do p=0,xmap%npes-1
              if(grid%tile(p) == tile2) then
                 if (in_box(i2_side1(l), j2_side1(l), grid%is(p), grid%ie(p), &
                      grid%js(p), grid%je(p))) then
                    grid%x_repro(ll_repro)%pe = p + xmap%root_pe
                 end if
              end if
           end do
        end if ! make_exchange_reproduce
     end do
  end if

  deallocate(i1, j1, i2, j2, area)
  if(use_higher_order) deallocate(di, dj)
  if(scale_exist) deallocate(scale)
  if(is_distribute) then
     deallocate(i1_side1, j1_side1, i2_side1, j2_side1, area_side1)
     if(use_higher_order) deallocate(di_side1, dj_side1)
  endif  

  i1=>NULL(); j1=>NULL(); i2=>NULL(); j2=>NULL()
  call mpp_clock_end(id_load_xgrid5)



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
  integer                  :: start(4), nread(4), isc2, iec2, jsc2, jec2

  d2r = PI/180.0

  call mpp_get_compute_domain(grid%domain, is, ie, js, je)

  select case(grid_version)
  case(VERSION1)
     allocate(grid%lon(grid%im), grid%lat(grid%jm))
     if(grid_id == 'ATM') then
        call read_data(grid_file, 'xta', lonb)
        call read_data(grid_file, 'yta', latb)

        if(.not. allocated(AREA_ATM_MODEL)) then
           allocate(AREA_ATM_MODEL(is:ie, js:je))
           call get_area_elements(grid_file, 'AREA_ATM_MODEL', grid%domain, AREA_ATM_MODEL)
        endif
        if(.not. allocated(AREA_ATM_SPHERE)) then
           allocate(AREA_ATM_SPHERE(is:ie, js:je))
           call get_area_elements(grid_file, 'AREA_ATM', grid%domain, AREA_ATM_SPHERE)
        endif
     else if(grid_id == 'LND') then
        call read_data(grid_file, 'xtl', lonb)
        call read_data(grid_file, 'ytl', latb)
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

     if( grid_id == 'LND' .or. grid_id == 'ATM'  .or. grid_id == 'WAV' ) then
        isc2 = 2*grid%is_me-1; iec2 = 2*grid%ie_me+1
        jsc2 = 2*grid%js_me-1; jec2 = 2*grid%je_me+1
        allocate(tmpx(isc2:iec2, jsc2:jec2) )
        allocate(tmpy(isc2:iec2, jsc2:jec2) )   
        start = 1; nread = 1          
        start(1) = isc2; nread(1) = iec2 - isc2 + 1
        start(2) = jsc2; nread(2) = jec2 - jsc2 + 1 
        call read_data(grid_file, 'x', tmpx, start, nread, no_domain=.TRUE.)
        call read_data(grid_file, 'y', tmpy, start, nread, no_domain=.TRUE.)      
        if(is_lat_lon(tmpx, tmpy) ) then
           deallocate(tmpx, tmpy)
           start = 1; nread = 1
           start(2) = 2; nread(1) = nlon*2+1
           allocate(tmpx(nlon*2+1, 1), tmpy(1, nlat*2+1))
           call read_data(grid_file, "x", tmpx, start, nread, no_domain=.TRUE.)    
           allocate(grid%lon(grid%im), grid%lat(grid%jm))
           do i = 1, grid%im
              grid%lon(i) = tmpx(2*i,1) * d2r
           end do
           start = 1; nread = 1
           start(1) = 2; nread(2) = nlat*2+1
           call read_data(grid_file, "y", tmpy, start, nread, no_domain=.TRUE.)    
           do j = 1, grid%jm
              grid%lat(j) = tmpy(1, 2*j) * d2r
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
        deallocate(tmpx, tmpy)
     end if     
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
     data = -1.0
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
!     call setup_xmap(xmap, grid_ids, grid_domains, grid_file, atm_grid)
!   </TEMPLATE>

!   <IN NAME="grid_ids" TYPE="character(len=3)" DIM="(:)"> </IN>
!   <IN NAME="grid_domains" TYPE="type(Domain2d)" DIM="(:)"> </IN>
!   <IN NAME="grid_file" TYPE="character(len=*)" > </IN>
!   <IN NAME="atmos_grid" TYPE="type(grid_box_type),optional" > </IN>
!   <OUT NAME="xmap" TYPE="xmap_type"  > </OUT>

subroutine setup_xmap(xmap, grid_ids, grid_domains, grid_file, atm_grid)
  type (xmap_type),                        intent(inout) :: xmap
  character(len=3), dimension(:),            intent(in ) :: grid_ids
  type(Domain2d), dimension(:),              intent(in ) :: grid_domains
  character(len=*),                          intent(in ) :: grid_file
  type(grid_box_type), optional,             intent(in ) :: atm_grid

  integer :: g,     p, send_size, recv_size, i, siz(4)
  integer :: unit, nxgrid_file, i1, i2, i3, tile1, tile2, j
  integer :: nxc, nyc, out_unit
  type (grid_type), pointer, save :: grid =>NULL(), grid1 =>NULL()
  real, dimension(3) :: xxx
  real, dimension(:,:), allocatable   :: check_data
  real, dimension(:,:,:), allocatable :: check_data_3D
  character(len=256)                  :: xgrid_file, xgrid_name
  character(len=256)                  :: tile_file, mosaic_file
  character(len=256)                  :: mosaic1, mosaic2, contact
  character(len=256)                  :: tile1_name, tile2_name
  character(len=256),     allocatable :: tile1_list(:), tile2_list(:)
  integer                             :: npes, npes2
  integer,                allocatable :: pelist(:)
  type(domain2d), save                :: domain2
  logical :: use_higher_order = .false.

  call mpp_clock_begin(id_setup_xmap)

  if(interp_method .ne. 'first_order')  use_higher_order = .true.

  out_unit = stdout()
  xmap%me   = mpp_pe  ()
  xmap%npes = mpp_npes()
  xmap%root_pe = mpp_root_pe()

  allocate( xmap%grids(1:size(grid_ids(:))) )

  allocate ( xmap%your1my2(0:xmap%npes-1), xmap%your2my1(0:xmap%npes-1) )
  allocate ( xmap%your2my1_size(0:xmap%npes-1) )

  xmap%your1my2 = .false.; xmap%your2my1 = .false.;
  xmap%your2my1_size = 0

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

  call mpp_clock_begin(id_load_xgrid)
  do g=1,size(grid_ids(:))
     grid => xmap%grids(g)
     if (g==1) grid1 => xmap%grids(g)
     grid%id     = grid_ids    (g)
     grid%domain = grid_domains(g)
     grid%on_this_pe = mpp_domain_is_initialized(grid_domains(g))
     allocate ( grid%is(0:xmap%npes-1), grid%ie(0:xmap%npes-1) )
     allocate ( grid%js(0:xmap%npes-1), grid%je(0:xmap%npes-1) )
     allocate ( grid%tile(0:xmap%npes-1) )
     grid%npes   = 0
     grid%ni = 0
     grid%nj = 0
     grid%is = 0
     grid%ie = -1
     grid%js = 0
     grid%je = -1
     grid%tile = -1

     select case(xmap%version)
     case(VERSION1)
        grid%ntile = 1
     case(VERSION2)
        call read_data(grid_file, lowercase(grid_ids(g))//'_mosaic_file', mosaic_file)
        grid%ntile = get_mosaic_ntiles('INPUT/'//trim(mosaic_file))
     end select

     if( g == 1 .AND. grid_ids(1) == 'ATM' ) then 
        if( .NOT. grid%on_this_pe ) call error_mesg('xgrid_mod', 'ATM domain is not defined on some processor' ,FATAL)
     endif
     grid%npes =  mpp_get_domain_npes(grid%domain)
     if( xmap%npes > grid%npes .AND. g == 1 .AND. grid_ids(1) == 'ATM' ) then
        call mpp_broadcast_domain(grid%domain, domain2)
     else if(xmap%npes > grid%npes) then
        call mpp_broadcast_domain(grid%domain)
        grid%npes =  mpp_get_domain_npes(grid%domain)
     endif

     npes = grid%npes
     allocate(grid%pelist(0:npes-1))
     call mpp_get_domain_pelist(grid%domain, grid%pelist)
     grid%root_pe = mpp_get_domain_root_pe(grid%domain)

     call mpp_get_data_domain(grid%domain, grid%isd_me, grid%ied_me, grid%jsd_me, grid%jed_me, &
                              xsize=grid%nxd_me, ysize=grid%nyd_me)
     call mpp_get_global_domain(grid%domain, xsize=grid%ni, ysize=grid%nj)

     if( grid%root_pe == xmap%root_pe ) then
        call mpp_get_compute_domains(grid%domain,  xbegin=grid%is(0:npes-1), xend=grid%ie(0:npes-1), &
                                     ybegin=grid%js(0:npes-1), yend=grid%je(0:npes-1) )
        call mpp_get_tile_list(grid%domain, grid%tile(0:npes-1))
        if( xmap%npes > npes .AND. g == 1 .AND. grid_ids(1) == 'ATM' ) then
           call mpp_get_compute_domains(domain2, xbegin=grid%is(npes:xmap%npes-1), xend=grid%ie(npes:xmap%npes-1), &
                                        ybegin=grid%js(npes:xmap%npes-1), yend=grid%je(npes:xmap%npes-1) )
           call mpp_get_tile_list(domain2, grid%tile(npes:xmap%npes-1))
        endif
     else
        npes2 = xmap%npes-npes
        call mpp_get_compute_domains(domain2,  xbegin=grid%is(0:npes2-1), xend=grid%ie(0:npes2-1), &
                                     ybegin=grid%js(0:npes2-1), yend=grid%je(0:npes2-1) )
        call mpp_get_compute_domains(grid%domain, xbegin=grid%is(npes2:xmap%npes-1), xend=grid%ie(npes2:xmap%npes-1), &
                                     ybegin=grid%js(npes2:xmap%npes-1), yend=grid%je(npes2:xmap%npes-1) )
        call mpp_get_tile_list(domain2, grid%tile(0:npes2-1))
        call mpp_get_tile_list(grid%domain, grid%tile(npes2:xmap%npes-1))
     endif
     if( xmap%npes > grid%npes .AND. g == 1 .AND. grid_ids(1) == 'ATM' ) then
        call mpp_deallocate_domain(domain2)
     endif
     npes = grid%npes
     if(  g == 1 .AND. grid_ids(1) == 'ATM' ) npes = xmap%npes
     do p = 0, npes-1
        if(grid%tile(p) > grid%ntile .or. grid%tile(p) < 1) call error_mesg('xgrid_mod', &
                 'tile id should between 1 and ntile', FATAL)
     end do

     grid%im = grid%ni
     grid%jm = grid%nj
     call mpp_max(grid%ni)
     call mpp_max(grid%nj)
    
     grid%is_me => grid%is(xmap%me-xmap%root_pe); grid%ie_me => grid%ie(xmap%me-xmap%root_pe)
     grid%js_me => grid%js(xmap%me-xmap%root_pe); grid%je_me => grid%je(xmap%me-xmap%root_pe)
     grid%nxc_me = grid%ie_me - grid%is_me + 1
     grid%nyc_me = grid%je_me - grid%js_me + 1
     grid%tile_me => grid%tile(xmap%me-xmap%root_pe)

     grid%km = 1

     if( grid%on_this_pe ) then
        allocate( grid%area    (grid%is_me:grid%ie_me, grid%js_me:grid%je_me) )
        allocate( grid%area_inv(grid%is_me:grid%ie_me, grid%js_me:grid%je_me) )
        grid%area       = 0.0
        grid%size       = 0
        grid%size_repro = 0
     endif
    
     ! get the center point of the grid box
     select case(xmap%version)
     case(VERSION1)
        if( grid%npes .NE. xmap%npes ) then
           call error_mesg('xgrid_mod', ' grid%npes .NE. xmap%npes ', FATAL)
        endif
        call get_grid(grid, grid_ids(g), grid_file, xmap%version)
     case(VERSION2)
        allocate(pelist(0:xmap%npes-1))
        call mpp_get_current_pelist(pelist)
        if( grid%on_this_pe ) then
           call mpp_set_current_pelist(grid%pelist)
           call get_mosaic_tile_grid(tile_file, 'INPUT/'//trim(mosaic_file), grid%domain)
           call get_grid(grid, grid_ids(g), tile_file, xmap%version)
        endif
        call mpp_set_current_pelist(pelist)
        deallocate(pelist) 
        ! read the contact information from mosaic_file to check if atmosphere is nested model
        if( g == 1 .AND. grid_ids(1) == 'ATM' ) then
           nnest = get_nest_contact('INPUT/'//trim(mosaic_file), tile_nest, tile_parent, is_nest, &
                ie_nest, js_nest, je_nest, is_parent, ie_parent, js_parent, je_parent)         
        endif       
     end select

     if( use_higher_order .AND. grid%id == 'ATM') then
        if( nnest > 0 ) call error_mesg('xgrid_mod', 'second_order is not supported for nested coupler', FATAL)
        if( grid%is_latlon ) then
           call mpp_modify_domain(grid%domain, grid%domain_with_halo, whalo=1, ehalo=1, shalo=1, nhalo=1)
           call mpp_get_data_domain(grid%domain_with_halo, grid%isd_me, grid%ied_me, grid%jsd_me, grid%jed_me, &
                              xsize=grid%nxd_me, ysize=grid%nyd_me) 
        else
           if(.NOT. present(atm_grid)) call error_mesg('xgrid_mod', 'when first grid is "ATM", atm_grid should be present', FATAL)
           if(grid%is_me-grid%isd_me .NE. 1 .or. grid%ied_me-grid%ie_me .NE. 1 .or.               &
              grid%js_me-grid%jsd_me .NE. 1 .or. grid%jed_me-grid%je_me .NE. 1 ) call error_mesg( &
              'xgrid_mod', 'for non-latlon grid (cubic grid), the halo size should be 1 in all four direction', FATAL)
           if(.NOT.( ASSOCIATED(atm_grid%dx) .AND. ASSOCIATED(atm_grid%dy) .AND. ASSOCIATED(atm_grid%edge_w) .AND.    &
                ASSOCIATED(atm_grid%edge_e) .AND. ASSOCIATED(atm_grid%edge_s) .AND. ASSOCIATED(atm_grid%edge_n) .AND. &
                ASSOCIATED(atm_grid%en1) .AND. ASSOCIATED(atm_grid%en2) .AND. ASSOCIATED(atm_grid%vlon) .AND.         &
                ASSOCIATED(atm_grid%vlat) ) )  call error_mesg( &
            'xgrid_mod', 'for non-latlon grid (cubic grid), all the fields in atm_grid data type should be allocated', FATAL)
           nxc = grid%ie_me  - grid%is_me  + 1
           nyc = grid%je_me  - grid%js_me  + 1
           if(size(atm_grid%dx,1) .NE. nxc .OR. size(atm_grid%dx,2) .NE. nyc+1)               &
                call error_mesg('xgrid_mod', 'incorrect dimension size of atm_grid%dx', FATAL)
           if(size(atm_grid%dy,1) .NE. nxc+1 .OR. size(atm_grid%dy,2) .NE. nyc)               &
                call error_mesg('xgrid_mod', 'incorrect dimension sizeof atm_grid%dy', FATAL)
           if(size(atm_grid%area,1) .NE. nxc .OR. size(atm_grid%area,2) .NE. nyc)             &
                call error_mesg('xgrid_mod', 'incorrect dimension size of atm_grid%area', FATAL)
           if(size(atm_grid%edge_w(:)) .NE. nyc+1 .OR. size(atm_grid%edge_e(:)) .NE. nyc+1)    &
                call error_mesg('xgrid_mod', 'incorrect dimension size of atm_grid%edge_w/edge_e', FATAL)
           if(size(atm_grid%edge_s(:)) .NE. nxc+1 .OR. size(atm_grid%edge_n(:)) .NE. nxc+1)    &
                call error_mesg('xgrid_mod', 'incorrect dimension size of atm_grid%edge_s/edge_n', FATAL)
           if(size(atm_grid%en1,1) .NE. 3 .OR. size(atm_grid%en1,2) .NE. nxc .OR. size(atm_grid%en1,3) .NE. nyc+1) & 
                call error_mesg( 'xgrid_mod', 'incorrect dimension size of atm_grid%en1', FATAL)
           if(size(atm_grid%en2,1) .NE. 3 .OR. size(atm_grid%en2,2) .NE. nxc+1 .OR. size(atm_grid%en2,3) .NE. nyc) &
                call error_mesg( 'xgrid_mod', 'incorrect dimension size of atm_grid%en2', FATAL)
           if(size(atm_grid%vlon,1) .NE. 3 .OR. size(atm_grid%vlon,2) .NE. nxc .OR. size(atm_grid%vlon,3) .NE. nyc)   &
                call error_mesg('xgrid_mod', 'incorrect dimension size of atm_grid%vlon', FATAL)
           if(size(atm_grid%vlat,1) .NE. 3 .OR. size(atm_grid%vlat,2) .NE. nxc .OR. size(atm_grid%vlat,3) .NE. nyc)   &
                call error_mesg('xgrid_mod', 'incorrect dimension size of atm_grid%vlat', FATAL)
           allocate(grid%box%dx    (grid%is_me:grid%ie_me,   grid%js_me:grid%je_me+1 ))
           allocate(grid%box%dy    (grid%is_me:grid%ie_me+1, grid%js_me:grid%je_me   ))
           allocate(grid%box%area  (grid%is_me:grid%ie_me,   grid%js_me:grid%je_me   ))
           allocate(grid%box%edge_w(grid%js_me:grid%je_me+1))
           allocate(grid%box%edge_e(grid%js_me:grid%je_me+1))
           allocate(grid%box%edge_s(grid%is_me:grid%ie_me+1))
           allocate(grid%box%edge_n(grid%is_me:grid%ie_me+1))
           allocate(grid%box%en1   (3, grid%is_me:grid%ie_me,   grid%js_me:grid%je_me+1 ))
           allocate(grid%box%en2   (3, grid%is_me:grid%ie_me+1, grid%js_me:grid%je_me   ))
           allocate(grid%box%vlon  (3, grid%is_me:grid%ie_me,   grid%js_me:grid%je_me   ))
           allocate(grid%box%vlat  (3, grid%is_me:grid%ie_me,   grid%js_me:grid%je_me   ))
           grid%box%dx     = atm_grid%dx
           grid%box%dy     = atm_grid%dy
           grid%box%area   = atm_grid%area
           grid%box%edge_w = atm_grid%edge_w
           grid%box%edge_e = atm_grid%edge_e
           grid%box%edge_s = atm_grid%edge_s
           grid%box%edge_n = atm_grid%edge_n
           grid%box%en1    = atm_grid%en1
           grid%box%en2    = atm_grid%en2
           grid%box%vlon   = atm_grid%vlon
           grid%box%vlat   = atm_grid%vlat
        end if
     end if

     if (g>1) then
        if(grid%on_this_pe) then
           allocate( grid%frac_area(grid%is_me:grid%ie_me, grid%js_me:grid%je_me, grid%km) )
           grid%frac_area = 1.0
        endif

        ! load exchange cells, sum grid cell areas, set your1my2/your2my1
        select case(xmap%version)
        case(VERSION1)
           call load_xgrid (xmap, grid, grid_file, grid_ids(1), grid_ids(g), 1, 1, use_higher_order)
        case(VERSION2)
           select case(grid_ids(1))
           case( 'ATM' )
              xgrid_name = 'a'
           case( 'LND' )
              xgrid_name = 'l'
           case( 'WAV' )
              xgrid_name = 'w'
           case default 
              call error_mesg('xgrid_mod', 'grid_ids(1) should be ATM, LND or WAV', FATAL)
           end select
           select case(grid_ids(g))
           case( 'LND' )
              xgrid_name = trim(xgrid_name)//'Xl_file'
           case( 'OCN' )
              xgrid_name = trim(xgrid_name)//'Xo_file'
           case( 'WAV' )
              xgrid_name = trim(xgrid_name)//'Xw_file'
           case default 
              call error_mesg('xgrid_mod', 'grid_ids(g) should be LND, OCN or WAV', FATAL)
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
           if(field_exist(grid_file, xgrid_name)) then
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

                 call load_xgrid (xmap, grid, xgrid_file, grid_ids(1), grid_ids(g), tile1, tile2, &
                                  use_higher_order)
              end do
           endif
           deallocate(tile1_list, tile2_list)
        end select
        if(grid%on_this_pe) then
           grid%area_inv = 0.0;
           where (grid%area>0.0) grid%area_inv = 1.0/grid%area
        endif
     end if
  end do

  call mpp_clock_end(id_load_xgrid)

  grid1%area_inv = 0.0;
  where (grid1%area>0.0)
     grid1%area_inv = 1.0/grid1%area
  end where

  xmap%your1my2(xmap%me-xmap%root_pe) = .false. ! this is not necessarily true but keeps
  xmap%your2my1(xmap%me-xmap%root_pe) = .false. ! a PE from communicating with itself

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
     xmap%send_count_repro_tot = sum(xmap%send_count_repro)
     xmap%recv_count_repro_tot = sum(xmap%recv_count_repro)
  else
     xmap%send_count_repro_tot = 0
     xmap%recv_count_repro_tot = 0
  end if

  if (xgrid_log) then
    call mpp_open( unit, 'xgrid.out', action=MPP_OVERWR, threading=MPP_MULTI, &
         fileset=MPP_MULTI, nohdrs=.TRUE. )  

    write( unit,* )xmap%grids(:)%id, ' GRID: PE ', xmap%me, ' #XCELLS=', &
       xmap%grids(2:size(xmap%grids(:)))%size, ' #COMM. PARTNERS=', &
       count(xmap%your1my2), '/', count(xmap%your2my1), &
       pack((/(p+xmap%root_pe,p=0,xmap%npes-1)/), xmap%your1my2),  &
       '/', pack((/(p+xmap%root_pe,p=0,xmap%npes-1)/), xmap%your2my1)
    call close_file (unit)
  endif

  allocate( xmap%x1(1:sum(xmap%grids(2:size(xmap%grids(:)))%size)) )
  allocate( xmap%x2(1:sum(xmap%grids(2:size(xmap%grids(:)))%size)) )
  allocate( xmap%x1_put(1:sum(xmap%grids(2:size(xmap%grids(:)))%size)) )
  allocate( xmap%x2_get(1:sum(xmap%grids(2:size(xmap%grids(:)))%size)) )

  !--- The following will setup indx to be used in regen
  allocate(xmap%get1, xmap%put1)
  call mpp_clock_begin(id_set_comm)

  call set_comm_get1(xmap)

  call set_comm_put1(xmap)

  if(make_exchange_reproduce) then
    allocate(xmap%get1_repro)
    call set_comm_get1_repro(xmap)
  endif

  call mpp_clock_end(id_set_comm)

  call mpp_clock_begin(id_regen)
  call regen(xmap)
  call mpp_clock_end(id_regen)

  call mpp_clock_begin(id_conservation_check)

  xxx = conservation_check(grid1%area*0.0+1.0, grid1%id, xmap)
  write(out_unit,* )"Checked data is array of constant 1"
  write(out_unit,* )grid1%id,'(',xmap%grids(:)%id,')=', xxx 

  do g=2,size(xmap%grids(:))
     xxx = conservation_check(xmap%grids(g)%frac_area*0.0+1.0, xmap%grids(g)%id, xmap )
     write( out_unit,* )xmap%grids(g)%id,'(',xmap%grids(:)%id,')=', xxx 
  enddo
  ! create an random number 2d array
  if(grid1%id == "ATM") then
     allocate(check_data(size(grid1%area,1), size(grid1%area,2)))
     call random_number(check_data)

     !--- second order along both zonal and meridinal direction
     xxx = conservation_check(check_data, grid1%id, xmap,  remap_method = remapping_method )
     write( out_unit,* ) &
          "Checked data is array of random number between 0 and 1 using "//trim(interp_method)
     write( out_unit,* )grid1%id,'(',xmap%grids(:)%id,')=', xxx 

     deallocate(check_data)
     do g=2,size(xmap%grids(:))
        allocate(check_data_3d(size(xmap%grids(g)%frac_area,1),size(xmap%grids(g)%frac_area,2), &
             size(xmap%grids(g)%frac_area,3) )) 
        call random_number(check_data_3d)
        xxx = conservation_check(check_data_3d, xmap%grids(g)%id, xmap,  remap_method = remapping_method )
        write( out_unit,* )xmap%grids(g)%id,'(',xmap%grids(:)%id,')=', xxx
        deallocate( check_data_3d)
     end do
  endif
  call mpp_clock_end(id_conservation_check)

  call mpp_clock_end(id_setup_xmap)

end subroutine setup_xmap
! </SUBROUTINE>

!----------------------------------------------------------------------------
! currently we are assuming there is only one nest region
function get_nest_contact(mosaic_file, tile_nest_out, tile_parent_out, is_nest_out, &
                          ie_nest_out, js_nest_out, je_nest_out, is_parent_out, &
                          ie_parent_out, js_parent_out, je_parent_out)
character(len=*), intent(in) :: mosaic_file
integer,         intent(out) :: tile_nest_out, tile_parent_out
integer,         intent(out) :: is_nest_out, ie_nest_out
integer,         intent(out) :: js_nest_out, je_nest_out
integer,         intent(out) :: is_parent_out, ie_parent_out
integer,         intent(out) :: js_parent_out, je_parent_out
integer                      :: get_nest_contact
!--- local variables
integer                            :: ntiles, ncontacts, n, t1, t2
integer                            :: nx1_contact, ny1_contact
integer                            :: nx2_contact, ny2_contact
integer, allocatable, dimension(:) :: nx, ny
integer, allocatable, dimension(:) :: tile1, tile2
integer, allocatable, dimension(:) :: istart1, iend1, jstart1, jend1
integer, allocatable, dimension(:) :: istart2, iend2, jstart2, jend2

  tile_nest_out = 0; tile_parent_out = 0
  is_nest_out   = 0; ie_nest_out     = 0
  js_nest_out   = 0; je_nest_out     = 0
  is_parent_out = 0; ie_parent_out   = 0
  js_parent_out = 0; je_parent_out   = 0
  get_nest_contact = 0

  ! first read the contact information
  ntiles = get_mosaic_ntiles(mosaic_file)
  if( ntiles == 1 ) return
  allocate(nx(ntiles), ny(ntiles))  
  call get_mosaic_grid_sizes(mosaic_file, nx, ny)

  ncontacts = get_mosaic_ncontacts(mosaic_file)
  if(ncontacts == 0) return
  allocate(tile1(ncontacts), tile2(ncontacts))
  allocate(istart1(ncontacts), iend1(ncontacts))
  allocate(jstart1(ncontacts), jend1(ncontacts))
  allocate(istart2(ncontacts), iend2(ncontacts))
  allocate(jstart2(ncontacts), jend2(ncontacts))

  call get_mosaic_contact( mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1, &
                           istart2, iend2, jstart2, jend2)

  do n = 1, ncontacts
    if( tile1(n) == tile2(n) ) cycle ! same tile could not be nested 

    nx1_contact = iend1(n)-istart1(n)+1
    ny1_contact = jend1(n)-jstart1(n)+1
    nx2_contact = iend2(n)-istart2(n)+1
    ny2_contact = jend2(n)-jstart2(n)+1
    t1 = tile1(n);
    t2 = tile2(n);
    ! For nesting, the contact index of one tile must match its global domain 
    if( (nx(t1) .NE. nx1_contact .OR. ny(t1) .NE. ny1_contact ) .AND. &
        (nx(t2) .NE. nx2_contact .OR. ny(t2) .NE. ny2_contact ) ) cycle
    if(nx1_contact == nx2_contact .AND. ny1_contact == ny2_contact) then
      call error_mesg('xgrid_mod', 'There is no refinement for the overlapping region', FATAL)
    endif

    get_nest_contact = get_nest_contact + 1
    if(get_nest_contact>1) then
       call error_mesg('xgrid_mod', 'only support one nest region, contact developer' ,FATAL)
    endif
    if(nx2_contact*ny2_contact > nx1_contact*ny1_contact) then
      is_nest_out     = istart2(n);
      ie_nest_out     = iend2  (n);
      js_nest_out     = jstart2(n);
      je_nest_out     = jend2  (n);
      tile_nest_out   = tile2  (n);
      is_parent_out   = istart1(n);
      ie_parent_out   = iend1  (n);
      js_parent_out   = jstart1(n);
      je_parent_out   = jend1  (n);
      tile_parent_out = tile1  (n);
    else 
      is_nest_out     = istart1(n);
      ie_nest_out     = iend1  (n);
      js_nest_out     = jstart1(n);
      je_nest_out     = jend1  (n);
      tile_nest_out   = tile1  (n);
      is_parent_out   = istart2(n);
      ie_parent_out   = iend2  (n);
      js_parent_out   = jstart2(n);
      je_parent_out   = jend2  (n);
      tile_parent_out = tile2  (n);
    endif
  enddo

  deallocate(nx, ny, tile1, tile2)
  deallocate(istart1, iend1, jstart1, jend1)
  deallocate(istart2, iend2, jstart2, jend2)


  return
  
end function get_nest_contact

!#######################################################################
subroutine set_comm_get1_repro(xmap)
  type (xmap_type), intent(inout) :: xmap
  integer, dimension(xmap%npes) :: pe_ind, cnt
  integer, dimension(0:xmap%npes-1) :: send_ind, recv_ind, pl
  integer :: npes, nsend, nrecv, mypos
  integer :: m, p, pos, n, g, l
  type(comm_type), pointer, save :: comm => NULL()

  comm => xmap%get1_repro
  npes = xmap%npes

  nrecv = 0
  mypos = mpp_pe() - mpp_root_pe()
  do m=0,npes-1 
    p = mod(mypos+npes-m, npes)
    if( xmap%recv_count_repro(p) > 0 ) then
      nrecv = nrecv + 1
      pe_ind(nrecv) = p
    endif
  enddo

  comm%nrecv = nrecv
  if( nrecv > 0 ) then
    allocate(comm%recv(nrecv))
    pos = 0
    do n = 1, nrecv
      p = pe_ind(n)
      comm%recv(n)%count = xmap%recv_count_repro(p)
      comm%recv(n)%pe = p + xmap%root_pe
      comm%recv(n)%buffer_pos = pos 
      pos = pos + comm%recv(n)%count
    enddo
  endif


  ! send information
  nsend = 0
  mypos = mpp_pe() - mpp_root_pe()
  do m=0,xmap%npes-1 
    p = mod(mypos+m, npes)
    if( xmap%send_count_repro(p) > 0 ) then
      nsend = nsend + 1
      pe_ind(nsend) = p
      send_ind(p) = nsend
    endif
  enddo

  comm%nsend = nsend
  if( nsend > 0 ) then
     allocate(comm%send(nsend))
     pos = 0
     cnt(:) = 0
     do n = 1, nsend
        p = pe_ind(n)
        comm%send(n)%count = xmap%send_count_repro(p)
        comm%send(n)%pe = p + xmap%root_pe 
        comm%send(n)%buffer_pos = pos
        pos = pos + comm%send(n)%count
        allocate(comm%send(n)%i(comm%send(n)%count))
        allocate(comm%send(n)%j(comm%send(n)%count))
        allocate(comm%send(n)%g(comm%send(n)%count))
        allocate(comm%send(n)%xLoc(comm%send(n)%count))
     enddo

     do g=2,size(xmap%grids(:))
        do l=1,xmap%grids(g)%size ! index into this side 2 grid's patterns
           p = xmap%grids(g)%x(l)%pe-xmap%root_pe 
           n = send_ind(p)
           cnt(n) = cnt(n) + 1
           pos = cnt(n)
           comm%send(n)%i(pos) = xmap%grids(g)%x(l)%i2
           comm%send(n)%j(pos) = xmap%grids(g)%x(l)%j2
           comm%send(n)%g(pos) = g
        enddo
     enddo
     !--- make sure the count is correct
     do n = 1, nsend
        if( comm%send(n)%count .NE. cnt(n) ) call error_mesg('xgrid_mod', &
             'comm%send(n)%count .NE. cnt(n)', FATAL)
     enddo
   endif

   !--- set up the recv_pos for unpack the data.
   pl(:) = 1
   do g=2,size(xmap%grids(:))
      do l=1,xmap%grids(g)%size_repro ! index into side1 grid's patterns
         p = xmap%grids(g)%x_repro(l)%pe-xmap%root_pe
         xmap%grids(g)%x_repro(l)%recv_pos = pl(p)
         pl(p) = pl(p) + 1
      end do
   end do



end subroutine set_comm_get1_repro

!#######################################################################
subroutine set_comm_get1(xmap)
  type (xmap_type), intent(inout) :: xmap
  type (grid_type), pointer, save :: grid1 =>NULL()
  integer, allocatable :: send_size(:)
  integer, allocatable :: recv_size(:)
  integer              :: max_size, g, npes, l, ll, nset, m
  integer              :: i1, j1, tile1, p, n, pos, buffer_pos, mypos
  integer              :: nsend, nrecv, rbuf_size, sbuf_size, msgsize
  logical              :: found
  real,    allocatable :: recv_buf(:), send_buf(:)
  real,    allocatable :: diarray(:), djarray(:)
  integer, allocatable :: iarray(:), jarray(:), tarray(:)
  integer, allocatable :: pos_x(:), pelist(:), size_pe(:), pe_side1(:)
  integer              :: recv_buffer_pos(0:xmap%npes)
  integer              :: send_buffer_pos(0:xmap%npes)
  type(comm_type), pointer, save :: comm => NULL()

  max_size = 0
  do g=2,size(xmap%grids(:))
    max_size = max_size + xmap%grids(g)%size
  enddo
  comm => xmap%get1
  grid1 => xmap%grids(1)
  comm%nsend = 0
  comm%nrecv = 0
  npes = xmap%npes

  allocate(pelist(0:npes-1))
  call mpp_get_current_pelist(pelist)
  allocate(send_size(0:npes-1))
  allocate(recv_size(0:npes-1))
  allocate(size_pe(0:npes-1))
  allocate(pos_x(0:npes-1))
  size_pe = 0
  send_size = 0
  recv_size = 0

  if(max_size > 0) then
     allocate(pe_side1(max_size))
     allocate(xmap%ind_get1(max_size))

     !--- find the recv_indx
     ll = 0 
     do g=2,size(xmap%grids(:))
        do l=1,xmap%grids(g)%size
           i1 = xmap%grids(g)%x(l)%i1
           j1 = xmap%grids(g)%x(l)%j1
           tile1 = xmap%grids(g)%x(l)%tile
           do p=0,npes-1
              if(grid1%tile(p) == tile1) then
                 if(in_box(i1, j1, grid1%is(p), grid1%ie(p), grid1%js(p), grid1%je(p))) then
                    size_pe(p) = size_pe(p) + 1     
                    exit
                 endif
              endif
           enddo
           if( p == npes ) then
              call error_mesg('xgrid_mod', 'tile is not in grid1%tile(:)', FATAL)
           endif
           ll = ll + 1
           pe_side1(ll) = p
        enddo
     enddo

     pos_x = 0
     do p = 1, npes-1
        pos_x(p) = pos_x(p-1) + size_pe(p-1)
     enddo

     !---find the send size for get_1_from_xgrid
     allocate(iarray(max_size))
     allocate(jarray(max_size))
     allocate(tarray(max_size))
     if(monotonic_exchange) then
        allocate(diarray(max_size))
        allocate(djarray(max_size))        
     endif

     ll = 0

     do g=2,size(xmap%grids(:))
        do l=1,xmap%grids(g)%size
           i1 = xmap%grids(g)%x(l)%i1
           j1 = xmap%grids(g)%x(l)%j1
           tile1 = xmap%grids(g)%x(l)%tile
           ll = ll + 1
           p = pe_side1(ll)

           found = .false.
           if(send_size(p) > 0) then
              if( i1 == iarray(pos_x(p)+send_size(p)) .AND. j1 == jarray(pos_x(p)+send_size(p)) &
                   .AND. tile1 == tarray(pos_x(p)+send_size(p))) then
                 found = .true.
                 n = send_size(p)
              else
                 !---may need to replace with a fast search algorithm
                 do n = 1, send_size(p)
                    if(i1 == iarray(pos_x(p)+n) .AND. j1 == jarray(pos_x(p)+n) .AND. tile1 == tarray(pos_x(p)+n)) then
                       found = .true.
                       exit
                    endif
                 enddo
              endif
           endif
           if( (.NOT. found) .OR. monotonic_exchange ) then
              send_size(p) = send_size(p)+1
              pos = pos_x(p)+send_size(p)
              iarray(pos) = i1
              jarray(pos) = j1
              tarray(pos) = tile1
              if(monotonic_exchange) then
                 diarray(pos) = xmap%grids(g)%x(l)%di
                 djarray(pos) = xmap%grids(g)%x(l)%dj
              endif
              n = send_size(p)
           endif
           xmap%ind_get1(ll) = n
        enddo
     enddo

     pos_x = 0
     do p = 1, npes-1
        pos_x(p) = pos_x(p-1) + send_size(p-1)
     enddo

     ll = 0
     do g=2,size(xmap%grids(:))
        do l=1,xmap%grids(g)%size
           ll = ll + 1
           p = pe_side1(ll)  
           xmap%ind_get1(ll) = pos_x(p) + xmap%ind_get1(ll)
        enddo
     enddo
  endif

  mypos = mpp_pe()-mpp_root_pe()
 
  ! send/recv for get_1_from_xgrid_recv
  recv_size(:) = xmap%your2my1_size(:)
  nsend = count( send_size> 0)  
  comm%nsend = nsend
  if(nsend>0) then
     allocate(comm%send(nsend))
     comm%send(:)%count = 0
  endif

  pos = 0
  do p = 0, npes-1
     send_buffer_pos(p) = pos
     pos = pos + send_size(p)
  enddo

  pos = 0
  comm%sendsize = 0
  do n = 0, npes-1
     p = mod(mypos+n, npes)
     if(send_size(p)>0) then
        pos = pos + 1
        allocate(comm%send(pos)%i(send_size(p)))
        comm%send(pos)%buffer_pos = send_buffer_pos(p)
        comm%send(pos)%count = send_size(p)
        comm%send(pos)%pe = pelist(p)
        comm%sendsize = comm%sendsize + send_size(p)
     endif
  enddo

  nset = 3
  if(monotonic_exchange) nset = 5
  rbuf_size = sum(recv_size)*nset
  sbuf_size = sum(send_size)*nset
  if(rbuf_size>0) allocate(recv_buf(rbuf_size))
  if(sbuf_size>0) allocate(send_buf(sbuf_size))

  pos = 0
  do n = 0, npes-1
     p = mod(mypos+npes-n, npes)
     if(recv_size(p) ==0) cycle
     msgsize = recv_size(p)*nset
     call mpp_recv(recv_buf(pos+1), glen=msgsize, from_pe=pelist(p), block=.false., tag=COMM_TAG_4)
     pos = pos + msgsize
  enddo

  pos_x = 0
  do p = 1, npes-1
     pos_x(p) = pos_x(p-1) + size_pe(p-1)
  enddo
  ll = 0
  pos = 0
  do n = 0, npes-1
     p = mod(mypos+n, npes)
     do l = 1, send_size(p)
        send_buf(pos+1) = iarray(pos_x(p)+l)
        send_buf(pos+2) = jarray(pos_x(p)+l)
        send_buf(pos+3) = tarray(pos_x(p)+l)
        if(monotonic_exchange) then
           send_buf(pos+4) = diarray(pos_x(p)+l)
           send_buf(pos+5) = djarray(pos_x(p)+l)
        endif
        pos = pos + nset
     enddo
  enddo

  pos = 0
  do n = 0, npes-1
     p = mod(mypos+n, npes)
     if(send_size(p) ==0) cycle
     msgsize = send_size(p)*nset
     call mpp_send(send_buf(pos+1), plen=msgsize, to_pe=pelist(p), tag=COMM_TAG_4 )
     pos = pos + msgsize
  enddo

  call mpp_sync_self(check=EVENT_RECV)
  nrecv = count(recv_size>0)
  comm%nrecv = nrecv
  comm%recvsize = 0

  if(nrecv >0) then
     allocate(comm%recv(nrecv))
     comm%recv(:)%count = 0
     !--- set up the buffer pos for each receiving
     buffer_pos = 0
     do p = 0, npes-1
        recv_buffer_pos(p) = buffer_pos
        buffer_pos = buffer_pos +  recv_size(p)
     enddo
     pos = 0
     buffer_pos = 0        
     do m=0,npes-1
        p = mod(mypos+npes-m, npes)
        if(recv_size(p)>0) then
           pos = pos + 1
           allocate(comm%recv(pos)%i(recv_size(p)))
           allocate(comm%recv(pos)%j(recv_size(p)))
           allocate(comm%recv(pos)%tile(recv_size(p)))
           comm%recv(pos)%buffer_pos = recv_buffer_pos(p)
           comm%recv(pos)%pe = pelist(p)
           comm%recv(pos)%count = recv_size(p)
           comm%recvsize = comm%recvsize + recv_size(p)
           if(monotonic_exchange) then
              allocate(comm%recv(pos)%di(recv_size(p)))
              allocate(comm%recv(pos)%dj(recv_size(p)))
           endif
           do n = 1, recv_size(p)
              comm%recv(pos)%i(n) = recv_buf(buffer_pos+1) - grid1%is_me + 1
              comm%recv(pos)%j(n) = recv_buf(buffer_pos+2) - grid1%js_me + 1
              comm%recv(pos)%tile(n) = recv_buf(buffer_pos+3) 
              if(monotonic_exchange) then
                 comm%recv(pos)%di(n) = recv_buf(buffer_pos+4) 
                 comm%recv(pos)%dj(n) = recv_buf(buffer_pos+5) 
              endif
              buffer_pos = buffer_pos + nset
           enddo
        endif
     enddo
     allocate(comm%unpack_ind(nrecv))
     pos = 0
     do p = 0, npes-1
        if(recv_size(p)>0) then
           pos = pos + 1
           do m = 1, nrecv
              if(comm%recv(m)%pe == pelist(p)) then
                 comm%unpack_ind(pos) = m
                 exit
              endif
           enddo
        endif
     enddo
  endif
  call mpp_sync_self()

  if(allocated(send_buf) ) deallocate(send_buf)
  if(allocated(recv_buf) ) deallocate(recv_buf)
  if(allocated(pelist)   ) deallocate(pelist)
  if(allocated(pos_x)    ) deallocate(pos_x)
  if(allocated(pelist)   ) deallocate(pelist)
  if(allocated(iarray)   ) deallocate(iarray)
  if(allocated(jarray)   ) deallocate(jarray)
  if(allocated(tarray)   ) deallocate(tarray)
  if(allocated(size_pe)  ) deallocate(size_pe)

end subroutine set_comm_get1

!###############################################################################
subroutine set_comm_put1(xmap)
  type (xmap_type), intent(inout) :: xmap
  type (grid_type), pointer, save :: grid1 =>NULL()
  integer, allocatable :: send_size(:)
  integer, allocatable :: recv_size(:)
  integer              :: max_size, g, npes, l, ll, m, mypos
  integer              :: i1, j1, tile1, p, n, pos, buffer_pos
  integer              :: nsend, nrecv, msgsize, nset, rbuf_size, sbuf_size
  logical              :: found
  real,    allocatable :: recv_buf(:), send_buf(:)
  real,    allocatable :: diarray(:), djarray(:)
  integer, allocatable :: iarray(:), jarray(:), tarray(:)
  integer, allocatable :: pos_x(:), pelist(:), size_pe(:), pe_put1(:)
  integer              :: root_pe, recvsize, sendsize
  integer              :: recv_buffer_pos(0:xmap%npes)
  type(comm_type), pointer, save :: comm => NULL()


  comm => xmap%put1
  if(nnest == 0 .OR. xmap%grids(1)%id .NE. 'ATM' ) then
     comm%nsend    = xmap%get1%nrecv
     comm%nrecv    = xmap%get1%nsend   
     comm%sendsize = xmap%get1%recvsize
     comm%recvsize = xmap%get1%sendsize    
     comm%send     => xmap%get1%recv
     comm%recv     => xmap%get1%send
     xmap%ind_put1 => xmap%ind_get1 
    return
  endif

  max_size = 0
  do g=2,size(xmap%grids(:))
     max_size = max_size + xmap%grids(g)%size
  enddo
  grid1 => xmap%grids(1)
  comm%nsend = 0
  comm%nrecv = 0
  npes = xmap%npes
  allocate(pelist(0:npes-1))
  call mpp_get_current_pelist(pelist)
  allocate(send_size(0:npes-1))
  allocate(recv_size(0:npes-1))
  allocate(size_pe(0:npes-1))
  allocate(pos_x(0:npes-1))
  size_pe = 0
  send_size = 0
  recv_size = 0

  if(max_size > 0) then
     allocate(pe_put1(max_size))
     allocate(xmap%ind_put1(max_size))

     !--- find the recv_indx
     ll = 0 
     do g=2,size(xmap%grids(:))
        do l=1,xmap%grids(g)%size
           i1 = xmap%grids(g)%x(l)%i1
           j1 = xmap%grids(g)%x(l)%j1
           tile1 = xmap%grids(g)%x(l)%tile
           do p=0,npes-1
              if(grid1%tile(p) == tile1) then
                 if(in_box(i1, j1, grid1%is(p), grid1%ie(p), grid1%js(p), grid1%je(p))) then
                    size_pe(p) = size_pe(p) + 1     
                    exit
                 endif
              endif
           enddo
           ll = ll + 1
           pe_put1(ll) = p
        enddo
     enddo

     pos_x = 0
     do p = 1, npes-1
        pos_x(p) = pos_x(p-1) + size_pe(p-1)
     enddo

     !---find the send size for get_1_from_xgrid
     allocate(iarray(max_size))
     allocate(jarray(max_size))
     allocate(tarray(max_size))
     if(monotonic_exchange) then
        allocate(diarray(max_size))
        allocate(djarray(max_size))        
     endif

     ll = 0

     do g=2,size(xmap%grids(:))
        do l=1,xmap%grids(g)%size
           i1 = xmap%grids(g)%x(l)%i1
           j1 = xmap%grids(g)%x(l)%j1
           tile1 = xmap%grids(g)%x(l)%tile
           ll = ll + 1
           p = pe_put1(ll)

           found = .false.
           if(send_size(p) > 0) then
              if( i1 == iarray(pos_x(p)+send_size(p)) .AND. j1 == jarray(pos_x(p)+send_size(p)) &
                   .AND. tile1 == tarray(pos_x(p)+send_size(p))) then
                 found = .true.
                 n = send_size(p)
              else
                 !---may need to replace with a fast search algorithm
                 do n = 1, send_size(p)
                    if(i1 == iarray(pos_x(p)+n) .AND. j1 == jarray(pos_x(p)+n) .AND. tile1 == tarray(pos_x(p)+n)) then
                       found = .true.
                       exit
                    endif
                 enddo
              endif
           endif
           if( (.NOT. found) .OR. monotonic_exchange ) then
              send_size(p) = send_size(p)+1
              pos = pos_x(p)+send_size(p)
              iarray(pos) = i1
              jarray(pos) = j1
              tarray(pos) = tile1
              if(monotonic_exchange) then
                 diarray(pos) = xmap%grids(g)%x(l)%di
                 djarray(pos) = xmap%grids(g)%x(l)%dj
              endif
              n = send_size(p)
           endif
           xmap%ind_put1(ll) = n
        enddo
     enddo

     pos_x = 0
     do p = 1, npes-1
        pos_x(p) = pos_x(p-1) + send_size(p-1)
     enddo

     ll = 0
     do g=2,size(xmap%grids(:))
        do l=1,xmap%grids(g)%size
           i1 = xmap%grids(g)%x(l)%i1
           j1 = xmap%grids(g)%x(l)%j1
           tile1 = xmap%grids(g)%x(l)%tile
           ll = ll + 1
           p = pe_put1(ll)  
           xmap%ind_put1(ll) = pos_x(p) + xmap%ind_put1(ll)
        enddo
     enddo
  endif

  mypos = mpp_pe()-mpp_root_pe()
  do n = 0, npes-1
     p = mod(mypos+npes-n, npes)
     call mpp_recv(recv_size(p), glen=1, from_pe=pelist(p), block=.false., tag=COMM_TAG_5)
  enddo

  !--- send data
  do n = 0, npes-1
     p = mod(mypos+n, npes)
     call mpp_send(send_size(p), plen=1, to_pe=pelist(p), tag=COMM_TAG_5)
  enddo

  call mpp_sync_self(check=EVENT_RECV)
  call mpp_sync_self()

  !--- recv for put_1_to_xgrid
  nrecv = count( send_size> 0)  
  comm%nrecv = nrecv
  if(nrecv>0) then
     allocate(comm%recv(nrecv))
     comm%recv(:)%count = 0
  endif
  pos = 0
  comm%recvsize = 0
  do p = 0, npes-1
     recv_buffer_pos(p) = pos
     pos = pos + send_size(p)
  enddo

  pos = 0
  do n = 0, npes-1
     p = mod(mypos+npes-n, npes)
     if(send_size(p)>0) then
        pos = pos + 1
        allocate(comm%recv(pos)%i(send_size(p)))
        comm%recv(pos)%buffer_pos = recv_buffer_pos(p)
        comm%recv(pos)%count = send_size(p)
        comm%recv(pos)%pe = pelist(p)
        comm%recvsize = comm%recvsize + send_size(p)
     endif
  enddo

  nset = 3
  if(monotonic_exchange) nset = 5
  rbuf_size = sum(recv_size)*nset
  sbuf_size = sum(send_size)*nset
  if(rbuf_size>0) allocate(recv_buf(rbuf_size))
  if(sbuf_size>0) allocate(send_buf(sbuf_size))

  pos = 0
  do n = 0, npes-1
     p = mod(mypos+npes-n, npes)
     if(recv_size(p) ==0) cycle
     msgsize = recv_size(p)*nset
     call mpp_recv(recv_buf(pos+1), glen=msgsize, from_pe=pelist(p), block=.false., tag=COMM_TAG_6)
     pos = pos + msgsize
  enddo

  pos_x = 0
  do p = 1, npes-1
     pos_x(p) = pos_x(p-1) + size_pe(p-1)
  enddo
  ll = 0
  pos = 0
  do n = 0, npes-1
     p = mod(mypos+n, npes)
     do l = 1, send_size(p)
        send_buf(pos+1) = iarray(pos_x(p)+l)
        send_buf(pos+2) = jarray(pos_x(p)+l)
        send_buf(pos+3) = tarray(pos_x(p)+l)
        if(monotonic_exchange) then
           send_buf(pos+4) = diarray(pos_x(p)+l)
           send_buf(pos+5) = djarray(pos_x(p)+l)
        endif
        pos = pos + nset
     enddo
  enddo

  pos = 0
  do n = 0, npes-1
     p = mod(mypos+n, npes)
     if(send_size(p) ==0) cycle
     msgsize = send_size(p)*nset
     call mpp_send(send_buf(pos+1), plen=msgsize, to_pe=pelist(p), tag=COMM_TAG_6 )
     pos = pos + msgsize
  enddo

  call mpp_sync_self(check=EVENT_RECV)
  nsend = count(recv_size>0)
  comm%nsend = nsend
  comm%sendsize = 0

  if(nsend >0) then
     allocate(comm%send(nsend))
     comm%send(:)%count = 0
     pos = 0
     buffer_pos = 0
     do m=0,npes-1
        p = mod(mypos+npes-m, npes)
        if(recv_size(p)>0) then
           pos = pos + 1
           allocate(comm%send(pos)%i(recv_size(p)))
           allocate(comm%send(pos)%j(recv_size(p)))
           allocate(comm%send(pos)%tile(recv_size(p)))
           comm%send(pos)%pe = pelist(p)
           comm%send(pos)%count = recv_size(p)
           comm%sendsize = comm%sendsize + recv_size(p)
           if(monotonic_exchange) then
              allocate(comm%send(pos)%di(recv_size(p)))
              allocate(comm%send(pos)%dj(recv_size(p)))
           endif
           do n = 1, recv_size(p)
              comm%send(pos)%i(n) = recv_buf(buffer_pos+1) - grid1%is_me + 1
              comm%send(pos)%j(n) = recv_buf(buffer_pos+2) - grid1%js_me + 1
              comm%send(pos)%tile(n) = recv_buf(buffer_pos+3) 
              if(monotonic_exchange) then
                 comm%send(pos)%di(n) = recv_buf(buffer_pos+4) 
                 comm%send(pos)%dj(n) = recv_buf(buffer_pos+5) 
              endif
              buffer_pos = buffer_pos + nset
           enddo
        endif
     enddo
  endif

  call mpp_sync_self()
  if(allocated(send_buf) ) deallocate(send_buf)
  if(allocated(recv_buf) ) deallocate(recv_buf)
  if(allocated(pelist)   ) deallocate(pelist)
  if(allocated(pos_x)    ) deallocate(pos_x)
  if(allocated(pelist)   ) deallocate(pelist)
  if(allocated(iarray)   ) deallocate(iarray)
  if(allocated(jarray)   ) deallocate(jarray)
  if(allocated(tarray)   ) deallocate(tarray)
  if(allocated(size_pe)  ) deallocate(size_pe)

end subroutine set_comm_put1


!###############################################################################
subroutine regen(xmap)
type (xmap_type), intent(inout) :: xmap

  integer              :: g, l, k, max_size
  integer              :: i1, j1, i2, j2, p
  integer              :: tile1
  integer              :: ll
  logical              :: overlap_with_nest
  integer              :: cnt(xmap%get1%nsend)
  integer              :: i,j,n,xloc,pos,nsend,m,npes, mypos
  integer              :: send_ind(0:xmap%npes-1)

  max_size = 0

  do g=2,size(xmap%grids(:))
    max_size = max_size + xmap%grids(g)%size * xmap%grids(g)%km
  end do

  if (max_size>size(xmap%x1(:))) then
    deallocate(xmap%x1)
    deallocate(xmap%x2)
    allocate( xmap%x1(1:max_size) )
    allocate( xmap%x2(1:max_size) )
  endif


  do g=2,size(xmap%grids(:))
    xmap%grids(g)%first = 1
    xmap%grids(g)%last  = 0
  end do  

  xmap%size = 0
  ll = 0
  do g=2,size(xmap%grids(:))
     xmap%grids(g)%first = xmap%size + 1;
     
     do l=1,xmap%grids(g)%size
        i1 = xmap%grids(g)%x(l)%i1
        j1 = xmap%grids(g)%x(l)%j1
        i2 = xmap%grids(g)%x(l)%i2
        j2 = xmap%grids(g)%x(l)%j2
        tile1 = xmap%grids(g)%x(l)%tile
        ll = ll + 1
        do k=1,xmap%grids(g)%km
           if (xmap%grids(g)%frac_area(i2,j2,k)/=0.0) then
              xmap%size = xmap%size+1
              xmap%x1(xmap%size)%pos  = xmap%ind_get1(ll)
              xmap%x1(xmap%size)%i    = xmap%grids(g)%x(l)%i1
              xmap%x1(xmap%size)%j    = xmap%grids(g)%x(l)%j1
              xmap%x1(xmap%size)%tile = xmap%grids(g)%x(l)%tile
              xmap%x1(xmap%size)%area = xmap%grids(g)%x(l)%area &
                   *xmap%grids(g)%frac_area(i2,j2,k)
              xmap%x1(xmap%size)%di   = xmap%grids(g)%x(l)%di 
              xmap%x1(xmap%size)%dj   = xmap%grids(g)%x(l)%dj 
              xmap%x2(xmap%size)%i    = xmap%grids(g)%x(l)%i2
              xmap%x2(xmap%size)%j    = xmap%grids(g)%x(l)%j2
              xmap%x2(xmap%size)%k    = k
              xmap%x2(xmap%size)%area = xmap%grids(g)%x(l)%area * xmap%grids(g)%x(l)%scale 
           end if
        end do
     end do
     xmap%grids(g)%last = xmap%size
  end do


  if (max_size>size(xmap%x1_put(:))) then
    deallocate(xmap%x1_put)
    allocate( xmap%x1_put(1:max_size) )
  endif
  if (max_size>size(xmap%x2_get(:))) then
    deallocate(xmap%x2_get)
    allocate( xmap%x2_get(1:max_size) )
  endif

  do g=2,size(xmap%grids(:))
    xmap%grids(g)%first_get = 1
    xmap%grids(g)%last_get  = 0
  end do  

  xmap%size_put1 = 0
  xmap%size_get2 = 0
  ll = 0
  do g=2,size(xmap%grids(:))
     xmap%grids(g)%first_get = xmap%size_get2 + 1;
     
     do l=1,xmap%grids(g)%size
        i1 = xmap%grids(g)%x(l)%i1
        j1 = xmap%grids(g)%x(l)%j1
        i2 = xmap%grids(g)%x(l)%i2
        j2 = xmap%grids(g)%x(l)%j2
        tile1 = xmap%grids(g)%x(l)%tile
        ll = ll + 1
        overlap_with_nest = .false.
        if(  xmap%grids(1)%id == "ATM" .AND. tile1 == tile_parent .AND. &
             in_box(i1, j1, is_parent, ie_parent, js_parent, je_parent) ) overlap_with_nest = .true.
        do k=1,xmap%grids(g)%km
           if (xmap%grids(g)%frac_area(i2,j2,k)/=0.0) then
              xmap%size_put1 = xmap%size_put1+1
              xmap%x1_put(xmap%size_put1)%pos  = xmap%ind_put1(ll)
              xmap%x1_put(xmap%size_put1)%i    = xmap%grids(g)%x(l)%i1
              xmap%x1_put(xmap%size_put1)%j    = xmap%grids(g)%x(l)%j1
              xmap%x1_put(xmap%size_put1)%tile = xmap%grids(g)%x(l)%tile
              xmap%x1_put(xmap%size_put1)%area = xmap%grids(g)%x(l)%area &
                   *xmap%grids(g)%frac_area(i2,j2,k)
              xmap%x1_put(xmap%size_put1)%di   = xmap%grids(g)%x(l)%di 
              xmap%x1_put(xmap%size_put1)%dj   = xmap%grids(g)%x(l)%dj 
              if( .not. overlap_with_nest) then
                 xmap%size_get2 = xmap%size_get2+1
                 xmap%x2_get(xmap%size_get2)%i    = xmap%grids(g)%x(l)%i2
                 xmap%x2_get(xmap%size_get2)%j    = xmap%grids(g)%x(l)%j2
                 xmap%x2_get(xmap%size_get2)%k    = k
                 xmap%x2_get(xmap%size_get2)%area = xmap%grids(g)%x(l)%area * xmap%grids(g)%x(l)%scale 
                 xmap%x2_get(xmap%size_get2)%pos  = xmap%size_put1
              endif
           end if
        end do
     end do
     xmap%grids(g)%last_get = xmap%size_get2
  end do

  !---set up information for get_1_from_xgrid_repro
  if (make_exchange_reproduce) then
  if (xmap%get1_repro%nsend > 0) then
     xloc = 0
     nsend = 0
     npes = xmap%npes
     mypos = mpp_pe() - mpp_root_pe()
     cnt(:) = 0
     do m=0,npes-1 
        p = mod(mypos+m, npes)
        if( xmap%send_count_repro(p) > 0 ) then
          nsend = nsend + 1
          send_ind(p) = nsend
        endif
     enddo
     do g=2,size(xmap%grids(:))
        do l=1,xmap%grids(g)%size ! index into this side 2 grid's patterns
           i = xmap%grids(g)%x(l)%i2
           j = xmap%grids(g)%x(l)%j2
           p = xmap%grids(g)%x(l)%pe-xmap%root_pe 
           n = send_ind(p) 
           cnt(n) = cnt(n) + 1
           pos = cnt(n) 
           xmap%get1_repro%send(n)%xLoc(pos) = xloc
           xloc = xloc + count(xmap%grids(g)%frac_area(i,j,:)/=0.0)
        enddo
     enddo
  endif
  endif

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

subroutine put_side1_to_xgrid(d, grid_id, x, xmap, remap_method, complete)
  real, dimension(:,:), intent(in   )    :: d
  character(len=3),     intent(in   )    :: grid_id
  real, dimension(:),   intent(inout)    :: x
  type (xmap_type),     intent(inout)    :: xmap
  integer, intent(in), optional          :: remap_method
  logical, intent(in), optional          :: complete

  logical                                         :: is_complete, set_mismatch
  integer                                         :: g, method
  character(len=2)                                :: text
  integer,                                   save :: isize=0
  integer,                                   save :: jsize=0
  integer,                                   save :: lsize=0
  integer,                                   save :: xsize=0
  integer,                                   save :: method_saved=0
  character(len=3),                          save :: grid_id_saved=""
  integer(LONG_KIND), dimension(MAX_FIELDS), save :: d_addrs=-9999
  integer(LONG_KIND), dimension(MAX_FIELDS), save :: x_addrs=-9999
  
  if (grid_id==xmap%grids(1)%id) then
     method = FIRST_ORDER      ! default
     if(present(remap_method)) method = remap_method
     is_complete = .true.
     if(present(complete)) is_complete=complete
     lsize = lsize + 1
     if( lsize > MAX_FIELDS ) then
        write( text,'(i2)' ) MAX_FIELDS
        call error_mesg ('xgrid_mod',  'MAX_FIELDS='//trim(text)//' exceeded for group put_side1_to_xgrid', FATAL)
     endif
     d_addrs(lsize) = LOC(d)
     x_addrs(lsize) = LOC(x)  
   
     if(lsize == 1) then
        isize = size(d,1)
        jsize = size(d,2)
        xsize = size(x(:))
        method_saved = method
        grid_id_saved = grid_id
     else
        set_mismatch = .false.
        set_mismatch = set_mismatch .OR. (isize /= size(d,1))
        set_mismatch = set_mismatch .OR. (jsize /= size(d,2))
        set_mismatch = set_mismatch .OR. (xsize /= size(x(:)))
        set_mismatch = set_mismatch .OR. (method_saved /= method)
        set_mismatch = set_mismatch .OR. (grid_id_saved /= grid_id)
        if(set_mismatch)then
           write( text,'(i2)' ) lsize
           call error_mesg ('xgrid_mod', 'Incompatible field at count '//text//' for group put_side1_to_xgrid', FATAL )
        endif
     endif

     if(is_complete) then
        !--- when exchange_monotonic is true and the side 1 ia atm, will always use monotonic second order conservative.
        if(monotonic_exchange .AND. grid_id == 'ATM') then
           call put_1_to_xgrid_order_2(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)
        else if(method == FIRST_ORDER) then
           call put_1_to_xgrid_order_1(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)
        else 
           if(grid_id .NE. 'ATM') call error_mesg ('xgrid_mod',  &
                "second order put_to_xgrid should only be applied to 'ATM' model, "//&
                "contact developer", FATAL)
           call put_1_to_xgrid_order_2(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)
        endif

        d_addrs = -9999
        x_addrs = -9999
        isize   = 0
        jsize   = 0
        xsize   = 0
        lsize   = 0
        method_saved = 0
        grid_id_saved = ""
     endif
     return
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

subroutine get_side1_from_xgrid(d, grid_id, x, xmap, complete)
  real, dimension(:,:), intent(  out) :: d
  character(len=3),     intent(in   ) :: grid_id
  real, dimension(:),   intent(in   ) :: x
  type (xmap_type),     intent(inout) :: xmap
  logical, intent(in), optional     :: complete

  logical                                         :: is_complete, set_mismatch
  integer                                         :: g
  character(len=2)                                :: text
  integer,                                   save :: isize=0
  integer,                                   save :: jsize=0
  integer,                                   save :: lsize=0
  integer,                                   save :: xsize=0
  character(len=3),                          save :: grid_id_saved=""
  integer(LONG_KIND), dimension(MAX_FIELDS), save :: d_addrs=-9999
  integer(LONG_KIND), dimension(MAX_FIELDS), save :: x_addrs=-9999

  if (grid_id==xmap%grids(1)%id) then
     is_complete = .true.
     if(present(complete)) is_complete=complete
     lsize = lsize + 1
     if( lsize > MAX_FIELDS ) then
        write( text,'(i2)' ) MAX_FIELDS
        call error_mesg ('xgrid_mod',  'MAX_FIELDS='//trim(text)//' exceeded for group get_side1_from_xgrid', FATAL)
     endif
     d_addrs(lsize) = LOC(d)
     x_addrs(lsize) = LOC(x)  

     if(lsize == 1) then
        isize = size(d,1)
        jsize = size(d,2)
        xsize = size(x(:))
        grid_id_saved = grid_id
     else
        set_mismatch = .false.
        set_mismatch = set_mismatch .OR. (isize /= size(d,1))
        set_mismatch = set_mismatch .OR. (jsize /= size(d,2))
        set_mismatch = set_mismatch .OR. (xsize /= size(x(:)))
        set_mismatch = set_mismatch .OR. (grid_id_saved /= grid_id)
        if(set_mismatch)then
           write( text,'(i2)' ) lsize
           call error_mesg ('xgrid_mod', 'Incompatible field at count '//text//' for group get_side1_from_xgrid', FATAL )
        endif
     endif

     if(is_complete) then
        if (make_exchange_reproduce) then
           call get_1_from_xgrid_repro(d_addrs, x_addrs, xmap, xsize, lsize)
        else
           call get_1_from_xgrid(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)
        end if
        d_addrs = -9999
        x_addrs = -9999
        isize   = 0
        jsize   = 0
        xsize   = 0
        lsize   = 0
        grid_id_saved = ""
     endif
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
logical, dimension(:), intent(out) :: some_arr

  integer :: g

  if (.not.present(grid_id)) then

    if(xmap%size > 0) then
       some_arr = .true.
    else 
       some_arr = .false.
    end if
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
  call mpp_clock_begin(id_put_2_to_xgrid)

  do l=grid%first,grid%last
    x(l) = d(xmap%x2(l)%i,xmap%x2(l)%j,xmap%x2(l)%k)
  end do

  call mpp_clock_end(id_put_2_to_xgrid)
end subroutine put_2_to_xgrid

!#######################################################################

subroutine get_2_from_xgrid(d, grid, x, xmap)
type (grid_type),                                intent(in ) :: grid
real, dimension(grid%is_me:grid%ie_me, &
                grid%js_me:grid%je_me, grid%km), intent(out) :: d
real, dimension(:),     intent(in   ) :: x
type (xmap_type),       intent(in   ) :: xmap

  integer                 :: l, k

  call mpp_clock_begin(id_get_2_from_xgrid)

  d = 0.0
  do l=grid%first_get,grid%last_get
    d(xmap%x2_get(l)%i,xmap%x2_get(l)%j,xmap%x2_get(l)%k) = &
            d(xmap%x2_get(l)%i,xmap%x2_get(l)%j,xmap%x2_get(l)%k) + xmap%x2_get(l)%area*x(xmap%x2_get(l)%pos)
  end do
  !
  !  normalize with side 2 grid cell areas
  !
  do k=1,size(d,3)
    d(:,:,k) = d(:,:,k) * grid%area_inv
  end do

  call mpp_clock_end(id_get_2_from_xgrid)

end subroutine get_2_from_xgrid

!#######################################################################

subroutine put_1_to_xgrid_order_1(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)
  integer(LONG_KIND), dimension(:), intent(in) :: d_addrs
  integer(LONG_KIND), dimension(:), intent(in) :: x_addrs
  type (xmap_type),              intent(inout) :: xmap
  integer,                          intent(in) :: isize, jsize, xsize, lsize

  integer                         :: i, j, p, buffer_pos, msgsize
  integer                         :: from_pe, to_pe, pos, n, l, count
  integer                         :: ibegin, istart, iend, start_pos
  type (comm_type), pointer, save :: comm =>NULL()
  real                            :: recv_buffer(xmap%put1%recvsize*lsize)
  real                            :: send_buffer(xmap%put1%sendsize*lsize)
  real                            :: unpack_buffer(xmap%put1%recvsize)

  real, dimension(isize, jsize)   :: d
  real, dimension(xsize)          :: x
  pointer(ptr_d, d)
  pointer(ptr_x, x)

  call mpp_clock_begin(id_put_1_to_xgrid_order_1)

  !--- pre-post receiving
  comm => xmap%put1
  do p = 1, comm%nrecv
     msgsize = comm%recv(p)%count*lsize
     from_pe = comm%recv(p)%pe
     buffer_pos = comm%recv(p)%buffer_pos*lsize
     call mpp_recv(recv_buffer(buffer_pos+1), glen=msgsize, from_pe = from_pe, block=.false., tag=COMM_TAG_7)
  enddo

  !--- send the data
  buffer_pos = 0
  do p = 1, comm%nsend
     msgsize = comm%send(p)%count*lsize
     to_pe = comm%send(p)%pe
     pos = buffer_pos
     do l = 1, lsize
        ptr_d = d_addrs(l)
        do n = 1, comm%send(p)%count
           pos = pos + 1
           i = comm%send(p)%i(n)
           j = comm%send(p)%j(n)
           send_buffer(pos) = d(i,j)
        enddo
     enddo
     call mpp_send(send_buffer(buffer_pos+1), plen=msgsize, to_pe = to_pe, tag=COMM_TAG_7 )
     buffer_pos = buffer_pos + msgsize
  enddo

  call mpp_sync_self(check=EVENT_RECV)

  !--- unpack the buffer
  if( lsize == 1) then
     ptr_x = x_addrs(1)
     do l=1,xmap%size_put1
        x(l) =  recv_buffer(xmap%x1_put(l)%pos)
     end do
  else
     start_pos = 0
!$OMP parallel do default(none) shared(lsize,x_addrs,comm,recv_buffer,xmap) &
!$OMP                          private(ptr_x,count,ibegin,istart,iend,pos,unpack_buffer)
     do l = 1, lsize
        ptr_x = x_addrs(l)
        do p = 1, comm%nrecv
           count = comm%recv(p)%count
           ibegin = comm%recv(p)%buffer_pos*lsize + 1
           istart = ibegin + (l-1)*count
           iend = istart + count - 1
           pos = comm%recv(p)%buffer_pos
           do n = istart, iend
              pos = pos + 1
              unpack_buffer(pos) = recv_buffer(n)
           enddo
        enddo
        do i=1,xmap%size_put1
           x(i) =  unpack_buffer(xmap%x1_put(i)%pos)
        end do    
     enddo
  endif

  call mpp_sync_self()

  call mpp_clock_end(id_put_1_to_xgrid_order_1)

end subroutine put_1_to_xgrid_order_1

!#######################################################################


subroutine put_1_to_xgrid_order_2(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)
  integer(LONG_KIND), dimension(:), intent(in) :: d_addrs
  integer(LONG_KIND), dimension(:), intent(in) :: x_addrs
  type (xmap_type),              intent(inout) :: xmap
  integer,                          intent(in) :: isize, jsize, xsize, lsize

  !: NOTE: halo size is assumed to be 1 in setup_xmap
  real, dimension(0:isize+1, 0:jsize+1, lsize) :: tmp
  real, dimension(isize,     jsize,     lsize) :: tmpx, tmpy
  real, dimension(isize,     jsize,     lsize) :: d_bar_max, d_bar_min
  real, dimension(isize,     jsize,     lsize) :: d_max, d_min
  real                            :: d_bar
  integer                         :: i, is, ie, im, j, js, je, jm, ii, jj
  integer                         :: p, l, ioff, joff, isd, jsd
  type (grid_type), pointer, save :: grid1 =>NULL()
  type (comm_type), pointer, save :: comm  =>NULL()
  integer                         :: buffer_pos, msgsize, from_pe, to_pe, pos, n
  integer                         :: ibegin, count, istart, iend
  real                            :: recv_buffer(xmap%put1%recvsize*lsize*3)
  real                            :: send_buffer(xmap%put1%sendsize*lsize*3)
  real                            :: unpack_buffer(xmap%put1%recvsize*3)
  logical                         :: on_west_edge, on_east_edge, on_south_edge, on_north_edge
  real, dimension(isize, jsize)   :: d
  real, dimension(xsize)          :: x
  pointer(ptr_d, d)
  pointer(ptr_x, x)

  call mpp_clock_begin(id_put_1_to_xgrid_order_2)
  grid1 => xmap%grids(1)

  is = grid1%is_me;   ie = grid1%ie_me
  js = grid1%js_me;   je = grid1%je_me
  isd = grid1%isd_me
  jsd = grid1%jsd_me

!$OMP parallel do default(none) shared(lsize,tmp,d_addrs,isize,jsize) private(ptr_d)
  do l = 1, lsize
     tmp(:,:,l) = LARGE_NUMBER
     ptr_d = d_addrs(l)
     tmp(1:isize,1:jsize,l) = d(:,:)
  enddo

  if(grid1%is_latlon) then
     call mpp_update_domains(tmp,grid1%domain_with_halo)
!$OMP parallel do default(none) shared(lsize,tmp,grid1,is,ie,js,je,isd,jsd,tmpx,tmpy)
     do l = 1, lsize
        tmpy(:,:,l) = grad_merid_latlon(tmp(:,:,l), grid1%lat, is, ie, js, je, isd, jsd)
        tmpx(:,:,l) = grad_zonal_latlon(tmp(:,:,l), grid1%lon, grid1%lat, is, ie, js, je, isd, jsd)
     enddo
  else
     call mpp_update_domains(tmp,grid1%domain)
     on_west_edge  = (is==1)
     on_east_edge  = (ie==grid1%im)
     on_south_edge = (js==1)
     on_north_edge = (je==grid1%jm)
!$OMP parallel do default(none) shared(lsize,tmp,grid1,tmpx,tmpy, &
!$OMP                                  on_west_edge,on_east_edge,on_south_edge,on_north_edge)
     do l = 1, lsize
        call gradient_cubic(tmp(:,:,l), grid1%box%dx, grid1%box%dy, grid1%box%area,   &
                            grid1%box%edge_w, grid1%box%edge_e, grid1%box%edge_s,     &
                            grid1%box%edge_n, grid1%box%en1, grid1%box%en2,           &
                            grid1%box%vlon, grid1%box%vlat, tmpx(:,:,l), tmpy(:,:,l), &
                            on_west_edge, on_east_edge, on_south_edge, on_north_edge)
     enddo
  end if     

  !--- pre-post receiving
  buffer_pos = 0
  comm => xmap%put1
  do p = 1, comm%nrecv
     msgsize = comm%recv(p)%count*lsize
     buffer_pos = comm%recv(p)%buffer_pos*lsize
     if(.NOT. monotonic_exchange) then
        msgsize = msgsize*3
        buffer_pos = buffer_pos*3
     endif
     from_pe = comm%recv(p)%pe
     call mpp_recv(recv_buffer(buffer_pos+1), glen=msgsize, from_pe = from_pe, block=.false., tag=COMM_TAG_8)
  enddo

  !--- compute d_bar_max and d_bar_min.
  if(monotonic_exchange) then
!$OMP parallel do default(none) shared(lsize,isize,jsize,d_bar_max,d_bar_min,d_max,d_min,tmp)
     do l = 1, lsize
        do j = 1, jsize
           do i = 1, isize
              d_bar_max(i,j,l) = -LARGE_NUMBER
              d_bar_min(i,j,l) =  LARGE_NUMBER
              d_max    (i,j,l) = -LARGE_NUMBER
              d_min    (i,j,l) =  LARGE_NUMBER          
              do jj = j-1, j+1
                 do ii = i-1, i+1
                    if(tmp(i,j,l) .NE. LARGE_NUMBER) then
                       if(tmp(i,j,l) > d_bar_max(i,j,l)) d_bar_max(i,j,l) = tmp(i,j,l)
                       if(tmp(i,j,l) < d_bar_min(i,j,l)) d_bar_min(i,j,l) = tmp(i,j,l)
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  !--- send the data
  buffer_pos = 0
  if(monotonic_exchange) then
     pos = 0
     do p = 1, comm%nsend
        msgsize = comm%send(p)%count*lsize
        to_pe = comm%send(p)%pe
        do l = 1, lsize
           ptr_d = d_addrs(l)
           do n = 1, comm%send(p)%count
              pos = pos + 1
              i = comm%send(p)%i(n)
              j = comm%send(p)%j(n)
              send_buffer(pos) = d(i,j) + tmpy(i,j,l)*comm%send(p)%dj(n) + tmpx(i,j,l)*comm%send(p)%di(n)
              if(send_buffer(pos) > d_max(i,j,l)) d_max(i,j,l) = send_buffer(pos)
              if(send_buffer(pos) < d_min(i,j,l)) d_min(i,j,l) = send_buffer(pos)              
           enddo
        enddo
     enddo

     do p = 1, comm%nsend
        msgsize = comm%send(p)%count*lsize
        to_pe = comm%send(p)%pe
        pos = buffer_pos
        do l = 1, lsize
           ptr_d = d_addrs(l)
           do n = 1, comm%send(p)%count
              pos = pos + 1
              i = comm%send(p)%i(n)
              j = comm%send(p)%j(n)           
              d_bar = d(i,j)
              if( d_max(i,j,l) > d_bar_max(i,j,l) ) then
                 send_buffer(pos) = d_bar + ((send_buffer(pos)-d_bar)/(d_max(i,j,l)-d_bar)) * (d_bar_max(i,j,l)-d_bar)
              else if( d_min(i,j,l) < d_bar_min(i,j,l) ) then
                 send_buffer(pos) = d_bar + ((send_buffer(pos)-d_bar)/(d_min(i,j,l)-d_bar)) * (d_bar_min(i,j,l)-d_bar)
              endif
           enddo
        enddo        
        call mpp_send(send_buffer(buffer_pos+1), plen=msgsize, to_pe = to_pe, tag=COMM_TAG_8 )
        buffer_pos = buffer_pos + msgsize
     enddo
  else
     do p = 1, comm%nsend
        msgsize = comm%send(p)%count*lsize*3
        to_pe = comm%send(p)%pe
        pos = buffer_pos
        do l = 1, lsize
           ptr_d = d_addrs(l)
           do n = 1, comm%send(p)%count
              pos = pos + 3
              i = comm%send(p)%i(n)
              j = comm%send(p)%j(n)
              send_buffer(pos-2) = d(i,j)
              send_buffer(pos-1) = tmpy(i,j,l)
              send_buffer(pos  ) = tmpx(i,j,l)
           enddo
        enddo
        call mpp_send(send_buffer(buffer_pos+1), plen=msgsize, to_pe = to_pe, tag=COMM_TAG_8 )
        buffer_pos = buffer_pos + msgsize
     enddo
  endif

  call mpp_sync_self(check=EVENT_RECV)

  !--- unpack the buffer
  if(monotonic_exchange) then
     if( lsize == 1) then
        ptr_x = x_addrs(1)
        do l=1,xmap%size_put1
           pos = xmap%x1_put(l)%pos
           x(l) =  recv_buffer(pos)
        end do
     else
        do l = 1, lsize
           ptr_x = x_addrs(l)
           pos = 0
           do p = 1, comm%nsend
              count = comm%send(p)%count
              ibegin = comm%recv(p)%buffer_pos*lsize + 1
              istart = ibegin + (l-1)*count
              iend = istart + count - 1
              pos = comm%recv(p)%buffer_pos
              do n = istart, iend
                 pos = pos + 1
                 unpack_buffer(pos) = recv_buffer(n)
              enddo
           enddo
           do i=1,xmap%size_put1
              pos = xmap%x1_put(i)%pos
              x(i) =  unpack_buffer(pos)
           end do
        enddo
     endif
  else
     if( lsize == 1) then
        ptr_x = x_addrs(1)
!$OMP parallel do default(none) shared(xmap,recv_buffer,ptr_x) private(pos)
        do l=1,xmap%size_put1
           pos = xmap%x1_put(l)%pos
           x(l) = recv_buffer(3*pos-2) + recv_buffer(3*pos-1)*xmap%x1_put(l)%dj + recv_buffer(3*pos)*xmap%x1_put(l)%di
        end do
     else
!$OMP parallel do default(none) shared(lsize,comm,xmap,recv_buffer,x_addrs) &
!$OMP                          private(ptr_x,pos,ibegin,istart,iend,count,unpack_buffer)
        do l = 1, lsize
           ptr_x = x_addrs(l)
           pos = 0
           ibegin = 1
           do p = 1, comm%nrecv
              count = comm%recv(p)%count*3
              ibegin = comm%recv(p)%buffer_pos*lsize*3 + 1
              istart = ibegin + (l-1)*count
              iend = istart + count - 1
              pos =  comm%recv(p)%buffer_pos*3
              do n = istart, iend
                 pos = pos + 1
                 unpack_buffer(pos) = recv_buffer(n)
              enddo
           enddo
           do i=1,xmap%size_put1
              pos = xmap%x1_put(i)%pos
              x(i) = unpack_buffer(3*pos-2) + unpack_buffer(3*pos-1)*xmap%x1_put(i)%dj + unpack_buffer(3*pos)*xmap%x1_put(i)%di
           end do
        enddo
     endif
  endif

  call mpp_sync_self()
  call mpp_clock_end(id_put_1_to_xgrid_order_2)

end subroutine put_1_to_xgrid_order_2

!#######################################################################

subroutine get_1_from_xgrid(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)
  integer(LONG_KIND), dimension(:), intent(in) :: d_addrs
  integer(LONG_KIND), dimension(:), intent(in) :: x_addrs
  type (xmap_type),              intent(inout) :: xmap
  integer,                          intent(in) :: isize, jsize, xsize, lsize

  real, dimension(xmap%size), target :: dg(xmap%size, lsize)
  integer                            :: i, j, l, p, n, m
  integer                            :: msgsize, buffer_pos, pos
  integer                            :: istart, iend, count
  real              , pointer, save  :: dgp =>NULL()
  type  (grid_type) , pointer, save  :: grid1 =>NULL()
  type  (comm_type) , pointer, save  :: comm  =>NULL()
  type(overlap_type), pointer, save  :: send => NULL()
  type(overlap_type), pointer, save  :: recv => NULL()
  real                               :: recv_buffer(xmap%get1%recvsize*lsize*3)
  real                               :: send_buffer(xmap%get1%sendsize*lsize*3)
  real                               :: unpack_buffer(xmap%get1%recvsize*3)
  real                               :: d(isize,jsize)
  real, dimension(xsize)             :: x
  pointer(ptr_d, d)
  pointer(ptr_x, x)

  call mpp_clock_begin(id_get_1_from_xgrid)

  comm => xmap%get1
  grid1 => xmap%grids(1)

  do p = 1, comm%nrecv
     recv => comm%recv(p)
     msgsize = recv%count*lsize
     buffer_pos = recv%buffer_pos*lsize
     call mpp_recv(recv_buffer(buffer_pos+1), glen=msgsize, from_pe = recv%pe, block=.false., tag=COMM_TAG_9)
  enddo

  dg = 0.0;
!$OMP parallel do default(none) shared(lsize,xmap,dg,x_addrs) private(dgp,ptr_x)
  do l = 1, lsize
     ptr_x = x_addrs(l)
     do i=1,xmap%size
        dgp => dg(xmap%x1(i)%pos,l)
        dgp =  dgp + xmap%x1(i)%area*x(i)
     enddo
  enddo


  !--- send the data
  buffer_pos = 0
  istart = 1
  do p = 1, comm%nsend
     send => comm%send(p)
     msgsize = send%count*lsize
     pos = buffer_pos
     istart = send%buffer_pos+1
     iend = istart + send%count - 1
     do l = 1, lsize
        do n = istart, iend
           pos = pos + 1
           send_buffer(pos) = dg(n,l)
        enddo
     enddo
     call mpp_send(send_buffer(buffer_pos+1), plen=msgsize, to_pe = send%pe, tag=COMM_TAG_9 )
     buffer_pos = buffer_pos + msgsize
     istart = iend + 1
  enddo

  call mpp_sync_self(check=EVENT_RECV)

  !--- unpack the buffer
  do l = 1, lsize
     ptr_d = d_addrs(l)
     d = 0.0
  enddo
  !--- To bitwise reproduce old results, first copy the data onto its own pe.

  do p = 1, comm%nrecv
     recv => comm%recv(p)
     count = recv%count
     buffer_pos = recv%buffer_pos*lsize
     if( recv%pe == xmap%me ) then
!$OMP parallel do default(none) shared(lsize,recv,recv_buffer,buffer_pos,d_addrs,count) &
!$OMP                          private(ptr_d,i,j,pos)
        do l = 1, lsize
           pos = buffer_pos + (l-1)*count
           ptr_d = d_addrs(l)
           do n = 1,count
              i = recv%i(n)
              j = recv%j(n)
              pos = pos + 1
              d(i,j) = recv_buffer(pos) 
           enddo
        enddo
        exit
     endif
  enddo

  pos = 0
  do m = 1, comm%nrecv
     p = comm%unpack_ind(m)
     recv => comm%recv(p)
     if( recv%pe == xmap%me ) then
        cycle
     endif
     buffer_pos = recv%buffer_pos*lsize
!$OMP parallel do default(none) shared(lsize,recv,recv_buffer,buffer_pos,d_addrs) &
!$OMP                          private(ptr_d,i,j,pos)
     do l = 1, lsize
        pos = buffer_pos + (l-1)*recv%count
        ptr_d = d_addrs(l)
        do n = 1, recv%count
           i = recv%i(n)
           j = recv%j(n)
           pos = pos + 1
           d(i,j) = d(i,j) + recv_buffer(pos) 
        enddo
     enddo
  enddo

  !
  ! normalize with side 1 grid cell areas
  !
!$OMP parallel do default(none) shared(lsize,d_addrs,grid1) private(ptr_d)
  do l = 1, lsize
     ptr_d = d_addrs(l)
     d = d * grid1%area_inv
  enddo
  call mpp_sync_self()
  call mpp_clock_end(id_get_1_from_xgrid)

end subroutine get_1_from_xgrid

!#######################################################################

subroutine get_1_from_xgrid_repro(d_addrs, x_addrs, xmap, xsize, lsize)
  integer(LONG_KIND), dimension(:), intent(in) :: d_addrs
  integer(LONG_KIND), dimension(:), intent(in) :: x_addrs
  type (xmap_type),              intent(inout) :: xmap
  integer,                          intent(in) :: xsize, lsize

  integer                            :: g, i, j, k, p, l, n, l2, m, l3
  integer                            :: msgsize, buffer_pos, pos
  type (grid_type), pointer, save :: grid =>NULL()
  type(comm_type),  pointer, save :: comm => NULL()
  type(overlap_type), pointer, save  :: send => NULL()
  type(overlap_type), pointer, save  :: recv => NULL()
    integer,  dimension(0:xmap%npes-1) :: pl, ml
  real                               :: recv_buffer(xmap%recv_count_repro_tot*lsize)
  real                               :: send_buffer(xmap%send_count_repro_tot*lsize)
  real                               :: d(xmap%grids(1)%is_me:xmap%grids(1)%ie_me, &
                                          xmap%grids(1)%js_me:xmap%grids(1)%je_me)
  real, dimension(xsize)             :: x
  pointer(ptr_d, d)
  pointer(ptr_x, x)

  call mpp_clock_begin(id_get_1_from_xgrid_repro)
  comm => xmap%get1_repro
  !--- pre-post receiving
  do p = 1, comm%nrecv
     recv => comm%recv(p)
     msgsize = recv%count*lsize
     buffer_pos = recv%buffer_pos*lsize
     call mpp_recv(recv_buffer(buffer_pos+1), glen=msgsize, from_pe = recv%pe, block=.false., tag=COMM_TAG_10)
     n = recv%pe -xmap%root_pe
     pl(n) = buffer_pos 
     ml(n) = recv%count
  enddo

  !pack the data
  send_buffer(:) = 0.0
!$OMP parallel do default(none) shared(lsize,x_addrs,comm,xmap,send_buffer) &
!$OMP                          private(ptr_x,i,j,g,l2,pos,send)
  do p = 1, comm%nsend
     pos = comm%send(p)%buffer_pos*lsize
     send => comm%send(p)
     do l = 1,lsize
        ptr_x = x_addrs(l)
        do n = 1, send%count
           i = send%i(n)
           j = send%j(n)
           g = send%g(n)
           l2 = send%xloc(n)
           pos = pos + 1
           do k =1, xmap%grids(g)%km
             if(xmap%grids(g)%frac_area(i,j,k)/=0.0) then
              l2 = l2+1
              send_buffer(pos) = send_buffer(pos) + xmap%x1(l2)%area *x(l2)
             endif
           enddo
         enddo
      enddo
   enddo

  do p =1, comm%nsend
     buffer_pos = comm%send(p)%buffer_pos*lsize
     msgsize = comm%send(p)%count*lsize
     call mpp_send(send_buffer(buffer_pos+1), plen=msgsize, to_pe=comm%send(p)%pe, tag=COMM_TAG_10)
  enddo

  do l = 1, lsize
     ptr_d = d_addrs(l)
     d = 0
  enddo

  call mpp_sync_self(check=EVENT_RECV)

!$OMP parallel do default(none) shared(lsize,d_addrs,xmap,recv_buffer,pl,ml) &
!$OMP                          private(ptr_d,grid,i,j,p,pos)
  do l = 1, lsize
     ptr_d = d_addrs(l)
     do g=2,size(xmap%grids(:))
        grid => xmap%grids(g)
        do l3=1,grid%size_repro ! index into side1 grid's patterns
           i = grid%x_repro(l3)%i1
           j = grid%x_repro(l3)%j1
           p = grid%x_repro(l3)%pe-xmap%root_pe
           pos = pl(p) + (l-1)*ml(p) + grid%x_repro(l3)%recv_pos
           d(i,j) = d(i,j) + recv_buffer(pos)
        end do
     end do
     ! normalize with side 1 grid cell areas
     d = d * xmap%grids(1)%area_inv
  enddo

  call mpp_sync_self()       

  call mpp_clock_end(id_get_1_from_xgrid_repro)

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
  conservation_check_side1 = 0.0
  if(grid1%tile_me .NE. tile_nest) conservation_check_side1(1) = sum(grid1%area*d)
!  if(grid1%tile_me .NE. tile_parent .OR. grid1%id .NE. "ATM") &
!      conservation_check_side1(1) = sum(grid1%area*d) 

  call put_to_xgrid (d, grid1%id, x_over, xmap, remap_method)    ! put from side 1
  do g=2,size(xmap%grids(:))
    grid2 => xmap%grids(g)
    if(grid2%on_this_pe) then
      allocate (d2 (grid2%is_me:grid2%ie_me, grid2%js_me:grid2%je_me,  grid2%km) )
    endif
    call get_from_xgrid (d2, grid2%id, x_over, xmap) ! get onto side 2's
    if(grid2%on_this_pe) then
       conservation_check_side1(2) = conservation_check_side1(2) + sum( grid2%area * sum(grid2%frac_area*d2,DIM=3) )
    endif
    call put_to_xgrid (d2, grid2%id, x_back, xmap) ! put from side 2's
    if(allocated(d2))deallocate (d2)
  end do
  call get_from_xgrid(d1, grid1%id, x_back, xmap)  ! get onto side 1
  if(grid1%tile_me .NE. tile_nest) conservation_check_side1(3) = sum(grid1%area*d1)
!  if(grid1%tile_me .NE. tile_parent .OR. grid1%id .NE. "ATM") &
!     conservation_check_side1(3) = sum(grid1%area*d1)
  call mpp_sum(conservation_check_side1,3)

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
  conservation_check_side2 = 0.0
  do g = 2,size(xmap%grids(:))
    grid2 => xmap%grids(g)
    if (grid_id==grid2%id) then
      if(grid2%on_this_pe) then
         conservation_check_side2(1) = sum( grid2%area * sum(grid2%frac_area*d,DIM=3) )
      endif
      call put_to_xgrid(d, grid_id, x_over, xmap)  ! put from this side 2
    else
      call put_to_xgrid(0.0 * grid2%frac_area, grid2%id, x_over, xmap) ! zero rest
    end if
  end do

  allocate ( d1(size(grid1%area,1),size(grid1%area,2)) )
  call get_from_xgrid(d1, grid1%id, x_over, xmap)  ! get onto side 1
  if(grid1%tile_me .NE. tile_nest)conservation_check_side2(2) = sum(grid1%area*d1)
  call put_to_xgrid(d1,  grid1%id, x_back, xmap,remap_method)   ! put from side 1
  deallocate ( d1 )

  conservation_check_side2(3) = 0.0;
  do g = 2,size(xmap%grids(:))
    grid2 => xmap%grids(g)
    if(grid2%on_this_pe) then
       allocate ( d2 ( size(grid2%frac_area, 1), size(grid2%frac_area, 2),  &
                                                 size(grid2%frac_area, 3) ) )
    endif
    call get_from_xgrid(d2,  grid2%id, x_back, xmap) ! get onto side 2's
    conservation_check_side2(3) = conservation_check_side2(3) + sum( grid2%area * sum(grid2%frac_area*d2,DIM=3) )
    if(allocated(d2) )deallocate ( d2 )
  end do
  call mpp_sum(conservation_check_side2, 3)  

end function conservation_check_side2
! </FUNCTION>

!******************************************************************************
! This routine is used to get the grid area of component model with id.
subroutine get_xmap_grid_area(id, xmap, area)
  character(len=3),     intent(in   ) :: id
  type (xmap_type),     intent(inout) :: xmap
  real, dimension(:,:), intent(out  ) :: area
  integer                             :: g
  logical                             :: found

   found = .false.
   do g = 1, size(xmap%grids(:))
      if (id==xmap%grids(g)%id ) then
         if(size(area,1) .NE. size(xmap%grids(g)%area,1) .OR. size(area,2) .NE. size(xmap%grids(g)%area,2) ) &
           call error_mesg("xgrid_mod", "size mismatch between area and xmap%grids(g)%area", FATAL)
         area = xmap%grids(g)%area
         found = .true.
         exit
      end if
   end do

   if(.not. found) call error_mesg("xgrid_mod", id//" is not found in xmap%grids id", FATAL)

end subroutine get_xmap_grid_area

!#######################################################################

! This function is used to calculate the gradient along zonal direction.
! Maybe need to setup a limit for the gradient. The grid is assumeed 
! to be regular lat-lon grid

function grad_zonal_latlon(d, lon, lat, is, ie, js, je, isd, jsd) 

  integer,                    intent(in) :: isd, jsd
  real, dimension(isd:,jsd:), intent(in) :: d
  real, dimension(:),         intent(in) :: lon
  real, dimension(:),         intent(in) :: lat
  integer,                    intent(in) :: is, ie, js, je 
  real, dimension(is:ie,js:je)           :: grad_zonal_latlon
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
     if(abs(dx).lt.EPS )  call error_mesg('xgrids_mod(grad_zonal_latlon)', 'Improper grid size in lontitude', FATAL)
     if(dx .gt. PI)  dx = dx - 2.0* PI
     if(dx .lt. -PI) dx = dx + 2.0* PI
     do j = js, je
        costheta = cos(lat(j))
        if(abs(costheta) .lt. EPS) call error_mesg('xgrids_mod(grad_zonal_latlon)', 'Improper latitude grid', FATAL)
        grad_zonal_latlon(i,j) = (d(ip1,j)-d(im1,j))/(dx*costheta)
     enddo
  enddo

  return

end function grad_zonal_latlon

!#######################################################################

! This function is used to calculate the gradient along meridinal direction.
! Maybe need to setup a limit for the gradient. regular lat-lon grid are assumed

function grad_merid_latlon(d, lat, is, ie, js, je, isd, jsd) 
  integer,                    intent(in) :: isd, jsd
  real, dimension(isd:,jsd:), intent(in) :: d
  real, dimension(:),         intent(in) :: lat
  integer,                    intent(in) :: is, ie, js, je 
  real, dimension(is:ie,js:je)           :: grad_merid_latlon
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
     if(abs(dy).lt.EPS) call error_mesg('xgrids_mod(grad_merid_latlon)', 'Improper grid size in latitude', FATAL)

     do i = is, ie
        grad_merid_latlon(i,j) = (d(i,jp1) - d(i,jm1))/dy
     enddo
  enddo

  return
end function grad_merid_latlon

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

     from_dq = delta_t * 4.0*PI*radius**2 * sum( sum(xmap%grids(grid_index)%area * &
          & sum(xmap%grids(grid_index)%frac_area * data, DIM=3), DIM=1))
     to_dq = from_dq

  ! update only if argument is present.
  if(present(to  )) to   % dq(  to_side) = to   % dq(  to_side) + to_dq
  if(present(from)) from % dq(from_side) = from % dq(from_side) - from_dq

  if(present(verbose).and.debug_stocks) then
     call mpp_sum(from_dq)
     call mpp_sum(to_dq)
     from_dq = from_dq/(4.0*PI*radius**2)
     to_dq   = to_dq  /(4.0*PI*radius**2)
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
     from_dq = delta_t * 4.0*PI*radius**2 * sum(sum(xmap%grids(1)%area * data, DIM=1))
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
     from_dq = from_dq/(4.0*PI*radius**2)
     to_dq   = to_dq  /(4.0*PI*radius**2)
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

  ier = 0
  res = 0.0

  if(.not. associated(xmap%grids) ) then
     ier = 6
     return
  endif

  res = delta_t * 4.0*PI*radius**2 * sum(sum(xmap%grids(1)%area * data, DIM=1))

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

  integer, parameter :: initID = -2 ! initial value for diag IDs. Must not be equal to the value 
  ! that register_diag_field returns when it can't register the filed -- otherwise the registration 
  ! is attempted every time this subroutine is called

  real :: f_value, c_value, planet_area
  character(len=80) :: formatString
  integer :: iday, isec, hours
  integer :: diagID, compInd
  integer, dimension(NELEMS,4), save :: f_valueDiagID = initID
  integer, dimension(NELEMS,4), save :: c_valueDiagID = initID
  integer, dimension(NELEMS,4), save :: fmc_valueDiagID = initID

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
     planet_area = 4.0*PI*radius**2
     f_value       = f_value     / planet_area
     c_value       = c_value     / planet_area

     if(comp_name == 'ATM') compInd = 1
     if(comp_name == 'LND') compInd = 2
     if(comp_name == 'ICE') compInd = 3
     if(comp_name == 'OCN') compInd = 4


     if(f_valueDiagID(index,compInd) == initID) then
        field_name = trim(comp_name) // trim(STOCK_NAMES(index))
        field_name  = trim(field_name) // 'StocksChange_Flux'
        units = trim(STOCK_UNITS(index))
        f_valueDiagID(index,compInd) = register_diag_field('stock_print', field_name, Time, &
             units=units)
     endif

     if(c_valueDiagID(index,compInd) == initID) then
        field_name = trim(comp_name) // trim(STOCK_NAMES(index))
        field_name = trim(field_name) // 'StocksChange_Comp'
        units = trim(STOCK_UNITS(index))
        c_valueDiagID(index,compInd) = register_diag_field('stock_print', field_name, Time, &
             units=units)
     endif

     if(fmc_valueDiagID(index,compInd) == initID) then
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
    integer                          :: i, j, nlon, nlat, num

    is_lat_lon = .true.
    nlon = size(lon,1)
    nlat = size(lon,2)
    LOOP_LAT: do j = 1, nlat
       do i = 2, nlon
          if(lat(i,j) .NE. lat(1,j)) then
             is_lat_lon = .false.
             exit LOOP_LAT
          end if
       end do
    end do LOOP_LAT

    if(is_lat_lon) then
       LOOP_LON: do i = 1, nlon
          do j = 2, nlat
             if(lon(i,j) .NE. lon(i,1)) then
                is_lat_lon = .false.
                exit LOOP_LON
             end if
          end do
       end do LOOP_LON
    end if

    num = 0
    if(is_lat_lon) num = 1
    call mpp_min(num)
    if(num == 1) then
       is_lat_lon = .true.
    else
       is_lat_lon = .false.
    end if

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
