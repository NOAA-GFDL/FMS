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
!                                    Michael Winton (mw@gfdl.noaa.gov) Oct 2001
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module xgrid_mod

! <CONTACT EMAIL="mw@gfdl.noaa.gov">
!   Michael Winton
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
                           error_mesg, close_file, FATAL, stdlog,            &
                           write_version_number 
use mpp_mod,         only: mpp_npes, mpp_pe, mpp_root_pe, mpp_send, mpp_recv, &
                           mpp_sync_self
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                           Domain2d, mpp_global_sum, mpp_update_domains,    &
                           mpp_copy_domains, mpp_get_data_domain
use mpp_io_mod,      only: mpp_open, MPP_MULTI, MPP_OVERWR
use constants_mod,   only: PI


implicit none
include 'netcdf.inc'
private

public xmap_type, setup_xmap, set_frac_area, put_to_xgrid, get_from_xgrid, &
       xgrid_count, some, conservation_check, xgrid_init

! <NAMELIST NAME="xgrid_nml">
!   <DATA NAME="make_exchange_reproduce" TYPE="logical"  DEFAULT=".false.">
!     Set to .true. to make <TT>xgrid_mod</TT> reproduce answers on different
!     numbers of PEs.  This option has a considerable performance impact.
!   </DATA>
!   <DATA NAME="use_higher_order" TYPE="logical"  DEFAULT=".false.">
!     set to .ture. to use second order remapping scheme, also need to set
!     remap_order of flux_exchange_nml to 2. 
!   </DATA>
logical :: make_exchange_reproduce = .false. ! exactly same on different # PEs
logical :: use_higher_order = .false.
logical :: debug = .false.
namelist /xgrid_nml/ make_exchange_reproduce, use_higher_order, debug
! </NAMELIST>
logical :: init = .true.


! <INTERFACE NAME="put_to_xgrid">

!   <OVERVIEW>
!     Scatters data from model grid onto exchange grid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Scatters data from model grid onto exchange grid.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call put_to_xgrid(d, grid_id, x, xmap, remapping_order)
!   </TEMPLATE>
!   <IN NAME="d"  TYPE="real"  > </IN>
!   <IN NAME="grid_id"  TYPE=" character(len=3)"  > </IN>
!   <INOUT NAME="x"  TYPE="real"  > </INOUT>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>
!   <IN NAME="remapping_order" TYPE="integer,optional">
!     default value is 1 (linear interpolation). =2 use second order
!     interpolation.
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
!     call conservation_check(d, grid_id, xmap,remapping_order)
!   </TEMPLATE>
!   <IN NAME="d"  TYPE="real" DIM="(:,:)" > </IN>
!   <IN NAME="grid_id"  TYPE="character(len=3)"  > </IN>
!   <INOUT NAME="xmap"  TYPE="xmap_type"  > </INOUT>
!   <OUT NAME="" TYPE="real" DIM="3">The global sum of a variable.</OUT>
!   <IN NAME="remapping_order" TYPE="integer,optional">
!     default value is 1 (linear interpolation). =2 use second order
!     interpolation.
!   </IN>
interface conservation_check
  module procedure conservation_check_side1
  module procedure conservation_check_side2
end interface
! </INTERFACE>

type xcell_type
  integer :: i1, j1, i2, j2 ! indices of cell in model arrays on both sides
  integer :: pe             ! other side pe that has this cell
  real    :: area           ! geographic area of exchange cell
  real    :: di, dj         ! Weight for the gradient of flux
end type xcell_type

type grid_type
  character(len=3)                :: id           ! grid identifier
  integer, pointer, dimension(:)  :: is, ie       ! domain - i-range (pe index)
  integer, pointer, dimension(:)  :: js, je       ! domain - j-range (pe index)
  integer, pointer                :: is_me, ie_me ! my domain - i-range
  integer, pointer                :: js_me, je_me ! my domain - j-range
  integer                         :: im, jm, km   ! global domain range
  real, pointer, dimension(:,:)   :: lon, lat     ! center of global grids
  real, pointer, dimension(:,:,:) :: frac_area    ! partition fractions
  real, pointer, dimension(:,:)   :: area         ! cell area
  real, pointer, dimension(:,:)   :: area_inv     ! 1 / area for normalization
  integer                         :: first, last  ! xgrid index range
  integer                         :: size         ! # xcell patterns
  type(xcell_type), pointer, dimension(:) :: x    ! xcell patterns
  integer                         :: size_repro   ! # side 1 patterns for repro
  type(xcell_type), pointer, dimension(:) :: x_repro ! side 1 patterns for repro
  type(Domain2d) :: domain                        ! used for conservation checks
  type(Domain2d) :: domain_with_halo              ! used for getting data halo
end type grid_type

type x1_type
  integer :: i, j
  real    :: area   ! (= geographic area * frac_area)
  real    :: di, dj ! weight for the gradient of flux
end type x1_type

type x2_type
  integer :: i, j, k
  real    :: area   ! geographic area of exchange cell
end type x2_type

type xmap_type
  private
  integer :: size            ! # of exchange grid cells with area > 0 on this pe

  integer :: me, npes, root_pe
  logical, pointer, dimension(:) :: your1my2 ! true if side 1 domain on
                                             ! indexed pe overlaps side 2
                                             ! domain on this pe
  logical, pointer, dimension(:) :: your2my1 ! true if a side 2 domain on
                                             ! indexed pe overlaps side 1
                                             ! domain on this pe

  type (grid_type), pointer, dimension(:) :: grids ! 1st grid is side 1;
                                                   ! rest on side 2
  !
  ! Description of the individual exchange grid cells (index is cell #)
  !
  type(x1_type), pointer, dimension(:) :: x1 ! side 1 info
  type(x2_type), pointer, dimension(:) :: x2 ! side 2 info

  real, pointer, dimension(:) :: send_buffer ! for non-blocking sends
  integer, pointer, dimension(:) :: send_count_repro, recv_count_repro
end type xmap_type

!-----------------------------------------------------------------------
 character(len=128) :: version = '$Id: xgrid.f90,v 1.4 2003/04/09 21:16:19 fms Exp $'
 character(len=128) :: tagname = '$Name: inchon $'

 logical :: module_is_initialized = .FALSE.

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

subroutine xgrid_init 

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

  
end subroutine xgrid_init
! </SUBROUTINE>

!#######################################################################

subroutine load_xgrid (xmap, grid, domain, ncid, id_i1, id_j1, id_i2, id_j2, &
                           id_area, n_areas, id_di, id_dj  )
type(xmap_type), intent(inout)  :: xmap
type(grid_type), intent(inout)  :: grid
type(Domain2d), intent(inout)   :: domain
integer, intent(in) :: ncid, id_i1, id_j1, id_i2, id_j2, id_area, n_areas
integer, intent(in), optional :: id_di, id_dj 

  integer, dimension(n_areas) :: i1, j1, i2, j2 ! xgrid quintuples
  real,    dimension(n_areas) :: area, di, dj     ! from grid file
  type (grid_type), pointer   :: grid1
  integer, dimension(0:xmap%npes-1) :: is_2, ie_2, js_2, je_2 ! side 2 decomp.
  integer :: start(4), nread(4), rcode, l, ll, ll_repro, p

  grid1 => xmap%grids(1)
  start = 1; nread = 1; nread(1) = n_areas
  rcode = nf_get_vara_int(ncid, id_i1, start, nread, i1)
  rcode = nf_get_vara_int(ncid, id_j1, start, nread, j1)
  rcode = nf_get_vara_int(ncid, id_i2, start, nread, i2)
  rcode = nf_get_vara_int(ncid, id_j2, start, nread, j2)
  rcode = nf_get_vara_double(ncid, id_area, start, nread, area)
  di = 0.0;  dj = 0.0
  if(present(id_di))  rcode = nf_get_vara_double(ncid, id_di, start, nread, di)
  if(present(id_dj))  rcode = nf_get_vara_double(ncid, id_dj, start, nread, dj)
  

  do l=1,n_areas
    if (in_box(i1(l), j1(l), grid1%is_me, grid1%ie_me, &
                             grid1%js_me, grid1%je_me) ) then
      grid1%area(i1(l),j1(l)) = grid1%area(i1(l),j1(l))+area(l)
      do p=0,xmap%npes-1
        if (in_box(i2(l), j2(l), grid%is(p), grid%ie(p), &
                                 grid%js(p), grid%je(p)))  then
          xmap%your2my1(p) = .true.
        end if
      end do
      grid%size_repro = grid%size_repro + 1
    end if
    if (in_box(i2(l), j2(l), grid%is_me, grid%ie_me, &
                             grid%js_me, grid%je_me) ) then
      grid%size = grid%size + 1
      grid%area(i2(l),j2(l)) = grid%area(i2(l),j2(l))+area(l)
      do p=0,xmap%npes-1
        if (in_box(i1(l), j1(l), grid1%is(p), grid1%ie(p), &
                                 grid1%js(p), grid1%je(p))) then
          xmap%your1my2(p) = .true.
        end if
      end do
    end if
  end do

  allocate( grid%x( grid%size ) )
  if (make_exchange_reproduce) allocate ( grid%x_repro(grid%size_repro) )
  ll = 0
  ll_repro = 0
  do l=1,n_areas
    if (in_box(i2(l), j2(l), grid%is_me, grid%ie_me, grid%js_me, grid%je_me)) then
      ! insert in this grids cell pattern list and add area to side 2 area
      ll = ll + 1
      grid%x(ll)%i1   = i1(l); grid%x(ll)%i2   = i2(l)
      grid%x(ll)%j1   = j1(l); grid%x(ll)%j2   = j2(l)
      grid%x(ll)%area = area(l)
      grid%x(ll)%di  = di(l)
      grid%x(ll)%dj  = dj(l)

      if (make_exchange_reproduce) then
        do p=0,xmap%npes-1
          if (in_box(i1(l), j1(l), grid1%is(p), grid1%ie(p), &
                                   grid1%js(p), grid1%je(p))) then
            grid%x(ll)%pe = p + xmap%root_pe
          end if
        end do
      end if ! make_exchange reproduce
    end if
    if (in_box(i1(l),j1(l), grid1%is_me,grid1%ie_me, grid1%js_me,grid1%je_me) &
        .and. make_exchange_reproduce                                     ) then
      ll_repro = ll_repro + 1
      grid%x_repro(ll_repro)%i1   = i1(l); grid%x_repro(ll_repro)%i2   = i2(l)
      grid%x_repro(ll_repro)%j1   = j1(l); grid%x_repro(ll_repro)%j2   = j2(l)
      grid%x_repro(ll_repro)%area = area(l)
      grid%x_repro(ll_repro)%di  = di(l)
      grid%x_repro(ll_repro)%dj  = dj(l)
      do p=0,xmap%npes-1
        if (in_box(i2(l), j2(l), grid%is(p), grid%ie(p), &
                                 grid%js(p), grid%je(p))) then
          grid%x_repro(ll_repro)%pe = p + xmap%root_pe
        end if
      end do
    end if ! make_exchange_reproduce
  end do
  grid%area_inv = 0.0;
  where (grid%area>0.0) grid%area_inv = 1.0/grid%area
end subroutine load_xgrid

!#######################################################################
!
! get_grid - read the center point of the grid from grid_spec.nc.
!
!

subroutine get_grid(grid, grid_id, ncid)
type(grid_type), intent(inout) :: grid
integer,         intent(in)    :: ncid
character(len=3), intent(in)   :: grid_id

   integer :: rcode, start(4), nread(4), id_lon, id_lat, i, j
   real, dimension(grid%im) :: lonb
   real, dimension(grid%jm) :: latb
   real, dimension(grid%im,grid%jm) :: lon, lat
   real ::  d2r

   d2r = PI/180.0

   start = 1; nread = 1  

   if(grid_id == 'ATM') then
      nread(1) =  grid%im     
      rcode = nf_inq_varid(ncid,'xta',id_lon)
      if (rcode/=0) call error_mesg('xgrid_mod', 'cannot find grid file field xta', FATAL)
      rcode = nf_get_vara_double(ncid, id_lon, start, nread,lonb)

      rcode = nf_inq_varid(ncid,'yta',id_lat)
      if (rcode/=0) call error_mesg('xgrid_mod', 'cannot find grid file field yta', FATAL)
      nread(1) =  grid%jm
      rcode = nf_get_vara_double(ncid, id_lat, start, nread,latb)

      do j= 1, grid%jm
         grid%lon(:,j) = lonb * d2r
      enddo
      do i=1, grid%im
         grid%lat(i,:) = latb * d2r
      enddo

   else if(grid_id == 'LND') then
      nread(1) =  grid%im     
      rcode = nf_inq_varid(ncid,'xtl',id_lon)
      if (rcode/=0) call error_mesg('xgrid_mod', 'cannot find grid file field xtl', FATAL)
      rcode = nf_get_vara_double(ncid, id_lon, start, nread,lonb)

      nread(1) =  grid%jm
      rcode = nf_inq_varid(ncid,'ytl',id_lat)
      if (rcode/=0) call error_mesg('xgrid_mod', 'cannot find grid file field ytl', FATAL)
      rcode = nf_get_vara_double(ncid, id_lat, start, nread,latb)

      do j= 1, grid%jm
         grid%lon(:,j) = lonb * d2r
      enddo

      do i=1, grid%im
         grid%lat(i,:) = latb * d2r
      enddo

   else if(grid_id == 'OCN') then

      nread(1) =  grid%im;  nread(2) =  grid%jm      
      rcode = nf_inq_varid(ncid,'geolon_t',id_lon)      
      if (rcode/=0) call error_mesg('xgrid_mod', 'cannot find grid file field geolon_t', FATAL)
      rcode = nf_get_vara_double(ncid, id_lon, start, nread,lon)
      grid%lon = lon * d2r

      rcode = nf_inq_varid(ncid,'geolat_t',id_lat) 
      if (rcode/=0) call error_mesg('xgrid_mod', 'cannot find grid file field geolat_t', FATAL)
      rcode = nf_get_vara_double(ncid, id_lat, start, nread,lat)
      grid%lat = lat * d2r
   else
      call error_mesg('xgrid_mod','setup_xmap: could not find grid id', FATAL)
   endif

return

end subroutine get_grid
  


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

subroutine setup_xmap(xmap, grid_ids, grid_domains, grid_file)
type (xmap_type),                          intent(out) :: xmap
character(len=3), dimension(:),            intent(in ) :: grid_ids
type(Domain2d), dimension(size(grid_ids)), intent(in ) :: grid_domains
character(len=*)                         , intent(in ) :: grid_file

  integer :: g, l, ll, p, n_areas, send_size
  integer :: ncid, i1_id, j1_id, i2_id, j2_id, area_id, di_id, dj_id
  integer :: dims(4), rcode
  integer :: unit
  type (grid_type), pointer :: grid, grid1
  real, dimension(3) :: xxx
  real, dimension(:,:), allocatable :: check_data
  real, dimension(:,:,:), allocatable :: check_data_3D
  integer :: i, j

  xmap%me   = mpp_pe  ()
  xmap%npes = mpp_npes()
  xmap%root_pe = mpp_root_pe()

  allocate( xmap%grids(1:size(grid_ids)) )

  allocate ( xmap%your1my2(0:xmap%npes-1), xmap%your2my1(0:xmap%npes-1) )

  xmap%your1my2 = .false.; xmap%your2my1 = .false.;

  rcode = nf_open(grid_file,0,ncid)
  if (rcode/=0) call error_mesg ('xgrid_mod', 'cannot open grid file', FATAL)

  do g=1,size(grid_ids)
    grid => xmap%grids(g)
    if (g==1) grid1 => xmap%grids(g)
    grid%id     = grid_ids    (g)
    grid%domain = grid_domains(g)
    call mpp_copy_domains(grid%domain, grid%domain_with_halo, xhalo = 1, yhalo =1)

    allocate ( grid%is(0:xmap%npes-1), grid%ie(0:xmap%npes-1) )
    allocate ( grid%js(0:xmap%npes-1), grid%je(0:xmap%npes-1) )
    call mpp_get_compute_domains(grid%domain, xbegin=grid%is, xend=grid%ie, &
                                              ybegin=grid%js, yend=grid%je  )
    grid%is_me => grid%is(xmap%me-xmap%root_pe); grid%ie_me => grid%ie(xmap%me-xmap%root_pe)
    grid%js_me => grid%js(xmap%me-xmap%root_pe); grid%je_me => grid%je(xmap%me-xmap%root_pe)
    grid%im = maxval(grid%ie)
    grid%jm = maxval(grid%je)
    grid%km = 1

    allocate( grid%area    (grid%is_me:grid%ie_me, grid%js_me:grid%je_me) )
    allocate( grid%area_inv(grid%is_me:grid%ie_me, grid%js_me:grid%je_me) )
    grid%area       = 0.0
    grid%size       = 0
    grid%size_repro = 0
    if (g>1) then
      allocate( grid%frac_area(grid%is_me:grid%ie_me, grid%js_me:grid%je_me, &
                                                      grid%km              ) )
      grid%frac_area = 1.0

      rcode = nf_inq_varid(ncid, &
                         'I_'//grid_ids(1)//'_'//grid_ids(1)//'x'//grid_ids(g),&
                         i1_id)
      if (rcode/=0) &
        call error_mesg('xgrid_mod', 'cannot find grid file field', FATAL)
      rcode = nf_inq_varid(ncid, &
                         'J_'//grid_ids(1)//'_'//grid_ids(1)//'x'//grid_ids(g),&
                         j1_id)
      if (rcode/=0) &
        call error_mesg('xgrid_mod', 'cannot find grid file field', FATAL)
      rcode = nf_inq_varid(ncid, &
                         'I_'//grid_ids(g)//'_'//grid_ids(1)//'x'//grid_ids(g),&
                         i2_id)
      if (rcode/=0) &
        call error_mesg('xgrid_mod', 'cannot find grid file field', FATAL)
      rcode = nf_inq_varid(ncid, &
                         'J_'//grid_ids(g)//'_'//grid_ids(1)//'x'//grid_ids(g),&
                         j2_id)
      if (rcode/=0) &
        call error_mesg('xgrid_mod', 'cannot find grid file field', FATAL)

      rcode = nf_inq_varid(ncid, 'AREA_'//grid_ids(1)//'x'//grid_ids(g), area_id)
      if (rcode/=0) &
        call error_mesg('xgrid_mod', 'cannot find grid file field', FATAL)
      if(use_higher_order) then
         rcode = nf_inq_varid(ncid, 'DI_'//grid_ids(1)//'x'//grid_ids(g), di_id)
         if (rcode/=0) &
           call error_mesg('xgrid_mod', 'cannot find grid file field DI1_'//grid_ids(1)//'x'//grid_ids(g), FATAL)

         rcode = nf_inq_varid(ncid, 'DJ_'//grid_ids(1)//'x'//grid_ids(g), dj_id)
         if (rcode/=0) &
           call error_mesg('xgrid_mod', 'cannot find grid file field DJ1'//grid_ids(1)//'x'//grid_ids(g), FATAL)
      endif
      rcode = nf_inq_vardimid(ncid, area_id, dims)
      rcode = nf_inq_dimlen(ncid, dims(1), n_areas)

      ! load exchange cells, sum grid cell areas, set your1my2/your2my1
      if(use_higher_order) then
         call load_xgrid (xmap, grid, grid%domain, ncid, i1_id, j1_id, i2_id, j2_id, &
                       area_id, n_areas, di_id, dj_id)
      else 
         call load_xgrid (xmap, grid, grid%domain, ncid, i1_id, j1_id, i2_id, j2_id, &
                       area_id, n_areas )
      endif
    end if

      ! get the center point of the grid box

      allocate(grid%lon(grid%im,grid%jm), grid%lat(grid%im,grid%jm))
      call get_grid(grid, grid_ids(g), ncid)

  end do
  rcode = nf_close(ncid)

  grid1%area_inv = 0.0;
  where (grid1%area>0.0)
    grid1%area_inv = 1.0/grid1%area
  end where

  xmap%your1my2(xmap%me-xmap%root_pe) = .false. ! this is not necessarily true but keeps
  xmap%your2my1(xmap%me-xmap%root_pe) = .false. ! a PE from communicating with itself

  send_size = grid1%im*grid1%jm
  if (make_exchange_reproduce) then
    allocate( xmap%send_count_repro(0:xmap%npes-1) )
    allocate( xmap%recv_count_repro(0:xmap%npes-1) )
    xmap%send_count_repro = 0
    xmap%recv_count_repro = 0
    do g=2,size(xmap%grids)
      do p=0,xmap%npes-1
        xmap%send_count_repro(p) = xmap%send_count_repro(p) &
                                  +count(xmap%grids(g)%x      (:)%pe==p+xmap%root_pe)
        xmap%recv_count_repro(p) = xmap%recv_count_repro(p) &
                                  +count(xmap%grids(g)%x_repro(:)%pe==p+xmap%root_pe)
      end do
    end do
    send_size = max(send_size, sum(xmap%send_count_repro))
  end if
  allocate (xmap%send_buffer(send_size))

  call mpp_open( unit, 'xgrid.out', action= MPP_OVERWR,threading= MPP_MULTI )  

  write( unit,* )xmap%grids(:)%id, ' GRID: PE ', xmap%me, ' #XCELLS=', &
           xmap%grids(2:size(xmap%grids))%size, ' #COMM. PARTNERS=', &
           count(xmap%your1my2), '/', count(xmap%your2my1), &
           pack((/(p+xmap%root_pe,p=0,xmap%npes-1)/), xmap%your1my2),  &
           '/', pack((/(p+xmap%root_pe,p=0,xmap%npes-1)/), xmap%your2my1)

  allocate( xmap%x1(1:sum(xmap%grids(2:size(xmap%grids))%size)) )
  allocate( xmap%x2(1:sum(xmap%grids(2:size(xmap%grids))%size)) )

  call regen(xmap)

  if ( debug ) then
     xxx = conservation_check(grid1%area*0+1.0, grid1%id, xmap)
     if (xmap%me==xmap%root_pe)write( *,* )"Checked data is array of constant 1"
     if (xmap%me==xmap%root_pe)write( *,* )grid1%id,'(',xmap%grids(:)%id,')=', xxx 

     do g=2,size(xmap%grids)
       xxx = conservation_check(xmap%grids(g)%frac_area*0+1.0, xmap%grids(g)%id, &
                                                                           xmap )
       if (xmap%me==xmap%root_pe) &
         write( *,* )xmap%grids(g)%id,'(',xmap%grids(:)%id,')=', xxx 
     end do

   ! create an random number 2d array
     if(grid1%id == "ATM") then
        allocate(check_data(size(grid1%area,1), size(grid1%area,2)))
        call random_number(check_data)
        xxx = conservation_check(check_data, grid1%id, xmap,  remapping_order =2)
        deallocate(check_data)
        if (xmap%me==xmap%root_pe)write( *,* )"Checked data is array of random number between 0 and 1"
        if (xmap%me==xmap%root_pe)write( *,* )grid1%id,'(',xmap%grids(:)%id,')=', xxx  

        do g=2,size(xmap%grids)
          allocate(check_data_3d(size(xmap%grids(g)%frac_area,1),size(xmap%grids(g)%frac_area,2), &
                                size(xmap%grids(g)%frac_area,3) )) 
          call random_number(check_data_3d)
          xxx = conservation_check(check_data_3d, xmap%grids(g)%id, xmap,  remapping_order =2 )
          if (xmap%me==xmap%root_pe) write( *,* )xmap%grids(g)%id,'(',xmap%grids(:)%id,')=', xxx
          deallocate( check_data_3d)
        end do
    endif
  endif
  call close_file (unit)

end subroutine setup_xmap
! </SUBROUTINE>

!#######################################################################


subroutine regen(xmap)
type (xmap_type), intent(inout) :: xmap

  integer :: g, l, i, j, k, max_size

  max_size = 0;
  do g=2,size(xmap%grids)
    max_size = max_size + xmap%grids(g)%size * xmap%grids(g)%km
  end do
  if (max_size>size(xmap%x1)) then
    deallocate(xmap%x1)
    deallocate(xmap%x2)
    allocate( xmap%x1(1:max_size) )
    allocate( xmap%x2(1:max_size) )
  end if

  xmap%size = 0
  do g=2,size(xmap%grids)
    xmap%grids(g)%first = xmap%size + 1;
    do l=1,xmap%grids(g)%size
      i = xmap%grids(g)%x(l)%i2
      j = xmap%grids(g)%x(l)%j2
      do k=1,xmap%grids(g)%km
        if (xmap%grids(g)%frac_area(i,j,k)/=0.0) then
          xmap%size = xmap%size+1
          xmap%x1(xmap%size)%i    = xmap%grids(g)%x(l)%i1
          xmap%x1(xmap%size)%j    = xmap%grids(g)%x(l)%j1
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
  type(grid_type), pointer :: grid

  if (grid_id==xmap%grids(1)%id) call error_mesg ('xgrid_mod',  &
                                   'set_frac_area called on side 1 grid', FATAL)
  do g=2,size(xmap%grids)
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
!   <IN NAME="remapping_order" TYPE="integer,optional"></IN>

subroutine put_side1_to_xgrid(d, grid_id, x, xmap, remapping_order)
real, dimension(:,:), intent(in   ) :: d
character(len=3),     intent(in   ) :: grid_id
real, dimension(:),   intent(inout) :: x
type (xmap_type),     intent(inout) :: xmap
integer, intent(in), optional       :: remapping_order

  integer :: g
  integer :: order 

  order = 1
  if(present(remapping_order)) order = remapping_order

  if ((order .ne. 1) .and. (order .ne. 2)) &
    call error_mesg ('xgrid_mod', 'xgrid only support first and second order interpolation', FATAL) 

  if (grid_id==xmap%grids(1)%id) then
       call put_1_to_xgrid(d, x, xmap, order)
    return;
  end if

  do g=2,size(xmap%grids)
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

  do g=2,size(xmap%grids)
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
  
  do g=2,size(xmap%grids)
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
  
  do g=2,size(xmap%grids)
    if (grid_id==xmap%grids(g)%id) then
      call get_2_from_xgrid(d, xmap%grids(g), x, xmap)
      return;
    end if
  end do
  
  call error_mesg ('xgrid_mod', 'get_from_xgrid: could not find grid id', FATAL)

end subroutine get_side2_from_xgrid
! </SUBROUTINE>

!#######################################################################

! <FUNCTION NAME="some">

!   <OVERVIEW>
!     Returns logical associating exchange grid cells with given side two grid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns logical associating exchange grid cells with given side two grid.
!   </DESCRIPTION>
!   <TEMPLATE>
!     some(xmap, grid_id)
!   </TEMPLATE>

!   <IN NAME="xmap"  TYPE="xmap_type"  ></IN>
!   <IN NAME="grid_id"  TYPE="character(len=3)"  ></IN>
!   <OUT NAME="some"  TYPE="logical" DIM="(xmap%size)" >
!     logical associating exchange grid cells with given side 2 grid.
!   </OUT>

function some(xmap, grid_id)
type (xmap_type),           intent(in) :: xmap
character(len=3), optional, intent(in) :: grid_id
logical, dimension(xmap%size)          :: some

  integer :: g, l

  if (.not.present(grid_id)) then
    some = .true.
    return;
  end if

  if (grid_id==xmap%grids(1)%id) &
    call error_mesg ('xgrid_mod', 'some expects a side 2 grid id', FATAL)
  
  do g=2,size(xmap%grids)
    if (grid_id==xmap%grids(g)%id) then
      some = .false.
      some(xmap%grids(g)%first:xmap%grids(g)%last) = .true.;
      return;
    end if
  end do
  
  call error_mesg ('xgrid_mod', 'some could not find grid id', FATAL)

end function some
! </FUNCTION>

!#######################################################################

subroutine put_2_to_xgrid(d, grid, x, xmap)
type (grid_type),                                intent(in) :: grid
real, dimension(grid%is_me:grid%ie_me, &
                grid%js_me:grid%je_me, grid%km), intent(in) :: d
real, dimension(:    ), intent(inout) :: x
type (xmap_type),       intent(in   ) :: xmap

  integer                 :: k, l
  type (x2_type), pointer :: c

  do l=grid%first,grid%last
    c => xmap%x2(l)
    x(l) = d(c%i,c%j,c%k)
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
  type (x2_type), pointer :: c

  d = 0.0
  do l=grid%first,grid%last
    c => xmap%x2(l)
    d(c%i,c%j,c%k) = d(c%i,c%j,c%k) + c%area*x(l)
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

  real, dimension(im*jm) :: buf
  integer :: i, j, l

!  call mpp_recv(buf, im*jm, pe)
!  l = 0
!  do j=1,jm; do i=1,im;
!    l = l + 1
!    get_side_1(i,j) = buf(l)
!  end do; end do
  call mpp_recv( get_side_1, im*jm, pe )
end function get_side_1

!#######################################################################

subroutine put_1_to_xgrid(d, x, xmap, remapping_order)
real, dimension(:,:), intent(in   ) :: d
real, dimension(:  ), intent(inout) :: x
type (xmap_type),     intent(inout) :: xmap
integer,              intent(in   ) :: remapping_order

  integer :: i, is, ie, im, j, js, je, jm, p, l
  type (x1_type), pointer :: c
  real, dimension(xmap%grids(1)%im,xmap%grids(1)%jm) :: dg
  type (grid_type), pointer :: grid1

  if(remapping_order .eq. 2) then
     call put_1_to_xgrid_order_2(d, x, xmap)
     return
  endif

  grid1 => xmap%grids(1)
  is = grid1%is_me; ie = grid1%ie_me;
  js = grid1%js_me; je = grid1%je_me;
  dg(is:ie,js:je) = d;

  im = ie-is+1; jm = je-js+1;
  l = 0
  call mpp_sync_self()          !Balaji
  do j=1,jm; do i=1,im;
    l = l + 1;
    xmap%send_buffer(l) =  d(i,j)
  end do; end do;
  do p=0,xmap%npes-1
    if (xmap%your2my1(p)) then
      call mpp_send(xmap%send_buffer, im*jm, p+xmap%root_pe);
    end if
  end do
  do p=0,xmap%npes-1
    if (xmap%your1my2(p)) then
      is = grid1%is(p); ie = grid1%ie(p);
      js = grid1%js(p); je = grid1%je(p);
      dg(is:ie,js:je) = get_side_1(p+xmap%root_pe,ie-is+1,je-js+1);
    end if
  end do
  do l=1,xmap%size
    c    => xmap%x1(l)
    x(l) =  dg(c%i,c%j)
  end do

!  call mpp_sync_self
end subroutine put_1_to_xgrid

!#######################################################################


subroutine put_1_to_xgrid_order_2(d, x, xmap)
real, dimension(:,:), intent(in   ) :: d
real, dimension(:  ), intent(inout) :: x
type (xmap_type),     intent(inout) :: xmap

  integer :: i, is, ie, im, j, js, je, jm, p, l, isd, ied, jsd, jed
  type (x1_type), pointer :: c
  real, dimension(xmap%grids(1)%im,xmap%grids(1)%jm) :: dg
  real, dimension(xmap%grids(1)%im,xmap%grids(1)%jm,2) :: grad_d
  real, dimension(:,:), allocatable :: d_halo  ! store the data with halo
  type (grid_type), pointer :: grid1
  real           :: dlon, dlat

  grid1 => xmap%grids(1)
  is = grid1%is_me; ie = grid1%ie_me
  js = grid1%js_me; je = grid1%je_me
  dg(is:ie,js:je) = d

! first get the halo of data
  call mpp_get_compute_domain( grid1%domain_with_halo, is,  ie,  js,  je  )
  call mpp_get_data_domain   ( grid1%domain_with_halo, isd, ied, jsd, jed )
  allocate( d_halo(isd:ied,jsd:jed))
  d_halo = 0.0
  d_halo(is:ie,js:je) = d(:,:)
  call mpp_update_domains(d_halo,grid1%domain_with_halo)      
  grad_d(is:ie,js:je,:) = grad(d_halo, grid1, xmap,xmap%me-xmap%root_pe)

  im = ie-is+1; jm = je-js+1
  l = 0
  call mpp_sync_self()          !Balaji

! if size of send_buffer is not enough, need to reallocate send_buffer
    if(size(xmap%send_buffer) .lt. size(d_halo)) then
       deallocate(xmap%send_buffer)
       allocate(xmap%send_buffer(size(d_halo)))
     endif

    do j=js-1,je+1; do i=is-1,ie+1
      l = l + 1
      xmap%send_buffer(l) =  d_halo(i,j)
    end do; end do

    do p=0,xmap%npes-1
      if (xmap%your2my1(p)) then
        call mpp_send(xmap%send_buffer, (im+2)*(jm+2), p+xmap%root_pe);
      end if
    end do

  do p=0,xmap%npes-1
    if (xmap%your1my2(p)) then
      is = grid1%is(p); ie = grid1%ie(p)
      js = grid1%js(p); je = grid1%je(p)
      isd = is-1; ied = ie + 1
      jsd = js-1; jed = je + 1
      if(allocated(d_halo)) deallocate(d_halo)
      allocate(d_halo(isd:ied,jsd:jed))
      d_halo(isd:ied,jsd:jed) = get_side_1(p+xmap%root_pe,ie-is+3,je-js+3)
      dg(is:ie,js:je) = d_halo(is:ie,js:je)
      grad_d(is:ie,js:je,:) = grad(d_halo, grid1, xmap,p)
    end if
  end do

  do l=1,xmap%size
    c    => xmap%x1(l)
    x(l) =  dg(c%i,c%j)
    x(l) = x(l) + grad_d(c%i,c%j,1)*c%di + grad_d(c%i,c%j,2)*c%dj
  end do


! free the memory of the local data
  deallocate (d_halo)

end subroutine put_1_to_xgrid_order_2

!#######################################################################

subroutine get_1_from_xgrid(d, x, xmap)
real, dimension(:,:), intent(out)   :: d
real, dimension(:  ), intent(in )   :: x
type (xmap_type),     intent(inout) :: xmap

  real, dimension(xmap%grids(1)%im,xmap%grids(1)%jm), target :: dg
  integer :: i, is, ie, im, j, js, je, jm, l, le, p
  real             , pointer :: dgp
  type (x1_type)   , pointer :: c
  type (grid_type) , pointer :: grid1

  grid1 => xmap%grids(1)

  dg = 0.0;
  do l=1,xmap%size
    c   => xmap%x1(l)
    dgp => dg(c%i,c%j)
    dgp =  dgp + c%area*x(l)
  end do

  le = 0;
  call mpp_sync_self()          !Balaji
  do p=0,xmap%npes-1
    if (xmap%your1my2(p)) then
      l = le + 1;
      is = grid1%is(p); ie = grid1%ie(p);
      js = grid1%js(p); je = grid1%je(p);
      do j=js,je; do i=is,ie;
        le = le + 1
        xmap%send_buffer(le) = dg(i,j)
      end do; end do;
      call mpp_send(xmap%send_buffer(l:le), le-l+1, p+xmap%root_pe);
    end if
  end do
  d = dg(grid1%is_me:grid1%ie_me,grid1%js_me:grid1%je_me);
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
  type (grid_type), pointer :: grid

  allocate ( x_psum  (sum(xmap%send_count_repro)) )
  allocate ( pe_psum (sum(xmap%send_count_repro)) )
  x_psum = 0.0
  l1 = 0 ! index into partition summed exchange grid variable
  l2 = 0 ! index into exchange grid variable
  do g=2,size(xmap%grids)
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
      call mpp_send(xmap%send_buffer(l1:l2), l2-l1+1, p+xmap%root_pe);
    end if
  end do
  deallocate ( x_psum, pe_psum)
  allocate ( x_psum (sum(xmap%recv_count_repro)) )
  l2 = 0;
  do p=0,xmap%npes-1
    l1 = l2 + 1
    l2 = l2 + xmap%recv_count_repro(p)
    if (xmap%recv_count_repro(p)>0) then ! can receive from myself
      call mpp_recv(x_psum(l1:l2), l2-l1+1, p+xmap%root_pe);
      pl(p) = l1
    end if
  end do
  d = 0.0
  do g=2,size(xmap%grids)
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
!   <IN NAME="remapping_order" TYPE="integer,optional"></IN>
! conservation_check - returns three numbers which are the global sum of a
! variable (1) on its home model grid, (2) after interpolation to the other
! side grid(s), and (3) after re_interpolation back onto its home side grid(s).
!
function conservation_check_side1(d, grid_id, xmap,remapping_order) ! this one for 1->2->1
real, dimension(:,:), intent(in   ) :: d
character(len=3),     intent(in   ) :: grid_id
type (xmap_type),     intent(inout) :: xmap
real, dimension(3)                  :: conservation_check_side1
integer,              intent(in), optional :: remapping_order

  real                       :: gsum
  real, dimension(xmap%size) :: x_over, x_back
  real, dimension(size(d,1),size(d,2)) :: d1
  real, dimension(:,:,:), allocatable  :: d2
  integer                              :: g
  type (grid_type), pointer            :: grid1, grid2
  integer                              :: order

  order = 1
  if(present(remapping_order)) order = remapping_order
  grid1 => xmap%grids(1)
  conservation_check_side1(1) = mpp_global_sum(grid1%domain, grid1%area*d)
  conservation_check_side1(2) = 0.0
  call put_to_xgrid (d, grid1%id, x_over, xmap, remapping_order = order)    ! put from side 1
  do g=2,size(xmap%grids)
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

function conservation_check_side2(d, grid_id, xmap,remapping_order) ! this one for 2->1->2
real, dimension(:,:,:), intent(in   ) :: d
character(len=3),       intent(in   ) :: grid_id
type (xmap_type),       intent(inout) :: xmap
real, dimension(3)                    :: conservation_check_side2
integer,              intent(in), optional :: remapping_order

  real                       :: gsum
  real, dimension(xmap%size) :: x_over, x_back
  real, dimension(:,:  ), allocatable :: d1
  real, dimension(:,:,:), allocatable :: d2
  integer                             :: g
  type (grid_type), pointer           :: grid1, grid2
  integer                              :: order

  order = 1
  if(present(remapping_order)) order = remapping_order
  grid1 => xmap%grids(1)
  do g = 2,size(xmap%grids)
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
  call put_to_xgrid(d1,  grid1%id, x_back, xmap,remapping_order = order)   ! put from side 1
  deallocate ( d1 )

  conservation_check_side2(3) = 0.0;
  do g = 2,size(xmap%grids)
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

! This function is used to calculate the gradient of a given data.
! Maybe need to setup a limit for the gradient. 

function grad(d, grid, xmap,p) 

type(grid_type),      intent(in) :: grid   
real, dimension(:,:), intent(in) :: d
type(xmap_type),      intent(in) :: xmap
integer, intent(in)              :: p

real, dimension(grid%is(p):grid%ie(p),grid%js(p):grid%je(p),2) :: grad
real, dimension(:,:), allocatable :: d_local, di, dj
integer :: im, jm, i, j, k, is, ie, js, je
real    :: dx1, dx2, dy1, dy2, costheta
real    :: EPS = 1.0e-10
real    :: data(4), lon(4), lat(4)
integer ::  ilon(4), jlat(4)

   is = grid%is(p); ie = grid%ie(p)
   js = grid%js(p); je = grid%je(p)

   im = ie - is + 1
   jm = je - js + 1

   if(im /= size(d,1)-2 .or. jm /= size(d,2)-2) call error_mesg('xgrids_mod', 'The size of data is not &
                                                   &consistent with the sie of model grid', FATAL) 
   allocate(d_local(is-1:ie+1,js-1:je+1), di(is:ie,js:je), dj(is:ie,js:je))
   d_local(is-1:ie+1,js-1:je+1) = d

!  calculate the gradient of the data on each grid
   do i = is, ie
   do j = js, je
      costheta = cos(grid%lat(i,j))
      if(abs(costheta) .lt. EPS)  &
            call error_mesg('xgrids_mod', 'Impropere latitude grid', FATAL)

      if(i == grid%im) then
         if(j == 1) then
            ilon(1) = i-1; ilon(2) = i;   ilon(3) = i;   ilon(4) = i-1
            jlat(1) = j;   jlat(2) = j;   jlat(3) = j+1; jlat(4) = j+1
         else if( j == grid%jm) then
            ilon(1) = i-1; ilon(2) = i;   ilon(3) = i;   ilon(4) = i-1
            jlat(1) = j-1; jlat(2) = j-1; jlat(3) = j;   jlat(4) = j    
         else
            ilon(1) = i-1; ilon(2) = i;   ilon(3) = i;   ilon(4) = i-1
            jlat(1) = j-1; jlat(2) = j-1; jlat(3) = j+1; jlat(4) = j+1
         endif
      else if(i==1) then
         if(j == 1) then
            ilon(1) = i;   ilon(2) = i+1; ilon(3) = i+1; ilon(4) = i
            jlat(1) = j;   jlat(2) = j;   jlat(3) = j+1; jlat(4) = j+1
         else if( j == grid%jm) then
            ilon(1) = i;   ilon(2) = i+1; ilon(3) = i+1; ilon(4) = i
            jlat(1) = j-1; jlat(2) = j-1; jlat(3) = j;   jlat(4) = j    
         else
            ilon(1) = i;   ilon(2) = i+1; ilon(3) = i+1; ilon(4) = i
            jlat(1) = j-1; jlat(2) = j-1; jlat(3) = j+1; jlat(4) = j+1
         endif
      else 
         if(j == 1) then
            ilon(1) = i-1; ilon(2) = i+1; ilon(3) = i+1; ilon(4) = i-1
            jlat(1) = j;   jlat(2) = j;   jlat(3) = j+1; jlat(4) = j+1
         else if( j == grid%jm) then
            ilon(1) = i-1; ilon(2) = i+1; ilon(3) = i+1; ilon(4) = i-1
            jlat(1) = j-1; jlat(2) = j-1; jlat(3) = j;   jlat(4) = j    
         else
            ilon(1) = i-1; ilon(2) = i+1; ilon(3) = i+1; ilon(4) = i-1
            jlat(1) = j-1; jlat(2) = j-1; jlat(3) = j+1; jlat(4) = j+1
         endif
      endif
         do k = 1, 4
            lon(k) = grid%lon(ilon(k),jlat(k))
            lat(k) = grid%lat(ilon(k),jlat(k))            
            data(k)= d_local(ilon(k),jlat(k)) 
         enddo
         dx1 = lon(2) - lon(1)
         dx2 = lon(3) - lon(4)
         dy1 = lat(4) - lat(1)
         dy2 = lat(3) - lat(2)

         if(abs(dx1).lt.EPS .or. abs(dx2) .lt. EPS)  &
            call error_mesg('xgrids_mod', 'Improper grid size in lontitude', FATAL)
         if(dx1 .gt. PI) dx1 = dx1 - 2.0* PI
         if(dx1 .lt. -PI) dx1 = dx1 + 2.0* PI
         if(dx2 .gt. PI) dx2 = dx2 - 2.0* PI
         if(dx2 .lt. -PI) dx2 = dx2 + 2.0* PI
         di(i,j) = ((data(2)-data(1))/dx1 + (data(3) -data(4))/dx2)/(2.0*costheta)

         if(abs(dy1).lt.EPS .or. abs(dy2).lt.EPS )  &
            call error_mesg('xgrids_mod', 'Improper grid size in latitude', FATAL)
         dj(i,j) = ((data(4)-data(1))/dy1 + (data(3) -data(2))/dy2)/2.0

   enddo
   enddo

   grad(:,:,1) = di
   grad(:,:,2) = dj

   deallocate(d_local)
   deallocate(di, dj)

   return

end function grad


!#######################################################################


end module xgrid_mod

! <INFO>

!   <REFERENCE>   
!      A <LINK SRC="http://www.gfdl.noaa.gov/~mw/docs/grid_coupling.html"> guide </LINK>to grid coupling in FMS.
!   </REFERENCE>
!   <REFERENCE>
!      A simple xgrid <LINK SRC="http://www.gfdl.gov/~mw/docs/xgrid_example.f90.txt"> example. </LINK>
!   </REFERENCE>

! </INFO>
