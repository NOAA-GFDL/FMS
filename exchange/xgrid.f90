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

use       fms_mod,   only: file_exist, open_namelist_file, check_nml_error,  &
                           error_mesg, close_file, FATAL, stdlog,            &
                           write_version_number 
use mpp_mod,         only: mpp_npes, mpp_pe, mpp_root_pe, mpp_send, mpp_recv, &
                           mpp_sync_self
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                           Domain2d, mpp_global_sum
use mpp_io_mod,      only: mpp_open, MPP_MULTI, MPP_OVERWR

implicit none
include 'netcdf.inc'
private

public xmap_type, setup_xmap, set_frac_area, put_to_xgrid, get_from_xgrid, &
       xgrid_count, some, conservation_check, xgrid_init

logical :: make_exchange_reproduce = .false. ! exactly same on different # PEs
namelist /xgrid_nml/ make_exchange_reproduce
logical :: init = .true.

interface put_to_xgrid
  module procedure put_side1_to_xgrid
  module procedure put_side2_to_xgrid
end interface

interface get_from_xgrid
  module procedure get_side1_from_xgrid
  module procedure get_side2_from_xgrid
end interface

interface conservation_check
  module procedure conservation_check_side1
  module procedure conservation_check_side2
end interface

type xcell_type
  integer :: i1, j1, i2, j2 ! indices of cell in model arrays on both sides
  integer :: pe             ! other side pe that has this cell
  real    :: area           ! geographic area of exchange cell
end type xcell_type

type grid_type
  character(len=3)                :: id           ! grid identifier
  integer, pointer, dimension(:)  :: is, ie       ! domain - i-range (pe index)
  integer, pointer, dimension(:)  :: js, je       ! domain - j-range (pe index)
  integer, pointer                :: is_me, ie_me ! my domain - i-range
  integer, pointer                :: js_me, je_me ! my domain - j-range
  integer                         :: im, jm, km   ! global domain range
  real, pointer, dimension(:,:,:) :: frac_area    ! partition fractions
  real, pointer, dimension(:,:)   :: area         ! cell area
  real, pointer, dimension(:,:)   :: area_inv     ! 1 / area for normalization
  integer                         :: first, last  ! xgrid index range
  integer                         :: size         ! # xcell patterns
  type(xcell_type), pointer, dimension(:) :: x    ! xcell patterns
  integer                         :: size_repro   ! # side 1 patterns for repro
  type(xcell_type), pointer, dimension(:) :: x_repro ! side 1 patterns for repro
  type(Domain2d) :: domain                        ! used for conservation checks
end type grid_type

type x1_type
  integer :: i, j
  real    :: area ! (= geographic area * frac_area)
end type x1_type

type x2_type
  integer :: i, j, k
  real    :: area ! geographic area of exchange cell
end type x2_type

type xmap_type
  private
  integer :: size            ! # of exchange grid cells with area > 0 on this pe

  integer :: me, npes
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
 character(len=128) :: version = '$Id: xgrid.f90,v 1.3 2002/07/16 22:55:10 fms Exp $'
 character(len=128) :: tagname = '$Name: havana $'

 logical :: module_is_initialized = .FALSE.

contains

!#######################################################################

logical function in_box(i, j, is, ie, js, je)
integer :: i, j, is, ie, js, je

  in_box = (i>=is) .and. (i<=ie) .and. (j>=js) .and. (j<=je)
end function in_box

!#######################################################################

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

!#######################################################################

subroutine load_xgrid (xmap, grid, domain, ncid, id_i1, id_j1, id_i2, id_j2, &
                                                           id_area, n_areas  )
type(xmap_type), intent(inout)  :: xmap
type(grid_type), intent(inout)  :: grid
type(Domain2d), intent(inout)   :: domain
integer, intent(in) :: ncid, id_i1, id_j1, id_i2, id_j2, id_area, n_areas

  integer, dimension(n_areas) :: i1, j1, i2, j2 ! xgrid quintuples
  real,    dimension(n_areas) :: area           ! from grid file
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
      if (make_exchange_reproduce) then
        do p=0,xmap%npes-1
          if (in_box(i1(l), j1(l), grid1%is(p), grid1%ie(p), &
                                   grid1%js(p), grid1%je(p))) then
            grid%x(ll)%pe = p
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
      do p=0,xmap%npes-1
        if (in_box(i2(l), j2(l), grid%is(p), grid%ie(p), &
                                 grid%js(p), grid%je(p))) then
          grid%x_repro(ll_repro)%pe = p
        end if
      end do
    end if ! make_exchange_reproduce
  end do
  grid%area_inv = 0.0;
  where (grid%area>0.0) grid%area_inv = 1.0/grid%area
end subroutine load_xgrid

!#######################################################################
!
! setup_xmap - sets up exchange grid connectivity using grid specification file
!              and processor domain decomposition.  initializes xmap.
!

subroutine setup_xmap(xmap, grid_ids, grid_domains, grid_file)
type (xmap_type),                          intent(out) :: xmap
character(len=3), dimension(:),            intent(in ) :: grid_ids
type(Domain2d), dimension(size(grid_ids)), intent(in ) :: grid_domains
character(len=*)                         , intent(in ) :: grid_file

  integer :: g, l, ll, p, n_areas, send_size
  integer :: ncid, i1_id, j1_id, i2_id, j2_id, area_id
  integer :: dims(4), rcode
  integer :: unit
  type (grid_type), pointer :: grid, grid1
  real, dimension(3) :: xxx

  xmap%me   = mpp_pe  ()
  xmap%npes = mpp_npes()

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

    allocate ( grid%is(0:xmap%npes-1), grid%ie(0:xmap%npes-1) )
    allocate ( grid%js(0:xmap%npes-1), grid%je(0:xmap%npes-1) )
    call mpp_get_compute_domains(grid%domain, xbegin=grid%is, xend=grid%ie, &
                                              ybegin=grid%js, yend=grid%je  )
    grid%is_me => grid%is(xmap%me); grid%ie_me => grid%ie(xmap%me)
    grid%js_me => grid%js(xmap%me); grid%je_me => grid%je(xmap%me)
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
      rcode = nf_inq_vardimid(ncid, area_id, dims)
      rcode = nf_inq_dimlen(ncid, dims(1), n_areas)

      ! load exchange cells, sum grid cell areas, set your1my2/your2my1
      call load_xgrid (xmap, grid, grid%domain, ncid, &
                       i1_id, j1_id, i2_id, j2_id, area_id, n_areas)
    end if
  end do
  rcode = nf_close(ncid)

  grid1%area_inv = 0.0;
  where (grid1%area>0.0)
    grid1%area_inv = 1.0/grid1%area
  end where

  xmap%your1my2(xmap%me) = .false. ! this is not necessarily true but keeps
  xmap%your2my1(xmap%me) = .false. ! a PE from communicating with itself

  send_size = grid1%im*grid1%jm
  if (make_exchange_reproduce) then
    allocate( xmap%send_count_repro(0:xmap%npes-1) )
    allocate( xmap%recv_count_repro(0:xmap%npes-1) )
    xmap%send_count_repro = 0
    xmap%recv_count_repro = 0
    do g=2,size(xmap%grids)
      do p=0,xmap%npes-1
        xmap%send_count_repro(p) = xmap%send_count_repro(p) &
                                  +count(xmap%grids(g)%x      (:)%pe==p)
        xmap%recv_count_repro(p) = xmap%recv_count_repro(p) &
                                  +count(xmap%grids(g)%x_repro(:)%pe==p)
      end do
    end do
    send_size = max(send_size, sum(xmap%send_count_repro))
  end if
  allocate (xmap%send_buffer(send_size))

  call mpp_open( unit, 'xgrid.out', action= MPP_OVERWR,threading= MPP_MULTI )  

  write( unit,* )xmap%grids(:)%id, ' GRID: PE ', xmap%me, ' #XCELLS=', &
           xmap%grids(2:size(xmap%grids))%size, ' #COMM. PARTNERS=', &
           count(xmap%your1my2), '/', count(xmap%your2my1), &
           pack((/(p,p=0,xmap%npes-1)/), xmap%your1my2),  &
           '/', pack((/(p,p=0,xmap%npes-1)/), xmap%your2my1)

  allocate( xmap%x1(1:sum(xmap%grids(2:size(xmap%grids))%size)) )
  allocate( xmap%x2(1:sum(xmap%grids(2:size(xmap%grids))%size)) )

  call regen(xmap)

  xxx = conservation_check(grid1%area*0+1.0, grid1%id, xmap)
  if (xmap%me==mpp_root_pe())write( unit,* )grid1%id,'(',xmap%grids(:)%id,')=', xxx 
  do g=2,size(xmap%grids)
    xxx = conservation_check(xmap%grids(g)%frac_area*0+1.0, xmap%grids(g)%id, &
                                                                        xmap  )
    if (xmap%me==mpp_root_pe()) &
      write( unit,* )xmap%grids(g)%id,'(',xmap%grids(:)%id,')=', xxx 
  end do
  call close_file (unit)
end subroutine setup_xmap

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
!
! set_frac_area - changes sub-grid partion areas and/or number.
!
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

!#######################################################################
!
! xgrid_count - returns current size of exchange grid variables.
!
integer function xgrid_count(xmap)
type (xmap_type), intent(inout) :: xmap

  xgrid_count = xmap%size
end function xgrid_count

!#######################################################################

subroutine put_side1_to_xgrid(d, grid_id, x, xmap)
real, dimension(:,:), intent(in   ) :: d
character(len=3),     intent(in   ) :: grid_id
real, dimension(:),   intent(inout) :: x
type (xmap_type),     intent(inout) :: xmap

  integer :: g

  if (grid_id==xmap%grids(1)%id) then
    call put_1_to_xgrid(d, x, xmap)
    return;
  end if

  do g=2,size(xmap%grids)
    if (grid_id==xmap%grids(g)%id)    &
      call error_mesg ('xgrid_mod',  &
                       'put_to_xgrid expects a 3D side 2 grid', FATAL)
  end do

  call error_mesg ('xgrid_mod', 'put_to_xgrid: could not find grid id', FATAL)

end subroutine put_side1_to_xgrid

!#######################################################################

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

!#######################################################################

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

!#######################################################################

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

!#######################################################################
!
! some - returns logical associating exchange grid cells with given side 2 grid
!
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

!#######################################################################

subroutine put_2_to_xgrid(d, grid, x, xmap)
type (grid_type),                                intent(in) :: grid
real, dimension(grid%is_me:grid%ie_me, &
                grid%js_me:grid%je_me, grid%km), intent(in) :: d
real, dimension(:    ), intent(inout) :: x
type (xmap_type),       intent(in   ) :: xmap

  integer                 :: i, j, k, l
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

subroutine put_1_to_xgrid(d, x, xmap)
real, dimension(:,:), intent(in   ) :: d
real, dimension(:  ), intent(inout) :: x
type (xmap_type),     intent(inout) :: xmap

  integer :: i, is, ie, im, j, js, je, jm, p, l
  type (x1_type), pointer :: c
  real, dimension(xmap%grids(1)%im,xmap%grids(1)%jm) :: dg
  type (grid_type), pointer :: grid1

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
      call mpp_send(xmap%send_buffer, im*jm, p);
    end if
  end do
  do p=0,xmap%npes-1
    if (xmap%your1my2(p)) then
      is = grid1%is(p); ie = grid1%ie(p);
      js = grid1%js(p); je = grid1%je(p);
      dg(is:ie,js:je) = get_side_1(p,ie-is+1,je-js+1);
    end if
  end do
  do l=1,xmap%size
    c    => xmap%x1(l)
    x(l) =  dg(c%i,c%j)
  end do

!  call mpp_sync_self
end subroutine put_1_to_xgrid

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
      call mpp_send(xmap%send_buffer(l:le), le-l+1, p);
    end if
  end do
  d = dg(grid1%is_me:grid1%ie_me,grid1%js_me:grid1%je_me);
  im = grid1%ie_me-grid1%is_me+1;
  jm = grid1%je_me-grid1%js_me+1;
  do p=0,xmap%npes-1
    if (xmap%your2my1(p)) d = d + get_side_1(p,im,jm)
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
      xmap%send_buffer(l1:l2) = pack(x_psum, pe_psum==p)
      call mpp_send(xmap%send_buffer(l1:l2), l2-l1+1, p);
    end if
  end do
  deallocate ( x_psum, pe_psum)
  allocate ( x_psum (sum(xmap%recv_count_repro)) )
  l2 = 0;
  do p=0,xmap%npes-1
    l1 = l2 + 1
    l2 = l2 + xmap%recv_count_repro(p)
    if (xmap%recv_count_repro(p)>0) then ! can receive from myself
      call mpp_recv(x_psum(l1:l2), l2-l1+1, p);
      pl(p) = l1
    end if
  end do
  d = 0.0
  do g=2,size(xmap%grids)
    grid => xmap%grids(g)
    do l3=1,grid%size_repro ! index into side1 grid's patterns
      i = grid%x_repro(l3)%i1
      j = grid%x_repro(l3)%j1
      d(i,j) = d(i,j) + x_psum(pl(grid%x_repro(l3)%pe))
      pl(grid%x_repro(l3)%pe) = pl(grid%x_repro(l3)%pe) + 1
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
!
! conservation_check - returns three numbers which are the global sum of a
! variable (1) on its home model grid, (2) after interpolation to the other
! side grid(s), and (3) after re_interpolation back onto its home side grid(s).
!
function conservation_check_side1(d, grid_id, xmap) ! this one for 1->2->1
real, dimension(:,:), intent(in   ) :: d
character(len=3),     intent(in   ) :: grid_id
type (xmap_type),     intent(inout) :: xmap
real, dimension(3)                  :: conservation_check_side1

  real                       :: gsum
  real, dimension(xmap%size) :: x_over, x_back
  real, dimension(size(d,1),size(d,2)) :: d1
  real, dimension(:,:,:), allocatable  :: d2
  integer                              :: g
  type (grid_type), pointer            :: grid1, grid2

  grid1 => xmap%grids(1)
  conservation_check_side1(1) = mpp_global_sum(grid1%domain, grid1%area*d)
  conservation_check_side1(2) = 0.0
  call put_to_xgrid (d, grid1%id, x_over, xmap)    ! put from side 1
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

!#######################################################################
!
! conservation_check - returns three numbers which are the global sum of a
! variable (1) on its home model grid, (2) after interpolation to the other
! side grid(s), and (3) after re_interpolation back onto its home side grid(s).
!
function conservation_check_side2(d, grid_id, xmap) ! this one for 2->1->2
real, dimension(:,:,:), intent(in   ) :: d
character(len=3),       intent(in   ) :: grid_id
type (xmap_type),       intent(inout) :: xmap
real, dimension(3)                    :: conservation_check_side2

  real                       :: gsum
  real, dimension(xmap%size) :: x_over, x_back
  real, dimension(:,:  ), allocatable :: d1
  real, dimension(:,:,:), allocatable :: d2
  integer                             :: g
  type (grid_type), pointer           :: grid1, grid2

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
  call put_to_xgrid(d1,  grid1%id, x_back, xmap)   ! put from side 1
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

!#######################################################################


end module xgrid_mod
