module mpp_domains_define_mod

use mpp_mod,              only : mpp_init, mpp_root_pe, mpp_npes, stdlog, stderr, mpp_error, FATAL,   &
                                 NULL_PE, mpp_sync, mpp_sync_self, mpp_get_current_pelist,            &
                                 mpp_transmit, mpp_clock_end, mpp_clock_begin
use mpp_parameter_mod,    only : GLOBAL, CYCLIC, WEST, EAST, SOUTH, NORTH, CYCLIC_GLOBAL_DOMAIN,      &
                                 GLOBAL_DATA_DOMAIN, FOLD_WEST_EDGE, FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, &
                                 FOLD_NORTH_EDGE, DOMAIN_ID_BASE
use mpp_datatype_mod,     only : domain1D, domain2D
use mpp_data_mod,         only : module_is_initialized=>mpp_domains_is_initialized, pe,               &
                                 debug=>debug_mpp_domains, domain_info_buf, ptr_info, NULL_DOMAIN2D
use mpp_domains_util_mod, only : mpp_get_compute_domains, mpp_get_compute_domain, mpp_get_data_domain, &
                                 mpp_get_global_domain, mpp_get_layout, mpp_get_pelist,                &
                                 mpp_domains_set_stack_size, compute_overlaps

  implicit none
  private

#include <fms_platform.h>

  character(len=128), public :: version= &
       '$Id: mpp_domains_define.F90,v 12.0 2005/04/14 17:58:04 fms Exp $'
  character(len=128), public :: tagname= &
       '$Name: lima $'

  public :: mpp_define_layout, mpp_define_domains, mpp_modify_domain

  ! <INTERFACE NAME="mpp_define_layout">
  !  <OVERVIEW>
  !    Retrieve layout associated with a domain decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Given a global 2D domain and the number of divisions in the
  !    decomposition (<TT>ndivs</TT>: usually the PE count unless some
  !    domains are masked) this calls returns a 2D domain layout.
  !    
  !    By default, <TT>mpp_define_layout</TT> will attempt to divide the
  !    2D index space into domains that maintain the aspect ratio of the
  !    global domain. If this cannot be done, the algorithm favours domains
  !    that are longer in <TT>x</TT> than <TT>y</TT>, a preference that could
  !    improve vector performance.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_define_layout( global_indices, ndivs, layout )
  !  </TEMPLATE>
  !  <IN NAME="global_indices"></IN>
  !  <IN NAME="ndivs"></IN>
  !  <OUT NAME="layout"></OUT>
  ! </INTERFACE>

  interface mpp_define_layout
     module procedure mpp_define_layout2D
  end interface


  ! <INTERFACE NAME="mpp_define_domains">

  !   <OVERVIEW>
  !     Set up a domain decomposition.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     There are two forms for the <TT>mpp_define_domains</TT> call. The 2D
  !     version is generally to be used but is built by repeated calls to the
  !     1D version, also provided.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call mpp_define_domains( global_indices, ndivs, domain, &
  !                                   pelist, flags, halo, extent, maskmap )
  !   </TEMPLATE>
  !  <TEMPLATE>
  !    call mpp_define_domains( global_indices, layout, domain, pelist, &
  !                                   xflags, yflags, xhalo, yhalo,           &
  !                                   xextent, yextent, maskmap, name )
  !  </TEMPLATE>
  !   <IN NAME="global_indices" >
  !     Defines the global domain.
  !   </IN>
  !   <IN NAME="ndivs">
  !     Is the number of domain divisions required.
  !   </IN>
  !   <INOUT NAME="domain">
  !     Holds the resulting domain decomposition.
  !   </INOUT>
  !   <IN NAME="pelist">
  !     List of PEs to which the domains are to be assigned.
  !   </IN>
  !   <IN NAME="flags">
  !      An optional flag to pass additional information
  !      about the desired domain topology. Useful flags in a 1D decomposition
  !      include <TT>GLOBAL_DATA_DOMAIN</TT> and
  !      <TT>CYCLIC_GLOBAL_DOMAIN</TT>. Flags are integers: multiple flags may
  !      be added together. The flag values are public parameters available by
  !      use association.
  !   </IN>
  !   <IN NAME="halo">
  !     Width of the halo.
  !   </IN>
  !   <IN NAME="extent">
  !      Normally <TT>mpp_define_domains</TT> attempts
  !      an even division of the global domain across <TT>ndivs</TT>
  !      domains. The <TT>extent</TT> array can be used by the user to pass a
  !      custom domain division. The <TT>extent</TT> array has <TT>ndivs</TT>
  !      elements and holds the compute domain widths, which should add up to
  !      cover the global domain exactly.
  !   </IN>
  !   <IN NAME="maskmap">
  !     Some divisions may be masked
  !     (<TT>maskmap=.FALSE.</TT>) to exclude them from the computation (e.g
  !     for ocean model domains that are all land). The <TT>maskmap</TT> array
  !     is dimensioned <TT>ndivs</TT> and contains <TT>.TRUE.</TT> values for
  !     any domain that must be <I>included</I> in the computation (default
  !     all). The <TT>pelist</TT> array length should match the number of
  !     domains included in the computation.
  !    </IN>   

  !  <IN NAME="layout"></IN>
  !  <IN NAME="xflags, yflags"></IN>
  !  <IN NAME="xhalo, yhalo"></IN>
  !  <IN NAME="xextent, yextent"></IN>
  !  <IN NAME="name" ></IN>

  !  <NOTE>    
  !    For example:
  !    
  !    <PRE>
  !    call mpp_define_domains( (/1,100/), 10, domain, &
  !         flags=GLOBAL_DATA_DOMAIN+CYCLIC_GLOBAL_DOMAIN, halo=2 )
  !    </PRE>
  !    
  !    defines 10 compute domains spanning the range [1,100] of the global
  !    domain. The compute domains are non-overlapping blocks of 10. All the data
  !    domains are global, and with a halo of 2 span the range [-1:102]. And
  !    since the global domain has been declared to be cyclic,
  !    <TT>domain(9)%next => domain(0)</TT> and <TT>domain(0)%prev =>
  !    domain(9)</TT>. A field is allocated on the data domain, and computations proceed on
  !    the compute domain. A call to <LINK
  !    SRC="#mpp_update_domains"><TT>mpp_update_domains</TT></LINK> would fill in
  !    the values in the halo region:

  !    <PRE>
  !    call mpp_get_data_domain( domain, isd, ied ) !returns -1 and 102
  !    call mpp_get_compute_domain( domain, is, ie ) !returns (1,10) on PE 0 ...
  !    allocate( a(isd:ied) )
  !    do i = is,ie
  !       a(i) = &lt;perform computations&gt;
  !    end do
  !    call mpp_update_domains( a, domain )
  !    </PRE>

  !    The call to <TT>mpp_update_domains</TT> fills in the regions outside
  !    the compute domain. Since the global domain is cyclic, the values at
  !    <TT>i=(-1,0)</TT> are the same as at <TT>i=(99,100)</TT>; and
  !    <TT>i=(101,102)</TT> are the same as <TT>i=(1,2)</TT>.
  !    
  !    The 2D version is just an extension of this syntax to two
  !    dimensions.
  !
  !    The 2D version of the above should generally be used in
  !    codes, including 1D-decomposed ones, if there is a possibility of
  !    future evolution toward 2D decomposition. The arguments are similar to
  !    the 1D case, except that now we have optional arguments
  !    <TT>flags</TT>, <TT>halo</TT>, <TT>extent</TT> and <TT>maskmap</TT>
  !    along two axes.
  !    
  !    <TT>flags</TT> can now take an additional possible value to fold
  !    one or more edges. This is done by using flags
  !    <TT>FOLD_WEST_EDGE</TT>, <TT>FOLD_EAST_EDGE</TT>,
  !    <TT>FOLD_SOUTH_EDGE</TT> or <TT>FOLD_NORTH_EDGE</TT>. When a fold
  !    exists (e.g cylindrical domain), vector fields reverse sign upon
  !    crossing the fold. This parity reversal is performed only in the
  !    vector version of <LINK
  !    SRC="#mpp_update_domains"><TT>mpp_update_domains</TT></LINK>. In
  !    addition, shift operations may need to be applied to vector fields on
  !    staggered grids, also described in the vector interface to
  !    <TT>mpp_update_domains</TT>.
  !    
  !    <TT>name</TT> is the name associated with the decomposition,
  !    e.g <TT>'Ocean model'</TT>. If this argument is present,
  !    <TT>mpp_define_domains</TT> will print the domain decomposition
  !    generated to <TT>stdlog</TT>.
  !    
  !    Examples:
  !    
  !    <PRE>
  !    call mpp_define_domains( (/1,100,1,100/), (/2,2/), domain, xhalo=1 )
  !    </PRE>
  !    
  !    will create the following domain layout:
  !    <PRE>
  !                   |---------|-----------|-----------|-------------|
  !                   |domain(1)|domain(2)  |domain(3)  |domain(4)    |
  !    |--------------|---------|-----------|-----------|-------------|
  !    |Compute domain|1,50,1,50|51,100,1,50|1,50,51,100|51,100,51,100|
  !    |--------------|---------|-----------|-----------|-------------|
  !    |Data domain   |0,51,1,50|50,101,1,50|0,51,51,100|50,101,51,100|
  !    |--------------|---------|-----------|-----------|-------------|
  !    </PRE>
  !    
  !    Again, we allocate arrays on the data domain, perform computations
  !    on the compute domain, and call <TT>mpp_update_domains</TT> to update
  !    the halo region.
  !    
  !    If we wished to perfom a 1D decomposition along <TT>Y</TT>
  !    on the same global domain, we could use:

  !    <PRE>
  !    call mpp_define_domains( (/1,100,1,100/), layout=(/4,1/), domain, xhalo=1 )
  !    </PRE>

  !    This will create the following domain layout:
  !    <PRE>
  !                   |----------|-----------|-----------|------------|
  !                   |domain(1) |domain(2)  |domain(3)  |domain(4)   |
  !    |--------------|----------|-----------|-----------|------------|
  !    |Compute domain|1,100,1,25|1,100,26,50|1,100,51,75|1,100,76,100|
  !    |--------------|----------|-----------|-----------|------------|
  !    |Data domain   |0,101,1,25|0,101,26,50|0,101,51,75|1,101,76,100|
  !    |--------------|----------|-----------|-----------|------------|
  !    </PRE>
  !   </NOTE>
  ! </INTERFACE>
  interface mpp_define_domains
     module procedure mpp_define_domains1D
     module procedure mpp_define_domains2D
  end interface

! <INTERFACE NAME="mpp_modify_domain">
!   <OVERVIEW>
!     modifies the extents (compute, data and global) of domain
!   </OVERVIEW>
!   <IN NAME="domain_in">
!     The source domain.
!   </IN>
!   <IN NAME="halo">
!     Halo size of the returned 1D doamin. Default value is 0.
!   </IN>
!   <IN NAME="cbegin,cend">
!    Axis specifications associated with the compute domain of the returned 1D domain.
!   </IN>
!   <IN NAME="gbegin,gend">
!    Axis specifications associated with the global domain of the returned 1D domain.
!   </IN>
!   <IN NAME="isc,iec">
!    Zonal axis specifications associated with the compute domain of the returned 2D domain.
!   </IN>
!   <IN NAME="jsc,jec">
!    Meridinal axis specifications associated with the compute domain of the returned 2D domain.
!   </IN>
!   <IN NAME="isg,ieg">
!    Zonal axis specifications associated with the global domain of the returned 2D domain.
!   </IN>
!   <IN NAME="jsg,jeg">
!    Meridinal axis specifications associated with the global domain of the returned 2D domain.
!   </IN>
!   <IN NAME="xhalo,yhalo">
!     Halo size of the returned 2D doamin. Default value is 0.
!   </IN>
!   <INOUT NAME="domain_out">
!     The returned domain.
!   </INOUT>

! </INTERFACE>

  interface mpp_modify_domain
     module procedure mpp_modify_domain1D
     module procedure mpp_modify_domain2D
  end interface

contains


  ! <SUBROUTINE NAME="mpp_define_layout2D" INTERFACE="mpp_define_layout">
  !  <IN NAME="global_indices" TYPE="integer" DIM="(4)"></IN>
  !  <IN NAME="ndivs" TYPE="integer"></IN>
  !  <OUT NAME="layout" TYPE="integer" DIM="(2)"></OUT>
  ! </SUBROUTINE>
  subroutine mpp_define_layout2D( global_indices, ndivs, layout )
    integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
    integer, intent(in) :: ndivs !number of divisions to divide global domain
    integer, intent(out) :: layout(2)

    integer :: isg, ieg, jsg, jeg, isz, jsz, idiv, jdiv

    isg = global_indices(1)
    ieg = global_indices(2)
    jsg = global_indices(3)
    jeg = global_indices(4)

    isz = ieg - isg + 1
    jsz = jeg - jsg + 1
    !first try to divide ndivs in the domain aspect ratio: if imperfect aspect, reduce idiv till it divides ndivs
    idiv = nint( sqrt(float(ndivs*isz)/jsz) )
    idiv = max(idiv,1) !for isz=1 line above can give 0
    do while( mod(ndivs,idiv).NE.0 )
       idiv = idiv - 1
    end do                 !will terminate at idiv=1 if not before
    jdiv = ndivs/idiv

    layout = (/ idiv, jdiv /)
    return
  end subroutine mpp_define_layout2D

  !#####################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !              MPP_DEFINE_DOMAINS: define layout and decomposition            !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! <SUBROUTINE NAME="mpp_define_domains1D" INTERFACE="mpp_define_domains">>
  !   <IN NAME="global_indices" TYPE="integer" DIM="(2)"> </IN>
  !   <IN NAME="ndivs" TYPE="integer">  </IN>
  !   <INOUT NAME="domain" TYPE="type(domain1D)"> </INOUT>
  !   <IN NAME="pelist" TYPE="integer" DIM="(0:)">  </IN>
  !   <IN NAME="flags" TYPE="integer">  </IN>
  !   <IN NAME="halo" TYPE="integer">  </IN>
  !   <IN NAME="extent" TYPE="integer" DIM="(0:)">  </IN>
  !   <IN NAME="maskmap" TYPE="logical" DIM="(0:)"> </IN>
  ! </SUBROUTINE>
  subroutine mpp_define_domains1D( global_indices, ndivs, domain, pelist, flags, halo, extent, maskmap )
    !routine to divide global array indices among domains, and assign domains to PEs
    !domain is of type domain1D
    !ARGUMENTS:
    !      global_indices(2)=(isg,ieg) gives the extent of global domain
    !      ndivs is number of divisions of domain: even divisions unless extent is present.
    !      domain is the returned domain1D
    !      pelist (optional) list of PEs to which domains are to be assigned (default 0...npes-1)
    !                 size of pelist must correspond to number of mask=.TRUE. divisions
    !      flags define whether compute and data domains are global (undecomposed) and whether global domain has periodic boundaries
    !      halo (optional) defines halo width (currently the same on both sides)
    !      extent (optional) array defines width of each division (used for non-uniform domain decomp, for e.g load-balancing)
    !      maskmap (optional) a division whose maskmap=.FALSE. is not assigned to any domain
    !  By default we assume decomposition of compute and data domains, non-periodic boundaries, no halo, as close to uniform extents
    !  as the input parameters permit
    integer, intent(in) :: global_indices(2) !(/ isg, ieg /)
    integer, intent(in) :: ndivs
    type(domain1D), intent(inout) :: domain !declared inout so that existing links, if any, can be nullified
    integer, intent(in), optional :: pelist(0:)
    integer, intent(in), optional :: flags, halo
    integer, intent(in), optional :: extent(0:)
    logical, intent(in), optional :: maskmap(0:)

    logical :: compute_domain_is_global, data_domain_is_global
    integer :: ndiv, n, isg, ieg, is, ie, i, off, pos, hs, he
    integer, allocatable :: pes(:)
    logical, allocatable :: mask(:)
    integer :: halosz
    !used by symmetry algorithm
    integer :: imax, ndmax, ndmirror
    logical :: symmetrize
    !statement functions
    logical :: even, odd
    even(n) = (mod(n,2).EQ.0)
    odd (n) = (mod(n,2).EQ.1)

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: You must first call mpp_domains_init.' )
    !get global indices
    isg = global_indices(1)
    ieg = global_indices(2)
    if( ndivs.GT.ieg-isg+1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: more divisions requested than rows available.' )
    !get the list of PEs on which to assign domains; if pelist is absent use 0..npes-1
    if( PRESENT(pelist) )then
       if( .NOT.any(pelist.EQ.pe) )then
          write( stderr(),* )'pe=', pe, ' pelist=', pelist
          call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: pe must be in pelist.' )
       end if
       allocate( pes(0:size(pelist(:))-1) )
       pes(:) = pelist(:)
    else
       allocate( pes(0:mpp_npes()-1) )
       pes(:) = (/ (i,i=0,mpp_npes()-1) /)
    end if

    !get number of real domains: 1 mask domain per PE in pes
    allocate( mask(0:ndivs-1) )
    mask = .TRUE.                 !default mask
    if( PRESENT(maskmap) )then
       if( size(maskmap(:)).NE.ndivs ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: maskmap array size must equal number of domain divisions.' )
       mask(:) = maskmap(:)
    end if
    if( count(mask).NE.size(pes(:)) ) &
         call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: number of TRUEs in maskmap array must match PE count.' )
    if( PRESENT(extent) )then
       if( size(extent(:)).NE.ndivs ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: extent array size must equal number of domain divisions.' )
    end if
    !get halosize
    halosz = 0
    if( PRESENT(halo) )halosz = halo

    !get flags
    compute_domain_is_global = .FALSE.
    data_domain_is_global    = .FALSE.
    domain%cyclic = .FALSE.
    if( PRESENT(flags) )then
       !NEW: obsolete flag global_compute_domain, since ndivs is non-optional and you cannot have global compute and ndivs.NE.1 
       compute_domain_is_global = ndivs.EQ.1
       !if compute domain is global, data domain must also be
       data_domain_is_global    = BTEST(flags,GLOBAL) .OR. compute_domain_is_global
       domain%cyclic  = BTEST(flags,CYCLIC) .AND. halosz.NE.0
    end if

    !set up links list
    allocate( domain%list(0:ndivs-1) )

    !set global domain
    domain%list(:)%global%begin     = isg
    domain%list(:)%global%end       = ieg
    domain%list(:)%global%size      = ieg-isg+1
    domain%list(:)%global%max_size  = ieg-isg+1
    domain%list(:)%global%is_global = .TRUE. !always

    !get compute domain
    if( compute_domain_is_global )then
       domain%list(:)%compute%begin = isg
       domain%list(:)%compute%end   = ieg
       domain%list(:)%compute%is_global = .TRUE.
       domain%list(:)%pe = pes(:)
       domain%pos = 0
    else
       domain%list(:)%compute%is_global = .FALSE.
       is = isg
       n = 0
       do ndiv=0,ndivs-1
          if( PRESENT(extent) )then
             ie = is + extent(ndiv) - 1
             if( ndiv.EQ.ndivs-1 .AND. ie.NE.ieg ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: extent array limits do not match global domain.' )
          else
             !modified for mirror-symmetry
             !original line
             !                 ie = is + CEILING( float(ieg-is+1)/(ndivs-ndiv) ) - 1

             !problem of dividing nx points into n domains maintaining symmetry
             !i.e nx=18 n=4 4554 and 5445 are solutions but 4455 is not.
             !this will always work for nx even n even or odd
             !this will always work for nx odd, n odd
             !this will never  work for nx odd, n even: for this case we supersede the mirror calculation
             !                 symmetrize = .NOT. ( mod(ndivs,2).EQ.0 .AND. mod(ieg-isg+1,2).EQ.1 )
             !nx even n odd fails if n>nx/2
             symmetrize = ( even(ndivs) .AND. even(ieg-isg+1) ) .OR. &
                  (  odd(ndivs) .AND.  odd(ieg-isg+1) ) .OR. &
                  (  odd(ndivs) .AND. even(ieg-isg+1) .AND. ndivs.LT.(ieg-isg+1)/2 )

             !mirror domains are stored in the list and retrieved if required.
             if( ndiv.EQ.0 )then
                !initialize max points and max domains
                imax = ieg
                ndmax = ndivs
             end if
             !do bottom half of decomposition, going over the midpoint for odd ndivs
             if( ndiv.LT.(ndivs-1)/2+1 )then
                !domain is sized by dividing remaining points by remaining domains
                ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                ndmirror = (ndivs-1) - ndiv !mirror domain
                if( ndmirror.GT.ndiv .AND. symmetrize )then !only for domains over the midpoint
                   !mirror extents, the max(,) is to eliminate overlaps
                   domain%list(ndmirror)%compute%begin = max( isg+ieg-ie, ie+1 )
                   domain%list(ndmirror)%compute%end   = max( isg+ieg-is, ie+1 )
                   imax = domain%list(ndmirror)%compute%begin - 1
                   ndmax = ndmax - 1
                end if
             else
                if( symmetrize )then
                   !do top half of decomposition by retrieving saved values
                   is = domain%list(ndiv)%compute%begin
                   ie = domain%list(ndiv)%compute%end
                else
                   ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                end if
             end if
          end if
          if( ie.LT.is )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents must be positive definite.' )
          domain%list(ndiv)%compute%begin = is
          domain%list(ndiv)%compute%end   = ie
          if( ndiv.GT.0 ) then
            if( is.NE.domain%list(ndiv-1)%compute%end+1 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents do not span space completely.' )
          endif
          if( ndiv.EQ.ndivs-1 .AND. domain%list(ndiv)%compute%end.NE.ieg ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents do not span space completely.' )
          if( mask(ndiv) )then
             domain%list(ndiv)%pe = pes(n)
             if( pe.EQ.pes(n) )domain%pos = ndiv
             n = n + 1
          end if
          is = ie + 1
       end do
    end if

    domain%list(:)%compute%size  = domain%list(:)%compute%end - domain%list(:)%compute%begin + 1

    !get data domain
    !data domain is at least equal to compute domain
    domain%list(:)%data%begin = domain%list(:)%compute%begin
    domain%list(:)%data%end   = domain%list(:)%compute%end
    domain%list(:)%data%is_global = .FALSE.
    !apply global flags
    if( data_domain_is_global )then
       domain%list(:)%data%begin  = isg
       domain%list(:)%data%end    = ieg
       domain%list(:)%data%is_global = .TRUE.
    end if
    !apply margins
    domain%list(:)%data%begin = domain%list(:)%data%begin - halosz
    domain%list(:)%data%end   = domain%list(:)%data%end   + halosz  
    domain%list(:)%data%size  = domain%list(:)%data%end - domain%list(:)%data%begin + 1

    !      domain = domain%list(pos) !load domain from domain%list(pos)
    domain%compute = domain%list(domain%pos)%compute
    domain%data = domain%list(domain%pos)%data
    domain%global = domain%list(domain%pos)%global
    domain%compute%max_size = MAXVAL( domain%list(:)%compute%size )
    domain%data%max_size    = MAXVAL( domain%list(:)%data%size )
    domain%global%max_size  = domain%global%size

    !PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
    deallocate( pes, mask )
    return

  contains

    function if_overlap( hs, he, cs, ce, os, oe )
      !function to look if a portion of the halo [hs,he] lies with in the compute region [cs,ce]
      !if yes, if_overlap returns true, and the overlap region is returned in [os,oe]
      logical :: if_overlap
      integer, intent(in) :: hs, he, cs, ce
      integer, intent(out) :: os, oe
      os = max(hs,cs)
      oe = min(he,ce)
      if( debug )write( stderr(),'(a,7i4)' ) &
           'MPP_DEFINE_DOMAINS1D: pe, hs, he, cs, ce, os, oe=', pe, hs, he, cs, ce, os, oe
      if_overlap = oe.GE.os
      return
    end function if_overlap

  end subroutine mpp_define_domains1D
  ! <SUBROUTINE NAME="mpp_define_domains2D" INTERFACE="mpp_define_domains">
  !  <IN NAME="global_indices" TYPE="integer" DIM="(4)"> </IN>
  !  <IN NAME="layout" TYPE="integer" DIM="(2)"></IN>
  !  <INOUT NAME="domain" TYPE="type(domain2D)"></INOUT>
  !  <IN NAME="pelist" TYPE="integer" DIM="(0:)"></IN>
  !  <IN NAME="xflags, yflags" TYPE="integer"></IN>
  !  <IN NAME="xhalo, yhalo" TYPE="integer"></IN>
  !  <IN NAME="xextent, yextent" TYPE="integer" DIM="(0:)"></IN>
  !  <IN NAME="maskmap" TYPE="logical" DIM="(:,:)"></IN>
  !  <IN NAME="name" TYPE="character(len=*)"></IN>
  ! </SUBROUTINE>
  subroutine mpp_define_domains2D( global_indices, layout, domain, pelist, &
       xflags, yflags, xhalo, yhalo, xextent, yextent, maskmap, name )
    !define 2D data and computational domain on global rectilinear cartesian domain (isg:ieg,jsg:jeg) and assign them to PEs
    integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
    integer, intent(in) :: layout(2)
    type(domain2D), intent(inout) :: domain
    integer, intent(in), optional :: pelist(0:)
    integer, intent(in), optional :: xflags, yflags, xhalo, yhalo
    integer, intent(in), optional :: xextent(0:), yextent(0:)
    logical, intent(in), optional :: maskmap(0:,0:)
    character(len=*), intent(in), optional :: name
    integer :: i, j, m, n
    integer :: ipos, jpos, pos
    integer :: ndivx, ndivy, isg, ieg, jsg, jeg, isd, ied, jsd, jed
    integer(LONG_KIND),save :: domain_cnt=0
    integer(LONG_KIND)      :: d_base

    logical, allocatable :: mask(:,:)
    integer, allocatable :: pes(:), pearray(:,:)
    character(len=8) :: text

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: You must first call mpp_domains_init.' )
    ndivx = layout(1); ndivy = layout(2)
    isg = global_indices(1); ieg = global_indices(2); jsg = global_indices(3); jeg = global_indices(4)

    if( PRESENT(pelist) )then
       allocate( pes(0:size(pelist(:))-1) )
       pes = pelist
    else
       allocate( pes(0:mpp_npes()-1) )
       call mpp_get_current_pelist(pes)
       !          pes = (/ (i,i=0,mpp_npes()-1) /)
    end if

    allocate( mask(0:ndivx-1,0:ndivy-1) )
    mask = .TRUE.
    if( PRESENT(maskmap) )then
       if( size(maskmap,1).NE.ndivx .OR. size(maskmap,2).NE.ndivy ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: maskmap array does not match layout.' )
       mask(:,:) = maskmap(:,:)
    end if
    !number of unmask domains in layout must equal number of PEs assigned
    n = count(mask)
    if( n.NE.size(pes(:)) )then
       write( text,'(i8)' )n
       call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: incorrect number of PEs assigned for this layout and maskmap. Use ' &
            //text//' PEs for this domain decomposition.' )
    end if

    !place on PE array; need flag to assign them to j first and then i
    allocate( pearray(0:ndivx-1,0:ndivy-1) )
    pearray(:,:) = NULL_PE
    ipos = NULL_PE; jpos = NULL_PE; pos = NULL_PE
    n = 0
    do j = 0,ndivy-1
       do i = 0,ndivx-1
          if( mask(i,j) )then
             pearray(i,j) = pes(n)
             if( pes(n).EQ.pe )then
                pos = n
                ipos = i
                jpos = j
             end if
             n = n + 1
          end if
       end do
    end do
    if( ipos.EQ.NULL_PE .OR. jpos.EQ.NULL_PE .or. pos.EQ.NULL_PE ) &
         call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: pelist must include this PE.' )
    if( debug )write( stderr(), * )'pe, ipos, jpos=', pe, ipos, jpos, ' pearray(:,jpos)=', pearray(:,jpos), &
         ' pearray(ipos,:)=', pearray(ipos,:)

    !do domain decomposition using 1D versions in X and Y
    call mpp_define_domains( global_indices(1:2), ndivx, domain%x, &
         pack(pearray(:,jpos),mask(:,jpos)), xflags, xhalo, xextent, mask(:,jpos) )
    call mpp_define_domains( global_indices(3:4), ndivy, domain%y, &
         pack(pearray(ipos,:),mask(ipos,:)), yflags, yhalo, yextent, mask(ipos,:) )
    if( domain%x%list(domain%x%pos)%pe.NE.domain%y%list(domain%y%pos)%pe ) &
         call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: domain%x%list(ipos)%pe.NE.domain%y%list(jpos)%pe.' ) 
    domain%pos = pos
    domain%pe  = pe

    domain_cnt = domain_cnt + INT(1,KIND=LONG_KIND)
    d_base = DOMAIN_ID_BASE
    domain%id = domain_cnt*DOMAIN_ID_BASE  ! Must be LONG_KIND arithmetic

    !set up fold
    domain%fold = 0
    if( PRESENT(xflags) )then
       if( BTEST(xflags,WEST) )domain%fold = domain%fold + FOLD_WEST_EDGE
       if( BTEST(xflags,EAST) )domain%fold = domain%fold + FOLD_EAST_EDGE
    end if
    if( PRESENT(yflags) )then
       if( BTEST(yflags,SOUTH) )domain%fold = domain%fold + FOLD_SOUTH_EDGE
       if( BTEST(yflags,NORTH) )domain%fold = domain%fold + FOLD_NORTH_EDGE
    end if

    if( BTEST(domain%fold,SOUTH) .OR. BTEST(domain%fold,NORTH) )then
       if( domain%y%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
       if( modulo(domain%x%global%size,2).NE.0 ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in X must be even when there is a fold in Y.' )
       !check if folded domain boundaries line up in X: compute domains lining up is a sufficient condition for symmetry
       n = ndivx - 1
       do i = 0,n/2
          if( domain%x%list(i)%compute%size.NE.domain%x%list(n-i)%compute%size ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
       end do
    end if
    if( BTEST(domain%fold,WEST) .OR. BTEST(domain%fold,EAST) )then
       if( domain%x%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
       if( modulo(domain%y%global%size,2).NE.0 ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in Y must be even when there is a fold in X.' )
       !check if folded domain boundaries line up in Y: compute domains lining up is a sufficient condition for symmetry
       n = ndivy - 1
       do i = 0,n/2
          if( domain%y%list(i)%compute%size.NE.domain%y%list(n-i)%compute%size ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
       end do
    end if

    !set up domain%list
    if( debug )write( stderr(),'(a,9i4)' )'pe, domain=', pe, domain_info_buf(1:8)
    if( pe.EQ.mpp_root_pe() .AND. PRESENT(name) )then
       write( stdlog(), '(/a,i3,a,i3)' )trim(name)//' domain decomposition: ', ndivx, ' X', ndivy
       write( stdlog(), '(3x,a)' )'pe,   is,  ie,  js,  je,    isd, ied, jsd, jed'
    end if
    call mpp_sync()
    call mpp_get_compute_domain( domain, domain_info_buf(1), domain_info_buf(2), domain_info_buf(3), domain_info_buf(4) )
    call mpp_get_data_domain   ( domain, domain_info_buf(5), domain_info_buf(6), domain_info_buf(7), domain_info_buf(8) )
    n = size(pes(:))
    allocate( domain%list(0:n-1) ) !this is only used for storage of remote compute and data domain info
    do i = 0,n-1
       m = mod(pos+i,n)
       domain%list(m)%pe = pes(m)
     ! Force use of "scalar", integer pointer mpp interface
       call mpp_transmit( put_data=domain_info_buf(1), plen=8, to_pe=pes(mod(pos+n-i,n)), &
                          get_data=domain_info_buf(9), glen=8, from_pe=pes(m) )
       domain%list(m)%x%compute%begin = domain_info_buf(9)
       domain%list(m)%x%compute%end   = domain_info_buf(10)
       domain%list(m)%y%compute%begin = domain_info_buf(11)
       domain%list(m)%y%compute%end   = domain_info_buf(12)
       domain%list(m)%x%data%begin = domain_info_buf(13)
       domain%list(m)%x%data%end   = domain_info_buf(14)
       domain%list(m)%y%data%begin = domain_info_buf(15)
       domain%list(m)%y%data%end   = domain_info_buf(16)
       if( pe.EQ.mpp_root_pe() .AND. PRESENT(name) )write( stdlog(), '(2x,i3,x,4i5,3x,4i5)' )pes(m), domain_info_buf(9:)
    end do
    call mpp_sync_self(pes)
    domain%list(:)%x%compute%size = domain%list(:)%x%compute%end - domain%list(:)%x%compute%begin + 1
    domain%list(:)%y%compute%size = domain%list(:)%y%compute%end - domain%list(:)%y%compute%begin + 1
    domain%list(:)%x%data%size = domain%list(:)%x%data%end - domain%list(:)%x%data%begin + 1
    domain%list(:)%y%data%size = domain%list(:)%y%data%end - domain%list(:)%y%data%begin + 1

    domain%remote_domains_initialized = .FALSE.
    domain%remote_off_domains_initialized = .FALSE.
    call compute_overlaps(domain)
    !PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)

    !print out decomposition
    if( pe.EQ.mpp_root_pe() .AND. PRESENT(name) )then
       write(*,*) trim(name)//' domain decomposition'
       write (*,110) (domain%list(i)%x%compute%size, i= 0, layout(1)-1)
       write (*,120) (domain%list(i)%y%compute%size, i= 0, layout(2)-1)
110    format ('  X-AXIS = ',24i4,/,(11x,24i4))
120    format ('  Y-AXIS = ',24i4,/,(11x,24i4))
    endif

    deallocate( pes, mask, pearray )

    return
  end subroutine mpp_define_domains2D

  !#####################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_MODIFY_DOMAIN: modify extent of domain                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <SUBROUTINE NAME="mpp_modify_domain1D" INTERFACE="mpp_modify_domain">
!   <IN NAME="domain_in" TYPE="type(domain1D)" > </IN>
!   <IN NAME="halo" TYPE="integer,optional" > </IN>
!   <IN NAME="cbegin,cend" TYPE="integer,optional" > </IN>
!   <IN NAME="gbegin,gend" TYPE="integer,optional" > </IN>
!   <INOUT NAME="domain_out" TYPE="type(domain1D)" > </INOUT>

! <PUBLICROUTINE>
    subroutine mpp_modify_domain1D(domain_in,domain_out,cbegin,cend,gbegin,gend, halo)
      ! </PUBLICROUTINE>
      type(domain1D), intent(in)    :: domain_in
      type(domain1D), intent(inout) :: domain_out
      integer, intent(in), optional :: halo
      integer, intent(in), optional :: cbegin, cend             ! extent of compute_domain
      integer, intent(in), optional :: gbegin, gend             ! extent of global domain
      integer :: ndivs, global_indices(2) !(/ isg, ieg /)
      integer                       :: halosz, flag
! get the global indices of the input domain
      global_indices(1) = domain_in%global%begin;  global_indices(2) = domain_in%global%end

! get the halo
      halosz = 0
      if(present(halo)) halosz = halo

      ! get the layout
      ndivs = size(domain_in%list(:))

! get the flag
      flag = 0
      if(domain_in%cyclic) flag = flag + CYCLIC_GLOBAL_DOMAIN
      if(domain_in%data%is_global) flag = flag + GLOBAL_DATA_DOMAIN

      call mpp_define_domains( global_indices, ndivs, domain_out, pelist = domain_in%list(:)%pe, &
           flags = flag, halo = halosz, extent = domain_in%list(:)%compute%size )
           
      if(present(cbegin)) domain_out%compute%begin = cbegin
      if(present(cend))   domain_out%compute%end = cend
      domain_out%compute%size = domain_out%compute%end - domain_out%compute%begin + 1
      if(present(gbegin)) domain_out%global%begin = gbegin
      if(present(gend))   domain_out%global%end = gend
      domain_out%global%size = domain_out%global%end - domain_out%global%begin + 1
      
    end subroutine mpp_modify_domain1D
! </SUBROUTINE>

  !#######################################################################
!----------------------------------------------------------------------------------
! <SUBROUTINE NAME="mpp_modify_domain2D" INTERFACE="mpp_modify_domain">
!   <IN NAME="domain_in" TYPE="type(domain2D)" > </IN>
!   <IN NAME="isc,iec" TYPE="integer,optional" > </IN>
!   <IN NAME="jsc,jec" TYPE="integer,optional" > </IN>
!   <IN NAME="isg,ieg" TYPE="integer,optional" > </IN>
!   <IN NAME="jsg,jeg" TYPE="integer,optional" > </IN>
!   <IN NAME="xhalo,yhalo" TYPE="integer,optional" > </IN>
!   <INOUT NAME="domain_out" TYPE="type(domain2D)" > </INOUT>

! <PUBLICROUTINE>
    subroutine mpp_modify_domain2D(domain_in, domain_out, isc, iec, jsc, jec, isg, ieg, jsg, jeg, xhalo, yhalo)
      ! </PUBLICROUTINE>
      type(domain2D), intent(in)    :: domain_in
      type(domain2D), intent(inout) :: domain_out
      integer, intent(in), optional :: isc, iec, jsc, jec
      integer, intent(in), optional :: isg, ieg, jsg, jeg
      integer, intent(in), optional :: xhalo, yhalo
      integer                       :: ndivx, ndivy, global_indices(4), layout(2)
      integer                       :: xhalosz, yhalosz, xflag, yflag
      
      if(present(xhalo) .or. present(yhalo)) then
! get the global indices of the input domain
         global_indices(1) = domain_in%x%global%begin;  global_indices(2) = domain_in%x%global%end
         global_indices(3) = domain_in%y%global%begin;  global_indices(4) = domain_in%y%global%end
         
! get the halo
      xhalosz = 0; yhalosz = 0
      if(present(xhalo)) xhalosz = xhalo
      if(present(yhalo)) yhalosz = yhalo
      
         ! get the layout
         layout(1) = size(domain_in%x%list(:)); layout(2) = size(domain_in%y%list(:))
         
! get the flag
      xflag = 0; yflag = 0
      if(domain_in%x%cyclic) xflag = xflag + CYCLIC_GLOBAL_DOMAIN
      if(domain_in%x%data%is_global) xflag = xflag + GLOBAL_DATA_DOMAIN
      if(domain_in%y%cyclic) yflag = yflag + CYCLIC_GLOBAL_DOMAIN
      if(domain_in%y%data%is_global) yflag = yflag + GLOBAL_DATA_DOMAIN
      
         call mpp_define_domains( global_indices, layout, domain_out, pelist = domain_in%list(:)%pe, &
                               xflags = xflag, yflags = yflag,  xhalo = xhalosz,    &
              yhalo = yhalosz, xextent = domain_in%x%list(:)%compute%size, yextent = domain_in%y%list(:)%compute%size)
                               
      else    
  domain_out = NULL_DOMAIN2D
         call mpp_modify_domain(domain_in%x, domain_out%x, isc, iec, isg, ieg)
         call mpp_modify_domain(domain_in%y, domain_out%y, jsc, jec, jsg, jeg)
      endif
         
    end subroutine mpp_modify_domain2D
! </SUBROUTINE>

  !#####################################################################

end module mpp_domains_define_mod
