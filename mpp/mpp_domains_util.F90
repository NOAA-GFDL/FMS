module mpp_domains_util_mod

  use mpp_mod,           only : mpp_root_pe, mpp_malloc, mpp_error, FATAL, NOTE, NULL_PE
  use mpp_parameter_mod, only : WEST, EAST, SOUTH, NORTH, AGRID
  use mpp_datatype_mod,  only : domain1D, domain2D
  use mpp_data_mod,      only : mpp_domains_stack, ptr_domains_stack, pe, mpp_domains_stack_size
  use mpp_data_mod,      only : grid_offset_type, module_is_initialized=>mpp_domains_is_initialized

  implicit none
  private

  character(len=128), public :: version= &
       '$Id mpp_domains_util.F90 $'
  character(len=128), public :: tagname= &
       '$Name: lima $'

  public :: mpp_domains_set_stack_size, mpp_get_compute_domain, mpp_get_compute_domains
  public :: mpp_get_data_domain, mpp_get_global_domain, mpp_get_domain_components
  public :: mpp_get_layout, mpp_get_pelist, operator(.EQ.), operator(.NE.), compute_overlaps

  !--- public interface

  ! <INTERFACE NAME="operator">
  !  <OVERVIEW>
  !    Equality/inequality operators for domaintypes.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The module provides public operators to check for
  !    equality/inequality of domaintypes, e.g:
  !    
  !    <PRE>
  !    type(domain1D) :: a, b
  !    type(domain2D) :: c, d
  !    ...
  !    if( a.NE.b )then
  !        ...
  !    end if
  !    if( c==d )then
  !        ...
  !    end if
  !    </PRE>
  !    
  !    Domains are considered equal if and only if the start and end
  !    indices of each of their component global, data and compute domains
  !    are equal.
  !  </DESCRIPTION>
  ! </INTERFACE>
  interface operator(.EQ.)
     module procedure mpp_domain1D_eq
     module procedure mpp_domain2D_eq
  end interface

  interface operator(.NE.)
     module procedure mpp_domain1D_ne
     module procedure mpp_domain2D_ne
  end interface

  ! <INTERFACE NAME="mpp_get_compute_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the compute domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the compute domains
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_compute_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_compute_domain
     module procedure mpp_get_compute_domain1D
     module procedure mpp_get_compute_domain2D
  end interface

  ! <INTERFACE NAME="mpp_get_compute_domains">
  !  <OVERVIEW>
  !    Retrieve the entire array of compute domain extents associated with a decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Retrieve the entire array of compute domain extents associated with a decomposition.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_compute_domains( domain, xbegin, xend, xsize, &
  !                                                ybegin, yend, ysize )
  !  </TEMPLATE>
  !  <IN NAME="domain" TYPE="type(domain2D)"></IN>
  !  <OUT NAME="xbegin,ybegin" TYPE="integer" DIM="(:)"></OUT>
  !  <OUT NAME="xend,yend" TYPE="integer" DIM="(:)"></OUT>
  !  <OUT NAME="xsize,ysize" TYPE="integer" DIM="(:)"></OUT>
  ! </INTERFACE>
  interface mpp_get_compute_domains
     module procedure mpp_get_compute_domains1D
     module procedure mpp_get_compute_domains2D
  end interface

  ! <INTERFACE NAME="mpp_get_data_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the data domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the data domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_data_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_data_domain
     module procedure mpp_get_data_domain1D
     module procedure mpp_get_data_domain2D
  end interface

  ! <INTERFACE NAME="mpp_get_global_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the global domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the global domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_global_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_global_domain
     module procedure mpp_get_global_domain1D
     module procedure mpp_get_global_domain2D
  end interface

  ! <INTERFACE NAME="mpp_get_pelist">
  !  <OVERVIEW>
  !    Retrieve list of PEs associated with a domain decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The 1D version of this call returns an array of the PEs assigned to this 1D domain
  !    decomposition. In addition the optional argument <TT>pos</TT> may be
  !    used to retrieve the 0-based position of the domain local to the
  !    calling PE, i.e <TT>domain%list(pos)%pe</TT> is the local PE,
  !    as returned by <LINK SRC="mpp.html#mpp_pe"><TT>mpp_pe()</TT></LINK>.
  !    The 2D version of this call is identical to 1D version.
  !  </DESCRIPTION>
  !  <IN NAME="domain"></IN>
  !  <OUT NAME="pelist"></OUT>
  !  <OUT NAME="pos"></OUT>
  ! </INTERFACE>
  interface mpp_get_pelist
     module procedure mpp_get_pelist1D
     module procedure mpp_get_pelist2D
  end interface

  ! <INTERFACE NAME="mpp_get_layout">
  !  <OVERVIEW>
  !    Retrieve layout associated with a domain decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The 1D version of this call returns the number of divisions that was assigned to this
  !    decomposition axis. The 2D version of this call returns an array of
  !    dimension 2 holding the results on two axes.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_layout( domain, layout )
  !  </TEMPLATE>
  !  <IN NAME="domain"></IN>
  !  <OUT NAME="layout"></OUT>
  ! </INTERFACE>
  interface mpp_get_layout
     module procedure mpp_get_layout1D
     module procedure mpp_get_layout2D
  end interface

contains

  ! <SUBROUTINE NAME="mpp_domains_set_stack_size">
  !  <OVERVIEW>
  !    Set user stack size.
  ! </OVERVIEW>
  ! <DESCRIPTION>
  !    This sets the size of an array that is used for internal storage by
  !    <TT>mpp_domains</TT>. This array is used, for instance, to buffer the
  !    data sent and received in halo updates.
  !    
  !    This call has implied global synchronization. It should be
  !    placed somewhere where all PEs can call it.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_domains_set_stack_size(n)
  !  </TEMPLATE>
  !  <IN NAME="n" TYPE="integer"></IN>
  ! </SUBROUTINE>
  subroutine mpp_domains_set_stack_size(n)
    !set the mpp_domains_stack variable to be at least n LONG words long
    integer, intent(in) :: n
    character(len=8) :: text

    if( n.LE.mpp_domains_stack_size )return
#ifdef use_libSMA
    call mpp_malloc( ptr_domains_stack, n, mpp_domains_stack_size )
#else
    if( allocated(mpp_domains_stack) )deallocate(mpp_domains_stack)
    allocate( mpp_domains_stack(n) )
    mpp_domains_stack_size = n
#endif
    write( text,'(i8)' )n
    if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, 'MPP_DOMAINS_SET_STACK_SIZE: stack size set to '//text//'.' )

    return
  end subroutine mpp_domains_set_stack_size


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !                MPP_DOMAINS: overloaded operators (==, /=)                   !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function mpp_domain1D_eq( a, b )
    logical                    :: mpp_domain1D_eq
    type(domain1D), intent(in) :: a, b

    mpp_domain1D_eq = ( a%compute%begin.EQ.b%compute%begin .AND. &
         a%compute%end  .EQ.b%compute%end   .AND. &
         a%data%begin   .EQ.b%data%begin    .AND. &
         a%data%end     .EQ.b%data%end      .AND. & 
         a%global%begin .EQ.b%global%begin  .AND. &
         a%global%end   .EQ.b%global%end    )
    !compare pelists
    !      if( mpp_domain1D_eq )mpp_domain1D_eq = ASSOCIATED(a%list) .AND. ASSOCIATED(b%list)
    !      if( mpp_domain1D_eq )mpp_domain1D_eq = size(a%list(:)).EQ.size(b%list(:))
    !      if( mpp_domain1D_eq )mpp_domain1D_eq = ALL(a%list%pe.EQ.b%list%pe)

    return
  end function mpp_domain1D_eq

  function mpp_domain1D_ne( a, b )
    logical                    :: mpp_domain1D_ne
    type(domain1D), intent(in) :: a, b

    mpp_domain1D_ne = .NOT. ( a.EQ.b )
    return
  end function mpp_domain1D_ne

  function mpp_domain2D_eq( a, b )
    logical                    :: mpp_domain2D_eq
    type(domain2D), intent(in) :: a, b

    mpp_domain2D_eq = a%x.EQ.b%x .AND. a%y.EQ.b%y
    if( mpp_domain2D_eq .AND. ((a%pe.EQ.NULL_PE).OR.(b%pe.EQ.NULL_PE)) )return !NULL_DOMAIN2D
    !compare pelists
    if( mpp_domain2D_eq )mpp_domain2D_eq = ASSOCIATED(a%list) .AND. ASSOCIATED(b%list)
    if( mpp_domain2D_eq )mpp_domain2D_eq = size(a%list(:)).EQ.size(b%list(:))
    if( mpp_domain2D_eq )mpp_domain2D_eq = ALL(a%list%pe.EQ.b%list%pe)
    return
  end function mpp_domain2D_eq

  !#####################################################################

  function mpp_domain2D_ne( a, b )
    logical                    :: mpp_domain2D_ne
    type(domain2D), intent(in) :: a, b

    mpp_domain2D_ne = .NOT. ( a.EQ.b )
    return
  end function mpp_domain2D_ne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !     MPP_GET and SET routiness: retrieve various components of domains       !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpp_get_compute_domain1D( domain, begin, end, size, max_size, is_global )
    type(domain1D),     intent(in) :: domain
    integer, intent(out), optional :: begin, end, size, max_size
    logical, intent(out), optional :: is_global

    if( PRESENT(begin)     )begin     = domain%compute%begin
    if( PRESENT(end)       )end       = domain%compute%end
    if( PRESENT(size)      )size      = domain%compute%size
    if( PRESENT(max_size)  )max_size  = domain%compute%max_size
    if( PRESENT(is_global) )is_global = domain%compute%is_global
    return
  end subroutine mpp_get_compute_domain1D

  !#####################################################################
  subroutine mpp_get_data_domain1D( domain, begin, end, size, max_size, is_global )
    type(domain1D),     intent(in) :: domain
    integer, intent(out), optional :: begin, end, size, max_size
    logical, intent(out), optional :: is_global

    if( PRESENT(begin)     )begin     = domain%data%begin
    if( PRESENT(end)       )end       = domain%data%end
    if( PRESENT(size)      )size      = domain%data%size
    if( PRESENT(max_size)  )max_size  = domain%data%max_size
    if( PRESENT(is_global) )is_global = domain%data%is_global
    return
  end subroutine mpp_get_data_domain1D

  !#####################################################################
  subroutine mpp_get_global_domain1D( domain, begin, end, size, max_size )
    type(domain1D),     intent(in) :: domain
    integer, intent(out), optional :: begin, end, size, max_size

    if( PRESENT(begin)    )begin    = domain%global%begin
    if( PRESENT(end)      )end      = domain%global%end
    if( PRESENT(size)     )size     = domain%global%size
    if( PRESENT(max_size) )max_size = domain%global%max_size
    return
  end subroutine mpp_get_global_domain1D

  !#####################################################################
  subroutine mpp_get_compute_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
       x_is_global, y_is_global )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    call mpp_get_compute_domain( domain%x, xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_compute_domain( domain%y, ybegin, yend, ysize, ymax_size, y_is_global )
    return
  end subroutine mpp_get_compute_domain2D

  !#####################################################################
  subroutine mpp_get_data_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
       x_is_global, y_is_global )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    call mpp_get_data_domain( domain%x, xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_data_domain( domain%y, ybegin, yend, ysize, ymax_size, y_is_global )
    return
  end subroutine mpp_get_data_domain2D

  !#####################################################################
  subroutine mpp_get_global_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    call mpp_get_global_domain( domain%x, xbegin, xend, xsize, xmax_size )
    call mpp_get_global_domain( domain%y, ybegin, yend, ysize, ymax_size )
    return
  end subroutine mpp_get_global_domain2D

  !#####################################################################
  ! <SUBROUTINE NAME="mpp_get_domain_components">
  !  <OVERVIEW>
  !    Retrieve 1D components of 2D decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    It is sometime necessary to have direct recourse to the domain1D types
  !    that compose a domain2D object. This call retrieves them.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_domain_components( domain, x, y )
  !  </TEMPLATE>
  !  <IN NAME="domain" TYPE="type(domain2D)"></IN>
  !  <OUT NAME="x,y"  TYPE="type(domain1D)"></OUT>
  ! </SUBROUTINE>
  subroutine mpp_get_domain_components( domain, x, y )
    type(domain2D),            intent(in) :: domain
    type(domain1D), intent(inout), optional :: x, y
    if( PRESENT(x) )x = domain%x
    if( PRESENT(y) )y = domain%y
    return
  end subroutine mpp_get_domain_components

  !#####################################################################
  subroutine mpp_get_compute_domains1D( domain, begin, end, size )
    type(domain1D),                   intent(in) :: domain
    integer, intent(out), optional, dimension(:) :: begin, end, size 

    if( .NOT.module_is_initialized ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: must first call mpp_domains_init.' )
    !we use shape instead of size for error checks because size is used as an argument
    if( PRESENT(begin) )then
       if( any(shape(begin).NE.shape(domain%list)) ) &
            call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: begin array size does not match domain.' )
       begin(:) = domain%list(:)%compute%begin
    end if
    if( PRESENT(end) )then
       if( any(shape(end).NE.shape(domain%list)) ) &
            call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: end array size does not match domain.' )
            end(:) = domain%list(:)%compute%end
    end if
    if( PRESENT(size) )then
       if( any(shape(size).NE.shape(domain%list)) ) &
           call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: size array size does not match domain.' )
           size(:) = domain%list(:)%compute%size
    end if
    return
end subroutine mpp_get_compute_domains1D


!#####################################################################
subroutine mpp_get_compute_domains2D( domain, xbegin, xend, xsize, ybegin, yend, ysize )
 type(domain2D),                   intent(in) :: domain
 integer, intent(out), optional, dimension(:) :: xbegin, xend, xsize, ybegin, yend, ysize

 if( .NOT.module_is_initialized ) &
      call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: must first call mpp_domains_init.' )

 if( PRESENT(xbegin) )then
    if( size(xbegin(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xbegin array size does not match domain.' )
    xbegin(:) = domain%list(:)%x%compute%begin
 end if
 if( PRESENT(xend) )then
    if( size(xend(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xend array size does not match domain.' )
    xend(:) = domain%list(:)%x%compute%end
 end if
 if( PRESENT(xsize) )then
    if( size(xsize(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xsize array size does not match domain.' )
    xsize(:) = domain%list(:)%x%compute%size
 end if
 if( PRESENT(ybegin) )then
    if( size(ybegin(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: ybegin array size does not match domain.' )
    ybegin(:) = domain%list(:)%y%compute%begin
 end if
 if( PRESENT(yend) )then
    if( size(yend(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: yend array size does not match domain.' )
    yend(:) = domain%list(:)%y%compute%end
 end if
 if( PRESENT(ysize) )then
    if( size(ysize(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: ysize array size does not match domain.' )
    ysize(:) = domain%list(:)%y%compute%size
 end if
 return
end subroutine mpp_get_compute_domains2D

!#####################################################################
! <SUBROUTINE NAME="mpp_get_pelist1D" INTERFACE="mpp_get_pelist">
!  <IN NAME="domain" TYPE="type(domain1D)"></IN>
!  <OUT NAME="pelist" TYPE="integer" DIM="(:)"></OUT>
!  <OUT NAME="pos" TYPE="integer"></OUT>
! </SUBROUTINE>
subroutine mpp_get_pelist1D( domain, pelist, pos )
 type(domain1D),     intent(in) :: domain
 integer,           intent(out) :: pelist(:)
 integer, intent(out), optional :: pos
 integer                        :: ndivs

 if( .NOT.module_is_initialized ) &
      call mpp_error( FATAL, 'MPP_GET_PELIST: must first call mpp_domains_init.' )
 ndivs = size(domain%list(:))

 if( size(pelist(:)).NE.ndivs ) &
      call mpp_error( FATAL, 'MPP_GET_PELIST: pelist array size does not match domain.' )

 pelist(:) = domain%list(0:ndivs-1)%pe
 if( PRESENT(pos) )pos = domain%pos
 return
end subroutine mpp_get_pelist1D

!#####################################################################
! <SUBROUTINE NAME="mpp_get_pelist2D" INTERFACE="mpp_get_pelist">
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <OUT NAME="pelist" TYPE="integer" DIM="(:)"></OUT>
!  <OUT NAME="pos" TYPE="integer"></OUT>
! </SUBROUTINE>
subroutine mpp_get_pelist2D( domain, pelist, pos )
 type(domain2D),     intent(in) :: domain
 integer,           intent(out) :: pelist(:)
 integer, intent(out), optional :: pos

 if( .NOT.module_is_initialized ) &
      call mpp_error( FATAL, 'MPP_GET_PELIST: must first call mpp_domains_init.' )
 if( size(pelist(:)).NE.size(domain%list(:)) ) &
      call mpp_error( FATAL, 'MPP_GET_PELIST: pelist array size does not match domain.' )

 pelist(:) = domain%list(:)%pe
 if( PRESENT(pos) )pos = domain%pos
 return
end subroutine mpp_get_pelist2D

!#####################################################################
! <SUBROUTINE NAME="mpp_get_layout1D" INTERFACE="mpp_get_layout">
!  <IN NAME="domain" TYPE="type(domain1D)"></IN>
!  <OUT NAME="layout" TYPE="integer"></OUT>
! </SUBROUTINE>
subroutine mpp_get_layout1D( domain, layout )
 type(domain1D), intent(in) :: domain
 integer,       intent(out) :: layout

 if( .NOT.module_is_initialized ) &
      call mpp_error( FATAL, 'MPP_GET_LAYOUT: must first call mpp_domains_init.' )

 layout = size(domain%list(:))
 return
end subroutine mpp_get_layout1D

!#####################################################################
! <SUBROUTINE NAME="mpp_get_layout2D" INTERFACE="mpp_get_layout">
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <OUT NAME="layout" TYPE="integer" DIM="(2)"></OUT>
! </SUBROUTINE>
subroutine mpp_get_layout2D( domain, layout )
 type(domain2D), intent(in) :: domain
 integer,       intent(out) :: layout(2)

 if( .NOT.module_is_initialized ) &
      call mpp_error( FATAL, 'MPP_GET_LAYOUT: must first call mpp_domains_init.' )

 layout(1) = size(domain%x%list(:))
 layout(2) = size(domain%y%list(:))
 return
end subroutine mpp_get_layout2D

!#####################################################################

    subroutine compute_overlaps( domain )
!computes remote domain overlaps
!assumes only one in each direction
      type(domain2D), intent(inout) :: domain
      integer :: i, j, k, m, n
      integer :: is, ie, js, je, isc, iec, jsc, jec, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ioff, joff
      integer :: list

      if( grid_offset_type.EQ.AGRID .AND. domain%remote_domains_initialized     )return
      if( grid_offset_type.NE.AGRID .AND. domain%remote_off_domains_initialized )return
      domain%gridtype = grid_offset_type
      n = size(domain%list(:))
!send
      call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=ioff, ysize=joff ) !cyclic offsets
      domain%list(:)%overlap = .FALSE.
      do list = 0,n-1
         m = mod( domain%pos+list, n )
!to_pe's eastern halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%compute%begin; je = domain%list(m)%y%compute%end
         if( ie.GT.ieg )then
             if( domain%x%cyclic .AND. iec.LT.is )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if    
         end if    
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_w_off%overlap = .TRUE.
                 domain%list(m)%send_w_off%is = is
                 domain%list(m)%send_w_off%ie = ie
                 domain%list(m)%send_w_off%js = js
                 domain%list(m)%send_w_off%je = je
             else
                 domain%list(m)%send_w%overlap = .TRUE.
                 domain%list(m)%send_w%is = is
                 domain%list(m)%send_w%ie = ie
                 domain%list(m)%send_w%js = js
                 domain%list(m)%send_w%je = je
             end if
         else    
             domain%list(m)%send_w%overlap = .FALSE.
             domain%list(m)%send_w_off%overlap = .FALSE.
         end if
!to_pe's SE halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( ie.GT.ieg )then
             if( domain%x%cyclic .AND. iec.LT.is )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if  
         end if  
         if( jsg.GT.js )then
             if( domain%y%cyclic .AND. je.LT.jsc )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if  
         end if  
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_nw_off%overlap = .TRUE.
                 domain%list(m)%send_nw_off%is = is
                 domain%list(m)%send_nw_off%ie = ie
                 domain%list(m)%send_nw_off%js = js
                 domain%list(m)%send_nw_off%je = je
             else
                 domain%list(m)%send_nw%overlap = .TRUE.
                 domain%list(m)%send_nw%is = is
                 domain%list(m)%send_nw%ie = ie
                 domain%list(m)%send_nw%js = js
                 domain%list(m)%send_nw%je = je
             end if
         else    
             domain%list(m)%send_nw%overlap = .FALSE.
             domain%list(m)%send_nw_off%overlap = .FALSE.
         end if
!to_pe's southern halo
         is = domain%list(m)%x%compute%begin; ie = domain%list(m)%x%compute%end
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( jsg.GT.js )then
             if( domain%y%cyclic .AND. je.LT.jsc )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if  
         end if  
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_n_off%overlap = .TRUE.
                 domain%list(m)%send_n_off%is = is
                 domain%list(m)%send_n_off%ie = ie
                 domain%list(m)%send_n_off%js = js
                 domain%list(m)%send_n_off%je = je
             else
                 domain%list(m)%send_n%overlap = .TRUE.
                 domain%list(m)%send_n%is = is
                 domain%list(m)%send_n%ie = ie
                 domain%list(m)%send_n%js = js
                 domain%list(m)%send_n%je = je
             end if
         else    
             domain%list(m)%send_n%overlap = .FALSE.
             domain%list(m)%send_n_off%overlap = .FALSE.
         end if
!to_pe's SW halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( isg.GT.is )then
             if( domain%x%cyclic .AND. ie.LT.isc )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if  
         end if  
         if( jsg.GT.js )then
             if( domain%y%cyclic .AND. je.LT.jsc )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if  
         end if  
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_ne_off%overlap = .TRUE.
                 domain%list(m)%send_ne_off%is = is
                 domain%list(m)%send_ne_off%ie = ie
                 domain%list(m)%send_ne_off%js = js
                 domain%list(m)%send_ne_off%je = je
             else
                 domain%list(m)%send_ne%overlap = .TRUE.
                 domain%list(m)%send_ne%is = is
                 domain%list(m)%send_ne%ie = ie
                 domain%list(m)%send_ne%js = js
                 domain%list(m)%send_ne%je = je
             end if
         else    
             domain%list(m)%send_ne%overlap = .FALSE.
             domain%list(m)%send_ne_off%overlap = .FALSE.
         end if
!to_pe's western halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%compute%begin; je = domain%list(m)%y%compute%end
         if( isg.GT.is )then
             if( domain%x%cyclic .AND. ie.LT.isc )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if  
         end if  
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_e_off%overlap = .TRUE.
                 domain%list(m)%send_e_off%is = is
                 domain%list(m)%send_e_off%ie = ie
                 domain%list(m)%send_e_off%js = js
                 domain%list(m)%send_e_off%je = je
             else
                 domain%list(m)%send_e%overlap = .TRUE.
                 domain%list(m)%send_e%is = is
                 domain%list(m)%send_e%ie = ie
                 domain%list(m)%send_e%js = js
                 domain%list(m)%send_e%je = je
             end if
         else    
             domain%list(m)%send_e%overlap = .FALSE.
             domain%list(m)%send_e_off%overlap = .FALSE.
         end if
!to_pe's NW halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( isg.GT.is )then
             if( domain%x%cyclic .AND. ie.LT.isc )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if  
         end if  
         if( je.GT.jeg )then
             if( domain%y%cyclic .AND. jec.LT.js )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if  
         end if  
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_se_off%overlap = .TRUE.
                 domain%list(m)%send_se_off%is = is
                 domain%list(m)%send_se_off%ie = ie
                 domain%list(m)%send_se_off%js = js
                 domain%list(m)%send_se_off%je = je
             else
                 domain%list(m)%send_se%overlap = .TRUE.
                 domain%list(m)%send_se%is = is
                 domain%list(m)%send_se%ie = ie
                 domain%list(m)%send_se%js = js
                 domain%list(m)%send_se%je = je
             end if
         else    
             domain%list(m)%send_se%overlap = .FALSE.
             domain%list(m)%send_se_off%overlap = .FALSE.
         end if
!to_pe's northern halo
         is = domain%list(m)%x%compute%begin; ie = domain%list(m)%x%compute%end
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( je.GT.jeg )then
             if( domain%y%cyclic .AND. jec.LT.js )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if  
         end if  
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_s_off%overlap = .TRUE.
                 domain%list(m)%send_s_off%is = is
                 domain%list(m)%send_s_off%ie = ie
                 domain%list(m)%send_s_off%js = js
                 domain%list(m)%send_s_off%je = je
             else
                 domain%list(m)%send_s%overlap = .TRUE.
                 domain%list(m)%send_s%is = is
                 domain%list(m)%send_s%ie = ie
                 domain%list(m)%send_s%js = js
                 domain%list(m)%send_s%je = je
             end if
         else    
             domain%list(m)%send_s%overlap = .FALSE.
             domain%list(m)%send_s_off%overlap = .FALSE.
         end if
!to_pe's NE halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( ie.GT.ieg )then
             if( domain%x%cyclic .AND. iec.LT.is )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if
         end if  
         if( je.GT.jeg )then
             if( domain%y%cyclic .AND. jec.LT.js )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if  
         end if  
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_sw_off%overlap = .TRUE.
                 domain%list(m)%send_sw_off%is = is
                 domain%list(m)%send_sw_off%ie = ie
                 domain%list(m)%send_sw_off%js = js
                 domain%list(m)%send_sw_off%je = je
             else
                 domain%list(m)%send_sw%overlap = .TRUE.
                 domain%list(m)%send_sw%is = is
                 domain%list(m)%send_sw%ie = ie
                 domain%list(m)%send_sw%js = js
                 domain%list(m)%send_sw%je = je
             end if
         else    
             domain%list(m)%send_sw%overlap = .FALSE.
             domain%list(m)%send_sw_off%overlap = .FALSE.
         end if
      end do 
            
!recv       
      do list = 0,n-1
         m = mod( domain%pos+n-list, n )
         call mpp_get_compute_domain( domain%list(m), isc, iec, jsc, jec )
!recv_e  
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%compute%begin; jed = domain%y%compute%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_e%folded = .FALSE.
         if( ied.GT.ieg )then
             if( domain%x%cyclic .AND. ie.LT.isd )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 domain%list(m)%recv_e%folded = .TRUE.
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if  
         end if  
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_e_off%overlap = .TRUE.
                 domain%list(m)%recv_e_off%is = is
                 domain%list(m)%recv_e_off%ie = ie
                 domain%list(m)%recv_e_off%js = js
                 domain%list(m)%recv_e_off%je = je
             else
                 domain%list(m)%recv_e%overlap = .TRUE.
                 domain%list(m)%recv_e%is = is
                 domain%list(m)%recv_e%ie = ie
                 domain%list(m)%recv_e%js = js
                 domain%list(m)%recv_e%je = je
             endif
         else    
             domain%list(m)%recv_e%overlap = .FALSE.
             domain%list(m)%recv_e_off%overlap = .FALSE.
         end if
!recv_se     
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_se%folded = .FALSE.
         if( jsd.LT.jsg )then
             if( domain%y%cyclic .AND. js.GT.jed )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 domain%list(m)%recv_se%folded = .TRUE.
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if  
         end if  
         if( ied.GT.ieg )then
             if( domain%x%cyclic .AND. ie.LT.isd )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_se%folded = .TRUE.
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if  
         end if  
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_se_off%overlap = .TRUE.
                 domain%list(m)%recv_se_off%is = is
                 domain%list(m)%recv_se_off%ie = ie
                 domain%list(m)%recv_se_off%js = js
                 domain%list(m)%recv_se_off%je = je
             else
                 domain%list(m)%recv_se%overlap = .TRUE.
                 domain%list(m)%recv_se%is = is
                 domain%list(m)%recv_se%ie = ie
                 domain%list(m)%recv_se%js = js
                 domain%list(m)%recv_se%je = je
             endif
         else    
             domain%list(m)%recv_se%overlap = .FALSE.
             domain%list(m)%recv_se_off%overlap = .FALSE.
         end if
!recv_s      
         isd = domain%x%compute%begin; ied = domain%x%compute%end
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_s%folded = .FALSE.
         if( jsd.LT.jsg )then
             if( domain%y%cyclic .AND. js.GT.jed )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 domain%list(m)%recv_s%folded = .TRUE.
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if  
         end if  
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_s_off%overlap = .TRUE.
                 domain%list(m)%recv_s_off%is = is
                 domain%list(m)%recv_s_off%ie = ie
                 domain%list(m)%recv_s_off%js = js
                 domain%list(m)%recv_s_off%je = je
             else
                 domain%list(m)%recv_s%overlap = .TRUE.
                 domain%list(m)%recv_s%is = is
                 domain%list(m)%recv_s%ie = ie
                 domain%list(m)%recv_s%js = js
                 domain%list(m)%recv_s%je = je
             endif
         else    
             domain%list(m)%recv_s%overlap = .FALSE.
             domain%list(m)%recv_s_off%overlap = .FALSE.
             
         end if
!recv_sw 
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_sw%folded = .FALSE.
         if( jsd.LT.jsg )then
             if( domain%y%cyclic .AND. js.GT.jed )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 domain%list(m)%recv_sw%folded = .TRUE.
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if  
         end if  
         if( isd.LT.isg )then
             if( domain%x%cyclic .AND. is.GT.ied )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_sw%folded = .TRUE.
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if  
         end if  
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_sw_off%overlap = .TRUE.
                 domain%list(m)%recv_sw_off%is = is
                 domain%list(m)%recv_sw_off%ie = ie
                 domain%list(m)%recv_sw_off%js = js
                 domain%list(m)%recv_sw_off%je = je
             else
                 domain%list(m)%recv_sw%overlap = .TRUE.
                 domain%list(m)%recv_sw%is = is
                 domain%list(m)%recv_sw%ie = ie
                 domain%list(m)%recv_sw%js = js
                 domain%list(m)%recv_sw%je = je
             endif
         else    
             domain%list(m)%recv_sw%overlap = .FALSE.
             domain%list(m)%recv_sw_off%overlap = .FALSE.
         end if
!recv_w      
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%compute%begin; jed = domain%y%compute%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_w%folded = .FALSE.
         if( isd.LT.isg )then
             if( domain%x%cyclic .AND. is.GT.ied )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_w%folded = .TRUE.
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if  
         end if  
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_w_off%overlap = .TRUE.
                 domain%list(m)%recv_w_off%is = is
                 domain%list(m)%recv_w_off%ie = ie
                 domain%list(m)%recv_w_off%js = js
                 domain%list(m)%recv_w_off%je = je
             else
                 domain%list(m)%recv_w%overlap = .TRUE.
                 domain%list(m)%recv_w%is = is
                 domain%list(m)%recv_w%ie = ie
                 domain%list(m)%recv_w%js = js
                 domain%list(m)%recv_w%je = je
             endif
         else    
             domain%list(m)%recv_w%overlap = .FALSE.
             domain%list(m)%recv_w_off%overlap = .FALSE.
         end if
!recv_nw     
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_nw%folded = .FALSE.
         if( jed.GT.jeg )then
             if( domain%y%cyclic .AND. je.LT.jsd )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 domain%list(m)%recv_nw%folded = .TRUE.
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if  
         end if  
         if( isd.LT.isg )then
             if( domain%x%cyclic .AND. is.GT.ied )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_nw%folded = .TRUE.
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if  
         end if  
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_nw_off%overlap = .TRUE.
                 domain%list(m)%recv_nw_off%is = is
                 domain%list(m)%recv_nw_off%ie = ie
                 domain%list(m)%recv_nw_off%js = js
                 domain%list(m)%recv_nw_off%je = je
             else
                 domain%list(m)%recv_nw%overlap = .TRUE.
                 domain%list(m)%recv_nw%is = is
                 domain%list(m)%recv_nw%ie = ie
                 domain%list(m)%recv_nw%js = js
                 domain%list(m)%recv_nw%je = je
             endif
         else    
             domain%list(m)%recv_nw%overlap = .FALSE.
             domain%list(m)%recv_nw_off%overlap = .FALSE.
         end if
!recv_n      
         isd = domain%x%compute%begin; ied = domain%x%compute%end
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_n%folded = .FALSE.
         if( jed.GT.jeg )then
             if( domain%y%cyclic .AND. je.LT.jsd )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 domain%list(m)%recv_n%folded = .TRUE.
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if  
         end if  
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_n_off%overlap = .TRUE.
                 domain%list(m)%recv_n_off%is = is
                 domain%list(m)%recv_n_off%ie = ie
                 domain%list(m)%recv_n_off%js = js
                 domain%list(m)%recv_n_off%je = je
             else
                 domain%list(m)%recv_n%overlap = .TRUE.
                 domain%list(m)%recv_n%is = is
                 domain%list(m)%recv_n%ie = ie
                 domain%list(m)%recv_n%js = js
                 domain%list(m)%recv_n%je = je
             end if
         else    
             domain%list(m)%recv_n%overlap = .FALSE.
             domain%list(m)%recv_n_off%overlap = .FALSE.
         end if
!recv_ne     
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_ne%folded = .FALSE.
         if( jed.GT.jeg )then
             if( domain%y%cyclic .AND. je.LT.jsd )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 domain%list(m)%recv_ne%folded = .TRUE.
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if  
         end if  
         if( ied.GT.ieg )then
             if( domain%x%cyclic .AND. ie.LT.isd )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_ne%folded = .TRUE.
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if  
         end if  
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_ne_off%overlap = .TRUE.
                 domain%list(m)%recv_ne_off%is = is
                 domain%list(m)%recv_ne_off%ie = ie
                 domain%list(m)%recv_ne_off%js = js
                 domain%list(m)%recv_ne_off%je = je
             else
                 domain%list(m)%recv_ne%overlap = .TRUE.
                 domain%list(m)%recv_ne%is = is
                 domain%list(m)%recv_ne%ie = ie
                 domain%list(m)%recv_ne%js = js
                 domain%list(m)%recv_ne%je = je
             end if
         else    
             domain%list(m)%recv_ne%overlap = .FALSE.
             domain%list(m)%recv_ne_off%overlap = .FALSE.
         end if
      end do 
      if( grid_offset_type.EQ.AGRID )domain%remote_domains_initialized = .TRUE.
      if( grid_offset_type.NE.AGRID )domain%remote_off_domains_initialized = .TRUE.
      return
    end subroutine compute_overlaps

!#######################################################################

end module mpp_domains_util_mod
