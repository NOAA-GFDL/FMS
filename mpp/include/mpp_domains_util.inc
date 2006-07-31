
! -*-f90-*-
! $Id: mpp_domains_util.inc,v 13.0.2.3 2006/05/05 19:51:14 z1l Exp $

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
    if( mpp_pe().EQ.mpp_root_pe() )call mpp_error( NOTE, 'MPP_DOMAINS_SET_STACK_SIZE: stack size set to '//text//'.' )

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
    integer                    :: nt, n
   
    mpp_domain2d_eq = size(a%x(:)) .EQ. size(b%x(:))
    nt = size(a%x(:))
    do n = 1, nt
       if(mpp_domain2d_eq) mpp_domain2D_eq = a%x(n).EQ.b%x(n) .AND. a%y(n).EQ.b%y(n)
    end do

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
  subroutine mpp_get_memory_domain1D( domain, begin, end, size, max_size, is_global )
    type(domain1D),     intent(in) :: domain
    integer, intent(out), optional :: begin, end, size, max_size
    logical, intent(out), optional :: is_global

    if( PRESENT(begin)     )begin     = domain%memory%begin
    if( PRESENT(end)       )end       = domain%memory%end
    if( PRESENT(size)      )size      = domain%memory%size
    if( PRESENT(max_size)  )max_size  = domain%memory%max_size
    if( PRESENT(is_global) )is_global = domain%memory%is_global
    return
  end subroutine mpp_get_memory_domain1D

  !#####################################################################
  subroutine mpp_get_memory_compute_domain1D( domain, begin, end, size, max_size, is_global )
    type(domain1D),     intent(in) :: domain
    integer, intent(out), optional :: begin, end, size, max_size
    logical, intent(out), optional :: is_global

    if( PRESENT(begin)     )begin     = domain%memory%begin + domain%compute%begin - domain%data%begin
    if( PRESENT(end)       )end       = domain%memory%begin + domain%compute%end - domain%data%begin
    if( PRESENT(size)      )size      = domain%compute%size
    if( PRESENT(max_size)  )max_size  = domain%compute%max_size
    if( PRESENT(is_global) )is_global = domain%compute%is_global
    return
  end subroutine mpp_get_memory_compute_domain1D

  !#####################################################################
  subroutine mpp_get_memory_data_domain1D( domain, begin, end, size, max_size, is_global )
    type(domain1D),     intent(in) :: domain
    integer, intent(out), optional :: begin, end, size, max_size
    logical, intent(out), optional :: is_global

    if( PRESENT(begin)     )begin     = domain%memory%begin
    if( PRESENT(end)       )end       = domain%memory%begin + domain%data%end - domain%data%begin
    if( PRESENT(size)      )size      = domain%data%size
    if( PRESENT(max_size)  )max_size  = domain%data%max_size
    if( PRESENT(is_global) )is_global = domain%data%is_global
    return
  end subroutine mpp_get_memory_data_domain1D

  !#####################################################################
  subroutine mpp_get_compute_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
       x_is_global, y_is_global, tile_number )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_number
    integer                        :: tile

    tile = 1
    if(present(tile_number)) tile = tile_number

    call mpp_get_compute_domain( domain%x(tile), xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_compute_domain( domain%y(tile), ybegin, yend, ysize, ymax_size, y_is_global )

    return
  end subroutine mpp_get_compute_domain2D

  !#####################################################################
  subroutine mpp_get_data_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
       x_is_global, y_is_global, tile_number )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_number
    integer                        :: tile

    tile = 1
    if(present(tile_number)) tile = tile_number

    call mpp_get_data_domain( domain%x(tile), xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_data_domain( domain%y(tile), ybegin, yend, ysize, ymax_size, y_is_global )

    return
  end subroutine mpp_get_data_domain2D

  !#####################################################################
  subroutine mpp_get_global_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, tile_number )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    integer, intent(in),  optional :: tile_number
    integer                        :: tile

    tile = 1
    if(present(tile_number)) tile = tile_number

    call mpp_get_global_domain( domain%x(tile), xbegin, xend, xsize, xmax_size )
    call mpp_get_global_domain( domain%y(tile), ybegin, yend, ysize, ymax_size )

    return
  end subroutine mpp_get_global_domain2D

  !#####################################################################
  subroutine mpp_get_memory_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
       x_is_global, y_is_global)
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    integer                        :: tile

    tile = 1

    call mpp_get_memory_domain( domain%x(tile), xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_memory_domain( domain%y(tile), ybegin, yend, ysize, ymax_size, y_is_global )

    return
  end subroutine mpp_get_memory_domain2D

  !#####################################################################
  subroutine mpp_get_memory_compute_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, &
       ysize, ymax_size, x_is_global, y_is_global, tile_number )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_number
    integer                        :: tile

    tile = 1
    if(present(tile_number)) tile = tile_number

    call mpp_get_memory_compute_domain( domain%x(tile), xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_memory_compute_domain( domain%y(tile), ybegin, yend, ysize, ymax_size, y_is_global )

    return
  end subroutine mpp_get_memory_compute_domain2D

  !#####################################################################
  subroutine mpp_get_memory_data_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, &
       ysize, ymax_size, x_is_global, y_is_global, tile_number )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_number
    integer                        :: tile

    tile = 1
    if(present(tile_number)) tile = tile_number

    call mpp_get_memory_data_domain( domain%x(tile), xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_memory_data_domain( domain%y(tile), ybegin, yend, ysize, ymax_size, y_is_global )

    return
  end subroutine mpp_get_memory_data_domain2D

  !#####################################################################
  subroutine mpp_set_compute_domain1D( domain, begin, end, size, is_global )
    type(domain1D), intent(inout) :: domain
    integer, intent(in), optional :: begin, end, size
    logical, intent(in), optional :: is_global

    if(present(begin)) domain%compute%begin = begin
    if(present(end))   domain%compute%end   = end
    if(present(size))  domain%compute%size  = size
    if(present(is_global)) domain%compute%is_global = is_global

  end subroutine mpp_set_compute_domain1D

  !#####################################################################
  subroutine mpp_set_compute_domain2D( domain, xbegin, xend, ybegin, yend, xsize, ysize, &
                                       x_is_global, y_is_global, tile_number )
    type(domain2D), intent(inout) :: domain
    integer, intent(in), optional :: xbegin, xend, ybegin, yend, xsize, ysize
    logical, intent(in), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_number
    integer                        :: tile

    tile = 1
    if(present(tile_number)) tile = tile_number

    call mpp_set_compute_domain(domain%x(tile), xbegin, xend, xsize, x_is_global)
    call mpp_set_compute_domain(domain%y(tile), ybegin, yend, ysize, y_is_global)

  end subroutine mpp_set_compute_domain2D

  !#####################################################################
  subroutine mpp_set_data_domain1D( domain, begin, end, size, is_global )
    type(domain1D), intent(inout) :: domain
    integer, intent(in), optional :: begin, end, size
    logical, intent(in), optional :: is_global

    if(present(begin)) domain%data%begin = begin
    if(present(end))   domain%data%end   = end
    if(present(size))  domain%data%size  = size
    if(present(is_global)) domain%data%is_global = is_global

  end subroutine mpp_set_data_domain1D

  !#####################################################################
  subroutine mpp_set_data_domain2D( domain, xbegin, xend, ybegin, yend, xsize, ysize, &
                                    x_is_global, y_is_global, tile_number )
    type(domain2D), intent(inout) :: domain
    integer, intent(in), optional :: xbegin, xend, ybegin, yend, xsize, ysize
    logical, intent(in), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_number
    integer                        :: tile

    tile = 1
    if(present(tile_number)) tile = tile_number

    call mpp_set_data_domain(domain%x(tile), xbegin, xend, xsize, x_is_global)
    call mpp_set_data_domain(domain%y(tile), ybegin, yend, ysize, y_is_global)

  end subroutine mpp_set_data_domain2D

  !#####################################################################
  subroutine mpp_set_global_domain1D( domain, begin, end, size)
    type(domain1D), intent(inout) :: domain
    integer, intent(in), optional :: begin, end, size

    if(present(begin)) domain%global%begin = begin
    if(present(end))   domain%global%end   = end
    if(present(size))  domain%global%size  = size

  end subroutine mpp_set_global_domain1D

  !#####################################################################
  subroutine mpp_set_global_domain2D( domain, xbegin, xend, ybegin, yend, xsize, ysize, tile_number )
    type(domain2D), intent(inout) :: domain
    integer, intent(in), optional :: xbegin, xend, ybegin, yend, xsize, ysize
    integer, intent(in),  optional :: tile_number
    integer                        :: tile

    tile = 1
    if(present(tile_number)) tile = tile_number
    call mpp_set_global_domain(domain%x(tile), xbegin, xend, xsize)
    call mpp_set_global_domain(domain%y(tile), ybegin, yend, ysize)

  end subroutine mpp_set_global_domain2D

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
  subroutine mpp_get_domain_components( domain, x, y, tile_number )
    type(domain2D),            intent(in) :: domain
    type(domain1D), intent(inout), optional :: x, y
    integer, intent(in),  optional :: tile_number
    integer                        :: tile

    tile = 1
    if(present(tile_number)) tile = tile_number
    if( PRESENT(x) )x = domain%x(tile)
    if( PRESENT(y) )y = domain%y(tile)
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
integer :: i
 if( .NOT.module_is_initialized ) &
      call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: must first call mpp_domains_init.' )

 if( PRESENT(xbegin) )then
    if( size(xbegin(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xbegin array size does not match domain.' )
    do i = 1, size(xbegin(:))
       xbegin(i) = domain%list(i-1)%x(1)%compute%begin
    end do
 end if
 if( PRESENT(xend) )then
    if( size(xend(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xend array size does not match domain.' )
    do i = 1, size(xend(:))
       xend(i) = domain%list(i-1)%x(1)%compute%end
    end do
 end if
 if( PRESENT(xsize) )then
    if( size(xsize(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xsize array size does not match domain.' )
    do i = 1, size(xsize(:))
       xsize(i) = domain%list(i-1)%x(1)%compute%size
    end do
 end if
 if( PRESENT(ybegin) )then
    if( size(ybegin(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: ybegin array size does not match domain.' )
    do i = 1, size(ybegin(:))
       ybegin(i) = domain%list(i-1)%y(1)%compute%begin
    end do
 end if
 if( PRESENT(yend) )then
    if( size(yend(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: yend array size does not match domain.' )
    do i = 1, size(yend(:))
       yend(i) = domain%list(i-1)%y(1)%compute%end
    end do
 end if
 if( PRESENT(ysize) )then
    if( size(ysize(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: ysize array size does not match domain.' )
    do i = 1, size(ysize(:))
       ysize(i) = domain%list(i-1)%y(1)%compute%size
    end do
 end if
 return
end subroutine mpp_get_compute_domains2D


!#####################################################################
! This will return xextent and yextent for each tile
subroutine mpp_get_domain_extents(domain, xextent, yextent)
   type(domain2d),             intent(in) :: domain
   integer, dimension(:,:), intent(inout) :: xextent, yextent
   integer                                :: ntile, nlist, n, m, ndivx, ndivy, tile, pos

   ntile = domain%ntiles
   nlist = size(domain%list(:))
   if(size(xextent,2) .ne. ntile .or. size(yextent,2) .ne. ntile) call mpp_error(FATAL, &
       "mpp_domains_utile.inc: the second dimension size of xextent/yextent is not correct")
   ndivx = size(xextent,1); ndivy = size(yextent,1)
   do n = 0, nlist-1
      if(ANY(domain%list(n)%x(:)%pos>ndivx-1) ) call mpp_error(FATAL, &
         "mpp_domains_utile.inc: first dimension size of xextent is less than the x-layout in some tile")
      if(ANY(domain%list(n)%y(:)%pos>ndivy-1) ) call mpp_error(FATAL, &
         "mpp_domains_utile.inc: first dimension size of yextent is less than the y-layout in some tile")
   end do

   xextent = 0; yextent=0

   do n = 0, nlist-1
      do m = 1, size(domain%list(n)%tile_id(:))
         tile = domain%list(n)%tile_id(m)
         pos = domain%list(n)%x(m)%pos+1
         if(xextent(pos, tile) == 0) xextent(pos,tile) = domain%list(n)%x(m)%compute%size
         pos = domain%list(n)%y(m)%pos+1
         if(yextent(pos, tile) == 0) yextent(pos,tile) = domain%list(n)%y(m)%compute%size
      end do
   end do


end subroutine mpp_get_domain_extents

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

 layout(1) = size(domain%x(1)%list(:))
 layout(2) = size(domain%y(1)%list(:))
 return
end subroutine mpp_get_layout2D

!#####################################################################
  ! <SUBROUTINE NAME="mpp_get_domain_shift">
  !  <OVERVIEW>
  !    Returns the shift value in x and y-direction according to domain position.. 
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    When domain is symmetry, one extra point maybe needed in
  !    x- and/or y-direction. This routine will return the shift value based
  !    on the position
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_domain_shift( domain, ishift, jshift, position )
  !  </TEMPLATE>
  !  <IN NAME="domain" TYPE="type(domain2D)">
  !    predefined data contains 2-d domain decomposition.
  !  </IN>
  !  <OUT NAME="ishift, jshift"  TYPE="integer">
  !    return value will be 0 or 1.
  !  </OUT>
  !  <IN NAME="position" TYPE="integer">
  !   position of data. Its value can be CENTER, EAST, NORTH or CORNER.
  !  </OUT>
  ! </SUBROUTINE>
subroutine mpp_get_domain_shift(domain, ishift, jshift, position)
  type(domain2D),     intent(in) :: domain
  integer,           intent(out) :: ishift, jshift
  integer, optional,  intent(in) :: position
  integer :: pos

  ishift = 0 ; jshift = 0
  pos = CENTER
  if(present(position)) pos = position

  if(domain%symmetry) then ! shift is non-zero only when the domain is symmetry.
     select case(pos)
     case(CORNER)
        ishift = 1; jshift = 1
     case(EAST)
        ishift = 1
     case(NORTH)
        jshift = 1
     end select
  end if

end subroutine mpp_get_domain_shift

!#####################################################################

    subroutine mpp_get_neighbor_pe_1d(domain, direction, pe)

      ! Return PE to the righ/left of this PE-domain.    
  
      type(domain1D), intent(inout) :: domain
      integer, intent(in)           :: direction
      integer, intent(out)          :: pe
      
      integer ipos, ipos2, npx

      pe   = NULL_PE
      npx  = size(domain%list(:)) ! 0..npx-1
      ipos = domain%pos

      select case (direction)

      case (:-1)
         ! neighbor on the left
         ipos2 = ipos - 1
         if(ipos2 <    0) then
            if(domain%cyclic) then 
               ipos2 = npx-1
            else
               ipos2 = -999
            endif
         endif
             
      case (0)
         ! identity
         ipos2 = ipos

      case (1:)
         ! neighbor on the right
         ipos2 = ipos + 1
         if(ipos2 > npx-1) then
            if(domain%cyclic) then
               ipos2 = 0
            else
               ipos2 = -999
            endif
         endif

      end select

      if(ipos2 >= 0) pe = domain%list(ipos2)%pe
         
    end subroutine mpp_get_neighbor_pe_1d
!#####################################################################

    subroutine mpp_get_neighbor_pe_2d(domain, direction, pe)

      ! Return PE North/South/East/West of this PE-domain.
      ! direction must be NORTH, SOUTH, EAST or WEST.

      type(domain2D), intent(inout) :: domain
      integer, intent(in)           :: direction
      integer, intent(out)          :: pe

      integer ipos, jpos, npx, npy, ix, iy, ipos0, jpos0

      pe   = NULL_PE
      npx  = size(domain%x(1)%list(:)) ! 0..npx-1
      npy  = size(domain%y(1)%list(:)) ! 0..npy-1
      ipos0 = domain%x(1)%pos
      jpos0 = domain%y(1)%pos

      select case (direction)
      case (NORTH)
         ix = 0
         iy = 1
      case (NORTH_EAST)
         ix = 1
         iy = 1
      case (EAST)
         ix = 1
         iy = 0
      case (SOUTH_EAST)
         ix = 1
         iy =-1
      case (SOUTH)
         ix = 0
         iy =-1
      case (SOUTH_WEST)
         ix =-1
         iy =-1
      case (WEST)
         ix =-1
         iy = 0
      case (NORTH_WEST)
         ix =-1
         iy = 1

      case default
         call mpp_error( FATAL, &
 & 'MPP_GET_NEIGHBOR_PE_2D: direction must be either NORTH, ' &
 & // 'SOUTH, EAST, WEST, NORTH_EAST, SOUTH_EAST, SOUTH_WEST or NORTH_WEST')
      end select

      ipos = ipos0 + ix
      jpos = jpos0 + iy

      
      if( (ipos < 0 .or. ipos > npx-1) .and. domain%x(1)%cyclic ) then
         ! E/W cyclic domain
         ipos = modulo(ipos, npx)
      endif

      if(    (ipos < 0     .and. btest(domain%fold,WEST)) .or. &
           & (ipos > npx-1 .and. btest(domain%fold,EAST)) ) then  
         ! E or W folded domain
           ipos = ipos0
           jpos = npy-jpos-1
      endif

      if( (jpos < 0 .or. jpos > npy-1) .and. domain%y(1)%cyclic ) then
         ! N/S cyclic
         jpos = modulo(jpos, npy)
      endif

      if(    (jpos < 0     .and. btest(domain%fold,SOUTH)) .or. &
           & (jpos > npy-1 .and. btest(domain%fold,NORTH)) ) then         
         ! N or S folded
           ipos = npx-ipos-1
           jpos = jpos0
      endif

      ! get the PE number
      pe = NULL_PE
      if(ipos >= 0 .and. ipos <= npx-1 .and. jpos >= 0 .and. jpos <= npy-1) then
         pe = domain%pearray(ipos, jpos)
      endif


    end subroutine mpp_get_neighbor_pe_2d
      

!#######################################################################

  subroutine nullify_domain2d_list(domain)
     type(domain2d), intent(inout) :: domain

     domain%list =>NULL()

  end subroutine nullify_domain2d_list

!#######################################################################
  function mpp_domain_is_symmetry(domain)
    type(domain2d), intent(in) :: domain
    logical                    :: mpp_domain_is_symmetry

    mpp_domain_is_symmetry = domain%symmetry
    return

  end function mpp_domain_is_symmetry

!#######################################################################
  !--- private routine used only for mpp_update_domains. This routine will
  !--- compare whalo, ehalo, shalo, nhalo with the halo size when defining "domain"
  !--- to decide if update is needed. Also it check the sign of whalo, ehalo, shalo and nhalo.
  function domain_update_is_needed(domain, whalo, ehalo, shalo, nhalo)
    type(domain2d), intent(in) :: domain
    integer,        intent(in) :: whalo, ehalo, shalo, nhalo
    logical                    :: domain_update_is_needed

    domain_update_is_needed = .true.

    if(whalo == 0 .AND. ehalo==0 .AND. shalo == 0 .AND. nhalo==0 ) then
       domain_update_is_needed = .false.
       call mpp_error(NOTE, 'mpp_domains_util.inc: halo size to be updated are all zero, no update will be done')
       return
    end if
    if(  (whalo == -domain%whalo .AND. domain%whalo .NE. 0) .or.  & 
         (ehalo == -domain%ehalo .AND. domain%ehalo .NE. 0) .or.  &
         (shalo == -domain%shalo .AND. domain%shalo .NE. 0) .or.  & 
         (nhalo == -domain%nhalo .AND. domain%nhalo .NE. 0) ) then
       domain_update_is_needed = .false.
       call mpp_error(NOTE, 'mpp_domains_util.inc: at least one of w/e/s/n halo size to be updated '// &
            'is the inverse of the original halo when defining domain, no update will be done')
       return
    end if

  end function domain_update_is_needed
!#######################################################################
  ! this routine found the domain has the same halo size with the input
  ! whalo, ehalo, 
  function search_domain(domain, whalo, ehalo, shalo, nhalo, position, gridtype, direction)
    type(domain2d),          intent(inout) :: domain
    integer,                    intent(in) :: whalo, ehalo, shalo, nhalo
    integer, optional,          intent(in) :: position, gridtype
    character(len=1), optional, intent(in) :: direction
    type(domain2d),             pointer    :: search_domain
    type(domain2d),             pointer    :: Dom_ref => NULL()
    integer                                :: update_position

    ! select the domain
    if(present(position)) then
       update_position = position
    else if(present(gridtype)) then
       if(.NOT. present(direction)) call mpp_error(FATAL,  &
            "mpp_omains_util.inc: when gridtype is present, direction must present.")
       select case(gridtype)
       case(AGRID)
          update_position = CENTER
       case(BGRID_NE, BGRID_SW)
          update_position = CORNER
       case(CGRID_NE, CGRID_SW)
          if(direction == 'x' ) then
             update_position = EAST
          else
             update_position = NORTH
          end if
       case default
          call mpp_error(FATAL,"mpp_domains_util.inc: gridtype should be AGRID|BGRID_NE|BGRID_SW|CGRID_NE|CGRID_SW")
       end select
    else
       update_position = CENTER
    end if


    select case(update_position)
    case (CENTER)
       Dom_ref => domain%T
    case (CORNER)
       Dom_ref => domain%C
    case (NORTH)
       Dom_ref => domain%N
    case (EAST)
       Dom_ref => domain%E
    case default
       call mpp_error(FATAL,"mpp_domains_util.inc: position should be CENTER|CORNER|EAST|NORTH")
    end select

    if( whalo == Dom_ref%whalo .AND. ehalo == Dom_ref%ehalo .AND. shalo == Dom_ref%shalo &
         .AND. nhalo == Dom_ref%nhalo ) then
       search_domain => Dom_ref
    else     
       search_domain => domain%next

       do 
          if(.NOT. ASSOCIATED(search_domain)) then
             allocate(search_domain)
             search_domain%initialized = .false.
             call set_overlaps(Dom_ref, search_domain, whalo, ehalo, shalo, nhalo )
          end if

          if(whalo == search_domain%whalo .AND. ehalo == search_domain%ehalo .AND. &
               shalo == search_domain%shalo .AND. nhalo == search_domain%nhalo ) exit ! found domain

          search_domain => search_domain%next
       end do
    end if

  end function search_domain

  !########################################################################
  ! return the domain list correponding to position
  function get_domain(domain, position)
     type(domain2d),    intent(in) :: domain
     integer, optional, intent(in) :: position
     integer                       :: update_position
     type(domain2d), pointer       :: get_domain

     update_position = CENTER
     if(present(position)) update_position = position
     select case (update_position)
     case ( CENTER )
        get_domain => domain%T
     case ( CORNER )
        get_domain => domain%C
     case ( EAST )
        get_domain => domain%E
     case ( NORTH )
        get_domain => domain%N     
     end select

     return

  end function get_domain

  !########################################################################
  ! return the tile_id on current pe
  function mpp_get_tile_id(domain)
     type(domain2d),                  intent(in) :: domain
     integer, dimension(size(domain%tile_id(:))) :: mpp_get_tile_id

     mpp_get_tile_id = domain%tile_id
     return

  end function mpp_get_tile_id

  !########################################################################
  ! return number of tile on current pe
  function mpp_get_current_ntile(domain)
     type(domain2d),                  intent(in) :: domain
     integer                                     :: mpp_get_current_ntile

     mpp_get_current_ntile = size(domain%tile_id(:))
     return

  end function mpp_get_current_ntile

  !#######################################################################
  ! return if current pe is the root pe of the tile, if number of tiles on current pe 
  ! is greater than 1, will return true, if isc==isg and jsc==jsg also will return true,
  ! otherwise false will be returned.
  function mpp_domain_is_root_pe(domain)
     type(domain2d), intent(in) :: domain
     logical                    :: mpp_domain_is_root_pe

     if(size(domain%tile_id(:)) > 1) then
        mpp_domain_is_root_pe = .true.
     else if(domain%x(1)%compute%begin == domain%x(1)%global%begin .AND.  & 
             domain%y(1)%compute%begin == domain%y(1)%global%begin ) then
        mpp_domain_is_root_pe = .true.
     else
        mpp_domain_is_root_pe = .false.
     end if

  end function
