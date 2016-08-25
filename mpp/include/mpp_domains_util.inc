! -*-f90-*-

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
    call mpp_malloc( ptr_domains_stack_nonblock, n, mpp_domains_stack_size )
#else
    if( allocated(mpp_domains_stack) )deallocate(mpp_domains_stack)
    allocate( mpp_domains_stack(n) )
    if( allocated(mpp_domains_stack_nonblock) )deallocate(mpp_domains_stack_nonblock)
    allocate( mpp_domains_stack_nonblock(n) )

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
    if( mpp_domain2D_eq )mpp_domain2D_eq = ALL(a%io_layout .EQ. b%io_layout)

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
  subroutine mpp_get_compute_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
       x_is_global, y_is_global, tile_count, position )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_count, position
    integer                        :: tile, ishift, jshift

    tile = 1
    if(present(tile_count)) tile = tile_count

    call mpp_get_compute_domain( domain%x(tile), xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_compute_domain( domain%y(tile), ybegin, yend, ysize, ymax_size, y_is_global )
    call mpp_get_domain_shift( domain, ishift, jshift, position )
    if( PRESENT(xend) ) xend  = xend + ishift
    if( PRESENT(yend) ) yend  = yend + jshift
    if( PRESENT(xsize)) xsize = xsize + ishift
    if( PRESENT(ysize)) ysize = ysize + jshift
    if(PRESENT(xmax_size))xmax_size = xmax_size + ishift
    if(PRESENT(ymax_size))ymax_size = ymax_size + jshift

    return
  end subroutine mpp_get_compute_domain2D

  !#####################################################################
  subroutine mpp_get_data_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
       x_is_global, y_is_global, tile_count, position )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_count, position
    integer                        :: tile, ishift, jshift

    tile = 1
    if(present(tile_count)) tile = tile_count

    call mpp_get_data_domain( domain%x(tile), xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_data_domain( domain%y(tile), ybegin, yend, ysize, ymax_size, y_is_global )
    call mpp_get_domain_shift( domain, ishift, jshift, position )
    if( PRESENT(xend) ) xend  = xend + ishift
    if( PRESENT(yend) ) yend  = yend + jshift
    if( PRESENT(xsize)) xsize = xsize + ishift
    if( PRESENT(ysize)) ysize = ysize + jshift
    if(PRESENT(xmax_size))xmax_size = xmax_size + ishift
    if(PRESENT(ymax_size))ymax_size = ymax_size + jshift

    return
  end subroutine mpp_get_data_domain2D

  !#####################################################################
  subroutine mpp_get_global_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
                                      tile_count, position )
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    integer, intent(in),  optional :: tile_count, position
    integer                        :: tile, ishift, jshift

    tile = 1
    if(present(tile_count)) tile = tile_count

    call mpp_get_global_domain( domain%x(tile), xbegin, xend, xsize, xmax_size )
    call mpp_get_global_domain( domain%y(tile), ybegin, yend, ysize, ymax_size )
    call mpp_get_domain_shift( domain, ishift, jshift, position )
    if( PRESENT(xend) ) xend  = xend + ishift
    if( PRESENT(yend) ) yend  = yend + jshift
    if( PRESENT(xsize)) xsize = xsize + ishift
    if( PRESENT(ysize)) ysize = ysize + jshift
    if(PRESENT(xmax_size))xmax_size = xmax_size + ishift
    if(PRESENT(ymax_size))ymax_size = ymax_size + jshift

    return
  end subroutine mpp_get_global_domain2D

  !#####################################################################
  subroutine mpp_get_memory_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
       x_is_global, y_is_global, position)
    type(domain2D),     intent(in) :: domain
    integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
    logical, intent(out), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: position
    integer                        :: tile, ishift, jshift

    tile = 1

    call mpp_get_memory_domain( domain%x(tile), xbegin, xend, xsize, xmax_size, x_is_global )
    call mpp_get_memory_domain( domain%y(tile), ybegin, yend, ysize, ymax_size, y_is_global )
    call mpp_get_domain_shift( domain, ishift, jshift, position )
    if( PRESENT(xend) ) xend  = xend + ishift
    if( PRESENT(yend) ) yend  = yend + jshift
    if( PRESENT(xsize)) xsize = xsize + ishift
    if( PRESENT(ysize)) ysize = ysize + jshift
    if(PRESENT(xmax_size))xmax_size = xmax_size + ishift
    if(PRESENT(ymax_size))ymax_size = ymax_size + jshift

    return
  end subroutine mpp_get_memory_domain2D

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
                                       x_is_global, y_is_global, tile_count )
    type(domain2D), intent(inout) :: domain
    integer, intent(in), optional :: xbegin, xend, ybegin, yend, xsize, ysize
    logical, intent(in), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_count
    integer                        :: tile

    tile = 1
    if(present(tile_count)) tile = tile_count

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
                                    x_is_global, y_is_global, tile_count )
    type(domain2D), intent(inout) :: domain
    integer, intent(in), optional :: xbegin, xend, ybegin, yend, xsize, ysize
    logical, intent(in), optional :: x_is_global, y_is_global
    integer, intent(in),  optional :: tile_count
    integer                        :: tile

    tile = 1
    if(present(tile_count)) tile = tile_count

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
  subroutine mpp_set_global_domain2D( domain, xbegin, xend, ybegin, yend, xsize, ysize, tile_count )
    type(domain2D), intent(inout) :: domain
    integer, intent(in), optional :: xbegin, xend, ybegin, yend, xsize, ysize
    integer, intent(in),  optional :: tile_count
    integer                        :: tile

    tile = 1
    if(present(tile_count)) tile = tile_count
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
  subroutine mpp_get_domain_components( domain, x, y, tile_count )
    type(domain2D),            intent(in) :: domain
    type(domain1D), intent(inout), optional :: x, y
    integer, intent(in),  optional :: tile_count
    integer                        :: tile

    tile = 1
    if(present(tile_count)) tile = tile_count
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
subroutine mpp_get_compute_domains2D( domain, xbegin, xend, xsize, ybegin, yend, ysize, position )
 type(domain2D),                   intent(in) :: domain
 integer, intent(out), optional, dimension(:) :: xbegin, xend, xsize, ybegin, yend, ysize
 integer, intent(in ), optional               :: position

 integer :: i, ishift, jshift

 call mpp_get_domain_shift( domain, ishift, jshift, position )


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
       xend(i) = domain%list(i-1)%x(1)%compute%end + ishift
    end do
 end if
 if( PRESENT(xsize) )then
    if( size(xsize(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xsize array size does not match domain.' )
    do i = 1, size(xsize(:))
       xsize(i) = domain%list(i-1)%x(1)%compute%size + ishift
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
       yend(i) = domain%list(i-1)%y(1)%compute%end + jshift
    end do
 end if
 if( PRESENT(ysize) )then
    if( size(ysize(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: ysize array size does not match domain.' )
    do i = 1, size(ysize(:))
       ysize(i) = domain%list(i-1)%y(1)%compute%size + jshift
    end do
 end if
 return
end subroutine mpp_get_compute_domains2D

!#####################################################################
subroutine mpp_get_domain_extents1D(domain, xextent, yextent)
   type(domain2d),            intent(in) :: domain
   integer, dimension(0:), intent(inout) :: xextent, yextent
   integer                               :: n

   if(domain%ntiles .NE. 1) call mpp_error(FATAL,"mpp_domains_util.inc(mpp_get_domain_extents1D): "// &
       "ntiles is more than 1, please use mpp_get_domain_extents2D")
   if(size(xextent) .NE. size(domain%x(1)%list(:)))  call mpp_error(FATAL,"mpp_domains_util.inc(mpp_get_domain_extents1D): "// &
       "size(xextent) does not equal to size(domain%x(1)%list(:)))")
   if(size(yextent) .NE. size(domain%y(1)%list(:)))  call mpp_error(FATAL,"mpp_domains_util.inc(mpp_get_domain_extents1D): "// &
       "size(yextent) does not equal to size(domain%y(1)%list(:)))")
   do n = 0, size(domain%x(1)%list(:))-1
      xextent(n) = domain%x(1)%list(n)%compute%size
   enddo
   do n = 0, size(domain%y(1)%list(:))-1
      yextent(n) = domain%y(1)%list(n)%compute%size
   enddo

end subroutine mpp_get_domain_extents1D

!#####################################################################
! This will return xextent and yextent for each tile
subroutine mpp_get_domain_extents2D(domain, xextent, yextent)
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


end subroutine mpp_get_domain_extents2D
   
!#####################################################################
function mpp_get_domain_pe(domain)
   type(domain2d), intent(in) :: domain
   integer                    :: mpp_get_domain_pe

   mpp_get_domain_pe = domain%pe


end function mpp_get_domain_pe


function mpp_get_domain_tile_root_pe(domain)
   type(domain2d), intent(in) :: domain
   integer                    :: mpp_get_domain_tile_root_pe

   mpp_get_domain_tile_root_pe = domain%tile_root_pe

end function mpp_get_domain_tile_root_pe

function mpp_get_io_domain(domain)
   type(domain2d), intent(in) :: domain
   type(domain2d), pointer    :: mpp_get_io_domain

   if(ASSOCIATED(domain%io_domain)) then
      mpp_get_io_domain => domain%io_domain
   else
      mpp_get_io_domain => NULL()
   endif

end function mpp_get_io_domain

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
  function mpp_domain_is_initialized(domain)
    type(domain2d), intent(in) :: domain
    logical                    :: mpp_domain_is_initialized

    mpp_domain_is_initialized = domain%initialized

    return

  end function mpp_domain_is_initialized

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
       if( debug )call mpp_error(NOTE, &
          'mpp_domains_util.inc: halo size to be updated are all zero, no update will be done')
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
  function search_update_overlap(domain, whalo, ehalo, shalo, nhalo, position)
    type(domain2d),          intent(inout) :: domain
    integer,                    intent(in) :: whalo, ehalo, shalo, nhalo
    integer,                    intent(in) :: position
    type(overlapSpec),          pointer    :: search_update_overlap
    type(overlapSpec),          pointer    :: update_ref

    select case(position)
    case (CENTER)
       update_ref => domain%update_T
    case (CORNER)
       update_ref => domain%update_C
    case (NORTH)
       update_ref => domain%update_N
    case (EAST)
       update_ref => domain%update_E
    case default
       call mpp_error(FATAL,"mpp_domains_util.inc(search_update_overlap): position should be CENTER|CORNER|EAST|NORTH")
    end select

    search_update_overlap => update_ref

    do 
       if(whalo == search_update_overlap%whalo .AND. ehalo == search_update_overlap%ehalo .AND. &
            shalo == search_update_overlap%shalo .AND. nhalo == search_update_overlap%nhalo ) then
            exit ! found domain
       endif
       !--- if not found, switch to next
       if(.NOT. ASSOCIATED(search_update_overlap%next)) then
          allocate(search_update_overlap%next)
          search_update_overlap => search_update_overlap%next
          call set_overlaps(domain, update_ref, search_update_overlap, whalo, ehalo, shalo, nhalo )
          exit
       else
          search_update_overlap => search_update_overlap%next
       end if
  
    end do

    update_ref => NULL()

  end function search_update_overlap

!#######################################################################
  ! this routine found the check at certain position 
  function search_check_overlap(domain, position)
    type(domain2d),             intent(in) :: domain
    integer,                    intent(in) :: position
    type(overlapSpec),          pointer    :: search_check_overlap

    select case(position)
    case (CENTER)
       search_check_overlap => NULL()
    case (CORNER)
       search_check_overlap => domain%check_C
    case (NORTH)
       search_check_overlap => domain%check_N
    case (EAST)
       search_check_overlap => domain%check_E
    case default
       call mpp_error(FATAL,"mpp_domains_util.inc(search_check_overlap): position should be CENTER|CORNER|EAST|NORTH")
    end select

  end function search_check_overlap

!#######################################################################
  ! this routine found the bound at certain position 
  function search_bound_overlap(domain, position)
    type(domain2d),             intent(in) :: domain
    integer,                    intent(in) :: position
    type(overlapSpec),          pointer    :: search_bound_overlap

    select case(position)
    case (CENTER)
       search_bound_overlap => NULL()
    case (CORNER)
       search_bound_overlap => domain%bound_C
    case (NORTH)
       search_bound_overlap => domain%bound_N
    case (EAST)
       search_bound_overlap => domain%bound_E
    case default
       call mpp_error(FATAL,"mpp_domains_util.inc(search_bound_overlap): position should be CENTER|CORNER|EAST|NORTH")
    end select

  end function search_bound_overlap

  !########################################################################
  ! return the tile_id on current pe
  function mpp_get_tile_id(domain)
     type(domain2d),                  intent(in) :: domain
     integer, dimension(size(domain%tile_id(:))) :: mpp_get_tile_id

     mpp_get_tile_id = domain%tile_id
     return

  end function mpp_get_tile_id

  !#######################################################################
  ! return the tile_id on current pelist. one-tile-per-pe is assumed.
  subroutine mpp_get_tile_list(domain, tiles)
     type(domain2d), intent(in) :: domain
     integer,     intent(inout) :: tiles(:)
     integer                    :: i

    if( size(tiles(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'mpp_get_tile_list: tiles array size does not match domain.' )
    do i = 1, size(tiles(:))
       if(size(domain%list(i-1)%tile_id(:)) > 1) call mpp_error( FATAL,  &
               'mpp_get_tile_list: only support one-tile-per-pe now, contact developer');
       tiles(i) = domain%list(i-1)%tile_id(1)
    end do     

  end subroutine mpp_get_tile_list

  !########################################################################
  ! return number of tiles in mosaic
  function mpp_get_ntile_count(domain)
     type(domain2d),                  intent(in) :: domain
     integer                                     :: mpp_get_ntile_count

     mpp_get_ntile_count = domain%ntiles
     return

  end function mpp_get_ntile_count

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
  function mpp_domain_is_tile_root_pe(domain)
     type(domain2d), intent(in) :: domain
     logical                    :: mpp_domain_is_tile_root_pe

     mpp_domain_is_tile_root_pe = domain%pe == domain%tile_root_pe;

  end function mpp_domain_is_tile_root_pe

  !#########################################################################
  ! return number of processors used on current tile.
  function mpp_get_tile_npes(domain)
    type(domain2d), intent(in)  :: domain
    integer                     :: mpp_get_tile_npes
    integer                     :: i, tile

    !--- When there is more than one tile on this pe, we assume each tile will be 
    !--- limited to this pe.
    if(size(domain%tile_id(:)) > 1) then
       mpp_get_tile_npes = 1
    else
      mpp_get_tile_npes = 0
      tile = domain%tile_id(1)
      do i = 0, size(domain%list(:))-1
         if(tile == domain%list(i)%tile_id(1) )  mpp_get_tile_npes =  mpp_get_tile_npes + 1
      end do     
    endif

  end function mpp_get_tile_npes

  !########################################################################
  ! get the processors list used on current tile.
  subroutine mpp_get_tile_pelist(domain, pelist)
    type(domain2d), intent(in)  :: domain
    integer,      intent(inout) :: pelist(:)
    integer                     :: npes_on_tile
    integer                     :: i, tile, pos

    npes_on_tile = mpp_get_tile_npes(domain)
    if(size(pelist(:)) .NE. npes_on_tile) call mpp_error(FATAL, &
           "mpp_domains_util.inc(mpp_get_tile_pelist): size(pelist) does not equal npes on current tile")
    tile = domain%tile_id(1)
    pos = 0
    do i = 0, size(domain%list(:))-1
       if(tile == domain%list(i)%tile_id(1)) then
          pos = pos+1
          pelist(pos) = domain%list(i)%pe
       endif
    enddo

    return

  end subroutine mpp_get_tile_pelist

!#####################################################################
subroutine mpp_get_tile_compute_domains( domain, xbegin, xend, ybegin, yend, position )
 type(domain2D),         intent(in) :: domain
 integer, intent(out), dimension(:) :: xbegin, xend, ybegin, yend
 integer, intent(in ), optional     :: position

 integer :: i, ishift, jshift
 integer :: npes_on_tile, pos, tile

 call mpp_get_domain_shift( domain, ishift, jshift, position )


 if( .NOT.module_is_initialized ) &
      call mpp_error( FATAL, 'mpp_get_compute_domains2D: must first call mpp_domains_init.' )

 npes_on_tile = mpp_get_tile_npes(domain)
 if(size(xbegin(:)) .NE. npes_on_tile) call mpp_error(FATAL, &
       "mpp_domains_util.inc(mpp_get_compute_domains2D): size(xbegin) does not equal npes on current tile")
 if(size(xend(:)) .NE. npes_on_tile) call mpp_error(FATAL, &
       "mpp_domains_util.inc(mpp_get_compute_domains2D): size(xend) does not equal npes on current tile")
 if(size(ybegin(:)) .NE. npes_on_tile) call mpp_error(FATAL, &
       "mpp_domains_util.inc(mpp_get_compute_domains2D): size(ybegin) does not equal npes on current tile")
 if(size(yend(:)) .NE. npes_on_tile) call mpp_error(FATAL, &
       "mpp_domains_util.inc(mpp_get_compute_domains2D): size(yend) does not equal npes on current tile")

 tile = domain%tile_id(1)
 pos = 0
 do i = 0, size(domain%list(:))-1
    if(tile == domain%list(i)%tile_id(1)) then
       pos = pos+1
       xbegin(pos) = domain%list(i)%x(1)%compute%begin
       xend  (pos) = domain%list(i)%x(1)%compute%end + ishift
       ybegin(pos) = domain%list(i)%y(1)%compute%begin
       yend  (pos) = domain%list(i)%y(1)%compute%end + jshift
    endif
 enddo

 return

end subroutine mpp_get_tile_compute_domains



  !#############################################################################
  function mpp_get_num_overlap(domain, action, p, position)
    type(domain2d),    intent(in) :: domain
    integer,           intent(in) :: action
    integer,           intent(in) :: p
    integer, optional, intent(in) :: position
    integer                       :: mpp_get_num_overlap
    type(overlapSpec), pointer    :: update => NULL()
    integer                       :: pos

    pos = CENTER
    if(present(position)) pos = position
    select case(pos)
    case (CENTER)
       update => domain%update_T
    case (CORNER)
       update => domain%update_C
    case (EAST)
       update => domain%update_E
    case (NORTH)
       update => domain%update_N
    case default
       call mpp_error( FATAL, "mpp_domains_mod(mpp_get_num_overlap): invalid option of position")
    end select

    if(action == EVENT_SEND) then
       if(p< 1 .OR. p > update%nsend) call mpp_error( FATAL, &
                "mpp_domains_mod(mpp_get_num_overlap): p should be between 1 and update%nsend")
       mpp_get_num_overlap = update%send(p)%count
    else if(action == EVENT_RECV) then
       if(p< 1 .OR. p > update%nrecv) call mpp_error( FATAL, &
                "mpp_domains_mod(mpp_get_num_overlap): p should be between 1 and update%nrecv")
       mpp_get_num_overlap = update%recv(p)%count
    else
       call mpp_error( FATAL, "mpp_domains_mod(mpp_get_num_overlap): invalid option of action")
    end if

  end function mpp_get_num_overlap

  !#############################################################################
  subroutine mpp_get_update_size(domain, nsend, nrecv, position)
    type(domain2d),      intent(in) :: domain
    integer,            intent(out) :: nsend, nrecv
    integer, optional,   intent(in) :: position
    integer                         :: pos

    pos = CENTER
    if(present(position)) pos = position
    select case(pos)
    case (CENTER)
       nsend = domain%update_T%nsend
       nrecv = domain%update_T%nrecv
    case (CORNER)
       nsend = domain%update_C%nsend
       nrecv = domain%update_C%nrecv
    case (EAST)
       nsend = domain%update_E%nsend
       nrecv = domain%update_E%nrecv
    case (NORTH)
       nsend = domain%update_N%nsend
       nrecv = domain%update_N%nrecv
    case default
       call mpp_error( FATAL, "mpp_domains_mod(mpp_get_update_size): invalid option of position")
    end select

  end subroutine mpp_get_update_size

  !#############################################################################
  subroutine mpp_get_update_pelist(domain, action, pelist, position)
    type(domain2d),       intent(in) :: domain
    integer,              intent(in) :: action
    integer,           intent(inout) :: pelist(:)
    integer, optional,    intent(in) :: position
    type(overlapSpec), pointer    :: update => NULL()
    integer                       :: pos, p

    pos = CENTER
    if(present(position)) pos = position
    select case(pos)
    case (CENTER)
       update => domain%update_T
    case (CORNER)
       update => domain%update_C
    case (EAST)
       update => domain%update_E
    case (NORTH)
       update => domain%update_N
    case default
       call mpp_error( FATAL, "mpp_domains_mod(mpp_get_update_pelist): invalid option of position")
    end select

    if(action == EVENT_SEND) then
       if(size(pelist) .NE. update%nsend) call mpp_error( FATAL,  &
             "mpp_domains_mod(mpp_get_update_pelist): size of pelist does not match update%nsend")
       do p = 1, update%nsend
          pelist(p) = update%send(p)%pe
       enddo
    else if(action == EVENT_RECV) then
       if(size(pelist) .NE. update%nrecv) call mpp_error( FATAL,  &
             "mpp_domains_mod(mpp_get_update_pelist): size of pelist does not match update%nrecv")
       do p = 1, update%nrecv
          pelist(p) = update%recv(p)%pe
       enddo
    else
       call mpp_error( FATAL, "mpp_domains_mod(mpp_get_update_pelist): invalid option of action")
    end if

  end subroutine mpp_get_update_pelist

  !#############################################################################
  subroutine mpp_get_overlap(domain, action, p, is, ie, js, je, dir, rot, position)
    type(domain2d),         intent(in) :: domain
    integer,                intent(in) :: action
    integer,                intent(in) :: p
    integer, dimension(:), intent(out) :: is, ie, js, je
    integer, dimension(:), intent(out) :: dir, rot    
    integer, optional,      intent(in) :: position
    type(overlapSpec),      pointer    :: update => NULL()
    type(overlap_type),     pointer    :: overlap => NULL()
    integer                            :: count, pos

    pos = CENTER
    if(present(position)) pos = position
    select case(pos)
    case (CENTER)
       update => domain%update_T
    case (CORNER)
       update => domain%update_C
    case (EAST)
       update => domain%update_E
    case (NORTH)
       update => domain%update_N
    case default
       call mpp_error( FATAL, "mpp_domains_mod(mpp_get_overlap): invalid option of position")
    end select

    if(action == EVENT_SEND) then
       overlap => update%send(p)
    else if(action == EVENT_RECV) then
       overlap => update%recv(p)
    else
       call mpp_error( FATAL, "mpp_domains_mod(mpp_get_overlap): invalid option of action")
    end if

    count = overlap%count
    if(size(is(:)) .NE. count .OR. size(ie (:)) .NE. count .OR. size(js (:)) .NE. count .OR.  &
       size(je(:)) .NE. count .OR. size(dir(:)) .NE. count .OR. size(rot(:)) .NE. count ) &
       call mpp_error( FATAL, "mpp_domains_mod(mpp_get_overlap): size mismatch between number of overlap and array size")

    is  = overlap%is      (1:count)
    ie  = overlap%ie      (1:count)
    js  = overlap%js      (1:count)
    je  = overlap%je      (1:count)
    dir = overlap%dir     (1:count)
    rot = overlap%rotation(1:count)

    update  => NULL()
    overlap => NULL()

  end subroutine mpp_get_overlap

  !##################################################################
  function mpp_get_domain_name(domain)
     type(domain2d),    intent(in) :: domain
     character(len=NAME_LENGTH)    :: mpp_get_domain_name

     mpp_get_domain_name = domain%name

  end function mpp_get_domain_name

  !#################################################################
  function mpp_get_domain_root_pe(domain)
     type(domain2d),    intent(in) :: domain
     integer                       :: mpp_get_domain_root_pe
     
     mpp_get_domain_root_pe = domain%list(0)%pe

  end function mpp_get_domain_root_pe
  !#################################################################
  function mpp_get_domain_npes(domain)
     type(domain2d),    intent(in) :: domain
     integer                       :: mpp_get_domain_npes  

     mpp_get_domain_npes = size(domain%list(:))

     return

  end function mpp_get_domain_npes

  !################################################################
  subroutine mpp_get_domain_pelist(domain, pelist)
     type(domain2d),    intent(in) :: domain
     integer,          intent(out) :: pelist(:)
     integer :: p

     if(size(pelist(:)) .NE. size(domain%list(:)) ) then
        call mpp_error(FATAL, "mpp_get_domain_pelist: size(pelist(:)) .NE. size(domain%list(:)) ")
     endif

     do p = 0, size(domain%list(:))-1     
        pelist(p+1) = domain%list(p)%pe 
     enddo

     return

  end subroutine mpp_get_domain_pelist

  !#################################################################
  function mpp_get_io_domain_layout(domain)
     type(domain2d),    intent(in) :: domain
     integer, dimension(2)         :: mpp_get_io_domain_layout

     mpp_get_io_domain_layout = domain%io_layout

  end function mpp_get_io_domain_layout

  !################################################################
  function get_rank_send(domain, overlap_x, overlap_y, rank_x, rank_y, ind_x, ind_y)
    type(domain2D),    intent(in) :: domain
    type(overlapSpec), intent(in) :: overlap_x, overlap_y
    integer,          intent(out) :: rank_x, rank_y, ind_x, ind_y
    integer                       :: get_rank_send
    integer                       :: nlist, nsend_x, nsend_y

    nlist = size(domain%list(:))
    nsend_x = overlap_x%nsend
    nsend_y = overlap_y%nsend
    rank_x = nlist+1
    rank_y = nlist+1
    if(nsend_x>0) rank_x = overlap_x%send(1)%pe - domain%pe
    if(nsend_y>0) rank_y = overlap_y%send(1)%pe - domain%pe
    if(rank_x .LT. 0) rank_x = rank_x + nlist
    if(rank_y .LT. 0) rank_y = rank_y + nlist
    get_rank_send = min(rank_x, rank_y)
    ind_x = nsend_x + 1
    ind_y = nsend_y + 1
    if(get_rank_send < nlist+1) then
       if(nsend_x>0) ind_x = 1
       if(nsend_y>0) ind_y = 1
    endif

  end function get_rank_send

  !############################################################################
  function get_rank_recv(domain, overlap_x, overlap_y, rank_x, rank_y, ind_x, ind_y)
    type(domain2D),    intent(in) :: domain
    type(overlapSpec), intent(in) :: overlap_x, overlap_y
    integer,          intent(out) :: rank_x, rank_y, ind_x, ind_y
    integer                       :: get_rank_recv
    integer                       :: nlist, nrecv_x, nrecv_y

    nlist = size(domain%list(:))
    nrecv_x = overlap_x%nrecv
    nrecv_y = overlap_y%nrecv
    rank_x = -1
    rank_y = -1
    if(nrecv_x>0) then
       rank_x = overlap_x%recv(1)%pe - domain%pe
       if(rank_x .LE. 0) rank_x = rank_x + nlist
    endif
    if(nrecv_y>0) then
       rank_y = overlap_y%recv(1)%pe - domain%pe    
       if(rank_y .LE. 0) rank_y = rank_y + nlist
    endif
    get_rank_recv = max(rank_x, rank_y)
    ind_x = nrecv_x + 1
    ind_y = nrecv_y + 1
    if(get_rank_recv < nlist+1) then
       if(nrecv_x>0) ind_x = 1
       if(nrecv_y>0) ind_y = 1
    endif

  end function get_rank_recv

  function get_vector_recv(domain, update_x, update_y, ind_x, ind_y, start_pos, pelist)
    type(domain2D),    intent(in) :: domain
    type(overlapSpec), intent(in) :: update_x, update_y
    integer,          intent(out) :: ind_x(:), ind_y(:)
    integer,          intent(out) :: start_pos(:)
    integer,          intent(out) :: pelist(:)
    integer                       :: nlist, nrecv_x, nrecv_y, ntot, n
    integer                       :: ix, iy, rank_x, rank_y, cur_pos
    integer                       :: get_vector_recv

    nlist = size(domain%list(:))
    nrecv_x = update_x%nrecv
    nrecv_y = update_y%nrecv
 
    ntot = nrecv_x + nrecv_y

    n  = 1
    ix = 1
    iy = 1
    ind_x = -1
    ind_y = -1
    get_vector_recv = 0
    cur_pos = 0
    do while (n<=ntot)
       if(ix <= nrecv_x ) then
          rank_x = update_x%recv(ix)%pe-domain%pe
          if(rank_x .LE. 0) rank_x = rank_x + nlist
       else
          rank_x = -1
       endif
       if(iy <= nrecv_y ) then
          rank_y = update_y%recv(iy)%pe-domain%pe
          if(rank_y .LE. 0) rank_y = rank_y + nlist
       else
          rank_y = -1
       endif
       get_vector_recv = get_vector_recv + 1
       start_pos(get_vector_recv) = cur_pos
       if( rank_x == rank_y ) then
          n = n+2
          ind_x (get_vector_recv) = ix
          ind_y (get_vector_recv) = iy
          cur_pos = cur_pos + update_x%recv(ix)%totsize + update_y%recv(iy)%totsize
          pelist(get_vector_recv) = update_x%recv(ix)%pe
          ix = ix + 1
          iy = iy + 1
       else if ( rank_x > rank_y ) then
          n = n+1
          ind_x (get_vector_recv) = ix
          ind_y (get_vector_recv) = -1
          cur_pos = cur_pos + update_x%recv(ix)%totsize
          pelist(get_vector_recv) = update_x%recv(ix)%pe
          ix = ix + 1
       else if ( rank_y > rank_x ) then
          n = n+1
          ind_x (get_vector_recv) = -1
          ind_y (get_vector_recv) = iy
          cur_pos = cur_pos + update_y%recv(iy)%totsize
          pelist(get_vector_recv) = update_y%recv(iy)%pe
          iy = iy+1
       endif
    end do


  end function get_vector_recv

  function get_vector_send(domain, update_x, update_y, ind_x, ind_y, start_pos, pelist)
    type(domain2D),    intent(in) :: domain
    type(overlapSpec), intent(in) :: update_x, update_y
    integer,          intent(out) :: ind_x(:), ind_y(:)
    integer,          intent(out) :: start_pos(:)
    integer,          intent(out) :: pelist(:)
    integer                       :: nlist, nsend_x, nsend_y, ntot, n
    integer                       :: ix, iy, rank_x, rank_y, cur_pos
    integer                       :: get_vector_send

    nlist = size(domain%list(:))
    nsend_x = update_x%nsend
    nsend_y = update_y%nsend
 
    ntot = nsend_x + nsend_y
    n  = 1
    ix = 1
    iy = 1
    ind_x = -1
    ind_y = -1
    get_vector_send = 0
    cur_pos = 0
    do while (n<=ntot)
       if(ix <= nsend_x ) then
          rank_x = update_x%send(ix)%pe-domain%pe
          if(rank_x .LT. 0) rank_x = rank_x + nlist
       else
          rank_x = nlist+1
       endif
       if(iy <= nsend_y ) then
          rank_y = update_y%send(iy)%pe-domain%pe
          if(rank_y .LT. 0) rank_y = rank_y + nlist
       else
          rank_y = nlist+1
       endif
       get_vector_send = get_vector_send + 1
       start_pos(get_vector_send) = cur_pos

       if( rank_x == rank_y ) then
          n = n+2
          ind_x  (get_vector_send) = ix
          ind_y  (get_vector_send) = iy
          cur_pos = cur_pos + update_x%send(ix)%totsize + update_y%send(iy)%totsize
          pelist (get_vector_send) = update_x%send(ix)%pe
          ix = ix + 1
          iy = iy + 1
       else if ( rank_x < rank_y ) then
          n = n+1
          ind_x  (get_vector_send) = ix
          ind_y  (get_vector_send) = -1
          cur_pos = cur_pos + update_x%send(ix)%totsize
          pelist (get_vector_send) = update_x%send(ix)%pe
          ix = ix + 1
       else if ( rank_y < rank_x ) then
          n = n+1
          ind_x  (get_vector_send) = -1
          ind_y  (get_vector_send) = iy
          cur_pos = cur_pos + update_y%send(iy)%totsize
          pelist (get_vector_send) = update_y%send(iy)%pe
          iy = iy+1
       endif
    end do


  end function get_vector_send


  !############################################################################  
  function get_rank_unpack(domain, overlap_x, overlap_y, rank_x, rank_y, ind_x, ind_y)
    type(domain2D),    intent(in) :: domain
    type(overlapSpec), intent(in) :: overlap_x, overlap_y
    integer,          intent(out) :: rank_x, rank_y, ind_x, ind_y
    integer                       :: get_rank_unpack
    integer                       :: nlist, nrecv_x, nrecv_y

    nlist = size(domain%list(:))
    nrecv_x = overlap_x%nrecv
    nrecv_y = overlap_y%nrecv

    rank_x = nlist+1
    rank_y = nlist+1
    if(nrecv_x>0) rank_x = overlap_x%recv(nrecv_x)%pe - domain%pe
    if(nrecv_y>0) rank_y = overlap_y%recv(nrecv_y)%pe - domain%pe
    if(rank_x .LE.0) rank_x = rank_x + nlist
    if(rank_y .LE.0) rank_y = rank_y + nlist

    get_rank_unpack = min(rank_x, rank_y)
    ind_x = 0
    ind_y = 0
    if(get_rank_unpack < nlist+1) then
       if(nrecv_x >0) ind_x = nrecv_x
       if(nrecv_y >0) ind_y = nrecv_y
    endif

  end function get_rank_unpack

  function get_mesgsize(overlap, do_dir)
   type(overlap_type), intent(in) :: overlap
   logical,            intent(in) :: do_dir(:)
   integer                       :: get_mesgsize
   integer                       :: n, dir

   get_mesgsize = 0
   do n = 1, overlap%count
      dir = overlap%dir(n)
      if(do_dir(dir)) then
         get_mesgsize = get_mesgsize + overlap%msgsize(n)
      end if
   end do

  end function get_mesgsize

  !#############################################################################
  subroutine mpp_set_domain_symmetry(domain, symmetry)
    type(domain2D), intent(inout) :: domain
    logical,        intent(in   ) :: symmetry

    domain%symmetry = symmetry

  end subroutine mpp_set_domain_symmetry
  

  subroutine mpp_copy_domain1D(domain_in, domain_out)
     type(domain1D),    intent(in) :: domain_in
     type(domain1D), intent(inout) :: domain_out

     domain_out%compute = domain_in%compute
     domain_out%data    = domain_in%data
     domain_out%global  = domain_in%global
     domain_out%memory  = domain_in%memory
     domain_out%cyclic  = domain_in%cyclic
     domain_out%pe      = domain_in%pe 
     domain_out%pos     = domain_in%pos 

  end subroutine mpp_copy_domain1D

  !#################################################################
  !z1l: This is not fully implemented. The current purpose is to make
  !     it work in read_data.
  subroutine mpp_copy_domain2D(domain_in, domain_out)
     type(domain2D),    intent(in) :: domain_in
     type(domain2D), intent(inout) :: domain_out

     integer :: n, ntiles

     domain_out%id             = domain_in%id
     domain_out%pe             = domain_in%pe
     domain_out%fold           = domain_in%fold
     domain_out%pos            = domain_in%pos
     domain_out%symmetry       = domain_in%symmetry
     domain_out%whalo          = domain_in%whalo
     domain_out%ehalo          = domain_in%ehalo
     domain_out%shalo          = domain_in%shalo
     domain_out%nhalo          = domain_in%nhalo
     domain_out%ntiles         = domain_in%ntiles
     domain_out%max_ntile_pe   = domain_in%max_ntile_pe
     domain_out%ncontacts      = domain_in%ncontacts
     domain_out%rotated_ninety = domain_in%rotated_ninety
     domain_out%initialized    = domain_in%initialized
     domain_out%tile_root_pe   = domain_in%tile_root_pe
     domain_out%io_layout      = domain_in%io_layout
     domain_out%name           = domain_in%name

     ntiles = size(domain_in%x(:))
     allocate(domain_out%x(ntiles), domain_out%y(ntiles), domain_out%tile_id(ntiles) )
     do n = 1, ntiles
        call mpp_copy_domain1D(domain_in%x(n), domain_out%x(n))
        call mpp_copy_domain1D(domain_in%y(n), domain_out%y(n))
     enddo
     domain_out%tile_id = domain_in%tile_id

     return

  end subroutine mpp_copy_domain2D

  !######################################################################
  subroutine set_group_update(group, domain)
     type(mpp_group_update_type), intent(inout) :: group  
     type(domain2D),              intent(inout) :: domain
     integer   :: nscalar, nvector, nlist
     integer   :: nsend, nrecv, nsend_old, nrecv_old
     integer   :: nsend_s, nsend_x, nsend_y
     integer   :: nrecv_s, nrecv_x, nrecv_y
     integer   :: update_buffer_pos, tot_recv_size, tot_send_size
     integer   :: msgsize_s, msgsize_x, msgsize_y, msgsize
     logical   :: recv_s(8), send_s(8), recv_v(8), send_v(8)
     integer   :: ntot, n, l, m, ksize
     integer   :: i_s, i_x, i_y, rank_s, rank_x, rank_y, rank
     integer   :: ind_s(3*MAXOVERLAP)
     integer   :: ind_x(3*MAXOVERLAP)
     integer   :: ind_y(3*MAXOVERLAP)
     integer   :: pelist(3*MAXOVERLAP), to_pe_list(3*MAXOVERLAP)
     integer   :: buffer_pos_recv(3*MAXOVERLAP), buffer_pos_send(3*MAXOVERLAP)
     integer   :: recv_size(3*MAXOVERLAP), send_size(3*MAXOVERLAP)
     integer   :: position_x, position_y, npack, nunpack, dir
     integer   :: pack_buffer_pos, unpack_buffer_pos
     integer   :: omp_get_num_threads, nthreads
     character(len=8)            :: text
     type(overlap_type), pointer :: overPtr => NULL()
     type(overlapSpec), pointer  :: update_s => NULL()
     type(overlapSpec), pointer  :: update_x => NULL()
     type(overlapSpec), pointer  :: update_y => NULL()

     nscalar = group%nscalar
     nvector = group%nvector

     !--- get the overlap data type
     select case(group%gridtype)
     case (AGRID)
        position_x = CENTER
        position_y = CENTER
     case (BGRID_NE, BGRID_SW)
        position_x = CORNER
        position_y = CORNER
     case (CGRID_NE, CGRID_SW)
        position_x = EAST
        position_y = NORTH
     case (DGRID_NE, DGRID_SW)
        position_x = NORTH
        position_y = EAST
     case default
        call mpp_error(FATAL, "set_group_update: invalid value of gridtype")
     end select
     if(nscalar>0) then
        update_s => search_update_overlap(domain, group%whalo_s, group%ehalo_s, &
             group%shalo_s, group%nhalo_s, group%position)
     endif
     if(nvector>0) then
        update_x => search_update_overlap(domain, group%whalo_v, group%ehalo_v, &
             group%shalo_v, group%nhalo_v, position_x)
        update_y => search_update_overlap(domain, group%whalo_v, group%ehalo_v, &
             group%shalo_v, group%nhalo_v, position_y)
     endif

     if(nscalar > 0) then
        recv_s = group%recv_s
        send_s = recv_s
     endif
     if(nvector > 0) then
        recv_v = group%recv_v
        send_v = recv_v
     end if
     nlist   = size(domain%list(:))
     group%initialized = .true.
     nsend_s = 0; nsend_x = 0; nsend_y = 0
     nrecv_s = 0; nrecv_x = 0; nrecv_y = 0

     if(nscalar > 0) then
         !--- This check could not be done because of memory domain
!        if( group%isize_s .NE. (group%ie_s-group%is_s+1) .OR. group%jsize_s .NE. (group%je_s-group%js_s+1))  &
!           call mpp_error(FATAL, "set_group_update: mismatch of size of the field and domain memory domain")
        nsend_s = update_s%nsend
        nrecv_s = update_s%nrecv
     endif

     !--- ksize_s must equal ksize_v
     if(nvector > 0 .AND. nscalar > 0) then
        if(group%ksize_s .NE. group%ksize_v) then
           call mpp_error(FATAL, "set_group_update: ksize_s and ksize_v are not equal")
        endif
        ksize = group%ksize_s
     else if (nscalar > 0) then
        ksize = group%ksize_s
     else if (nvector > 0) then
        ksize = group%ksize_v
     else
        call mpp_error(FATAL, "set_group_update: nscalar and nvector are all 0")
     endif

    nthreads = 1
!$OMP PARALLEL
!$   nthreads = omp_get_num_threads()
!$OMP END PARALLEL
     if( nthreads > nthread_control_loop ) then
        group%k_loop_inside = .FALSE.
     else
        group%k_loop_inside = .TRUE.
     endif

     if(nvector > 0) then
        !--- This check could not be done because of memory domain
!        if( group%isize_x .NE. (group%ie_x-group%is_x+1) .OR. group%jsize_x .NE. (group%je_x-group%js_x+1))  &
!           call mpp_error(FATAL, "set_group_update: mismatch of size of the fieldx and domain memory domain")
!        if( group%isize_y .NE. (group%ie_y-group%is_y+1) .OR. group%jsize_y .NE. (group%je_y-group%js_y+1))  &
!           call mpp_error(FATAL, "set_group_update: mismatch of size of the fieldy and domain memory domain")
        nsend_x = update_x%nsend
        nrecv_x = update_x%nrecv
        nsend_y = update_y%nsend
        nrecv_y = update_y%nrecv
     endif

     !figure out message size for each processor.
     ntot = nrecv_s + nrecv_x + nrecv_y
     if(ntot > 3*MAXOVERLAP) call mpp_error(FATAL, "set_group_update: ntot is greater than 3*MAXOVERLAP")
     n = 1
     i_s = 1
     i_x = 1
     i_y = 1
     ind_s = -1
     ind_x = -1
     ind_y = -1
     nrecv = 0
     do while(n<=ntot)
        if( i_s <= nrecv_s ) then
           rank_s = update_s%recv(i_s)%pe-domain%pe
           if(rank_s .LE. 0) rank_s = rank_s + nlist
        else
           rank_s = -1
        endif
        if( i_x <= nrecv_x ) then
           rank_x = update_x%recv(i_x)%pe-domain%pe
           if(rank_x .LE. 0) rank_x = rank_x + nlist
        else
           rank_x = -1
        endif
        if( i_y <= nrecv_y ) then
           rank_y = update_y%recv(i_y)%pe-domain%pe
           if(rank_y .LE. 0) rank_y = rank_y + nlist
        else
           rank_y = -1
        endif
        nrecv = nrecv + 1   
        rank = maxval((/rank_s, rank_x, rank_y/))
        if(rank == rank_s) then
           n = n + 1
           ind_s(nrecv) = i_s
           pelist(nrecv) = update_s%recv(i_s)%pe
           i_s = i_s + 1
        endif
        if(rank == rank_x) then
           n = n + 1
           ind_x(nrecv) = i_x
           pelist(nrecv) = update_x%recv(i_x)%pe
           i_x = i_x + 1
        endif
        if(rank == rank_y) then
           n = n + 1
           ind_y(nrecv) = i_y
           pelist(nrecv) = update_y%recv(i_y)%pe
           i_y = i_y + 1
        endif
     enddo

     nrecv_old = nrecv
     nrecv     = 0
     update_buffer_pos = 0
     tot_recv_size = 0

     !--- setup for recv
     do l = 1, nrecv_old
        msgsize_s = 0
        msgsize_x = 0
        msgsize_y = 0
        m = ind_s(l)
        if(m>0) msgsize_s = get_mesgsize(update_s%recv(m), recv_s)*ksize*nscalar
        m = ind_x(l)
        if(m>0) msgsize_x = get_mesgsize(update_x%recv(m), recv_v)*ksize*nvector
        m = ind_y(l)
        if(m>0) msgsize_y = get_mesgsize(update_y%recv(m), recv_v)*ksize*nvector
        msgsize = msgsize_s + msgsize_x + msgsize_y
        if( msgsize.GT.0 )then
           tot_recv_size = tot_recv_size + msgsize
           nrecv = nrecv + 1
           if(nrecv > MAXOVERLAP) then
              call mpp_error(FATAL, "set_group_update: nrecv is greater than MAXOVERLAP, increase MAXOVERLAP")
           endif
           group%from_pe(nrecv) = pelist(l)
           group%recv_size(nrecv) = msgsize
           group%buffer_pos_recv(nrecv) = update_buffer_pos
           update_buffer_pos = update_buffer_pos + msgsize
        end if
     end do
     group%nrecv = nrecv

     !--- setup for unpack
     nunpack = 0
     unpack_buffer_pos = 0
     do l = 1, nrecv_old
        m = ind_s(l)
        if(m>0) then
           overptr => update_s%recv(m)
           do n = 1, overptr%count
              dir = overptr%dir(n)
              if(recv_s(dir)) then
                 nunpack = nunpack + 1
                 if(nunpack > MAXOVERLAP) call mpp_error(FATAL,  &
                   "set_group_update: nunpack is greater than MAXOVERLAP, increase MAXOVERLAP 1")
                 group%unpack_type(nunpack) = FIELD_S
                 group%unpack_buffer_pos(nunpack) = unpack_buffer_pos
                 group%unpack_rotation(nunpack) = overptr%rotation(n) 
                 group%unpack_is(nunpack) = overptr%is(n)
                 group%unpack_ie(nunpack) = overptr%ie(n)
                 group%unpack_js(nunpack) = overptr%js(n)
                 group%unpack_je(nunpack) = overptr%je(n) 
                 group%unpack_size(nunpack) = overptr%msgsize(n)*nscalar
                 unpack_buffer_pos = unpack_buffer_pos + group%unpack_size(nunpack)*ksize
              end if
           end do   
        end if

        m = ind_x(l)
        if(m>0) then
           overptr => update_x%recv(m)
           do n = 1, overptr%count
              dir = overptr%dir(n)
              if(recv_v(dir)) then
                 nunpack = nunpack + 1
                 if(nunpack > MAXOVERLAP) call mpp_error(FATAL,  &
                   "set_group_update: nunpack is greater than MAXOVERLAP, increase MAXOVERLAP 2")
                 group%unpack_type(nunpack) = FIELD_X
                 group%unpack_buffer_pos(nunpack) = unpack_buffer_pos
                 group%unpack_rotation(nunpack) = overptr%rotation(n) 
                 group%unpack_is(nunpack) = overptr%is(n)
                 group%unpack_ie(nunpack) = overptr%ie(n)
                 group%unpack_js(nunpack) = overptr%js(n)
                 group%unpack_je(nunpack) = overptr%je(n) 
                 group%unpack_size(nunpack) = overptr%msgsize(n)*nvector
                 unpack_buffer_pos = unpack_buffer_pos + group%unpack_size(nunpack)*ksize
              end if
           end do   
        end if

        m = ind_y(l)
        if(m>0) then
           overptr => update_y%recv(m)
           do n = 1, overptr%count
              dir = overptr%dir(n)
              if(recv_v(dir)) then
                 nunpack = nunpack + 1
                 if(nunpack > MAXOVERLAP) call mpp_error(FATAL,  &
                   "set_group_update: nunpack is greater than MAXOVERLAP, increase MAXOVERLAP 3")
                 group%unpack_type(nunpack) = FIELD_Y
                 group%unpack_buffer_pos(nunpack) = unpack_buffer_pos
                 group%unpack_rotation(nunpack) = overptr%rotation(n) 
                 group%unpack_is(nunpack) = overptr%is(n)
                 group%unpack_ie(nunpack) = overptr%ie(n)
                 group%unpack_js(nunpack) = overptr%js(n)
                 group%unpack_je(nunpack) = overptr%je(n) 
                 group%unpack_size(nunpack) = overptr%msgsize(n)*nvector
                 unpack_buffer_pos = unpack_buffer_pos + group%unpack_size(nunpack)*ksize
              end if
           end do   
        end if
     end do
     group%nunpack = nunpack

     if(update_buffer_pos .NE. unpack_buffer_pos ) call mpp_error(FATAL,  &
                   "set_group_update: update_buffer_pos .NE. unpack_buffer_pos")

     !figure out message size for each processor.
     ntot = nsend_s + nsend_x + nsend_y
     n = 1
     i_s = 1
     i_x = 1
     i_y = 1
     ind_s = -1
     ind_x = -1
     ind_y = -1
     nsend = 0
     do while(n<=ntot)
        if( i_s <= nsend_s ) then
           rank_s = update_s%send(i_s)%pe-domain%pe
           if(rank_s .LT. 0) rank_s = rank_s + nlist
        else
           rank_s = nlist+1
        endif
        if( i_x <= nsend_x ) then
           rank_x = update_x%send(i_x)%pe-domain%pe
           if(rank_x .LT. 0) rank_x = rank_x + nlist
        else
           rank_x = nlist+1
        endif
        if( i_y <= nsend_y ) then
           rank_y = update_y%send(i_y)%pe-domain%pe
           if(rank_y .LT. 0) rank_y = rank_y + nlist
        else
           rank_y = nlist+1
        endif
        nsend = nsend + 1   
        rank = minval((/rank_s, rank_x, rank_y/))
        if(rank == rank_s) then
           n = n + 1
           ind_s(nsend) = i_s
           pelist(nsend) = update_s%send(i_s)%pe
           i_s = i_s + 1
        endif
        if(rank == rank_x) then
           n = n + 1
           ind_x(nsend) = i_x
           pelist(nsend) = update_x%send(i_x)%pe
           i_x = i_x + 1
        endif
        if(rank == rank_y) then
           n = n + 1
           ind_y(nsend) = i_y
           pelist(nsend) = update_y%send(i_y)%pe
           i_y = i_y + 1
        endif
     enddo

     nsend_old = nsend
     nsend     = 0
     tot_send_size = 0
     do l = 1, nsend_old
        msgsize_s = 0
        msgsize_x = 0
        msgsize_y = 0
        m = ind_s(l)
        if(m>0) msgsize_s = get_mesgsize(update_s%send(m), send_s)*ksize*nscalar
        m = ind_x(l)
        if(m>0) msgsize_x = get_mesgsize(update_x%send(m), send_v)*ksize*nvector
        m = ind_y(l)
        if(m>0) msgsize_y = get_mesgsize(update_y%send(m), send_v)*ksize*nvector
        msgsize = msgsize_s + msgsize_x + msgsize_y
        if( msgsize.GT.0 )then
           tot_send_size = tot_send_size + msgsize
           nsend = nsend + 1
           if(nsend > MAXOVERLAP) then
              call mpp_error(FATAL, "set_group_update: nsend is greater than MAXOVERLAP, increase MAXOVERLAP")
           endif
           send_size(nsend) = msgsize
           group%to_pe(nsend) = pelist(l)
           group%buffer_pos_send(nsend) = update_buffer_pos
           group%send_size(nsend) = msgsize
           update_buffer_pos = update_buffer_pos + msgsize
        end if
     end do
     group%nsend = nsend

     !--- setup for pack
     npack = 0
     pack_buffer_pos = unpack_buffer_pos
     do l = 1, nsend_old
        m = ind_s(l)
        if(m>0) then
           overptr => update_s%send(m)
           do n = 1, overptr%count
              dir = overptr%dir(n)
              if(send_s(dir)) then
                 npack = npack + 1
                 if(npack > MAXOVERLAP) call mpp_error(FATAL,  &
                   "set_group_update: npack is greater than MAXOVERLAP, increase MAXOVERLAP 1")
                 group%pack_type(npack) = FIELD_S
                 group%pack_buffer_pos(npack) = pack_buffer_pos
                 group%pack_rotation(npack) = overptr%rotation(n) 
                 group%pack_is(npack) = overptr%is(n)
                 group%pack_ie(npack) = overptr%ie(n)
                 group%pack_js(npack) = overptr%js(n)
                 group%pack_je(npack) = overptr%je(n) 
                 group%pack_size(npack) = overptr%msgsize(n)*nscalar
                 pack_buffer_pos = pack_buffer_pos + group%pack_size(npack)*ksize
              end if
           end do   
        end if

        m = ind_x(l)
        if(m>0) then
           overptr => update_x%send(m)
           do n = 1, overptr%count
              dir = overptr%dir(n)
              if(send_v(dir)) then
                 npack = npack + 1
                 if(npack > MAXOVERLAP) call mpp_error(FATAL,  &
                   "set_group_update: npack is greater than MAXOVERLAP, increase MAXOVERLAP 2")
                 group%pack_type(npack) = FIELD_X
                 group%pack_buffer_pos(npack) = pack_buffer_pos
                 group%pack_rotation(npack) = overptr%rotation(n) 
                 group%pack_is(npack) = overptr%is(n)
                 group%pack_ie(npack) = overptr%ie(n)
                 group%pack_js(npack) = overptr%js(n)
                 group%pack_je(npack) = overptr%je(n) 
                 group%pack_size(npack) = overptr%msgsize(n)*nvector
                 pack_buffer_pos = pack_buffer_pos + group%pack_size(npack)*ksize
              end if
           end do   
        end if

        m = ind_y(l)
        if(m>0) then
           overptr => update_y%send(m)
           do n = 1, overptr%count
              dir = overptr%dir(n)
              if(send_v(dir)) then
                 npack = npack + 1
                 if(npack > MAXOVERLAP) call mpp_error(FATAL,  &
                   "set_group_update: npack is greater than MAXOVERLAP, increase MAXOVERLAP 3")
                 group%pack_type(npack) = FIELD_Y
                 group%pack_buffer_pos(npack) = pack_buffer_pos
                 group%pack_rotation(npack) = overptr%rotation(n) 
                 group%pack_is(npack) = overptr%is(n)
                 group%pack_ie(npack) = overptr%ie(n)
                 group%pack_js(npack) = overptr%js(n)
                 group%pack_je(npack) = overptr%je(n) 
                 group%pack_size(npack) = overptr%msgsize(n)*nvector
                 pack_buffer_pos = pack_buffer_pos + group%pack_size(npack)*ksize
              end if
           end do   
        end if
     end do
     group%npack = npack
     if(update_buffer_pos .NE. pack_buffer_pos ) call mpp_error(FATAL,  &
                   "set_group_update: update_buffer_pos .NE. pack_buffer_pos")

     !--- make sure the buffer is large enough
     mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, tot_recv_size+tot_send_size )

     if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
        write( text,'(i8)' )mpp_domains_stack_hwm
        call mpp_error( FATAL, 'set_group_update: mpp_domains_stack overflow, '// &
                     'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
     end if

     group%tot_msgsize = tot_recv_size+tot_send_size

end subroutine set_group_update


!######################################################################
  subroutine mpp_clear_group_update(group)
     type(mpp_group_update_type), intent(inout) :: group

     group%nscalar = 0
     group%nvector = 0
     group%nsend   = 0
     group%nrecv   = 0
     group%npack   = 0
     group%nunpack = 0
     group%initialized = .false.

  end subroutine mpp_clear_group_update

!#####################################################################
  function mpp_group_update_initialized(group)
     type(mpp_group_update_type), intent(in) :: group
     logical :: mpp_group_update_initialized

     mpp_group_update_initialized = group%initialized

  end function mpp_group_update_initialized

!#####################################################################
  function mpp_group_update_is_set(group)
     type(mpp_group_update_type), intent(in) :: group
     logical :: mpp_group_update_is_set

     mpp_group_update_is_set = (group%nscalar > 0 .OR. group%nvector > 0)

  end function mpp_group_update_is_set



