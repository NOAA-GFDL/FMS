  subroutine mpp_pset_create(npset,pset,stacksize,pelist)
!create PSETs
!  called by all PEs in parent pelist
!  mpset must be exact divisor of npes
    integer, intent(in) :: npset !number of PSETs per set
    type(mpp_pset_type), intent(inout) :: pset
    integer, intent(in), optional :: stacksize
    integer, intent(in), optional :: pelist(:)

    integer :: npes
    integer :: i, j, k, out_unit, errunit
    integer, allocatable :: my_pelist(:), root_pelist(:)

    call mpp_init()
    call mpp_pset_init()

#ifdef PSET_DEBUG
    verbose=.TRUE.
#endif
    out_unit = stdout()
    errunit  = stderr()
    pe = mpp_pe()
    if(present(pelist)) then
       npes = size(pelist(:))
    else
       npes = mpp_npes()
    endif
    if( mod(npes,npset).NE.0 )then
        write( text,'(a,2i6)' ) &
             'MPP_PSET_CREATE: PSET size (npset) must divide npes exactly:'// &
             ' npset, npes=', npset, npes
        call mpp_error( FATAL, text )
    end if

    !configure out root_pelist
    allocate(my_pelist(0:npes-1) )
    allocate(root_pelist(0:npes/npset-1) )
    if(present(pelist)) then
       my_pelist = pelist
    else
       call mpp_get_current_pelist(my_pelist)
    endif
    do i = 0,npes/npset-1
       root_pelist(i) = my_pelist(npset*i)
    enddo
    write( out_unit,'(a,i6)' )'MPP_PSET_CREATE creating PSETs... npset=', npset
    if(ANY(my_pelist == pe) ) then
    if( pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_CREATE: PSET already initialized!' )
    pset%npset = npset
    allocate( pset%pelist(0:npes-1) )
    allocate( pset%root_pelist(0:npes/npset-1) )
    pset%pelist = my_pelist
!create the root PElist
    pset%root_pelist = root_pelist
    allocate( pset%pset(0:npset-1) )
    do i = 0,npes/npset-1
       k = npset*i
!designate the root PE, next PE, prev PE
       do j = 0,npset-1
          if( pe.EQ.pset%pelist(k+j) )then
              pset%pset(:) =  pset%pelist(k:k+npset-1)
              pset%pos = j
              pset%root_in_pset = pset%root_pelist(i)
              if( j.EQ.0 )then
                  pset%prev_in_pset = pset%pelist(k+npset-1)
              else
                  pset%prev_in_pset = pset%pelist(k+j-1)
              end if
              if( j.EQ.npset-1 )then
                  pset%next_in_pset = pset%pelist(k)
              else
                  pset%next_in_pset = pset%pelist(k+j+1)
              end if
          end if
       end do
    end do

    pset%root = pe.EQ.pset%root_in_pset

!stack
    pset%hiWM = 0 !initialize hi-water-mark
    pset%maxstack = 1000000 !default
    if( PRESENT(stacksize) )pset%maxstack = stacksize
    write( out_unit,'(a,i8)' ) &
         'MPP_PSET_CREATE: setting stacksize=', pset%maxstack
    if( pset%root )then
        allocate( pset%stack(pset%maxstack) )
#ifdef use_CRI_pointers
        pset%p_stack = LOC(pset%stack)
#endif
    end if
    pset%initialized = .TRUE. !must be called before using pset
    call mpp_pset_broadcast_ptr(pset,pset%p_stack)
    endif

    call mpp_declare_pelist(root_pelist)

    if( verbose )then
        write( errunit,'(a,4i6)' )'MPP_PSET_CREATE: pe, root, next, prev=', &
             pe, pset%root_in_pset, pset%next_in_pset, pset%prev_in_pset
        write( errunit,* )'PE ', pe, ' pset=', pset%pset(:)
        write( out_unit,* )'root pelist=', pset%root_pelist(:)
    end if
  end subroutine mpp_pset_create

  subroutine mpp_pset_get_root_pelist(pset,pelist)
    type(mpp_pset_type), intent(in) :: pset
    integer, intent(out) :: pelist(:)

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_GET_ROOT_PELIST: called with uninitialized PSET.' )
    if( size(pelist).NE.size(pset%root_pelist) )then
        write( text,'(a,2i6)' ) &
             'pelist argument has wrong size: requested, actual=', &
             size(pelist), size(pset%root_pelist)
        call mpp_error( FATAL, 'MPP_PSET_GET_ROOT_PELIST: '//text )
    end if
    pelist(:) = pset%root_pelist(:)
  end subroutine mpp_pset_get_root_pelist
