
! -*-f90-*-
! $Id: mpp_domains_comm.inc,v 13.0 2006/03/28 21:40:34 fms Exp $

#ifdef use_GSM
#include "mpp/shmem.fh"
#endif

  ! <INTERFACE NAME="mpp_domains_init_comm">
  !  <OVERVIEW>
  !    Instantiate mpp_domains_communicator derived type
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The call returns the allocated and filled mpp_communictor derived type
  !    This data object is designed to allows all data start and end points
  !    to be calculated once and re-used in an domains update, update_V and
  !    redistribute operation
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_domains_update_init_comm(d_comm,domain,flags )
  !  </TEMPLATE>
  !  <OUT NAME="d_comm"></OUT>
  !  <IN NAME="domain"></IN>
  !  <IN NAME="flags"></IN>
  ! </INTERFACE>

  !--- we may provide seperate interface for scalar and vector field.
  function mpp_update_init_comm(domain,l_addr,isize,jsize,ksize,l_addr2,flags, position, gridtype) RESULT(d_comm)
      type(DomainCommunicator2D), pointer       :: d_comm
      type(domain2D),target,      intent(inout) :: domain  ! use of compute_overlaps requires intent(inout)
      integer(LONG_KIND),         intent(in)    :: l_addr(:)
      integer,                    intent(in)    :: isize(:)
      integer,                    intent(in)    :: jsize(:)
      integer,                    intent(in)    :: ksize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2(:)
      integer, optional,          intent(in)    :: flags
      integer, optional,          intent(in)    :: position
      integer, optional,          intent(in)    :: gridtype

      integer(LONG_KIND) :: domain_id
      integer :: m, n, list, l
      integer :: update_flags
      integer :: is, ie, js, je, ioff, joff, list_size
      integer :: i, lsize,rsize,msgsize,to_pe,from_pe
      integer :: num_component, x_position, y_position
      integer, allocatable :: isL(:,:), jsL(:,:)
      logical :: recv(8), send(8)
      logical :: folded
      real,dimension(*) :: f, f2
      pointer(ptr,f);  pointer(ptr2,f2)
      real,dimension(*) :: f_r, f2_r
      pointer(ptr_r,f_r);  pointer(ptr2_r,f2_r)
      logical :: use_shmem_ptr
      character(len=8) :: text
      type(DomainCommunicator2D), pointer :: d_comm_local => NULL()
      type(domain2d),             pointer :: Dom          => NULL()

      use_shmem_ptr=.false.   ! False for any but GSM
#ifdef use_GSM
      call check_sma_env()
      ptr=l_addr(1); if(PRESENT(l_addr2))ptr2=l_addr2(1)  ! Set up target for shmem_ptr
      use_shmem_ptr=.false.; ptr_r=shmem_ptr(f,mpp_pe()); if(ptr_r>0)use_shmem_ptr=.true.
#endif

!for all gridtypes
      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) ) then
          update_flags = flags
          ! The following test is so that SCALAR_PAIR can be used alone with the
          ! same default update pattern as without.
          if (BTEST(update_flags,SCALAR_BIT)) then
            if (.NOT.(BTEST(update_flags,WEST) .OR. BTEST(update_flags,EAST) &
                 .OR. BTEST(update_flags,NORTH) .OR. BTEST(update_flags,SOUTH))) &
              update_flags = update_flags + XUPDATE+YUPDATE   !default with SCALAR_PAIR
          end if
      end if

    ! Create unique domain identifier
    list_size = size(l_addr(:))
    domain_id=set_domain_id(domain%id,ksize+list_size,update_flags, gridtype, position )

    if(PRESENT(l_addr2))then
       d_comm =>get_comm(domain_id,l_addr(1),l_addr2(1))
    else
       d_comm =>get_comm(domain_id,l_addr(1))
    endif

    if(d_comm%initialized)return   ! Found existing field/domain communicator

    !--- position and gridtype can not be presented at the same time
    if(PRESENT(position) .AND. PRESENT(gridtype)) then
       call mpp_error( FATAL, 'mpp_update_comm_mod: position and gridtype can not both be presented')
    else if (PRESENT(position)) then
       x_position = position; y_position = position
       !--- for cubic grid, can not use scalar version update_domain for data on E or N-cell
       if( domain%topology_type == CUBIC_GRID .AND. ( position == EAST .OR. position == NORTH ) ) then
          call mpp_error(FATAL, 'mpp_domains_comm_mod: for cubic grid, ' // &
               'can not use scalar version update_domain for data on E or N-cell' )
       end if
    else if ( PRESENT(gridtype)) then
       if( BTEST(update_flags,WEST) .AND. BTEST(domain%fold,WEST)  .AND. BTEST(gridtype,EAST)  ) &
            call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
       if( BTEST(update_flags,SOUTH) .AND. BTEST(domain%fold,SOUTH) .AND. BTEST(gridtype,NORTH) ) &
            call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
       if( BTEST(update_flags,EAST) .AND. BTEST(domain%fold,EAST)  .AND. BTEST(gridtype,WEST)  ) &
            call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
       if( BTEST(update_flags,NORTH) .AND. BTEST(domain%fold,NORTH) .AND. BTEST(gridtype,SOUTH) ) &
            call mpp_error( FATAL, 'Incompatible grid offset and fold.' )

       !--- For cubic grid, velocity is not allowed on the CORNER.
       if( domain%topology_type == CUBIC_GRID .AND. ( .NOT. BTEST(update_flags,SCALAR_BIT)) .AND. &
            ( gridtype == BGRID_NE .OR. gridtype == BGRID_SW ) ) then
          call mpp_error(FATAL,"mpp_update_comm_mod: velocity is not allowed on the C-cell for cubic grid")
       end if

       select case(gridtype)
       case (AGRID)
          x_position = CENTER; y_position = CENTER
       case (BGRID_NE, BGRID_SW)  ! acutally not quite right for BGRID_SW, BGRID_SW is not well implemented.
          x_position = CORNER; y_position = CORNER
       case (CGRID_NE, CGRID_SW)  ! acutally not quite right for CGRID_SW, CGRID_SW is not well implemented.
          x_position = EAST; y_position = NORTH
       case default
          call mpp_error( FATAL, 'mpp_update_comm_mod: gridtype must be one of ' //&
               'AGRID|BGRID_NE|BGRID_SW|CGRID_NE|CGRID_SW.' )
       end select
    else
       call mpp_error( FATAL, 'mpp_update_comm_mod: one of position and gridtype should be presented')
    end if

    recv(1) = BTEST(update_flags,EAST)
    recv(3) = BTEST(update_flags,SOUTH)
    recv(5) = BTEST(update_flags,WEST)
    recv(7) = BTEST(update_flags,NORTH)
    recv(2) = recv(1) .AND. recv(3)
    recv(4) = recv(3) .AND. recv(5)
    recv(6) = recv(5) .AND. recv(7)
    recv(8) = recv(7) .AND. recv(1)
    send    = recv

    if(x_position == y_position) then
      num_component = 1
      d_comm%staggered = .FALSE.
    else 
         num_component = 2
         d_comm%staggered = .TRUE.
      endif

      do n = 1, num_component
         if( n == 1 ) then
            d_comm_local => d_comm
            d_comm_local%position = x_position
         else
            allocate(d_comm%y_comm)
            d_comm_local => d_comm%y_comm
            d_comm_local%position = y_position
         endif
       select case(d_comm_local%position)
       case (CENTER)
          Dom => domain%T
       case (EAST)
          Dom => domain%E
       case (NORTH)
          Dom => domain%N
       case (CORNER)
          Dom => domain%C
       end select

         d_comm_local%l_addr = l_addr(1)
         d_comm_local%domain =>domain

       d_comm_local%Slist_size = size(Dom%list(:))
       d_comm_local%isize = isize(n)
       d_comm_local%jsize = jsize(n)
       d_comm_local%ke = ksize

!send
         lsize = d_comm_local%Slist_size-1
       allocate( d_comm_local%S_msize(8,0:lsize),   d_comm_local%S_msize2(8,0:lsize), &
                 d_comm_local%do_thisS(8,0:lsize),  d_comm_local%S_do_buf(0:lsize),   &
                 d_comm_local%do_thisS2(8,0:lsize),                                   &
                 d_comm_local%do_thisS3(8,0:lsize), d_comm_local%S_do_buf3(0:lsize),  &
                 d_comm_local%do_thisS4(8,0:lsize), d_comm_local%S_do_buf4(0:lsize),  &
                 d_comm_local%send(8,0:lsize),      d_comm_local%recv(8,0:lsize),     &
                 d_comm_local%cto_pe(0:lsize),      isL(8,0:lsize),jsL(8,0:lsize))
#ifdef use_CAF
         allocate(d_comm_local%Rcaf_idx(0:lsize))
         d_comm_local%Rcaf_idx = 0
#endif
      isL=0;jsL=0
         d_comm_local%cto_pe=-1
         d_comm_local%S_msize=0
         d_comm_local%S_do_buf=.false.
         d_comm_local%do_thisS=.false.
       d_comm_local%S_msize2=0
       d_comm_local%do_thisS2=.false.
       d_comm_local%S_do_buf3=.false.
       d_comm_local%do_thisS3=.false.
       d_comm_local%S_do_buf4=.false.
       d_comm_local%do_thisS4=.false.
      ioff = domain%x%data%begin
      joff = domain%y%data%begin

      do list = 0,lsize
         m = mod( domain%pos+list, lsize+1 )
          d_comm_local%cto_pe(list) = Dom%list(m)%pe
          if( .NOT. ANY(Dom%list(m)%overlap) )cycle

            to_pe = d_comm_local%cto_pe(list)
#ifdef use_CAF
            d_comm_local%Rcaf_idx(to_pe) = list   ! local CAF pe needs to know it's position in remote CAF list
#endif
            do l = 1, 8  ! loop over 8 direction
             if(send(l)) then
                if( Dom%list(m)%send(l)%overlap(1) ) then
                  d_comm_local%do_thisS(l,list)=.true.
                   is = Dom%list(m)%send(l)%is;       ie = Dom%list(m)%send(l)%ie
                   js = Dom%list(m)%send(l)%js;       je = Dom%list(m)%send(l)%je
                   d_comm_local%send(l,list)%is = is; d_comm_local%send(l,list)%ie = ie
                   d_comm_local%send(l,list)%js = js; d_comm_local%send(l,list)%je = je
                   d_comm_local%S_msize(l,list) = (ie-is+1)*(je-js+1)*ksize
                   d_comm_local%send(l,list)%rotation = Dom%list(m)%send(l)%rotation 
                end if
                if( Dom%list(m)%send(l)%overlap(2) ) then
                   d_comm_local%do_thisS2(l,list)=.true.
                   d_comm_local%send(l,list)%n = Dom%list(m)%send(l)%n
                   d_comm_local%send(l,list)%i => Dom%list(m)%send(l)%i
                   d_comm_local%send(l,list)%j => Dom%list(m)%send(l)%j
                   d_comm_local%S_msize(l,list) = d_comm_local%S_msize(l,list)+d_comm_local%send(l,list)%n*ksize 
                end if
                if( Dom%list(m)%send(l)%overlap(3) ) then
                   is = Dom%list(m)%send(l)%is; ie = Dom%list(m)%send(l)%ie
                   js = Dom%list(m)%send(l)%js; je = Dom%list(m)%send(l)%je
                   d_comm_local%do_thisS3(l,list)=.true.
                   d_comm_local%send(l,list)%is = is; d_comm_local%send(l,list)%ie = ie
                   d_comm_local%send(l,list)%js = js; d_comm_local%send(l,list)%je = je
                   d_comm_local%S_msize2(l,list) = (ie-is+1)*(je-js+1)*ksize
                   d_comm_local%send(l,list)%rotation = Dom%list(m)%send(l)%rotation
                   !               isL(l,list) = is-ioff+1; jsL(l,list) = js-joff+1
                end if
                if( Dom%list(m)%send(l)%overlap(4) ) then
                   d_comm_local%do_thisS4(l,list)=.true.
                   d_comm_local%send(l,list)%i2 = Dom%list(m)%send(l)%i2
                   d_comm_local%send(l,list)%j2 = Dom%list(m)%send(l)%j2
                end if
             end if
            enddo
          if(ANY(d_comm_local%do_thisS (:,list)))d_comm_local%S_do_buf (list) = .true.
          if(ANY(d_comm_local%do_thisS2(:,list)))d_comm_local%S_do_buf (list) = .true.
          if(ANY(d_comm_local%do_thisS3(:,list)))d_comm_local%S_do_buf3(list) = .true.
          if(ANY(d_comm_local%do_thisS4(:,list)))d_comm_local%S_do_buf4(list) = .true.
#ifdef use_GSM
            if(d_comm_local%S_do_buf(list))then
           if(.not.use_shmem_ptr)then  ! address is from GLOBAL_ALLOC
             if(.not.PRESENT(l_addr2))then
               call mpp_send(l_addr(1), plen=list_size, to_pe=to_pe )
             else
               call mpp_send(l_addr(1), plen=list_size, to_pe=to_pe )
               call mpp_send(l_addr2(1), plen=list_size, to_pe=to_pe )
             endif
           endif
           call mpp_send(isL(1,list), plen=8, to_pe=to_pe )
           call mpp_send(jsL(1,list), plen=8, to_pe=to_pe )
             call mpp_send(isize(n), plen=1, to_pe=to_pe )
             call mpp_send(jsize(n), plen=1, to_pe=to_pe )
         endif
#endif
      end do

      call mpp_sync_self()
!recv
       d_comm_local%Rlist_size = size(Dom%list(:))

         rsize = d_comm_local%Rlist_size-1
       allocate(d_comm_local%R_msize(8,0:rsize),   d_comm_local%R_msize2(8,0:rsize), &
                d_comm_local%do_thisR(8,0:rsize),  d_comm_local%R_do_buf(0:rsize),   &
                d_comm_local%do_thisR2(8,0:rsize),                                   &
                d_comm_local%do_thisR3(8,0:rsize), d_comm_local%R_do_buf3(0:rsize),  &
                d_comm_local%do_thisR4(8,0:rsize), d_comm_local%R_do_buf4(0:rsize),  &
                d_comm_local%isizeR(0:rsize),      d_comm_local%jsizeR(0:rsize),     &
                d_comm_local%sendisR(8,0:rsize),   d_comm_local%sendjsR(8,0:rsize),  &
                d_comm_local%cfrom_pe(0:rsize) )
#ifdef use_GSM
      if(.not.PRESENT(l_addr2))then                                 
            allocate(d_comm_local%rem_addrl(list_size,0:rsize))
            d_comm_local%rem_addrl=-9999
      else
            allocate(d_comm_local%rem_addrlx(list_size,0:rsize),d_comm_local%rem_addrly(list_size,0:rsize))
            d_comm_local%rem_addrlx=-9999; d_comm_local%rem_addrly=-9999
      endif
#endif

         d_comm_local%cfrom_pe=-1
         d_comm_local%R_msize=0
         d_comm_local%R_do_buf=.false.
         d_comm_local%do_thisR=.false.
       d_comm_local%do_thisR2 = .false.
       d_comm_local%R_do_buf3 = .false.
       d_comm_local%do_thisR3 = .false.
       d_comm_local%R_do_buf4 = .false.
       d_comm_local%do_thisR4 = .false.
       d_comm_local%R_msize2=0
         d_comm_local%isizeR=0; d_comm_local%jsizeR=0
         d_comm_local%sendisR=0; d_comm_local%sendjsR=0
      do list = 0,rsize
          m = mod( Dom%pos+rsize+1-list, rsize+1 )
          d_comm_local%cfrom_pe(list) = Dom%list(m)%pe
          if( .NOT. ANY(Dom%list(m)%overlap) )cycle
          from_pe = d_comm_local%cfrom_pe(list)

            do l = 1, 8  ! loop over 8 directions
             if(recv(l)) then
                if( Dom%list(m)%recv(l)%overlap(1) ) then
                   d_comm_local%do_thisR(l,list)=.true.
                   is = Dom%list(m)%recv(l)%is; ie = Dom%list(m)%recv(l)%ie
                   js = Dom%list(m)%recv(l)%js; je = Dom%list(m)%recv(l)%je
                   d_comm_local%recv(l,list)%is = is; d_comm_local%recv(l,list)%ie = ie
                   d_comm_local%recv(l,list)%js = js; d_comm_local%recv(l,list)%je = je
                   d_comm_local%R_msize(l,list) = (ie-is+1)*(je-js+1)*ksize
                   d_comm_local%recv(l,list)%rotation = Dom%list(m)%recv(l)%rotation
                end if
                if( Dom%list(m)%recv(l)%overlap(2) ) then
                   d_comm_local%do_thisR2(l,list)=.true.
                   d_comm_local%recv(l,list)%n = Dom%list(m)%recv(l)%n
                   d_comm_local%recv(l,list)%i => Dom%list(m)%recv(l)%i
                   d_comm_local%recv(l,list)%j => Dom%list(m)%recv(l)%j
                   d_comm_local%R_msize(l,list) = d_comm_local%R_msize(l,list)+d_comm_local%recv(l,list)%n*ksize
                end if
                if( Dom%list(m)%recv(l)%overlap(3) ) then
                   is = Dom%list(m)%recv(l)%is; ie = Dom%list(m)%recv(l)%ie
                   js = Dom%list(m)%recv(l)%js; je = Dom%list(m)%recv(l)%je
                   d_comm_local%do_thisR3(l,list)=.true.
                   d_comm_local%recv(l,list)%is = is; d_comm_local%recv(l,list)%ie = ie
                   d_comm_local%recv(l,list)%js = js; d_comm_local%recv(l,list)%je = je
                   d_comm_local%R_msize2(l,list) = (ie-is+1)*(je-js+1)*ksize
                   d_comm_local%recv(l,list)%rotation = Dom%list(m)%recv(l)%rotation
                end if
                if( Dom%list(m)%recv(l)%overlap(4) ) then
                   d_comm_local%do_thisR4(l,list)=.true.
                   d_comm_local%recv(l,list)%i2 = Dom%list(m)%recv(l)%i2
                   d_comm_local%recv(l,list)%j2 = Dom%list(m)%recv(l)%j2
                end if
             end if
            enddo

         if(ANY(d_comm_local%do_thisR (:,list)))d_comm_local%R_do_buf (list) = .true.
         if(ANY(d_comm_local%do_thisR2(:,list)))d_comm_local%R_do_buf (list) = .true.
         if(ANY(d_comm_local%do_thisR3(:,list)))d_comm_local%R_do_buf3(list) = .true.
         if(ANY(d_comm_local%do_thisR4(:,list)))d_comm_local%R_do_buf4(list) = .true.
#ifdef use_GSM
            if(d_comm_local%R_do_buf(list))then
           if(.not.PRESENT(l_addr2))then
             if(use_shmem_ptr)then                          
               do i=1,list_size   
                 ptr = l_addr(i)
                        d_comm_local%rem_addrl(i,list) = shmem_ptr(f,from_pe)
               end do
             else
                     call mpp_recv(d_comm_local%rem_addrl(1,list), glen=list_size, from_pe=from_pe )
             endif
           else
             if(use_shmem_ptr)then
               do i=1,list_size
                 ptr  = l_addr(i)
                 ptr2 = l_addr2(i)
                        d_comm_local%rem_addrlx(i,list) = shmem_ptr(f,from_pe)
                        d_comm_local%rem_addrly(i,list) = shmem_ptr(f2,from_pe)
               end do
             else    
                     call mpp_recv(d_comm_local%rem_addrlx(1,list), glen=list_size, from_pe=from_pe )
                     call mpp_recv(d_comm_local%rem_addrly(1,list), glen=list_size, from_pe=from_pe )
             endif
           endif
               call mpp_recv(d_comm_local%sendisR(1,list), glen=8, from_pe=from_pe )
               call mpp_recv(d_comm_local%sendjsR(1,list), glen=8, from_pe=from_pe )
               call mpp_recv(d_comm_local%isizeR(list), glen=1, from_pe=from_pe )
               call mpp_recv(d_comm_local%jsizeR(list), glen=1, from_pe=from_pe )
         endif
#endif
      end do

      d_comm_local%isize_max = isize(n); call mpp_max(d_comm_local%isize_max)
      d_comm_local%jsize_max = jsize(n); call mpp_max(d_comm_local%jsize_max)

      DEALLOCATE(isL,jsL)
      d_comm_local%initialized = .true.
      d_comm_local => NULL()
   enddo

!--- check the stack size
#if (!defined(use_GSM)) && (!defined(use_CAF))
    ! Handles case where S_msize and/or R_msize are 0 size array
      msgsize = ( MAXVAL( (/0,sum(d_comm%S_msize(:,:))/) ) + MAXVAL( (/0,sum(d_comm%R_msize(:,:))/) ) ) * list_size
      if( PRESENT(l_addr2) ) then
         if( Associated(d_comm%y_comm) ) then
            msgsize = msgsize + ( MAXVAL( (/0,sum(d_comm%y_comm%S_msize(:,:))/) ) + &
                      MAXVAL( (/0,sum(d_comm%y_comm%R_msize(:,:))/) ) ) * list_size
         else
            msgsize = msgsize*2
         end if
      end if
      if(msgsize>0)then
         mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, msgsize )
         if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
             write( text,'(i8)' )mpp_domains_stack_hwm
             call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                  //trim(text)//') from all PEs.' )
         end if
      end if
#endif

    end function mpp_update_init_comm


    function mpp_redistribute_init_comm(domain_in,l_addrs_in, domain_out,l_addrs_out, &
                                        isize_in,jsize_in,ksize_in,isize_out,jsize_out,ksize_out) RESULT(d_comm)
      type(DomainCommunicator2D), pointer       :: d_comm
      type(domain2D),target,      intent(in)    :: domain_in
      integer(LONG_KIND),         intent(in)    :: l_addrs_in(:)
      type(domain2D),target,      intent(in)    :: domain_out
      integer(LONG_KIND),         intent(in)    :: l_addrs_out(:)
      integer,                    intent(in)    :: isize_in
      integer,                    intent(in)    :: jsize_in
      integer,                    intent(in)    :: ksize_in
      integer,                    intent(in)    :: isize_out 
      integer,                    intent(in)    :: jsize_out
      integer,                    intent(in)    :: ksize_out 

      integer(LONG_KIND) :: domain_id
      integer :: m, list
      integer :: is, ie, js, je, ke, ioff, joff, list_size
      integer :: isc, iec, jsc, jec
      integer :: i, lsize,rsize,msgsize,to_pe,from_pe
      integer,           allocatable,dimension(:) :: isL, jsL
      integer(LONG_KIND),allocatable,dimension(:,:) :: slist_addr
      real,dimension(*) :: f
      pointer(ptr,f)
      real,dimension(*) :: f_r
      pointer(ptr_r,f_r)
      logical :: use_shmem_ptr
      character(len=8) :: text


    ! This test determines whether input fields are from allocated memory (LOC gets global
    ! address) or "static" memory (need shmem_ptr). This probably needs to be generalized
    ! to determine appropriate mechanism for each incoming address.

    ! "Concurrent" run mode may leave field_in or field_out unassociated if pe does not
    ! contain in/out data. Use of STATIC option for ocean complicates this as ocean component
    ! always defined. Field_out is always a boundary structure and so is always allocated or
    ! not depending on whether it's used. If field out is defined (>0), then it is used otherwise
    ! field in must be defined.

      use_shmem_ptr=.false.   ! False for any but use_GSM
#ifdef use_GSM
      call check_sma_env()
      if(l_addrs_in(1)>0)then  ! Set up target for shmem_ptr
        ptr=l_addrs_in(1)
        ptr_r=shmem_ptr(f,mpp_pe())
        if(ptr_r>0)use_shmem_ptr=.true.
      endif
#endif

!fix ke
      ke = 0
      if( domain_in%pe /= NULL_PE )ke = ksize_in
      if( domain_out%pe /= NULL_PE )then
          if( ke /= 0 .AND. ke /= ksize_out ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: mismatch between field_in and field_out.' )
          ke = ksize_out
      end if
      if( ke == 0 )call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: either domain_in or domain_out must be native.' )
!check sizes
      if( domain_in%pe /= NULL_PE )then
          if( isize_in /= domain_in%x%data%size .OR. jsize_in /= domain_in%y%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: field_in must be on data domain of domain_in.' )
      end if
      if( domain_out%pe /= NULL_PE )then
          if( isize_out /= domain_out%x%data%size .OR. jsize_out /= domain_out%y%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: field_out must be on data domain of domain_out.' )
      end if


    ! Create unique domain identifier
      list_size = size(l_addrs_in(:))
      if(l_addrs_out(1) > 0)then
        domain_id = set_domain_id(domain_out%id,ke+list_size)
      else
        domain_id = set_domain_id(domain_in%id,ke+list_size)
      endif

      d_comm =>get_comm(domain_id,l_addrs_in(1),l_addrs_out(1))

      if(d_comm%initialized)return   ! Found existing field/domain communicator

      d_comm%l_addr = l_addrs_in(1)
      d_comm%domain_in =>domain_in
      d_comm%Slist_size = size(domain_out%list(:))
      d_comm%isize_in = isize_in
      d_comm%jsize_in = jsize_in
      d_comm%ke = ke

!send 
      lsize = d_comm%Slist_size-1
      allocate(d_comm%sendis(1,0:lsize), d_comm%sendie(1,0:lsize), &
               d_comm%sendjs(1,0:lsize), d_comm%sendje(1,0:lsize), &
               d_comm%S_msize(1,0:lsize),isL(0:lsize),jsL(0:lsize))
      allocate(slist_addr(list_size,0:lsize))
      allocate(d_comm%cto_pe(0:lsize), d_comm%S_do_buf(0:lsize))
#ifdef use_CAF
      allocate(d_comm%Rcaf_idx(0:lsize))
      d_comm%Rcaf_idx = 0
#endif
      isL=0;jsL=0
      slist_addr = -9999
      d_comm%cto_pe=-1
      d_comm%sendis=0; d_comm%sendie=0
      d_comm%sendjs=0; d_comm%sendje=0;
      d_comm%S_msize=0
      d_comm%S_do_buf=.false.

      ioff = domain_in%x%data%begin
      joff = domain_in%y%data%begin

      call mpp_get_compute_domain( domain_in, isc, iec, jsc, jec )
      do list = 0,lsize
         m = mod( domain_out%pos+list+lsize+1, lsize+1 )
         d_comm%cto_pe(list) = domain_out%list(m)%pe
         to_pe =  d_comm%cto_pe(list)
#ifdef use_CAF
         d_comm%Rcaf_idx(to_pe) = list   ! local CAF pe needs to know it's position in remote CAF list
#endif
         call mpp_get_compute_domain( domain_out%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie >= is .AND. je >= js )then
           d_comm%S_do_buf(list) = .true.
           d_comm%sendis(1,list)=is; d_comm%sendie(1,list)=ie
           d_comm%sendjs(1,list)=js; d_comm%sendje(1,list)=je
           d_comm%S_msize(1,list) = (ie-is+1)*(je-js+1)*ke
           isL(list) = is-ioff+1; jsL(list) = js-joff+1
#ifdef use_GSM
           if(.not.use_shmem_ptr)then  ! address is from GLOBAL_ALLOC
             slist_addr(1:list_size,list) = l_addrs_in(1:list_size)
             call mpp_send(slist_addr(1,list), plen=list_size, to_pe=to_pe )
           endif
           call mpp_send(isL(list), plen=1, to_pe=to_pe )
           call mpp_send(jsL(list), plen=1, to_pe=to_pe )
           call mpp_send(isize_in, plen=1, to_pe=to_pe )   ! I extent of remote array
           call mpp_send(jsize_in, plen=1, to_pe=to_pe )   ! J extent of remote array
#endif
         end if
      end do 

      call mpp_sync_self()
!recv 
      d_comm%domain_out =>domain_out
      d_comm%Rlist_size = size(domain_in%list(:))
      d_comm%isize_out = isize_out
      d_comm%jsize_out = jsize_out

      rsize = d_comm%Rlist_size-1
      allocate(d_comm%recvis(1,0:rsize), d_comm%recvie(1,0:rsize), &
               d_comm%recvjs(1,0:rsize), d_comm%recvje(1,0:rsize), &
               d_comm%R_msize(1,0:rsize))
      allocate(d_comm%cfrom_pe(0:rsize), d_comm%R_do_buf(0:rsize))
      allocate(d_comm%isizeR(0:rsize), d_comm%jsizeR(0:rsize))
      allocate(d_comm%sendisR(1,0:rsize), d_comm%sendjsR(1,0:rsize))
      allocate(d_comm%rem_addrl(list_size,0:rsize))
      d_comm%rem_addrl=-9999
      d_comm%cfrom_pe=-1
      d_comm%recvis=0; d_comm%recvie=0
      d_comm%recvjs=0; d_comm%recvje=0;
      d_comm%R_msize=0
      d_comm%R_do_buf=.false.
      d_comm%isizeR=0; d_comm%jsizeR=0
      d_comm%sendisR=0; d_comm%sendjsR=0

      call mpp_get_compute_domain( domain_out, isc, iec, jsc, jec )
      do list = 0,rsize
         m = mod( domain_in%pos+rsize+1-list, rsize+1 )
         d_comm%cfrom_pe(list) = domain_in%list(m)%pe
         from_pe = d_comm%cfrom_pe(list)
         call mpp_get_compute_domain( domain_in%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie >= is .AND. je >= js )then
           d_comm%R_do_buf(list) = .true.
           d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
           d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
           d_comm%R_msize(1,list) = (ie-is+1)*(je-js+1)*ke
#ifdef use_GSM
           if(use_shmem_ptr)then
             do i=1,list_size
               ptr = l_addrs_in(i)
               d_comm%rem_addrl(i,list) = shmem_ptr(f,from_pe)
             end do
           else
             call mpp_recv(d_comm%rem_addrl(1,list), glen=list_size, from_pe=from_pe )
           endif
           call mpp_recv(d_comm%sendisR(1,list), glen=1, from_pe=from_pe )
           call mpp_recv(d_comm%sendjsR(1,list), glen=1, from_pe=from_pe )
           call mpp_recv(d_comm%isizeR(list), glen=1, from_pe=from_pe )
           call mpp_recv(d_comm%jsizeR(list), glen=1, from_pe=from_pe )
#endif
         end if
      end do 

#ifdef use_CAF
    ! just need address of local input array for CAF
      d_comm%rem_addrl(:,1) = l_addrs_in(1:list_size)
#endif

      d_comm%isize_max = isize_in; call mpp_max(d_comm%isize_max)
      d_comm%jsize_max = jsize_in; call mpp_max(d_comm%jsize_max)

#if (!defined(use_GSM)) && (!defined(use_CAF))
    ! Handles case where S_msize and/or R_msize are 0 size array
      msgsize = ( MAXVAL( (/0,sum(d_comm%S_msize(:,:))/) ) + MAXVAL( (/0,sum(d_comm%R_msize(:,:))/) ) ) * list_size
      if(msgsize>0)then
         mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, msgsize )
         if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
             write( text,'(i8)' )mpp_domains_stack_hwm
             call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                  //trim(text)//') from all PEs.' )
         end if
      end if
#endif

      DEALLOCATE(slist_addr,isL,jsL)

      d_comm%initialized = .true.

    end function mpp_redistribute_init_comm


    function mpp_global_field_init_comm(domain,l_addr,isize_g,jsize_g,isize_l, &
                                        jsize_l, ksize,l_addr2,flags) RESULT(d_comm)
      type(DomainCommunicator2D), pointer       :: d_comm
      type(domain2D),target,      intent(in)    :: domain
      integer(LONG_KIND),         intent(in)    :: l_addr
      integer,                    intent(in)    :: isize_g
      integer,                    intent(in)    :: jsize_g
      integer,                    intent(in)    :: isize_l
      integer,                    intent(in)    :: jsize_l 
      integer,                    intent(in)    :: ksize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2
      integer, optional,          intent(in)    :: flags

      integer(LONG_KIND) :: domain_id
      integer :: n, lpos, rpos, list, nlist, tile_id
      integer :: update_flags
      logical :: xonly, yonly
      integer :: is, ie, js, je, ioff, joff, ishift, jshift
      integer :: lsize,rsize,msgsize,to_pe,from_pe
      integer,           allocatable,dimension(:) :: isL, jsL
      integer(LONG_KIND),allocatable,dimension(:,:) :: slist_addr
      integer(LONG_KIND),save       ,dimension(2)   :: rem_addr
      real,dimension(*) :: f, f2
      pointer(ptr,f);  pointer(ptr2,f2)
      real,dimension(*) :: f_r, f2_r
      pointer(ptr_r,f_r);  pointer(ptr2_r,f2_r)
      logical :: use_shmem_ptr
      character(len=8) :: text


      use_shmem_ptr=.false.  ! False for any but use_GSM
#ifdef use_GSM
      call check_sma_env()
      ptr=l_addr; if(PRESENT(l_addr2))ptr2=l_addr2  ! Set up target for shmem_ptr
      use_shmem_ptr=.false.  ! False for CAF ; ptr_r=shmem_ptr(f,mpp_pe()); if(ptr_r>0)use_shmem_ptr=.true.
#endif

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: must first call mpp_domains_init.' )
      update_flags=0; xonly = .FALSE.;  yonly = .FALSE.
      if( PRESENT(flags) )then
          update_flags = flags
          xonly = flags == XUPDATE
          yonly = flags == YUPDATE
          if( .NOT.xonly .AND. .NOT.yonly )call mpp_error( WARNING, 'MPP_GLOBAL_FIELD: you must have flags=XUPDATE or YUPDATE.' )
      end if

      !--- get shift size for global data, when domain is non-symmetry, ishift and jshift will be 0.
      call mpp_get_global_shift(domain, isize_l, jsize_l, ishift, jshift)

      if( isize_g /= (domain%x%global%size+ishift) .OR. jsize_g /= (domain%y%global%size+jshift) ) &
           call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_INIT_COMM: incoming arrays do not match domain.' )

      if( isize_l == (domain%x%compute%size+ishift) .AND. jsize_l == (domain%y%compute%size+jshift) )then
!local is on compute domain                                                                      
        ioff = -domain%x%compute%begin + 1                                                       
        joff = -domain%y%compute%begin + 1
      elseif( isize_l == (domain%x%data%size+ishift) .AND. jsize_l == (domain%y%data%size+jshift) )then
!local is on data domain                                                                        
        ioff = -domain%x%data%begin + 1                                                         
        joff = -domain%y%data%begin + 1
      else
        call mpp_error(FATAL,'MPP_GLOBAL_FIELD_INIT_COMM: incoming field array must match either compute domain or data domain.')
      endif

    ! Create unique domain identifier
      domain_id=set_domain_id(domain%id,ksize,update_flags)
      d_comm =>get_comm(domain_id,l_addr,l_addr2)
      
      if(d_comm%initialized)return   ! Found existing field/domain communicator

      d_comm%domain =>domain
      d_comm%isize_in = isize_l; d_comm%isize_out = isize_g
      d_comm%jsize_in = jsize_l; d_comm%jsize_out = jsize_g
      d_comm%ke = ksize
      d_comm%gf_ioff=ioff; d_comm%gf_joff=joff

!fill off-domains (note loops begin at an offset of 1)
      if( xonly )then
!send
        d_comm%Slist_size = size(domain%x%list(:))
        lsize = d_comm%Slist_size-1
        allocate(d_comm%cto_pe(0:lsize))
        d_comm%cto_pe=-1
        do list = 0,lsize
          lpos = mod(domain%x%pos+lsize+1-list,lsize+1)
          d_comm%cto_pe(list) = domain%x%list(lpos)%pe
        end do
!recv
        d_comm%Rlist_size = d_comm%Slist_size
        rsize = d_comm%Rlist_size-1
        allocate(d_comm%cfrom_pe(0:rsize))
        allocate(d_comm%recvis(1,0:rsize), d_comm%recvie(1,0:rsize), &
                 d_comm%recvjs(1,0:rsize), d_comm%recvje(1,0:rsize), &
                 d_comm%R_msize(1,0:rsize))                          
        d_comm%cfrom_pe=-1
        d_comm%recvis=0; d_comm%recvie=0
        d_comm%recvjs=0; d_comm%recvje=0;
        d_comm%R_msize=0
        do list = 0,rsize
          rpos = mod(domain%x%pos+list,rsize+1)
          from_pe = domain%x%list(rpos)%pe
          d_comm%cfrom_pe(list) = from_pe
          is = domain%list(from_pe)%x%compute%begin; ie = domain%list(from_pe)%x%compute%end + ishift
          js = domain%y%compute%begin; je = domain%y%compute%end + jshift
          d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
          d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
          d_comm%R_msize(1,list) = (ie-is+1) * (je-js+1) * ksize
        end do

      elseif( yonly )then
!send  
        d_comm%Slist_size = size(domain%y%list(:))
        lsize = d_comm%Slist_size-1
        allocate(d_comm%cto_pe(0:lsize))
        d_comm%cto_pe=-1
        do list = 0,lsize
          lpos = mod(domain%y%pos+lsize+1-list,lsize+1)
          d_comm%cto_pe(list) = domain%y%list(lpos)%pe
        end do
!recv    
        d_comm%Rlist_size = d_comm%Slist_size
        rsize = d_comm%Rlist_size-1
        allocate(d_comm%cfrom_pe(0:rsize))
        allocate(d_comm%recvis(1,0:rsize), d_comm%recvie(1,0:rsize), &
                 d_comm%recvjs(1,0:rsize), d_comm%recvje(1,0:rsize), &
                 d_comm%R_msize(1,0:rsize))                          
        d_comm%cfrom_pe=-1
        d_comm%recvis=0; d_comm%recvie=0
        d_comm%recvjs=0; d_comm%recvje=0;
        d_comm%R_msize=0
        do list = 0,rsize
          rpos = mod(domain%y%pos+list,rsize+1)
          from_pe = domain%y%list(rpos)%pe
          d_comm%cfrom_pe(list) = from_pe
          is = domain%x%compute%begin; ie = domain%x%compute%end + ishift
          js = domain%list(from_pe)%y%compute%begin; je = domain%list(from_pe)%y%compute%end + jshift
          d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
          d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
          d_comm%R_msize(1,list) = (ie-is+1) * (je-js+1) * ksize
        end do

      else

!send
        tile_id = domain%tile_id
        d_comm%Slist_size =  count(domain%list(:)%tile_id == tile_id)
        lsize = d_comm%Slist_size-1
        allocate(d_comm%cto_pe(0:lsize))
        d_comm%cto_pe=-1
        n = 0
        nlist = size(domain%list(:))        
        do list = 0,nlist-1
          lpos = mod(domain%pos+nlist-list,nlist)
          if( domain%list(lpos)%tile_id .NE. tile_id ) cycle
          d_comm%cto_pe(n) = domain%list(lpos)%pe
          n = n+1
        end do
!recv
        d_comm%Rlist_size = d_comm%Slist_size
        rsize = d_comm%Rlist_size-1
        allocate(d_comm%cfrom_pe(0:rsize))
        allocate(d_comm%recvis(1,0:rsize), d_comm%recvie(1,0:rsize), &
                 d_comm%recvjs(1,0:rsize), d_comm%recvje(1,0:rsize), &
                 d_comm%R_msize(1,0:rsize))
        d_comm%cfrom_pe=-1
        d_comm%recvis=0; d_comm%recvie=0
        d_comm%recvjs=0; d_comm%recvje=0;
        d_comm%R_msize=0
        n = 0
        do list = 0,nlist-1
          rpos = mod(domain%pos+list,nlist)
          if( domain%list(rpos)%tile_id .NE. tile_id ) cycle
          d_comm%cfrom_pe(n) = domain%list(rpos)%pe
          is = domain%list(rpos)%x%compute%begin; ie = domain%list(rpos)%x%compute%end + ishift
          js = domain%list(rpos)%y%compute%begin; je = domain%list(rpos)%y%compute%end + jshift
          d_comm%recvis(1,n)=is; d_comm%recvie(1,n)=ie
          d_comm%recvjs(1,n)=js; d_comm%recvje(1,n)=je
          d_comm%R_msize(1,n) = (je-js+1) * (ie-is+1) * ksize
          n = n+1
        end do

      endif

!send 

      allocate(d_comm%sendis(1,0:lsize), d_comm%sendie(1,0:lsize), &
               d_comm%sendjs(1,0:lsize), d_comm%sendje(1,0:lsize), &
               d_comm%S_msize(1,0:lsize),isL(0:lsize),jsL(0:lsize))
      allocate(slist_addr(2,0:lsize))
      isL=0; jsL=0
      slist_addr = -9999
      d_comm%sendis=0; d_comm%sendie=0
      d_comm%sendjs=0; d_comm%sendje=0;
      d_comm%S_msize=0
      do list = 0,lsize
        to_pe = d_comm%cto_pe(list)
        is=domain%x%compute%begin; ie=domain%x%compute%end + ishift
        js=domain%y%compute%begin; je=domain%y%compute%end + jshift
        d_comm%sendis(1,list)=is; d_comm%sendie(1,list)=ie
        d_comm%sendjs(1,list)=js; d_comm%sendje(1,list)=je
        d_comm%S_msize(1,list) = (je-js+1) * (ie-is+1) * ksize
        isL(list) = ioff+domain%x%compute%begin; jsL(list) = joff+domain%y%compute%begin
#ifdef use_GSM
        if(.not.use_shmem_ptr)then  ! address is from GLOBAL_ALLOC
          if(.not.PRESENT(l_addr2))then
            slist_addr(1,list) = d_comm%l_addr
            call mpp_send(slist_addr(1,list), plen=1, to_pe=to_pe )
          else
            slist_addr(1,list) = d_comm%l_addrx
            slist_addr(2,list) = d_comm%l_addry
            call mpp_send(slist_addr(1,list), plen=2, to_pe=to_pe )
          endif
        endif
        call mpp_send(isL(list), plen=1, to_pe=to_pe )
        call mpp_send(jsL(list), plen=1, to_pe=to_pe )
        call mpp_send(isize_l, plen=1, to_pe=to_pe )
        call mpp_send(jsize_l, plen=1, to_pe=to_pe )
#endif
      end do

      if( xonly ) then
        call mpp_sync_self(domain%x%list(:)%pe)
      elseif( yonly ) then
        call mpp_sync_self(domain%y%list(:)%pe)
      else
        call mpp_sync_self(domain%list(:)%pe)
      endif

!recv
      allocate(d_comm%isizeR(0:rsize), d_comm%jsizeR(0:rsize))
      allocate(d_comm%sendisR(1,0:rsize), d_comm%sendjsR(1,0:rsize))
      if(.not.PRESENT(l_addr2))then
        allocate(d_comm%rem_addr(0:rsize))
        d_comm%rem_addr=-9999
      else
        allocate(d_comm%rem_addrx(0:rsize),d_comm%rem_addry(0:rsize))
        d_comm%rem_addrx=-9999; d_comm%rem_addry=-9999
      endif
      d_comm%isizeR=0; d_comm%jsizeR=0
      d_comm%sendisR=0; d_comm%sendjsR=0
      rem_addr = -9999
      do list = 0,rsize
        from_pe = d_comm%cfrom_pe(list)
#ifdef use_GSM
        if(.not.PRESENT(l_addr2))then
          if(use_shmem_ptr)then
            d_comm%rem_addr(list) = shmem_ptr(f,from_pe)
          else
            call mpp_recv(rem_addr(1), glen=1, from_pe=from_pe )
            d_comm%rem_addr(list) = rem_addr(1)
          endif
        else   
          if(use_shmem_ptr)then
            d_comm%rem_addrx(list) = shmem_ptr(f,from_pe)
            d_comm%rem_addry(list) = shmem_ptr(f2,from_pe)
          else
            call mpp_recv(rem_addr(1), glen=2, from_pe=from_pe )
            d_comm%rem_addrx(list) = rem_addr(1)
            d_comm%rem_addry(list) = rem_addr(2)
          endif
        endif  
        call mpp_recv(d_comm%sendisR(1,list), glen=1, from_pe=from_pe )
        call mpp_recv(d_comm%sendjsR(1,list), glen=1, from_pe=from_pe )
        call mpp_recv(d_comm%isizeR(list), glen=1, from_pe=from_pe )
        call mpp_recv(d_comm%jsizeR(list), glen=1, from_pe=from_pe )
#endif
      end do

#ifdef use_CAF
      if(.not.PRESENT(l_addr2))then
        d_comm%rem_addr(0) = d_comm%l_addr
      else
        d_comm%rem_addrx(0) = d_comm%l_addrx
        d_comm%rem_addry(0) = d_comm%l_addry
      endif
#endif

#ifdef use_GSM
      d_comm%isize_max = isize_l; call mpp_max(d_comm%isize_max)
      d_comm%jsize_max = jsize_l; call mpp_max(d_comm%jsize_max)
#endif

#if (!defined(use_GSM)) && (!defined(use_CAF))
    ! Handles case where S_msize and/or R_msize are 0 size array
      msgsize = MAXVAL( (/0,sum(d_comm%S_msize(:,:))/) ) + MAXVAL( (/0,sum(d_comm%R_msize(:,:))/) )
      if(msgsize>0)then
         mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, msgsize )
         if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
             write( text,'(i8)' )mpp_domains_stack_hwm
             call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_INIT_COMM: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                  //trim(text)//') from all PEs.' )
         end if
      end if
#endif
            
      DEALLOCATE(slist_addr,isL,jsL)

      d_comm%initialized = .true.

    end function mpp_global_field_init_comm


    subroutine mpp_update_free_comm(domain,l_addr,ksize,lsize,l_addr2,flags,gridtype)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      type(domain2D),target,      intent(in)    :: domain
      integer(LONG_KIND),         intent(in)    :: l_addr
      integer,                    intent(in)    :: ksize
      integer,                    intent(in)    :: lsize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2
      integer, optional,          intent(in)    :: flags
      integer, optional,          intent(in)    :: gridtype

      integer(LONG_KIND) :: domain_id
      integer :: update_flags

      update_flags = XUPDATE+YUPDATE; if( PRESENT(flags) )update_flags = flags
      domain_id=set_domain_id(domain%id,ksize+lsize,update_flags,gridtype)
      call free_comm(domain_id,l_addr,l_addr2)
    end subroutine mpp_update_free_comm


    subroutine mpp_redistribute_free_comm(domain_in,l_addr,domain_out,l_addr2,ksize,lsize)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      type(domain2D),             intent(in)    :: domain_in
      integer(LONG_KIND),         intent(in)    :: l_addr
      type(domain2D),             intent(in)    :: domain_out  
      integer(LONG_KIND),         intent(in)    :: l_addr2
      integer,                    intent(in)    :: ksize,lsize

      integer(LONG_KIND) :: domain_id

      if(l_addr2 > 0)then
        domain_id = set_domain_id(domain_out%id,ksize+lsize)
      else
        domain_id = set_domain_id(domain_in%id,ksize+lsize)
      endif
      call free_comm(domain_id,l_addr,l_addr2)
    end subroutine mpp_redistribute_free_comm


    subroutine mpp_global_field_free_comm(domain,l_addr,ksize,l_addr2,flags)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      type(domain2D),             intent(in)    :: domain
      integer(LONG_KIND),         intent(in)    :: l_addr
      integer,                    intent(in)    :: ksize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2
      integer,           optional,intent(in)    :: flags

      integer :: update_flags
      integer(LONG_KIND) :: domain_id

      update_flags=0; if(PRESENT(flags))update_flags=flags
      domain_id=set_domain_id(domain%id,ksize,update_flags)
      call free_comm(domain_id,l_addr,l_addr2)
    end subroutine mpp_global_field_free_comm


    subroutine free_comm(domain_id,l_addr,l_addr2)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      integer(LONG_KIND),         intent(in) :: domain_id
      integer(LONG_KIND),         intent(in) :: l_addr
      integer(LONG_KIND),optional,intent(in) :: l_addr2

      integer(LONG_KIND) :: dc_key,a_key
      integer :: dc_idx,a_idx,i_idx,insert,insert_a,insert_i
      integer :: a2_idx,insert_a2


      i_idx = find_key(domain_id,ids_sorted(1:n_ids),insert_i)
      a_idx = find_key(l_addr,addrs_sorted(1:a_sort_len),insert_a)
      a_key = int(addrs_idx(a_idx),KIND(LONG_KIND))
      if(PRESENT(l_addr2))then
        a2_idx = find_key(l_addr2,addrs2_sorted(1:a2_sort_len),insert_a2)
        a_key = a_key + ADDR2_BASE*int(addrs2_idx(a2_idx),KIND(LONG_KIND))
      endif
      dc_key = DOMAIN_ID_BASE*int(ids_idx(i_idx),KIND(LONG_KIND)) + a_key
      dc_idx = find_key(dc_key,dcKey_sorted(1:dc_sort_len),insert)

      if(dc_idx < 0)then
        call mpp_error(FATAL,'FREE_COMM: attempt to remove nonexistent domains communicator key')
      endif
      call deallocate_comm(d_comm(dc_idx))
      call pop_key(dcKey_sorted,d_comm_idx,dc_sort_len,dc_idx)
      call pop_key(addrs_sorted,addrs_idx,a_sort_len,a_idx)
      if(PRESENT(l_addr2))call pop_key(addrs2_sorted,addrs2_idx,a2_sort_len,a2_idx)
    end subroutine free_comm


    function get_comm(domain_id,l_addr,l_addr2)
      integer(LONG_KIND),intent(in)       :: domain_id
      integer(LONG_KIND),intent(in)       :: l_addr
      integer(LONG_KIND),intent(in),optional :: l_addr2
      type(DomainCommunicator2D), pointer :: get_comm

      integer(LONG_KIND) :: dc_key,a_key
      integer :: i,dc_idx,a_idx,i_idx,insert,insert_a,insert_i
      integer :: a2_idx,insert_a2

      if(.not.ALLOCATED(d_comm))ALLOCATE(d_comm(MAX_FIELDS))
      i_idx   = find_key(domain_id,ids_sorted(1:n_ids),insert_i)
      a_idx = find_key(l_addr,addrs_sorted(1:a_sort_len),insert_a)
      a_key = int(addrs_idx(a_idx),KIND(LONG_KIND))
      if(PRESENT(l_addr2))then
        a2_idx = find_key(l_addr2,addrs2_sorted(1:a2_sort_len),insert_a2)
        a_key = a_key + ADDR2_BASE*int(addrs2_idx(a2_idx),KIND(LONG_KIND))
      endif
      dc_key   = DOMAIN_ID_BASE*int(ids_idx(i_idx),KIND(LONG_KIND)) + a_key
      dc_idx   = find_key(dc_key,dcKey_sorted(1:dc_sort_len),insert)
      if(dc_idx > 0)then
        get_comm =>d_comm(d_comm_idx(dc_idx))
      else
        if(i_idx<0)then
          if(n_ids == MAX_DOM_IDS)then
            call mpp_error(FATAL,'GET_COMM: Maximum number of domains exceeded')
          endif
          n_ids = n_ids+1
          i_idx = push_key(ids_sorted,ids_idx,i_sort_len,insert_i,domain_id,n_ids)
        endif
        if(a_idx<0)then
          if(n_addrs == MAX_ADDRS)then
             call mpp_error(FATAL,'GET_COMM: Maximum number of memory addresses exceeded')
          endif
          n_addrs = n_addrs + 1
          a_idx = push_key(addrs_sorted,addrs_idx,a_sort_len,insert_a,l_addr,n_addrs)
        endif
        if(PRESENT(l_addr2))then
          if(a2_idx<0)then
            if(n_addrs2 == MAX_ADDRS2)then
               call mpp_error(FATAL,'GET_COMM: Maximum number of 2nd memory addresses exceeded')
            endif
            n_addrs2 = n_addrs2 + 1
            a2_idx = push_key(addrs2_sorted,addrs2_idx,a2_sort_len,insert_a2,l_addr2,n_addrs2)
          endif
        endif
        if(n_comm == MAX_FIELDS)then
          call mpp_error(FATAL,'GET_COMM: Maximum number of fields exceeded')
        endif
        a_key = int(addrs_idx(a_idx),KIND(8))
        if(PRESENT(l_addr2))a_key = a_key + ADDR2_BASE*int(addrs2_idx(a2_idx),KIND(8))
        dc_key   = DOMAIN_ID_BASE*int(ids_idx(i_idx),KIND(LONG_KIND)) + a_key
        dc_idx   = find_key(dc_key,dcKey_sorted(1:dc_sort_len),insert)
        if(dc_idx /= -1)call mpp_error(FATAL,'GET_COMM: attempt to insert existing key')
        n_comm = n_comm + 1
        i = push_key(dcKey_sorted,d_comm_idx,dc_sort_len,insert,dc_key,n_comm)
        d_comm_idx(insert) = n_comm
        if(PRESENT(l_addr2))then
          d_comm(n_comm)%l_addrx = l_addr
          d_comm(n_comm)%l_addry = l_addr2
        else
          d_comm(n_comm)%l_addr = l_addr
        endif
        get_comm =>d_comm(n_comm)
      endif
    end function get_comm


    function push_key(sorted,idx,n_idx,insert,key,ival)
      integer(LONG_KIND),intent(inout),dimension(:)   :: sorted
      integer,           intent(inout),dimension(-1:) :: idx  ! Start -1 to simplify first call logic in get_comm
      integer,           intent(inout)                :: n_idx
      integer,           intent(in)                   :: insert
      integer(LONG_KIND),intent(in)                   :: key
      integer,           intent(in)                   :: ival

      integer                                         :: push_key,i

      do i=n_idx,insert,-1
        sorted(i+1) = sorted(i)
        idx(i+1)    = idx(i)
      end do
      sorted(insert) = key
      n_idx = n_idx + 1
      idx(insert) = ival
      push_key = insert
    end function push_key


    subroutine pop_key(sorted,idx,n_idx,key_idx)
      integer(LONG_KIND),intent(inout),dimension(:)   :: sorted
      integer,           intent(inout),dimension(-1:) :: idx  ! Start -1 to simplify first call logic in get_comm
      integer,           intent(inout)                :: n_idx
      integer,           intent(in)                   :: key_idx

      integer                                         :: i

      do i=key_idx,n_idx-1
        sorted(i) = sorted(i+1)
        idx(i)    = idx(i+1)  
      end do                  
      sorted(n_idx) = -9999 
      idx(n_idx) = -9999
      n_idx = n_idx - 1  
    end subroutine pop_key


    function find_key(key,sorted,insert) RESULT(n)
    ! The algorithm used here requires monotonic keys w/out repetition.
      integer(LONG_KIND),intent(in)              :: key        ! new address to be found in list
      integer(LONG_KIND),dimension(:),intent(in) :: sorted  ! list of sorted local addrs
      integer,                        intent(out) :: insert
      integer :: n, n_max, n_min, n_key
      logical :: not_found

      n_key = size(sorted(:))
      insert = 1
      n = -1  ! value not in list
      if(n_key == 0)return  ! first call
      
      if(key < sorted(1))then
        insert = 1; return
      elseif(key > sorted(n_key))then
        insert = n_key+1; return
      endif
        
      if(key == sorted(1))then
        n = 1; return
      elseif(key == sorted(n_key))then
        n = n_key; return
      endif
        
      not_found = .true.
      n = n_key/2 + 1
      n_min=1; n_max=n_key
      do while(not_found)
        if(key == sorted(n))then
          not_found = .false.
        elseif(key > sorted(n))then
          if(key < sorted(n+1))then
            insert = n+1; exit
          endif
          n_min = n
          n = (n+1+n_max)/2
        else
          if(key > sorted(n-1))then
            insert = n; exit
          endif
          n_max = n
          n = (n+n_min)/2
        endif
        if(n==1 .or. n==n_key)exit
      end do
      if(not_found)n = -1  ! value not in list
    end function find_key


    subroutine deallocate_comm(d_comm)
      type(DomainCommunicator2D), intent(inout) :: d_comm

      d_comm%domain  =>NULL()
      d_comm%domain_in  =>NULL()
      d_comm%domain_out =>NULL()

      d_comm%initialized=.false.
      d_comm%id=-9999
      d_comm%l_addr  =-9999
      d_comm%l_addrx =-9999
      d_comm%l_addry =-9999

      if( _ALLOCATED(d_comm%sendis) )   DEALLOCATE(d_comm%sendis);    !!d_comm%sendis =>NULL()
      if( _ALLOCATED(d_comm%sendie) )   DEALLOCATE(d_comm%sendie);    !!d_comm%sendie =>NULL()
      if( _ALLOCATED(d_comm%sendjs) )   DEALLOCATE(d_comm%sendjs);    !!d_comm%sendjs =>NULL()
      if( _ALLOCATED(d_comm%sendje) )   DEALLOCATE(d_comm%sendje);    !!d_comm%sendje =>NULL()
      if( _ALLOCATED(d_comm%S_msize) )  DEALLOCATE(d_comm%S_msize);   !!d_comm%S_msize =>NULL()
      if( _ALLOCATED(d_comm%do_thisS) ) DEALLOCATE(d_comm%do_thisS);  !!d_comm%do_thisS =>NULL()
      if( _ALLOCATED(d_comm%S_do_buf) ) DEALLOCATE(d_comm%S_do_buf);  !!d_comm%S_do_buf =>NULL()
      if( _ALLOCATED(d_comm%cto_pe) )   DEALLOCATE(d_comm%cto_pe);    !!d_comm%cto_pe  =>NULL()
      if( _ALLOCATED(d_comm%Rcaf_idx) ) DEALLOCATE(d_comm%Rcaf_idx);  !!d_comm%Rcaf_idx  =>NULL()
      if( _ALLOCATED(d_comm%recvis) )   DEALLOCATE(d_comm%recvis);    !!d_comm%recvis =>NULL()
      if( _ALLOCATED(d_comm%recvie) )   DEALLOCATE(d_comm%recvie);    !!d_comm%recvie =>NULL()
      if( _ALLOCATED(d_comm%recvjs) )   DEALLOCATE(d_comm%recvjs);    !!d_comm%recvjs =>NULL()
      if( _ALLOCATED(d_comm%recvje) )   DEALLOCATE(d_comm%recvje);    !!d_comm%recvje =>NULL()
      if( _ALLOCATED(d_comm%R_msize) )  DEALLOCATE(d_comm%R_msize);   !!d_comm%R_msize =>NULL()
      if( _ALLOCATED(d_comm%do_thisR) ) DEALLOCATE(d_comm%do_thisR);  !!d_comm%do_thisR =>NULL()
      if( _ALLOCATED(d_comm%R_do_buf) ) DEALLOCATE(d_comm%R_do_buf);  !!d_comm%R_do_buf =>NULL()
      if( _ALLOCATED(d_comm%cfrom_pe) ) DEALLOCATE(d_comm%cfrom_pe);  !!d_comm%cfrom_pe  =>NULL()
      d_comm%Slist_size=0; d_comm%Rlist_size=0
      d_comm%isize=0; d_comm%jsize=0; d_comm%ke=0
      d_comm%isize_in=0; d_comm%jsize_in=0
      d_comm%isize_out=0; d_comm%jsize_out=0
      d_comm%isize_max=0; d_comm%jsize_max=0
      d_comm%gf_ioff=0; d_comm%gf_joff=0
    ! Remote data
      if( _ALLOCATED(d_comm%isizeR) )   DEALLOCATE(d_comm%isizeR);    !!dd_comm%isizeR =>NULL()
      if( _ALLOCATED(d_comm%jsizeR) )   DEALLOCATE(d_comm%jsizeR);    !!dd_comm%jsizeR =>NULL()
      if( _ALLOCATED(d_comm%sendisR) )  DEALLOCATE(d_comm%sendisR);   !!dd_comm%sendisR =>NULL()
      if( _ALLOCATED(d_comm%sendjsR) )  DEALLOCATE(d_comm%sendjsR);   !!dd_comm%sendjsR =>NULL()
      if( _ALLOCATED(d_comm%rem_addr) ) DEALLOCATE(d_comm%rem_addr);  !!dd_comm%rem_addr  =>NULL()
      if( _ALLOCATED(d_comm%rem_addrx) )DEALLOCATE(d_comm%rem_addrx); !!dd_comm%rem_addrx =>NULL()
      if( _ALLOCATED(d_comm%rem_addry) )DEALLOCATE(d_comm%rem_addry); !!dd_comm%rem_addry =>NULL()
      if( _ALLOCATED(d_comm%rem_addrl) )DEALLOCATE(d_comm%rem_addrl); !!dd_comm%rem_addrl =>NULL()
    end subroutine deallocate_comm


    function set_domain_id(d_id,ksize,flags,gtype, position)
      integer(LONG_KIND), intent(in) :: d_id
      integer           , intent(in) :: ksize
      integer           , optional, intent(in) :: flags
      integer           , optional, intent(in) :: gtype
      integer           , optional, intent(in) :: position

      integer(LONG_KIND)             :: set_domain_id
      
      set_domain_id=d_id + KE_BASE*int(ksize,KIND(d_id))
      if(PRESENT(flags))set_domain_id=set_domain_id+int(flags,KIND(d_id))
      if(PRESENT(gtype))set_domain_id=set_domain_id+GT_BASE*int(gtype,KIND(d_id))  ! Must be LONG_KIND arithmetic
      !--- gtype is never been used to set id. we need to add position to calculate id to seperate
      !--- BGRID and CGRID or scalar variable.
      if(present(position)) then
         set_domain_id=set_domain_id+GT_BASE*int(2**position, KIND(d_id))
      endif

    end function set_domain_id


#ifdef use_GSM
    subroutine check_sma_env()
      character(len=15) :: sma_env
      call getenv('SMA_GLOBAL_HEAP_SIZE',sma_env)
      if(sma_env == '               ') &
         call mpp_error( FATAL, 'Environment variables SMA_GLOBAL_ALLOC and SMA_GLOBAL_HEAP_SIZE must be set.' )
    end subroutine check_sma_env
#endif

!#######################################################################


#ifdef use_CAF
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_r8_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_r8_1d
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_r8_3d_type
#undef MPP_ASSOC_CAF_FIELD_3D_
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_r8_3d
#include <mpp_domains_comm.h>

#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_c8_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_c8_1d
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_c8_3d_type
#undef MPP_ASSOC_CAF_FIELD_3D_
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_c8_3d
#include <mpp_domains_comm.h>

#ifndef no_8byte_integers
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_i8_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_i8_1d
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_i8_3d_type
#undef MPP_ASSOC_CAF_FIELD_3D_
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_i8_3d
#include <mpp_domains_comm.h>

#undef MPP_TYPE_
#define MPP_TYPE_ logical(LONG_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_l8_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_l8_1d
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_l8_3d_type
#undef MPP_ASSOC_CAF_FIELD_3D_
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_l8_3d
#include <mpp_domains_comm.h>
#endif

#ifndef no_4byte_reals
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_r4_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_r4_1d
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_r4_3d_type
#undef MPP_ASSOC_CAF_FIELD_3D_
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_r4_3d
#include <mpp_domains_comm.h>

#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_c4_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_c4_1d
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_c4_3d_type
#undef MPP_ASSOC_CAF_FIELD_3D_
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_c4_3d
#include <mpp_domains_comm.h>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_i4_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_i4_1d
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_i4_3d_type
#undef MPP_ASSOC_CAF_FIELD_3D_
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_i4_3d
#include <mpp_domains_comm.h>

#undef MPP_TYPE_
#define MPP_TYPE_ logical(INT_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_l4_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_l4_1d
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_l4_3d_type
#undef MPP_ASSOC_CAF_FIELD_3D_
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_l4_3d
#include <mpp_domains_comm.h>                                                                      
#endif

