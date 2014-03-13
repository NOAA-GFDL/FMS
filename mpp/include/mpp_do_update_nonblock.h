! -*-f90-*- 
subroutine MPP_START_DO_UPDATE_3D_(id_update, f_addrs, domain, update, d_type, ke_max, ke_list, flags, reuse_id_update, name)
  integer,                    intent(in) :: id_update
  integer(LONG_KIND),         intent(in) :: f_addrs(:,:)
  type(domain2D),             intent(in) :: domain
  type(overlapSpec),          intent(in) :: update
  MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
  integer,                    intent(in) :: ke_max
  integer,                    intent(in) :: ke_list(:,:)
  logical,                    intent(in) :: reuse_id_update
  character(len=*),           intent(in) :: name
  integer,                    intent(in) :: flags
  
  !--- local variables
  integer                     :: i, j, k, m, n, l, dir, tMe
  integer                     :: buffer_pos, msgsize, from_pe, to_pe, pos
  integer                     :: is, ie, js, je, sendsize, recvsize
  logical                     :: send(8), recv(8), update_edge_only
  integer                     :: l_size, ke_sum, my_id_update
  integer                     :: request
  integer                     :: send_msgsize(MAXLIST)
  character(len=128)          :: text
  MPP_TYPE_                   :: buffer(size(mpp_domains_stack_nonblock(:)))
  MPP_TYPE_                   :: field(update%xbegin:update%xend, update%ybegin:update%yend,ke_max)
  pointer( ptr, buffer )
  pointer(ptr_field, field)

  update_edge_only = BTEST(flags, EDGEONLY)
  recv(1) = BTEST(flags,EAST)
  recv(3) = BTEST(flags,SOUTH)
  recv(5) = BTEST(flags,WEST)
  recv(7) = BTEST(flags,NORTH)
  if( update_edge_only ) then
     if( .NOT. (recv(1) .OR. recv(3) .OR. recv(5) .OR. recv(7)) ) then
        recv(1) = .true.
        recv(3) = .true.
        recv(5) = .true.
        recv(7) = .true.
     endif
  else
     recv(2) = recv(1) .AND. recv(3)
     recv(4) = recv(3) .AND. recv(5)
     recv(6) = recv(5) .AND. recv(7)
     recv(8) = recv(7) .AND. recv(1)
  endif
  send    = recv

  l_size = size(f_addrs,1)
  ke_sum = sum(ke_list)
  ptr = LOC(mpp_domains_stack_nonblock)

  buffer_pos = nonblock_data(id_update)%recv_pos

  if( update%nrecv > MAX_REQUEST ) then
     write( text,'(a,i8,a,i8)' ) 'update%nrecv =', update%nrecv, ' greater than MAX_REQEUST =', MAX_REQUEST
     call mpp_error(FATAL,'MPP_START_DO_UPDATE: '//trim(text))
  endif
  if( update%nsend > MAX_REQUEST ) then
     write( text,'(a,i8,a,i8)' ) 'update%nsend =', update%nsend, ' greater than MAX_REQEUST =', MAX_REQUEST
     call mpp_error(FATAL,'MPP_START_DO_UPDATE: '//trim(text))
  endif

  ! pre-postrecv
  !--- make sure the domain stack size is big enough.
  recvsize = 0
  do m = 1, update%nrecv
     nonblock_data(id_update)%size_recv(m) = 0
     if( update%recv(m)%count == 0 )cycle
     msgsize = 0
     do n = 1, update%recv(m)%count
        dir = update%recv(m)%dir(n)
        if(recv(dir)) then
           msgsize = msgsize + update%recv(m)%msgsize(n)           
        end if
     end do
     if( msgsize.GT.0 )then
        msgsize = msgsize*ke_sum
        recvsize = recvsize + msgsize
        nonblock_data(id_update)%size_recv(m) = msgsize
        nonblock_data(id_update)%buffer_pos_recv(m) = buffer_pos
        buffer_pos = buffer_pos + msgsize
     end if
  end do
     
  sendsize = 0
  do m = 1, update%nsend
     if( update%send(m)%count == 0 )cycle

     ! make sure the stacksize is big enough
     msgsize = 0
     do n = 1, update%send(m)%count
        dir = update%send(m)%dir(n)
        if( send(dir) )  msgsize = msgsize + update%send(m)%msgsize(n)
     enddo
     if( msgsize.GT.0 )then
        msgsize = msgsize*ke_sum
        sendsize = sendsize + msgsize
        nonblock_data(id_update)%buffer_pos_send(m) = buffer_pos
        buffer_pos = buffer_pos + msgsize        
     end if
  end do

  mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, &
      nonblock_data(id_update)%recv_pos+recvsize+sendsize )
  if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
     write( text,'(i8)' )mpp_domains_stack_hwm
     call mpp_error( FATAL, 'MPP_START_DO_UPDATE: mpp_domains_stack overflow, ' // &
          'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
  end if

  if( reuse_id_update ) then
     if(recvsize .NE. nonblock_data(id_update)%recv_msgsize) then
        call mpp_error(FATAL,'MPP_START_DO_UPDATE: mismatch of recv msgsize for field '//trim(name) )
     endif
     if(sendsize .NE. nonblock_data(id_update)%send_msgsize) then
        call mpp_error(FATAL,'MPP_START_DO_UPDATE: mismatch of send msgsize for field '//trim(name) )
     endif
  else
     nonblock_data(id_update)%recv_msgsize = recvsize
     nonblock_data(id_update)%send_msgsize = sendsize
     nonblock_data(id_update)%send_pos = nonblock_data(id_update)%recv_pos + recvsize
     nonblock_buffer_pos = nonblock_buffer_pos + recvsize + sendsize
  endif

  ! pre-postrecv
  call mpp_clock_begin(recv_clock_nonblock)
  do m = 1, update%nrecv
     msgsize = nonblock_data(id_update)%size_recv(m)
     if( msgsize.GT.0 )then
        from_pe =  update%recv(m)%pe
        buffer_pos = nonblock_data(id_update)%buffer_pos_recv(m)
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., &
             tag=id_update, request=request)
        nonblock_data(id_update)%request_recv(m) = request

#ifdef use_libMPI
        nonblock_data(id_update)%type_recv(m) = MPI_TYPE_
#endif
     end if
  end do ! end do m = 1, update%nrecv  

  call mpp_clock_end(recv_clock_nonblock)

  ! send
  call mpp_clock_begin(send_pack_clock_nonblock)
!$OMP parallel do schedule(dynamic) default(shared) private(buffer_pos,pos,dir,tMe,is,ie,js,je,ptr_field,to_pe, &
!$OMP msgsize,request)
  do m = 1, update%nsend
     send_msgsize(m) = 0
     if( update%send(m)%count == 0 )cycle
     
     buffer_pos = nonblock_data(id_update)%buffer_pos_send(m)
     pos = buffer_pos

     do n = 1, update%send(m)%count
        dir = update%send(m)%dir(n)
        if( send(dir) ) then
           tMe = update%send(m)%tileMe(n)
           is = update%send(m)%is(n); ie = update%send(m)%ie(n)
           js = update%send(m)%js(n); je = update%send(m)%je(n)
           select case( update%send(m)%rotation(n) )
           case(ZERO)
              do l=1,l_size  ! loop over number of fields
                 ptr_field = f_addrs(l, tMe)
                 do k = 1,ke_list(l,tMe)  
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              enddo
           case( MINUS_NINETY ) 
              do l=1,l_size  ! loop over number of fields
                 ptr_field = f_addrs(l, tMe)
                 do k = 1,ke_list(l,tMe)  
                    do i = is, ie
                       do j = je, js, -1
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case( NINETY ) 
              do l=1,l_size  ! loop over number of fields
                 ptr_field = f_addrs(l, tMe)

                 do k = 1,ke_list(l,tMe)  
                    do i = ie, is, -1
                       do j = js, je
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case( ONE_HUNDRED_EIGHTY ) 
              do l=1,l_size  ! loop over number of fields
                 ptr_field = f_addrs(l, tMe)
                 do k = 1,ke_list(l,tMe)  
                    do j = je, js, -1
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           end select
       endif
    end do ! do n = 1, update%send(m)%count
    send_msgsize(m) = pos - buffer_pos
 enddo
 !$OMP end parallel do

  do m = 1, update%nsend
    msgsize = send_msgsize(m)
     if( msgsize .GT.0 )then
        buffer_pos = nonblock_data(id_update)%buffer_pos_send(m)
        to_pe = update%send(m)%pe
        call mpp_send( buffer(buffer_pos+1), plen=msgsize , to_pe=to_pe, &
                       tag=id_update, request=request)
        nonblock_data(id_update)%request_send(m) = request
     end if
  end do ! end do ist = 0,nlist-1

  call mpp_clock_end(send_pack_clock_nonblock)

  return


end subroutine MPP_START_DO_UPDATE_3D_

!###############################################################################

subroutine MPP_COMPLETE_DO_UPDATE_3D_(id_update, f_addrs, domain, update, d_type, ke_max, ke_list, flags) 
  integer,             intent(in) :: id_update
  integer(LONG_KIND),  intent(in) :: f_addrs(:,:)
  type(domain2d),      intent(in) :: domain
  type(overlapSpec),   intent(in) :: update
  integer,             intent(in) :: ke_max
  integer,             intent(in) :: ke_list(:,:)
  MPP_TYPE_,           intent(in) :: d_type  ! creates unique interface
  integer,             intent(in) :: flags

  !--- local variables
  integer                     :: i, j, k, m, n, l, dir, count, tMe, tNbr
  integer                     :: buffer_pos, msgsize, from_pe, pos
  integer                     :: is, ie, js, je
  logical                     :: send(8), recv(8), update_edge_only
  integer                     :: l_size, ke_sum, sendsize, recvsize
  character(len=128)          :: text
  MPP_TYPE_                   :: recv_buffer(size(mpp_domains_stack_nonblock(:)))
  MPP_TYPE_                   :: field(update%xbegin:update%xend, update%ybegin:update%yend,ke_max)
  pointer( ptr, recv_buffer )
  pointer(ptr_field, field)

  update_edge_only = BTEST(flags, EDGEONLY)
  recv(1) = BTEST(flags,EAST)
  recv(3) = BTEST(flags,SOUTH)
  recv(5) = BTEST(flags,WEST)
  recv(7) = BTEST(flags,NORTH)
  if( update_edge_only ) then
     if( .NOT. (recv(1) .OR. recv(3) .OR. recv(5) .OR. recv(7)) ) then
        recv(1) = .true.
        recv(3) = .true.
        recv(5) = .true.
        recv(7) = .true.
     endif
  else
     recv(2) = recv(1) .AND. recv(3)
     recv(4) = recv(3) .AND. recv(5)
     recv(6) = recv(5) .AND. recv(7)
     recv(8) = recv(7) .AND. recv(1)
  endif
  send    = recv

  ke_sum = sum(ke_list)
  l_size = size(f_addrs,1)
  ptr = LOC(mpp_domains_stack_nonblock)

  count = update%nrecv
  if(count > 0) then
     call mpp_clock_begin(wait_clock_nonblock)
     call mpp_sync_self(check=EVENT_RECV, request=nonblock_data(id_update)%request_recv(1:count), &
                        msg_size=nonblock_data(id_update)%size_recv(1:count),                     &
                        msg_type=nonblock_data(id_update)%type_recv(1:count) )
     call mpp_clock_end(wait_clock_nonblock)
#ifdef use_libMPI
     nonblock_data(id_update)%request_recv(:)    = MPI_REQUEST_NULL
#else
     nonblock_data(id_update)%request_recv(:)    = 0
#endif
     nonblock_data(id_update)%type_recv(:) = 0
  endif 

  !--unpack the data
  call mpp_clock_begin(unpk_clock_nonblock)
!$OMP parallel do schedule(dynamic) default(shared) private(dir,buffer_pos,pos,tMe,is,ie,js,je,msgsize, &
!$OMP          ptr_field)
  do m = update%nrecv, 1, -1
     if( update%recv(m)%count == 0 )cycle
     buffer_pos = nonblock_data(id_update)%buffer_pos_recv(m) + nonblock_data(id_update)%size_recv(m)

     pos = buffer_pos
     do n = update%recv(m)%count, 1, -1
        dir = update%recv(m)%dir(n)
        if( recv(dir) ) then
           tMe = update%recv(m)%tileMe(n)
           is = update%recv(m)%is(n); ie = update%recv(m)%ie(n)
           js = update%recv(m)%js(n); je = update%recv(m)%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke_sum
           pos = buffer_pos - msgsize
           buffer_pos = pos
           do l=1, l_size  ! loop over number of fields
              ptr_field = f_addrs(l, tMe)
              do k = 1,ke_list(l,tMe)
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       field(i,j,k) = recv_buffer(pos)
                    end do
                 end do
              end do
           end do
        end if
     end do ! do n = 1, update%recv(m)%count
  end do
!$OMP end parallel do
  call mpp_clock_end(unpk_clock_nonblock)

  count = update%nrecv
  if(count > 0) then
     nonblock_data(id_update)%size_recv(:) = 0
  endif

  count = update%nsend
  if(count > 0) then
     call mpp_clock_begin(wait_clock_nonblock)
     call mpp_sync_self(check=EVENT_SEND, request=nonblock_data(id_update)%request_send(1:count))
     call mpp_clock_end(wait_clock_nonblock)
     nonblock_data(id_update)%request_send_count = 0
#ifdef use_libMPI
     nonblock_data(id_update)%request_send(:)    = MPI_REQUEST_NULL
#else
     nonblock_data(id_update)%request_send(:)    = 0
#endif
  endif 
  
!  call init_nonblock_type(nonblock_data(id_update))

  return

end subroutine MPP_COMPLETE_DO_UPDATE_3D_
