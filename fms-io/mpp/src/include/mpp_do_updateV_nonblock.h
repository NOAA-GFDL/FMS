! -*-f90-*- 
!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
subroutine MPP_START_DO_UPDATE_3D_V_(id_update, f_addrsx, f_addrsy, domain, update_x, update_y,     &
                                     d_type, ke_max, ke_list, gridtype, flags, reuse_id_update, name)
  integer,             intent(in) :: id_update
  integer(LONG_KIND),  intent(in) :: f_addrsx(:,:), f_addrsy(:,:)
  type(domain2d),      intent(in) :: domain
  type(overlapSpec),   intent(in) :: update_x, update_y
  integer,             intent(in) :: ke_max
  integer,             intent(in) :: ke_list(:,:)
  MPP_TYPE_,           intent(in) :: d_type  ! creates unique interface
  integer,             intent(in) :: gridtype
  logical,             intent(in) :: reuse_id_update
  character(len=*),    intent(in) :: name
  integer,             intent(in) :: flags

  !---local variable ------------------------------------------
  integer            :: i, j, k, l, is, ie, js, je, n, m
  integer            :: pos, nlist, msgsize, tile, l_size
  integer            :: to_pe, from_pe, buffer_pos
  integer            :: tMe, dir, ke_sum
  logical            :: send(8), recv(8), update_edge_only
  character(len=128) :: text
  integer            :: ind_x, ind_y
  integer            :: nsend, nrecv, sendsize, recvsize
  integer            :: request
  integer            :: send_msgsize(update_x%nsend+update_y%nsend)
  integer            :: ind_send_x(update_x%nsend+update_y%nsend), ind_send_y(update_x%nsend+update_y%nsend)
  integer            :: ind_recv_x(update_x%nrecv+update_y%nrecv), ind_recv_y(update_x%nrecv+update_y%nrecv)
  integer            :: from_pe_list(update_x%nrecv+update_y%nrecv), to_pe_list(update_x%nsend+update_y%nsend)
  integer            :: start_pos_recv(update_x%nrecv+update_y%nrecv), start_pos_send(update_x%nsend+update_y%nsend)
  MPP_TYPE_          :: fieldx(update_x%xbegin:update_x%xend, update_x%ybegin:update_x%yend,ke_max)
  MPP_TYPE_          :: fieldy(update_y%xbegin:update_y%xend, update_y%ybegin:update_y%yend,ke_max)
  MPP_TYPE_          :: buffer(size(mpp_domains_stack_nonblock(:)))

  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)
  pointer( ptr, buffer )

  update_edge_only = BTEST(flags, EDGEONLY)
  recv = .false.
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
  l_size = size(f_addrsx,1)
  nlist  = size(domain%list(:))
  ptr    = LOC(mpp_domains_stack_nonblock)

  nrecv = get_vector_recv(domain, update_x, update_y, ind_recv_x, ind_recv_y, start_pos_recv, from_pe_list)
  nsend = get_vector_send(domain, update_x, update_y, ind_send_x, ind_send_y, start_pos_send, to_pe_list)
  if( nrecv > MAX_REQUEST  ) then
     write( text,'(a,i8,a,i8)' ) 'nrecv =', nrecv, ' greater than MAX_REQEUST =', MAX_REQUEST
     call mpp_error(FATAL,'MPP_START_DO_UPDATE_V: '//trim(text))
  endif
  if( nsend > MAX_REQUEST ) then
     write( text,'(a,i8,a,i8)' ) 'nsend =', nsend, ' greater than MAX_REQEUST =', MAX_REQUEST
     call mpp_error(FATAL,'MPP_START_DO_UPDATE_V: '//trim(text))
  endif
  !--- make sure the domain stack size is big enough.
  buffer_pos = nonblock_data(id_update)%recv_pos
  recvsize = 0
  do m = 1, nrecv
     msgsize = 0
     nonblock_data(id_update)%size_recv(m) = 0
     ind_x = ind_recv_x(m)
     ind_y = ind_recv_y(m)
     if(ind_x >= 0) then
        do n = 1, update_x%recv(ind_x)%count
           dir = update_x%recv(ind_x)%dir(n)
           if(recv(dir)) then
              msgsize = msgsize + update_x%recv(ind_x)%msgsize(n)
           end if
        end do
     endif
     if(ind_y >= 0) then
        do n = 1, update_y%recv(ind_y)%count
           dir = update_y%recv(ind_y)%dir(n)
           if(recv(dir)) then
              msgsize = msgsize + update_y%recv(ind_y)%msgsize(n)
           end if
        end do
     endif
     if( msgsize.GT.0 )then
        msgsize = msgsize*ke_sum
        recvsize = recvsize + msgsize
        nonblock_data(id_update)%size_recv(m) = msgsize
        nonblock_data(id_update)%buffer_pos_recv(m) = buffer_pos
        buffer_pos = buffer_pos + msgsize
     end if
  end do   

  sendsize = 0
  do m = 1, nsend
     msgsize = 0
     ind_x = ind_send_x(m)
     ind_y = ind_send_y(m)
     if(ind_x >= 0) then
        do n = 1, update_x%send(ind_x)%count
           dir = update_x%send(ind_x)%dir(n)
           if(send(dir)) then
              msgsize = msgsize + update_x%send(ind_x)%msgsize(n)
           end if
        end do
     endif
     if(ind_y >= 0) then
        do n = 1, update_y%send(ind_y)%count
           dir = update_y%send(ind_y)%dir(n)
           if(send(dir)) then
              msgsize = msgsize + update_y%send(ind_y)%msgsize(n)
           end if
        end do
     endif
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
     call mpp_error( FATAL, 'MPP_START_DO_UPDATE_V: mpp_domains_stack overflow, '// &
          'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
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

  !--- recv
  call mpp_clock_begin(recv_clock_nonblock)
  do m = 1, nrecv
     msgsize = nonblock_data(id_update)%size_recv(m)
     from_pe = from_pe_list(m)
     if( msgsize .GT. 0 )then
        buffer_pos = nonblock_data(id_update)%buffer_pos_recv(m)
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., &
                   tag=id_update, request=request)
        nonblock_data(id_update)%request_recv(m) = request
#ifdef use_libMPI
        nonblock_data(id_update)%type_recv(m) = MPI_TYPE_
#endif
     end if
  end do

  call mpp_clock_end(recv_clock_nonblock)

  !--- send

  call mpp_clock_begin(send_pack_clock_nonblock)

!$OMP parallel do schedule(dynamic) default(shared) private(ind_x,ind_y,buffer_pos,pos,dir,tMe, &
!$OMP          is,ie,js,je,ptr_fieldx,ptr_fieldy)
  do m = 1, nsend
     send_msgsize(m) = 0
     ind_x = ind_send_x(m)
     ind_y = ind_send_y(m)
     buffer_pos = nonblock_data(id_update)%buffer_pos_send(m)
     pos = buffer_pos
     
     select case( gridtype )
     case(BGRID_NE, BGRID_SW, AGRID)
        if(ind_x >= 0) then
           do n = 1, update_x%send(ind_x)%count
              dir = update_x%send(ind_x)%dir(n)
              if( send(dir) ) then 
                 tMe = update_x%send(ind_x)%tileMe(n)
                 is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                 js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                 select case( update_x%send(ind_x)%rotation(n) )
                 case(ZERO)
                    do l=1,l_size  ! loop over number of fields
                       ptr_fieldx = f_addrsx(l,tMe)
                       ptr_fieldy = f_addrsy(l,tMe)
                       do k = 1,ke_list(l,tMe)
                          do j = js, je
                             do i = is, ie
                                pos = pos + 2
                                buffer(pos-1) = fieldx(i,j,k)
                                buffer(pos)   = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 case( MINUS_NINETY ) 
                    if( BTEST(flags,SCALAR_BIT) ) then
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do i = is, ie
                                do j = je, js, -1
                                   pos = pos + 2
                                   buffer(pos-1) = fieldy(i,j,k)
                                   buffer(pos)   = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    else
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do i = is, ie
                                do j = je, js, -1
                                   pos = pos + 2
                                   buffer(pos-1) = -fieldy(i,j,k)
                                   buffer(pos)   =  fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    end if
                 case( NINETY )
                    if( BTEST(flags,SCALAR_BIT) ) then
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do i = ie, is, -1
                                do j = js, je
                                   pos = pos + 2
                                   buffer(pos-1) = fieldy(i,j,k)
                                   buffer(pos)   = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    else
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do i = ie, is, -1
                                do j = js, je
                                   pos = pos + 2
                                   buffer(pos-1) = fieldy(i,j,k)
                                   buffer(pos)   = -fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    end if
                 case( ONE_HUNDRED_EIGHTY )
                    if( BTEST(flags,SCALAR_BIT) ) then
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 2
                                   buffer(pos-1) =  fieldx(i,j,k)
                                   buffer(pos)   =  fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    else
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 2
                                   buffer(pos-1) =  -fieldx(i,j,k)
                                   buffer(pos)   =  -fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    end if
                 end select ! select case( rotation(n) )
              end if ! if( send(dir) ) 
           end do ! do n = 1, update_x%send(ind_x)%count
        endif
     case(CGRID_NE, CGRID_SW)
        if(ind_x>=0) then
           do n = 1, update_x%send(ind_x)%count
              dir = update_x%send(ind_x)%dir(n)
              if( send(dir) ) then
                 tMe = update_x%send(ind_x)%tileMe(n)
                 is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                 js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                 select case( update_x%send(ind_x)%rotation(n) )
                 case(ZERO)
                    do l=1,l_size  ! loop over number of fields
                       ptr_fieldx = f_addrsx(l,tMe)
                       ptr_fieldy = f_addrsy(l,tMe)
                       do k = 1,ke_list(l,tMe)
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldx(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 case(MINUS_NINETY)
                    if( BTEST(flags,SCALAR_BIT) ) then
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do i = is, ie
                                do j = je, js, -1
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    else
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do i = is, ie
                                do j = je, js, -1
                                   pos = pos + 1
                                   buffer(pos) = -fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    end if
                 case(NINETY)
                    do l=1,l_size  ! loop over number of fields
                       ptr_fieldx = f_addrsx(l,tMe)
                       ptr_fieldy = f_addrsy(l,tMe)
                       do k = 1, ke_list(l,tMe)
                          do i = ie, is, -1
                             do j = js, je
                                pos = pos + 1
                                buffer(pos) = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 case(ONE_HUNDRED_EIGHTY)
                    if( BTEST(flags,SCALAR_BIT) ) then
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    else
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 1
                                   buffer(pos) = -fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    end if
                 end select
              end if
           end do
        endif
        if(ind_y>=0) then
           do n = 1, update_y%send(ind_y)%count
              dir = update_y%send(ind_y)%dir(n)
              if( send(dir) ) then
                 tMe = update_y%send(ind_y)%tileMe(n)
                 is = update_y%send(ind_y)%is(n); ie = update_y%send(ind_y)%ie(n)
                 js = update_y%send(ind_y)%js(n); je = update_y%send(ind_y)%je(n)
                 select case( update_y%send(ind_y)%rotation(n) )
                 case(ZERO)
                    do l=1,l_size  ! loop over number of fields
                       ptr_fieldx = f_addrsx(l,tMe)
                       ptr_fieldy = f_addrsy(l,tMe)
                       do k = 1,ke_list(l,tMe)
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 case(MINUS_NINETY)
                    do l=1,l_size  ! loop over number of fields
                       ptr_fieldx = f_addrsx(l,tMe)
                       ptr_fieldy = f_addrsy(l,tMe)
                       do k = 1,ke_list(l,tMe)
                          do i = is, ie
                             do j = je, js, -1
                                pos = pos + 1
                                buffer(pos) = fieldx(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 case(NINETY)
                    if( BTEST(flags,SCALAR_BIT) ) then
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do i = ie, is, -1
                                do j = js, je
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    else
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do i = ie, is, -1
                                do j = js, je
                                   pos = pos + 1
                                   buffer(pos) = -fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    end if
                 case(ONE_HUNDRED_EIGHTY)
                    if( BTEST(flags,SCALAR_BIT) ) then
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    else
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 1
                                   buffer(pos) = -fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    end if
                 end select
              endif
           enddo
        endif
     end select
     send_msgsize(m) = pos - buffer_pos
  enddo
!$OMP end parallel do
  do m = 1, nsend
     msgsize = send_msgsize(m)
     to_pe = to_pe_list(m)
     buffer_pos = nonblock_data(id_update)%buffer_pos_send(m)
     if( msgsize .GT.0 )then
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, &
                       tag=id_update, request=request )
        nonblock_data(id_update)%request_send(m) = request
     end if
  end do

  call mpp_clock_end(send_pack_clock_nonblock)


end subroutine MPP_START_DO_UPDATE_3D_V_

!###############################################################################
subroutine MPP_COMPLETE_DO_UPDATE_3D_V_(id_update, f_addrsx, f_addrsy, domain, update_x, update_y,     &
                                        d_type, ke_max, ke_list, gridtype, flags) 
  integer,             intent(in) :: id_update
  integer(LONG_KIND),  intent(in) :: f_addrsx(:,:), f_addrsy(:,:)
  type(domain2d),      intent(in) :: domain
  type(overlapSpec),   intent(in) :: update_x, update_y
  integer,             intent(in) :: ke_max
  integer,             intent(in) :: ke_list(:,:)
  MPP_TYPE_,           intent(in) :: d_type  ! creates unique interface
  integer,             intent(in) :: gridtype
  integer,             intent(in) :: flags


  !--- local variables
  MPP_TYPE_ :: fieldx(update_x%xbegin:update_x%xend, update_x%ybegin:update_x%yend,ke_max)
  MPP_TYPE_ :: fieldy(update_y%xbegin:update_y%xend, update_y%ybegin:update_y%yend,ke_max)
  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)

  MPP_TYPE_ :: recv_buffer(size(mpp_domains_stack_nonblock(:)))
  pointer( ptr, recv_buffer )

  integer :: i, j, k, l, is, ie, js, je, n, ke_sum, l_size, m
  integer :: pos, nlist, msgsize, tile, buffer_pos
  integer :: ind_x, ind_y, nrecv, nsend
  integer :: ind_recv_x(update_x%nrecv+update_y%nrecv), ind_recv_y(update_x%nrecv+update_y%nrecv)
  integer :: start_pos_recv(update_x%nrecv+update_y%nrecv)
  integer :: from_pe_list(update_x%nrecv+update_y%nrecv)
  logical :: recv(8), send(8), update_edge_only
  integer :: shift, midpoint
  integer :: tMe, dir

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
  l_size = size(f_addrsx,1)
  nlist  = size(domain%list(:))
  ptr = LOC(mpp_domains_stack_nonblock)

  nrecv = get_vector_recv(domain, update_x, update_y, ind_recv_x, ind_recv_y, start_pos_recv, from_pe_list)

  if(nrecv > 0) then
     call mpp_clock_begin(wait_clock_nonblock)
     call mpp_sync_self(check=EVENT_RECV, request=nonblock_data(id_update)%request_recv(1:nrecv), &
                        msg_size=nonblock_data(id_update)%size_recv(1:nrecv),                     &
                        msg_type=nonblock_data(id_update)%type_recv(1:nrecv) )
     call mpp_clock_end(wait_clock_nonblock)
#ifdef use_libMPI
     nonblock_data(id_update)%request_recv(:)    = MPI_REQUEST_NULL
#else
     nonblock_data(id_update)%request_recv(:)    = 0
#endif
     nonblock_data(id_update)%type_recv(:)    = 0
  endif 

  call mpp_clock_begin(unpk_clock_nonblock)
!$OMP parallel do schedule(dynamic) default(shared) private(ind_x,ind_y,buffer_pos,pos,dir,tMe,is,ie,js,je, &
!$OMP          msgsize,ptr_fieldx,ptr_fieldy)
  do m = nrecv,1,-1
     ind_x = ind_recv_x(m)
     ind_y = ind_recv_y(m)
     buffer_pos = nonblock_data(id_update)%buffer_pos_recv(m)+nonblock_data(id_update)%size_recv(m) 
     pos = buffer_pos
     select case ( gridtype )
     case(BGRID_NE, BGRID_SW, AGRID)
        if(ind_x>=0) then
           do n = update_x%recv(ind_x)%count, 1, -1    
              dir = update_x%recv(ind_x)%dir(n)
              if( recv(dir) ) then
                 tMe = update_x%recv(ind_x)%tileMe(n)
                 is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                 js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n) 
                 msgsize = (ie-is+1)*(je-js+1)*ke_sum*2
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 do l=1, l_size  ! loop over number of fields
                    ptr_fieldx = f_addrsx(l, tMe)
                    ptr_fieldy = f_addrsy(l, tMe)
                    do k = 1,ke_list(l,tMe)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 2
                             fieldx(i,j,k) = recv_buffer(pos-1)
                             fieldy(i,j,k) = recv_buffer(pos)
                          end do
                       end do
                    enddo
                 end do
              end if ! end if( recv(dir) )
           end do  ! do dir=8,1,-1 
        endif
     case(CGRID_NE, CGRID_SW)
        if(ind_y>=0) then
           do n = update_y%recv(ind_y)%count, 1, -1

              dir = update_y%recv(ind_y)%dir(n)
              if( recv(dir) ) then
                 tMe = update_y%recv(ind_y)%tileMe(n)
                 is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
                 js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)
                 msgsize = (ie-is+1)*(je-js+1)*ke_sum
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 do l=1, l_size  ! loop over number of fields
                    ptr_fieldy = f_addrsy(l, tMe)
                    do k = 1,ke_list(l,tMe)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             fieldy(i,j,k) = recv_buffer(pos)
                          end do
                       end do
                    end do
                 end do
              end if
           end do
        endif
        if(ind_x>=0) then
           do n = update_x%recv(ind_x)%count, 1, -1
              dir = update_x%recv(ind_x)%dir(n)
              if( recv(dir) ) then
                 tMe = update_x%recv(ind_x)%tileMe(n)
                 is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                 js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n) 
                 msgsize = (ie-is+1)*(je-js+1)*ke_sum
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 do l=1, l_size  ! loop over number of fields
                    ptr_fieldx = f_addrsx(l, tMe)
                    do k = 1,ke_list(l,tMe)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             fieldx(i,j,k) = recv_buffer(pos)
                          end do
                       end do
                    end do
                 end do
              end if
           end do
        endif
     end select
  end do
!$OMP end parallel do
  call mpp_clock_end(unpk_clock_nonblock)
  ! ---northern boundary fold
  shift = 0
  tMe = 1
  if(domain%symmetry) shift = 1
  if( BTEST(domain%fold,NORTH) .AND. (.NOT.BTEST(flags,SCALAR_BIT)) )then
     j = domain%y(1)%global%end+shift
     if( domain%y(1)%data%begin.LE.j .AND. j.LE.domain%y(1)%data%end+shift )then !fold is within domain
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
           j  = domain%y(1)%global%end+shift
           is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
           if( .NOT. domain%symmetry ) is = is - 1
           do i = is ,ie, midpoint
              if( domain%x(1)%data%begin.LE.i .AND. i.LE. domain%x(1)%data%end+shift )then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       fieldx(i,j,k) = 0.
                       fieldy(i,j,k) = 0.
                    end do
                 enddo
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off west edge, when update north or west direction
        j = domain%y(1)%global%end+shift 
        if ( recv(7) .OR. recv(5) ) then
           select case(gridtype)
           case(BGRID_NE)
              if(domain%symmetry) then
                 is = domain%x(1)%global%begin
              else
                 is = domain%x(1)%global%begin - 1
              end if
              if( is.GT.domain%x(1)%data%begin )then

                 if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-north BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)

                        do k = 1,ke_list(l,tMe)
                           do i = domain%x(1)%data%begin,is-1
                              fieldx(i,j,k) = fieldx(2*is-i,j,k)
                              fieldy(i,j,k) = fieldy(2*is-i,j,k)
                           end do
                        end do
                     end do
              end if
           case(CGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-north CGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do i = domain%x(1)%data%begin,is-1
                          fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                       end do
                    end do
                 end do
              end if
           end select
        end if

        !off east edge
        is = domain%x(1)%global%end
        if(domain%x(1)%cyclic .AND. is.LT.domain%x(1)%data%end )then
           ie = domain%x(1)%data%end
           is = is + 1
           select case(gridtype)
           case(BGRID_NE)
              is = is + shift
              ie = ie + shift
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do i = is,ie
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           case(CGRID_NE)
              do l=1,l_size
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do i = is, ie
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,SOUTH) .AND. (.NOT.BTEST(flags,SCALAR_BIT)) )then      ! ---southern boundary fold
     ! NOTE: symmetry is assumed for fold-south boundary
     j = domain%y(1)%global%begin
     if( domain%y(1)%data%begin.LE.j .AND. j.LE.domain%y(1)%data%end+shift )then !fold is within domain
        midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           j  = domain%y(1)%global%begin
           is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
           do i = is ,ie, midpoint
              if( domain%x(1)%data%begin.LE.i .AND. i.LE. domain%x(1)%data%end+shift )then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       fieldx(i,j,k) = 0.
                       fieldy(i,j,k) = 0.
                    end do
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off west edge, when update north or west direction
        j = domain%y(1)%global%begin
        if ( recv(3) .OR. recv(5) ) then
           select case(gridtype)
           case(BGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-south BGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do i = domain%x(1)%data%begin,is-1
                          fieldx(i,j,k) = fieldx(2*is-i,j,k)
                          fieldy(i,j,k) = fieldy(2*is-i,j,k)
                       end do
                    end do
                 end do
              end if
           case(CGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-south CGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do i = domain%x(1)%data%begin,is-1
                          fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                       end do
                    end do
                 end do
              end if
           end select
        end if

        !off east edge
        is = domain%x(1)%global%end
        if(domain%x(1)%cyclic .AND. is.LT.domain%x(1)%data%end )then
           ie = domain%x(1)%data%end
           is = is + 1
           select case(gridtype)
           case(BGRID_NE)
              is = is + shift
              ie = ie + shift
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do i = is,ie
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           case(CGRID_NE)
              do l=1,l_size
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do i = is, ie
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,WEST) .AND. (.NOT.BTEST(flags,SCALAR_BIT)) )then      ! ---eastern boundary fold
     ! NOTE: symmetry is assumed for fold-west boundary
     i = domain%x(1)%global%begin
     if( domain%x(1)%data%begin.LE.i .AND. i.LE.domain%x(1)%data%end+shift )then !fold is within domain
        midpoint = (domain%y(1)%global%begin+domain%y(1)%global%end-1+shift)/2
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           i  = domain%x(1)%global%begin
           js = domain%y(1)%global%begin; je = domain%y(1)%global%end+shift
           do j = js ,je, midpoint
              if( domain%y(1)%data%begin.LE.j .AND. j.LE. domain%y(1)%data%end+shift )then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       fieldx(i,j,k) = 0.
                       fieldy(i,j,k) = 0.
                    end do
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off south edge, when update south or west direction
        i = domain%x(1)%global%begin
        if ( recv(3) .OR. recv(5) ) then
           select case(gridtype)
           case(BGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then

                 if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-west BGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do j = domain%y(1)%data%begin,js-1
                          fieldx(i,j,k) = fieldx(i,2*js-j,k)
                          fieldy(i,j,k) = fieldy(i,2*js-j,k)
                       end do
                    end do
                 end do
              end if
           case(CGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then
                 if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-west CGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do j = domain%y(1)%data%begin,js-1
                          fieldx(i,j,k) = fieldx(i, 2*js-j-1,k)
                       end do
                    end do
                 end do
              end if
           end select
        end if

        !off north edge
        js = domain%y(1)%global%end
        if(domain%y(1)%cyclic .AND. js.LT.domain%y(1)%data%end )then
           je = domain%y(1)%data%end
           js = js + 1
           select case(gridtype)
           case(BGRID_NE)
              js = js + shift
              je = je + shift
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do j = js,je
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           case(CGRID_NE)
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do j = js, je
                       fieldx(i,j,k) = -fieldx(i,j,k)
                    end do
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,EAST) .AND. (.NOT.BTEST(flags,SCALAR_BIT)) )then      ! ---eastern boundary fold
     ! NOTE: symmetry is assumed for fold-west boundary
     i = domain%x(1)%global%end+shift
     if( domain%x(1)%data%begin.LE.i .AND. i.LE.domain%x(1)%data%end+shift )then !fold is within domain
        midpoint = (domain%y(1)%global%begin+domain%y(1)%global%end-1+shift)/2
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           i  = domain%x(1)%global%end+shift
           js = domain%y(1)%global%begin; je = domain%y(1)%global%end+shift
           do j = js ,je, midpoint
              if( domain%y(1)%data%begin.LE.j .AND. j.LE. domain%y(1)%data%end+shift )then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       fieldx(i,j,k) = 0.
                       fieldy(i,j,k) = 0.
                    end do
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off south edge, when update south or west direction
        i = domain%x(1)%global%end+shift
        if ( recv(3) .OR. recv(1) ) then
           select case(gridtype)
           case(BGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then

                 if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-east BGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do j = domain%y(1)%data%begin,js-1
                          fieldx(i,j,k) = fieldx(i,2*js-j,k)
                          fieldy(i,j,k) = fieldy(i,2*js-j,k)
                       end do
                    end do
                 end do
              end if
           case(CGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then
                 if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-east CGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do j = domain%y(1)%data%begin,js-1
                          fieldx(i,j,k) = fieldx(i, 2*js-j-1,k)
                       end do
                    end do
                 end do
              end if
           end select
        end if

        !off north edge
        js = domain%y(1)%global%end
        if(domain%y(1)%cyclic .AND. js.LT.domain%y(1)%data%end )then
           je = domain%y(1)%data%end
           js = js + 1
           select case(gridtype)
           case(BGRID_NE)
              js = js + shift
              je = je + shift
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do j = js,je
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           case(CGRID_NE)
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do j = js, je
                       fieldx(i,j,k) = -fieldx(i,j,k)
                    end do
                 end do
              end do
           end select
        end if
     end if
  end if

  
  if(nrecv>0) then
    nonblock_data(id_update)%size_recv(:) = 0
  endif

  nsend = update_x%nsend+update_y%nsend
  if(nsend > 0) then
     call mpp_clock_begin(wait_clock_nonblock)
     call mpp_sync_self(check=EVENT_SEND, request=nonblock_data(id_update)%request_send(1:nsend))
     call mpp_clock_end(wait_clock_nonblock)
     nonblock_data(id_update)%request_send_count = 0
#ifdef use_libMPI
     nonblock_data(id_update)%request_send(:)    = MPI_REQUEST_NULL
#else
     nonblock_data(id_update)%request_send(:)    = 0
#endif
  endif 

  return

end subroutine MPP_COMPLETE_DO_UPDATE_3D_V_
