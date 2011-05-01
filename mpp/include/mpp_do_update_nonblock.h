! -*-f90-*- 
subroutine MPP_START_DO_UPDATE_3D_V_(id_update, fieldx, fieldy, domain, update_x, update_y, gridtype, update_flags, &
                                  tile_count)
  integer,             intent(in)        :: id_update
  type(domain2d),      intent(in)        :: domain
  MPP_TYPE_,        intent(inout)        :: fieldx(domain%x(1)%data%begin:,domain%y(1)%data%begin:,:)
  MPP_TYPE_,        intent(inout)        :: fieldy(domain%x(1)%data%begin:,domain%y(1)%data%begin:,:)
  type(overlapSpec),   intent(in)        :: update_x, update_y
  integer, intent(in)                    :: gridtype
  integer, intent(in)                    :: update_flags
  integer, intent(in),          optional :: tile_count

  MPP_TYPE_ :: buffer(size(mpp_domains_stack_nonblock(:)))
  pointer( ptr, buffer )

  integer :: i, j, k, is, ie, js, je, ke, n
  integer :: pos, nlist, msgsize, tile
  integer :: to_pe, from_pe
  integer :: tMe, dir, count
  logical :: send(8), recv(8)
  character(len=8) :: text
  integer :: rank_x, rank_y, ind_x, ind_y, cur_rank
  integer :: nsend_x, nsend_y, nrecv_x, nrecv_y

  recv(1) = BTEST(update_flags,EAST)
  recv(3) = BTEST(update_flags,SOUTH)
  recv(5) = BTEST(update_flags,WEST)
  recv(7) = BTEST(update_flags,NORTH)
  recv(2) = recv(1) .AND. recv(3)
  recv(4) = recv(3) .AND. recv(5)
  recv(6) = recv(5) .AND. recv(7)
  recv(8) = recv(7) .AND. recv(1)
  send    = recv


  ke = size(fieldx,3)
  tile = 1
  if(PRESENT(tile_count)) tile = tile_count

  nlist = size(domain%list(:))
  ptr = LOC(mpp_domains_stack_nonblock)

  !recv
  nsend_x = update_x%nsend
  nsend_y = update_y%nsend
  nrecv_x = update_x%nrecv
  nrecv_y = update_y%nrecv

  !--- recv
  cur_rank = get_rank_recv(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y) 

  call mpp_clock_begin(recv_clock_nonblock)
  do while (ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y)
     msgsize = 0
     select case(gridtype)
     case(BGRID_NE, BGRID_SW, AGRID)
        if(cur_rank == rank_x) then
           from_pe = update_x%recv(ind_x)%pe
           do n = 1, update_x%recv(ind_x)%count
              tMe = update_x%recv(ind_x)%tileMe(n)
              if( tMe .NE. tile ) cycle
              dir = update_x%recv(ind_x)%dir(n)
              if(recv(dir)) then
                 msgsize = msgsize + update_x%recv(ind_x)%msgsize(n)
              end if
           end do
           msgsize = msgsize*2
           ind_x = ind_x+1
           ind_y = ind_x
           if(ind_x .LE. nrecv_x) then
              rank_x = update_x%recv(ind_x)%pe - domain%pe 
              if(rank_x .LE.0) rank_x = rank_x + nlist
           else
              rank_x = -1
           endif
           rank_y = rank_x
        endif
     case(CGRID_NE, CGRID_SW)
        if(cur_rank == rank_x) then
           from_pe = update_x%recv(ind_x)%pe
           do n = 1, update_x%recv(ind_x)%count
              tMe = update_x%recv(ind_x)%tileMe(n)
              if( tMe .NE. tile ) cycle
              dir = update_x%recv(ind_x)%dir(n)
              if(recv(dir)) then
                 msgsize = msgsize + update_x%recv(ind_x)%msgsize(n)
              end if
           end do
           ind_x = ind_x+1
           if(ind_x .LE. nrecv_x) then
              rank_x = update_x%recv(ind_x)%pe - domain%pe 
              if(rank_x .LE.0) rank_x = rank_x + nlist
           else
              rank_x = -1
           endif
        endif
        if(cur_rank == rank_y) then
           from_pe = update_y%recv(ind_y)%pe
           do n = 1, update_y%recv(ind_y)%count
              tMe = update_y%recv(ind_y)%tileMe(n)
              if( tMe .NE. tile ) cycle
              dir = update_y%recv(ind_y)%dir(n)
              if(recv(dir)) then
                 msgsize = msgsize + update_y%recv(ind_y)%msgsize(n)
              end if
           end do
           ind_y = ind_y+1
           if(ind_y .LE. nrecv_y) then
              rank_y = update_y%recv(ind_y)%pe - domain%pe 
              if(rank_y .LE.0) rank_y = rank_y + nlist
           else
              rank_y = -1
           endif
        endif
     end select
     cur_rank = max(rank_x, rank_y)
     msgsize = msgsize*ke

     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, nonblock_buffer_pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_START_DO_UPDATE_V: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        count = request_recv(id_update)%count + 1
        if( count > MAX_REQUEST ) then
           write( text,'(a,i,a,i)' ) 'request count =', count, ' greater than MAX_REQEUST =', MAX_REQUEST
           call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: '//trim(text))
        endif
        request_recv(id_update)%count = count

        call mpp_recv( buffer(nonblock_buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., tag=NONBLOCK_UPDATE_TAG, &
                   request=request_recv(id_update)%request(count))
        nonblock_buffer_pos = nonblock_buffer_pos + msgsize
     end if
  end do
  call mpp_clock_end(recv_clock_nonblock)

  recv_pos_list(id_update) = nonblock_buffer_pos  

  !--- send
  cur_rank = get_rank_send(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y) 

  do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
     call mpp_clock_begin(pack_clock_nonblock)
     pos = nonblock_buffer_pos
     !--- make sure the domain stack size is big enough
     msgsize = 0
     if(cur_rank == rank_x) then
        do n = 1, update_x%send(ind_x)%count
           tMe = update_x%send(ind_x)%tileMe(n)
           if( tMe .NE. tile ) cycle
           dir = update_x%send(ind_x)%dir(n)
           if( send(dir) ) msgsize = msgsize +  update_x%send(ind_x)%msgsize(n)
        enddo
     endif
     if(cur_rank == rank_y) then
        do n = 1, update_y%send(ind_y)%count
           tMe = update_y%send(ind_y)%tileMe(n)
           if( tMe .NE. tile ) cycle
           dir = update_y%send(ind_y)%dir(n)
           if( send(dir) ) msgsize = msgsize +  update_y%send(ind_y)%msgsize(n)
        enddo
     endif

     if( msgsize.GT.0 )then
        msgsize = msgsize*ke
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_START_DO_UPDATE_V: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if
     
     select case( gridtype )
     case(BGRID_NE, BGRID_SW, AGRID)
        if(cur_rank == rank_x) then
           to_pe = update_x%send(ind_x)%pe
           do n = 1, update_x%send(ind_x)%count
              tMe = update_x%send(ind_x)%tileMe(n)
              if( tMe .NE. tile ) cycle
              dir = update_x%send(ind_x)%dir(n)
              if( send(dir) ) then   
                 is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                 js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                 if( update_x%send(ind_x)%is_refined(n) ) then
                    select case( update_x%send(ind_x)%rotation(n) )
                    case(ZERO)
                       do k = 1,ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 2
                                buffer(pos-1) = fieldx(i,j,k)
                                buffer(pos)   = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    case(MINUS_NINETY)
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 2
                                   buffer(pos-1) =  fieldy(i,j,k)
                                   buffer(pos)   =  fieldx(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 2
                                   buffer(pos-1) = -fieldy(i,j,k)
                                   buffer(pos)   =  fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    case(NINETY)
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 2
                                   buffer(pos-1) =  fieldy(i,j,k)
                                   buffer(pos)   =  fieldx(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 2
                                   buffer(pos-1) =  fieldy(i,j,k)
                                   buffer(pos)   = -fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    case(ONE_HUNDRED_EIGHTY)
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 2
                                   buffer(pos-1) = fieldx(i,j,k)
                                   buffer(pos)   = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 2
                                   buffer(pos-1) = -fieldx(i,j,k)
                                   buffer(pos)   = -fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    end select  ! select case( rotation(n) )
                 else   ! if( is_refined(n) )
                    select case( update_x%send(ind_x)%rotation(n) )
                    case(ZERO)
                       do k = 1,ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 2
                                buffer(pos-1) = fieldx(i,j,k)
                                buffer(pos)   = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    case( MINUS_NINETY ) 
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do i = is, ie
                                do j = je, js, -1
                                   pos = pos + 2
                                   buffer(pos-1) = fieldy(i,j,k)
                                   buffer(pos)   = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do i = is, ie
                                do j = je, js, -1
                                   pos = pos + 2
                                   buffer(pos-1) = -fieldy(i,j,k)
                                   buffer(pos)   =  fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    case( NINETY )
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do i = ie, is, -1
                                do j = js, je
                                   pos = pos + 2
                                   buffer(pos-1) = fieldy(i,j,k)
                                   buffer(pos)   = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do i = ie, is, -1
                                do j = js, je
                                   pos = pos + 2
                                   buffer(pos-1) = fieldy(i,j,k)
                                   buffer(pos)   = -fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    case( ONE_HUNDRED_EIGHTY )
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 2
                                   buffer(pos-1) =  fieldx(i,j,k)
                                   buffer(pos)   =  fieldy(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 2
                                   buffer(pos-1) =  -fieldx(i,j,k)
                                   buffer(pos)   =  -fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    end select ! select case( rotation(n) )
                 end if ! if( is_refined(n) )
              end if ! if( send(dir) ) 
           end do ! do n = 1, update_x%send(ind_x)%count
           ind_x = ind_x+1
           ind_y = ind_x
           if(ind_x .LE. nsend_x) then
              rank_x = update_x%send(ind_x)%pe - domain%pe 
              if(rank_x .LT.0) rank_x = rank_x + nlist
           else
              rank_x = nlist+1
           endif
           rank_y = rank_x
        endif
     case(CGRID_NE, CGRID_SW)
        if(cur_rank == rank_x) then
           to_pe = update_x%send(ind_x)%pe
           do n = 1, update_x%send(ind_x)%count
              tMe = update_x%send(ind_x)%tileMe(n)
              if( tMe .NE. tile ) cycle
              dir = update_x%send(ind_x)%dir(n)
              if( send(dir) ) then
                 is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                 js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                 if( update_x%send(ind_x)%is_refined(n) ) then
                    select case( update_x%send(ind_x)%rotation(n) )
                    case(ZERO)
                       do k = 1,ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldx(i,j,k)
                             end do
                          end do
                       end do
                    case(MINUS_NINETY)
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = -fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    case(NINETY)
                       do k = 1, ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    case(ONE_HUNDRED_EIGHTY)
                       do k = 1,ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldx(i,j,k)
                             end do
                          end do
                       end do
                    end select
                 else
                    select case( update_x%send(ind_x)%rotation(n) )
                    case(ZERO)
                       do k = 1,ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldx(i,j,k)
                             end do
                          end do
                       end do
                    case(MINUS_NINETY)
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do i = is, ie
                                do j = je, js, -1
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do i = is, ie
                                do j = je, js, -1
                                   pos = pos + 1
                                   buffer(pos) = -fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    case(NINETY)
                       do k = 1, ke
                          do i = ie, is, -1
                             do j = js, je
                                pos = pos + 1
                                buffer(pos) = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    case(ONE_HUNDRED_EIGHTY)
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 1
                                   buffer(pos) = -fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    end select
                 end if
              end if
           end do
           ind_x = ind_x+1
           if(ind_x .LE. nsend_x) then
              rank_x = update_x%send(ind_x)%pe - domain%pe 
              if(rank_x .LT.0) rank_x = rank_x + nlist
           else
              rank_x = nlist+1 
           endif
        endif
        if(cur_rank == rank_y) then
           to_pe = update_y%send(ind_y)%pe
           do n = 1, update_y%send(ind_y)%count
              dir = update_y%send(ind_y)%dir(n)
              if( send(dir) ) then
                 tMe = update_y%send(ind_y)%tileMe(n)
                 is = update_y%send(ind_y)%is(n); ie = update_y%send(ind_y)%ie(n)
                 js = update_y%send(ind_y)%js(n); je = update_y%send(ind_y)%je(n)
                 if( update_y%send(ind_y)%is_refined(n) ) then
                    select case( update_y%send(ind_y)%rotation(n) )
                    case(ZERO)
                       do k = 1,ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    case(MINUS_NINETY)
                       do k = 1,ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldx(i,j,k)
                             end do
                          end do
                       end do
                    case(NINETY)
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = -fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    case(ONE_HUNDRED_EIGHTY)
                       do k = 1,ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    end select
                 else
                    select case( update_y%send(ind_y)%rotation(n) )
                    case(ZERO)
                       do k = 1,ke
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    case(MINUS_NINETY)
                       do k = 1,ke
                          do i = is, ie
                             do j = je, js, -1
                                pos = pos + 1
                                buffer(pos) = fieldx(i,j,k)
                             end do
                          end do
                       end do
                    case(NINETY)
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do i = ie, is, -1
                                do j = js, je
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do i = ie, is, -1
                                do j = js, je
                                   pos = pos + 1
                                   buffer(pos) = -fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    case(ONE_HUNDRED_EIGHTY)
                       if( BTEST(update_flags,SCALAR_BIT) ) then
                          do k = 1,ke
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       else
                          do k = 1,ke
                             do j = je, js, -1
                                do i = ie, is, -1
                                   pos = pos + 1
                                   buffer(pos) = -fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end if
                    end select
                 end if
              endif
           enddo
           ind_y = ind_y+1
           if(ind_y .LE. nsend_y) then
              rank_y = update_y%send(ind_y)%pe - domain%pe 
              if(rank_y .LT.0) rank_y = rank_y + nlist
           else
              rank_y = nlist+1
           endif
        endif
     end select
     call mpp_clock_end(pack_clock_nonblock)
     call mpp_clock_begin(send_clock_nonblock)
     cur_rank = min(rank_x, rank_y)
     msgsize = pos - nonblock_buffer_pos
     if( msgsize.GT.0 )then
        call mpp_send( buffer(nonblock_buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=NONBLOCK_UPDATE_TAG )
        nonblock_buffer_pos = pos
     end if
     call mpp_clock_end(send_clock_nonblock)
  end do


end subroutine MPP_START_DO_UPDATE_3D_V_

!###############################################################################
subroutine MPP_COMPLETE_DO_UPDATE_3D_V_ (id_update, fieldx, fieldy, domain, update_x, update_y, gridtype, update_flags, &
                                      tile_count, bufferx, buffery )
  integer,             intent(in)        :: id_update
  type(domain2d),      intent(in)        :: domain
  MPP_TYPE_,        intent(inout)        :: fieldx(domain%x(1)%data%begin:,domain%y(1)%data%begin:,:)
  MPP_TYPE_,        intent(inout)        :: fieldy(domain%x(1)%data%begin:,domain%y(1)%data%begin:,:)
  type(overlapSpec),   intent(in)        :: update_x, update_y
  integer, intent(in)                    :: gridtype
  integer, intent(in)                    :: update_flags
  integer, intent(in),          optional :: tile_count
  MPP_TYPE_,     intent(inout), optional :: bufferx(:)  
  MPP_TYPE_,     intent(inout), optional :: buffery(:)  

  !--- local variables
  MPP_TYPE_ :: recv_buffer(size(mpp_domains_stack_nonblock(:)))
  pointer( ptr, recv_buffer )

  integer :: i, j, k, is, ie, js, je, ke, n
  integer :: pos, nlist, msgsize, tile, buffer_pos
  integer :: rank_x, rank_y, ind_x, ind_y, cur_rank
  integer :: index, is1, ie1, js1, je1, ni, nj, total, start1, start, start2
  logical :: recv(8)
  integer :: shift, midpoint
  integer :: tMe, dir


  ke = size(fieldx,3)
  tile = 1
  if(PRESENT(tile_count)) tile = tile_count

  buffer_pos = recv_pos_list(id_update)
  cur_rank = get_rank_unpack(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y) 

  recv(1) = BTEST(update_flags,EAST)
  recv(3) = BTEST(update_flags,SOUTH)
  recv(5) = BTEST(update_flags,WEST)
  recv(7) = BTEST(update_flags,NORTH)
  recv(2) = recv(1) .AND. recv(3)
  recv(4) = recv(3) .AND. recv(5)
  recv(6) = recv(5) .AND. recv(7)
  recv(8) = recv(7) .AND. recv(1)

  nlist = size(domain%list(:))
  ptr = LOC(mpp_domains_stack_nonblock)

  do while (ind_x > 0 .OR. ind_y > 0)
     call mpp_clock_begin(unpk_clock_nonblock)
     pos = buffer_pos
     select case ( gridtype )
     case(BGRID_NE, BGRID_SW, AGRID)
        if(cur_rank == rank_x) then
           do n = update_x%recv(ind_x)%count, 1, -1    
              tMe = update_x%recv(ind_x)%tileMe(n)
              if( tMe .NE. tile ) cycle
              dir = update_x%recv(ind_x)%dir(n)
              if( recv(dir) ) then
                 is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                 js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n) 
                 msgsize = (ie-is+1)*(je-js+1)*ke*2
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 if(update_x%recv(ind_x)%is_refined(n)) then
                 if(.not. present(bufferx) .or. .not. present(buffery)) then
                    call mpp_error(FATAL,  &
                      "MPP_COMPLETE_DO_UPDATE_V: refinement exists but argument bufferx or buffery is not present")
                 endif
                 index = update_x%recv(ind_x)%index(n)
                 is1 = update_x%rSpec(tMe)%isNbr(index); ie1 = update_x%rSpec(tMe)%ieNbr(index)
                 js1 = update_x%rSpec(tMe)%jsNbr(index); je1 = update_x%rSpec(tMe)%jeNbr(index)
                 ni = ie1 - is1 + 1
                 nj = je1 - js1 + 1
                 total = ni*nj
                 start = (update_x%rSpec(tMe)%start(index)-1)*ke

                 if(start+total*ke>size(bufferx) .or. start+total*ke>size(buffery) ) call mpp_error(FATAL, &
                      "MPP_COMPLETE_DO_UPDATE_3D_V: size of bufferx or buffery is less than the size of the data to be filled.")
                 msgsize = ie - is + 1
                 start1 = start + (js-js1)*ni + is - is1 
                 do k = 1, ke
                    start2 = start1
                    do j = js, je
                       do i = start2+1, start2+msgsize
                          pos = pos + 2
                          bufferx(i) = recv_buffer(pos-1)
                          buffery(i) = recv_buffer(pos)
                       end do
                       start2 = start2 + ni
                    end do
                    start1 = start1 + total
                 end do
                 else
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 2
                             fieldx(i,j,k) = recv_buffer(pos-1)
                             fieldy(i,j,k) = recv_buffer(pos)
                          end do
                       end do
                    end do
                 end if
              end if ! end if( recv(dir) )
           end do  ! do dir=8,1,-1 
           ind_x = ind_x-1
           ind_y = ind_x
           if(ind_x .GT. 0) then
              rank_x = update_x%recv(ind_x)%pe - domain%pe 
              if(rank_x .LE.0) rank_x = rank_x + nlist
           else
              rank_x = nlist+1
           endif
           rank_y = rank_x
        endif
     case(CGRID_NE, CGRID_SW)
        if(cur_rank == rank_y) then
           do n = update_y%recv(ind_y)%count, 1, -1
              tMe = update_y%recv(ind_y)%tileMe(n)
              if( tMe .NE. tile ) cycle
              dir = update_y%recv(ind_y)%dir(n)
              if( recv(dir) ) then
                 is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
                 js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)
                 msgsize = (ie-is+1)*(je-js+1)*ke
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 if(update_y%recv(ind_y)%is_refined(n)) then
                    if(.not. present(buffery)) then
                       call mpp_error(FATAL,  &
                            "MPP_COMPLETE_DO_UPDATE_V: refinement exists but argument buffery is not present")
                    endif
                    index = update_y%recv(ind_y)%index(n)
                    is1 = update_y%rSpec(tMe)%isNbr(index); ie1 = update_y%rSpec(tMe)%ieNbr(index)
                    js1 = update_y%rSpec(tMe)%jsNbr(index); je1 = update_y%rSpec(tMe)%jeNbr(index)
                    ni = ie1 - is1 + 1
                    nj = je1 - js1 + 1
                    total = ni*nj
                    start = (update_y%rSpec(tMe)%start(index)-1)*ke
                    if(start+total*ke>size(buffery) ) call mpp_error(FATAL, &
                      "MPP_COMPLETE_DO_UPDATE_3D_V: size of buffery is less than the size of the data to be filled.")

                    msgsize = ie - is + 1
                    start1 = start + (js-js1)*ni + is - is1 
                    do k = 1, ke
                       start2 = start1
                       do j = js, je
                          do i = start2+1, start2+msgsize
                             pos = pos + 1
                             buffery(i) = recv_buffer(pos)
                          end do
                          start2 = start2 + ni
                       end do
                       start1 = start1 + total
                    end do
                 else
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             fieldy(i,j,k) = recv_buffer(pos)
                          end do
                       end do
                    end do
                 end if
              end if
           end do
           ind_y = ind_y-1
           if(ind_y .GT. 0) then
              rank_y = update_y%recv(ind_y)%pe - domain%pe
              if(rank_y .LE.0) rank_y = rank_y + nlist
           else
              rank_y = nlist+1
           endif
        endif
        if(cur_rank == rank_x) then
           do n = update_x%recv(ind_x)%count, 1, -1
              tMe = update_x%recv(ind_x)%tileMe(n)
              if( tMe .NE. tile ) cycle
              dir = update_x%recv(ind_x)%dir(n)
              if( recv(dir) ) then
                 is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                 js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n) 
                 msgsize = (ie-is+1)*(je-js+1)*ke
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 if(update_x%recv(ind_x)%is_refined(n)) then
                    if(.not. present(bufferx)) then
                       call mpp_error(FATAL,  &
                            "MPP_COMPLETE_DO_UPDATE_V: refinement exists but argument bufferx is not present")
                    endif
                    index = update_x%recv(ind_x)%index(n)
                    is1 = update_x%rSpec(tMe)%isNbr(index); ie1 = update_x%rSpec(tMe)%ieNbr(index)
                    js1 = update_x%rSpec(tMe)%jsNbr(index); je1 = update_x%rSpec(tMe)%jeNbr(index)
                    ni = ie1 - is1 + 1
                    nj = je1 - js1 + 1
                    total = ni*nj
                    start = (update_x%rSpec(tMe)%start(index)-1)*ke
                    if(start+total*ke>size(bufferx) ) call mpp_error(FATAL, &
                      "MPP_COMPLETE_DO_UPDATE_3D_V: size of bufferx is less than the size of the data to be filled.")

                    msgsize = ie - is + 1
                    start1 = start + (js-js1)*ni + is - is1 
                    do k = 1, ke
                       start2 = start1
                       do j = js, je
                          do i = start2+1, start2+msgsize
                             pos = pos + 1
                             bufferx(i) = recv_buffer(pos)
                          end do
                          start2 = start2 + ni
                       end do
                       start1 = start1 + total
                    end do
                 else
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             fieldx(i,j,k) = recv_buffer(pos)
                          end do
                       end do
                    end do
                 end if
              end if
           end do
           ind_x = ind_x-1
           if(ind_x .GT. 0) then
              rank_x = update_x%recv(ind_x)%pe - domain%pe 
              if(rank_x .LE.0) rank_x = rank_x + nlist
           else
              rank_x = nlist+1
           endif
        endif
     end select
     cur_rank = min(rank_x, rank_y)     
     call mpp_clock_end(unpk_clock_nonblock)
  end do

  ! ---northern boundary fold
  shift = 0
  if(domain%symmetry) shift = 1
  if( BTEST(domain%fold,NORTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then
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
                 do k = 1,ke
                    fieldx(i,j,k) = 0.
                    fieldy(i,j,k) = 0.
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off west edge, when update north or west direction
        j = domain%y(1)%global%end+shift 
        if ( BTEST(update_flags,NORTH) .OR. BTEST(update_flags,WEST) ) then
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
                 do k = 1,ke
                    do i = domain%x(1)%data%begin,is-1
                       fieldx(i,j,k) = fieldx(2*is-i,j,k)
                       fieldy(i,j,k) = fieldy(2*is-i,j,k)
                    end do
                 end do
              end if
           case(CGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-north CGRID_NE west edge ubound error.' )
                 do k = 1,ke
                    do i = domain%x(1)%data%begin,is-1
                       fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
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
              do k = 1,ke
                 do i = is,ie
                    fieldx(i,j,k) = -fieldx(i,j,k)
                    fieldy(i,j,k) = -fieldy(i,j,k)
                 end do
              end do
           case(CGRID_NE)
              do k = 1,ke
                 do i = is, ie
                    fieldy(i,j,k) = -fieldy(i,j,k)
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,SOUTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then      ! ---southern boundary fold
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
                 do k = 1,ke
                    fieldx(i,j,k) = 0.
                    fieldy(i,j,k) = 0.
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off west edge, when update north or west direction
        j = domain%y(1)%global%begin
        if ( BTEST(update_flags,SOUTH) .OR. BTEST(update_flags,WEST) ) then
           select case(gridtype)
           case(BGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 
                 if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-south BGRID_NE west edge ubound error.' )
                 do k = 1,ke
                    do i = domain%x(1)%data%begin,is-1
                       fieldx(i,j,k) = fieldx(2*is-i,j,k)
                       fieldy(i,j,k) = fieldy(2*is-i,j,k)
                    end do
                 end do
              end if
           case(CGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-south CGRID_NE west edge ubound error.' )
                 do k = 1,ke
                    do i = domain%x(1)%data%begin,is-1
                       fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
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
              do k = 1,ke
                 do i = is,ie
                    fieldx(i,j,k) = -fieldx(i,j,k)
                    fieldy(i,j,k) = -fieldy(i,j,k)
                 end do
              end do
           case(CGRID_NE)
              do k = 1,ke
                 do i = is, ie
                    fieldy(i,j,k) = -fieldy(i,j,k)
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,WEST) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then      ! ---eastern boundary fold
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
                 do k = 1,ke
                    fieldx(i,j,k) = 0.
                    fieldy(i,j,k) = 0.
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off south edge, when update south or west direction
        i = domain%x(1)%global%begin
        if ( BTEST(update_flags,SOUTH) .OR. BTEST(update_flags,WEST) ) then
           select case(gridtype)
           case(BGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then

                 if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-west BGRID_NE west edge ubound error.' )
                 do k = 1,ke
                    do j = domain%y(1)%data%begin,js-1
                       fieldx(i,j,k) = fieldx(i,2*js-j,k)
                       fieldy(i,j,k) = fieldy(i,2*js-j,k)
                    end do
                 end do
              end if
           case(CGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then
                 if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-west CGRID_NE west edge ubound error.' )
                 do k = 1,ke
                    do j = domain%y(1)%data%begin,js-1
                       fieldx(i,j,k) = fieldx(i, 2*js-j-1,k)
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
              do k = 1,ke
                 do j = js,je
                    fieldx(i,j,k) = -fieldx(i,j,k)
                    fieldy(i,j,k) = -fieldy(i,j,k)
                 end do
              end do
           case(CGRID_NE)
              do k = 1,ke
                 do j = js, je
                    fieldx(i,j,k) = -fieldx(i,j,k)
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,EAST) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then      ! ---eastern boundary fold
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
                 do k = 1,ke
                    fieldx(i,j,k) = 0.
                    fieldy(i,j,k) = 0.
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off south edge, when update south or west direction
        i = domain%x(1)%global%end+shift
        if ( BTEST(update_flags,SOUTH) .OR. BTEST(update_flags,EAST) ) then
           select case(gridtype)
           case(BGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then

                 if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-east BGRID_NE west edge ubound error.' )
                 do k = 1,ke
                    do j = domain%y(1)%data%begin,js-1
                       fieldx(i,j,k) = fieldx(i,2*js-j,k)
                       fieldy(i,j,k) = fieldy(i,2*js-j,k)
                    end do
                 end do
              end if
           case(CGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then
                 if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-east CGRID_NE west edge ubound error.' )
                 do k = 1,ke
                    do j = domain%y(1)%data%begin,js-1
                       fieldx(i,j,k) = fieldx(i, 2*js-j-1,k)
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
              do k = 1,ke
                 do j = js,je
                    fieldx(i,j,k) = -fieldx(i,j,k)
                    fieldy(i,j,k) = -fieldy(i,j,k)
                 end do
              end do
           case(CGRID_NE)
              do k = 1,ke
                 do j = js, je
                    fieldx(i,j,k) = -fieldx(i,j,k)
                 end do
              end do
           end select
        end if
     end if
  end if

  return

end subroutine MPP_COMPLETE_DO_UPDATE_3D_V_
