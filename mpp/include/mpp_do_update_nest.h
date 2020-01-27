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
subroutine MPP_DO_UPDATE_NEST_FINE_3D_(f_addrs, nest_domain, update, d_type, ke, wb_addrs, eb_addrs, &
                                   sb_addrs, nb_addrs, flags, xbegin, xend, ybegin, yend)
  integer(LONG_KIND),         intent(in) :: f_addrs(:)
  type(nest_level_type),      intent(in) :: nest_domain
  type(nestSpec),          intent(in) :: update
  MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
  integer,                    intent(in) :: ke
  integer(LONG_KIND),         intent(in) :: wb_addrs(:)
  integer(LONG_KIND),         intent(in) :: eb_addrs(:)
  integer(LONG_KIND),         intent(in) :: sb_addrs(:)
  integer(LONG_KIND),         intent(in) :: nb_addrs(:)
  integer,                    intent(in) :: flags
  integer,                    intent(in) :: xbegin, xend, ybegin, yend

  character(len=8)            :: text
  type(overlap_type), pointer :: overPtr => NULL()
  logical   :: send(8), recv(8)
  integer   :: from_pe, to_pe, dir
  integer   :: m, n, l, i, j, k
  integer   :: is, ie, js, je, l_size
  integer   :: buffer_pos, msgsize
  integer   :: buffer_recv_size, pos
  MPP_TYPE_ :: field(xbegin:xend, ybegin:yend,ke)
  MPP_TYPE_ :: wbuffer(update%west%is_you :update%west%ie_you,  update%west%js_you :update%west%je_you, ke)
  MPP_TYPE_ :: ebuffer(update%east%is_you :update%east%ie_you,  update%east%js_you :update%east%je_you, ke)
  MPP_TYPE_ :: sbuffer(update%south%is_you:update%south%ie_you, update%south%js_you:update%south%je_you,ke)
  MPP_TYPE_ :: nbuffer(update%north%is_you:update%north%ie_you, update%north%js_you:update%north%je_you,ke)
  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))

  pointer(ptr_field, field)
  pointer(ptr_buffer, buffer )
  pointer(ptr_wbuffer, wbuffer)
  pointer(ptr_ebuffer, ebuffer)
  pointer(ptr_sbuffer, sbuffer)
  pointer(ptr_nbuffer, nbuffer)

  recv(1) = BTEST(flags,EAST)
  recv(3) = BTEST(flags,SOUTH)
  recv(5) = BTEST(flags,WEST)
  recv(7) = BTEST(flags,NORTH)

  send    = recv

  ptr_buffer = LOC(mpp_domains_stack)
  l_size = size(f_addrs(:))

  !--- pre-post receiving
  buffer_pos = 0
  do m = 1, update%nrecv
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(nest_recv_clock)
     msgsize = 0
     do n = 1, overPtr%count
        dir = overPtr%dir(n)
        if(recv(dir)) then
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           msgsize = msgsize + (ie-is+1)*(je-js+1)
        end if
     end do
     msgsize = msgsize*ke*l_size
     if( msgsize.GT.0 )then
        from_pe = overPtr%pe
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_FINE_3D_: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
     call mpp_clock_end(nest_recv_clock)
  end do ! end do m = 1, update%nrecv
  buffer_recv_size = buffer_pos

  !--- pack and send the data
  do m = 1, update%nsend
     overPtr => update%send(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(nest_pack_clock)
     pos = buffer_pos
     msgsize = 0
     do n = 1, overPtr%count
        dir = overPtr%dir(n)
        if( send(dir) )  msgsize = msgsize + overPtr%msgsize(n)
     enddo
     if( msgsize.GT.0 )then
        msgsize = msgsize*ke*l_size
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_FINE_3D_: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if

     do n = 1, overPtr%count
        dir = overPtr%dir(n)
        if( send(dir) ) then
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           select case(overPtr%rotation(n))
           case(ZERO)
           do l=1,l_size  ! loop over number of fields
              ptr_field = f_addrs(l)
              do k = 1,ke
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              end do
           end do
           case(MINUS_NINETY)
              do l=1,l_size  ! loop over number of fields
                 ptr_field = f_addrs(l)
                 do k = 1,ke
                    do i = is, ie
                       do j = je, js, -1
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case(NINETY)
              do l=1,l_size  ! loop over number of fields
                 ptr_field = f_addrs(l)
                 do k = 1,ke
                    do i = ie, is, -1
                       do j = js, je
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case default
              call mpp_error(FATAL, "MPP_DO_UPDATE_NEST_FINE_3D: pack rotate must be ZERO, MINUS_NINETY or NINETY")
           end select
        endif
     end do ! do n = 1, overPtr%count

     call mpp_clock_end(nest_pack_clock)
     call mpp_clock_begin(nest_send_clock)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then
        to_pe = overPtr%pe
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1 )
        buffer_pos = pos
     end if
     call mpp_clock_end(nest_send_clock)
  end do ! end do list = 0,nlist-1

  !unpack buffer
  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_clock_end(nest_wait_clock)

  buffer_pos = buffer_recv_size

  call mpp_clock_begin(nest_unpk_clock)
  do m = update%nrecv, 1, -1
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle

     pos = buffer_pos
     do n = overPtr%count, 1, -1
        dir = overPtr%dir(n)
        if( recv(dir) ) then
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           select case (dir)
           case ( 1 ) ! east
              do l=1,l_size  ! loop over number of fields
                 ptr_ebuffer = eb_addrs(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          ebuffer(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
              end do
           case ( 3 ) ! south
              do l=1,l_size  ! loop over number of fields
                 ptr_sbuffer = sb_addrs(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          sbuffer(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
              end do
           case ( 5 ) ! west
              do l=1,l_size  ! loop over number of fields
                 ptr_wbuffer = wb_addrs(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          wbuffer(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
              end do
           case ( 7 ) ! north
              do l=1,l_size  ! loop over number of fields
                 ptr_nbuffer = nb_addrs(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          nbuffer(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
              end do
           end select
        endif
     end do ! do n = 1, overPtr%count
  end do
  call mpp_clock_end(nest_unpk_clock)

  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self( )
  call mpp_clock_end(nest_wait_clock)
      return

end subroutine MPP_DO_UPDATE_NEST_FINE_3D_

#ifdef VECTOR_FIELD_
subroutine MPP_DO_UPDATE_NEST_FINE_3D_V_(f_addrsx, f_addrsy, nest_domain, update_x, update_y, d_type, ke, wb_addrsx, wb_addrsy, &
                                   eb_addrsx, eb_addrsy, sb_addrsx, sb_addrsy, nb_addrsx, nb_addrsy, flags)
  integer(LONG_KIND),         intent(in) :: f_addrsx(:), f_addrsy(:)
  type(nest_level_type),      intent(in) :: nest_domain
  type(nestSpec),             intent(in) :: update_x, update_y
  MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
  integer,                    intent(in) :: ke
  integer(LONG_KIND),         intent(in) :: wb_addrsx(:), wb_addrsy(:)
  integer(LONG_KIND),         intent(in) :: eb_addrsx(:), eb_addrsy(:)
  integer(LONG_KIND),         intent(in) :: sb_addrsx(:), sb_addrsy(:)
  integer(LONG_KIND),         intent(in) :: nb_addrsx(:), nb_addrsy(:)
  integer,                    intent(in) :: flags

  character(len=8)            :: text
  logical   :: send(8), recv(8)
  integer   :: from_pe, to_pe, dir
  integer   :: m, n, l, i, j, k
  integer   :: is, ie, js, je, l_size
  integer   :: buffer_pos, msgsize
  integer   :: buffer_recv_size, pos
  integer   :: nrecv, nsend, ind_x, ind_y
  integer   :: ind_send_x(update_x%nsend+update_y%nsend), ind_send_y(update_x%nsend+update_y%nsend)
  integer   :: ind_recv_x(update_x%nrecv+update_y%nrecv), ind_recv_y(update_x%nrecv+update_y%nrecv)
  integer   :: from_pelist(update_x%nrecv+update_y%nrecv), to_pelist(update_x%nsend+update_y%nsend)
  integer   :: start_pos_recv(update_x%nrecv+update_y%nrecv), start_pos_send(update_x%nsend+update_y%nsend)

  MPP_TYPE_ :: fieldx(update_x%xbegin:update_x%xend, update_x%ybegin:update_x%yend,ke)
  MPP_TYPE_ :: fieldy(update_y%xbegin:update_y%xend, update_y%ybegin:update_y%yend,ke)
  MPP_TYPE_ :: wbufferx(update_x%west%is_you :update_x%west%ie_you,  update_x%west%js_you :update_x%west%je_you, ke)
  MPP_TYPE_ :: ebufferx(update_x%east%is_you :update_x%east%ie_you,  update_x%east%js_you :update_x%east%je_you, ke)
  MPP_TYPE_ :: sbufferx(update_x%south%is_you:update_x%south%ie_you, update_x%south%js_you:update_x%south%je_you,ke)
  MPP_TYPE_ :: nbufferx(update_x%north%is_you:update_x%north%ie_you, update_x%north%js_you:update_x%north%je_you,ke)
  MPP_TYPE_ :: wbuffery(update_y%west%is_you :update_y%west%ie_you,  update_y%west%js_you :update_y%west%je_you, ke)
  MPP_TYPE_ :: ebuffery(update_y%east%is_you :update_y%east%ie_you,  update_y%east%js_you :update_y%east%je_you, ke)
  MPP_TYPE_ :: sbuffery(update_y%south%is_you:update_y%south%ie_you, update_y%south%js_you:update_y%south%je_you,ke)
  MPP_TYPE_ :: nbuffery(update_y%north%is_you:update_y%north%ie_you, update_y%north%js_you:update_y%north%je_you,ke)
  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))

  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)
  pointer(ptr_buffer, buffer )
  pointer(ptr_wbufferx, wbufferx)
  pointer(ptr_ebufferx, ebufferx)
  pointer(ptr_sbufferx, sbufferx)
  pointer(ptr_nbufferx, nbufferx)
  pointer(ptr_wbuffery, wbuffery)
  pointer(ptr_ebuffery, ebuffery)
  pointer(ptr_sbuffery, sbuffery)
  pointer(ptr_nbuffery, nbuffery)

  recv(1) = BTEST(flags,EAST)
  recv(3) = BTEST(flags,SOUTH)
  recv(5) = BTEST(flags,WEST)
  recv(7) = BTEST(flags,NORTH)

  send    = recv

  ptr_buffer = LOC(mpp_domains_stack)
  l_size = size(f_addrsx(:))

  nrecv = get_nest_vector_recv(nest_domain, update_x, update_y, ind_recv_x, ind_recv_y, start_pos_recv, from_pelist)
  nsend = get_nest_vector_send(nest_domain, update_x, update_y, ind_send_x, ind_send_y, start_pos_send, to_pelist)
  !--- pre-post receiving
  buffer_pos = 0
  call mpp_clock_begin(nest_recv_clock)
  do m = 1, nrecv
     msgsize = 0
     ind_x = ind_recv_x(m)
     ind_y = ind_recv_y(m)
     if(ind_x >= 0) then
        do n = 1, update_x%recv(ind_x)%count
           msgsize = msgsize + update_x%recv(ind_x)%msgsize(n)
        end do
     endif
     if(ind_y >= 0) then
        do n = 1, update_y%recv(ind_y)%count
           msgsize = msgsize + update_y%recv(ind_y)%msgsize(n)
        end do
     endif
     if( msgsize.GT.0 )then
     msgsize = msgsize*ke*l_size
        from_pe = from_pelist(m)

        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_FINE_3D_V: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do ! end do m = 1, update%nrecv
  call mpp_clock_end(nest_recv_clock)
  buffer_recv_size = buffer_pos

  !--- pack and send the data
  do m = 1, nsend
     ind_x = ind_send_x(m)
     ind_y = ind_send_y(m)

     call mpp_clock_begin(nest_pack_clock)
     pos = buffer_pos
     msgsize = 0

     if(ind_x >= 0) then
        do n = 1, update_x%send(ind_x)%count
           msgsize = msgsize + update_x%send(ind_x)%msgsize(n)
        end do
     endif
     if(ind_y >= 0) then
        do n = 1, update_y%send(ind_y)%count
           msgsize = msgsize + update_y%send(ind_y)%msgsize(n)
        end do
     endif

     if( msgsize.GT.0 )then
        msgsize = msgsize*ke*l_size
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_FINE_3D_V: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if

     if(ind_x>=0) then
        do n = 1, update_x%send(ind_x)%count
           dir = update_x%send(ind_x)%dir(n)
           if( send(dir) ) then
              is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
              js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
              select case(update_x%send(ind_x)%rotation(n))
              case(ZERO)
                 do l=1,l_size  ! loop over number of fields
                    ptr_fieldx = f_addrsx(l)
                    do k = 1,ke
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
                       ptr_fieldy = f_addrsy(l)
                       do k = 1,ke
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
                       ptr_fieldy = f_addrsy(l)
                       do k = 1,ke
                          do i = is, ie
                             do j = je, js, -1
                                pos = pos + 1
                                buffer(pos) = -fieldy(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 endif
              case(NINETY)
                 do l=1,l_size  ! loop over number of fields
                    ptr_fieldy = f_addrsy(l)
                    do k = 1,ke
                       do i = ie, is, -1
                          do j = js, je
                             pos = pos + 1
                             buffer(pos) = fieldy(i,j,k)
                          end do
                       end do
                    end do
                 end do
              case default
                 call mpp_error(FATAL, "MPP_DO_UPDATE_NEST_FINE_3D_V X: pack rotate must be ZERO, MINUS_NINETY or NINETY")
              end select
           endif
        end do ! do n = 1, update_x%send(ind_x)%count
     endif

     if(ind_y>=0) then
        do n = 1, update_y%send(ind_y)%count
           dir = update_y%send(ind_y)%dir(n)
           if( send(dir) ) then
              is = update_y%send(ind_y)%is(n); ie = update_y%send(ind_y)%ie(n)
              js = update_y%send(ind_y)%js(n); je = update_y%send(ind_y)%je(n)

              select case(update_y%send(ind_y)%rotation(n))
              case(ZERO)
                 do l=1,l_size  ! loop over number of fields
                    ptr_fieldy = f_addrsy(l)
                    do k = 1,ke
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
                    ptr_fieldx = f_addrsx(l)
                    do k = 1,ke
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
                       ptr_fieldx = f_addrsx(l)
                       do k = 1,ke
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
                       ptr_fieldx = f_addrsx(l)
                       do k = 1,ke
                          do i = ie, is, -1
                             do j = js, je
                                pos = pos + 1
                                buffer(pos) = -fieldx(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 endif
              case default
                 call mpp_error(FATAL, "MPP_DO_UPDATE_NEST_FINE_3D_V Y: pack rotate must be ZERO, MINUS_NINETY or NINETY")
              end select
           endif
        end do ! do n = 1, update_x%send(ind_x)%count
     endif

     call mpp_clock_end(nest_pack_clock)
     call mpp_clock_begin(nest_send_clock)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then
        to_pe = to_pelist(m)
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1 )
        buffer_pos = pos
     end if
     call mpp_clock_end(nest_send_clock)
  end do ! end do list = 0,nlist-1

  !unpack buffer
  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_clock_end(nest_wait_clock)

  buffer_pos = buffer_recv_size

  call mpp_clock_begin(nest_unpk_clock)
  do m = nrecv, 1, -1
     ind_x = ind_recv_x(m)
     ind_y = ind_recv_y(m)

     pos = buffer_pos
     if(ind_y>=0) then
        do n = update_y%recv(ind_y)%count, 1, -1
           dir = update_y%recv(ind_y)%dir(n)
           if( recv(dir) ) then
              is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
              js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)
              msgsize = (ie-is+1)*(je-js+1)*ke*l_size
              pos = buffer_pos - msgsize
              buffer_pos = pos
              select case (dir)
              case ( 1 ) ! east
                 do l=1,l_size  ! loop over number of fields
                    ptr_ebuffery = eb_addrsy(l)
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             ebuffery(i,j,k) = buffer(pos)
                          end do
                       end do
                    end do
                 end do
              case ( 3 ) ! south
                 do l=1,l_size  ! loop over number of fields
                    ptr_sbuffery = sb_addrsy(l)
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             sbuffery(i,j,k) = buffer(pos)
                          end do
                       end do
                    end do
                 end do
              case ( 5 ) ! west
                 do l=1,l_size  ! loop over number of fields
                    ptr_wbuffery = wb_addrsy(l)
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             wbuffery(i,j,k) = buffer(pos)
                          end do
                       end do
                    end do
                 end do
              case ( 7 ) ! north
                 do l=1,l_size  ! loop over number of fields
                    ptr_nbuffery = nb_addrsy(l)
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             nbuffery(i,j,k) = buffer(pos)
                          end do
                       end do
                    end do
                 end do
              end select
           endif
        end do ! do n = update_y%recv(ind_y)%count, 1, -1
     endif

     if(ind_x>=0) then
        do n = update_x%recv(ind_x)%count, 1, -1
           dir = update_x%recv(ind_x)%dir(n)
           if( recv(dir) ) then
              is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
              js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n)
              msgsize = (ie-is+1)*(je-js+1)*ke*l_size
              pos = buffer_pos - msgsize
              buffer_pos = pos
              select case (dir)
              case ( 1 ) ! east
                 do l=1,l_size  ! loop over number of fields
                    ptr_ebufferx = eb_addrsx(l)
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             ebufferx(i,j,k) = buffer(pos)
                          end do
                       end do
                    end do
                 end do
              case ( 3 ) ! south
                 do l=1,l_size  ! loop over number of fields
                    ptr_sbufferx = sb_addrsx(l)
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             sbufferx(i,j,k) = buffer(pos)
                          end do
                       end do
                    end do
                 end do
              case ( 5 ) ! west
                 do l=1,l_size  ! loop over number of fields
                    ptr_wbufferx = wb_addrsx(l)
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             wbufferx(i,j,k) = buffer(pos)
                          end do
                       end do
                    end do
                 end do
              case ( 7 ) ! north
                 do l=1,l_size  ! loop over number of fields
                    ptr_nbufferx = nb_addrsx(l)
                    do k = 1,ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             nbufferx(i,j,k) = buffer(pos)
                          end do
                       end do
                    end do
                 end do
              end select
           endif
        end do ! do n = update_x%recv(ind_x)%count, 1, -1
     endif
  end do
  call mpp_clock_end(nest_unpk_clock)

  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self( )
  call mpp_clock_end(nest_wait_clock)
  return

end subroutine MPP_DO_UPDATE_NEST_FINE_3D_V_

#endif   !VECTOR_FIELD_


!###############################################################################
subroutine MPP_DO_UPDATE_NEST_COARSE_3D_(f_addrs_in, f_addrs_out, nest_domain, update, d_type, ke)
  integer(LONG_KIND),         intent(in) :: f_addrs_in(:)
  integer(LONG_KIND),         intent(in) :: f_addrs_out(:)
  type(nest_domain_type),     intent(in) :: nest_domain
  type(nestSpec),             intent(in) :: update
  MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
  integer,                    intent(in) :: ke

  character(len=8)            :: text
  type(overlap_type), pointer :: overPtr => NULL()
  integer   :: from_pe, to_pe
  integer   :: m, n, l, i, j, k
  integer   :: is, ie, js, je, l_size
  integer   :: buffer_pos, msgsize
  integer   :: buffer_recv_size, pos
  MPP_TYPE_ :: field_in(update%xbegin_c:update%xend_c, update%ybegin_c:update%yend_c,ke)
  MPP_TYPE_ :: field_out(update%xbegin:update%xend, update%ybegin:update%yend,ke)
  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))

  pointer(ptr_field_in, field_in)
  pointer(ptr_field_out, field_out)
  pointer(ptr_buffer, buffer )

  ptr_buffer = LOC(mpp_domains_stack)
  l_size = size(f_addrs_in(:))

  !--- pre-post receiving
  buffer_pos = 0
  call mpp_clock_begin(nest_recv_clock)
  do m = 1, update%nrecv
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle
     msgsize = 0
     do n = 1, overPtr%count
        is = overPtr%is(n); ie = overPtr%ie(n)
        js = overPtr%js(n); je = overPtr%je(n)
        msgsize = msgsize + (ie-is+1)*(je-js+1)
     end do

     msgsize = msgsize*ke*l_size
     if( msgsize.GT.0 )then
        from_pe = overPtr%pe
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_COARSE_3D_: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_2 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do ! end do m = 1, update%nrecv
  call mpp_clock_end(nest_recv_clock)
  buffer_recv_size = buffer_pos

  !--- pack and send the data
  do m = 1, update%nsend
     overPtr => update%send(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(nest_pack_clock)
     pos = buffer_pos
     msgsize = 0
     do n = 1, overPtr%count
        msgsize = msgsize + overPtr%msgsize(n)
     enddo
     if( msgsize.GT.0 )then
        msgsize = msgsize*ke*l_size
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_COARSE_3D_: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if

     do n = 1, overPtr%count
        is = overPtr%is(n); ie = overPtr%ie(n)
        js = overPtr%js(n); je = overPtr%je(n)

        select case(overPtr%rotation(n))
        case(ZERO)
        do l=1,l_size  ! loop over number of fields
              ptr_field_in = f_addrs_in(l)
           do k = 1,ke
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                       buffer(pos) = field_in(i,j,k)
                 end do
              end do
           end do
        end do
        case(MINUS_NINETY)
           do l=1,l_size  ! loop over number of fields
              ptr_field_in = f_addrs_in(l)
              do k = 1,ke
                 do i = is, ie
                    do j = je, js, -1
                       pos = pos + 1
                       buffer(pos) = field_in(i,j,k)
                    end do
                 end do
              end do
           end do
        case(NINETY)
           do l=1,l_size  ! loop over number of fields
              ptr_field_in = f_addrs_in(l)
              do k = 1,ke
                 do i = ie, is, -1
                    do j = js, je
                       pos = pos + 1
                       buffer(pos) = field_in(i,j,k)
                    end do
                 end do
              end do
           end do
        case default
           call mpp_error(FATAL, "MPP_DO_UPDATE_NEST_COARSE_3D: pack rotate must be ZERO, MINUS_NINETY or NINETY")
        end select
     end do ! do n = 1, overPtr%count

     call mpp_clock_end(nest_pack_clock)
     call mpp_clock_begin(nest_send_clock)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then
        to_pe = overPtr%pe
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_2 )
        buffer_pos = pos
     end if
     call mpp_clock_end(nest_send_clock)
  end do ! end do list = 0,nlist-1

  !unpack buffer
  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_clock_end(nest_wait_clock)

  buffer_pos = buffer_recv_size

  call mpp_clock_begin(nest_unpk_clock)
  do m = update%nrecv, 1, -1
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle

     pos = buffer_pos
     do n = overPtr%count, 1, -1
        is = overPtr%is(n); ie = overPtr%ie(n)
        js = overPtr%js(n); je = overPtr%je(n)
        msgsize = (ie-is+1)*(je-js+1)*ke*l_size
        pos = buffer_pos - msgsize
        buffer_pos = pos
        do l=1,l_size  ! loop over number of fields
           ptr_field_out = f_addrs_out(l)
           do k = 1,ke
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                    field_out(i,j,k) = buffer(pos)
                 end do
              end do
           end do
        end do
     end do ! do n = 1, overPtr%count
  end do
  call mpp_clock_end(nest_unpk_clock)

  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self( )
  call mpp_clock_end(nest_wait_clock)
      return

end subroutine MPP_DO_UPDATE_NEST_COARSE_3D_


#ifdef VECTOR_FIELD_
!###############################################################################
subroutine MPP_DO_UPDATE_NEST_COARSE_3D_V_(f_addrsx_in, f_addrsy_in, f_addrsx_out, f_addrsy_out, &
                                           nest_domain, nest, update_x, update_y, d_type, ke, flags)
  integer(LONG_KIND),         intent(in) :: f_addrsx_in(:), f_addrsy_in(:)
  integer(LONG_KIND),         intent(in) :: f_addrsx_out(:), f_addrsy_out(:)
  type(nest_domain_type),     intent(in) :: nest_domain
  type(nest_level_type),      intent(in) :: nest
  type(nestSpec),             intent(in) :: update_x, update_y
  MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
  integer,                    intent(in) :: ke
  integer,                    intent(in) :: flags

  character(len=8)            :: text
  type(overlap_type), pointer :: overPtr => NULL()
  integer   :: from_pe, to_pe
  integer   :: m, n, l, i, j, k
  integer   :: is, ie, js, je, l_size
  integer   :: buffer_pos, msgsize
  integer   :: buffer_recv_size, pos
  integer   :: nrecv, nsend, ind_x, ind_y
  integer   :: ind_send_x(update_x%nsend+update_y%nsend), ind_send_y(update_x%nsend+update_y%nsend)
  integer   :: ind_recv_x(update_x%nrecv+update_y%nrecv), ind_recv_y(update_x%nrecv+update_y%nrecv)
  integer   :: from_pelist(update_x%nrecv+update_y%nrecv), to_pelist(update_x%nsend+update_y%nsend)
  integer   :: start_pos_recv(update_x%nrecv+update_y%nrecv), start_pos_send(update_x%nsend+update_y%nsend)

  MPP_TYPE_ :: fieldx_in(update_x%xbegin_c:update_x%xend_c, update_x%ybegin_c:update_x%yend_c,ke)
  MPP_TYPE_ :: fieldy_in(update_y%xbegin_c:update_y%xend_c, update_y%ybegin_c:update_y%yend_c,ke)
  MPP_TYPE_ :: fieldx_out(update_x%xbegin:update_x%xend, update_x%ybegin:update_x%yend,ke)
  MPP_TYPE_ :: fieldy_out(update_y%xbegin:update_y%xend, update_y%ybegin:update_y%yend,ke)
  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))

  pointer(ptr_fieldx_in, fieldx_in)
  pointer(ptr_fieldy_in, fieldy_in)
  pointer(ptr_fieldx_out, fieldx_out)
  pointer(ptr_fieldy_out, fieldy_out)
  pointer(ptr_buffer, buffer )

  ptr_buffer = LOC(mpp_domains_stack)
  l_size = size(f_addrsx_in(:))

  nrecv = get_nest_vector_recv(nest, update_x, update_y, ind_recv_x, ind_recv_y, start_pos_recv, from_pelist)
  nsend = get_nest_vector_send(nest, update_x, update_y, ind_send_x, ind_send_y, start_pos_send, to_pelist)
  !--- pre-post receiving
  buffer_pos = 0
  call mpp_clock_begin(nest_recv_clock)
  do m = 1, nrecv
     msgsize = 0
     ind_x = ind_recv_x(m)
     ind_y = ind_recv_y(m)
     if(ind_x >= 0) then
        do n = 1, update_x%recv(ind_x)%count
           msgsize = msgsize + update_x%recv(ind_x)%msgsize(n)
        end do
     endif
     if(ind_y >= 0) then
        do n = 1, update_y%recv(ind_y)%count
           msgsize = msgsize + update_y%recv(ind_y)%msgsize(n)
        end do
     endif

     if( msgsize.GT.0 )then
        msgsize = msgsize*ke*l_size
        from_pe = from_pelist(m)
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_COARSE_3D_: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_2 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do ! end do m = 1, nrecv
  call mpp_clock_end(nest_recv_clock)
  buffer_recv_size = buffer_pos

  !--- pack and send the data
  do m = 1, nsend
     ind_x = ind_send_x(m)
     ind_y = ind_send_y(m)

     call mpp_clock_begin(nest_pack_clock)
     pos = buffer_pos
     msgsize = 0
     if(ind_x >= 0) then
        do n = 1, update_x%send(ind_x)%count
           msgsize = msgsize + update_x%send(ind_x)%msgsize(n)
        end do
     endif
     if(ind_y >= 0) then
        do n = 1, update_y%send(ind_y)%count
           msgsize = msgsize + update_y%send(ind_y)%msgsize(n)
        end do
     endif

     if( msgsize.GT.0 )then
        msgsize = msgsize*ke*l_size
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_COARSE_3D_: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if

     if(ind_x>=0) then
        do n = 1, update_x%send(ind_x)%count
           is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
           js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)

           select case(update_x%send(ind_x)%rotation(n) )
           case(ZERO)
              do l=1,l_size  ! loop over number of fields
                 ptr_fieldx_in = f_addrsx_in(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = fieldx_in(i,j,k)
                       end do
                    end do
                 end do
              end do
           case(MINUS_NINETY)
              if( BTEST(flags,SCALAR_BIT) ) then

                 do l=1,l_size  ! loop over number of fields
                    ptr_fieldy_in = f_addrsy_in(l)
                    do k = 1,ke
                       do i = is, ie
                          do j = je, js, -1
                             pos = pos + 1
                             buffer(pos) = fieldy_in(i,j,k)
                          end do
                       end do
                    end do
                 end do
              else
                 do l=1,l_size  ! loop over number of fields
                    ptr_fieldy_in = f_addrsy_in(l)
                    do k = 1,ke
                       do i = is, ie
                          do j = je, js, -1
                             pos = pos + 1
                             buffer(pos) = -fieldy_in(i,j,k)
                          end do
                       end do
                    end do
                 end do
              endif
           case(NINETY)
              do l=1,l_size  ! loop over number of fields
                 ptr_fieldy_in = f_addrsy_in(l)
                 do k = 1,ke
                    do i = ie, is, -1
                       do j = js, je
                          pos = pos + 1
                          buffer(pos) = fieldy_in(i,j,k)
                       end do
                    end do
                 end do
              end do
           case default
              call mpp_error(FATAL, "MPP_DO_UPDATE_NEST_COARSE_3D: pack rotate must be ZERO, MINUS_NINETY or NINETY")
           end select
        end do ! do n = 1, update_x%send(ind_x)%count
     end if

     if(ind_y>=0) then
        do n = 1, update_y%send(ind_y)%count
           is = update_y%send(ind_y)%is(n); ie = update_y%send(ind_y)%ie(n)
           js = update_y%send(ind_y)%js(n); je = update_y%send(ind_y)%je(n)

           select case(update_y%send(ind_y)%rotation(n) )
           case(ZERO)
              do l=1,l_size  ! loop over number of fields
                 ptr_fieldy_in = f_addrsy_in(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = fieldy_in(i,j,k)
                       end do
                    end do
                 end do
              end do
           case(MINUS_NINETY)
              do l=1,l_size  ! loop over number of fields
                 ptr_fieldx_in = f_addrsx_in(l)
                 do k = 1,ke
                    do i = is, ie
                       do j = je, js, -1
                          pos = pos + 1
                          buffer(pos) = fieldx_in(i,j,k)
                       end do
                    end do
                 end do
              end do
           case(NINETY)
              if( BTEST(flags,SCALAR_BIT) ) then
                 do l=1,l_size  ! loop over number of fields
                    ptr_fieldx_in = f_addrsx_in(l)
                    do k = 1,ke
                       do i = ie, is, -1
                          do j = js, je
                             pos = pos + 1
                             buffer(pos) = fieldx_in(i,j,k)
                          end do
                       end do
                    end do
                 end do
              else
                 do l=1,l_size  ! loop over number of fields
                    ptr_fieldx_in = f_addrsx_in(l)
                    do k = 1,ke
                       do i = ie, is, -1
                          do j = js, je
                             pos = pos + 1
                             buffer(pos) = -fieldx_in(i,j,k)
                          end do
                       end do
                    end do
                 end do
              end if
           case default
              call mpp_error(FATAL, "MPP_DO_UPDATE_NEST_COARSE_3D: pack rotate must be ZERO, MINUS_NINETY or NINETY")
           end select
        end do ! do n = 1, update_y%send(ind_y)%count
     end if

     call mpp_clock_end(nest_pack_clock)
     call mpp_clock_begin(nest_send_clock)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then
        to_pe = to_pelist(m)
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_2 )
        buffer_pos = pos
     end if
     call mpp_clock_end(nest_send_clock)
  end do ! end do list = 0,nsend

  !unpack buffer
  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_clock_end(nest_wait_clock)

  buffer_pos = buffer_recv_size

  call mpp_clock_begin(nest_unpk_clock)
  do m = nrecv, 1, -1
     ind_x = ind_recv_x(m)
     ind_y = ind_recv_y(m)
     pos = buffer_pos
     if(ind_y>=0) then
        do n = update_y%recv(ind_y)%count, 1, -1
           is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
           js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)

           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           do l=1,l_size  ! loop over number of fields
              ptr_fieldy_out = f_addrsy_out(l)
              do k = 1,ke
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       fieldy_out(i,j,k) = buffer(pos)
                    end do
                 end do
              end do
           end do
        end do
     end if
     if(ind_x>=0) then
        do n = update_x%recv(ind_x)%count, 1, -1
           is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
           js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n)

           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           do l=1,l_size  ! loop over number of fields
              ptr_fieldx_out = f_addrsx_out(l)
              do k = 1,ke
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       fieldx_out(i,j,k) = buffer(pos)
                    end do
                 end do
              end do
           end do
        end do
     end if
  end do
  call mpp_clock_end(nest_unpk_clock)

  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self( )
  call mpp_clock_end(nest_wait_clock)
  return

end subroutine MPP_DO_UPDATE_NEST_COARSE_3D_V_

#endif   !VECTOR_FIELD_
