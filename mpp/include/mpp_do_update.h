! -*-f90-*- 
    subroutine MPP_DO_UPDATE_3D_( f_addrs, domain, d_type, ke, b_addrs, b_size, flags, name)
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND), intent(in)         :: f_addrs(:,:)
      type(domain2D), intent(in)    :: domain
      MPP_TYPE_, intent(in)                  :: d_type  ! creates unique interface
      integer,   intent(in)                  :: ke
      integer(LONG_KIND), intent(in)         :: b_addrs(:,:)
      integer,            intent(in)         :: b_size
      integer, intent(in), optional :: flags
      character(len=*), intent(in), optional :: name

      MPP_TYPE_ :: field(domain%x(1)%memory%begin:domain%x(1)%memory%end, domain%y(1)%memory%begin:domain%y(1)%memory%end,ke)
      MPP_TYPE_ :: fillbuffer(b_size)
      pointer(ptr_field, field)
      pointer(ptr_buffer, fillbuffer)
      integer                    :: update_flags
      type(boundary),   pointer  :: check   => NULL()
      type(overlapSpec), pointer :: overPtr => NULL()      
      character(len=8)           :: text
      character(len=64)          :: field_name      

!equate to mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer( ptr, buffer )
      integer :: buffer_pos
     
!receive domains saved here for unpacking
!for non-blocking version, could be recomputed
      logical :: send(8), recv(8)
      integer :: to_pe, from_pe, list, pos, msgsize
      integer :: n, l_size, l, m, i, j, k, nlist
      integer :: is, ie, js, je, tMe, dir
      integer :: start, start1, start2, index, is1, ie1, js1, je1, ni, nj, total
      integer :: buffer_recv_size

      nlist = size(domain%list(:))
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
      ptr = LOC(mpp_domains_stack)
      l_size = size(f_addrs,1)

      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) )update_flags = flags

      recv(1) = BTEST(update_flags,EAST)
      recv(3) = BTEST(update_flags,SOUTH)
      recv(5) = BTEST(update_flags,WEST)
      recv(7) = BTEST(update_flags,NORTH)
      recv(2) = recv(1) .AND. recv(3)
      recv(4) = recv(3) .AND. recv(5)
      recv(6) = recv(5) .AND. recv(7)
      recv(8) = recv(7) .AND. recv(1)
      send    = recv

      !--- if debug_update_level is not NO_DEBUG, check the consistency on the bounds 
      !--- (domain is symmetry or folded north edge). North bound will be checked when north edge is folded.
      !--- when domain is symmetry, For data on T-cell, no check is needed; for data on E-cell, 
      !--- data on East and West boundary will be checked ; For data on N-cell, data on North and South 
      !--- boundary will be checked; For data on C-cell, data on West, East, South, North will be checked.
      !--- The check will be done in the following way: Western boundary data sent to Eastern boundary to check
      !--- and Southern boundary to check

      if( debug_update_level .NE. NO_CHECK ) then      
         if(present(name)) then
            field_name = name
         else
            field_name = "un-named"
         end if

         !--- pre-post recv the data 
         do list = 0,nlist-1
            m = mod( domain%pos+nlist-list, nlist )
            check=>domain%check%recv(m)
            msgsize = 0
            do n = 1, check%count
               is = check%is(n); ie = check%ie(n)
               js = check%js(n); je = check%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            msgsize = msgsize*ke*l_size

            if( msgsize.GT.0 )then
               from_pe = domain%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_DO_UPDATE: mpp_domains_stack overflow, '// &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
               end if
               call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE. )
               buffer_pos = buffer_pos + msgsize
            end if
         end do
         buffer_recv_size = buffer_pos

         !--- send the data
         do list = 0,nlist-1
            m = mod( domain%pos+list, nlist )
            check => domain%check%send(m)
            pos = buffer_pos
            do n = 1, check%count
               is = check%is(n); ie = check%ie(n)
               js = check%js(n); je = check%je(n)
               tMe = check%tileMe(n)
               select case( check%rotation(n) )
               case(ZERO)
                  do l = 1, l_size ! loop over number of fields
                     ptr_field = f_addrs(l, tMe)
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
                  do l = 1, l_size ! loop over number of fields
                     ptr_field = f_addrs(l, tMe)
                     do k = 1,ke  
                        do j = je, js, -1
                           do i = is, ie

                              pos = pos + 1
                              buffer(pos) = field(i,j,k)
                           end do
                        end do
                     end do
                  end do
               case(NINETY)
                  do l = 1, l_size ! loop over number of fields
                     ptr_field = f_addrs(l, tMe)
                     do k = 1,ke  
                        do j = js, je
                           do i = ie, is, -1

                              pos = pos + 1
                              buffer(pos) = field(i,j,k)
                           end do
                        end do
                     end do
                  end do
               case(ONE_HUNDRED_EIGHTY)
                  do l = 1, l_size ! loop over number of fields
                     ptr_field = f_addrs(l, tMe)
                     do k = 1,ke  
                        do j = je, js, -1
                           do i = ie, is, -1
                              pos = pos + 1
                              buffer(pos) = field(i,j,k)
                           end do
                        end do
                     end do
                  end do
               end select
            end do
            msgsize = pos - buffer_pos
            if( msgsize.GT.0 )then
               to_pe = domain%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_DO_UPDATE: mpp_domains_stack overflow, ' // &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
               end if
               call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
               buffer_pos = pos
            end if
         end do ! end do list = 0,nlist-1

         call mpp_sync_self(check=EVENT_RECV) ! To ensure recv is completed.
         buffer_pos = buffer_recv_size
         !--- compare the data in reverse order
         CHECK_LOOP: do list = nlist-1,0,-1
            m = mod( domain%pos+nlist-list, nlist )
            check=>domain%check%recv(m)
            do n = check%count, 1, -1
               is = check%is(n); ie = check%ie(n)
               js = check%js(n); je = check%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               tMe = check%tileMe(n)
               do l=1, l_size  ! loop over number of fields
                  ptr_field = f_addrs(l, tMe)
                  do k = 1,ke
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           if( field(i,j,k) .NE. buffer(pos) ) then
                              print*,"Error from MPP_DO_UPDATE on pe = ", mpp_pe(), ": field ", &
                                   trim(field_name), " at point (", i, ",", j, ",", k, ") = ", field(i,j,k), &
                                   " does not equal to the value = ", buffer(pos), " on pe ", domain%list(m)%pe
                              call mpp_error(debug_update_level, "MPP_DO_UPDATE: mismatch on the boundary for symmetry point")
                              exit CHECK_LOOP
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do CHECK_LOOP ! end do list = nlist-1,0,-1
         call mpp_sync_self()
      end if ! end if( check ... )

      buffer_pos = 0

      !recv
      do list = 0,nlist-1
         m = mod( domain%pos+nlist-list, nlist )
         overPtr => domain%update_recv(m)
         if( overPtr%count == 0 )cycle
         call mpp_clock_begin(recv_clock)
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
            from_pe = domain%list(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_UPDATE: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
            end if
            call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE. )
            buffer_pos = buffer_pos + msgsize
         end if
         call mpp_clock_end(recv_clock)
      end do ! end do list = 0,nlist-1
      buffer_recv_size = buffer_pos

      ! send
      do list = 0,nlist-1
         m = mod( domain%pos+list, nlist )
         overPtr => domain%update_send(m)
         if( overPtr%count == 0 )cycle
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos
         do n = 1, overPtr%count
            dir = overPtr%dir(n)
            if( send(dir) ) then
               tMe = overPtr%tileMe(n)
               is = overPtr%is(n); ie = overPtr%ie(n)
               js = overPtr%js(n); je = overPtr%je(n)
               if( overptr%is_refined(n) ) then
                  do l=1,l_size  ! loop over number of fields
                     ptr_field = f_addrs(l, tMe)
                     do k = 1,ke  
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = field(i,j,k)
                           end do
                        end do
                     end do
                  end do
               else
                  select case( overPtr%rotation(n) )
                  case(ZERO)
                     do l=1,l_size  ! loop over number of fields
                        ptr_field = f_addrs(l, tMe)
                        do k = 1,ke  
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 buffer(pos) = field(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  case( MINUS_NINETY ) 
                     do l=1,l_size  ! loop over number of fields
                        ptr_field = f_addrs(l, tMe)
                        do k = 1,ke  
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
                        do k = 1,ke  
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
                        do k = 1,ke  
                           do j = je, js, -1
                              do i = ie, is, -1
                                 pos = pos + 1
                                 buffer(pos) = field(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  end select
               end if
            endif
         end do ! do n = 1, overPtr%count

         call mpp_clock_end(pack_clock)
         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = domain%list(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_UPDATE: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
            buffer_pos = pos
         end if
         call mpp_clock_end(send_clock)
      end do ! end do ist = 0,nlist-1

      !unpack recv
      !unpack halos in reverse order
!     ptr_rfield = f_addrs(1)
      call mpp_clock_begin(wait_clock)
      call mpp_sync_self(check=EVENT_RECV)
      call mpp_clock_end(wait_clock)
      buffer_pos = buffer_recv_size      

      do list = nlist-1,0,-1
         m = mod( domain%pos+nlist-list, nlist )
         overPtr => domain%update_recv(m)
         if( overPtr%count == 0 )cycle
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         do n = overPtr%count, 1, -1
            dir = overPtr%dir(n)
            if( recv(dir) ) then
               tMe = overPtr%tileMe(n)
               is = overPtr%is(n); ie = overPtr%ie(n)
               js = overPtr%js(n); je = overPtr%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               if(OverPtr%is_refined(n)) then
                  index = overPtr%index(n)
                  is1 = domain%rSpec(tMe)%isNbr(index); ie1 = domain%rSpec(tMe)%ieNbr(index)
                  js1 = domain%rSpec(tMe)%jsNbr(index); je1 = domain%rSpec(tMe)%jeNbr(index)
                  ni = ie1 - is1 + 1
                  nj = je1 - js1 + 1
                  total = ni*nj
                  start = (domain%rSpec(tMe)%start(index)-1)*ke
                  if(start+total*ke>b_size ) call mpp_error(FATAL, &
                       "MPP_DO_UPDATE: b_size is less than the size of the data to be filled.")
                  msgsize = ie - is + 1
                  do l=1, l_size  ! loop over number of fields
                     ptr_buffer = b_addrs(l, tMe)
                     start1 = start + (js-js1)*ni + is - is1
                     do k = 1, ke
                        start2 = start1
                        do j = js, je
                           fillbuffer(start2+1:start2+msgsize) = buffer(pos+1:pos+msgsize)
                           start2 = start2 + ni
                           pos   = pos + msgsize
                        end do
                        start1 = start1 + total
                     end do
                  end do
               else
                  do l=1,l_size  ! loop over number of fields
                     ptr_field = f_addrs(l, tMe)
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              field(i,j,k) = buffer(pos)
                           end do
                        end do
                     end do
                  end do
               endif
            end if
         end do ! do n = 1, overPtr%count
         call mpp_clock_end(unpk_clock)
      end do

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)
      return
    end subroutine MPP_DO_UPDATE_3D_
