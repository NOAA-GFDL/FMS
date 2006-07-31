! -*-f90-*- 
    subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_( f_addrs, domain, d_type, ke, flags, name)
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND), intent(in)         :: f_addrs(:)
      type(domain2d),     intent(in)         :: domain
      MPP_TYPE_, intent(in)                  :: d_type  ! creates unique interface
      integer,   intent(in)                  :: ke
      integer, intent(in), optional          :: flags 
      character(len=*), intent(in), optional :: name

      MPP_TYPE_ :: field(domain%x(1)%memory%begin:domain%x(1)%memory%end, domain%y(1)%memory%begin:domain%y(1)%memory%end,ke)
      pointer(ptr_field, field)
      integer                    :: update_flags
      type(checkbound),  pointer :: bound   => NULL()
      type(overlapSpec), pointer :: overPtr => NULL()
      character(len=8)           :: text
      character(len=64)          :: field_name      

!equate to mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer(ptr,buffer )
      integer :: buffer_pos
!      real(kind(buffer)) :: rbuffer(size(mpp_domains_stack(:)))
!      pointer(r_ptr,rbuffer)

!receive domains saved here for unpacking
!for non-blocking version, could be recomputed
      logical :: send(8), recv(8)
      integer :: to_pe, from_pe, list, pos, msgsize
      integer :: n, l_size, l, m, i, j, k, nlist, nrecv, nsend
      integer :: is, ie, js, je, ntileMe, ntileNbr, ndir, tMe, tNbr, dir

      nlist = size(domain%list(:))
      ntileMe = size(domain%x(:))
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
      ptr = LOC(mpp_domains_stack)
      l_size = size(f_addrs(:))

      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) ) update_flags = flags

      recv(1) = BTEST(update_flags,EAST)
      recv(3) = BTEST(update_flags,SOUTH)
      recv(5) = BTEST(update_flags,WEST)
      recv(7) = BTEST(update_flags,NORTH)
      recv(2) = recv(1) .AND. recv(3)
      recv(4) = recv(3) .AND. recv(5)
      recv(6) = recv(5) .AND. recv(7)
      recv(8) = recv(7) .AND. recv(1)
      send    = recv

      !send
      do list = 0,nlist-1
         m = mod( domain%pos+list, nlist )
         if( .NOT. domain%list(m)%overlap )cycle
         ntileNbr = size(domain%list(m)%x(:) )
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos
         do dir = 1,8
            if( send(dir) ) then
               do l=1,l_size  ! loop over number of fields
                  tMe = min(l, ntileMe)
                  ptr_field = f_addrs(l)
                  do tNbr = 1, ntileNbr
                     overPtr => domain%list(m)%send(tNbr,tMe,dir)
                     do n = 1, 3
                        if( overPtr%overlap(n) ) then
                           is = overPtr%is(n); ie = overPtr%ie(n)
                           js = overPtr%js(n); je = overPtr%je(n)
                           do k = 1,ke
                              do j = js, je
                                 do i = is, ie
                                    pos = pos + 1
                                    buffer(pos) = field(i,j,k)
                                 end do
                              end do
                           end do
                        end if
                     end do ! end do n = 1, 3
                     nsend = overPtr%n
                     if(nsend>0) then
                        do n = 1, nsend
                           i = overPtr%i(n); j = overPtr%j(n)
                           do k = 1,ke  
                              pos = pos + 1
                              buffer(pos) = field(i,j,k)
                           end do
                        end do
                     end if
                  end do ! end do tNbr = 1, ntileNbr
               end do  ! do l=1,l_size
            end if ! end if ( send(dir) )
         end do ! end do dir = 1, 8

         call mpp_clock_end(pack_clock)
         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = domain%list(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
            buffer_pos = pos
         end if
         call mpp_clock_end(send_clock)
      end do ! end do ist = 0,nlist-1
             
!recv
      do list = 0,nlist-1
         m = mod( domain%pos+nlist-list, nlist )
         if( .NOT. domain%list(m)%overlap )cycle
         ntileNbr = size(domain%list(m)%x(:))
         call mpp_clock_begin(recv_clock)
         msgsize = 0
         do dir = 1, 8  ! loop over 8 direction
            if(recv(dir)) then
               do tNbr = 1, ntileNbr
                  do l = 1, l_size  ! loop over number of fields
                     tMe = min(l, ntileMe)
                     overPtr => domain%list(m)%recv(tMe,tNbr,dir)
                     do n = 1, 3
                        if( overPtr%overlap(n) ) then
                           is = overPtr%is(n); ie = overPtr%ie(n)
                           js = overPtr%js(n); je = overPtr%je(n)
                           msgsize = msgsize + (ie-is+1)*(je-js+1)
                        end if
                     end do
                     msgsize = msgsize + overPtr%n
                  end do
               end do
            end if
         end do

         msgsize = msgsize*ke
         if( msgsize.GT.0 )then
            from_pe = domain%list(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
            end if
            call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
            buffer_pos = buffer_pos + msgsize
         end if
         call mpp_clock_end(recv_clock)
      end do ! end do list = 0,nlist-1
             
!unpack recv
!unpack halos in reverse order
!     ptr_rfield = f_addrs(1)

      do list = nlist-1,0,-1
         m = mod( domain%pos+nlist-list, nlist )
         if( .NOT.domain%list(m)%overlap)cycle
         ntileNbr = size(domain%list(m)%x(:))         
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         do dir=8,1,-1
            if( recv(dir) ) then
               do tNbr = ntileNbr, 1, -1
                  do l=l_size,1,-1  ! loop over number of fields
                     tMe = min(l, ntileMe)
                     ptr_field = f_addrs(l)
                     overPtr => domain%list(m)%recv(tMe,tNbr,dir)
                     nrecv = overPtr%n
                     if( nrecv > 0 )then  ! direction
                        msgsize = nrecv*ke
                        pos = buffer_pos - msgsize
                        buffer_pos = pos
                        do n = 1, nrecv
                           i = overPtr%i(n)
                           j = overPtr%j(n)
                           do k = 1,ke
                              pos = pos + 1
                              field(i,j,k) = buffer(pos)
                           end do
                        end do
                     end if
                     do n = 3, 1, -1 
                        if( overPtr%overlap(n) )then  ! direction
                           is = overPtr%is(n); ie = overPtr%ie(n)
                           js = overPtr%js(n); je = overPtr%je(n)
                           msgsize = (ie-is+1)*(je-js+1)*ke
                           pos = buffer_pos - msgsize
                           buffer_pos = pos
                           select case ( overPtr%rotation )
                           case ( ZERO )
                              do k = 1,ke
                                 do j = js, je
                                    do i = is, ie
                                       pos = pos + 1
                                       field(i,j,k) = buffer(pos)
                                    end do
                                 end do
                              end do
                           case ( MINUS_NINETY )
                              do k = 1,ke
                                 do i = ie, is, -1
                                    do j = js, je
                                       pos = pos + 1
                                       field(i,j,k) = buffer(pos)
                                    end do
                                 end do
                              end do
                           case ( NINETY )
                              do k = 1,ke
                                 do i = is, ie
                                    do j = je, js, -1
                                       pos = pos + 1
                                       field(i,j,k) = buffer(pos)
                                    end do
                                 end do
                              end do
                           case ( ONE_HUNDRED_EIGHTY )
                              do k = 1,ke
                                 do j = je, js, -1
                                    do i = ie, is, -1
                                       pos = pos + 1
                                       field(i,j,k) = buffer(pos)
                                    end do
                                 end do
                              end do
                           end select
                        end if
                     end do ! end do n = 3, 1, -1 
                  end do  ! do l=l_size,1,-1
               end do ! end do tNbr = ntileNbr, 1, -1
            end if ! end if( recv(dir) )
         end do  ! do dir=8,1,-1

         call mpp_clock_end(unpk_clock)
      end do

      !--- if debug is true and domain is symmetry, check the consistency on the bounds between tiles.
      !--- For data on T-cell, no check is needed; for data on E-cell, data on East and West boundary
      !--- will be checked; For data on N-cell, data on North and South boundary will be checked;
      !--- For data on C-cell, data on West, East, South, North will be checked.
      !--- The check will be done in the following way: Western boundary data sent to Eastern boundary to check
      !--- and Southern boundary to check

      if( debug .AND. domain%symmetry .AND. domain%position .NE. CENTER ) then      
         if(present(name)) then
            field_name = name
         else
            field_name = "un-named"
         end if

         bound => domain%bound
         ndir = size(bound%list(0)%send,3)
         !--- send the data
         do list = 0,nlist-1
            m = mod( domain%pos+list, nlist )
            if(.NOT. bound%list(m)%overlap) cycle
            ntileNbr = size(domain%list(m)%x(:) )
            pos = buffer_pos
            do dir = 1, ndir  ! ndir = 1 for E or N-cell, ndir = 2 for C-cell
               do l = 1, l_size ! loop over number of fields
                  tMe = min(l, ntileMe)
                  ptr_field = f_addrs(l)
                  do tNbr = 1, ntileNbr
                     overPtr => bound%list(m)%send(tNbr,tMe,dir)
                     if( overPtr%overlap(1) ) then
                        is = overPtr%is(1); ie = overPtr%ie(1)
                        js = overPtr%js(1); je = overPtr%je(1)
                        do k = 1,ke  
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 buffer(pos) = field(i,j,k)
                              end do
                           end do
                        end do
                     end if
                  end do
               end do ! end do l = 1, l_size 
            end do ! end do dir = 1, ndir
            msgsize = pos - buffer_pos
            if( msgsize.GT.0 )then
               to_pe = domain%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, ' // &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
               end if
               call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
               buffer_pos = pos
            end if
         end do ! end do list = 0,nlist-1

         !--- recv the data 
         do list = 0,nlist-1
            m = mod( domain%pos+nlist-list, nlist )
            if( .NOT. bound%list(m)%overlap )cycle
            ntileNbr = size(domain%list(m)%x(:))            
            msgsize = 0
            do dir = 1, ndir  
               do tNbr = 1, ntileNbr
                  do l = 1, l_size  ! loop over number of fields
                     tMe = min(l, ntileMe)
                     overPtr => bound%list(m)%recv(tMe,tNbr,dir)
                     if( overPtr%overlap(1) ) then
                        is = overPtr%is(1); ie = overPtr%ie(1)
                        js = overPtr%js(1); je = overPtr%je(1)
                        msgsize = msgsize + (ie-is+1)*(je-js+1)
                     end if
                  end do
               end do
            end do
            msgsize = msgsize*ke

            if( msgsize.GT.0 )then
               from_pe = domain%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, '// &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
               end if
               call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
               buffer_pos = buffer_pos + msgsize
            end if
         end do

         !--- compare the data in reverse order
         do list = nlist-1,0,-1
            m = mod( domain%pos+nlist-list, nlist )
            if( .NOT.bound%list(m)%overlap)cycle
            pos = buffer_pos
            do dir = ndir, 1, -1
               do tNbr = ntileNbr, 1, -1
                  do l=l_size,1,-1  ! loop over number of fields
                     tMe = min(l, ntileMe)
                     ptr_field = f_addrs(l)
                     overPtr => bound%list(m)%recv(tMe,tNbr,dir)
                     if( overPtr%overlap(1) ) then
                        is = overPtr%is(1); ie = overPtr%ie(1)
                        js = overPtr%js(1); je = overPtr%je(1)
                        msgsize = (ie-is+1)*(je-js+1)*ke
                        pos = buffer_pos - msgsize
                        buffer_pos = pos
                        select case ( overPtr%rotation )
                        case ( ZERO )
                           do k = 1,ke
                              do j = js, je
                                 do i = is, ie
                                    pos = pos + 1
                                    if( field(i,j,k) .NE. buffer(pos) ) then
                                       print*,"Error from mpp_do_update_old.h on pe = ", mpp_pe(), ": field ", &
                                            trim(field_name), " at point (", i, ",", j, ",", k, ") = ", field(i,j,k), &
                                            " does not equal to the value = ", buffer(pos), " on pe ", domain%list(m)%pe
                                       call mpp_error(FATAL, "mpp_do_update_old.h: mismatch on the boundary for symmetry point")
                                    end if
                                 end do
                              end do
                           end do
                        case ( MINUS_NINETY )
                           do k = 1,ke
                              do i = ie, is, -1
                                 do j = js, je
                                    pos = pos + 1
                                    if( field(i,j,k) .NE. buffer(pos) ) then
                                       print*,"Error from mpp_do_update_old.h on pe = ", mpp_pe(), ": field ", &
                                            trim(field_name), " at point (", i, ",", j, ",", k, ") = ", field(i,j,k), &
                                            " does not equal to the value = ", buffer(pos), " on pe ", domain%list(m)%pe
                                       call mpp_error(FATAL, "mpp_do_update_old.h: mismatch on the boundary for symmetry point")
                                    end if
                                 end do
                              end do
                           end do
                        case ( NINETY )
                           do k = 1,ke
                              do i = is, ie
                                 do j = je, js, -1
                                    pos = pos + 1
                                    if( field(i,j,k) .NE. buffer(pos) ) then
                                       print*,"Error from mpp_do_update_old.h on pe = ", mpp_pe(), ": field ", &
                                            trim(field_name), " at point (", i, ",", j, ",", k, ") = ", field(i,j,k), &
                                            " does not equal to the value = ", buffer(pos), " on pe ", domain%list(m)%pe
                                       call mpp_error(FATAL, "mpp_do_update_old.h: mismatch on the boundary for symmetry point")
                                    end if
                                 end do
                              end do
                           end do
                        case ( ONE_HUNDRED_EIGHTY )
                           do k = 1,ke
                              do j = je, js, -1
                                 do i = ie, is, -1
                                    pos = pos + 1
                                    if( field(i,j,k) .NE. buffer(pos) ) then
                                       print*,"Error from mpp_do_update_old.h on pe = ", mpp_pe(), ": field ", &
                                            trim(field_name), " at point (", i, ",", j, ",", k, ") = ", field(i,j,k), &
                                            " does not equal to the value = ", buffer(pos), " on pe ", domain%list(m)%pe
                                       call mpp_error(FATAL, "mpp_do_update_old.h: mismatch on the boundary for symmetry point")
                                    end if
                                 end do
                              end do
                           end do
                        end select
                     end if ! end if ( overPtr%overlap(1) )
                  end do ! end do l=l_size,1,-1
               end do ! end do tNbr = ntileNbr, 1, -1
            end do ! end do dir = bound%ndir, 1, -1
         end do ! end do list = nlist-1,0,-1
      end if ! end if( check ... )

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)
      return
    end subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_
