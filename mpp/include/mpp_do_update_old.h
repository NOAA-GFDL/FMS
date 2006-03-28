    subroutine MPP_DOMAINS_DO_UPDATE_3Dold_( field, domain, flags, position)
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(in)    :: domain
      MPP_TYPE_, intent(inout)      :: field(domain%x%data%begin:,domain%y%data%begin:,:)
      integer, intent(in), optional :: flags
      integer, intent(in), optional :: position

      integer   :: update_flags, update_position
      type(domain2d), pointer :: Dom => NULL()

!equate to mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif
      integer :: buffer_pos
      integer :: i, j, k, m, n, l, nlist, nsend, nrecv
      integer :: is, ie, js, je, ke
!receive domains saved here for unpacking
!for non-blocking version, could be recomputed
      integer :: to_pe, from_pe, list, pos, msgsize
      logical :: send(8), recv(8)
      character(len=8) :: text

      update_position = CENTER
      if(PRESENT(position)) update_position = position

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

      nlist = size(domain%list(:))
      ke = size(field,3)
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking

#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
#endif

      ! select the domain      
      select case(update_position)
      case (CENTER)
         Dom => domain%T
      case (EAST)
         Dom => domain%E
      case (NORTH)
         Dom => domain%N
      case (CORNER)
         Dom => domain%C
      end select

      ! send
      do list = 0,nlist-1
         m = mod( Dom%pos+list, nlist )
         if( .NOT. Dom%list(m)%overlap(1) .AND. .NOT. Dom%list(m)%overlap(2) )cycle
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos

         do l = 1, 8  ! loop over 8 direction
            if( send(l) ) then
               if( Dom%list(m)%send(l)%overlap(1) ) then
   
!                  call mpp_clock_begin(pack_loop_clock)
                  is = Dom%list(m)%send(l)%is; ie = Dom%list(m)%send(l)%ie
                  js = Dom%list(m)%send(l)%js; je = Dom%list(m)%send(l)%je
                  do k = 1,ke  
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           buffer(pos) = field(i,j,k)
                        end do
                     end do
                  end do
!                  call mpp_clock_end(pack_loop_clock)
               endif
               if( Dom%list(m)%send(l)%overlap(2) ) then
!                  call mpp_clock_begin(pack_loop_clock)
                  nsend = Dom%list(m)%send(l)%n
                  do n = 1, nsend
                     i = Dom%list(m)%send(l)%i(n); j = Dom%list(m)%send(l)%j(n)
                     do k = 1,ke  
                        pos = pos + 1
                        buffer(pos) = field(i,j,k)
                     end do
                  end do
!                  call mpp_clock_end(pack_loop_clock)
               endif
            end if
         enddo
         call mpp_clock_end(pack_clock)

         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = Dom%list(m)%pe
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
      end do

      !recv
      do list = 0,nlist-1
         m = mod( Dom%pos+nlist-list, nlist )

         if( .NOT. Dom%list(m)%overlap(1) .AND. .NOT. Dom%list(m)%overlap(2) )cycle
         call mpp_clock_begin(recv_clock)
         msgsize = 0

         do l = 1, 8  ! loop over 8 direction
            if(recv(l)) then
              if( Dom%list(m)%recv(l)%overlap(1) ) then
                  is = Dom%list(m)%recv(l)%is; ie = Dom%list(m)%recv(l)%ie
                  js = Dom%list(m)%recv(l)%js; je = Dom%list(m)%recv(l)%je
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
              end if
              msgsize = msgsize + Dom%list(m)%recv(l)%n
            end if
         end do
         msgsize = msgsize*ke

         if( msgsize.GT.0 )then
            from_pe = Dom%list(m)%pe
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
      end do

      !unpack recv
      !unpack halos in reverse order
      do list = nlist-1,0,-1
         m = mod( Dom%pos+nlist-list, nlist )

         if( .NOT.Dom%list(m)%overlap(1) .AND. .NOT.Dom%list(m)%overlap(2))cycle
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         do l = 8,1,-1
            if( recv(l) ) then
               if( Dom%list(m)%recv(l)%overlap(2) ) then
                  nrecv = Dom%list(m)%recv(l)%n
                  msgsize = nrecv*ke
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  do n = 1, nrecv
                     i = Dom%list(m)%recv(l)%i(n); j = Dom%list(m)%recv(l)%j(n)
                     do k = 1,ke
                        pos = pos + 1
                        field(i,j,k) = buffer(pos)
                     end do
                  end do
               endif
               if( Dom%list(m)%recv(l)%overlap(1) ) then
                  is = Dom%list(m)%recv(l)%is; ie = Dom%list(m)%recv(l)%ie
                  js = Dom%list(m)%recv(l)%js; je = Dom%list(m)%recv(l)%je
                  msgsize = (ie-is+1)*(je-js+1)*ke
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  select case ( Dom%list(m)%recv(l)%rotation )
                  case ( ZERO )
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
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
               endif
            end if
         enddo
         call mpp_clock_end(unpk_clock)
      end do

      !------------------------------------------------------------------
      !  For the multiple-tile masaic, need to communication between tiles.
      !  Those overlapping regions are specified by rectangle. 
      !  The reason this step need to be seperated from the previous step,
      !  is because for cubic-grid, some points ( like at i=ieg+1 ), the overlapping
      !  is both within tile and between tiles. 
      !------------------------------------------------------------------
      !send
      if ( domain%ncontacts > 0 ) then
         !--- for cubic grid, can not use scalar version update_domain for data on E or N-cell
         if( domain%topology_type == CUBIC_GRID .AND. ( update_position == EAST .OR. update_position == NORTH ) ) then
            call mpp_error(FATAL, 'MPP_UPDATE_DOMAINS: for cubic grid, ' // &
                                  'can not use scalar version update_domain for data on E or N-cell' )
         end if

         do list = 0,nlist-1
            m = mod( Dom%pos+list, nlist )
            if( .NOT. Dom%list(m)%overlap(3) )cycle
            call mpp_clock_begin(pack_clock)
            pos = buffer_pos

            do l = 1, 8  ! loop over 8 direction
               if(send(l) .AND. Dom%list(m)%send(l)%overlap(3)) then
!                  call mpp_clock_begin(pack_loop_clock)
                  is = Dom%list(m)%send(l)%is; ie = Dom%list(m)%send(l)%ie
                  js = Dom%list(m)%send(l)%js; je = Dom%list(m)%send(l)%je
                  do k = 1, ke
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           buffer(pos) = field(i,j,k)
                        end do
                     end do
                  end do
!                  call mpp_clock_end(pack_loop_clock)
               endif
            enddo
            call mpp_clock_end(pack_clock)

            call mpp_clock_begin(send_clock)
            msgsize = pos - buffer_pos
            if( msgsize.GT.0 )then
               to_pe = Dom%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, ' // &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
               end if
               call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
               buffer_pos = pos
            end if
            call mpp_clock_end(send_clock)
         end do

         !recv
         do list = 0,nlist-1
            m = mod( Dom%pos+nlist-list, nlist )
            if( .NOT. Dom%list(m)%overlap(3) )cycle
            call mpp_clock_begin(recv_clock)
            msgsize = 0

            do l = 1, 8  ! loop over 8 direction
               if(recv(l) .AND. Dom%list(m)%recv(l)%overlap(3)) then
                  is = Dom%list(m)%recv(l)%is; ie = Dom%list(m)%recv(l)%ie
                  js = Dom%list(m)%recv(l)%js; je = Dom%list(m)%recv(l)%je
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               endif
            enddo
            msgsize = msgsize*ke

            if( msgsize.GT.0 )then
               from_pe = Dom%list(m)%pe
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
         end do

         !unpack recv
         !unpack halos in reverse order
         do list = nlist-1,0,-1
            m = mod( Dom%pos+nlist-list, nlist )

            if( .NOT.Dom%list(m)%overlap(3) )cycle
            call mpp_clock_begin(unpk_clock)
            pos = buffer_pos
            do l = 8,1,-1
               if(recv(l) .AND. Dom%list(m)%recv(l)%overlap(3)) then
                  is = Dom%list(m)%recv(l)%is; ie = Dom%list(m)%recv(l)%ie
                  js = Dom%list(m)%recv(l)%js; je = Dom%list(m)%recv(l)%je
                  msgsize = (ie-is+1)*(je-js+1)*ke
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  select case ( Dom%list(m)%recv(l)%rotation )
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
                  end select
               endif
            enddo
            call mpp_clock_end(unpk_clock)
         end do

         !-------------------------------------------------------------------
         ! for the data on corner, some extra communication is needed to 
         ! ensure the corner value on different tile is the same.
         !-------------------------------------------------------------------
         if( update_position == CORNER ) then
            do list = 0,nlist-1
               m = mod( Dom%pos+list, nlist )
               if( .NOT. Dom%list(m)%overlap(4) )cycle
               call mpp_clock_begin(pack_clock)
               pos = buffer_pos
               do l = 1, 8  ! loop over 8 direction
                  if(send(l) .AND. Dom%list(m)%send(l)%overlap(4)) then
!                     call mpp_clock_begin(pack_loop_clock)
                     i = Dom%list(m)%send(l)%i2; j = Dom%list(m)%send(l)%j2
                     do k = 1, ke
                        pos = pos + 1
                        buffer(pos) = field(i,j,k)
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  endif
               enddo
               call mpp_clock_end(pack_clock)

               call mpp_clock_begin(send_clock)
               msgsize = pos - buffer_pos
               if( msgsize.GT.0 )then
                  to_pe = Dom%list(m)%pe
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
            end do

            !recv
            do list = 0,nlist-1
               m = mod( Dom%pos+nlist-list, nlist )
               if( .NOT. Dom%list(m)%overlap(4) )cycle
               call mpp_clock_begin(recv_clock)
               msgsize = 0

               do l = 1, 8  ! loop over 8 direction
                  if(recv(l) .AND. Dom%list(m)%recv(l)%overlap(4)) then
                     msgsize = msgsize + ke
                  endif
               enddo

               if( msgsize.GT.0 )then
                  from_pe = Dom%list(m)%pe
                  mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, buffer_pos+msgsize )
                  if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                     write( text,'(i8)' )mpp_domains_stack_hwm
                     call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, '// &
                          'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
                  end if
                  call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
                  buffer_pos = buffer_pos + msgsize
               end if
               call mpp_clock_end(recv_clock)
            end do

            !unpack recv
            !unpack halos in reverse order
            do list = nlist-1,0,-1
               m = mod( Dom%pos+nlist-list, nlist )
               if( .NOT.Dom%list(m)%overlap(4) )cycle
               call mpp_clock_begin(unpk_clock)
               pos = buffer_pos
               do l = 8,1,-1
                  if(recv(l) .AND. Dom%list(m)%recv(l)%overlap(4)) then
                     msgsize = ke   
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     i = Dom%list(m)%recv(l)%i2; j = Dom%list(m)%recv(l)%j2
                     do k = 1,ke
                        pos = pos + 1
                        field(i,j,k) = buffer(pos)
                     end do
                  endif
               enddo
               call mpp_clock_end(unpk_clock)
            end do
         end if
      end if  ! if(ncontacts > 0)

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)
      return
    end subroutine MPP_DOMAINS_DO_UPDATE_3Dold_
