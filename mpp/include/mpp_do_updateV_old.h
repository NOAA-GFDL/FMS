    subroutine MPP_DOMAINS_DO_UPDATE_3Dold_V_( fieldx, fieldy, domain, flags, gridtype)
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(inout) :: domain
      MPP_TYPE_, intent(inout), dimension(domain%x%data%begin:,domain%y%data%begin:,:) :: fieldx, fieldy
      integer, intent(in), optional :: flags, gridtype

      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif

      integer :: update_flags, grid_offset_type
      integer :: buffer_pos, nsend, nrecv
      integer :: i, j, k, m, n, l, is, ie, js, je, ke, nlist, is2, ie2
      integer :: to_pe, from_pe, list, pos, msgsize
      type(domain2d), pointer :: Dom_x => NULL()
      type(domain2d), pointer :: Dom_y => NULL()
      logical :: send(8), recv(8)
      character(len=8) :: text
      logical :: do_flip, do_flip2, flip_n, flip_w, flip_e

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
      ke = size(fieldx,3)
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking

#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
#endif
      
      !grid_offset_type used by update domains to determine shifts.
      grid_offset_type = AGRID
      if( PRESENT(gridtype) ) grid_offset_type = gridtype

      select case(grid_offset_type)
      case(AGRID)
         Dom_x => domain%T
         Dom_y => domain%T
      case(BGRID_NE, BGRID_SW)
         Dom_x => domain%C
         Dom_y => domain%C
      case(CGRID_NE, CGRID_SW)
         Dom_x => domain%E
         Dom_y => domain%N
      case default
         call mpp_error( FATAL, 'MPP_DOMAINS_DO_UPDATE_3Dold_V: gridtype must be one of '// &
                                'AGRID|BGRID_NE|BGRID_SW|CGRID_NE|CGRID_SW.' )
      end select

      ! send
      do list = 0,nlist-1
         m = mod( Dom_x%pos+list, nlist )
         if( .NOT. Dom_x%list(m)%overlap(1) .AND. .NOT. Dom_y%list(m)%overlap(1) )cycle
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos
         select case ( grid_offset_type )
         case(BGRID_NE, BGRID_SW, AGRID)
            do l = 1, 8  ! loop over 8 direction
               if(send(l)) then
                  if( Dom_x%list(m)%send(l)%overlap(1) ) then
!                     call mpp_clock_begin(pack_loop_clock) ! z1l not necessary for this clock, no
                     is = Dom_x%list(m)%send(l)%is; ie = Dom_x%list(m)%send(l)%ie
                     js = Dom_x%list(m)%send(l)%js; je = Dom_x%list(m)%send(l)%je
                     do k = 1,ke  
                        do j = js, je
                           do i = is, ie
                              pos = pos + 2
                              buffer(pos-1) = fieldx(i,j,k)
                              buffer(pos)   = fieldy(i,j,k)
                           end do
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  end if
                  if( Dom_x%list(m)%send(l)%overlap(2) ) then
!                     call mpp_clock_begin(pack_loop_clock)
                     nsend = Dom_x%list(m)%send(l)%n
                     do n = 1, nsend
                        i = Dom_x%list(m)%send(l)%i(n); j = Dom_x%list(m)%send(l)%j(n)
                        do k = 1,ke  
                           pos = pos + 2
                           buffer(pos-1) = fieldx(i,j,k)
                           buffer(pos)   = fieldy(i,j,k)
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  end if
               end if
            enddo
         case(CGRID_NE, CGRID_SW)
            do l = 1, 8  ! loop over 8 direction
               if(send(l)) then
                  if( Dom_x%list(m)%send(l)%overlap(1) ) then
!                     call mpp_clock_begin(pack_loop_clock)
                     is = Dom_x%list(m)%send(l)%is; ie = Dom_x%list(m)%send(l)%ie
                     js = Dom_x%list(m)%send(l)%js; je = Dom_x%list(m)%send(l)%je
                     do k = 1,ke  
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = fieldx(i,j,k)
                           end do
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  end if
                  if( Dom_x%list(m)%send(l)%overlap(2) ) then
!                     call mpp_clock_begin(pack_loop_clock)
                     nsend = Dom_x%list(m)%send(l)%n
                     do n = 1, nsend
                        i = Dom_x%list(m)%send(l)%i(n); j = Dom_x%list(m)%send(l)%j(n)
                        do k = 1,ke  
                           pos = pos + 1
                           buffer(pos) = fieldx(i,j,k)
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  end if
                  if(Dom_y%list(m)%send(l)%overlap(1)) then
!                     call mpp_clock_begin(pack_loop_clock)
                     is = Dom_y%list(m)%send(l)%is; ie = Dom_y%list(m)%send(l)%ie
                     js = Dom_y%list(m)%send(l)%js; je = Dom_y%list(m)%send(l)%je
                     do k = 1,ke  
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = fieldy(i,j,k)
                           end do
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  endif
                  if(Dom_y%list(m)%send(l)%overlap(2)) then
!                     call mpp_clock_begin(pack_loop_clock)
                     nsend = Dom_y%list(m)%send(l)%n
                     do n = 1, nsend
                        i = Dom_y%list(m)%send(l)%i(n); j = Dom_y%list(m)%send(l)%j(n)
                        do k = 1,ke  
                           pos = pos + 1
                           buffer(pos) = fieldy(i,j,k)
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  endif
               end if
            enddo
         end select
         call mpp_clock_end(pack_clock)

         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = Dom_x%list(m)%pe
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
         m = mod( Dom_x%pos+nlist-list, nlist )
         if( .NOT. Dom_x%list(m)%overlap(1) .AND. .NOT. Dom_y%list(m)%overlap(1) )cycle
         call mpp_clock_begin(recv_clock)
         msgsize = 0

         select case (grid_offset_type)
         case(BGRID_NE, BGRID_SW, AGRID)
            do l = 1, 8  ! loop over 8 direction
               if(recv(l)) then
                  if( Dom_x%list(m)%recv(l)%overlap(1) ) then
                     is = Dom_x%list(m)%recv(l)%is; ie = Dom_x%list(m)%recv(l)%ie
                     js = Dom_x%list(m)%recv(l)%js; je = Dom_x%list(m)%recv(l)%je
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
                  msgsize = msgsize + Dom_x%list(m)%recv(l)%n
               end if
            end do
            msgsize = msgsize*2
         case(CGRID_NE, CGRID_SW)
            do l = 1, 8  ! loop over 8 direction
               if(recv(l)) then
                  if( Dom_x%list(m)%recv(l)%overlap(1) ) then
                     is = Dom_x%list(m)%recv(l)%is; ie = Dom_x%list(m)%recv(l)%ie
                     js = Dom_x%list(m)%recv(l)%js; je = Dom_x%list(m)%recv(l)%je
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
                  if( Dom_y%list(m)%recv(l)%overlap(1) ) then
                     is = Dom_y%list(m)%recv(l)%is; ie = Dom_y%list(m)%recv(l)%ie
                     js = Dom_y%list(m)%recv(l)%js; je = Dom_y%list(m)%recv(l)%je
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
                  msgsize = msgsize + Dom_x%list(m)%recv(l)%n + Dom_y%list(m)%recv(l)%n
               end if
            enddo
         end select
         msgsize = msgsize*ke

         if( msgsize.GT.0 )then
            from_pe = Dom_x%list(m)%pe
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
         m = mod( Dom_x%pos+nlist-list, nlist )
         if( .NOT.Dom_x%list(m)%overlap(1) .AND. .NOT.Dom_y%list(m)%overlap(1) )cycle
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         select case( grid_offset_type )
         case(BGRID_NE, BGRID_SW, AGRID)
            do l = 8,1,-1
               if(recv(l)) then
                  if( Dom_x%list(m)%recv(l)%overlap(2) ) then
                     nrecv = Dom_x%list(m)%recv(l)%n
                     msgsize = nrecv*ke*2
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = Dom_x%list(m)%recv(l)%i(n); j = Dom_x%list(m)%recv(l)%j(n)
                        do k = 1,ke
                           pos = pos + 2
                           fieldx(i,j,k) = buffer(pos-1)
                           fieldy(i,j,k) = buffer(pos)
                        end do
                     end do
                  endif
                  if( Dom_x%list(m)%recv(l)%overlap(1) ) then
                     is = Dom_x%list(m)%recv(l)%is; ie = Dom_x%list(m)%recv(l)%ie
                     js = Dom_x%list(m)%recv(l)%js; je = Dom_x%list(m)%recv(l)%je
                     msgsize = (ie-is+1)*(je-js+1)*ke*2
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     select case ( Dom_x%list(m)%recv(l)%rotation )
                     case ( ZERO )
                        do k = 1,ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 2
                                 fieldx(i,j,k) = buffer(pos-1)
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     case ( ONE_HUNDRED_EIGHTY )
                        do k = 1,ke
                           do j = je, js, -1
                              do i = ie, is, -1
                                 pos = pos + 2
                                 fieldx(i,j,k) = buffer(pos-1)
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end select
                  endif
               end if
            enddo
         case(CGRID_NE, CGRID_SW)
            do l = 8,1,-1
               if(recv(l)) then
                  if( Dom_y%list(m)%recv(l)%overlap(2) ) then
                     nrecv = Dom_y%list(m)%recv(l)%n
                     msgsize = nrecv*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = Dom_y%list(m)%recv(l)%i(n); j = Dom_y%list(m)%recv(l)%j(n)
                        do k = 1,ke
                           pos = pos + 1
                           fieldy(i,j,k) = buffer(pos)
                        end do
                     end do
                  endif
                  if( Dom_y%list(m)%recv(l)%overlap(1) ) then
                     is = Dom_y%list(m)%recv(l)%is; ie = Dom_y%list(m)%recv(l)%ie
                     js = Dom_y%list(m)%recv(l)%js; je = Dom_y%list(m)%recv(l)%je
                     msgsize = (ie-is+1)*(je-js+1)*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     select case ( Dom_y%list(m)%recv(l)%rotation )
                     case ( ZERO )
                        do k = 1,ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     case ( ONE_HUNDRED_EIGHTY )
                        do k = 1,ke
                           do j = je, js, -1
                              do i = ie, is, -1
                                 pos = pos + 1
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end select
                  endif
                  if(Dom_x%list(m)%recv(l)%overlap(2)) then
                     nrecv = Dom_x%list(m)%recv(l)%n
                     msgsize = nrecv*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = Dom_x%list(m)%recv(l)%i(n); j = Dom_x%list(m)%recv(l)%j(n)
                        do k = 1,ke
                           pos = pos + 1
                           fieldx(i,j,k) = buffer(pos)
                        end do
                     end do
                  end if
                  if(Dom_x%list(m)%recv(l)%overlap(1)) then
                     is = Dom_x%list(m)%recv(l)%is; ie = Dom_x%list(m)%recv(l)%ie
                     js = Dom_x%list(m)%recv(l)%js; je = Dom_x%list(m)%recv(l)%je
                     msgsize = (ie-is+1)*(je-js+1)*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     select case ( Dom_x%list(m)%recv(l)%rotation )
                     case ( ZERO )                  
                        do k = 1,ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 fieldx(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     case ( ONE_HUNDRED_EIGHTY )
                        do k = 1,ke
                           do j = je, js, -1
                              do i = ie, is, -1
                                 pos = pos + 1
                                 fieldx(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end select
                  endif
               end if
            enddo
         end select
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
         !--- For cubic grid, velocity is not allowed on the CORNER.
         if( domain%topology_type == CUBIC_GRID .AND. ( .NOT. BTEST(update_flags,SCALAR_BIT)) .AND. &
              ( grid_offset_type == BGRID_NE .OR. grid_offset_type == BGRID_SW ) ) then
            call mpp_error(FATAL,"MPP_DOMAINS_DO_UPDATE_3Dold_V: velocity is not allowed on the C-cell for cubic grid")
         end if
         do list = 0,nlist-1
            m = mod( Dom_x%pos+list, nlist )
            if( .NOT. Dom_x%list(m)%overlap(3) .AND. .NOT. Dom_y%list(m)%overlap(3) )cycle
            call mpp_clock_begin(pack_clock)
            pos = buffer_pos

            do l = 1, 8  ! loop over 8 direction
               if(send(l)) then
                  if( Dom_x%list(m)%send(l)%overlap(3) ) then
                     !                     call mpp_clock_begin(pack_loop_clock)
                     is = Dom_x%list(m)%send(l)%is; ie = Dom_x%list(m)%send(l)%ie
                     js = Dom_x%list(m)%send(l)%js; je = Dom_x%list(m)%send(l)%je
                     if( Dom_x%list(m)%send(l)%rotation == ZERO .OR. grid_offset_type == BGRID_NE &
                                                                .OR. grid_offset_type == BGRID_SW ) then
                        do k = 1, ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 buffer(pos) = fieldx(i,j,k)
                              end do
                           end do
                        end do
                     else
                        do k = 1, ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 buffer(pos) = fieldy(i,j,k)
                              end do
                           end do
                        end do
                     end if
                  !                     call mpp_clock_end(pack_loop_clock)
                  endif
                  if( Dom_y%list(m)%send(l)%overlap(3) ) then
                     !                     call mpp_clock_begin(pack_loop_clock)
                     is = Dom_y%list(m)%send(l)%is; ie = Dom_y%list(m)%send(l)%ie
                     js = Dom_y%list(m)%send(l)%js; je = Dom_y%list(m)%send(l)%je
                     if( Dom_y%list(m)%send(l)%rotation == ZERO .OR. grid_offset_type == BGRID_NE &
                                                                .OR. grid_offset_type == BGRID_SW ) then
                        do k = 1, ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 buffer(pos) = fieldy(i,j,k)
                              end do
                           end do
                        end do
                     else
                        do k = 1, ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 buffer(pos) = fieldx(i,j,k)
                              end do
                           end do
                        end do
                     end if
                     !                     call mpp_clock_end(pack_loop_clock)
                  endif
               end if
            enddo
            call mpp_clock_end(pack_clock)

            call mpp_clock_begin(send_clock)
            msgsize = pos - buffer_pos
            if( msgsize.GT.0 )then
               to_pe = Dom_x%list(m)%pe
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
            m = mod( Dom_x%pos+nlist-list, nlist )
            if( .NOT. Dom_x%list(m)%overlap(3) .AND. .NOT. Dom_y%list(m)%overlap(3) )cycle
            call mpp_clock_begin(recv_clock)
            msgsize = 0

            do l = 1, 8  ! loop over 8 direction
               if(recv(l)) then
                  if(Dom_x%list(m)%recv(l)%overlap(3)) then
                     is = Dom_x%list(m)%recv(l)%is; ie = Dom_x%list(m)%recv(l)%ie
                     js = Dom_x%list(m)%recv(l)%js; je = Dom_x%list(m)%recv(l)%je
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  endif
                  if(Dom_y%list(m)%recv(l)%overlap(3)) then
                     is = Dom_y%list(m)%recv(l)%is; ie = Dom_y%list(m)%recv(l)%ie
                     js = Dom_y%list(m)%recv(l)%js; je = Dom_y%list(m)%recv(l)%je
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  endif
               end if
            enddo
            msgsize = msgsize*ke

            if( msgsize.GT.0 )then
               from_pe = Dom_x%list(m)%pe
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
            m = mod( Dom_x%pos+nlist-list, nlist )
            if( .NOT.Dom_x%list(m)%overlap(3) .AND. .NOT.Dom_y%list(m)%overlap(3) )cycle
            call mpp_clock_begin(unpk_clock)
            pos = buffer_pos
            do l = 8,1,-1
               if(recv(l)) then
                  if(Dom_y%list(m)%recv(l)%overlap(3)) then
                     is = Dom_y%list(m)%recv(l)%is; ie = Dom_y%list(m)%recv(l)%ie
                     js = Dom_y%list(m)%recv(l)%js; je = Dom_y%list(m)%recv(l)%je
                     msgsize = (ie-is+1)*(je-js+1)*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     select case ( Dom_y%list(m)%recv(l)%rotation )
                     case ( ZERO )
                        do k = 1,ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     case ( MINUS_NINETY )
                        do k = 1,ke
                           do i = ie, is, -1
                              do j = js, je
                                 pos = pos + 1
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     case ( NINETY )
                        if( BTEST(update_flags,SCALAR_BIT) ) then
                           do k = 1,ke
                              do i = is, ie
                                 do j = je, js, -1
                                    pos = pos + 1
                                    fieldy(i,j,k) = buffer(pos)
                                 end do
                              end do
                           end do
                        else 
                           do k = 1,ke
                              do i = is, ie
                                 do j = je, js, -1
                                    pos = pos + 1
                                    fieldy(i,j,k) = -buffer(pos)
                                 end do
                              end do
                           end do
                        end if
                     end select
                  endif
                  if(Dom_x%list(m)%recv(l)%overlap(3)) then
                     is = Dom_x%list(m)%recv(l)%is; ie = Dom_x%list(m)%recv(l)%ie
                     js = Dom_x%list(m)%recv(l)%js; je = Dom_x%list(m)%recv(l)%je
                     msgsize = (ie-is+1)*(je-js+1)*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     select case ( Dom_x%list(m)%recv(l)%rotation )
                     case ( ZERO )
                        do k = 1,ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 fieldx(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     case ( MINUS_NINETY )
                        if( BTEST(update_flags,SCALAR_BIT) ) then
                           do k = 1,ke
                              do i = ie, is, -1
                                 do j = js, je
                                    pos = pos + 1
                                    fieldx(i,j,k) = buffer(pos)
                                 end do
                              end do
                           end do
                        else
                           do k = 1,ke
                              do i = ie, is, -1
                                 do j = js, je
                                    pos = pos + 1
                                    fieldx(i,j,k) = -buffer(pos)
                                 end do
                              end do
                           end do
                        end if
                     case ( NINETY )
                        do k = 1,ke
                           do i = is, ie
                              do j = je, js, -1
                                 pos = pos + 1
                                 fieldx(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end select
                  endif
               end if
            enddo
            call mpp_clock_end(unpk_clock)
         end do

         !-------------------------------------------------------------------
         ! for the data on corner, some extra communication is needed to 
         ! ensure the corner value on different tile is the same.
         ! the following may not be needed.
         !-------------------------------------------------------------------
         if( grid_offset_type == BGRID_NE ) then
            do list = 0,nlist-1
               m = mod( Dom_x%pos+list, nlist )
               if( .NOT. Dom_x%list(m)%overlap(4) )cycle
               call mpp_clock_begin(pack_clock)
               pos = buffer_pos
               do l = 1, 8  ! loop over 8 direction
                  if(send(l) .AND. Dom_x%list(m)%send(l)%overlap(4)) then
!                     call mpp_clock_begin(pack_loop_clock)
                     i = Dom_x%list(m)%send(l)%i2; j = Dom_x%list(m)%send(l)%j2
                     do k = 1, ke
                        pos = pos + 2
                        buffer(pos-1) = fieldx(i,j,k)
                        buffer(pos)   = fieldy(i,j,k)
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  endif
               enddo
               call mpp_clock_end(pack_clock)

               call mpp_clock_begin(send_clock)
               msgsize = pos - buffer_pos
               if( msgsize.GT.0 )then
                  to_pe = Dom_x%list(m)%pe
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
               m = mod( Dom_x%pos+nlist-list, nlist )
               if( .NOT. Dom_x%list(m)%overlap(4) )cycle
               call mpp_clock_begin(recv_clock)
               msgsize = 0

               do l = 1, 8  ! loop over 8 direction
                  if(recv(l) .AND. Dom_x%list(m)%recv(l)%overlap(4)) then
                     msgsize = msgsize + ke*2
                  endif
               enddo

               if( msgsize.GT.0 )then
                  from_pe = Dom_x%list(m)%pe
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
               m = mod( Dom_x%pos+nlist-list, nlist )
               if( .NOT.Dom_x%list(m)%overlap(4) )cycle
               call mpp_clock_begin(unpk_clock)
               pos = buffer_pos
               do l = 8,1,-1
                  if(recv(l) .AND. Dom_x%list(m)%recv(l)%overlap(4)) then
                     msgsize = 2*ke   
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     i = Dom_x%list(m)%recv(l)%i2; j = Dom_x%list(m)%recv(l)%j2
                     do k = 1,ke
                        pos = pos + 2
                        fieldx(i,j,k) = buffer(pos-1)
                        fieldy(i,j,k) = buffer(pos)
                     end do
                  endif
               enddo
               call mpp_clock_end(unpk_clock)
            end do
         end if
      end if  ! if(ntiles > 1)

      !--- for vector field flip the sign if necessary.            
      if( BTEST(domain%fold,NORTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then
         if( BTEST(update_flags,NORTH) ) then
            js = domain%y%global%end + 1
            je = domain%y%data%end
            if( je.GE.js )then
               is = domain%x%data%begin
               ie = domain%x%data%end
               !flip the sign if this is a vector
               select case(grid_offset_type)
               case (AGRID) 
                  do k = 1,ke
                     do j = js,je
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               case (BGRID_NE) 
                  if(domain%symmetry) then
                     do k = 1,ke
                        do j = js+1,je+1
                           do i = is,ie+1
                              fieldx(i,j,k) = -fieldx(i,j,k)
                              fieldy(i,j,k) = -fieldy(i,j,k)
                           end do
                        end do
                     end do
                  else
                     do k = 1,ke
                        do j = js,je
                           do i = is,ie
                              fieldx(i,j,k) = -fieldx(i,j,k)
                              fieldy(i,j,k) = -fieldy(i,j,k)
                           end do
                        end do
                     end do
                  end if
               case (CGRID_NE)
                  if(domain%symmetry) then
                     do k = 1,ke
                        do j = js,je
                           do i = is,ie+1
                              fieldx(i,j,k) = -fieldx(i,j,k)
                           end do
                        end do
                        do j = js+1,je+1
                           do i = is,ie
                              fieldy(i,j,k) = -fieldy(i,j,k)
                           end do
                        end do
                     end do
                  else
                     do k = 1,ke
                        do j = js,je
                           do i = is,ie
                              fieldx(i,j,k) = -fieldx(i,j,k)
                              fieldy(i,j,k) = -fieldy(i,j,k)
                           end do
                        end do
                     end do
                  endif
               end select
            end if
         endif

         !--- flip the sign at the folded north edge.
         j = domain%y%global%end
         if( domain%y%data%begin.LE.j .AND. j.LE.domain%y%data%end )then !fold is within domain
            if( domain%x%pos.GE.size(domain%x%list(:))/2 )then
               do_flip = .false.
               flip_n = BTEST(update_flags,NORTH)
               flip_e = BTEST(update_flags,EAST)
               flip_w = BTEST(update_flags,WEST)
               do_flip = flip_n .or. flip_w .or. flip_e
               do_flip2 = .FALSE.
               if( flip_n .AND. flip_e .AND. flip_w ) then
                  is = max(domain%x%data%begin,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie = domain%x%data%end
               else if (  flip_n .AND. flip_e ) then
                  is = max(domain%x%compute%begin,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie = domain%x%data%end
               else if (  flip_n .AND. flip_w ) then
                  is = max(domain%x%data%begin,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie = domain%x%compute%end
               else if ( flip_n ) then
                  is = max(domain%x%compute%begin,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie = domain%x%compute%end                 
               else if ( flip_w .AND. flip_e ) then  ! a little complicate
                  is = max(domain%x%compute%end+1,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie = domain%x%data%end
                  do_flip2 = .TRUE.
                  is2 = max(domain%x%data%begin,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie2 = domain%x%compute%begin-1
               else if ( flip_e ) then
                  is = max(domain%x%compute%end+1,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie = domain%x%data%end
               else if ( flip_w ) then
                  is = max(domain%x%data%begin,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie = domain%x%compute%begin-1
               end if
               if( do_flip ) then
                  select case(grid_offset_type)
                  case(BGRID_NE)
                     if(domain%symmetry) then
                        j = j + 1
                        ie = ie + 1
                     end if
                     do k = 1,ke
                        do i = is, ie
                           fieldx(i,j,k) = -fieldx(i,j,k) 
                           fieldy(i,j,k) = -fieldy(i,j,k) 
                        end do
                     end do
                  case(CGRID_NE)
                     if(domain%symmetry) j = j+1
                     do k = 1,ke
                        do i = is, ie
                           fieldy(i,j,k) = -fieldy(i,j,k) 
                        end do
                     end do
                  end select
                  if( do_flip2 ) then
                     select case(grid_offset_type)
                     case(BGRID_NE)
                        if(domain%symmetry) then
                           ie2 = ie2 + 1
                        end if
                        do k = 1,ke
                           do i = is2, ie2
                              fieldx(i,j,k) = -fieldx(i,j,k) 
                              fieldy(i,j,k) = -fieldy(i,j,k) 
                           end do
                        end do
                     case(CGRID_NE)
                        do k = 1,ke
                           do i = is2, ie2
                              fieldy(i,j,k) = -fieldy(i,j,k) 
                           end do
                        end do
                     end select
                  end if
               end if
            end if
            !poles set to 0: BGRID only
            if( grid_offset_type.EQ.BGRID_NE )then
               if(domain%symmetry) then
                  j = domain%y%global%end + 1
                  do i = domain%x%global%begin,domain%x%global%end+1,(domain%x%global%begin+domain%x%global%end)/2
                     if( domain%x%data%begin.LE.i .AND. i.LE.domain%x%data%end )then
                        do k = 1,ke
                           fieldx(i,j,k) = 0.
                           fieldy(i,j,k) = 0.
                        end do
                     end if
                  end do
               else
                  j = domain%y%global%end 
                  do i = domain%x%global%begin-1,domain%x%global%end,(domain%x%global%begin+domain%x%global%end)/2
                     if( domain%x%data%begin.LE.i .AND. i.LE.domain%x%data%end )then
                        do k = 1,ke
                           fieldx(i,j,k) = 0.
                           fieldy(i,j,k) = 0.
                        end do
                     end if
                  end do
               endif
            end if

            ! the last code block correct an error where the data in your halo coming from 
            ! other half may have the wrong sign
            !right of midpoint, when update north and east direction
            if ( BTEST(update_flags,NORTH) .OR. BTEST(update_flags,EAST) ) then
               j = domain%y%global%end 
               is = (domain%x%global%begin+domain%x%global%end)/2

               if( domain%x%compute%begin.LE.is .AND. is.LT.domain%x%data%end .AND. domain%x%pos < size(domain%x%list(:))/2 )then
                  ie = domain%x%data%end
                  if(domain%symmetry) j = j+1
                  select case(grid_offset_type)
                  case(BGRID_NE)
                     if(domain%symmetry) then
                        is = is+2
                        ie = ie+1
                     else
                        is = is + 1
                     endif
                     do k = 1,ke
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  case(CGRID_NE)
                     is = is+1
                     do k = 1,ke
                        do i = is,ie
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end select
               end if
            end if
         end if
      end if

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)
 
      return

    end subroutine MPP_DOMAINS_DO_UPDATE_3Dold_V_
