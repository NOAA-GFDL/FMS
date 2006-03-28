    subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_V_(f_addrsx,f_addrsy, d_comm,isize,jsize,ke, d_type,flags,gridtype)
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND), intent(in)         :: f_addrsx(:), f_addrsy(:)
      type(DomainCommunicator2D),    pointer :: d_comm
      integer,            intent(in)         :: isize(:), jsize(:), ke
      MPP_TYPE_, intent(in)                  :: d_type  ! creates unique interface
      integer, intent(in), optional :: flags, gridtype

      MPP_TYPE_ :: fieldx(d_comm%domain%x%data%begin:d_comm%domain%x%data%begin+isize(1)-1, &
                          d_comm%domain%y%data%begin:d_comm%domain%y%data%begin+jsize(1)-1,ke)
      pointer(ptr_fieldx, fieldx)
      MPP_TYPE_ :: fieldy(d_comm%domain%x%data%begin:d_comm%domain%x%data%begin+isize(2)-1, &
                          d_comm%domain%y%data%begin:d_comm%domain%y%data%begin+jsize(2)-1,ke)
      pointer(ptr_fieldy, fieldy)
      integer :: update_flags, grid_offset_type
      integer :: l_size, l, i, j, k, is, ie, js, je, m, n, is2, ie2
      integer :: pos, nlist, list, msgsize, nsend, nrecv
      integer :: to_pe, from_pe
      type(domain2d), pointer :: domain => NULL()
      type(DomainCommunicator2D), pointer :: d_commy => NULL()
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer(ptr,buffer )
      integer :: buffer_pos
      logical :: do_flip, do_flip2, flip_n, flip_w, flip_e
 
!--- the following will be moved into mpp_update_domain vector version
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

      grid_offset_type = AGRID
      if( PRESENT(gridtype) ) grid_offset_type = gridtype

      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
      l_size = size(f_addrsx(:))
      nlist = d_comm%Rlist_size
      ptr = LOC(mpp_domains_stack)

      if( d_comm%staggered ) then
         d_commy => d_comm%y_comm
      else
         d_commy => d_comm
      end if

      !--- send
      do list = 0,nlist-1
         if( .NOT.d_comm%S_do_buf(list) .AND. .NOT.d_commy%S_do_buf(list) )cycle
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos
         select case( grid_offset_type )
         case(BGRID_NE, BGRID_SW, AGRID)
            do l=1,l_size  ! loop over number of fields
               ptr_fieldx = f_addrsx(l)
               ptr_fieldy = f_addrsy(l)
               do m=1,8
                  if( d_comm%do_thisS(m,list) )then
!                     call mpp_clock_begin(pack_loop_clock) ! z1l: not necessary, not unpack_loop_clock
                     is = d_comm%send(m,list)%is; ie = d_comm%send(m,list)%ie
                     js = d_comm%send(m,list)%js; je = d_comm%send(m,list)%je
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

                  if( d_comm%do_thisS2(m,list) )then
!                     call mpp_clock_begin(pack_loop_clock)
                     nsend = d_comm%send(m,list)%n
                     do n = 1, nsend
                        i = d_comm%send(m,list)%i(n)
                        j = d_comm%send(m,list)%j(n)
                        do k = 1,ke
                           pos = pos + 2
                           buffer(pos-1) = fieldx(i,j,k)
                           buffer(pos) = fieldy(i,j,k)
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  end if

               end do  ! do m=1,8
            end do  ! do l=1,l_size
         case(CGRID_NE, CGRID_SW)
            do l=1,l_size  ! loop over number of fields
               ptr_fieldx = f_addrsx(l)
               ptr_fieldy = f_addrsy(l)
               do m=1,8
                  if( d_comm%do_thisS(m,list) )then
!                     call mpp_clock_begin(pack_loop_clock)
                     is = d_comm%send(m,list)%is; ie = d_comm%send(m,list)%ie
                     js = d_comm%send(m,list)%js; je = d_comm%send(m,list)%je
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
                  if( d_comm%do_thisS2(m,list) )then
!                     call mpp_clock_begin(pack_loop_clock)
                     nsend = d_comm%send(m,list)%n
                     do n = 1, nsend
                        i = d_comm%send(m,list)%i(n)
                        j = d_comm%send(m,list)%j(n)
                        do k = 1,ke
                           pos = pos + 1
                           buffer(pos) = fieldx(i,j,k)
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  end if
                  if( d_commy%do_thisS(m,list) )then
!                     call mpp_clock_begin(pack_loop_clock)
                     is = d_commy%send(m,list)%is; ie = d_commy%send(m,list)%ie
                     js = d_commy%send(m,list)%js; je = d_commy%send(m,list)%je
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = fieldy(i,j,k)
                           end do
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  end if
                  if( d_commy%do_thisS2(m,list) )then
!                     call mpp_clock_begin(pack_loop_clock)
                     nsend = d_commy%send(m,list)%n
                     do n = 1, nsend
                        i = d_commy%send(m,list)%i(n)
                        j = d_commy%send(m,list)%j(n)
                        do k = 1,ke
                           pos = pos + 1
                           buffer(pos) = fieldy(i,j,k)
                        end do
                     end do
!                     call mpp_clock_end(pack_loop_clock)
                  end if
               end do  ! do m=1,8
            end do  ! do l=1,l_size
         end select
         call mpp_clock_end(pack_clock)
         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = d_comm%cto_pe(list)
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
            buffer_pos = pos
         end if
         call mpp_clock_end(send_clock)
         !     write(stdout(),*) 'Update send checksum=',mpp_chksum(rbuffer(buffer_pos))
      end do

!recv
      nlist = d_comm%Rlist_size
      do list = 0,nlist-1
         if( .NOT. d_comm%R_do_buf(list) .AND. .NOT. d_commy%R_do_buf(list) )cycle
         call mpp_clock_begin(recv_clock)
         from_pe = d_comm%cfrom_pe(list)
         msgsize = 0
         select case(grid_offset_type)
         case(BGRID_NE, BGRID_SW, AGRID)
            do m=1,8
               msgsize = msgsize + d_comm%R_msize(m,list)
            end do
            msgsize = msgsize*2
         case(CGRID_NE, CGRID_SW)
            do m=1,8
               msgsize = msgsize + d_comm%R_msize(m,list)
               msgsize = msgsize + d_commy%R_msize(m,list)
            end do
         end select
         msgsize = msgsize*l_size
   
         if( msgsize.GT.0 )then
             call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
             buffer_pos = buffer_pos + msgsize
         end if
         call mpp_clock_end(recv_clock)
      end do

!unpack recv
!unpack halos in reverse order

      do list = nlist-1,0,-1
         if( .NOT.d_comm%R_do_buf(list) .AND. .NOT.d_commy%R_do_buf(list) )cycle
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         select case ( grid_offset_type )
         case(BGRID_NE, BGRID_SW, AGRID)
            do l=l_size,1,-1  ! loop over number of fields
               ptr_fieldx = f_addrsx(l)
               ptr_fieldy = f_addrsy(l)

               do m=8,1,-1
                  if( d_comm%do_thisR2(m,list) )then
                     nrecv = d_comm%recv(m,list)%n
                     msgsize = nrecv*ke*2
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = d_comm%recv(m,list)%i(n)
                        j = d_comm%recv(m,list)%j(n)
                        do k = 1,ke
                           pos = pos + 2
                           fieldx(i,j,k) = buffer(pos-1)
                           fieldy(i,j,k) = buffer(pos)
                        end do
                     end do
                  end if

                  if( d_comm%do_thisR(m,list) )then  ! direction
                     is = d_comm%recv(m,list)%is; ie = d_comm%recv(m,list)%ie
                     js = d_comm%recv(m,list)%js; je = d_comm%recv(m,list)%je 
                     msgsize = (ie-is+1)*(je-js+1)*ke*2
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     select case ( d_comm%recv(m,list)%rotation )
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
                           do j = je,js,-1
                              do i = ie,is,-1
                                 pos = pos + 2
                                 fieldx(i,j,k) = buffer(pos-1)
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end select
                  end if
               end do  ! do m=8,1,-1
            end do  ! do l=l_size,1,-1
         case(CGRID_NE, CGRID_SW)
            do l=l_size,1,-1  ! loop over number of fields
               ptr_fieldx = f_addrsx(l)
               ptr_fieldy = f_addrsy(l)

               do m=8,1,-1
                  if( d_commy%do_thisR2(m,list) )then  ! direction
                     nrecv = d_commy%recv(m,list)%n
                     msgsize = nrecv*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos

                     do n = 1, nrecv
                        i = d_commy%recv(m,list)%i(n)
                        j = d_commy%recv(m,list)%j(n)
                        do k = 1,ke
                           pos = pos + 1
                           fieldy(i,j,k) = buffer(pos)
                        end do
                     end do
                  end if

                  if( d_commy%do_thisR(m,list) )then  ! direction
                     is = d_commy%recv(m,list)%is; ie = d_commy%recv(m,list)%ie
                     js = d_commy%recv(m,list)%js; je = d_commy%recv(m,list)%je 
                     msgsize = (ie-is+1)*(je-js+1)*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     select case ( d_commy%recv(m,list)%rotation )
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
                           do j = je,js,-1
                              do i = ie,is,-1
                                 pos = pos + 1
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end select
                  end if

                  if( d_comm%do_thisR2(m,list) )then
                     nrecv = d_comm%recv(m,list)%n
                     msgsize = nrecv*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = d_comm%recv(m,list)%i(n)
                        j = d_comm%recv(m,list)%j(n)
                        do k = 1,ke
                           pos = pos + 1
                           fieldx(i,j,k) = buffer(pos)
                        end do
                     end do
                  end if

                  if( d_comm%do_thisR(m,list) )then  ! direction
                     is = d_comm%recv(m,list)%is; ie = d_comm%recv(m,list)%ie
                     js = d_comm%recv(m,list)%js; je = d_comm%recv(m,list)%je 
                     msgsize = (ie-is+1)*(je-js+1)*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     select case ( d_comm%recv(m,list)%rotation )
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
                           do j = je,js,-1
                              do i = ie,is,-1
                                 pos = pos + 1
                                 fieldx(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end select
                  end if
               end do  ! do m=8,1,-1
            end do  ! do l=l_size,1,-1
         end select
         call mpp_clock_end(unpk_clock)
!        write(stdout(),*) 'Update field checksum=',mpp_chksum(rfield)
      end do

      !------------------------------------------------------------------
      !  For the multiple-tile masaic, need to communication between tiles.
      !  Those overlapping regions are specified by rectangle. 
      !  The reason this step need to be seperated from the previous step,
      !  is because for cubic-grid, some points ( like at i=ieg+1 ), the overlapping
      !  is both within tile and between tiles. 
      !------------------------------------------------------------------

      domain => d_comm%domain
      !--- send
      if ( domain%ncontacts > 0 ) then
         do list = 0,nlist-1
            if( .NOT.d_comm%S_do_buf3(list) .AND. .NOT.d_commy%S_do_buf3(list) )cycle
            call mpp_clock_begin(pack_clock)
            to_pe = d_comm%cto_pe(list)
            pos = buffer_pos
            do l=1,l_size  ! loop over number of fields
               ptr_fieldx = f_addrsx(l)
               ptr_fieldy = f_addrsy(l)
               do m=1,8
                  if( d_comm%do_thisS3(m,list) )then
                     !                     call mpp_clock_begin(pack_loop_clock)
                     is = d_comm%send(m,list)%is; ie = d_comm%send(m,list)%ie
                     js = d_comm%send(m,list)%js; je = d_comm%send(m,list)%je
                     if( d_commy%send(m,list)%rotation == ZERO .OR. grid_offset_type == BGRID_NE &
                                                               .OR. grid_offset_type == BGRID_SW ) then
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
                                 buffer(pos) = fieldy(i,j,k)
                              end do
                           end do
                        end do
                     end if
!                     call mpp_clock_end(pack_loop_clock)
                  end if
                  if( d_commy%do_thisS3(m,list) )then
!                     call mpp_clock_begin(pack_loop_clock)
                     is = d_commy%send(m,list)%is; ie = d_commy%send(m,list)%ie
                     js = d_commy%send(m,list)%js; je = d_commy%send(m,list)%je
                     if( d_commy%send(m,list)%rotation == ZERO .OR. grid_offset_type == BGRID_NE &
                                                               .OR. grid_offset_type == BGRID_SW ) then
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
                                 buffer(pos) = fieldx(i,j,k)
                              end do
                           end do
                        end do
                     end if
!                     call mpp_clock_end(pack_loop_clock)
                  end if
               end do  ! do m=1,8
            end do  ! do l=1,l_size
            call mpp_clock_end(pack_clock)
            call mpp_clock_begin(send_clock)
            msgsize = pos - buffer_pos
            if( msgsize.GT.0 )then
               call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
               buffer_pos = pos
            end if
            call mpp_clock_end(send_clock)
            !     write(stdout(),*) 'Update send checksum=',mpp_chksum(rbuffer(buffer_pos))
         end do

         !recv
         nlist = d_comm%Rlist_size
         do list = 0,nlist-1
            if( .NOT. d_comm%R_do_buf3(list) .AND. .NOT. d_comm%R_do_buf3(list) )cycle
            call mpp_clock_begin(recv_clock)
            from_pe = d_comm%cfrom_pe(list)
            msgsize = 0
            do m=1,8
               if( d_comm%do_thisR3(m,list) ) msgsize = msgsize + d_comm%R_msize2(m,list)
               if( d_commy%do_thisR3(m,list) ) msgsize = msgsize + d_commy%R_msize2(m,list)
            end do
            msgsize = msgsize*l_size

            if( msgsize.GT.0 )then
               call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
               buffer_pos = buffer_pos + msgsize
            end if
            call mpp_clock_end(recv_clock)
         end do

         !unpack recv
         !unpack halos in reverse order

         do list = nlist-1,0,-1
            if( .NOT.d_comm%R_do_buf3(list) .AND. .NOT.d_commy%R_do_buf3(list) )cycle
            call mpp_clock_begin(unpk_clock)
            pos = buffer_pos

            do l=l_size,1,-1  ! loop over number of fields
               ptr_fieldx = f_addrsx(l)
               ptr_fieldy = f_addrsy(l)

               do m=8,1,-1
                  if( d_commy%do_thisR3(m,list) )then  ! direction
                     msgsize = d_commy%R_msize2(m,list)
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     is = d_commy%recv(m,list)%is; ie = d_commy%recv(m,list)%ie
                     js = d_commy%recv(m,list)%js; je = d_commy%recv(m,list)%je
                     select case ( d_commy%recv(m,list)%rotation )
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
                  end if

                  if( d_comm%do_thisR3(m,list) )then
                     msgsize = d_comm%R_msize2(m,list)
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     is = d_comm%recv(m,list)%is; ie = d_comm%recv(m,list)%ie
                     js = d_comm%recv(m,list)%js; je = d_comm%recv(m,list)%je
                     select case ( d_comm%recv(m,list)%rotation )
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
                  end if
               end do  ! do m=8,1,-1
            end do  ! do l=l_size,1,-1

            call mpp_clock_end(unpk_clock)
            !        write(stdout(),*) 'Update field checksum=',mpp_chksum(rfield)
         end do

         !--- to make sure the value on the corner is consistent between tiles.
         if( grid_offset_type == BGRID_NE ) then
            do list = 0,nlist-1
               if( .NOT.d_comm%S_do_buf4(list) )cycle
               call mpp_clock_begin(pack_clock)
               to_pe = d_comm%cto_pe(list)
               pos = buffer_pos
               do l=1,l_size  ! loop over number of fields
                  ptr_fieldx = f_addrsx(l)
                  ptr_fieldy = f_addrsy(l)
                  do m=1,8
                     if( d_comm%do_thisS4(m,list) )then
!                        call mpp_clock_begin(pack_loop_clock)
                        i = d_comm%send(m,list)%i2; j = d_comm%send(m,list)%j2
                        do k = 1,ke
                           pos = pos + 2
                           buffer(pos-1) = fieldx(i,j,k)
                           buffer(pos)   = fieldy(i,j,k)
                        end do
!                        call mpp_clock_end(pack_loop_clock)
                     end if
                  end do  ! do m=1,8
               end do  ! do l=1,l_size
               call mpp_clock_end(pack_clock)
               call mpp_clock_begin(send_clock)
               msgsize = pos - buffer_pos
               if( msgsize.GT.0 )then
                  call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
                  buffer_pos = pos
               end if
               call mpp_clock_end(send_clock)
               !     write(stdout(),*) 'Update send checksum=',mpp_chksum(rbuffer(buffer_pos))
            end do

            !recv
            nlist = d_comm%Rlist_size
            do list = 0,nlist-1
               if( .NOT. d_comm%R_do_buf4(list) )cycle
               call mpp_clock_begin(recv_clock)
               from_pe = d_comm%cfrom_pe(list)
               msgsize = 0
               do m=1,8
                  if( d_comm%do_thisR4(m,list) ) msgsize = msgsize + 2*ke
               end do
               msgsize = msgsize*l_size

               if( msgsize.GT.0 )then
                  call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
                  buffer_pos = buffer_pos + msgsize
               end if
               call mpp_clock_end(recv_clock)
            end do

            !unpack recv
            !unpack halos in reverse order

            do list = nlist-1,0,-1
               if( .NOT.d_comm%R_do_buf4(list) )cycle
               call mpp_clock_begin(unpk_clock)
               pos = buffer_pos

               do l=l_size,1,-1  ! loop over number of fields
                  ptr_fieldx = f_addrsx(l)
                  ptr_fieldy = f_addrsy(l)

                  do m=8,1,-1
                     if( d_commy%do_thisR4(m,list) )then  
                        msgsize = 2*ke
                        pos = buffer_pos - msgsize
                        buffer_pos = pos
                        i = d_commy%recv(m,list)%i2; j = d_commy%recv(m,list)%j2
                        do k = 1,ke
                           pos = pos + 2
                           fieldx(i,j,k) = buffer(pos-1)
                           fieldy(i,j,k) = buffer(pos)
                        end do
                     end if
                  end do  ! do m=8,1,-1
               end do  ! do l=l_size,1,-1

               call mpp_clock_end(unpk_clock)
               !        write(stdout(),*) 'Update field checksum=',mpp_chksum(rfield)
            end do
         end if
      end if  ! end if ( ncontacts > 0 )

     ! ---northern boundary fold
      if( BTEST(domain%fold,NORTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then
         if( BTEST(update_flags,NORTH) ) then
            js = domain%y%global%end + 1
            je = domain%y%data%end
            if( je.GE.js )then
               is = domain%x%data%begin
               ie = domain%x%data%end
               select case(grid_offset_type)
               case (AGRID) 
                  do l=1,l_size
                     ptr_fieldx = f_addrsx(l)
                     ptr_fieldy = f_addrsy(l)   
                     do k = 1,ke
                        do j = js,je
                           do i = is,ie
                              fieldx(i,j,k) = -fieldx(i,j,k)
                              fieldy(i,j,k) = -fieldy(i,j,k)
                           end do
                        end do
                     end do
                  end do
               case (BGRID_NE) 
                  if(domain%symmetry) then
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                           do j = js+1,je+1
                              do i = is,ie+1
                                 fieldx(i,j,k) = -fieldx(i,j,k)
                                 fieldy(i,j,k) = -fieldy(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  else
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                           do j = js,je
                              do i = is,ie
                                 fieldx(i,j,k) = -fieldx(i,j,k)
                                 fieldy(i,j,k) = -fieldy(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  end if
               case (CGRID_NE)
                  if(domain%symmetry) then
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)   
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
                     end do
                  else
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                           do j = js,je
                              do i = is,ie
                                 fieldx(i,j,k) = -fieldx(i,j,k)
                                 fieldy(i,j,k) = -fieldy(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  endif
               end select
            endif
         end if
         !--- flip the sign at the folded north edge.
         j = domain%y%global%end
         if( domain%y%data%begin.LE.j .AND. j.LE.domain%y%data%end )then !fold is within domain
            if ( domain%x%pos.GE.size(domain%x%list(:))/2 ) then
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
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                           do i = is, ie
                              fieldx(i,j,k) = -fieldx(i,j,k) 
                              fieldy(i,j,k) = -fieldy(i,j,k) 
                           end do
                        end do
                     end do
                  case(CGRID_NE)
                     if(domain%symmetry) j = j+1
                     do l=1,l_size
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                           do i = is, ie
                              fieldy(i,j,k) = -fieldy(i,j,k) 
                           end do
                        end do
                     end do
                  end select
                  if( do_flip2 ) then
                     select case(grid_offset_type)
                     case(BGRID_NE)
                        if(domain%symmetry) then
                           ie2 = ie2 + 1
                        end if
                        do l=1,l_size
                           ptr_fieldx = f_addrsx(l)
                           ptr_fieldy = f_addrsy(l)   
                           do k = 1,ke
                              do i = is2, ie2
                                 fieldx(i,j,k) = -fieldx(i,j,k) 
                                 fieldy(i,j,k) = -fieldy(i,j,k) 
                              end do
                           end do
                        end do
                     case(CGRID_NE)
                        do l=1,l_size
                           ptr_fieldy = f_addrsy(l)   
                           do k = 1,ke
                              do i = is2, ie2
                                 fieldy(i,j,k) = -fieldy(i,j,k) 
                              end do
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
                     if( domain%x%data%begin.LE.i .AND. i.LE.domain%x%data%end+1 )then
                        do l=1,l_size
                           ptr_fieldx = f_addrsx(l)
                           ptr_fieldy = f_addrsy(l)   
                           do k = 1,ke
                              fieldx(i,j,k) = 0.
                              fieldy(i,j,k) = 0.
                           end do
                        end do
                     end if
                  end do
               else
                  do i = domain%x%global%begin-1,domain%x%global%end,(domain%x%global%begin+domain%x%global%end)/2
                     if( domain%x%data%begin.LE.i .AND. i.LE.domain%x%data%end )then
                        do l=1,l_size
                           ptr_fieldx = f_addrsx(l)
                           ptr_fieldy = f_addrsy(l)   
                           do k = 1,ke
                              fieldx(i,j,k) = 0.
                              fieldy(i,j,k) = 0.
                           end do
                        end do
                     end if
                  end do
               endif
            end if

            ! the last code code block correct an error where the data in your halo coming from 
            ! other half may have the wrong sign

            !right of midpoint, when update north and east direction.
            if ( BTEST(update_flags,NORTH) .OR. BTEST(update_flags,EAST) ) then
               j = domain%y%global%end 
               is = (domain%x%global%begin+domain%x%global%end)/2

               if( domain%x%compute%begin.LE.is .AND. is.LT.domain%x%data%end &
                   .AND. domain%x%pos < size(domain%x%list(:))/2 )then
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
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                           do i = is,ie
                              fieldx(i,j,k) = -fieldx(i,j,k)
                              fieldy(i,j,k) = -fieldy(i,j,k)
                           end do
                        end do
                     end do
                  case(CGRID_NE)
                     do l=1,l_size
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                           do i = is+1,ie
                              fieldy(i,j,k) = -fieldy(2*is-i+1,j,k)
                           end do
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

    end subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_V_
