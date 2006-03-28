    subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_( f_addrs, d_comm, d_type)
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND), intent(in)         :: f_addrs(:)
      type(DomainCommunicator2D),    pointer :: d_comm
      MPP_TYPE_, intent(in)                  :: d_type  ! creates unique interface

      MPP_TYPE_ :: field(d_comm%domain%x%data%begin:d_comm%domain%x%data%begin+d_comm%isize-1, &
                         d_comm%domain%y%data%begin:d_comm%domain%y%data%begin+d_comm%jsize-1,d_comm%ke)
      pointer(ptr_field, field)
!     real(kind(field)) :: rfield(d_comm%domain%x%data%begin:d_comm%domain%x%data%begin+d_comm%isize-1, &
!                                 d_comm%domain%y%data%begin:d_comm%domain%y%data%begin+d_comm%jsize-1,d_comm%ke)
!     pointer (ptr_rfield, rfield)

!equate to mpp_domains_stack
      integer :: wordlen        !#words of MPP_TYPE_ fit in 1 word of mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer(ptr,buffer )
      integer :: buffer_pos
!      real(kind(buffer)) :: rbuffer(size(mpp_domains_stack(:)))
!      pointer(r_ptr,rbuffer)

!receive domains saved here for unpacking
!for non-blocking version, could be recomputed
      integer :: to_pe, from_pe, list, pos, msgsize
      integer :: n, l_size, l, ke, m, i, j, k, nlist, nrecv, nsend
      integer :: is, ie, js, je

      nlist = d_comm%Rlist_size
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
      ptr = LOC(mpp_domains_stack)
      wordlen = size(transfer(buffer(1),mpp_domains_stack))
      l_size = size(f_addrs(:))
      ke = d_comm%ke

!     r_ptr = LOC(buffer)

!send
      do list = 0,nlist-1
         if( .NOT.d_comm%S_do_buf(list) )cycle
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos
         do l=1,l_size  ! loop over number of fields
           ptr_field = f_addrs(l)
           do m=1,8
              if( d_comm%do_thisS(m,list) )then
!                 call mpp_clock_begin(pack_loop_clock)
                 is = d_comm%send(m,list)%is; ie = d_comm%send(m,list)%ie
                 js = d_comm%send(m,list)%js; je = d_comm%send(m,list)%je
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
!                 call mpp_clock_end(pack_loop_clock)
              end if
              if( d_comm%do_thisS2(m,list) )then
!                 call mpp_clock_begin(pack_loop_clock)
                 nsend = d_comm%send(m,list)%n
                 do n = 1, nsend
                    i = d_comm%send(m,list)%i(n)
                    j = d_comm%send(m,list)%j(n)
                    do k = 1,ke
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
!                 call mpp_clock_end(pack_loop_clock)
              end if
           end do  ! do m=1,8
         end do  ! do l=1,l_size
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
      do list = 0,nlist-1
         if( .NOT. d_comm%R_do_buf(list) )cycle
         call mpp_clock_begin(recv_clock)
         msgsize = 0
         do m=1,8
            msgsize = msgsize + d_comm%R_msize(m,list)
         end do
         msgsize = msgsize*l_size
         if( msgsize.GT.0 )then
             from_pe = d_comm%cfrom_pe(list)
             call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
             buffer_pos = buffer_pos + msgsize
         end if
         call mpp_clock_end(recv_clock)
      end do
             
!     write(stdout(),*) 'Update recv checksum=',mpp_chksum(rbuffer(buffer_pos))
             
!unpack recv
!unpack halos in reverse order
!     ptr_rfield = f_addrs(1)

      do list = nlist-1,0,-1
         if( .NOT.d_comm%R_do_buf(list) )cycle
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         do l=l_size,1,-1  ! loop over number of fields
           ptr_field = f_addrs(l)
           do m=8,1,-1
              if( d_comm%do_thisR2(m,list) )then  ! direction
                 nrecv = d_comm%recv(m,list)%n
                 msgsize = nrecv*ke
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 do n = 1, nrecv
                    i = d_comm%recv(m,list)%i(n)
                    j = d_comm%recv(m,list)%j(n)
                    do k = 1,ke
                       pos = pos + 1
                       field(i,j,k) = buffer(pos)
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
           end do  ! do m=8,1,-1
         end do  ! do l=l_size,1,-1
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
!send
      if ( d_comm%domain%ncontacts > 0 ) then
         do list = 0,nlist-1
            if( .NOT.d_comm%S_do_buf3(list) )cycle
            call mpp_clock_begin(pack_clock)
            to_pe = d_comm%cto_pe(list)
            pos = buffer_pos
            do l=1,l_size  ! loop over number of fields
               ptr_field = f_addrs(l)
               do m=1,8
                  if( d_comm%do_thisS3(m,list) )then
!                     call mpp_clock_begin(pack_loop_clock)
                     is = d_comm%send(m,list)%is; ie = d_comm%send(m,list)%ie
                     js = d_comm%send(m,list)%js; je = d_comm%send(m,list)%je
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = field(i,j,k)
                           end do
                        end do
                     end do
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
         do list = 0,nlist-1
            if( .NOT. d_comm%R_do_buf3(list) )cycle
            call mpp_clock_begin(recv_clock)
            from_pe = d_comm%cfrom_pe(list)
            msgsize = 0
            do l=1,l_size  ! loop over number of fields
               do m=1,8
                  if( d_comm%do_thisR3(m,list) )then
                     msgsize = msgsize + d_comm%R_msize2(m,list)
                  end if
               end do
            enddo
            if( msgsize.GT.0 )then
               call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
               buffer_pos = buffer_pos + msgsize
            end if
            call mpp_clock_end(recv_clock)
         end do

         !     write(stdout(),*) 'Update recv checksum=',mpp_chksum(rbuffer(buffer_pos))

         !unpack recv
         !unpack halos in reverse order
         !     ptr_rfield = f_addrs(1)

         do list = nlist-1,0,-1
            if( .NOT.d_comm%R_do_buf3(list) )cycle
            call mpp_clock_begin(unpk_clock)
            pos = buffer_pos
            do l=l_size,1,-1  ! loop over number of fields
               ptr_field = f_addrs(l)
               do m=8,1,-1
                  if( d_comm%do_thisR3(m,list) )then  ! direction
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
                  end if
               end do  ! do m=8,1,-1
            end do  ! do l=l_size,1,-1
            call mpp_clock_end(unpk_clock)
            !        write(stdout(),*) 'Update field checksum=',mpp_chksum(rfield)
         end do

         !--- for the data on corner, some extra communication is needed to ensure the corner value on different tile is the same.
         if( d_comm%position == CORNER ) then
            do list = 0,nlist-1
               if( .NOT.d_comm%S_do_buf4(list) )cycle
               call mpp_clock_begin(pack_clock)
               to_pe = d_comm%cto_pe(list)
               pos = buffer_pos
               do l=1,l_size  ! loop over number of fields
                  ptr_field = f_addrs(l)
                  do m=1,8
                     if( d_comm%do_thisS4(m,list) )then
!                        call mpp_clock_begin(pack_loop_clock)
                        i = d_comm%send(m,list)%i2; j = d_comm%send(m,list)%j2
                        do k = 1,ke
                           pos = pos + 1
                           buffer(pos) = field(i,j,k)
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
            do list = 0,nlist-1
               if( .NOT. d_comm%R_do_buf4(list) )cycle
               call mpp_clock_begin(recv_clock)
               from_pe = d_comm%cfrom_pe(list)
               msgsize = 0
               do l=1,l_size  ! loop over number of fields
                  do m=1,8
                     if( d_comm%do_thisR4(m,list) )then
                        msgsize = msgsize + ke
                     end if
                  end do
               enddo
               if( msgsize.GT.0 )then
                  call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
                  buffer_pos = buffer_pos + msgsize
               end if
               call mpp_clock_end(recv_clock)
            end do

            !unpack recv
            !unpack halos in reverse order
            !     ptr_rfield = f_addrs(1)

            do list = nlist-1,0,-1
               if( .NOT.d_comm%R_do_buf4(list) )cycle
               call mpp_clock_begin(unpk_clock)
               pos = buffer_pos
               do l=l_size,1,-1  ! loop over number of fields
                  ptr_field = f_addrs(l)
                  do m=8,1,-1
                     if( d_comm%do_thisR4(m,list) )then  ! direction
                        msgsize = ke
                        pos = buffer_pos - msgsize
                        buffer_pos = pos
                        i = d_comm%recv(m,list)%i2; j = d_comm%recv(m,list)%j2
                        do k = 1,ke
                           pos = pos + 1
                           field(i,j,k) = buffer(pos)
                        end do
                     end if
                  end do  ! do m=8,1,-1
               end do  ! do l=l_size,1,-1
               call mpp_clock_end(unpk_clock)
               !        write(stdout(),*) 'Update field checksum=',mpp_chksum(rfield)
            end do
         end if
      end if  ! if( domain%ncontacts > 1 )

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)
      return
    end subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_
