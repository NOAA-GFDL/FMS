    subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_( f_addrs, d_comm, d_type )
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND), intent(in)         :: f_addrs(:)
      type(DomainCommunicator2D), intent(in) :: d_comm
      MPP_TYPE_, intent(in)                  :: d_type  ! creates unique interface

      MPP_TYPE_ :: field(d_comm%domain%x%data%begin:d_comm%domain%x%data%begin+d_comm%isize-1, &
                         d_comm%domain%y%data%begin:d_comm%domain%y%data%begin+d_comm%jsize-1,d_comm%ke)
      pointer(ptr_field, field)
!     real(kind(field)) :: rfield(d_comm%domain%x%data%begin:d_comm%domain%x%data%begin+d_comm%isize-1, &
!                                 d_comm%domain%y%data%begin:d_comm%domain%y%data%begin+d_comm%jsize-1,d_comm%ke)
!     pointer (ptr_rfield, rfield)
      integer :: update_flags
!equate to mpp_domains_stack
      integer :: wordlen        !#words of MPP_TYPE_ fit in 1 word of mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer(ptr,buffer )
      integer :: buffer_pos
      real(kind(buffer)) :: rbuffer(size(mpp_domains_stack(:)))
      pointer(r_ptr,rbuffer)
      integer :: i, j, k, m, n, l, l_size
      integer :: is, ie, js, je, ke, isize, jsize,istat
!receive domains saved here for unpacking
!for non-blocking version, could be recomputed
      integer :: to_pe, from_pe, list, pos, msgsize
      character(len=8) :: text



      n = d_comm%Rlist_size
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
      ptr = LOC(mpp_domains_stack)
      wordlen = size(transfer(buffer(1),mpp_domains_stack))
      l_size = size(f_addrs(:))
      ke = d_comm%ke

!     r_ptr = LOC(buffer)

!send
      do list = 0,n-1
         if( .NOT.d_comm%S_do_buf(list) )cycle
!        if( d_comm%S_do_buf(list) )then
         call mpp_clock_begin(pack_clock)
         to_pe = d_comm%cto_pe(list)
         pos = buffer_pos
         do l=1,l_size  ! loop over number of fields
           ptr_field = f_addrs(l)
           do m=1,8
             if( d_comm%do_thisS(m,list) )then
               is=d_comm%sendis(m,list); ie=d_comm%sendie(m,list)
               js=d_comm%sendjs(m,list); je=d_comm%sendje(m,list)
               call mpp_clock_begin(pack_loop_clock)
               do k = 1,ke
                  do j = js,je
                     do i = is,ie
                        pos = pos + 1
                        buffer(pos) = field(i,j,k)
                     end do
                  end do
               end do
               call mpp_clock_end(pack_loop_clock)
             end if
           end do  ! do m=1,8
         end do  ! do l=1,l_size
         call mpp_clock_end(pack_clock)
         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
             buffer_pos = pos
         end if
         call mpp_clock_end(send_clock)
!        endif
!     write(stdout(),*) 'Update send checksum=',mpp_chksum(rbuffer(buffer_pos))
      end do
             
!recv
      n = d_comm%Rlist_size
      do list = 0,n-1
         if( .NOT. d_comm%R_do_buf(list) )cycle
         call mpp_clock_begin(recv_clock)
         from_pe = d_comm%cfrom_pe(list)
         msgsize = 0
         do m=1,8
           if( d_comm%do_thisR(1,m,list) )then
             is=d_comm%recvis(m,list); ie=d_comm%recvie(m,list)
             js=d_comm%recvjs(m,list); je=d_comm%recvje(m,list)
             msgsize = msgsize + d_comm%R_msize(m,list)
           end if
         end do
         msgsize = msgsize * l_size
         if( msgsize.GT.0 )then
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
             buffer_pos = buffer_pos + msgsize
         end if
         call mpp_clock_end(recv_clock)
      end do
             
!     write(stdout(),*) 'Update recv checksum=',mpp_chksum(rbuffer(buffer_pos))
             
!unpack recv
!unpack halos in reverse order
!     ptr_rfield = f_addrs(1)
      do list = n-1,0,-1
         if( .NOT.d_comm%R_do_buf(list) )cycle
!        if( d_comm%R_do_buf(list) )then
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         do l=l_size,1,-1  ! loop over number of fields
           ptr_field = f_addrs(l)
           do m=8,1,-1
             if( d_comm%do_thisR(1,m,list) )then  ! direction
               is=d_comm%recvis(m,list); ie=d_comm%recvie(m,list)
               js=d_comm%recvjs(m,list); je=d_comm%recvje(m,list)
               msgsize = d_comm%R_msize(m,list)
               pos = buffer_pos - msgsize
               buffer_pos = pos
               if( d_comm%do_thisR(2,m,list) )then  ! folded?
                 do k = 1,ke
                    do j = je,js,-1
                       do i = ie,is,-1
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
               else
                 do k = 1,ke
                    do j = js,je
                       do i = is,ie
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
               end if
             end if
           end do  ! do m=8,1,-1
         end do  ! do l=l_size,1,-1
         call mpp_clock_end(unpk_clock)
!        endif
!        write(stdout(),*) 'Update field checksum=',mpp_chksum(rfield)
      end do

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)
      return
    end subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_
