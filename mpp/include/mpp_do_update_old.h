    subroutine MPP_DOMAINS_DO_UPDATE_3Dold_( field, domain, flags )
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: field(domain%x%data%begin:,domain%y%data%begin:,:)
      integer, intent(in), optional :: flags
      integer :: update_flags
!equate to mpp_domains_stack
      integer :: wordlen        !#words of MPP_TYPE_ fit in 1 word of mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif
      integer :: buffer_pos
      integer :: i, j, k, m, n
      integer :: is, ie, js, je, ke
!receive domains saved here for unpacking
!for non-blocking version, could be recomputed
      integer, dimension(8) :: isr, ier, jsr, jer
      integer :: to_pe, from_pe, list, pos, msgsize
      logical :: recv_e, recv_se, recv_s, recv_sw, recv_w, recv_nw, recv_n, recv_ne
      logical :: send_e, send_se, send_s, send_sw, send_w, send_nw, send_n, send_ne
      logical :: folded
      character(len=8) :: text

      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) )update_flags = flags
      recv_w = BTEST(update_flags,WEST)
      recv_e = BTEST(update_flags,EAST)
      recv_s = BTEST(update_flags,SOUTH)
      recv_n = BTEST(update_flags,NORTH)
      recv_ne = recv_e .AND. recv_n
      recv_se = recv_e .AND. recv_s
      recv_sw = recv_w .AND. recv_s
      recv_nw = recv_w .AND. recv_n
      send_w = recv_e
      send_e = recv_w
      send_s = recv_n
      send_n = recv_s
      send_ne = send_e .AND. send_n
      send_se = send_e .AND. send_s
      send_sw = send_w .AND. send_s
      send_nw = send_w .AND. send_n
      if( recv_w .AND. BTEST(domain%fold,WEST)  .AND. BTEST(grid_offset_type,EAST)  ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
      if( recv_s .AND. BTEST(domain%fold,SOUTH) .AND. BTEST(grid_offset_type,NORTH) ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
      if( recv_e .AND. BTEST(domain%fold,EAST)  .AND. BTEST(grid_offset_type,WEST)  ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
      if( recv_n .AND. BTEST(domain%fold,NORTH) .AND. BTEST(grid_offset_type,SOUTH) ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )


      n = size(domain%list(:))
      ke = size(field,3)
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
      wordlen = size(transfer(buffer(1),mpp_domains_stack))
#endif
!send
      do list = 0,n-1
         m = mod( domain%pos+list, n )
         if( .NOT.domain%list(m)%overlap )cycle
         call mpp_clock_begin(pack_clock)
         to_pe = domain%list(m)%pe
         pos = buffer_pos
         if( send_w .AND. domain%list(m)%send_w%overlap )then
             is = domain%list(m)%send_w%is; ie = domain%list(m)%send_w%ie
             js = domain%list(m)%send_w%js; je = domain%list(m)%send_w%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_w_off%is; ie = domain%list(m)%send_w_off%ie
                 js = domain%list(m)%send_w_off%js; je = domain%list(m)%send_w_off%je
             end if
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
         if( send_nw .AND. domain%list(m)%send_nw%overlap )then
             is = domain%list(m)%send_nw%is; ie = domain%list(m)%send_nw%ie
             js = domain%list(m)%send_nw%js; je = domain%list(m)%send_nw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_nw_off%is; ie = domain%list(m)%send_nw_off%ie
                 js = domain%list(m)%send_nw_off%js; je = domain%list(m)%send_nw_off%je
             end if
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
         if( send_n .AND. domain%list(m)%send_n%overlap )then
             is = domain%list(m)%send_n%is; ie = domain%list(m)%send_n%ie
             js = domain%list(m)%send_n%js; je = domain%list(m)%send_n%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_n_off%is; ie = domain%list(m)%send_n_off%ie
                 js = domain%list(m)%send_n_off%js; je = domain%list(m)%send_n_off%je
             end if
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
         if( send_ne .AND. domain%list(m)%send_ne%overlap )then
             is = domain%list(m)%send_ne%is; ie = domain%list(m)%send_ne%ie
             js = domain%list(m)%send_ne%js; je = domain%list(m)%send_ne%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_ne_off%is; ie = domain%list(m)%send_ne_off%ie
                 js = domain%list(m)%send_ne_off%js; je = domain%list(m)%send_ne_off%je
             end if
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
         if( send_e .AND. domain%list(m)%send_e%overlap )then
             is = domain%list(m)%send_e%is; ie = domain%list(m)%send_e%ie
             js = domain%list(m)%send_e%js; je = domain%list(m)%send_e%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_e_off%is; ie = domain%list(m)%send_e_off%ie
                 js = domain%list(m)%send_e_off%js; je = domain%list(m)%send_e_off%je
             end if
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
         if( send_se .AND. domain%list(m)%send_se%overlap )then
             is = domain%list(m)%send_se%is; ie = domain%list(m)%send_se%ie
             js = domain%list(m)%send_se%js; je = domain%list(m)%send_se%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_se_off%is; ie = domain%list(m)%send_se_off%ie
                 js = domain%list(m)%send_se_off%js; je = domain%list(m)%send_se_off%je
             end if
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
         if( send_s .AND. domain%list(m)%send_s%overlap )then
             is = domain%list(m)%send_s%is; ie = domain%list(m)%send_s%ie
             js = domain%list(m)%send_s%js; je = domain%list(m)%send_s%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_s_off%is; ie = domain%list(m)%send_s_off%ie
                 js = domain%list(m)%send_s_off%js; je = domain%list(m)%send_s_off%je
             end if
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
         if( send_sw .AND. domain%list(m)%send_sw%overlap )then
             is = domain%list(m)%send_sw%is; ie = domain%list(m)%send_sw%ie
             js = domain%list(m)%send_sw%js; je = domain%list(m)%send_sw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_sw_off%is; ie = domain%list(m)%send_sw_off%ie
                 js = domain%list(m)%send_sw_off%js; je = domain%list(m)%send_sw_off%je
             end if
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
      end do
             
!recv
      do list = 0,n-1
         m = mod( domain%pos+n-list, n )
         if( .NOT.domain%list(m)%overlap )cycle
         call mpp_clock_begin(recv_clock)
         from_pe = domain%list(m)%pe
         msgsize = 0
         if( recv_e .AND. domain%list(m)%recv_e%overlap )then
             is = domain%list(m)%recv_e%is; ie = domain%list(m)%recv_e%ie
             js = domain%list(m)%recv_e%js; je = domain%list(m)%recv_e%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_e_off%is; ie = domain%list(m)%recv_e_off%ie
                 js = domain%list(m)%recv_e_off%js; je = domain%list(m)%recv_e_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_se .AND. domain%list(m)%recv_se%overlap )then
             is = domain%list(m)%recv_se%is; ie = domain%list(m)%recv_se%ie
             js = domain%list(m)%recv_se%js; je = domain%list(m)%recv_se%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_se_off%is; ie = domain%list(m)%recv_se_off%ie
                 js = domain%list(m)%recv_se_off%js; je = domain%list(m)%recv_se_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_s .AND. domain%list(m)%recv_s%overlap )then
             is = domain%list(m)%recv_s%is; ie = domain%list(m)%recv_s%ie
             js = domain%list(m)%recv_s%js; je = domain%list(m)%recv_s%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_s_off%is; ie = domain%list(m)%recv_s_off%ie
                 js = domain%list(m)%recv_s_off%js; je = domain%list(m)%recv_s_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_sw .AND. domain%list(m)%recv_sw%overlap )then
             is = domain%list(m)%recv_sw%is; ie = domain%list(m)%recv_sw%ie
             js = domain%list(m)%recv_sw%js; je = domain%list(m)%recv_sw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_sw_off%is; ie = domain%list(m)%recv_sw_off%ie
                 js = domain%list(m)%recv_sw_off%js; je = domain%list(m)%recv_sw_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_w .AND. domain%list(m)%recv_w%overlap )then
             is = domain%list(m)%recv_w%is; ie = domain%list(m)%recv_w%ie
             js = domain%list(m)%recv_w%js; je = domain%list(m)%recv_w%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_w_off%is; ie = domain%list(m)%recv_w_off%ie
                 js = domain%list(m)%recv_w_off%js; je = domain%list(m)%recv_w_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_nw .AND. domain%list(m)%recv_nw%overlap )then
             is = domain%list(m)%recv_nw%is; ie = domain%list(m)%recv_nw%ie
             js = domain%list(m)%recv_nw%js; je = domain%list(m)%recv_nw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_nw_off%is; ie = domain%list(m)%recv_nw_off%ie
                 js = domain%list(m)%recv_nw_off%js; je = domain%list(m)%recv_nw_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_n .AND. domain%list(m)%recv_n%overlap )then
             is = domain%list(m)%recv_n%is; ie = domain%list(m)%recv_n%ie
             js = domain%list(m)%recv_n%js; je = domain%list(m)%recv_n%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_n_off%is; ie = domain%list(m)%recv_n_off%ie
                 js = domain%list(m)%recv_n_off%js; je = domain%list(m)%recv_n_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_ne .AND. domain%list(m)%recv_ne%overlap )then
             is = domain%list(m)%recv_ne%is; ie = domain%list(m)%recv_ne%ie
             js = domain%list(m)%recv_ne%js; je = domain%list(m)%recv_ne%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_ne_off%is; ie = domain%list(m)%recv_ne_off%ie
                 js = domain%list(m)%recv_ne_off%js; je = domain%list(m)%recv_ne_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
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
             
!unpack recv
!unpack halos in reverse order
      do list = n-1,0,-1
         m = mod( domain%pos+n-list, n )
         if( .NOT.domain%list(m)%overlap )cycle
         call mpp_clock_begin(unpk_clock)
         from_pe = domain%list(m)%pe
         pos = buffer_pos
         if( recv_ne .AND. domain%list(m)%recv_ne%overlap )then
             is = domain%list(m)%recv_ne%is; ie = domain%list(m)%recv_ne%ie
             js = domain%list(m)%recv_ne%js; je = domain%list(m)%recv_ne%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_ne_off%is; ie = domain%list(m)%recv_ne_off%ie
                 js = domain%list(m)%recv_ne_off%js; je = domain%list(m)%recv_ne_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_ne%folded )then
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
         if( recv_n .AND. domain%list(m)%recv_n%overlap )then
             is = domain%list(m)%recv_n%is; ie = domain%list(m)%recv_n%ie
             js = domain%list(m)%recv_n%js; je = domain%list(m)%recv_n%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_n_off%is; ie = domain%list(m)%recv_n_off%ie
                 js = domain%list(m)%recv_n_off%js; je = domain%list(m)%recv_n_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_n%folded )then
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
         if( recv_nw .AND. domain%list(m)%recv_nw%overlap )then
             is = domain%list(m)%recv_nw%is; ie = domain%list(m)%recv_nw%ie
             js = domain%list(m)%recv_nw%js; je = domain%list(m)%recv_nw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_nw_off%is; ie = domain%list(m)%recv_nw_off%ie
                 js = domain%list(m)%recv_nw_off%js; je = domain%list(m)%recv_nw_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_nw%folded )then
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
         if( recv_w .AND. domain%list(m)%recv_w%overlap )then
             is = domain%list(m)%recv_w%is; ie = domain%list(m)%recv_w%ie
             js = domain%list(m)%recv_w%js; je = domain%list(m)%recv_w%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_w_off%is; ie = domain%list(m)%recv_w_off%ie
                 js = domain%list(m)%recv_w_off%js; je = domain%list(m)%recv_w_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_w%folded )then
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
         if( recv_sw .AND. domain%list(m)%recv_sw%overlap )then
             is = domain%list(m)%recv_sw%is; ie = domain%list(m)%recv_sw%ie
             js = domain%list(m)%recv_sw%js; je = domain%list(m)%recv_sw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_sw_off%is; ie = domain%list(m)%recv_sw_off%ie
                 js = domain%list(m)%recv_sw_off%js; je = domain%list(m)%recv_sw_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_sw%folded )then
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
         if( recv_s .AND. domain%list(m)%recv_s%overlap )then
             is = domain%list(m)%recv_s%is; ie = domain%list(m)%recv_s%ie
             js = domain%list(m)%recv_s%js; je = domain%list(m)%recv_s%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_s_off%is; ie = domain%list(m)%recv_s_off%ie
                 js = domain%list(m)%recv_s_off%js; je = domain%list(m)%recv_s_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_s%folded )then
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
         if( recv_se .AND. domain%list(m)%recv_se%overlap )then
             is = domain%list(m)%recv_se%is; ie = domain%list(m)%recv_se%ie
             js = domain%list(m)%recv_se%js; je = domain%list(m)%recv_se%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_se_off%is; ie = domain%list(m)%recv_se_off%ie
                 js = domain%list(m)%recv_se_off%js; je = domain%list(m)%recv_se_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_se%folded )then
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
         if( recv_e .AND. domain%list(m)%recv_e%overlap )then
             is = domain%list(m)%recv_e%is; ie = domain%list(m)%recv_e%ie
             js = domain%list(m)%recv_e%js; je = domain%list(m)%recv_e%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_e_off%is; ie = domain%list(m)%recv_e_off%ie
                 js = domain%list(m)%recv_e_off%js; je = domain%list(m)%recv_e_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_e%folded )then
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
         call mpp_clock_end(unpk_clock)
      end do

      call mpp_clock_begin(wait_clock)
!      call mpp_sync_self( domain%list(:)%pe )
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)
      return
    end subroutine MPP_DOMAINS_DO_UPDATE_3Dold_
