    subroutine MPP_UPDATE_DOMAINS_2D_( field, domain, flags )
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:)
      type(domain2D), intent(in) :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call mpp_update_domains( field3D, domain, flags )
      return
    end subroutine MPP_UPDATE_DOMAINS_2D_

    subroutine MPP_UPDATE_DOMAINS_4D_( field, domain, flags )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:)
      type(domain2D), intent(in) :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call mpp_update_domains( field3D, domain, flags )
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_

    subroutine MPP_UPDATE_DOMAINS_5D_( field, domain, flags )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:,:)
      type(domain2D), intent(in) :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call mpp_update_domains( field3D, domain, flags )
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_

    subroutine MPP_UPDATE_DOMAINS_3D_( field, domain, flags )
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: field(domain%x%data%begin:,domain%y%data%begin:,:)
      integer, intent(in), optional :: flags
      integer :: update_flags
!equate to mpp_domains_stack
      integer :: wordlen        !#words of MPP_TYPE_ fit in 1 word of mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack))
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


      n = size(domain%list)
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
             call mpp_send( buffer(buffer_pos+1), msgsize, to_pe )
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
             call mpp_recv( buffer(buffer_pos+1), msgsize, from_pe )
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
      call mpp_sync_self( domain%list(:)%pe )
      call mpp_clock_end(wait_clock)
      return
    end subroutine MPP_UPDATE_DOMAINS_3D_

    subroutine MPP_REDISTRIBUTE_2D_( domain_in, field_in, domain_out, field_out )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:)
      MPP_TYPE_ :: field3D_in (size(field_in, 1),size(field_in, 2),1)
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),1)
#ifdef use_CRI_pointers
      pointer( ptr_in,  field3D_in  )
      pointer( ptr_out, field3D_out )
      ptr_in  = LOC(field_in )
      ptr_out = LOC(field_out)
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out )
#else
      call mpp_error( FATAL, 'Requires Cray pointers.' )
#endif
      return
    end subroutine MPP_REDISTRIBUTE_2D_

    subroutine MPP_REDISTRIBUTE_3D_( domain_in, field_in, domain_out, field_out )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in ( domain_in%x%data%begin:, domain_in%y%data%begin:,:)
      MPP_TYPE_, intent(out) :: field_out(domain_out%x%data%begin:,domain_out%y%data%begin:,:)
      integer :: is, ie, js, je, ke, isc, iec, jsc, jec
      integer :: i, j, k
      integer :: list, m, n, pos, msgsize
      integer :: to_pe, from_pe
!      MPP_TYPE_, dimension(domain_in%x%compute%size*domain_in%y%compute%size*size(field_in,3)) :: send_buf, recv_buf
      MPP_TYPE_ :: buffer(size(mpp_domains_stack))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif
      integer :: buffer_pos, wordlen
      character(len=8) :: text

      ke = size(field_in,3)
      if( ke.NE.size(field_out,3) )call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mismatch between field_in and field_out.' )
      if( UBOUND(field_in,1).NE.domain_in%x%data%end .OR. &
          UBOUND(field_in,2).NE.domain_in%y%data%end ) &
          call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_in must be on data domain of domain_in.' )
      if( UBOUND(field_out,1).NE.domain_out%x%data%end .OR. &
          UBOUND(field_out,2).NE.domain_out%y%data%end ) &
          call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_out must be on data domain of domain_out.' )
      buffer_pos = 0
#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
      wordlen = size(TRANSFER(buffer(1),mpp_domains_stack))
#endif
!send
      call mpp_get_compute_domain( domain_in, isc, iec, jsc, jec )
      n = size(domain_out%list)
      do list = 0,n-1
         m = mod( domain_out%pos+list, n )
         to_pe = domain_out%list(m)%pe
         call mpp_get_compute_domain( domain_out%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             msgsize = (ie-is+1)*(je-js+1)*ke
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             pos = buffer_pos
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos+1
                      buffer(pos) = field_in(i,j,k)
                   end do
                end do
             end do
             call mpp_send( buffer(buffer_pos+1), msgsize, to_pe )
             buffer_pos = pos
         end if
      end do
!recv
      call mpp_get_compute_domain( domain_out, isc, iec, jsc, jec )
      n = size(domain_in%list)
      do list = 0,n-1
         m = mod( domain_in%pos+n-list, n )
         from_pe = domain_in%list(m)%pe
         call mpp_get_compute_domain( domain_in%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             msgsize = (ie-is+1)*(je-js+1)*ke
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             call mpp_recv( buffer(buffer_pos+1), msgsize, from_pe )
             pos = buffer_pos
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos+1
                      field_out(i,j,k) = buffer(pos)
                   end do
                end do
             end do
             buffer_pos = pos
         end if
      end do

      call mpp_sync_self( domain_out%list(:)%pe )
      return
    end subroutine MPP_REDISTRIBUTE_3D_

    subroutine MPP_REDISTRIBUTE_4D_( domain_in, field_in, domain_out, field_out )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:,:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:,:,:)
      MPP_TYPE_ :: field3D_in (size(field_in, 1),size(field_in, 2),size(field_in ,3)*size(field_in ,4))
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),size(field_out,3)*size(field_out,4))
#ifdef use_CRI_pointers
      pointer( ptr_in,  field3D_in  )
      pointer( ptr_out, field3D_out )
      ptr_in  = LOC(field_in )
      ptr_out = LOC(field_out)
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out )
#else
      call mpp_error( FATAL, 'Requires Cray pointers.' )
#endif
      return
    end subroutine MPP_REDISTRIBUTE_4D_

    subroutine MPP_REDISTRIBUTE_5D_( domain_in, field_in, domain_out, field_out )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:,:,:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:,:,:,:)
      MPP_TYPE_ :: field3D_in (size(field_in, 1),size(field_in, 2),size(field_in ,3)*size(field_in ,4)*size(field_in ,5))
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),size(field_out,3)*size(field_out,4)*size(field_out,5))
#ifdef use_CRI_pointers
      pointer( ptr_in,  field3D_in  )
      pointer( ptr_out, field3D_out )
      ptr_in  = LOC(field_in )
      ptr_out = LOC(field_out)
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out )
#else
      call mpp_error( FATAL, 'Requires Cray pointers.' )
#endif
      return
    end subroutine MPP_REDISTRIBUTE_5D_
#ifdef VECTOR_FIELD_
!VECTOR_FIELD_ is set to false for MPP_TYPE_ integer or logical.
!vector fields
    subroutine MPP_UPDATE_DOMAINS_2D_V_( fieldx, fieldy, domain, flags, gridtype )
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_, intent(inout), dimension(:,:) :: fieldx, fieldy
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)
#ifdef use_CRI_pointers
      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: requires Cray pointers.' )
#endif
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
      return
    end subroutine MPP_UPDATE_DOMAINS_2D_V_

    subroutine MPP_UPDATE_DOMAINS_3D_V_( fieldx, fieldy, domain, flags, gridtype )
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(inout) :: domain
      MPP_TYPE_, intent(inout), dimension(domain%x%data%begin:,domain%y%data%begin:,:) :: fieldx, fieldy
      integer, intent(in), optional :: flags, gridtype
      integer :: update_flags
      integer :: i,j,k,n, is, ie, js, je, ke, pos
      integer :: isg, ieg, jsg, jeg, ioff, joff
      MPP_TYPE_ :: buffer(size(mpp_domains_stack))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif
      integer :: buffer_pos, msgsize, wordlen
      character(len=8) :: text

!gridtype
      grid_offset_type = AGRID
      if( PRESENT(gridtype) )then
          if( gridtype.NE.AGRID .AND. &
              gridtype.NE.BGRID_NE .AND. gridtype.NE.BGRID_SW .AND. &
              gridtype.NE.CGRID_NE .AND. gridtype.NE.CGRID_SW ) &
               call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: gridtype must be one of AGRID|BGRID_NE|BGRID_SW|CGRID_NE|CGRID_SW.' )
          grid_offset_type = gridtype
          call compute_overlaps(domain)
          if( grid_offset_type.NE.domain%gridtype ) &
               call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: gridtype cannot be changed during run.' )
      end if
!grid_offset_type used by update domains to determine shifts.
      call mpp_update_domains( fieldx, domain, flags )
      call mpp_update_domains( fieldy, domain, flags )
#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
#endif
      wordlen = size(TRANSFER(buffer(1),mpp_domains_stack))
      buffer_pos = 0
!for all gridtypes
      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) )update_flags = flags
      ke = size(fieldx,3)
      call mpp_get_global_domain( domain, is, ie, js, je, xsize=ioff, ysize=joff )
!if vector data was moved over a northern boundary fold, then flip sign at the fold
!need to add code for EWS boundaries
      if( BTEST(domain%fold,NORTH) .AND. BTEST(update_flags,NORTH) )then
          js = max(domain%y%global%end+1,domain%y%compute%end+1)
          je = domain%y%data%end
          if( je.GE.js )then
!for BGRID_NE, we need to move all data leftward by 1 point (vertical shift already done by mpp_update_domains)
              if( grid_offset_type.EQ.BGRID_NE )then
                  pos = domain%x%pos - 1 !the one on your right
                  if( pos.GE.0 )then
                      is = domain%x%list(pos)%data%end+1; ie=is
                  else if( domain%x%cyclic )then
                      pos = pos + size(domain%x%list)
                      is = domain%x%list(pos)%data%end+1 - ioff; ie=is
                  else
                      is=1; ie=0
                  end if
                  n = buffer_pos
                  if( ie.GE.is .AND. je.GE.js )then
                      msgsize = (je-js+1)*ke
                      mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
                      if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                          write( text,'(i8)' )mpp_domains_stack_hwm
                          call mpp_error( FATAL, 'MPP_UPDATE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                               //trim(text)//') from all PEs.' )
                      end if
                      do k = 1,ke
                         do j = js,je
                            do i = is,ie
                               n = n + 2
                               buffer(n-1) = fieldx(i,j,k)
                               buffer(n  ) = fieldy(i,j,k)
                            end do
                         end do
                      end do
                      call mpp_send( buffer(buffer_pos+1), n, domain%x%list(pos)%pe )
                      buffer_pos = buffer_pos + n
                  end if
!shift and flip the sign in your own halo
                  is = domain%x%data%begin
                  ie = domain%x%data%end-1
                  do k = 1,ke
                     do j = js,je
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i+1,j,k)
                           fieldy(i,j,k) = -fieldy(i+1,j,k)
                        end do
                     end do
                  end do
!receive data at x%data%end
                  pos = domain%x%pos + 1
                  if( pos.LT.size(domain%x%list) )then
                      n = (je-js+1)*ke*2
                  else if( domain%x%cyclic )then
                      pos = pos - size(domain%x%list)
                      n = (je-js+1)*ke*2
                  else
                      n = 0
                  end if
                  if( n.GT.0 )then
                      mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+n)*wordlen )
                      if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                          write( text,'(i8)' )mpp_domains_stack_hwm
                          call mpp_error( FATAL, 'MPP_UPDATE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                               //trim(text)//') from all PEs.' )
                      end if
                      call mpp_recv( buffer(buffer_pos+1), n, domain%x%list(pos)%pe )
                      i = domain%x%data%end
                      n = buffer_pos
                      do k = 1,ke
                         do j = js,je
                            n = n + 2
                            fieldx(i,j,k) = -buffer(n-1)
                            fieldy(i,j,k) = -buffer(n  )
                         end do
                      end do
                  end if
              else
!flip the sign
                  is = domain%x%data%begin
                  ie = domain%x%data%end
                  do k = 1,ke
                     do j = js,je
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
              end if
          end if
      end if

      grid_offset_type = AGRID  !reset
      call mpp_sync_self( domain%list(:)%pe )
      return
    end subroutine MPP_UPDATE_DOMAINS_3D_V_

    subroutine MPP_UPDATE_DOMAINS_4D_V_( fieldx, fieldy, domain, flags, gridtype )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout), dimension(:,:,:,:) :: fieldx, fieldy
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
#ifdef use_CRI_pointers
      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: requires Cray pointers.' )
#endif
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_V_

    subroutine MPP_UPDATE_DOMAINS_5D_V_( fieldx, fieldy, domain, flags, gridtype )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout), dimension(:,:,:,:,:) :: fieldx, fieldy
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4)*size(fieldx,5))
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4)*size(fieldy,5))
#ifdef use_CRI_pointers
      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: requires Cray pointers.' )
#endif
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_V_
#endif VECTOR_FIELD_
