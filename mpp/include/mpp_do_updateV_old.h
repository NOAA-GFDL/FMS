! -*-f90-*- 
    subroutine MPP_DOMAINS_DO_UPDATE_3Dold_V_( fieldx, fieldy, domainx, domainy, flags, gridtype, name)
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(in)    :: domainx, domainy
      MPP_TYPE_, intent(inout), dimension(domainx%x(1)%memory%begin:,domainx%y(1)%memory%begin:,:) :: fieldx, fieldy
      integer, intent(in), optional :: flags, gridtype
      character(len=*), intent(in), optional :: name

      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif

      integer :: update_flags, grid_offset_type
      integer :: buffer_pos, nsend, nrecv, ndir
      integer :: sign, midpoint, ioff, joff, dir
      integer :: i, j, k, m, n, is, ie, js, je, ke, nlist, is2, ie2, js2, je2
      integer :: to_pe, from_pe, list, pos, msgsize
      logical :: send(8), recv(8)
      character(len=8) :: text
      logical :: do_flip, do_flip2, flip_n, flip_w, flip_e
      character(len=64) :: field_name      
      type(checkbound),  pointer :: bound_x  => NULL()
      type(checkbound),  pointer :: bound_y  => NULL()
      type(overlapSpec), pointer :: overPtrx => NULL()         
      type(overlapSpec), pointer :: overPtry => NULL()  

      !--- there should be only one domain on each pe for the old update
      if(size(domainx%x(:)) > 1 ) call mpp_error(FATAL, &
          "mpp_do_updateV_old: more than one domain on this pe, should use mpp_do_updateV_new")

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

      nlist = size(domainx%list(:))
      ke = size(fieldx,3)
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking

#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
#endif
      
      !grid_offset_type used by update domains to determine shifts.
      grid_offset_type = AGRID
      if( PRESENT(gridtype) ) grid_offset_type = gridtype

      ! send
      do list = 0,nlist-1
         m = mod( domainx%pos+list, nlist )
         if( .NOT. domainx%list(m)%overlap .AND. .NOT. domainy%list(m)%overlap )cycle
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos
         select case ( grid_offset_type )
         case(BGRID_NE, BGRID_SW, AGRID)
            do dir = 1, 8  ! loop over 8 direction
               overPtrx => domainx%list(m)%send(1,1,dir)
               if(send(dir)) then
                  do n = 1, 3
                     if( overPtrx%overlap(n) ) then
                        is = overPtrx%is(n); ie = overPtrx%ie(n)
                        js = overPtrx%js(n); je = overPtrx%je(n)
                        do k = 1,ke  
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 2
                                 buffer(pos-1) = fieldx(i,j,k)
                                 buffer(pos)   = fieldy(i,j,k)
                              end do
                           end do
                        end do
                     end if
                  end do
                  nsend = overPtrx%n
                  if( nsend > 0 ) then
                     do n = 1, nsend
                        i = overPtrx%i(n); j = overPtrx%j(n)
                        do k = 1,ke  
                           pos = pos + 2
                           buffer(pos-1) = fieldx(i,j,k)
                           buffer(pos)   = fieldy(i,j,k)
                        end do
                     end do
                  end if
               end if
            enddo
         case(CGRID_NE, CGRID_SW)
            do dir = 1, 8  ! loop over 8 direction
               overPtrx => domainx%list(m)%send(1,1,dir)
               overPtry => domainy%list(m)%send(1,1,dir)
               if(send(dir)) then
                  do n = 1, 3
                     if( overPtrx%overlap(n) ) then
                        is = overPtrx%is(n); ie = overPtrx%ie(n)
                        js = overPtrx%js(n); je = overPtrx%je(n)
                        if( overPtrx%rotation == ZERO .OR.              &
                             overPtrx%rotation == ONE_HUNDRED_EIGHTY ) then
                           do k = 1,ke  
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
                     end if
                     if(overPtry%overlap(n)) then
                        is = overPtry%is(n); ie = overPtry%ie(n)
                        js = overPtry%js(n); je = overPtry%je(n)
                        if( overPtry%rotation == ZERO .OR.              &
                             overPtry%rotation == ONE_HUNDRED_EIGHTY ) then
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
                     endif
                  end do

                  nsend = overPtrx%n
                  if( nsend>0 ) then
                     do n = 1, nsend
                        i = overPtrx%i(n); j = overPtrx%j(n)
                        do k = 1,ke  
                           pos = pos + 1
                           buffer(pos) = fieldx(i,j,k)
                        end do
                     end do
                  end if
                  nsend = overPtry%n
                  if( nsend > 0 ) then
                     do n = 1, nsend
                        i = overPtry%i(n); j = overPtry%j(n)
                        do k = 1,ke  
                           pos = pos + 1
                           buffer(pos) = fieldy(i,j,k)
                        end do
                     end do
                  endif
               end if
            enddo
         end select
         call mpp_clock_end(pack_clock)

         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = domainx%list(m)%pe
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
         m = mod( domainx%pos+nlist-list, nlist )
         if( .NOT. domainx%list(m)%overlap .AND. .NOT. domainy%list(m)%overlap )cycle
         call mpp_clock_begin(recv_clock)
         msgsize = 0

         select case (grid_offset_type)
         case(BGRID_NE, BGRID_SW, AGRID)
            do dir = 1, 8  ! loop over 8 direction
               if(recv(dir)) then
                  overPtrx => domainx%list(m)%recv(1,1,dir)
                  do n = 1, 3
                     if( overPtrx%overlap(n) ) then
                        is = overPtrx%is(n); ie = overPtrx%ie(n)
                        js = overPtrx%js(n); je = overPtrx%je(n)
                        msgsize = msgsize + (ie-is+1)*(je-js+1)
                     end if
                  end do
                  msgsize = msgsize + overPtrx%n
               end if
            end do
            msgsize = msgsize*2
         case(CGRID_NE, CGRID_SW)
            do dir = 1, 8  ! loop over 8 direction
               if(recv(dir)) then
                  overPtrx => domainx%list(m)%recv(1,1,dir)
                  overPtry => domainy%list(m)%recv(1,1,dir)
                  do n = 1, 3
                     if( overPtrx%overlap(n) ) then
                        is = overPtrx%is(n); ie = overPtrx%ie(n)
                        js = overPtrx%js(n); je = overPtrx%je(n)
                        msgsize = msgsize + (ie-is+1)*(je-js+1)
                     end if
                     if( overPtry%overlap(n) ) then
                        is = overPtry%is(n); ie = overPtry%ie(n)
                        js = overPtry%js(n); je = overPtry%je(n)
                        msgsize = msgsize + (ie-is+1)*(je-js+1)
                     end if
                  end do
                  msgsize = msgsize + overPtrx%n + overPtry%n
               end if
            enddo
         end select
         msgsize = msgsize*ke

         if( msgsize.GT.0 )then
            from_pe = domainx%list(m)%pe
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
         m = mod( domainx%pos+nlist-list, nlist )
         if( .NOT.domainx%list(m)%overlap .AND. .NOT.domainy%list(m)%overlap )cycle
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         select case( grid_offset_type )
         case(BGRID_NE, BGRID_SW, AGRID)
            do dir = 8,1,-1
               if(recv(dir)) then
                  overPtrx => domainx%list(m)%recv(1,1,dir)
                  nrecv = overPtrx%n
                  if( nrecv > 0 ) then
                     msgsize = nrecv*ke*2
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = overPtrx%i(n); j = overPtrx%j(n)
                        do k = 1,ke
                           pos = pos + 2
                           fieldx(i,j,k) = buffer(pos-1)
                           fieldy(i,j,k) = buffer(pos)
                        end do
                     end do
                  endif
                  do n = 3, 1, -1
                     if( overPtrx%overlap(n) ) then
                        is = overPtrx%is(n); ie = overPtrx%ie(n)
                        js = overPtrx%js(n); je = overPtrx%je(n)
                        msgsize = (ie-is+1)*(je-js+1)*ke*2
                        pos = buffer_pos - msgsize
                        buffer_pos = pos
                        select case ( overPtrx%rotation )
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
                        case ( MINUS_NINETY )
                           if( BTEST(update_flags,SCALAR_BIT) ) then
                              do k = 1,ke
                                 do i = ie, is, -1
                                    do j = js, je
                                       pos = pos + 2
                                       fieldx(i,j,k) = buffer(pos-1)
                                       fieldy(i,j,k) = buffer(pos)
                                    end do
                                 end do
                              end do
                           else
                              do k = 1,ke
                                 do i = ie, is, -1
                                    do j = js, je
                                       pos = pos + 2
                                       fieldy(i,j,k) = buffer(pos-1)
                                       fieldx(i,j,k) = -buffer(pos)
                                    end do
                                 end do
                              end do
                           end if
                        case ( NINETY )
                           if( BTEST(update_flags,SCALAR_BIT) ) then
                              do k = 1,ke
                                 do i = is, ie
                                    do j = je, js, -1
                                       pos = pos + 2
                                       fieldx(i,j,k) = buffer(pos-1)
                                       fieldy(i,j,k) = buffer(pos)
                                    end do
                                 end do
                              end do
                           else 
                              do k = 1,ke
                                 do i = is, ie
                                    do j = je, js, -1
                                       pos = pos + 2
                                       fieldy(i,j,k) = -buffer(pos-1)
                                       fieldx(i,j,k) = buffer(pos)
                                    end do
                                 end do
                              end do
                           end if
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
                  end do
               end if
            enddo
         case(CGRID_NE, CGRID_SW)
            do dir = 8,1,-1
               if(recv(dir)) then
                  overPtrx => domainx%list(m)%recv(1,1,dir)
                  overPtry => domainy%list(m)%recv(1,1,dir)
                  nrecv = overPtry%n
                  if( nrecv > 0 ) then
                     msgsize = nrecv*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = overPtry%i(n); j = overPtry%j(n)
                        do k = 1,ke
                           pos = pos + 1
                           fieldy(i,j,k) = buffer(pos)
                        end do
                     end do
                  endif
                  nrecv = overPtrx%n
                  if( nrecv > 0) then
                     msgsize = nrecv*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = overPtrx%i(n); j = overPtrx%j(n)
                        do k = 1,ke
                           pos = pos + 1
                           fieldx(i,j,k) = buffer(pos)
                        end do
                     end do
                  end if

                  do n = 3, 1, -1
                     if( overPtry%overlap(n) ) then
                        is = overPtry%is(n); ie = overPtry%ie(n)
                        js = overPtry%js(n); je = overPtry%je(n)
                        msgsize = (ie-is+1)*(je-js+1)*ke
                        pos = buffer_pos - msgsize
                        buffer_pos = pos
                        select case ( overPtry%rotation )
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
                     if(overPtrx%overlap(n)) then
                        is = overPtrx%is(n); ie = overPtrx%ie(n)
                        js = overPtrx%js(n); je = overPtrx%je(n)
                        msgsize = (ie-is+1)*(je-js+1)*ke
                        pos = buffer_pos - msgsize
                        buffer_pos = pos
                        select case ( overPtrx%rotation )
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
                  end do
               end if
            enddo
         end select
         call mpp_clock_end(unpk_clock)
      end do

      !--- for vector field flip the sign if necessary.            
      if( BTEST(domainx%fold,NORTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then
         ioff = domainx%x(1)%memory%begin - domainx%x(1)%data%begin
         joff = domainx%y(1)%memory%begin - domainx%y(1)%data%begin
         if( BTEST(update_flags,NORTH) ) then
            js = domainx%y(1)%global%end + 1 + joff
            je = domainx%y(1)%data%end + joff
            if( je.GE.js )then
               is = domainx%x(1)%data%begin + ioff
               ie = domainx%x(1)%data%end + ioff
               !flip the sign if this is a vector
               select case(grid_offset_type)
               case (AGRID, BGRID_NE) 
                  do k = 1,ke
                     do j = js,je
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               case (CGRID_NE)
                  if(domainx%symmetry) then
                     is2 = domainy%x(1)%data%begin + ioff;     ie2 = domainy%x(1)%data%end + ioff
                     js2 = domainy%y(1)%global%end + 1 + joff; je2 = domainy%y(1)%data%end + joff      
                     do k = 1,ke
                        do j = js,je
                           do i = is,ie
                              fieldx(i,j,k) = -fieldx(i,j,k)
                           end do
                        end do
                        do j = js2,je2
                           do i = is2,ie2
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
         j = domainy%y(1)%global%end
         if( domainy%y(1)%data%begin.LE.j .AND. j.LE.domainy%y(1)%data%end )then !fold is within domain
            midpoint = (domainy%x(1)%global%begin+domainy%x(1)%global%end-1)/2
            if( domainy%x(1)%pos.GE.size(domainy%x(1)%list(:))/2 )then
               do_flip = .false.
               flip_n = BTEST(update_flags,NORTH)
               flip_e = BTEST(update_flags,EAST)
               flip_w = BTEST(update_flags,WEST)
               do_flip = flip_n .or. flip_w .or. flip_e
               do_flip2 = .FALSE.
               if ( domainy%x(1)%pos.GE.size(domainy%x(1)%list(:))/2 ) then
                  do_flip = .false.
                  flip_n = BTEST(update_flags,NORTH)
                  flip_e = BTEST(update_flags,EAST)
                  flip_w = BTEST(update_flags,WEST)
                  do_flip = flip_n .or. flip_w .or. flip_e
                  do_flip2 = .FALSE.
                  if( flip_n .AND. flip_e .AND. flip_w ) then
                     is = max(domainy%x(1)%data%begin,midpoint+1)
                     ie = domainy%x(1)%data%end
                  else if (  flip_n .AND. flip_e ) then
                     is = max(domainy%x(1)%compute%begin,midpoint+1)
                     ie = domainy%x(1)%data%end
                  else if (  flip_n .AND. flip_w ) then
                     is = max(domainy%x(1)%data%begin,midpoint+1)
                     ie = domainy%x(1)%compute%end
                  else if ( flip_n ) then
                     is = max(domainy%x(1)%compute%begin,midpoint+1)
                     ie = domainy%x(1)%compute%end                 
                  else if ( flip_w .AND. flip_e ) then  ! a little complicate
                     is = max(domainy%x(1)%compute%end+1,midpoint+1)
                     ie = domainy%x(1)%data%end
                     do_flip2 = .TRUE.
                     is2 = max(domainy%x(1)%data%begin,midpoint+1)
                     ie2 = domainy%x(1)%compute%begin-1
                     if(domainy%symmetry) ie2 = ie2 + 1
                  else if ( flip_e ) then
                     is = max(domainy%x(1)%compute%end+1,midpoint+1)
                     ie = domainy%x(1)%data%end
                  else if ( flip_w ) then
                     is = max(domainy%x(1)%data%begin,midpoint+1)
                     ie = domainy%x(1)%compute%begin-1
                     if(domainy%symmetry) ie = ie + 1
                  end if
                  if( do_flip ) then
                     is = is + ioff; ie = ie + ioff; j = j + joff
                     select case(grid_offset_type)
                     case(BGRID_NE)
                        do k = 1,ke
                           do i = is, ie
                              fieldx(i,j,k) = -fieldx(i,j,k) 
                              fieldy(i,j,k) = -fieldy(i,j,k) 
                           end do
                        end do
                        if( do_flip2 ) then
                           is2 = is2 + ioff; ie2 = ie2 + ioff
                           do k = 1,ke
                              do i = is2, ie2
                                 fieldx(i,j,k) = -fieldx(i,j,k) 
                                 fieldy(i,j,k) = -fieldy(i,j,k) 
                              end do
                           end do
                        end if
                     case(CGRID_NE)
                        do k = 1,ke
                           do i = is, ie
                              fieldy(i,j,k) = -fieldy(i,j,k) 
                           end do
                        end do
                        if( do_flip2 ) then
                           is2 = is2 + ioff; ie2 = ie2 + ioff
                           do k = 1,ke
                              do i = is2, ie2
                                 fieldy(i,j,k) = -fieldy(i,j,k) 
                              end do
                           end do
                        end if
                     end select
                  end if
               end if
            end if
            !poles set to 0: BGRID only
            if( grid_offset_type.EQ.BGRID_NE )then
               j  = domainx%y(1)%global%end + joff
               is = domainx%x(1)%global%begin; ie = domainx%x(1)%global%end
               if( .NOT. domainx%symmetry ) is = is - 1
               do i = is ,ie, midpoint
                  if( domainx%x(1)%data%begin.LE.i .AND. i.LE. domainx%x(1)%data%end )then
                     do k = 1,ke
                        fieldx(i+ioff,j,k) = 0.
                        fieldy(i+ioff,j,k) = 0.
                     end do
                  end if
               end do
            end if

            ! the last code block correct an error where the data in your halo coming from 
            ! other half may have the wrong sign
            !right of midpoint, when update north and east direction
            if ( BTEST(update_flags,NORTH) .OR. BTEST(update_flags,EAST) ) then
               j = domainy%y(1)%global%end + joff
               is = midpoint

               if( domainy%x(1)%compute%begin.LE.is .AND. is.LT.domainy%x(1)%data%end &
                   .AND. domainy%x(1)%pos < size(domainy%x(1)%list(:))/2 )then
                  is = is + 1 + ioff
                  ie = domainy%x(1)%data%end + ioff
                  select case(grid_offset_type)
                  case(BGRID_NE)
                     if(domainx%symmetry) is = is + 1
                     do k = 1,ke
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  case(CGRID_NE)
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

      !--- if debug is true and domain is symmetry, check the consistency on the bounds between tiles.
      !--- For data on AGRID, no check is needed; for data on CGRID, u will be checked on east boundary
      !--- and v will be checked on north boundary; For data on BGRID, u and v will be checked on east 
      !--- and north boundary.
      !--- The check will be done in the following way: Western/southern boundary data sent to 
      !--- Eastern boundary to check,  Southern/western boundary sent to northern boundary to check
      !--- folded-north-edge will not be checked.

      if( debug .AND. domainx%symmetry .AND.  grid_offset_type .NE. AGRID ) then      
         if(present(name)) then
            field_name = name
         else
            field_name = "un-named"
         end if

         bound_x => domainx%bound
         bound_y => domainy%bound
         ndir    = size(bound_x%list(m)%send,3)

         !--- send the data
         do list = 0,nlist-1
            m = mod( domainx%pos+list, nlist )
            pos = buffer_pos
            do dir = 1, ndir  ! ndir = 1 for E or N-cell, ndir = 2 for C-cell
               overPtrx => bound_x%list(m)%send(1,1,dir)
               overPtry => bound_y%list(m)%send(1,1,dir)
               if( overPtrx%overlap(1) ) then
                  is = overPtrx%is(1); ie = overPtrx%ie(1)
                  js = overPtrx%js(1); je = overPtrx%je(1)
                  if( overPtrx%rotation == ZERO .OR. grid_offset_type == BGRID_NE &
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
                     do k = 1, ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = fieldy(i,j,k)
                           end do
                        end do
                     end do
                  end if
               endif
               if( overPtry%overlap(1) ) then
                  is = overPtry%is(1); ie = overPtry%ie(1)
                  js = overPtry%js(1); je = overPtry%je(1)
                  if( overPtry%rotation == ZERO .OR. grid_offset_type == BGRID_NE &
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
                     do k = 1, ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = fieldx(i,j,k)
                           end do
                        end do
                     end do
                  end if
               endif
            end do
            msgsize = pos - buffer_pos
            if( msgsize.GT.0 )then
               to_pe = domainx%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, ' // &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
               end if
               call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
               buffer_pos = pos
            end if
         end do

         !--- recv the data 
         do list = 0,nlist-1
            m = mod( domainx%pos+nlist-list, nlist )
            msgsize = 0

            do dir = 1, ndir
               overPtrx => bound_x%list(m)%recv(1,1,dir)
               overPtry => bound_y%list(m)%recv(1,1,dir)
               if( overPtrx%overlap(1) ) then
                  is = overPtrx%is(1); ie = overPtrx%ie(1)
                  js = overPtrx%js(1); je = overPtrx%je(1)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end if
               if( overPtry%overlap(1) ) then
                  is = overPtry%is(1); ie = overPtry%ie(1)
                  js = overPtry%js(1); je = overPtry%je(1)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end if
            end do
            msgsize = msgsize*ke

            if( msgsize.GT.0 )then
               from_pe = domainx%list(m)%pe
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

         sign = -1
         if(BTEST(update_flags,SCALAR_BIT)) sign = 1
         !--- compare the data in reverse order
         do list = nlist-1,0,-1
            m = mod( domainx%pos+nlist-list, nlist )
            pos = buffer_pos
            do dir = ndir, 1, -1
               overPtrx => bound_x%list(m)%recv(1,1,dir)
               overPtry => bound_y%list(m)%recv(1,1,dir)
               if( overPtry%overlap(1) ) then
                  is = overPtry%is(1); ie = overPtry%ie(1)
                  js = overPtry%js(1); je = overPtry%je(1)
                  msgsize = (ie-is+1)*(je-js+1)*ke
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  select case ( overPtry%rotation )
                  case ( ZERO )
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              if( fieldy(i,j,k) .NE. buffer(pos) ) then
                                 print*,"Error from mpp_do_updateV_old.h on pe = ", mpp_pe(), ": y component of vector ", &
                                        trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldy(i,j,k), &
                                        " does not equal to the value = ", buffer(pos), " on pe ", domainx%list(m)%pe
                                 call mpp_error(FATAL, "mpp_do_updateV_old.h: mismatch on the boundary for symmetry point")
                              end if
                           end do
                        end do
                     end do
                  case ( MINUS_NINETY )
                     do k = 1,ke
                        do i = ie, is, -1
                           do j = js, je
                              pos = pos + 1
                              if( fieldy(i,j,k) .NE. buffer(pos) ) then
                                 print*,"Error from mpp_do_updateV_old.h on pe = ", mpp_pe(), ": y component of vector ", &
                                        trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldy(i,j,k), &
                                        " does not equal to the value = ", buffer(pos), " on pe ", domainx%list(m)%pe
                                 call mpp_error(FATAL, "mpp_do_updateV_old.h: mismatch on the boundary for symmetry point")
                              end if
                           end do
                        end do
                     end do
                  case ( NINETY )
                     do k = 1,ke
                        do i = is, ie
                           do j = je, js, -1
                              pos = pos + 1
                              if( fieldy(i,j,k) .NE. sign*buffer(pos) ) then
                                 print*,"Error from mpp_do_updateV_old.h on pe = ", mpp_pe(), ": y component of vector ", &
                                        trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldy(i,j,k), &
                                        " does not equal to the value = ", sign*buffer(pos), " on pe ", domainx%list(m)%pe
                                 call mpp_error(FATAL, "mpp_do_updateV_old.h: mismatch on the boundary for symmetry point")
                              end if
                           end do
                        end do
                     end do
                  case ( ONE_HUNDRED_EIGHTY )
                     do k = 1,ke
                        do j = je, js, -1
                           do i = ie, is, -1
                              pos = pos + 1
                              if( fieldy(i,j,k) .NE. sign*buffer(pos) ) then
                                 print*,"Error from mpp_do_updateV_old.h on pe = ", mpp_pe(), ": y component of vector ", &
                                        trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldy(i,j,k), &
                                        " does not equal to the value = ", sign*buffer(pos), " on pe ", domainx%list(m)%pe
                                 call mpp_error(FATAL, "mpp_do_updateV_old.h: mismatch on the boundary for symmetry point")
                              end if
                           end do
                        end do
                     end do
                  end select
               end if
               if( overPtrx%overlap(1) ) then
                  is = overPtrx%is(1); ie = overPtrx%ie(1)
                  js = overPtrx%js(1); je = overPtrx%je(1)
                  msgsize = (ie-is+1)*(je-js+1)*ke
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  select case ( overPtrx%rotation )
                  case ( ZERO )
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              if( fieldx(i,j,k) .NE. buffer(pos) ) then
                                 print*,"Error from mpp_do_updateV_old.h on pe = ", mpp_pe(), ": x-component of vector ", &
                                        trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldx(i,j,k), &
                                        " does not equal to the value = ", buffer(pos), " on pe ", domainx%list(m)%pe
                                 call mpp_error(FATAL, "mpp_do_updateV_old.h: mismatch on the boundary for symmetry point")
                              end if
                           end do
                        end do
                     end do
                  case ( MINUS_NINETY )
                     do k = 1,ke
                        do i = ie, is, -1
                           do j = js, je
                              pos = pos + 1
                              if( fieldx(i,j,k) .NE. sign*buffer(pos) ) then
                                 print*,"Error from mpp_do_updateV_old.h on pe = ", mpp_pe(), ": x-component of vector ", &
                                        trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldx(i,j,k), &
                                        " does not equal to the value = ", sign*buffer(pos), " on pe ", domainx%list(m)%pe
                                 call mpp_error(FATAL, "mpp_do_updateV_old.h: mismatch on the boundary for symmetry point")
                              end if
                           end do
                        end do
                     end do
                  case ( NINETY )
                     do k = 1,ke
                        do i = is, ie
                           do j = je, js, -1
                              pos = pos + 1
                              if( fieldx(i,j,k) .NE. buffer(pos) ) then
                                 print*,"Error from mpp_do_updateV_old.h on pe = ", mpp_pe(), ": x-component of vector ", &
                                        trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldx(i,j,k), &
                                        " does not equal to the value = ", buffer(pos), " on pe ", domainx%list(m)%pe
                                 call mpp_error(FATAL, "mpp_do_updateV_old.h: mismatch on the boundary for symmetry point")
                              end if
                           end do
                        end do
                     end do
                  case ( ONE_HUNDRED_EIGHTY )
                     do k = 1,ke
                        do j = je, js, -1
                           do i = ie, is, -1
                              pos = pos + 1
                              if( fieldx(i,j,k) .NE. sign*buffer(pos) ) then
                                 print*,"Error from mpp_do_updateV_old.h on pe = ", mpp_pe(), ": x component of vector ", &
                                        trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldx(i,j,k), &
                                        " does not equal to the value = ", sign*buffer(pos), " on pe ", domainx%list(m)%pe
                                 call mpp_error(FATAL, "mpp_do_updateV_old.h: mismatch on the boundary for symmetry point")
                              end if
                           end do
                        end do
                     end do
                  end select
               end if
            end do
         end do
         write(stdout(),*) "NOTE from mpp_do_updateV_old.h: the data on the boundary is consistent for field " &
                           //trim(field_name)
      end if

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)
 
      return

    end subroutine MPP_DOMAINS_DO_UPDATE_3Dold_V_
