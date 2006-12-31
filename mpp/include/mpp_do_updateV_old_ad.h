! -*-f90-*- 
    subroutine MPP_DO_UPDATE_AD_3Dold_V_( fieldx, fieldy, domainx, domainy, gridtype, flags, name)
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(in)    :: domainx, domainy
      MPP_TYPE_, intent(inout), dimension(domainx%x(1)%memory%begin:,domainx%y(1)%memory%begin:,:) :: fieldx, fieldy
      integer, intent(in), optional :: gridtype
      integer, intent(in), optional :: flags
      character(len=*), intent(in), optional :: name

      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif

      integer :: update_flags
      integer :: buffer_pos, nsend, nrecv
      integer :: midpoint, dir
      integer :: i, j, k, m, n, is, ie, js, je, ke, nlist
      integer :: to_pe, from_pe, list, pos, msgsize
      logical :: send(8), recv(8)
      character(len=8) :: text
      character(len=64) :: field_name      
      type(boundary),    pointer :: check_x  => NULL()
      type(boundary),    pointer :: check_y  => NULL()
      type(overlapSpec), pointer :: overPtrx => NULL()         
      type(overlapSpec), pointer :: overPtry => NULL()  

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

      if( BTEST(update_flags,NORTH) .AND. BTEST(domainx%fold,NORTH) .AND. BTEST(gridtype,SOUTH) ) &
           call mpp_error( FATAL, 'MPP_DO_UPDATE_V_OLd_ad: Incompatible grid offset and fold.' )

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
      
      ! send
      do list = 0,nlist-1
         m = mod( domainx%pos+list, nlist )
         if( .NOT. domainx%list(m)%overlap .AND. .NOT. domainy%list(m)%overlap )cycle
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos
         select case ( gridtype )
         case(BGRID_NE, BGRID_SW, AGRID)
            do dir = 1, 8  ! loop over 8 direction
               if(recv(dir)) then
                  overPtrx => domainx%list(m)%recv(1,1,dir)
                  do n = 1, 3
                     if( overPtrx%overlap(n) ) then
                        is = overPtrx%is(n); ie = overPtrx%ie(n)
                        js = overPtrx%js(n); je = overPtrx%je(n)
                           select case(overPtrx%rotation)
                           case(ZERO)
                              do k = 1,ke  
                                 do j = js, je
                                    do i = is, ie
                                       pos = pos + 2
                                       buffer(pos-1) = fieldx(i,j,k)
                                       buffer(pos)   = fieldy(i,j,k)
                                    end do
                                 end do
                              end do
                           case(MINUS_NINETY)
                              if( BTEST(update_flags,SCALAR_BIT) ) then
                                 do k = 1,ke
                                       do i = is, ie
                                    do j = je, js, -1
                                          pos = pos + 2
                                          buffer(pos-1) = fieldx(i,j,k)
                                          buffer(pos)   = fieldy(i,j,k)
                                       end do
                                    end do
                                 end do
                              else
                                 do k = 1,ke
                                    do i = is, ie
                                       do j = je, js, -1
                                          pos = pos + 2
                                          buffer(pos-1) = -fieldy(i,j,k)
                                          buffer(pos)   =  fieldx(i,j,k)
                                       end do
                                    end do
                                 end do
                              end if
                           case( NINETY )
                              if( BTEST(update_flags,SCALAR_BIT) ) then
                                 do k = 1,ke
                                    do i = ie, is, -1
                                       do j = js, je
                                          pos = pos + 2
                                       buffer(pos-1) = fieldx(i,j,k)
                                       buffer(pos)   = fieldy(i,j,k)
                                       end do
                                    end do
                                 end do
                              else
                                 do k = 1,ke
                                    do i = ie, is, -1
                                       do j = js, je
                                          pos = pos + 2
                                          buffer(pos-1) = fieldy(i,j,k)
                                          buffer(pos)   = -fieldx(i,j,k)
                                       end do
                                    end do
                                 end do
                              end if
                           case (ONE_HUNDRED_EIGHTY)  ! no refinement is considered yet
                             if( BTEST(update_flags,SCALAR_BIT) ) then
                                do k = 1,ke
                                   do j = je, js, -1
                                      do i = ie, is, -1
                                         pos = pos + 2
                                         buffer(pos-1) = fieldx(i,j,k)
                                         buffer(pos)   = fieldy(i,j,k)
                                      end do
                                   end do
                                end do
                             else
                                do k = 1,ke
                                   do j = je, js, -1
                                      do i = ie, is, -1
                                         pos = pos + 2
                                         buffer(pos-1) = -fieldx(i,j,k)
                                         buffer(pos)   = -fieldy(i,j,k)
                                      end do
                                   end do
                                end do
                             end if
                           end select
                        end if
                  end do
                  nsend = overPtrx%n
                  if( nsend > 0 ) then
                     if( BTEST(update_flags,SCALAR_BIT) ) then
                        do n = 1, nsend
                           i = overPtrx%i(n); j = overPtrx%j(n)
                           do k = 1,ke  
                              pos = pos + 2
                              buffer(pos-1) = fieldx(i,j,k)
                              buffer(pos)   = fieldy(i,j,k)
                           end do
                        end do
                     else
                        do n = 1, nsend
                           i = overPtrx%i(n); j = overPtrx%j(n)
                           do k = 1,ke  
                              pos = pos + 2
                              buffer(pos-1) = -fieldx(i,j,k)
                              buffer(pos)   = -fieldy(i,j,k)
                           end do
                        end do
                     end if
                  end if
               end if
            enddo
         case(CGRID_NE, CGRID_SW)
            do dir = 1, 8  ! loop over 8 direction
               if(send(dir)) then
                  overPtrx => domainx%list(m)%send(1,1,dir)
                  overPtry => domainy%list(m)%send(1,1,dir)
                  do n = 1, 3
                     if( overPtrx%overlap(n) ) then
                        is = overPtrx%is(n); ie = overPtrx%ie(n)
                        js = overPtrx%js(n); je = overPtrx%je(n)
                           select case( overPtrx%rotation )
                           case( ZERO )
                              do k = 1,ke  
                                 do j = js, je
                                    do i = is, ie
                                       pos = pos + 1
                                       buffer(pos) = fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           case( MINUS_NINETY ) 
                              if( BTEST(update_flags,SCALAR_BIT) ) then
                                 do k = 1,ke
                                    do i = is, ie
                                       do j = je, js, -1
                                          pos = pos + 1
                                          buffer(pos)   = fieldy(i,j,k)
                                       end do
                                    end do
                                 end do
                              else
                                 do k = 1,ke
                                    do i = is, ie
                                       do j = je, js, -1
                                          pos = pos + 1
                                          buffer(pos) = -fieldy(i,j,k)
                                       end do
                                    end do
                                 end do
                              end if
                           case( NINETY )
                              do k = 1,ke
                                 do i = ie, is, -1
                                    do j = js, je
                                       pos = pos + 1
                                       buffer(pos) = fieldy(i,j,k)
                                    end do
                                 end do
                              end do
                           case (ONE_HUNDRED_EIGHTY)  ! no refinement is considered yet
                              if( BTEST(update_flags,SCALAR_BIT) ) then
                                 do k = 1,ke
                                    do j = je, js, -1
                                       do i = ie, is, -1
                                          pos = pos + 1
                                          buffer(pos) = fieldx(i,j,k)
                                       end do
                                    end do
                                 end do
                              else
                                 do k = 1,ke
                                    do j = je, js, -1
                                       do i = ie, is, -1
                                          pos = pos + 1
                                          buffer(pos) = -fieldx(i,j,k)
                                       end do
                                    end do
                                 end do
                              end if
                           end select
                        end if
                     if(overPtry%overlap(n)) then
                        is = overPtry%is(n); ie = overPtry%ie(n)
                        js = overPtry%js(n); je = overPtry%je(n)
                           select case( overPtry%rotation )
                           case(ZERO)
                              do k = 1,ke  
                                 do j = js, je
                                    do i = is, ie
                                       pos = pos + 1
                                       buffer(pos)   = fieldy(i,j,k)
                                    end do
                                 end do
                              end do
                           case( MINUS_NINETY ) 
                              do k = 1,ke
                                 do i = is, ie
                                    do j = je, js, -1
                                       pos = pos + 1
                                       buffer(pos) = fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           case( NINETY )
                              if( BTEST(update_flags,SCALAR_BIT) ) then
                                 do k = 1,ke
                                    do i = ie, is, -1
                                       do j = js, je
                                          pos = pos + 1
                                          buffer(pos) = fieldx(i,j,k)
                                       end do
                                    end do
                                 end do
                              else
                                 do k = 1,ke
                                    do i = ie, is, -1
                                       do j = js, je
                                          pos = pos + 1
                                          buffer(pos)   = -fieldx(i,j,k)
                                       end do
                                    end do
                                 end do
                              end if
                           case (ONE_HUNDRED_EIGHTY)  ! no refinement is considered yet
                              if( BTEST(update_flags,SCALAR_BIT) ) then
                                 do k = 1,ke
                                    do j = je, js, -1
                                       do i = ie, is, -1
                                          pos = pos + 1
                                          buffer(pos)   = fieldy(i,j,k)
                                       end do
                                    end do
                                 end do
                              else
                                 do k = 1,ke
                                    do j = je, js, -1
                                       do i = ie, is, -1
                                          pos = pos + 1
                                          buffer(pos)   = -fieldy(i,j,k)
                                       end do
                                    end do
                                 end do
                              end if
                           end select
                        end if
                  end do

                  nsend = overPtrx%n
                  if( nsend>0 ) then
                     if( BTEST(update_flags,SCALAR_BIT) ) then
                        do n = 1, nsend
                           i = overPtrx%i(n); j = overPtrx%j(n)
                           do k = 1,ke  
                              pos = pos + 1
                              buffer(pos) = fieldx(i,j,k)
                           end do
                        end do
                     else
                        do n = 1, nsend
                           i = overPtrx%i(n); j = overPtrx%j(n)
                           do k = 1,ke  
                              pos = pos + 1
                              buffer(pos) = -fieldx(i,j,k)
                           end do
                        end do
                     end if
                  end if
                  nsend = overPtry%n
                  if( nsend > 0 ) then
                     if( BTEST(update_flags,SCALAR_BIT) ) then
                        do n = 1, nsend
                           i = overPtry%i(n); j = overPtry%j(n)
                           do k = 1,ke  
                              pos = pos + 1
                              buffer(pos) = fieldy(i,j,k)
                           end do
                        end do
                     else
                        do n = 1, nsend
                           i = overPtry%i(n); j = overPtry%j(n)
                           do k = 1,ke  
                              pos = pos + 1
                              buffer(pos) = -fieldy(i,j,k)
                           end do
                        end do
                     end if
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
               call mpp_error( FATAL, 'MPP_DO_UPDATE_V_OLD_AD: mpp_domains_stack overflow, ' // &
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

         select case (gridtype)
         case(BGRID_NE, BGRID_SW, AGRID)
            do dir = 1, 8  ! loop over 8 direction
               if(send(dir)) then
                  overPtrx => domainx%list(m)%send(1,1,dir)
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
               if(send(dir)) then
                  overPtrx => domainx%list(m)%send(1,1,dir)
                  overPtry => domainy%list(m)%send(1,1,dir)
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
               call mpp_error( FATAL, 'MPP_DO_UPDATE_V_OLD_AD: mpp_domains_stack overflow, '// &
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
         select case( gridtype )
         case(BGRID_NE, BGRID_SW, AGRID)
            do dir = 8,1,-1
               if(send(dir)) then
                  overPtrx => domainx%list(m)%send(1,1,dir)
                  nrecv = overPtrx%n
                  if( nrecv > 0 ) then
                     msgsize = nrecv*ke*2
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = overPtrx%i(n); j = overPtrx%j(n)
                        do k = 1,ke
                           pos = pos + 2
                           fieldx(i,j,k) = fieldx(i,j,k) + buffer(pos-1)
                           fieldy(i,j,k) = fieldy(i,j,k) + buffer(pos)
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
                           do k = 1,ke
                              do j = js, je
                                 do i = is, ie
                                    pos = pos + 2
                                 fieldx(i,j,k) = fieldx(i,j,k) + buffer(pos-1)
                                 fieldy(i,j,k) = fieldy(i,j,k) + buffer(pos)
                                 end do
                              end do
                           end do
                        end if
                  end do
               end if
            enddo
         case(CGRID_NE, CGRID_SW)
            do dir = 8,1,-1
               if(send(dir)) then
                  overPtrx => domainx%list(m)%send(1,1,dir)
                  overPtry => domainy%list(m)%send(1,1,dir)
                  nrecv = overPtry%n
                  if( nrecv > 0 ) then
                     msgsize = nrecv*ke
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do n = 1, nrecv
                        i = overPtry%i(n); j = overPtry%j(n)
                        do k = 1,ke
                           pos = pos + 1
                           fieldy(i,j,k) = fieldy(i,j,k) + buffer(pos)
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
                           fieldx(i,j,k) = fieldx(i,j,k) + buffer(pos)
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
                           do k = 1,ke
                              do j = js, je
                                 do i = is, ie
                                    pos = pos + 1
                                 fieldy(i,j,k) = fieldy(i,j,k) + buffer(pos)
                                 end do
                              end do
                           end do
                        end if
                     if(overPtrx%overlap(n)) then
                        is = overPtrx%is(n); ie = overPtrx%ie(n)
                        js = overPtrx%js(n); je = overPtrx%je(n)
                        msgsize = (ie-is+1)*(je-js+1)*ke
                        pos = buffer_pos - msgsize
                        buffer_pos = pos
                           do k = 1,ke
                              do j = js, je
                                 do i = is, ie
                                    pos = pos + 1
                                 fieldx(i,j,k) = fieldx(i,j,k) + buffer(pos)
                                 end do
                              end do
                           end do
                        end if
                  end do
               end if
            enddo
         end select
         call mpp_clock_end(unpk_clock)
      end do

      !--- for vector field flip the sign if necessary.            
      if( BTEST(domainx%fold,NORTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then
         !--- flip the sign at the folded north edge.
         j = domainy%y(1)%global%end
         if( domainy%y(1)%data%begin.LE.j .AND. j.LE.domainy%y(1)%data%end )then !fold is within domain
            midpoint = (domainy%x(1)%global%begin+domainy%x(1)%global%end-1)/2
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               j  = domainx%y(1)%global%end 
               is = domainx%x(1)%global%begin; ie = domainx%x(1)%global%end
               if( .NOT. domainx%symmetry ) is = is - 1
               do i = is ,ie, midpoint
                  if( domainx%x(1)%data%begin.LE.i .AND. i.LE. domainx%x(1)%data%end )then
                     do k = 1,ke
                        fieldx(i,j,k) = 0.
                        fieldy(i,j,k) = 0.
                     end do
                  end if
               end do
            end if

            ! the following code block correct an error where the data in your halo coming from 
            ! other half may have the wrong sign

            !off west edge, when update north or west direction
            j = domainy%y(1)%global%end 
            if ( BTEST(update_flags,NORTH) .OR. BTEST(update_flags,WEST) ) then
               select case(gridtype)
               case(BGRID_NE)
                  is = domainx%x(1)%global%begin - 1
                  if(domainx%symmetry) then
                     is = domainx%x(1)%global%begin
                  else
                     is = domainx%x(1)%global%begin - 1
                  end if
                  if( is.GT.domainx%x(1)%data%begin )then

                     if( 2*is-domainx%x(1)%data%begin.GT.domainx%x(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V_OLD_AD: BGRID_NE west edge ubound error.' )
                     do k = 1,ke
                        do i = domainx%x(1)%data%begin,is-1
                           fieldx(i,j,k) = fieldx(2*is-i,j,k)
                           fieldy(i,j,k) = fieldy(2*is-i,j,k)
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  is = domainy%x(1)%global%begin
                  if( is.GT.domainy%x(1)%data%begin )then
                     if( 2*is-domainy%x(1)%data%begin-1.GT.domainy%x(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V_OLD_AD: CGRID_NE west edge ubound error.' )
                     do k = 1,ke
                        do i = domainy%x(1)%data%begin,is-1
                           fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                        end do
                     end do
                  end if
               end select
            end if

            !off east edge
            is = domainy%x(1)%global%end
            if(domainy%x(1)%cyclic .AND. is.LT.domainy%x(1)%data%end )then
               ie = domainy%x(1)%data%end
               is = is + 1
               select case(gridtype)
               case(BGRID_NE)
                  do k = 1,ke
                     do i = is,ie
                        fieldx(i,j,k) = -fieldx(i,j,k)
                        fieldy(i,j,k) = -fieldy(i,j,k)
                     end do
                  end do
               case(CGRID_NE)
                  do k = 1,ke
                     do i = is, ie
                        fieldy(i,j,k) = -fieldy(i,j,k)
                     end do
                  end do
               end select
            end if
         end if
      end if

      !--- if debug_update_domain is true and domain is symmetry, check the consistency on the bounds between tiles.
      !--- For data on AGRID, no check is needed; for data on CGRID, u will be checked on east boundary
      !--- and v will be checked on north boundary; For data on BGRID, u and v will be checked on east 
      !--- and north boundary.
      !--- The check will be done in the following way: Western/southern boundary data sent to 
      !--- Eastern boundary to check,  Southern/western boundary sent to northern boundary to check
      !--- folded-north-edge will not be checked.

      if( debug_update_domain .AND. domainx%symmetry .AND.  gridtype .NE. AGRID ) then      
         if(present(name)) then
            field_name = name
         else
            field_name = "un-named"
         end if

         !--- send the data
         do list = 0,nlist-1
            m = mod( domainx%pos+list, nlist )
            check_x => domainx%check%send(m)
            check_y => domainy%check%send(m)
            pos = buffer_pos
            do n = 1, check_x%count
               is = check_x%is(n); ie = check_x%ie(n)
               js = check_x%js(n); je = check_x%je(n)
               select case( check_x%rotation(n) )
               case(ZERO)
                  do k = 1,ke  
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           buffer(pos) = fieldx(i,j,k)
                        end do
                     end do
                  end do
               case(MINUS_NINETY)
                  if( BTEST(update_flags,SCALAR_BIT) ) then
                     do k = 1, ke
                        do j = je, js, -1
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = fieldy(i,j,k)
                           end do
                        end do
                     end do
                  else
                     do k = 1, ke
                        do j = je, js, -1
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = -fieldy(i,j,k)
                           end do
                        end do
                     end do
                  end if
               case(NINETY)
                  do k = 1, ke
                     do j = js, je
                        do i = ie, is, -1
                           pos = pos + 1
                           buffer(pos) = fieldy(i,j,k)
                        end do
                     end do
                  end do
               case(ONE_HUNDRED_EIGHTY) 
                  if( BTEST(update_flags,SCALAR_BIT) ) then
                     do k = 1, ke
                        do j = je, js, -1
                           do i = ie, is, -1
                              pos = pos + 1
                              buffer(pos) = fieldx(i,j,k)
                           end do
                        end do
                     end do
                  else
                     do k = 1, ke
                        do j = je, js, -1
                           do i = ie, is, -1
                              pos = pos + 1
                              buffer(pos) = -fieldx(i,j,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end do
            do n = 1, check_y%count
               is = check_y%is(n); ie = check_y%ie(n)
               js = check_y%js(n); je = check_y%je(n)
               select case( check_y%rotation(n) )
               case(ZERO)
                  do k = 1,ke  
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           buffer(pos) = fieldy(i,j,k)
                        end do
                     end do
                  end do
               case(MINUS_NINETY)
                  do k = 1, ke
                     do j = je, js, -1
                        do i = is, ie
                           pos = pos + 1
                           buffer(pos) = fieldx(i,j,k)
                        end do
                     end do
                  end do
               case(NINETY)
                  if( BTEST(update_flags,SCALAR_BIT) ) then
                     do k = 1, ke
                        do j = js, je
                           do i = ie, is, -1
                              pos = pos + 1
                              buffer(pos) = fieldx(i,j,k)
                           end do
                        end do
                     end do
                  else
                     do k = 1, ke
                        do j = js, je
                           do i = ie, is, -1
                              pos = pos + 1
                              buffer(pos) = -fieldx(i,j,k)
                           end do
                        end do
                     end do
                  end if
               case(ONE_HUNDRED_EIGHTY) 
                  if( BTEST(update_flags,SCALAR_BIT) ) then
                     do k = 1, ke
                        do j = je, js, -1
                           do i = ie, is, -1
                              pos = pos + 1
                              buffer(pos) = fieldy(i,j,k)
                           end do
                        end do
                     end do
                  else
                     do k = 1, ke
                        do j = je, js, -1
                           do i = ie, is, -1
                              pos = pos + 1
                              buffer(pos) = -fieldy(i,j,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end do
            msgsize = pos - buffer_pos
            if( msgsize.GT.0 )then
               to_pe = domainx%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_DO_UPDATE_V_OLD_AD: mpp_domains_stack overflow, ' // &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
               end if
               call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
               buffer_pos = pos
            end if
         end do

         !--- recv the data 
         do list = 0,nlist-1
            m = mod( domainx%pos+nlist-list, nlist )
            check_x=>domainx%check%recv(m)
            check_y=>domainy%check%recv(m)
            msgsize = 0
            do n = 1, check_x%count
               is = check_x%is(n); ie = check_x%ie(n)
               js = check_x%js(n); je = check_x%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            do n = 1, check_y%count
               is = check_y%is(n); ie = check_y%ie(n)
               js = check_y%js(n); je = check_y%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            msgsize = msgsize*ke

            if( msgsize.GT.0 )then
               from_pe = domainx%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_DO_UPDATE_V_OLD_AD: mpp_domains_stack overflow, '// &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
               end if
               call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
               buffer_pos = buffer_pos + msgsize
            end if
         end do

         !--- compare the data in reverse order
         do list = nlist-1,0,-1
            m = mod( domainx%pos+nlist-list, nlist )
            check_x=>domainx%check%recv(m)
            check_y=>domainy%check%recv(m)
            do n = check_y%count, 1, -1
               is = check_y%is(n); ie = check_y%ie(n)
               js = check_y%js(n); je = check_y%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke
               pos = buffer_pos - msgsize
               buffer_pos = pos
               do k = 1,ke
                  do j = js, je
                     do i = is, ie
                        pos = pos + 1
                        if( fieldy(i,j,k) .NE. buffer(pos) ) then
                              print*,"Error from MPP_DO_UPDATE_V_OLD_AD on pe = ", mpp_pe(), ": y component of vector ", &
                                trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldy(i,j,k), &
                                " does not equal to the value = ", buffer(pos), " on pe ", domainx%list(m)%pe
                              call mpp_error(FATAL, "MPP_DO_UPDATE_V_OLD_AD: mismatch on the boundary for symmetry point")
                        end if
                     end do
                  end do
               end do
            end do
            do n = check_x%count, 1, -1
               is = check_x%is(n); ie = check_x%ie(n)
               js = check_x%js(n); je = check_x%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke
               pos = buffer_pos - msgsize
               buffer_pos = pos
               do k = 1,ke
                  do j = js, je
                     do i = is, ie
                        pos = pos + 1
                        if( fieldx(i,j,k) .NE. buffer(pos) ) then
                              print*,"Error from MPP_DO_UPDATE_V_OLD_AD on pe = ", mpp_pe(), ": x-component of vector ", &
                                trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldx(i,j,k), &
                                " does not equal to the value = ", buffer(pos), " on pe ", domainx%list(m)%pe
                              call mpp_error(FATAL, "MPP_DO_UPDATE_V_OLD_AD: mismatch on the boundary for symmetry point")
                        end if
                     end do
                  end do
               end do
            end do
         end do
         write(stdout(),*) "NOTE from MPP_DO_UPDATE_V_OLD_AD: the data on the boundary is consistent for field " &
              //trim(field_name)
      end if

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)

      return

    end subroutine MPP_DO_UPDATE_AD_3Dold_V_
