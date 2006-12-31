! -*-f90-*- 
subroutine MPP_DO_GET_BOUNDARY_3Dnew_( f_addrs, domain, b_addrs, bsize, ke, d_type, flags)
  type(domain2D), intent(in)      :: domain
  integer(LONG_KIND), intent(in)  :: f_addrs(:,:)
  integer(LONG_KIND), intent(in)  :: b_addrs(:,:,:)
  integer,            intent(in)  :: bsize(:), ke
  MPP_TYPE_, intent(in)           :: d_type  ! creates unique interface
  integer, intent(in)             :: flags
  
  MPP_TYPE_ :: field(domain%x(1)%memory%begin:domain%x(1)%memory%end, domain%y(1)%memory%begin:domain%y(1)%memory%end,ke)
  MPP_TYPE_ :: ebuffer(bsize(1), ke), sbuffer(bsize(2), ke), wbuffer(bsize(3), ke), nbuffer(bsize(4), ke)
  pointer(ptr_field, field)
  pointer(ptr_ebuffer, ebuffer)
  pointer(ptr_sbuffer, sbuffer)  
  pointer(ptr_wbuffer, wbuffer)
  pointer(ptr_nbuffer, nbuffer)

  logical                 :: recv(4), send(4)
  integer                 :: nlist, buffer_pos, list, pos, tMe
  integer                 :: i, j, k, l, m, n, index
  integer                 :: is, ie, js, je, msgsize, l_size
  character(len=8)        :: text
  type(boundary), pointer :: bound => NULL()

  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))

  pointer( ptr, buffer )
  ptr = LOC(mpp_domains_stack)

  l_size = size(f_addrs,1)
  recv(1) = BTEST(flags,EAST)
  recv(2) = BTEST(flags,SOUTH)
  recv(3) = BTEST(flags,WEST)
  recv(4) = BTEST(flags,NORTH)

  send = recv

  nlist = size(domain%list(:))  
  buffer_pos = 0     

  ! send
  do list = 0,nlist-1
     m = mod( domain%pos+list, nlist )
     bound=>domain%bound%send(m)
     if( bound%count == 0) cycle
     pos = buffer_pos

     do n = 1, bound%count
        if(send(bound%dir(n))) then
           is = bound%is(n); ie = bound%ie(n)
           js = bound%js(n); je = bound%je(n)
           tMe = bound%tileMe(n)
           select case( bound%rotation(n) )
           case(ZERO)
              do l=1,l_size
                 ptr_field = f_addrs(l, tMe)
                 do k = 1, ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case( MINUS_NINETY )
              do l=1,l_size
                 ptr_field = f_addrs(l, tMe)
                 do k = 1, ke
                    do j = je, js, -1
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case( NINETY )
              do l=1,l_size
                 ptr_field = f_addrs(l, tMe)
                 do k = 1, ke
                    do j = js, je
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case (ONE_HUNDRED_EIGHTY) 
              do l=1,l_size
                 ptr_field = f_addrs(l, tMe)
                 do k = 1, ke
                    do j = je, js, -1
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           end select
        end if ! if(send(bound%dir(n)))
     end do ! do n = 1, bound%count
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then  
        !--- maybe we do not need the following stack size check.
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_OLD: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=domain%list(m)%pe )
        buffer_pos = pos
     end if
  end do

  !recv
  do list = 0,nlist-1
     m = mod( domain%pos+nlist-list, nlist )
     bound=>domain%bound%recv(m)
     if( bound%count == 0) cycle
     msgsize = 0
     do n = 1, bound%count
        if(recv(bound%dir(n))) then
           is = bound%is(n); ie = bound%ie(n)
           js = bound%js(n); je = bound%je(n)
           msgsize = msgsize + (ie-is+1)*(je-js+1)
        end if
     end do
     msgsize = msgsize*ke*l_size
     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_OLD: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=domain%list(m)%pe )
        buffer_pos = buffer_pos + msgsize
     end if
  end do

  !unpack recv
  !unpack buffer in reverse order.
  do list = nlist-1,0,-1
     m = mod( domain%pos+nlist-list, nlist )
     bound=>domain%bound%recv(m)
     if( bound%count == 0) cycle
     do n = bound%count, 1, -1
        if(recv(bound%dir(n))) then
           is = bound%is(n); ie = bound%ie(n)
           js = bound%js(n); je = bound%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           tMe = bound%tileMe(n)
           select case( bound%dir(n) )
           case ( 1 ) ! EAST
              do l=1,l_size
                 ptr_ebuffer = b_addrs(1, l, tMe)              
                 do k = 1, ke
                    index = bound%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          ebuffer(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 2 ) ! SOUTH
              do l=1,l_size
                 ptr_sbuffer = b_addrs(2, l, tMe)   
                 do k = 1, ke
                    index = bound%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          sbuffer(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 3 ) ! WEST
              do l=1,l_size
                 ptr_wbuffer = b_addrs(3, l, tMe)   
                 do k = 1, ke
                    index = bound%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          wbuffer(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 4 ) ! norTH
              do l=1,l_size
                 ptr_nbuffer = b_addrs(4, l, tMe)   
                 do k = 1, ke
                    index = bound%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          nbuffer(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           end select
        end if
     end do
  end do

  call mpp_sync_self( )


end subroutine MPP_DO_GET_BOUNDARY_3Dnew_


subroutine MPP_DO_GET_BOUNDARY_3Dnew_V_(f_addrsx, f_addrsy, domainx, domainy, b_addrsx, b_addrsy, &
                                        bsizex, bsizey, ke, d_type, flags)
  type(domain2D), intent(in)      :: domainx, domainy
  integer(LONG_KIND), intent(in)  :: f_addrsx(:,:), f_addrsy(:,:)
  integer(LONG_KIND), intent(in)  :: b_addrsx(:,:,:), b_addrsy(:,:,:)
  integer,            intent(in)  :: bsizex(:), bsizey(:), ke
  MPP_TYPE_, intent(in)           :: d_type  ! creates unique interface
  integer, intent(in)             :: flags

  MPP_TYPE_ :: fieldx(domainx%x(1)%memory%begin:domainx%x(1)%memory%end, domainx%y(1)%memory%begin:domainx%y(1)%memory%end,ke)
  MPP_TYPE_ :: fieldy(domainy%x(1)%memory%begin:domainy%x(1)%memory%end, domainy%y(1)%memory%begin:domainy%y(1)%memory%end,ke)
  MPP_TYPE_ :: ebufferx(bsizex(1), ke), sbufferx(bsizex(2), ke), wbufferx(bsizex(3), ke), nbufferx(bsizex(4), ke)
  MPP_TYPE_ :: ebuffery(bsizey(1), ke), sbuffery(bsizey(2), ke), wbuffery(bsizey(3), ke), nbuffery(bsizey(4), ke)
  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)
  pointer(ptr_ebufferx, ebufferx)
  pointer(ptr_sbufferx, sbufferx)  
  pointer(ptr_wbufferx, wbufferx)
  pointer(ptr_nbufferx, nbufferx)
  pointer(ptr_ebuffery, ebuffery)
  pointer(ptr_sbuffery, sbuffery)  
  pointer(ptr_wbuffery, wbuffery)
  pointer(ptr_nbuffery, nbuffery)

  logical                 :: recv(4), send(4)
  integer                 :: nlist, buffer_pos, list, pos, tMe
  integer                 :: is, ie, js, je, msgsize, l_size
  integer                 :: i, j, k, l, m, n, index
  character(len=8)        :: text
  type(boundary), pointer :: boundx => NULL()
  type(boundary), pointer :: boundy => NULL()

  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
#ifdef use_CRI_pointers
  pointer( ptr, buffer )
  ptr = LOC(mpp_domains_stack)
#endif

  l_size = size(f_addrsx,1)
  recv(1) = BTEST(flags,EAST)
  recv(2) = BTEST(flags,SOUTH)
  recv(3) = BTEST(flags,WEST)
  recv(4) = BTEST(flags,NORTH)
  send = recv

  nlist = size(domainx%list(:))  
  buffer_pos = 0     

  ! send
  do list = 0,nlist-1
     m = mod( domainx%pos+list, nlist )
     boundx=>domainx%bound%send(m)
     boundy=>domainy%bound%send(m)
     pos = buffer_pos
     do n = 1, boundx%count   
        if(send(boundx%dir(n))) then
           is = boundx%is(n); ie = boundx%ie(n)
           js = boundx%js(n); je = boundx%je(n)
           tMe = boundx%tileMe(n)
           select case( boundx%rotation(n) )
           case(ZERO)
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, tMe)
                 do k = 1, ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = fieldx(i,j,k)
                       end do
                    end do
                 end do
              end do
           case( MINUS_NINETY )
              if( BTEST(flags,SCALAR_BIT) ) then
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, tMe)
                    do k = 1, ke
                       do j = je, js, -1
                          do i = is, ie
                             pos = pos + 1
                             buffer(pos) = fieldy(i,j,k)
                          end do
                       end do
                    end do
                 end do
              else
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, tMe)
                    do k = 1, ke
                       do j = je, js, -1
                          do i = is, ie
                             pos = pos + 1
                             buffer(pos) = -fieldy(i,j,k)
                          end do
                       end do
                    end do
                 end do
              end if
           case( NINETY )
              do l=1,l_size
                 ptr_fieldy = f_addrsy(l, tMe)
                 do k = 1, ke
                    do j = js, je
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) = fieldy(i,j,k)
                       end do
                    end do
                 end do
              end do
           case (ONE_HUNDRED_EIGHTY) 
              if( BTEST(flags,SCALAR_BIT) ) then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, tMe)
                    do k = 1, ke
                       do j = je, js, -1
                          do i = ie, is, -1
                             pos = pos + 1
                             buffer(pos) = fieldx(i,j,k)
                          end do
                       end do
                    end do
                 end do
              else     
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, tMe)
                    do k = 1, ke
                       do j = je, js, -1
                          do i = ie, is, -1
                             pos = pos + 1
                             buffer(pos) = -fieldx(i,j,k)
                          end do
                       end do
                    end do
                 end do
              end if
           end select
        end if ! if(send(boundx%dir(n)))
     end do  !do n = 1, boundx%count

     do n = 1, boundy%count   
        if(send(boundy%dir(n))) then
           is = boundy%is(n); ie = boundy%ie(n)
           js = boundy%js(n); je = boundy%je(n)
           tMe = boundy%tileMe(n)
           select case( boundy%rotation(n) )
           case(ZERO)
              do l=1,l_size
                 ptr_fieldy = f_addrsy(l, tMe)
                 do k = 1, ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = fieldy(i,j,k)
                       end do
                    end do
                 end do
              end do
           case( MINUS_NINETY )
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, tMe)
                 do k = 1, ke
                    do j = je, js, -1
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = fieldx(i,j,k)
                       end do
                    end do
                 end do
              end do
           case( NINETY )
              if( BTEST(flags,SCALAR_BIT) ) then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, tMe)
                    do k = 1, ke
                       do j = js, je
                          do i = ie, is, -1
                             pos = pos + 1
                             buffer(pos) = fieldx(i,j,k)
                          end do
                       end do
                    end do
                 end do
              else
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, tMe)
                    do k = 1, ke
                       do j = js, je
                          do i = ie, is, -1
                             pos = pos + 1
                             buffer(pos) = -fieldx(i,j,k)
                          end do
                       end do
                    end do
                 end do
              end if
           case (ONE_HUNDRED_EIGHTY) 
              if( BTEST(flags,SCALAR_BIT) ) then
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, tMe)
                    do k = 1, ke
                       do j = je, js, -1
                          do i = ie, is, -1
                             pos = pos + 1
                             buffer(pos) = fieldy(i,j,k)
                          end do
                       end do
                    end do
                 end do
              else     
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, tMe)
                    do k = 1, ke
                       do j = je, js, -1
                          do i = ie, is, -1
                             pos = pos + 1
                             buffer(pos) = -fieldy(i,j,k)
                          end do
                       end do
                    end do
                 end do
              end if
           end select
        end if ! if(send(boundy%dir(n)))
     end do    ! do n = 1, boundy%count
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then  
        !--- maybe we do not need the following stack size check.
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_V_NEW: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=domainx%list(m)%pe )
        buffer_pos = pos
     end if

  end do       !  do list = 0,nlist-1

  !recv
  do list = 0,nlist-1
     m = mod( domainx%pos+nlist-list, nlist )
     boundx=>domainx%bound%recv(m)
     boundy=>domainy%bound%recv(m)
     msgsize = 0
     do n = 1, boundx%count
        if(recv(boundx%dir(n))) then
           is = boundx%is(n); ie = boundx%ie(n)
           js = boundx%js(n); je = boundx%je(n)
           msgsize = msgsize + (ie-is+1)*(je-js+1)
        end if
     end do
     do n = 1, boundy%count
        if(recv(boundy%dir(n))) then
           is = boundy%is(n); ie = boundy%ie(n)
           js = boundy%js(n); je = boundy%je(n)
           msgsize = msgsize + (ie-is+1)*(je-js+1)
        end if
     end do
     msgsize = msgsize*ke*l_size
     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_V_NEW: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=domainx%list(m)%pe )
        buffer_pos = buffer_pos + msgsize
     end if
  end do

  !unpack recv
  !unpack buffer in reverse order.
  do list = nlist-1,0,-1
     m = mod( domainx%pos+nlist-list, nlist )
     boundx=>domainx%bound%recv(m)
     boundy=>domainy%bound%recv(m)
     do n = boundy%count, 1, -1
        if(recv(boundy%dir(n))) then
           is = boundy%is(n); ie = boundy%ie(n)
           js = boundy%js(n); je = boundy%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           tMe = boundy%tileMe(n)
           select case( boundy%dir(n) )
           case ( 1 ) ! EAST
              do l=1,l_size
                 ptr_ebuffery = b_addrsy(1, l, tMe)     
                 do k = 1, ke
                    index = boundy%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          ebuffery(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 2 ) ! SOUTH
              do l=1,l_size
                 ptr_sbuffery = b_addrsy(2, l, tMe)     
                 do k = 1, ke
                    index = boundy%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          sbuffery(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 3 ) ! WEST
              do l=1,l_size
                 ptr_wbuffery = b_addrsy(3, l, tMe)     
                 do k = 1, ke
                    index = boundy%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          wbuffery(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 4 ) ! norTH
              do l=1,l_size
                 ptr_nbuffery = b_addrsy(4, l, tMe)     
                 do k = 1, ke
                    index = boundy%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          nbuffery(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           end select
        end if
     end do

     do n = boundx%count, 1, -1
        if(recv(boundx%dir(n))) then
           is = boundx%is(n); ie = boundx%ie(n)
           js = boundx%js(n); je = boundx%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           tMe = boundx%tileMe(n)
           select case( boundx%dir(n) )
           case ( 1 ) ! EAST
              do l=1,l_size
                 ptr_ebufferx = b_addrsx(1, l, tMe)     
                 do k = 1, ke
                    index = boundx%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          ebufferx(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 2 ) ! SOUTH
              do l=1,l_size
                 ptr_sbufferx = b_addrsx(2, l, tMe)     
                 do k = 1, ke
                    index = boundx%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          sbufferx(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 3 ) ! WEST
              do l=1,l_size
                 ptr_wbufferx = b_addrsx(3, l, tMe)     
                 do k = 1, ke
                    index = boundx%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          wbufferx(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 4 ) ! norTH
              do l=1,l_size
                 ptr_nbufferx = b_addrsx(4, l, tMe)     
                 do k = 1, ke
                    index = boundx%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          nbufferx(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           end select
        end if
     end do

  end do

  call mpp_sync_self( )



end subroutine MPP_DO_GET_BOUNDARY_3Dnew_V_
