! -*-f90-*- 
subroutine MPP_DO_GET_BOUNDARY_3Dold_( field, domain, flags, ebuffer, sbuffer, wbuffer, nbuffer)
  type(domain2D), intent(in)           :: domain
  integer, intent(in)                  :: flags
  MPP_TYPE_, intent(in)                :: field(domain%x(1)%memory%begin:,domain%y(1)%memory%begin:,:)
  MPP_TYPE_, intent(inout), optional   :: ebuffer(:,:), sbuffer(:,:), wbuffer(:,:), nbuffer(:,:)


  type(boundary), pointer :: bound => NULL()
  integer                 :: buffer_pos, pos, msgsize, index
  integer                 :: nlist, ke, list, m, n, k, i, j
  integer                 :: is, ie, js, je
  logical                 :: recv(4), send(4)
  character(len=8)        :: text

  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))

#ifdef use_CRI_pointers
  pointer( ptr, buffer )
  ptr = LOC(mpp_domains_stack)
#endif

  recv(1) = BTEST(flags,EAST)
  recv(2) = BTEST(flags,SOUTH)
  recv(3) = BTEST(flags,WEST)
  recv(4) = BTEST(flags,NORTH)

  send = recv

  nlist = size(domain%list(:))  
  ke = size(field,3)
  buffer_pos = 0     

  ! send
  do list = 0,nlist-1
     m = mod( domain%pos+list, nlist )
     bound=>domain%bound%send(m)
     pos = buffer_pos
     do n = 1, bound%count
        if(send(bound%dir(n))) then
           is = bound%is(n); ie = bound%ie(n)
           js = bound%js(n); je = bound%je(n)
           select case( bound%rotation(n) )
           case(ZERO)
              do k = 1, ke
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              end do
           case( MINUS_NINETY )
              do k = 1, ke
                 do j = je, js, -1
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              end do 
           case( NINETY )
              do k = 1, ke
                 do j = js, je
                    do i = ie, is, -1
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              end do 
           case (ONE_HUNDRED_EIGHTY) 
              do k = 1, ke
                 do j = je, js, -1
                    do i = ie, is, -1
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
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
     msgsize = 0
     do n = 1, bound%count
        if(recv(bound%dir(n))) then
           is = bound%is(n); ie = bound%ie(n)
           js = bound%js(n); je = bound%je(n)
           msgsize = msgsize + (ie-is+1)*(je-js+1)
        end if
     end do
     msgsize = msgsize*ke
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
     do n = bound%count, 1, -1
        if(recv(bound%dir(n))) then
           is = bound%is(n); ie = bound%ie(n)
           js = bound%js(n); je = bound%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke
           pos = buffer_pos - msgsize
           buffer_pos = pos
           select case( bound%dir(n) )
           case ( 1 ) ! EAST
              if(.NOT. present(ebuffer)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD: optional argument ebuffer should be presented")
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
           case ( 2 ) ! SOUTH
              if(.NOT. present(sbuffer)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD: optional argument sbuffer should be presented")
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
           case ( 3 ) ! WEST
              if(.NOT. present(wbuffer)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD: optional argument wbuffer should be presented")
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
           case ( 4 ) ! norTH
              if(.NOT. present(nbuffer)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD: optional argument nbuffer should be presented")
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
           end select
        end if
     end do
  end do

  call mpp_sync_self( )

end subroutine MPP_DO_GET_BOUNDARY_3Dold_

!##########################################################################
subroutine MPP_DO_GET_BOUNDARY_3Dold_V_(fieldx, fieldy, domainx, domainy, flags, ebufferx, sbufferx, wbufferx, nbufferx, &
                                        ebuffery, sbuffery, wbuffery, nbuffery)
  type(domain2D), intent(in)         :: domainx, domainy
  MPP_TYPE_, intent(in)              :: fieldx(domainx%x(1)%memory%begin:,domainx%y(1)%memory%begin:,:)
  MPP_TYPE_, intent(in)              :: fieldy(domainy%x(1)%memory%begin:,domainy%y(1)%memory%begin:,:)
  integer, intent(in)                :: flags
  MPP_TYPE_, intent(inout), optional :: ebufferx(:,:), sbufferx(:,:), wbufferx(:,:), nbufferx(:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffery(:,:), sbuffery(:,:), wbuffery(:,:), nbuffery(:,:)


  type(boundary), pointer :: boundx => NULL()
  type(boundary), pointer :: boundy => NULL()
  logical                 :: recv(4), send(4)
  integer                 :: nlist, ke, is, ie, js, je, msgsize, index
  integer                 :: i, j, k, m, n, list, pos, buffer_pos
  character(len=8)        :: text

  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
#ifdef use_CRI_pointers
  pointer( ptr, buffer )
  ptr = LOC(mpp_domains_stack)
#endif

  recv(1) = BTEST(flags,EAST)
  recv(2) = BTEST(flags,SOUTH)
  recv(3) = BTEST(flags,WEST)
  recv(4) = BTEST(flags,NORTH)
  send = recv

  nlist = size(domainx%list(:))  
  ke = size(fieldx,3)
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
           select case( boundx%rotation(n) )
           case(ZERO)
              do k = 1, ke
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = fieldx(i,j,k)
                    end do
                 end do
              end do
           case( MINUS_NINETY )
              if( BTEST(flags,SCALAR_BIT) ) then
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
           case( NINETY )
              do k = 1, ke
                 do j = js, je
                    do i = ie, is, -1
                       pos = pos + 1
                       buffer(pos) = fieldy(i,j,k)
                    end do
                 end do
              end do 
           case (ONE_HUNDRED_EIGHTY) 
              if( BTEST(flags,SCALAR_BIT) ) then
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
        end if ! if(send(boundx%dir(n)))
     end do  !do n = 1, boundx%count

     do n = 1, boundy%count   
        if(send(boundy%dir(n))) then
           is = boundy%is(n); ie = boundy%ie(n)
           js = boundy%js(n); je = boundy%je(n)
           select case( boundy%rotation(n) )
           case(ZERO)
              do k = 1, ke
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = fieldy(i,j,k)
                    end do
                 end do
              end do
           case( MINUS_NINETY )
              do k = 1, ke
                 do j = je, js, -1
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = fieldx(i,j,k)
                    end do
                 end do
              end do
           case( NINETY )
              if( BTEST(flags,SCALAR_BIT) ) then
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
           case (ONE_HUNDRED_EIGHTY) 
              if( BTEST(flags,SCALAR_BIT) ) then
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
        end if ! if(send(boundy%dir(n)))
     end do    ! do n = 1, boundy%count
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then  
        !--- maybe we do not need the following stack size check.
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_V_OLD: mpp_domains_stack overflow, ' // &
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
     msgsize = msgsize*ke
     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_V_OLD: mpp_domains_stack overflow, '// &
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
           msgsize = (ie-is+1)*(je-js+1)*ke
           pos = buffer_pos - msgsize
           buffer_pos = pos
           select case( boundy%dir(n) )
           case ( 1 ) ! EAST
              if(.NOT. present(ebuffery)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD_V: optional argument ebuffery should be presented")
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
           case ( 2 ) ! SOUTH
              if(.NOT. present(sbuffery)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD_V: optional argument sbuffery should be presented")
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
           case ( 3 ) ! WEST
              if(.NOT. present(wbuffery)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD_V: optional argument wbuffery should be presented")
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
           case ( 4 ) ! norTH
              if(.NOT. present(nbuffery)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD_V: optional argument nbuffery should be presented")
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
           end select
        end if
     end do

     do n = boundx%count, 1, -1
        if(recv(boundx%dir(n))) then
           is = boundx%is(n); ie = boundx%ie(n)
           js = boundx%js(n); je = boundx%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke
           pos = buffer_pos - msgsize
           buffer_pos = pos
           select case( boundx%dir(n) )
           case ( 1 ) ! EAST
              if(.NOT. present(ebufferx)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD_V: optional argument ebufferx should be presented")
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
           case ( 2 ) ! SOUTH
              if(.NOT. present(sbufferx)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD_V: optional argument sbufferx should be presented")
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
           case ( 3 ) ! WEST
              if(.NOT. present(wbufferx)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD_V: optional argument wbufferx should be presented")
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
           case ( 4 ) ! norTH
              if(.NOT. present(nbufferx)) call mpp_error(FATAL,  &
                     "MPP_DO_GET_BOUNDARY_OLD_V: optional argument nbufferx should be presented")
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
           end select
        end if
     end do

  end do

  call mpp_sync_self( )


end subroutine MPP_DO_GET_BOUNDARY_3Dold_V_
