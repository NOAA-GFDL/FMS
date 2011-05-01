function MPP_START_UPDATE_DOMAINS_2D_( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count)
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer                                :: MPP_START_UPDATE_DOMAINS_2D_

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
  pointer( ptr, field3D )
  ptr = LOC(field)

  MPP_START_UPDATE_DOMAINS_2D_ = mpp_start_update_domains(field3D, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count)
  return

end function MPP_START_UPDATE_DOMAINS_2D_

function MPP_START_UPDATE_DOMAINS_3D_( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count )

  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(domain%x(1)%data%begin:,domain%y(1)%data%begin:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer                                :: MPP_START_UPDATE_DOMAINS_3D_

  !--- local variables
  integer                                :: update_whalo, update_ehalo, update_shalo, update_nhalo
  integer                                :: tile, update_position, update_flags
  integer                                :: m, n, i, j, k, msgsize, dir, pos, tMe
  integer                                :: is, ie, js, je, ke
  logical                                :: send(8), recv(8)
  integer                                :: from_pe, to_pe, count
  character(len=128)                     :: text
  type(overlapSpec),  pointer            :: update => NULL()
  type(overlap_type), pointer            :: overPtr => NULL() 
  MPP_TYPE_ :: buffer(size(mpp_domains_stack_nonblock(:)))
  pointer( ptr, buffer )

  if(present(whalo)) then
     update_whalo = whalo
     if(abs(update_whalo) > domain%whalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D: "// &
          "optional argument whalo should not be larger than the whalo when define domain.")
  else
     update_whalo = domain%whalo
  end if
  if(present(ehalo)) then
     update_ehalo = ehalo
     if(abs(update_ehalo) > domain%ehalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D: "// &
          "optional argument ehalo should not be larger than the ehalo when define domain.")
  else
     update_ehalo = domain%ehalo
  end if
  if(present(shalo)) then
     update_shalo = shalo
     if(abs(update_shalo) > domain%shalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D: "// &
          "optional argument shalo should not be larger than the shalo when define domain.")
  else
     update_shalo = domain%shalo
  end if
  if(present(nhalo)) then
     update_nhalo = nhalo
     if(abs(update_nhalo) > domain%nhalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D: "// &
          "optional argument nhalo should not be larger than the nhalo when define domain.")
  else
     update_nhalo = domain%nhalo
  end if
  
  if( .not. domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) ) return

  current_id_update = current_id_update + 1
  num_update = num_update + 1
  if( current_id_update > MAX_DOMAIN_FIELDS ) then
     write( text,'(a,i,a,i)' ) 'num_fields =', current_id_update, ' greater than MAX_DOMAIN_FIELDS =', MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS: '//trim(text))
  endif
  MPP_START_UPDATE_DOMAINS_3D_ = current_id_update

  if( .not. start_update ) then
     call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS: mpp_start_update could not be called before all the update are completed')
  endif

  ke = size(field,3)

  tile = 1
  if(PRESENT(tile_count)) tile = tile_count
  update_position = CENTER
  if(present(position)) update_position = position  
  update_flags = XUPDATE+YUPDATE   !default
  if( PRESENT(flags) )update_flags = flags

  recv(1) = BTEST(update_flags,EAST)
  recv(3) = BTEST(update_flags,SOUTH)
  recv(5) = BTEST(update_flags,WEST)
  recv(7) = BTEST(update_flags,NORTH)
  recv(2) = recv(1) .AND. recv(3)
  recv(4) = recv(3) .AND. recv(5)
  recv(6) = recv(5) .AND. recv(7)
  recv(8) = recv(7) .AND. recv(1)
  send    = recv

  update_flags_list(current_id_update) = update_flags
  update_whalo_list(current_id_update) = update_whalo
  update_ehalo_list(current_id_update) = update_ehalo
  update_shalo_list(current_id_update) = update_shalo
  update_nhalo_list(current_id_update) = update_nhalo
  update_position_list(current_id_update) = update_position

  update => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, update_position)
  ptr = LOC(mpp_domains_stack_nonblock)

  ! pre-postrecv
  do m = 1, update%nrecv
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(recv_clock_nonblock)
     msgsize = 0
     !--- make sure the domain stack size is big enough.
     do n = 1, overPtr%count
        if( overPtr%tileMe(n) .NE. tile ) cycle
        dir = overPtr%dir(n)
        if(recv(dir)) then
           msgsize = msgsize + overPtr%msgsize(n)
        end if
     end do

     msgsize = msgsize*ke
     if( msgsize.GT.0 )then
        from_pe = overPtr%pe
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (nonblock_buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_START_UPDATE_DOMAINS: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        count = request_recv(current_id_update)%count + 1
        if( count > MAX_REQUEST ) then
           write( text,'(a,i,a,i)' ) 'request count =', count, ' greater than MAX_REQEUST =', MAX_REQUEST
           call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS: '//trim(text))
        endif
        request_recv(current_id_update)%count = count
        call mpp_recv( buffer(nonblock_buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=NONBLOCK_UPDATE_TAG, &
             request=request_recv(current_id_update)%request(count))
        nonblock_buffer_pos = nonblock_buffer_pos + msgsize
     end if
     call mpp_clock_end(recv_clock_nonblock)
  end do ! end do m = 1, update%nrecv  

  recv_pos_list(current_id_update) = nonblock_buffer_pos  

  ! send
  do m = 1, update%nsend
     overPtr => update%send(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(pack_clock_nonblock)
     pos = nonblock_buffer_pos

     ! make sure the stacksize is big enough
     msgsize = 0
     do n = 1, overPtr%count
        if( overPtr%tileMe(n) .NE. tile ) cycle
        dir = overPtr%dir(n)
        if( send(dir) )  msgsize = msgsize + overPtr%msgsize(n)
     enddo
     if( msgsize.GT.0 )then
        msgsize = msgsize*ke
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_START_UPDATE_DOMAINS: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if

     do n = 1, overPtr%count
        tMe = overPtr%tileMe(n)
        if( tMe .NE. tile ) cycle
        dir = overPtr%dir(n)
        if( send(dir) ) then
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           if( overptr%is_refined(n) ) then
              do k = 1,ke  
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              end do
           else
              select case( overPtr%rotation(n) )
              case(ZERO)
                 do k = 1,ke  
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              case( MINUS_NINETY ) 
                 do k = 1,ke  
                    do i = is, ie
                       do j = je, js, -1
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              case( NINETY ) 
                 do k = 1,ke  
                    do i = ie, is, -1
                       do j = js, je
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              case( ONE_HUNDRED_EIGHTY ) 
                 do k = 1,ke  
                    do j = je, js, -1
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end select
           end if
        endif
     end do ! do n = 1, overPtr%count

     call mpp_clock_end(pack_clock_nonblock)
     call mpp_clock_begin(send_clock_nonblock)
     msgsize = pos - nonblock_buffer_pos
     if( msgsize.GT.0 )then
        to_pe = overPtr%pe
        call mpp_send( buffer(nonblock_buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=NONBLOCK_UPDATE_TAG)
        nonblock_buffer_pos = pos
     end if
     call mpp_clock_end(send_clock_nonblock)
  end do ! end do ist = 0,nlist-1

  update => NULL()
  overPtr => NULL()

  return

end function MPP_START_UPDATE_DOMAINS_3D_

!##########################################################################################
function MPP_START_UPDATE_DOMAINS_4D_( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count )
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer                                :: MPP_START_UPDATE_DOMAINS_4D_

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
  pointer( ptr, field3D )
  ptr = LOC(field)

  MPP_START_UPDATE_DOMAINS_4D_ = mpp_start_update_domains(field3D, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count)
  return

end function MPP_START_UPDATE_DOMAINS_4D_

!##########################################################################################
function MPP_START_UPDATE_DOMAINS_5D_( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count)
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer                                :: MPP_START_UPDATE_DOMAINS_5D_

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
  pointer( ptr, field3D )
  ptr = LOC(field)

  MPP_START_UPDATE_DOMAINS_5D_ = mpp_start_update_domains(field3D, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count )
  return

end function MPP_START_UPDATE_DOMAINS_5D_

!##################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_2D_( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, buffer )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  MPP_TYPE_,     intent(inout), optional :: buffer(:)

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
  pointer( ptr, field3D )
  ptr = LOC(field)
  call mpp_complete_update_domains(id_update, field3D, domain, flags, position, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count, buffer )

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_2D_

!##################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_3D_( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, buffer )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(domain%x(1)%data%begin:,domain%y(1)%data%begin:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  MPP_TYPE_,     intent(inout), optional :: buffer(:)
  integer                                :: update_whalo, update_ehalo, update_shalo, update_nhalo
  integer                                :: tile, update_position, update_flags
  integer                                :: is, ie, js, je, ke, buffer_pos, pos
  integer                                :: is1, ie1, js1, je1, ni, nj, total
  integer                                :: start, start1, start2
  integer                                :: i, j, k, m, n, tMe, dir
  integer                                :: msgsize, index
  logical                                :: recv(8)
  type(overlapSpec),  pointer            :: update => NULL()
  type(overlap_type), pointer            :: overPtr => NULL() 
  MPP_TYPE_ :: recv_buffer(size(mpp_domains_stack_nonblock(:)))
  pointer( ptr, recv_buffer )

  if(present(whalo)) then
     update_whalo = whalo
     if(abs(update_whalo) > domain%whalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
          "optional argument whalo should not be larger than the whalo when define domain.")
  else
     update_whalo = domain%whalo
  end if
  if(present(ehalo)) then
     update_ehalo = ehalo
     if(abs(update_ehalo) > domain%ehalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
          "optional argument ehalo should not be larger than the ehalo when define domain.")
  else
     update_ehalo = domain%ehalo
  end if
  if(present(shalo)) then
     update_shalo = shalo
     if(abs(update_shalo) > domain%shalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
          "optional argument shalo should not be larger than the shalo when define domain.")
  else
     update_shalo = domain%shalo
  end if
  if(present(nhalo)) then
     update_nhalo = nhalo
     if(abs(update_nhalo) > domain%nhalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
          "optional argument nhalo should not be larger than the nhalo when define domain.")
  else
     update_nhalo = domain%nhalo
  end if
  if( .not. domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) ) return

  call mpp_clock_begin(wait_clock_nonblock)
  if(request_recv(id_update)%count > 0) then
     call mpp_sync_self(request=request_recv(id_update)%request(1:request_recv(id_update)%count))
     request_recv(id_update)%count = 0
  endif 
  call mpp_clock_end(wait_clock_nonblock)

  ke = size(field,3)
  tile = 1
  if(PRESENT(tile_count)) tile = tile_count
  update_position = CENTER
  if(present(position)) update_position = position  
  update_flags = XUPDATE+YUPDATE   !default
  if( PRESENT(flags) )update_flags = flags
  recv(1) = BTEST(update_flags,EAST)
  recv(3) = BTEST(update_flags,SOUTH)
  recv(5) = BTEST(update_flags,WEST)
  recv(7) = BTEST(update_flags,NORTH)
  recv(2) = recv(1) .AND. recv(3)
  recv(4) = recv(3) .AND. recv(5)
  recv(6) = recv(5) .AND. recv(7)
  recv(8) = recv(7) .AND. recv(1)

  !check to make sure the consistency of halo size, position and flags.
  if( update_flags_list(id_update) .NE. update_flags ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
       "mismatch of optional argument flag between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_whalo_list(id_update) .NE. update_whalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
       "mismatch of optional argument whalo between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_ehalo_list(id_update) .NE. update_ehalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
       "mismatch of optional argument ehalo between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_shalo_list(id_update) .NE. update_shalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
       "mismatch of optional argument shalo between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_nhalo_list(id_update) .NE. update_nhalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
       "mismatch of optional argument nhalo between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_position_list(id_update) .NE. update_position ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D: "// &
       "mismatch of optional argument position between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")

  update => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, update_position)
  ptr = LOC(mpp_domains_stack_nonblock)

  buffer_pos = recv_pos_list(id_update)
  !--unpack the data
  call mpp_clock_begin(unpk_clock_nonblock)
  do m = update%nrecv, 1, -1
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle

     pos = buffer_pos
     do n = overPtr%count, 1, -1
        tMe = overPtr%tileMe(n)
        if( tMe .NE. tile ) cycle
        dir = overPtr%dir(n)
        if( recv(dir) ) then
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke
           pos = buffer_pos - msgsize
           buffer_pos = pos
           if(OverPtr%is_refined(n)) then
              if(.not. present(buffer)) then
                 call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS: refinement exists but argument buffer is not present")
              endif
              index = overPtr%index(n)
              is1 = update%rSpec(tMe)%isNbr(index); ie1 = update%rSpec(tMe)%ieNbr(index)
              js1 = update%rSpec(tMe)%jsNbr(index); je1 = update%rSpec(tMe)%jeNbr(index)
              ni = ie1 - is1 + 1
              nj = je1 - js1 + 1
              total = ni*nj
              start = (update%rSpec(tMe)%start(index)-1)*ke
              if(start+total*ke>size(buffer) ) call mpp_error(FATAL, &
                   "MPP_COMPETE_UPDATE_DOMAINS: b_size is less than the size of the data to be filled.")
              msgsize = ie - is + 1
              start1 = start + (js-js1)*ni + is - is1
              do k = 1, ke
                 start2 = start1
                 do j = js, je
                    buffer(start2+1:start2+msgsize) = recv_buffer(pos+1:pos+msgsize)
                    start2 = start2 + ni
                    pos   = pos + msgsize
                 end do
                 start1 = start1 + total
              end do
           else
              do k = 1,ke
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       field(i,j,k) = recv_buffer(pos)
                    end do
                 end do
              end do
           endif
        end if
     end do ! do n = 1, overPtr%count
  end do

  call mpp_clock_end(unpk_clock_nonblock)

  num_update = num_update - 1

  !--- For the last call of mpp_complete_update_domains, need to call mpp_sync_self for send
  !--- also set start_update = .true., reset current_id_update to 0
  if( num_update == 0) then
     start_update      = .true.
     current_id_update = 0
     nonblock_buffer_pos   = 0
     call mpp_clock_begin(wait_clock_nonblock)
     call mpp_sync_self( )
     call mpp_clock_end(wait_clock_nonblock)
  endif

  update => NULL()
  overPtr => NULL()

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_3D_

!##################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_4D_( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, buffer )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  MPP_TYPE_,     intent(inout), optional :: buffer(:)

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
  pointer( ptr, field3D )
  ptr = LOC(field)
  call mpp_complete_update_domains(id_update, field3D, domain, flags, position, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count, buffer )

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_4D_

!##################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_5D_( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, buffer )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  MPP_TYPE_,     intent(inout), optional :: buffer(:)

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
  pointer( ptr, field3D )
  ptr = LOC(field)
  call mpp_complete_update_domains(id_update, field3D, domain, flags, position, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count, buffer )

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_5D_

#ifdef  VECTOR_FIELD_
function MPP_START_UPDATE_DOMAINS_2D_V_( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count )
  !updates data domain of 3D field whose computational domains have been computed
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:), fieldy(:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer                                :: MPP_START_UPDATE_DOMAINS_2D_V_
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  MPP_START_UPDATE_DOMAINS_2D_V_ = mpp_start_update_domains(field3Dx, field3Dy, domain, flags, gridtype, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count )

  return

end function MPP_START_UPDATE_DOMAINS_2D_V_

!###################################################################################
function MPP_START_UPDATE_DOMAINS_3D_V_( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count )
  !updates data domain of 3D field whose computational domains have been computed
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:), fieldy(:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  !--- local variables
  integer                                :: MPP_START_UPDATE_DOMAINS_3D_V_
  integer                                :: update_whalo, update_ehalo, update_shalo, update_nhalo
  integer                                :: grid_offset_type, position_x, position_y, update_flags
  logical                                :: exchange_uv
  character(len=128)                       :: text
  type(overlapSpec),  pointer            :: updatex => NULL()
  type(overlapSpec),  pointer            :: updatey => NULL()

  if(present(whalo)) then
     update_whalo = whalo
     if(abs(update_whalo) > domain%whalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D_V: "// &
          "optional argument whalo should not be larger than the whalo when define domain.")
  else
     update_whalo = domain%whalo
  end if
  if(present(ehalo)) then
     update_ehalo = ehalo
     if(abs(update_ehalo) > domain%ehalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D_V: "// &
          "optional argument ehalo should not be larger than the ehalo when define domain.")
  else
     update_ehalo = domain%ehalo
  end if
  if(present(shalo)) then
     update_shalo = shalo
     if(abs(update_shalo) > domain%shalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D_V: "// &
          "optional argument shalo should not be larger than the shalo when define domain.")
  else
     update_shalo = domain%shalo
  end if
  if(present(nhalo)) then
     update_nhalo = nhalo
     if(abs(update_nhalo) > domain%nhalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D_V: "// &
          "optional argument nhalo should not be larger than the nhalo when define domain.")
  else
     update_nhalo = domain%nhalo
  end if
  
  if( .not. domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) ) return
  current_id_update = current_id_update + 1
  num_update = num_update + 1
  if( current_id_update > MAX_DOMAIN_FIELDS ) then
     write( text,'(a,i,a,i)' ) 'num_fields =', current_id_update, ' greater than MAX_DOMAIN_FIELDS =', MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: '//trim(text) )
  endif
  MPP_START_UPDATE_DOMAINS_3D_V_ = current_id_update

  if( .not. start_update ) then
     call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: mpp_start_update could not be called before all the update are completed')
  endif

  grid_offset_type = AGRID
  if( PRESENT(gridtype) ) grid_offset_type = gridtype

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

  if( BTEST(update_flags,NORTH) .AND. BTEST(domain%fold,NORTH) .AND. BTEST(grid_offset_type,SOUTH) ) &
       call mpp_error( FATAL, 'MPP_START_UPDATE_DOMAINS_V: Incompatible grid offset and fold.' )

  update_flags_list(current_id_update) = update_flags
  update_whalo_list(current_id_update) = update_whalo
  update_ehalo_list(current_id_update) = update_ehalo
  update_shalo_list(current_id_update) = update_shalo
  update_nhalo_list(current_id_update) = update_nhalo
  update_gridtype_list(current_id_update) = grid_offset_type

  exchange_uv = .false.
  if(grid_offset_type == DGRID_NE) then
     exchange_uv = .true.
     grid_offset_type = CGRID_NE
  else if( grid_offset_type == DGRID_SW ) then
     exchange_uv = .true.
     grid_offset_type = CGRID_SW
  end if

  select case(grid_offset_type)
  case (AGRID)
     position_x = CENTER
     position_y = CENTER
  case (BGRID_NE, BGRID_SW)
     position_x = CORNER
     position_y = CORNER
  case (CGRID_NE, CGRID_SW)
     position_x = EAST
     position_y = NORTH
  case default
     call mpp_error(FATAL, "mpp_update_domains2D.h: invalid value of grid_offset_type")
  end select
  updatex => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, position_x)
  updatey => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, position_y)

  if(exchange_uv) then
     call mpp_start_do_update(current_id_update, fieldx, fieldy, domain, updatey, updatex, grid_offset_type, update_flags, &
                                  tile_count)
  else
     call mpp_start_do_update(current_id_update, fieldx, fieldy, domain, updatex, updatey, grid_offset_type, update_flags, &
                                  tile_count)     
  endif 

end function MPP_START_UPDATE_DOMAINS_3D_V_

function MPP_START_UPDATE_DOMAINS_4D_V_( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count )
  !updates data domain of 3D field whose computational domains have been computed
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer                                :: MPP_START_UPDATE_DOMAINS_4D_V_
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  MPP_START_UPDATE_DOMAINS_4D_V_ = mpp_start_update_domains(field3Dx, field3Dy, domain, flags, gridtype, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count )

  return

end function MPP_START_UPDATE_DOMAINS_4D_V_

function MPP_START_UPDATE_DOMAINS_5D_V_( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count )
  !updates data domain of 3D field whose computational domains have been computed
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:,:), fieldy(:,:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer                                :: MPP_START_UPDATE_DOMAINS_5D_V_
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4)*size(fieldx,5))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4)*size(fieldy,5))
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  MPP_START_UPDATE_DOMAINS_5D_V_ = mpp_start_update_domains(field3Dx, field3Dy, domain, flags, gridtype, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count )

  return

end function MPP_START_UPDATE_DOMAINS_5D_V_

!####################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_2D_V_( id_update, fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:), fieldy(:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  MPP_TYPE_,     intent(inout), optional :: bufferx(:), buffery(:)
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  call mpp_complete_update_domains(id_update, field3Dx, field3Dy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery )

  return

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_2D_V_

!####################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_3D_V_( id_update, fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:), fieldy(:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  MPP_TYPE_,     intent(inout), optional :: bufferx(:), buffery(:)

  integer                                :: update_whalo, update_ehalo, update_shalo, update_nhalo
  integer                                :: grid_offset_type, position_x, position_y, update_flags
  logical                                :: exchange_uv
  type(overlapSpec),  pointer            :: updatex => NULL()
  type(overlapSpec),  pointer            :: updatey => NULL()

  if(present(whalo)) then
     update_whalo = whalo
     if(abs(update_whalo) > domain%whalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D_V: "// &
          "optional argument whalo should not be larger than the whalo when define domain.")
  else
     update_whalo = domain%whalo
  end if
  if(present(ehalo)) then
     update_ehalo = ehalo
     if(abs(update_ehalo) > domain%ehalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D_V: "// &
          "optional argument ehalo should not be larger than the ehalo when define domain.")
  else
     update_ehalo = domain%ehalo
  end if
  if(present(shalo)) then
     update_shalo = shalo
     if(abs(update_shalo) > domain%shalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D_V: "// &
          "optional argument shalo should not be larger than the shalo when define domain.")
  else
     update_shalo = domain%shalo
  end if
  if(present(nhalo)) then
     update_nhalo = nhalo
     if(abs(update_nhalo) > domain%nhalo ) call mpp_error(FATAL, "MPP_START_UPDATE_DOMAINS_3D_V: "// &
          "optional argument nhalo should not be larger than the nhalo when define domain.")
  else
     update_nhalo = domain%nhalo
  end if
  
  if( .not. domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) ) return

  call mpp_clock_begin(wait_clock_nonblock)
  if(request_recv(id_update)%count > 0) then
     call mpp_sync_self(request=request_recv(id_update)%request(1:request_recv(id_update)%count))
     request_recv(id_update)%count = 0
  endif
  call mpp_clock_end(wait_clock_nonblock)

  grid_offset_type = AGRID
  if( PRESENT(gridtype) ) grid_offset_type = gridtype

  update_flags = XUPDATE+YUPDATE   !default
  if( PRESENT(flags) )update_flags = flags

  !check to make sure the consistency of halo size, position and flags.
  if( update_flags_list(id_update) .NE. update_flags ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D_V: "// &
       "mismatch of optional argument flag between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_whalo_list(id_update) .NE. update_whalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D_V: "// &
       "mismatch of optional argument whalo between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_ehalo_list(id_update) .NE. update_ehalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D_V: "// &
       "mismatch of optional argument ehalo between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_shalo_list(id_update) .NE. update_shalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D_V: "// &
       "mismatch of optional argument shalo between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_nhalo_list(id_update) .NE. update_nhalo ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D_V: "// &
       "mismatch of optional argument nhalo between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")
  if( update_gridtype_list(id_update) .NE. grid_offset_type ) call mpp_error(FATAL, "MPP_COMPLETE_UPDATE_DOMAINS_3D_V: "// &
       "mismatch of optional argument gridtype between MPP_COMPLETE_UPDATE_DOMAINS and MPP_START_UPDATE_DOMAINS")

  exchange_uv = .false.
  if(grid_offset_type == DGRID_NE) then
     exchange_uv = .true.
     grid_offset_type = CGRID_NE
  else if( grid_offset_type == DGRID_SW ) then
     exchange_uv = .true.
     grid_offset_type = CGRID_SW
  end if

  select case(grid_offset_type)
  case (AGRID)
     position_x = CENTER
     position_y = CENTER
  case (BGRID_NE, BGRID_SW)
     position_x = CORNER
     position_y = CORNER
  case (CGRID_NE, CGRID_SW)
     position_x = EAST
     position_y = NORTH
  case default
     call mpp_error(FATAL, "mpp_update_domains2D.h: invalid value of grid_offset_type")
  end select
  updatex => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, position_x)
  updatey => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, position_y)

  if(exchange_uv) then
     call mpp_complete_do_update(id_update, fieldx, fieldy, domain, updatey, updatex, grid_offset_type, update_flags, &
                                 tile_count, buffery, bufferx)
  else
     call mpp_complete_do_update(id_update, fieldx, fieldy, domain, updatex, updatey, grid_offset_type, update_flags, &
                                  tile_count, bufferx, buffery)     
  endif 

  num_update = num_update - 1

  !--- For the last call of mpp_complete_update_domains, need to call mpp_sync_self for send
  !--- also set start_update = .true., reset current_id_update to 0
  if( num_update == 0) then
     start_update      = .true.
     current_id_update = 0
     nonblock_buffer_pos   = 0
     call mpp_clock_begin(wait_clock_nonblock)
     call mpp_sync_self( )
     call mpp_clock_end(wait_clock_nonblock)
  endif

  updatex => NULL()
  updatey => NULL()

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_3D_V_

!####################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_4D_V_( id_update, fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  MPP_TYPE_,     intent(inout), optional :: bufferx(:), buffery(:)
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  call mpp_complete_update_domains(id_update, field3Dx, field3Dy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery )

  return

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_4D_V_

!####################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_5D_V_( id_update, fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:,:), fieldy(:,:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  MPP_TYPE_,     intent(inout), optional :: bufferx(:), buffery(:)
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4)*size(fieldx,5))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4)*size(fieldy,5))
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  call mpp_complete_update_domains(id_update, field3Dx, field3Dy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery )

  return

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_5D_V_

#endif
