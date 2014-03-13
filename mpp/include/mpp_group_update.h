! -*-f90-*-
subroutine MPP_CREATE_GROUP_UPDATE_2D_(group, field, domain, flags, position, &
     whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,            intent(inout)        :: field(:,:)
  type(domain2D),       intent(inout)        :: domain  
  integer,              intent(in), optional :: flags
  integer,              intent(in), optional :: position
  integer,              intent(in), optional :: whalo, ehalo, shalo, nhalo

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
  pointer( ptr, field3D )
  ptr = LOC(field)

  call mpp_create_group_update(group, field3D, domain, flags, position, whalo, ehalo, shalo, nhalo)

  return

end subroutine MPP_CREATE_GROUP_UPDATE_2D_

subroutine MPP_CREATE_GROUP_UPDATE_3D_(group, field, domain, flags, position, whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,        intent(inout)        :: field(:,:,:)
  type(domain2D),   intent(inout)        :: domain  
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.

  integer          :: update_position, update_whalo, update_ehalo, update_shalo, update_nhalo
  integer          :: update_flags, isize, jsize, ksize
  integer          :: nscalar
  character(len=3) :: text
  logical          :: set_mismatch, update_edge_only
  logical          :: recv(8)

  if(group%initialized) then
      call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE_3D: group is already initialized")
  endif

  if(present(whalo)) then
     update_whalo = whalo
     if(abs(update_whalo) > domain%whalo ) call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE: "// &
          "optional argument whalo should not be larger than the whalo when define domain.")
  else
     update_whalo = domain%whalo
  end if
  if(present(ehalo)) then
     update_ehalo = ehalo
     if(abs(update_ehalo) > domain%ehalo ) call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE: "// &
          "optional argument ehalo should not be larger than the ehalo when define domain.")
  else
     update_ehalo = domain%ehalo
  end if
  if(present(shalo)) then
     update_shalo = shalo
     if(abs(update_shalo) > domain%shalo ) call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE: "// &
          "optional argument shalo should not be larger than the shalo when define domain.")
  else
     update_shalo = domain%shalo
  end if
  if(present(nhalo)) then
     update_nhalo = nhalo
     if(abs(update_nhalo) > domain%nhalo ) call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE: "// &
          "optional argument nhalo should not be larger than the nhalo when define domain.")
  else
     update_nhalo = domain%nhalo
  end if
  update_position = CENTER
  !--- when there is NINETY or MINUS_NINETY rotation for some contact, the salar data can not be on E or N-cell,
  if(present(position)) then
     update_position = position 
     if(domain%rotated_ninety .AND. ( position == EAST .OR. position == NORTH ) )  &
          call mpp_error(FATAL, 'MPP_CREATE_GROUP_UPDATE_3D: hen there is NINETY or MINUS_NINETY rotation, ' // &
          'can not use scalar version update_domain for data on E or N-cell' )
  end if

  if( domain%max_ntile_pe > 1 ) then
     call mpp_error(FATAL,'MPP_CREATE_GROUP_UPDATE: do not support multiple tile per processor')
  endif

  update_flags = XUPDATE+YUPDATE 
  if(present(flags)) update_flags = flags

  group%nscalar = group%nscalar + 1
  nscalar = group%nscalar
  if( nscalar > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_CREATE_GROUP_UPDATE: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif

  isize = size(field,1); jsize=size(field,2); ksize = size(field,3)

  group%addrs_s(nscalar) = LOC(field)
  if( group%nscalar == 1 ) then
     group%flags_s  = update_flags
     group%whalo_s  = update_whalo
     group%ehalo_s  = update_ehalo
     group%shalo_s  = update_shalo
     group%nhalo_s  = update_nhalo
     group%position = update_position
     group%isize_s  = isize
     group%jsize_s  = jsize
     group%ksize_s  = ksize
     call mpp_get_memory_domain(domain, group%is_s, group%ie_s, group%js_s, group%je_s, position=position)

     update_edge_only = BTEST(update_flags, EDGEONLY)
     recv(1) = BTEST(update_flags,EAST)
     recv(3) = BTEST(update_flags,SOUTH)
     recv(5) = BTEST(update_flags,WEST)
     recv(7) = BTEST(update_flags,NORTH)
     if( update_edge_only ) then
        recv(2) = .false.
        recv(4) = .false.
        recv(6) = .false.
        recv(8) = .false.
        if( .NOT. (recv(1) .OR. recv(3) .OR. recv(5) .OR. recv(7)) ) then
           recv(1) = .true.
           recv(3) = .true.
           recv(5) = .true.
           recv(7) = .true.
        endif
     else
        recv(2) = recv(1) .AND. recv(3)
        recv(4) = recv(3) .AND. recv(5)
        recv(6) = recv(5) .AND. recv(7)
        recv(8) = recv(7) .AND. recv(1)
     endif
     group%recv_s = recv
     group%update_s => search_update_overlap(domain, group%whalo_s, group%ehalo_s, &
          group%shalo_s, group%nhalo_s, group%position)
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (group%flags_s  .NE. update_flags)
     set_mismatch = set_mismatch .OR. (group%whalo_s  .NE. update_whalo)
     set_mismatch = set_mismatch .OR. (group%ehalo_s  .NE. update_ehalo)
     set_mismatch = set_mismatch .OR. (group%shalo_s  .NE. update_shalo)
     set_mismatch = set_mismatch .OR. (group%nhalo_s  .NE. update_nhalo)
     set_mismatch = set_mismatch .OR. (group%position .NE. update_position)
     set_mismatch = set_mismatch .OR. (group%isize_s  .NE. isize)
     set_mismatch = set_mismatch .OR. (group%jsize_s  .NE. jsize)
     set_mismatch = set_mismatch .OR. (group%ksize_s  .NE. ksize)

     if(set_mismatch)then
        write( text,'(i2)' ) nscalar
        call mpp_error(FATAL,'MPP_CREATE_GROUP_UPDATE_3D: Incompatible field at count '//text//' for group update.' )
     endif
  endif

  return

end subroutine MPP_CREATE_GROUP_UPDATE_3D_


subroutine MPP_CREATE_GROUP_UPDATE_4D_(group, field, domain, flags, position, &
     whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,            intent(inout)        :: field(:,:,:,:)
  type(domain2D),       intent(inout)        :: domain  
  integer,              intent(in), optional :: flags
  integer,              intent(in), optional :: position
  integer,              intent(in), optional :: whalo, ehalo, shalo, nhalo

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
  pointer( ptr, field3D )
  ptr = LOC(field)

  call mpp_create_group_update(group, field3D, domain, flags, position, whalo, ehalo, shalo, nhalo)

  return

end subroutine MPP_CREATE_GROUP_UPDATE_4D_

subroutine MPP_CREATE_GROUP_UPDATE_2D_V_( group, fieldx, fieldy, domain, flags, gridtype, &
                                          whalo, ehalo, shalo, nhalo)

  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,        intent(inout)            :: fieldx(:,:), fieldy(:,:)
  type(domain2D),   intent(inout)            :: domain
  integer,          intent(in),     optional :: flags, gridtype
  integer,          intent(in),     optional :: whalo, ehalo, shalo, nhalo
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  call  mpp_create_group_update(group, field3Dx, field3Dy, domain, flags, gridtype, &
                                          whalo, ehalo, shalo, nhalo)

  return

end subroutine MPP_CREATE_GROUP_UPDATE_2D_V_



subroutine MPP_CREATE_GROUP_UPDATE_3D_V_( group, fieldx, fieldy, domain, flags, gridtype, &
                                          whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:), fieldy(:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo

  integer          :: update_whalo, update_ehalo, update_shalo, update_nhalo
  integer          :: update_flags, isize_x, jsize_x, ksize_x, isize_y, jsize_y, ksize_y
  integer          :: nvector, update_gridtype, position_x, position_y
  character(len=3) :: text
  logical          :: set_mismatch, update_edge_only
  logical          :: recv(8)


  if(group%initialized) then
     call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE_V: group is already initialized")
  endif

  if(present(whalo)) then
     update_whalo = whalo
     if(abs(update_whalo) > domain%whalo ) call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE_V: "// &
          "optional argument whalo should not be larger than the whalo when define domain.")
  else
     update_whalo = domain%whalo
  end if
  if(present(ehalo)) then
     update_ehalo = ehalo
     if(abs(update_ehalo) > domain%ehalo ) call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE_V: "// &
          "optional argument ehalo should not be larger than the ehalo when define domain.")
  else
     update_ehalo = domain%ehalo
  end if
  if(present(shalo)) then
     update_shalo = shalo
     if(abs(update_shalo) > domain%shalo ) call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE_V: "// &
          "optional argument shalo should not be larger than the shalo when define domain.")
  else
     update_shalo = domain%shalo
  end if
  if(present(nhalo)) then
     update_nhalo = nhalo
     if(abs(update_nhalo) > domain%nhalo ) call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE_V: "// &
          "optional argument nhalo should not be larger than the nhalo when define domain.")
  else
     update_nhalo = domain%nhalo
  end if

  update_gridtype = AGRID
  if(PRESENT(gridtype)) update_gridtype = gridtype

  if( domain%max_ntile_pe > 1 ) then
     call mpp_error(FATAL,'MPP_CREATE_GROUP_UPDATE_V: do not support multiple tile per processor')
  endif

  update_flags = XUPDATE+YUPDATE   !default
  if( PRESENT(flags) )update_flags = flags

  group%nvector = group%nvector + 1
  nvector = group%nvector
  if( nvector > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_CREATE_GROUP_UPDATE_V: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif

  isize_x = size(fieldx,1); jsize_x = size(fieldx,2); ksize_x = size(fieldx,3) 
  isize_y = size(fieldy,1); jsize_y = size(fieldy,2); ksize_y = size(fieldy,3) 
  
  if(ksize_x .NE. ksize_y) call mpp_error(FATAL,  &
           'MPP_CREATE_GROUP_UPDATE_V: mismatch of ksize between fieldx and fieldy')

  group%addrs_x(nvector) = LOC(fieldx)
  group%addrs_y(nvector) = LOC(fieldy)

  if( group%nvector == 1 ) then
     group%flags_v  = update_flags
     group%whalo_v  = update_whalo
     group%ehalo_v  = update_ehalo
     group%shalo_v  = update_shalo
     group%nhalo_v  = update_nhalo
     group%gridtype = update_gridtype
     group%isize_x  = isize_x
     group%jsize_x  = jsize_x
     group%isize_y  = isize_y
     group%jsize_y  = jsize_y
     group%ksize_v  = ksize_x
     update_edge_only = BTEST(update_flags, EDGEONLY)
     recv(1) = BTEST(update_flags,EAST)
     recv(3) = BTEST(update_flags,SOUTH)
     recv(5) = BTEST(update_flags,WEST)
     recv(7) = BTEST(update_flags,NORTH)
     if( update_edge_only ) then
        recv(2) = .false.
        recv(4) = .false.
        recv(6) = .false.
        recv(8) = .false.
        if( .NOT. (recv(1) .OR. recv(3) .OR. recv(5) .OR. recv(7)) ) then
           recv(1) = .true.
           recv(3) = .true.
           recv(5) = .true.
           recv(7) = .true.
        endif
     else
        recv(2) = recv(1) .AND. recv(3)
        recv(4) = recv(3) .AND. recv(5)
        recv(6) = recv(5) .AND. recv(7)
        recv(8) = recv(7) .AND. recv(1)
     endif
     group%recv_v = recv
     select case(group%gridtype)
     case (AGRID)
        position_x = CENTER
        position_y = CENTER
     case (BGRID_NE, BGRID_SW)
        position_x = CORNER
        position_y = CORNER
     case (CGRID_NE, CGRID_SW)
        position_x = EAST
        position_y = NORTH
     case (DGRID_NE, DGRID_SW)
        position_x = NORTH
        position_y = EAST
     case default
        call mpp_error(FATAL, "mpp_CREATE_GROUP_UPDATE_V: invalid value of gridtype")
     end select

     call mpp_get_memory_domain(domain, group%is_x, group%ie_x, group%js_x, group%je_x, position=position_x)
     call mpp_get_memory_domain(domain, group%is_y, group%ie_y, group%js_y, group%je_y, position=position_y)

     group%update_x => search_update_overlap(domain, group%whalo_v, group%ehalo_v, &
          group%shalo_v, group%nhalo_v, position_x)
     group%update_y => search_update_overlap(domain, group%whalo_v, group%ehalo_v, &
          group%shalo_v, group%nhalo_v, position_y)
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (group%flags_v  .NE. update_flags)
     set_mismatch = set_mismatch .OR. (group%whalo_v  .NE. update_whalo)
     set_mismatch = set_mismatch .OR. (group%ehalo_v  .NE. update_ehalo)
     set_mismatch = set_mismatch .OR. (group%shalo_v  .NE. update_shalo)
     set_mismatch = set_mismatch .OR. (group%nhalo_v  .NE. update_nhalo)
     set_mismatch = set_mismatch .OR. (group%gridtype .NE. update_gridtype)
     set_mismatch = set_mismatch .OR. (group%isize_x  .NE. isize_x)
     set_mismatch = set_mismatch .OR. (group%jsize_x  .NE. jsize_x)
     set_mismatch = set_mismatch .OR. (group%isize_y  .NE. isize_y)
     set_mismatch = set_mismatch .OR. (group%jsize_y  .NE. jsize_y)
     set_mismatch = set_mismatch .OR. (group%ksize_v  .NE. ksize_x)

     if(set_mismatch)then
        write( text,'(i2)' ) nvector
        call mpp_error(FATAL,'MPP_CREATE_GROUP_UPDATE_V: Incompatible field at count '//text//' for group update.' )
     endif
  endif

  return

end subroutine MPP_CREATE_GROUP_UPDATE_3D_V_

subroutine MPP_CREATE_GROUP_UPDATE_4D_V_( group, fieldx, fieldy, domain, flags, gridtype, &
                                          whalo, ehalo, shalo, nhalo)

  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,        intent(inout)            :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  type(domain2D),   intent(inout)            :: domain
  integer,          intent(in),     optional :: flags, gridtype
  integer,          intent(in),     optional :: whalo, ehalo, shalo, nhalo
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  call  mpp_create_group_update(group, field3Dx, field3Dy, domain, flags, gridtype, &
                                          whalo, ehalo, shalo, nhalo)

  return

end subroutine MPP_CREATE_GROUP_UPDATE_4D_V_


subroutine MPP_DO_GROUP_UPDATE_(group, domain, d_type)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D),              intent(inout) :: domain  
  MPP_TYPE_,                   intent(in)    :: d_type

  integer   :: nscalar, nvector, nlist
  integer   :: nsend, nrecv, flags_v, ntot, nsend_old, nrecv_old
  integer   :: nsend_s, nsend_x, nsend_y
  integer   :: nrecv_s, nrecv_x, nrecv_y
  integer   :: update_buffer_pos, tot_recv_size, tot_send_size
  integer   :: msgsize_s, msgsize_x, msgsize_y, msgsize
  logical   :: recv_s(8), send_s(8), recv_v(8), send_v(8)
  integer   :: i_s, i_x, i_y, rank_s, rank_x, rank_y, rank
  integer   :: from_pe, to_pe, buffer_pos, pos, dir
  integer   :: ksize, is, ie, js, je
  integer   :: n, l, m, t, i, j, k, tk
  integer   :: shift, gridtype, midpoint
  integer   :: ind_s(3*MAXOVERLAP), recv_ind_s(3*MAXOVERLAP), send_ind_s(3*MAXOVERLAP)
  integer   :: ind_x(3*MAXOVERLAP), recv_ind_x(3*MAXOVERLAP), send_ind_x(3*MAXOVERLAP)
  integer   :: ind_y(3*MAXOVERLAP), recv_ind_y(3*MAXOVERLAP), send_ind_y(3*MAXOVERLAP)
  integer   :: pelist(3*MAXOVERLAP), from_pe_list(3*MAXOVERLAP), to_pe_list(3*MAXOVERLAP)
  integer   :: recv_size(3*MAXOVERLAP), send_size(3*MAXOVERLAP)
  integer   :: buffer_pos_recv(3*MAXOVERLAP), buffer_pos_send(3*MAXOVERLAP)
  character(len=8)            :: text
  type(overlap_type), pointer :: overPtr => NULL()

  MPP_TYPE_ :: sbuffer(mpp_domains_stack_size)
  MPP_TYPE_ :: rbuffer(mpp_domains_stack_size)
  MPP_TYPE_ :: field (group%is_s:group%ie_s,group%js_s:group%je_s, group%ksize_s)
  MPP_TYPE_ :: fieldx(group%is_x:group%ie_x,group%js_x:group%je_x, group%ksize_v)
  MPP_TYPE_ :: fieldy(group%is_y:group%ie_y,group%js_y:group%je_y, group%ksize_v)
  pointer(ptr_field, field)
  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)

  nscalar = group%nscalar
  nvector = group%nvector
  nlist   = size(domain%list(:))
  gridtype = group%gridtype

  !--- ksize_s must equal ksize_v 
  if(nvector > 0 .AND. nscalar > 0) then
     if(group%ksize_s .NE. group%ksize_v) then
        call mpp_error(FATAL, "MPP_DO_GROUP_UPDATE: ksize_s and ksize_v are not equal")
     endif
     ksize = group%ksize_s
  else if (nscalar > 0) then
     ksize = group%ksize_s
  else if (nvector > 0) then
     ksize = group%ksize_v
  else
     call mpp_error(FATAL, "MPP_DO_GROUP_UPDATE: nscalar and nvector are all 0")
  endif

  if(nscalar > 0) then
     recv_s = group%recv_s
     send_s = recv_s
  endif
  if(nvector > 0) then
     recv_v = group%recv_v
     send_v = recv_v
  end if

  if(.not. group%initialized) then
     group%initialized = .true.
     nsend_s = 0; nsend_x = 0; nsend_y = 0
     nrecv_s = 0; nrecv_x = 0; nrecv_y = 0
     if(nscalar > 0) then
         !--- This check could not be done because of memory domain
!        if( group%isize_s .NE. (group%ie_s-group%is_s+1) .OR. group%jsize_s .NE. (group%je_s-group%js_s+1))  &
!           call mpp_error(FATAL, "MPP_DO_GROUP_UPDATE: mismatch of size of the field and domain memory domain")
        nsend_s = group%update_s%nsend
        nrecv_s = group%update_s%nrecv
     endif

     if(nvector > 0) then
        !--- This check could not be done because of memory domain
!        if( group%isize_x .NE. (group%ie_x-group%is_x+1) .OR. group%jsize_x .NE. (group%je_x-group%js_x+1))  &
!           call mpp_error(FATAL, "MPP_DO_GROUP_UPDATE: mismatch of size of the fieldx and domain memory domain")
!        if( group%isize_y .NE. (group%ie_y-group%is_y+1) .OR. group%jsize_y .NE. (group%je_y-group%js_y+1))  &
!           call mpp_error(FATAL, "MPP_DO_GROUP_UPDATE: mismatch of size of the fieldy and domain memory domain")
        nsend_x = group%update_x%nsend
        nrecv_x = group%update_x%nrecv
        nsend_y = group%update_y%nsend
        nrecv_y = group%update_y%nrecv
     endif

     !figure out message size for each processor.
     ntot = nrecv_s + nrecv_x + nrecv_y
     if(ntot > 3*MAXOVERLAP) call mpp_error(FATAL, "MPP_DO_GROUP_UPDATE: ntot is greater than 3*MAXOVERLAP")
     n = 1
     i_s = 1
     i_x = 1
     i_y = 1
     ind_s = -1
     ind_x = -1
     ind_y = -1
     nrecv = 0
     do while(n<=ntot)
        if( i_s <= nrecv_s ) then
           rank_s = group%update_s%recv(i_s)%pe-domain%pe
           if(rank_s .LE. 0) rank_s = rank_s + nlist
        else
           rank_s = -1
        endif
        if( i_x <= nrecv_x ) then
           rank_x = group%update_x%recv(i_x)%pe-domain%pe
           if(rank_x .LE. 0) rank_x = rank_x + nlist
        else
           rank_x = -1
        endif
        if( i_y <= nrecv_y ) then
           rank_y = group%update_y%recv(i_y)%pe-domain%pe
           if(rank_y .LE. 0) rank_y = rank_y + nlist
        else
           rank_y = -1
        endif
        nrecv = nrecv + 1   
        rank = maxval((/rank_s, rank_x, rank_y/))
        if(rank == rank_s) then
           n = n + 1
           ind_s(nrecv) = i_s
           pelist(nrecv) = group%update_s%recv(i_s)%pe
           i_s = i_s + 1
        endif
        if(rank == rank_x) then
           n = n + 1
           ind_x(nrecv) = i_x
           pelist(nrecv) = group%update_x%recv(i_x)%pe
           i_x = i_x + 1
        endif
        if(rank == rank_y) then
           n = n + 1
           ind_y(nrecv) = i_y
           pelist(nrecv) = group%update_y%recv(i_y)%pe
           i_y = i_y + 1
        endif
     enddo

     nrecv_old = nrecv
     nrecv     = 0
     update_buffer_pos = 0
     tot_recv_size = 0
     do l = 1, nrecv_old
        msgsize_s = 0
        msgsize_x = 0
        msgsize_y = 0
        m = ind_s(l)
        if(m>0) msgsize_s = get_mesgsize(group%update_s%recv(m), recv_s)*ksize*nscalar
        m = ind_x(l)
        if(m>0) msgsize_x = get_mesgsize(group%update_x%recv(m), recv_v)*ksize*nvector
        m = ind_y(l)
        if(m>0) msgsize_y = get_mesgsize(group%update_y%recv(m), recv_v)*ksize*nvector
        msgsize = msgsize_s + msgsize_x + msgsize_y
        if( msgsize.GT.0 )then
           tot_recv_size = tot_recv_size + msgsize
           nrecv = nrecv + 1
           from_pe_list(nrecv) = pelist(l)
           recv_ind_s(nrecv) = ind_s(l)
           recv_ind_x(nrecv) = ind_x(l)
           recv_ind_y(nrecv) = ind_y(l)
           recv_size(nrecv) = msgsize
           buffer_pos_recv(nrecv) = update_buffer_pos
           update_buffer_pos = update_buffer_pos + msgsize
        end if
     end do

     if(nrecv > MAXOVERLAP) then
        call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE: nrecv is greater than MAXOVERLAP, increase MAXOVERLAP")
     endif
     group%nrecv = nrecv
     if(nrecv > 0) then
        group%from_pe(1:nrecv) = from_pe_list(1:nrecv)
        group%buffer_pos_recv(1:nrecv) = buffer_pos_recv(1:nrecv)
        group%recv_size(1:nrecv) = recv_size(1:nrecv)
        group%msgsize_recv(1:nrecv) = recv_size(1:nrecv)/ksize
        group%recv_ind_s(1:nrecv) = ind_s(1:nrecv)
        group%recv_ind_x(1:nrecv) = ind_x(1:nrecv)
        group%recv_ind_y(1:nrecv) = ind_y(1:nrecv)
     endif

     !figure out message size for each processor.
     ntot = nsend_s + nsend_x + nsend_y
     n = 1
     i_s = 1
     i_x = 1
     i_y = 1
     ind_s = -1
     ind_x = -1
     ind_y = -1
     nsend = 0
     do while(n<=ntot)
        if( i_s <= nsend_s ) then
           rank_s = group%update_s%send(i_s)%pe-domain%pe
           if(rank_s .LT. 0) rank_s = rank_s + nlist
        else
           rank_s = nlist+1
        endif
        if( i_x <= nsend_x ) then
           rank_x = group%update_x%send(i_x)%pe-domain%pe
           if(rank_x .LT. 0) rank_x = rank_x + nlist
        else
           rank_x = nlist+1
        endif
        if( i_y <= nsend_y ) then
           rank_y = group%update_y%send(i_y)%pe-domain%pe
           if(rank_y .LT. 0) rank_y = rank_y + nlist
        else
           rank_y = nlist+1
        endif
        nsend = nsend + 1   
        rank = minval((/rank_s, rank_x, rank_y/))
        if(rank == rank_s) then
           n = n + 1
           ind_s(nsend) = i_s
           pelist(nsend) = group%update_s%send(i_s)%pe
           i_s = i_s + 1
        endif
        if(rank == rank_x) then
           n = n + 1
           ind_x(nsend) = i_x
           pelist(nsend) = group%update_x%send(i_x)%pe
           i_x = i_x + 1
        endif
        if(rank == rank_y) then
           n = n + 1
           ind_y(nsend) = i_y
           pelist(nsend) = group%update_y%send(i_y)%pe
           i_y = i_y + 1
        endif
     enddo

     nsend_old = nsend
     nsend     = 0
     update_buffer_pos = 0
     tot_send_size = 0
     do l = 1, nsend_old
        msgsize_s = 0
        msgsize_x = 0
        msgsize_y = 0
        m = ind_s(l)
        if(m>0) msgsize_s = get_mesgsize(group%update_s%send(m), send_s)*ksize*nscalar
        m = ind_x(l)
        if(m>0) msgsize_x = get_mesgsize(group%update_x%send(m), send_v)*ksize*nvector
        m = ind_y(l)
        if(m>0) msgsize_y = get_mesgsize(group%update_y%send(m), send_v)*ksize*nvector
        msgsize = msgsize_s + msgsize_x + msgsize_y
        if( msgsize.GT.0 )then
           tot_send_size = tot_send_size + msgsize
           nsend = nsend + 1
           send_size(nsend) = msgsize
           buffer_pos_send(nsend) = update_buffer_pos
           send_ind_s(nsend) = ind_s(l)
           send_ind_x(nsend) = ind_x(l)
           send_ind_y(nsend) = ind_y(l)
           to_pe_list(nsend) = pelist(l)
           update_buffer_pos = update_buffer_pos + msgsize
        end if
     end do

     if(nsend > MAXOVERLAP) then
        call mpp_error(FATAL, "MPP_CREATE_GROUP_UPDATE: nsend is greater than MAXOVERLAP, increase MAXOVERLAP")
     endif
     group%nsend = nsend
     if(nsend > 0) then
        group%to_pe(1:nsend) = to_pe_list(1:nsend)
        group%buffer_pos_send(1:nsend) = buffer_pos_send(1:nsend)
        group%send_size(1:nsend) = send_size(1:nsend)
        group%msgsize_send(1:nsend) = send_size(1:nsend)/ksize
        group%send_ind_s(1:nsend) = ind_s(1:nsend)
        group%send_ind_x(1:nsend) = ind_x(1:nsend)
        group%send_ind_y(1:nsend) = ind_y(1:nsend)
     endif
     !--- make sure the buffer is large enough
     mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, tot_recv_size )
     mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, tot_send_size )

     if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
        write( text,'(i8)' )mpp_domains_stack_hwm
        call mpp_error( FATAL, 'MPP_DO_GROUP_UPDATE: mpp_domains_stack overflow, '// &
                     'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
     end if
  else
     nrecv = group%nrecv
     nsend = group%nsend
  endif

  !---pre-post receive.
  call mpp_clock_begin(group_recv_clock)
  do m = 1, nrecv
     msgsize = group%recv_size(m)
     from_pe = group%from_pe(m)
     if( msgsize .GT. 0 )then
        buffer_pos = group%buffer_pos_recv(m)
        call mpp_recv( rbuffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., &
             tag=COMM_TAG_1)
     end if
  end do
  call mpp_clock_end(group_recv_clock)

  flags_v = group%flags_v

  call mpp_clock_begin(group_pack_clock)
  !pack the data
  !$OMP parallel do schedule(dynamic) default(shared) private(t, k, buffer_pos, msgsize, pos, m, &
  !$OMP                                                       overPtr, dir, is, ie, js, je,      &
  !$OMP                                                       ptr_field, ptr_fieldx, ptr_fieldy )
  do tk = 1, nsend*ksize
     t = (tk-1)/ksize + 1
     k = mod((tk-1), ksize) + 1
     buffer_pos = group%buffer_pos_send(t)
     msgsize    = group%msgsize_send(t)
     pos  = buffer_pos + (k-1)*msgsize
     m = group%send_ind_s(t)
     if( m > 0 ) then
        overPtr => group%update_s%send(m)
        do n = 1, overPtr%count
           dir = overPtr%dir(n)
           if( send_s(dir) ) then
              is = overPtr%is(n); ie = overPtr%ie(n)
              js = overPtr%js(n); je = overPtr%je(n)
              select case( overPtr%rotation(n) )
              case(ZERO)
                 do l=1, group%nscalar  ! loop over number of fields
                    ptr_field = group%addrs_s(l)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          sbuffer(pos) = field(i,j,k)
                       end do
                    end do
                 enddo
              case( MINUS_NINETY )
                 do l=1,group%nscalar  ! loop over number of fields
                    ptr_field = group%addrs_s(l)
                    do i = is, ie
                       do j = je, js, -1
                          pos = pos + 1
                          sbuffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              case( NINETY )
                 do l=1,group%nscalar  ! loop over number of fields
                    ptr_field = group%addrs_s(l)
                    do i = ie, is, -1
                       do j = js, je
                          pos = pos + 1
                          sbuffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              case( ONE_HUNDRED_EIGHTY )
                 do l=1,group%nscalar  ! loop over number of fields
                    ptr_field = group%addrs_s(l)
                    do j = je, js, -1
                       do i = ie, is, -1
                          pos = pos + 1
                          sbuffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end select
           endif
        end do
     endif

     m = group%send_ind_x(t)
     if( m > 0 ) then
        overPtr => group%update_x%send(m)
        do n = 1, overPtr%count
           dir = overPtr%dir(n)
           if( send_v(dir) ) then
              is = overPtr%is(n); ie = overPtr%ie(n)
              js = overPtr%js(n); je = overPtr%je(n)
              select case( overPtr%rotation(n) )
              case(ZERO)
                 do l=1, nvector  ! loop over number of fields
                    ptr_fieldx = group%addrs_x(l)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          sbuffer(pos) = fieldx(i,j,k)
                       end do
                    end do
                 end do
              case( MINUS_NINETY )
                 if( BTEST(group%flags_v,SCALAR_BIT) ) then
                    do l=1,nvector  ! loop over number of fields
                       ptr_fieldy = group%addrs_y(l)
                       do i = is, ie
                          do j = je, js, -1
                             pos = pos + 1
                             sbuffer(pos) = fieldy(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do l=1,nvector  ! loop over number of fields
                       ptr_fieldy = group%addrs_y(l)
                       do i = is, ie
                          do j = je, js, -1
                             pos = pos + 1
                             sbuffer(pos) = -fieldy(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case( NINETY )
                 do l=1, nvector  ! loop over number of fields
                    ptr_fieldy = group%addrs_y(l)
                    do i = ie, is, -1
                       do j = js, je
                          pos = pos + 1
                          sbuffer(pos) = fieldy(i,j,k)
                       end do
                    end do
                 end do
              case( ONE_HUNDRED_EIGHTY )
                 if( BTEST(group%flags_v,SCALAR_BIT) ) then
                    do l=1,nvector  ! loop over number of fields
                       ptr_fieldx = group%addrs_x(l)
                       do j = je, js, -1
                          do i = ie, is, -1
                             pos = pos + 1
                             sbuffer(pos) =  fieldx(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do l=1,nvector  ! loop over number of fields
                       ptr_fieldx = group%addrs_x(l)
                       do j = je, js, -1
                          do i = ie, is, -1
                             pos = pos + 1
                             sbuffer(pos) =  -fieldx(i,j,k)
                          end do
                       end do
                    end do
                 end if
              end select ! select case( rotation(n) )
           end if ! if( send(dir) )
        end do ! do n = 1, update_x%send(m)%count
     endif

     m = group%send_ind_y(t)
     if( m > 0 ) then
        overPtr => group%update_y%send(m)
        do n = 1, overPtr%count
           dir = overPtr%dir(n)
           if( send_v(dir) ) then
              is = overPtr%is(n); ie = overPtr%ie(n)
              js = overPtr%js(n); je = overPtr%je(n)
              select case( overPtr%rotation(n) )
              case(ZERO)
                 do l=1, nvector  ! loop over number of fields
                    ptr_fieldy = group%addrs_y(l)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          sbuffer(pos) = fieldy(i,j,k)
                       end do
                    end do
                 end do
              case( MINUS_NINETY )
                 do l=1,nvector  ! loop over number of fields
                    ptr_fieldx = group%addrs_x(l)
                    do i = is, ie
                       do j = je, js, -1
                          pos = pos + 1
                          sbuffer(pos) = fieldx(i,j,k)
                       end do
                    end do
                 end do
              case( NINETY )
                 if( BTEST(group%flags_v,SCALAR_BIT) ) then
                    do l=1, nvector  ! loop over number of fields
                       ptr_fieldx = group%addrs_x(l)
                       do i = ie, is, -1
                          do j = js, je
                             pos = pos + 1
                             sbuffer(pos) = fieldx(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do l=1,nvector  ! loop over number of fields
                       ptr_fieldx = group%addrs_x(l)
                       do i = ie, is, -1
                          do j = js, je
                             pos = pos + 1
                             sbuffer(pos) = -fieldx(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case( ONE_HUNDRED_EIGHTY )
                 if( BTEST(group%flags_v,SCALAR_BIT) ) then
                    do l=1,nvector  ! loop over number of fields
                       ptr_fieldy = group%addrs_y(l)
                       do j = je, js, -1
                          do i = ie, is, -1
                             pos = pos + 1
                             sbuffer(pos) =  fieldy(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do l=1,nvector  ! loop over number of fields
                       ptr_fieldy = group%addrs_y(l)
                       do j = je, js, -1
                          do i = ie, is, -1
                             pos = pos + 1
                             sbuffer(pos) =  -fieldy(i,j,k)
                          end do
                       end do
                    end do
                 end if
              end select ! select case( rotation(n) )
           end if ! if( send(dir) )
        end do ! do n = 1, overptr%count
     endif
  enddo
  call mpp_clock_end(group_pack_clock)

  call mpp_clock_begin(group_send_clock)
  do t = 1, nsend  
     msgsize = group%send_size(t)
     if( msgsize .GT. 0 )then
        buffer_pos = group%buffer_pos_send(t)
        to_pe = group%to_pe(t)
        call mpp_send( sbuffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1)
     endif
  enddo
  call mpp_clock_end(group_send_clock)

  if(nrecv>0) then
     call mpp_clock_begin(group_wait_clock)
     call mpp_sync_self(check=EVENT_RECV)
     call mpp_clock_end(group_wait_clock)
  endif

  !---unpack the buffer
  call mpp_clock_begin(group_unpk_clock)
  !$OMP parallel do schedule(dynamic) default(shared) private(t, k, buffer_pos, msgsize, pos, m, &
  !$OMP                                                       overPtr, dir, is, ie, js, je,      &
  !$OMP                                                       ptr_field, ptr_fieldx, ptr_fieldy )
  do tk = 1, nrecv*ksize
     t = (tk-1)/ksize + 1
     k = mod((tk-1), ksize) + 1
     msgsize = group%msgsize_recv(t)
     buffer_pos = group%buffer_pos_recv(t) 
     pos = buffer_pos + (k-1)*msgsize

     m = group%recv_ind_s(t)
     if( m > 0 ) then
        overPtr => group%update_s%recv(m)
        do n = 1, overPtr%count
           dir = overPtr%dir(n)
           if( recv_s(dir) ) then
              is = overPtr%is(n); ie = overPtr%ie(n)
              js = overPtr%js(n); je = overPtr%je(n)
              do l=1,nscalar  ! loop over number of fields
                 ptr_field = group%addrs_s(l)
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       field(i,j,k) = rbuffer(pos)
                    end do
                 end do
              end do
           end if
        end do
     end if

     m = group%recv_ind_x(t)
     if( m > 0 ) then
        overPtr => group%update_x%recv(m)
        do n = 1, overPtr%count
           dir = overPtr%dir(n)
           if( recv_v(dir) ) then
              is = overPtr%is(n); ie = overPtr%ie(n)
              js = overPtr%js(n); je = overPtr%je(n)
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldx = group%addrs_x(l)
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       fieldx(i,j,k) = rbuffer(pos)
                    end do
                 end do
              end do
           endif
        enddo
     endif

     m = group%recv_ind_y(t)
     if( m > 0 ) then
        overPtr => group%update_y%recv(m)
        do n = 1, overPtr%count
           dir = overPtr%dir(n)
           if( recv_v(dir) ) then
              is = overPtr%is(n); ie = overPtr%ie(n)
              js = overPtr%js(n); je = overPtr%je(n)
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldy = group%addrs_y(l)
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       fieldy(i,j,k) = rbuffer(pos)
                    end do
                 end do
              end do
           endif
        enddo
     endif
  enddo
  call mpp_clock_end(group_unpk_clock)

  ! ---northern boundary fold
  shift = 0
  if(domain%symmetry) shift = 1
  if( nvector >0 .AND. BTEST(domain%fold,NORTH) .AND. (.NOT.BTEST(flags_v,SCALAR_BIT)) )then
     j = domain%y(1)%global%end+shift
     if( domain%y(1)%data%begin.LE.j .AND. j.LE.domain%y(1)%data%end+shift )then !fold is within domain
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
           j  = domain%y(1)%global%end+shift
           is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
           if( .NOT. domain%symmetry ) is = is - 1
           do i = is ,ie, midpoint
              if( domain%x(1)%data%begin.LE.i .AND. i.LE. domain%x(1)%data%end+shift )then
                 do l=1,nvector
                    ptr_fieldx = group%addrs_x(l)
                    ptr_fieldy = group%addrs_y(l)
                    do k = 1,ksize
                       fieldx(i,j,k) = 0.
                       fieldy(i,j,k) = 0.
                    end do
                 end do
              end if
           end do
        endif
        ! the following code code block correct an error where the data in your halo coming from
        ! other half may have the wrong sign
        !off west edge, when update north or west direction
        j = domain%y(1)%global%end+shift
        if ( recv_v(7) .OR. recv_v(5) ) then
           select case(gridtype)
           case(BGRID_NE)
              if(domain%symmetry) then
                 is = domain%x(1)%global%begin
              else
                 is = domain%x(1)%global%begin - 1
              end if
              if( is.GT.domain%x(1)%data%begin )then

                 if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-north BGRID_NE west edge ubound error.' )
                 do l=1,nvector
                    ptr_fieldx = group%addrs_x(l)
                    ptr_fieldy = group%addrs_y(l)
                    do k = 1,ksize
                       do i = domain%x(1)%data%begin,is-1
                          fieldx(i,j,k) = fieldx(2*is-i,j,k)
                          fieldy(i,j,k) = fieldy(2*is-i,j,k)
                       end do
                    end do
                 end do
              end if
           case(CGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-north CGRID_NE west edge ubound error.' )
                 do l=1,nvector
                    ptr_fieldy = group%addrs_y(l)
                    do k = 1,ksize
                       do i = domain%x(1)%data%begin,is-1
                          fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                       end do
                    end do
                 end do
              end if
           end select
        end if
        !off east edge
        is = domain%x(1)%global%end
        if(domain%x(1)%cyclic .AND. is.LT.domain%x(1)%data%end )then
           ie = domain%x(1)%data%end
           is = is + 1
           select case(gridtype)
           case(BGRID_NE)
              is = is + shift
              ie = ie + shift
              do l=1,nvector
                 ptr_fieldx = group%addrs_x(l)
                 ptr_fieldy = group%addrs_y(l)
                 do k = 1,ksize
                    do i = is,ie
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           case(CGRID_NE)
              do l=1,nvector
                 ptr_fieldy = group%addrs_y(l)
                 do k = 1,ksize
                    do i = is, ie
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,SOUTH) .OR. BTEST(domain%fold,WEST) .OR. BTEST(domain%fold,EAST) ) then
     call mpp_error(FATAL, "MPP_DO_GROUP_UPDATE: this interface does not support folded_south, " // &
          "folded_west of folded_east, contact developer")
  endif

  if(nsend>0) then
     call mpp_clock_begin(group_wait_clock)
     call mpp_sync_self( )
     call mpp_clock_end(group_wait_clock)
  endif

end subroutine MPP_DO_GROUP_UPDATE_


