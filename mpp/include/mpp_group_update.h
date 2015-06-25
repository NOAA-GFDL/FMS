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
  logical   :: recv_v(8)
  integer   :: nsend, nrecv, flags_v
  integer   :: msgsize
  integer   :: from_pe, to_pe, buffer_pos, pos
  integer   :: ksize, is, ie, js, je
  integer   :: n, l, m, i, j, k, buffer_start_pos, nk
  integer   :: shift, gridtype, midpoint
  integer   :: npack, nunpack, rotation
  character(len=8)            :: text

  MPP_TYPE_ :: buffer(mpp_domains_stack_size)
  MPP_TYPE_ :: field (group%is_s:group%ie_s,group%js_s:group%je_s, group%ksize_s)
  MPP_TYPE_ :: fieldx(group%is_x:group%ie_x,group%js_x:group%je_x, group%ksize_v)
  MPP_TYPE_ :: fieldy(group%is_y:group%ie_y,group%js_y:group%je_y, group%ksize_v)
  pointer(ptr, buffer )
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
  if(nvector > 0) recv_v = group%recv_v

  ptr = LOC(mpp_domains_stack)

  !--- set reset_index_s and reset_index_v to 0 
  group%reset_index_s = 0
  group%reset_index_v = 0

  if(.not. group%initialized) call set_group_update(group,domain)

  nrecv = group%nrecv
  nsend = group%nsend

  !---pre-post receive.
  call mpp_clock_begin(group_recv_clock)
  do m = 1, nrecv
     msgsize = group%recv_size(m)
     from_pe = group%from_pe(m)
     if( msgsize .GT. 0 )then
        buffer_pos = group%buffer_pos_recv(m)
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., &
             tag=COMM_TAG_1)
     end if
  end do

 !pack the data
  call mpp_clock_end(group_recv_clock)

  flags_v = group%flags_v
  npack = group%npack

  call mpp_clock_begin(group_pack_clock)
  !pack the data
  buffer_start_pos = 0
#include <group_update_pack.inc>
  call mpp_clock_end(group_pack_clock)

  call mpp_clock_begin(group_send_clock)
  do n = 1, nsend  
     msgsize = group%send_size(n)
     if( msgsize .GT. 0 )then
        buffer_pos = group%buffer_pos_send(n)
        to_pe = group%to_pe(n)
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1)
     endif
  enddo
  call mpp_clock_end(group_send_clock)

  if(nrecv>0) then
     call mpp_clock_begin(group_wait_clock)
     call mpp_sync_self(check=EVENT_RECV)
     call mpp_clock_end(group_wait_clock)
  endif

  !---unpack the buffer
  nunpack = group%nunpack
  call mpp_clock_begin(group_unpk_clock)
#include <group_update_unpack.inc>
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


subroutine MPP_START_GROUP_UPDATE_(group, domain, d_type, reuse_buffer)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D),              intent(inout) :: domain  
  MPP_TYPE_,                   intent(in)    :: d_type
  logical,  optional,          intent(in)    :: reuse_buffer

  integer   :: nscalar, nvector
  integer   :: nsend, nrecv, flags_v
  integer   :: msgsize, npack, rotation
  integer   :: from_pe, to_pe, buffer_pos, pos
  integer   :: ksize, is, ie, js, je
  integer   :: n, l, m, i, j, k, buffer_start_pos, nk
  logical   :: reuse_buf_pos
  character(len=8)            :: text

  MPP_TYPE_ :: buffer(size(mpp_domains_stack_nonblock(:)))
  MPP_TYPE_ :: field (group%is_s:group%ie_s,group%js_s:group%je_s, group%ksize_s)
  MPP_TYPE_ :: fieldx(group%is_x:group%ie_x,group%js_x:group%je_x, group%ksize_v)
  MPP_TYPE_ :: fieldy(group%is_y:group%ie_y,group%js_y:group%je_y, group%ksize_v)
  pointer( ptr, buffer )
  pointer(ptr_field, field)
  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)

  nscalar = group%nscalar
  nvector = group%nvector

  if(nscalar>0) then
     ksize = group%ksize_s
  else
     ksize = group%ksize_v
  endif

  !--- set reset_index_s and reset_index_v to 0 
  group%reset_index_s = 0
  group%reset_index_v = 0

  reuse_buf_pos = .FALSE.
  if (PRESENT(reuse_buffer)) reuse_buf_pos = reuse_buffer

  if (.not. group%initialized) then
    call set_group_update(group,domain)
  endif

  if (.not. reuse_buf_pos) then
     group%buffer_start_pos = nonblock_group_buffer_pos
     nonblock_group_buffer_pos = nonblock_group_buffer_pos + group%tot_msgsize
     mpp_domains_stack_hwm = nonblock_group_buffer_pos + 1    
     if( mpp_domains_stack_hwm .GT. mpp_domains_stack_size )then
        write( text,'(i8)' )mpp_domains_stack_hwm
        call mpp_error( FATAL, 'set_group_update: mpp_domains_stack overflow, '// &
                     'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
     end if

  else if( group%buffer_start_pos < 0 ) then
     call mpp_error(FATAL, "MPP_START_GROUP_UPDATE: group%buffer_start_pos is not set")
  endif

  nrecv = group%nrecv
  nsend = group%nsend

  ptr = LOC(mpp_domains_stack_nonblock)

  ! Make sure it is not in the middle of the old version of non-blocking halo update.
  if(num_update>0) call mpp_error(FATAL, "MPP_START_GROUP_UPDATE: can not be called in the middle of "// &
              "mpp_start_update_domains/mpp_complete_update_domains call")

  num_nonblock_group_update = num_nonblock_group_update + 1

  !---pre-post receive.
  call mpp_clock_begin(nonblock_group_recv_clock)
  do m = 1, nrecv
     msgsize = group%recv_size(m)
     from_pe = group%from_pe(m)
     if( msgsize .GT. 0 )then
        buffer_pos = group%buffer_pos_recv(m) + group%buffer_start_pos
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., &
             tag=COMM_TAG_1, request=group%request_recv(m))
#ifdef use_libMPI
        group%type_recv(m) = MPI_TYPE_
#endif
     end if
  end do
  call mpp_clock_end(nonblock_group_recv_clock)

  flags_v = group%flags_v

  !pack the data
  call mpp_clock_begin(nonblock_group_pack_clock)
  npack = group%npack
  buffer_start_pos = group%buffer_start_pos
#include <group_update_pack.inc>
  call mpp_clock_end(nonblock_group_pack_clock)

  call mpp_clock_begin(nonblock_group_send_clock)
  do n = 1, nsend  
     msgsize = group%send_size(n)
     if( msgsize .GT. 0 )then
        buffer_pos = group%buffer_pos_send(n) + group%buffer_start_pos
        to_pe = group%to_pe(n)
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1, &
                       request=group%request_send(n))
     endif
  enddo
  call mpp_clock_end(nonblock_group_send_clock)

end subroutine MPP_START_GROUP_UPDATE_

subroutine MPP_COMPLETE_GROUP_UPDATE_(group, domain, d_type)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D),              intent(inout) :: domain  
  MPP_TYPE_,                   intent(in)    :: d_type

  integer   :: nsend, nrecv, nscalar, nvector
  integer   :: k, buffer_pos, msgsize, pos, m, n, l
  integer   :: is, ie, js, je, dir, ksize, i, j
  integer   :: shift, gridtype, midpoint, flags_v
  integer   :: nunpack, rotation, buffer_start_pos, nk
  logical   :: recv_v(8)
  MPP_TYPE_ :: buffer(size(mpp_domains_stack_nonblock(:)))
  MPP_TYPE_ :: field (group%is_s:group%ie_s,group%js_s:group%je_s, group%ksize_s)
  MPP_TYPE_ :: fieldx(group%is_x:group%ie_x,group%js_x:group%je_x, group%ksize_v)
  MPP_TYPE_ :: fieldy(group%is_y:group%ie_y,group%js_y:group%je_y, group%ksize_v)
  pointer(ptr, buffer )
  pointer(ptr_field, field)
  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)

  gridtype = group%gridtype
  flags_v = group%flags_v
  nscalar = group%nscalar
  nvector = group%nvector
  nrecv = group%nrecv
  nsend = group%nsend
  if(nscalar>0) then
     ksize = group%ksize_s
  else
     ksize = group%ksize_v
  endif
  if(nvector > 0) recv_v = group%recv_v
  ptr = LOC(mpp_domains_stack_nonblock)

  if(num_nonblock_group_update < 1) call mpp_error(FATAL, &
    'mpp_start_group_update must be called before calling mpp_end_group_update')  
  num_nonblock_group_update = num_nonblock_group_update - 1
  complete_group_update_on = .true.

  if(nrecv>0) then
     call mpp_clock_begin(nonblock_group_wait_clock)
     call mpp_sync_self(check=EVENT_RECV, request=group%request_recv(1:nrecv), &
                        msg_size=group%recv_size(1:nrecv), msg_type=group%type_recv(1:nrecv))
     call mpp_clock_end(nonblock_group_wait_clock)
  endif

  !---unpack the buffer
  nunpack = group%nunpack

  call mpp_clock_begin(nonblock_group_unpk_clock)
  buffer_start_pos = group%buffer_start_pos
#include <group_update_unpack.inc>
  call mpp_clock_end(nonblock_group_unpk_clock)

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
     call mpp_error(FATAL, "MPP_COMPLETE_GROUP_UPDATE: this interface does not support folded_south, " // &
          "folded_west of folded_east, contact developer")
  endif

  if(nsend>0) then
     call mpp_clock_begin(nonblock_group_wait_clock)
     call mpp_sync_self(check=EVENT_SEND, request=group%request_send(1:nsend) )
     call mpp_clock_end(nonblock_group_wait_clock)
  endif

  if( num_nonblock_group_update == 0) then
     nonblock_group_buffer_pos   = 0
  endif

end subroutine MPP_COMPLETE_GROUP_UPDATE_

subroutine MPP_RESET_GROUP_UPDATE_FIELD_2D_(group, field)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,                   intent(in)    :: field(:,:)

  group%reset_index_s = group%reset_index_s + 1

  if(group%reset_index_s > group%nscalar) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_2D_: group%reset_index_s > group%nscalar")
  if(size(field,1) .NE. group%isize_s .OR. size(field,2) .NE. group%jsize_s .OR. group%ksize_s .NE. 1) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_2D_: size of field does not match the size stored in group") 

  group%addrs_s(group%reset_index_s) = LOC(field)

end subroutine MPP_RESET_GROUP_UPDATE_FIELD_2D_

subroutine MPP_RESET_GROUP_UPDATE_FIELD_3D_(group, field)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,                   intent(in)    :: field(:,:,:)

  group%reset_index_s = group%reset_index_s + 1

  if(group%reset_index_s > group%nscalar) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_3D_: group%reset_index_s > group%nscalar")
  if(size(field,1) .NE. group%isize_s .OR. size(field,2) .NE. group%jsize_s .OR. size(field,3) .NE. group%ksize_s) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_3D_: size of field does not match the size stored in group") 

  group%addrs_s(group%reset_index_s) = LOC(field)

end subroutine MPP_RESET_GROUP_UPDATE_FIELD_3D_

subroutine MPP_RESET_GROUP_UPDATE_FIELD_4D_(group, field)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,                   intent(in)    :: field(:,:,:,:)

  group%reset_index_s = group%reset_index_s + 1

  if(group%reset_index_s > group%nscalar) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_4D_: group%reset_index_s > group%nscalar")
  if(size(field,1) .NE. group%isize_s .OR. size(field,2) .NE. group%jsize_s .OR. &
              size(field,3)*size(field,4) .NE. group%ksize_s) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_4D_: size of field does not match the size stored in group") 

  group%addrs_s(group%reset_index_s) = LOC(field)

end subroutine MPP_RESET_GROUP_UPDATE_FIELD_4D_


subroutine MPP_RESET_GROUP_UPDATE_FIELD_2D_V_(group, fieldx, fieldy)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,                   intent(in)    :: fieldx(:,:), fieldy(:,:)
  integer :: indx

  group%reset_index_v = group%reset_index_v + 1

  if(group%reset_index_v > group%nvector) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_2D_V_: group%reset_index_v > group%nvector")
  if(size(fieldx,1) .NE. group%isize_x .OR. size(fieldx,2) .NE. group%jsize_x .OR. group%ksize_v .NE. 1) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_2D_V_: size of fieldx does not match the size stored in group") 
  if(size(fieldy,1) .NE. group%isize_y .OR. size(fieldy,2) .NE. group%jsize_y ) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_2D_V_: size of fieldy does not match the size stored in group") 

  group%addrs_x(group%reset_index_v) = LOC(fieldx)
  group%addrs_y(group%reset_index_v) = LOC(fieldy)

end subroutine MPP_RESET_GROUP_UPDATE_FIELD_2D_V_


subroutine MPP_RESET_GROUP_UPDATE_FIELD_3D_V_(group, fieldx, fieldy)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,                   intent(in)    :: fieldx(:,:,:), fieldy(:,:,:)
  integer :: indx

  group%reset_index_v = group%reset_index_v + 1

  if(group%reset_index_v > group%nvector) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_3D_V_: group%reset_index_v > group%nvector")
  if(size(fieldx,1) .NE. group%isize_x .OR. size(fieldx,2) .NE. group%jsize_x .OR. size(fieldx,3) .NE. group%ksize_v) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_3D_V_: size of fieldx does not match the size stored in group") 
  if(size(fieldy,1) .NE. group%isize_y .OR. size(fieldy,2) .NE. group%jsize_y .OR. size(fieldy,3) .NE. group%ksize_v) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_3D_V_: size of fieldy does not match the size stored in group") 

  group%addrs_x(group%reset_index_v) = LOC(fieldx)
  group%addrs_y(group%reset_index_v) = LOC(fieldy)

end subroutine MPP_RESET_GROUP_UPDATE_FIELD_3D_V_


subroutine MPP_RESET_GROUP_UPDATE_FIELD_4D_V_(group, fieldx, fieldy)
  type(mpp_group_update_type), intent(inout) :: group
  MPP_TYPE_,                   intent(in)    :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  integer :: indx

  group%reset_index_v = group%reset_index_v + 1

  if(group%reset_index_v > group%nvector) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_4D_V_: group%reset_index_v > group%nvector")
  if(size(fieldx,1) .NE. group%isize_x .OR. size(fieldx,2) .NE. group%jsize_x .OR. &
              size(fieldx,3)*size(fieldx,4) .NE. group%ksize_v) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_4D_V_: size of fieldx does not match the size stored in group") 
  if(size(fieldy,1) .NE. group%isize_y .OR. size(fieldy,2) .NE. group%jsize_y .OR. &
              size(fieldy,3)*size(fieldy,4) .NE. group%ksize_v) &
     call mpp_error(FATAL, "MPP_RESET_GROUP_UPDATE_FIELD_4D_V_: size of fieldy does not match the size stored in group") 

  group%addrs_x(group%reset_index_v) = LOC(fieldx)
  group%addrs_y(group%reset_index_v) = LOC(fieldy)

end subroutine MPP_RESET_GROUP_UPDATE_FIELD_4D_V_


