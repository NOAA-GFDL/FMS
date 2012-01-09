function MPP_START_UPDATE_DOMAINS_2D_( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete)
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer,          intent(in), optional :: update_id
  logical,          intent(in), optional :: complete
  integer                                :: MPP_START_UPDATE_DOMAINS_2D_

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
  pointer( ptr, field3D )
  ptr = LOC(field)

  MPP_START_UPDATE_DOMAINS_2D_ = mpp_start_update_domains(field3D, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete)
  return

end function MPP_START_UPDATE_DOMAINS_2D_

function MPP_START_UPDATE_DOMAINS_3D_( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )

  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(domain%x(1)%data%begin:,domain%y(1)%data%begin:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer,          intent(in), optional :: update_id
  logical,          intent(in), optional :: complete
  integer                                :: MPP_START_UPDATE_DOMAINS_3D_

  !--- local variables
  integer                    :: current_id, ke_max
  integer                    :: update_whalo, update_ehalo, update_shalo, update_nhalo, update_flags, update_position
  integer                    :: tile, max_ntile, ntile
  logical                    :: set_mismatch, is_complete
  logical                    :: do_update, reuse_id_update
  integer, save              :: isize=0, jsize=0, l_size=0, list=0
  integer, save              :: pos, whalosz, ehalosz, shalosz, nhalosz, update_flags_saved
  character(len=128)         :: text, field_name
  integer, save              :: ke_list(MAX_DOMAIN_FIELDS, MAX_TILES)=0
  integer(LONG_KIND), save   :: f_addrs(MAX_DOMAIN_FIELDS, MAX_TILES)=-9999
  type(overlapSpec), pointer :: update => NULL()  
  MPP_TYPE_                  :: d_type

  field_name = "unknown"
  if(present(name)) field_name = name

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

  update_flags = XUPDATE+YUPDATE   !default
  if( PRESENT(flags) )update_flags = flags

  update_position = CENTER
  if(present(position)) then
     !--- when there is NINETY or MINUS_NINETY rotation for some contact, the salar data can not be on E or N-cell,
     if(domain%rotated_ninety .AND. ( position == EAST .OR. position == NORTH ) )  &
          call mpp_error(FATAL, 'MPP_START_UPDATE_DOMAINS_3D: hen there is NINETY or MINUS_NINETY rotation, ' // &
          'can not use scalar version update_domain for data on E or N-cell' )
     update_position = position    
  endif

  max_ntile = domain%max_ntile_pe
  ntile = size(domain%x(:))
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if
  tile = 1

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_3D: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_UPDATE_3D: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete

  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  f_addrs(list,tile) = LOC(field)
  ke_list(list,tile) = size(field,3)

  if(list == 1 .AND. tile == 1 )then
     isize=size(field,1); jsize=size(field,2); pos = update_position
     whalosz = update_whalo; ehalosz = update_ehalo; shalosz = update_shalo; nhalosz = update_nhalo
     update_flags_saved = update_flags
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize /= size(field,1))
     set_mismatch = set_mismatch .OR. (jsize /= size(field,2))
     set_mismatch = set_mismatch .OR. (update_position /= pos)
     set_mismatch = set_mismatch .OR. (update_whalo /= whalosz)
     set_mismatch = set_mismatch .OR. (update_ehalo /= ehalosz)
     set_mismatch = set_mismatch .OR. (update_shalo /= shalosz)
     set_mismatch = set_mismatch .OR. (update_nhalo /= nhalosz)
     set_mismatch = set_mismatch .OR. (update_flags_saved /= update_flags)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS: Incompatible field at count '//text//' for group update.' )
     endif
  endif

  if(is_complete) then
     l_size = list
     list = 0
  end if
  
  if(do_update) then
     num_update = num_update + 1

     if( PRESENT(update_id) ) then
        if( update_id < 1 .OR. update_id > MAX_NONBLOCK_UPDATE ) then
           write( text,'(a,i8,a,i8)' ) 'optional argument update_id =', update_id, &
                'is less than 1 or  greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
           call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS: '//trim(text))
        endif
        current_id = update_id
        reuse_id_update = .true.
        !--- when reuse the update_id, make sure update_flag, halo size and update_position are still the same
        if( update_flags_list(current_id) .NE. update_flags .OR. &
             update_whalo_list(current_id) .NE. update_whalo .OR. &
             update_ehalo_list(current_id) .NE. update_ehalo .OR. &
             update_shalo_list(current_id) .NE. update_shalo .OR. &
             update_nhalo_list(current_id) .NE. update_nhalo .OR. &
             update_position_list(current_id) .NE. update_position ) then
           call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS: mismatch for optional argument for field '//trim(field_name) )
        endif
     else
        reuse_id_update = .false.
        current_id_update = current_id_update + 1
        if( current_id_update > MAX_NONBLOCK_UPDATE ) then
           write( text,'(a,i8,a,i8)' ) 'num_fields =', current_id_update, ' greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
           call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS: '//trim(text))
        endif
        current_id = current_id_update
        update_flags_list(current_id) = update_flags
        update_whalo_list(current_id) = update_whalo
        update_ehalo_list(current_id) = update_ehalo
        update_shalo_list(current_id) = update_shalo
        update_nhalo_list(current_id) = update_nhalo
        update_position_list(current_id) = update_position
        recv_pos_list(current_id) = nonblock_buffer_pos  
     endif
     MPP_START_UPDATE_DOMAINS_3D_ = current_id

     ke_max = maxval(ke_list(1:l_size,1:ntile))
     if( domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) )then
        update => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, update_position)
        call mpp_start_do_update(current_id, f_addrs(1:l_size,1:ntile), domain, update, d_type, &
                                 ke_max, ke_list(1:l_size,1:ntile), update_flags, reuse_id_update, field_name )
     endif
     l_size=0; f_addrs=-9999; isize=0;  jsize=0;  ke_list=0
  else
     if(present(update_id)) then
        MPP_START_UPDATE_DOMAINS_3D_ = update_id
     else
        MPP_START_UPDATE_DOMAINS_3D_ = 0
     endif
  endif


end function MPP_START_UPDATE_DOMAINS_3D_

!##########################################################################################
function MPP_START_UPDATE_DOMAINS_4D_( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer,          intent(in), optional :: update_id
  logical,          intent(in), optional :: complete
  integer                                :: MPP_START_UPDATE_DOMAINS_4D_

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
  pointer( ptr, field3D )
  ptr = LOC(field)

  MPP_START_UPDATE_DOMAINS_4D_ = mpp_start_update_domains(field3D, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete)
  return

end function MPP_START_UPDATE_DOMAINS_4D_

!##########################################################################################
function MPP_START_UPDATE_DOMAINS_5D_( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete)
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer,          intent(in), optional :: update_id
  logical,          intent(in), optional :: complete
  integer                                :: MPP_START_UPDATE_DOMAINS_5D_

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
  pointer( ptr, field3D )
  ptr = LOC(field)

  MPP_START_UPDATE_DOMAINS_5D_ = mpp_start_update_domains(field3D, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )
  return

end function MPP_START_UPDATE_DOMAINS_5D_

!##################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_2D_( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, buffer, complete )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  MPP_TYPE_,     intent(inout), optional :: buffer(:)

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
  pointer( ptr, field3D )
  ptr = LOC(field)
  call mpp_complete_update_domains(id_update, field3D, domain, flags, position, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count, buffer, complete )

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_2D_

!##################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_3D_( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, buffer, complete )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(domain%x(1)%data%begin:,domain%y(1)%data%begin:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  MPP_TYPE_,     intent(inout), optional :: buffer(:)


  integer                    :: update_whalo, update_ehalo, update_shalo, update_nhalo
  integer                    :: update_position, update_flags
  type(overlapSpec), pointer :: update => NULL()
  integer                    :: tile, max_ntile, ntile
  logical                    :: is_complete
  logical                    :: do_update
  integer                    :: ke_max, buffer_size
  integer, save              :: list=0, bsize=0, l_size=0
  integer, save              :: ke_list(MAX_DOMAIN_FIELDS, MAX_TILES)=0
  integer(LONG_KIND), save   :: f_addrs(MAX_DOMAIN_FIELDS, MAX_TILES)=-9999
  integer(LONG_KIND), save   :: b_addrs(MAX_DOMAIN_FIELDS, MAX_TILES)=-9999
  character(len=128)         :: text

      MPP_TYPE_        :: d_type

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

  update_position = CENTER
  if(present(position)) update_position = position  
  update_flags = XUPDATE+YUPDATE   !default
  if( PRESENT(flags) )update_flags = flags

  max_ntile = domain%max_ntile_pe
  ntile = size(domain%x(:))
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if
  tile = 1

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_COMPLETE_UPDATE_DOMAINS_3D: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_UPDATE_3D: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if
  do_update = (tile == ntile) .AND. is_complete
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_COMPLETE_UPDATE_DOMAINS_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  f_addrs(list, tile) = LOC(field)
  buffer_size = 0
  if(present(buffer)) then
     buffer_size = size(buffer(:))
     b_addrs(list, tile) = LOC(buffer)
  end if

  ke_list(list,tile) = size(field,3)
  bsize = max(bsize, buffer_size)

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

  if(is_complete) then
     l_size = list
     list = 0
  end if

  if(do_update) then
     num_update = num_update - 1
     if( domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) ) then
        update => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, update_position)
        ke_max = maxval(ke_list(1:l_size,1:ntile))
        call mpp_complete_do_update(id_update, f_addrs(1:l_size,1:ntile), domain, update, d_type, &
                                    ke_max, ke_list(1:l_size,1:ntile), b_addrs(1:l_size,1:ntile), bsize, update_flags) 
     endif
     l_size=0; f_addrs=-9999; bsize=0; b_addrs=-9999; ke_list=0

     !--- For the last call of mpp_complete_update_domains
     !--- reset everything to init state
     if( num_update == 0) then
        recv_pos_list = 0
        send_pos_list = 0
        recv_msgsize_list = 0
        send_msgsize_list = 0
        update_flags_list = 0
        update_position_list = 0
        update_whalo_list = 0
        update_ehalo_list = 0
        update_shalo_list = 0
        update_nhalo_list = 0
        current_id_update = 0
        nonblock_buffer_pos   = 0
     endif
  endif

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_3D_

!##################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_4D_( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, buffer, complete )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  MPP_TYPE_,     intent(inout), optional :: buffer(:)

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
  pointer( ptr, field3D )
  ptr = LOC(field)
  call mpp_complete_update_domains(id_update, field3D, domain, flags, position, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count, buffer, complete )

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_4D_

!##################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_5D_( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, buffer, complete )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain  
  MPP_TYPE_,        intent(inout)        :: field(:,:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  MPP_TYPE_,     intent(inout), optional :: buffer(:)

  MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
  pointer( ptr, field3D )
  ptr = LOC(field)
  call mpp_complete_update_domains(id_update, field3D, domain, flags, position, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count, buffer, complete )

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_5D_

#ifdef  VECTOR_FIELD_
function MPP_START_UPDATE_DOMAINS_2D_V_( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )
  !updates data domain of 3D field whose computational domains have been computed
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:), fieldy(:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer,          intent(in), optional :: update_id
  logical,          intent(in), optional :: complete
  integer                                :: MPP_START_UPDATE_DOMAINS_2D_V_
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  MPP_START_UPDATE_DOMAINS_2D_V_ = mpp_start_update_domains(field3Dx, field3Dy, domain, flags, gridtype, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )

  return

end function MPP_START_UPDATE_DOMAINS_2D_V_

!###################################################################################
function MPP_START_UPDATE_DOMAINS_3D_V_( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )
  !updates data domain of 3D field whose computational domains have been computed
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:), fieldy(:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer,          intent(in), optional :: update_id
  logical,          intent(in), optional :: complete
  !--- local variables
  integer                     :: MPP_START_UPDATE_DOMAINS_3D_V_
  integer                     :: update_whalo, update_ehalo, update_shalo, update_nhalo
  integer                     :: grid_offset_type, position_x, position_y, update_flags, current_id
  logical                     :: do_update, is_complete, set_mismatch
  integer                     :: ntile, max_ntile, tile, ke_max
  logical                     :: exchange_uv, reuse_id_update
  character(len=128)          :: text, field_name
  integer,            save    :: whalosz, ehalosz, shalosz, nhalosz
  integer,            save    :: isize(2)=0,jsize(2)=0,l_size=0, offset_type=0, list=0
  integer,            save    :: ke_list (MAX_DOMAIN_FIELDS, MAX_TILES)=0
  integer(LONG_KIND), save    :: f_addrsx(MAX_DOMAIN_FIELDS, MAX_TILES)=-9999
  integer(LONG_KIND), save    :: f_addrsy(MAX_DOMAIN_FIELDS, MAX_TILES)=-9999
  type(overlapSpec),  pointer :: updatex => NULL()
  type(overlapSpec),  pointer :: updatey => NULL()
  MPP_TYPE_                   :: d_type

  field_name = "unknown"
  if(present(name)) field_name = name

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

  max_ntile = domain%max_ntile_pe
  ntile = size(domain%x(:))

  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if
  tile = 1

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_UPDATE_3D_V: "// &
          "optional argument tile_count should be present when number of tiles on some pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif

  f_addrsx(list, tile) = LOC(fieldx)
  f_addrsy(list, tile) = LOC(fieldy)

  ke_list(list, tile) = size(fieldx,3)

  if(list == 1 .AND. tile == 1)then
     isize(1)=size(fieldx,1); jsize(1)=size(fieldx,2)
     isize(2)=size(fieldy,1); jsize(2)=size(fieldy,2)
     offset_type = grid_offset_type
     whalosz = update_whalo; ehalosz = update_ehalo; shalosz = update_shalo; nhalosz = update_nhalo
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize(1) /= size(fieldx,1))
     set_mismatch = set_mismatch .OR. (jsize(1) /= size(fieldx,2))
     set_mismatch = set_mismatch .OR. (isize(2) /= size(fieldy,1))
     set_mismatch = set_mismatch .OR. (jsize(2) /= size(fieldy,2))
     set_mismatch = set_mismatch .OR. (grid_offset_type /= offset_type)
     set_mismatch = set_mismatch .OR. (update_whalo /= whalosz)
     set_mismatch = set_mismatch .OR. (update_ehalo /= ehalosz)
     set_mismatch = set_mismatch .OR. (update_shalo /= shalosz)
     set_mismatch = set_mismatch .OR. (update_nhalo /= nhalosz)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: Incompatible field at count '//text//' for group vector update.' )
     end if
  end if
  if(is_complete) then
     l_size = list
     list = 0
  end if
  if(do_update)then
     num_update = num_update + 1
     if( PRESENT(update_id) ) then
        reuse_id_update = .true.
        if( update_id < 1 .OR. update_id > MAX_NONBLOCK_UPDATE ) then
           write( text,'(a,i8,a,i8)' ) 'optional argument update_id =', update_id, &
                'is less than 1 or  greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
           call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: '//trim(text))
        endif
        current_id = update_id
        !--- when reuse the update_id, make sure update_flag, halo size and update_position are still the same
        if( update_flags_list(current_id) .NE. update_flags .OR. &
             update_whalo_list(current_id) .NE. update_whalo .OR. &
             update_ehalo_list(current_id) .NE. update_ehalo .OR. &
             update_shalo_list(current_id) .NE. update_shalo .OR. &
             update_nhalo_list(current_id) .NE. update_nhalo .OR. &
             update_gridtype_list(current_id) .NE. grid_offset_type ) then
           call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: mismatch for optional argument for field '//trim(field_name) )
        endif
     else
        reuse_id_update = .false.
        current_id_update = current_id_update + 1
        current_id = current_id_update
        if( current_id_update > MAX_NONBLOCK_UPDATE ) then
           write( text,'(a,i8,a,i8)' ) 'num_fields =', current_id_update, ' greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
           call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: '//trim(text))
        endif
        update_flags_list(current_id) = update_flags
        update_whalo_list(current_id) = update_whalo
        update_ehalo_list(current_id) = update_ehalo
        update_shalo_list(current_id) = update_shalo
        update_nhalo_list(current_id) = update_nhalo
        update_gridtype_list(current_id) = grid_offset_type
        recv_pos_list(current_id) = nonblock_buffer_pos  
     endif
     MPP_START_UPDATE_DOMAINS_3D_V_ = current_id
     if( domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) )then
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

        ke_max = maxval(ke_list(1:l_size,1:ntile))
        ke_max = maxval(ke_list(1:l_size,1:ntile))
        if(exchange_uv) then
           call mpp_start_do_update(current_id, f_addrsx(1:l_size,1:ntile), f_addrsy(1:l_size,1:ntile), domain, &
                                    updatey, updatex, d_type, ke_max, ke_list(1:l_size,1:ntile), grid_offset_type, &
                                    update_flags, reuse_id_update, field_name)
        else
           call mpp_start_do_update(current_id, f_addrsx(1:l_size,1:ntile), f_addrsy(1:l_size,1:ntile), domain, &
                                    updatex, updatey, d_type, ke_max, ke_list(1:l_size,1:ntile), gridtype,    & 
                                    update_flags, reuse_id_update, field_name)
        endif
     endif
     l_size=0; f_addrsx=-9999; f_addrsy=-9999; isize=0;  jsize=0;  ke_list=0
  else
     if(present(update_id)) then
        MPP_START_UPDATE_DOMAINS_3D_V_ = update_id
     else
        MPP_START_UPDATE_DOMAINS_3D_V_ = 0
     endif
  end if

  return

end function MPP_START_UPDATE_DOMAINS_3D_V_

function MPP_START_UPDATE_DOMAINS_4D_V_( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )
  !updates data domain of 3D field whose computational domains have been computed
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer,          intent(in), optional :: update_id
  logical,          intent(in), optional :: complete
  integer                                :: MPP_START_UPDATE_DOMAINS_4D_V_
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  MPP_START_UPDATE_DOMAINS_4D_V_ = mpp_start_update_domains(field3Dx, field3Dy, domain, flags, gridtype, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )

  return

end function MPP_START_UPDATE_DOMAINS_4D_V_

function MPP_START_UPDATE_DOMAINS_5D_V_( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )
  !updates data domain of 3D field whose computational domains have been computed
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:,:), fieldy(:,:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer,          intent(in), optional :: update_id
  logical,          intent(in), optional :: complete
  integer                                :: MPP_START_UPDATE_DOMAINS_5D_V_
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4)*size(fieldx,5))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4)*size(fieldy,5))
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  MPP_START_UPDATE_DOMAINS_5D_V_ = mpp_start_update_domains(field3Dx, field3Dy, domain, flags, gridtype, &
                                   whalo, ehalo, shalo, nhalo, name, tile_count, update_id, complete )

  return

end function MPP_START_UPDATE_DOMAINS_5D_V_

!####################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_2D_V_( id_update, fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery, complete )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:), fieldy(:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete

  MPP_TYPE_,     intent(inout), optional :: bufferx(:), buffery(:)
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  call mpp_complete_update_domains(id_update, field3Dx, field3Dy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery, complete )

  return

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_2D_V_

!####################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_3D_V_( id_update, fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery, complete )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:), fieldy(:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  MPP_TYPE_,     intent(inout), optional :: bufferx(:), buffery(:)
  logical,          intent(in), optional :: complete

  integer                     :: update_whalo, update_ehalo, update_shalo, update_nhalo
  integer                     :: grid_offset_type, position_x, position_y, update_flags
  logical                     :: do_update, is_complete
  integer                     :: ntile, max_ntile, tile, ke_max
  integer                     :: bufferx_size, buffery_size
  logical                     :: exchange_uv
  character(len=128)          :: text
  integer,            save    :: l_size=0, list=0, bsizex=1, bsizey=1
  integer,            save    :: ke_list (MAX_DOMAIN_FIELDS, MAX_TILES)=0
  integer(LONG_KIND), save    :: f_addrsx(MAX_DOMAIN_FIELDS, MAX_TILES)=-9999
  integer(LONG_KIND), save    :: f_addrsy(MAX_DOMAIN_FIELDS, MAX_TILES)=-9999
  integer(LONG_KIND) ,save    :: b_addrsx(MAX_DOMAIN_FIELDS, MAX_TILES)=-9999
  integer(LONG_KIND) ,save    :: b_addrsy(MAX_DOMAIN_FIELDS, MAX_TILES)=-9999
  type(overlapSpec),  pointer :: updatex => NULL()
  type(overlapSpec),  pointer :: updatey => NULL()
  MPP_TYPE_                   :: d_type

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

  max_ntile = domain%max_ntile_pe
  ntile = size(domain%x(:))

  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if
  tile = 1

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_UPDATE_3D_V: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_UPDATE_3D_V: "// &
          "optional argument tile_count should be present when number of tiles on some pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_UPDATE_3D_V: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif

  f_addrsx(list, tile) = LOC(fieldx)
  f_addrsy(list, tile) = LOC(fieldy)
  bufferx_size = 1; buffery_size = 1
  if(present(bufferx)) then
     bufferx_size = size(bufferx(:))
     buffery_size = size(buffery(:))
     b_addrsx(list, tile) = LOC(bufferx)
     b_addrsy(list, tile) = LOC(buffery)
  end if

  ke_list(list, tile) = size(fieldx,3)
  bsizex = max(bsizex, bufferx_size)
  bsizey = max(bsizey, buffery_size)

  if(is_complete) then
     l_size = list
     list = 0
  end if
  if(do_update)then
     num_update = num_update - 1
     if( domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) )then
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

        ke_max = maxval(ke_list(1:l_size,1:ntile))
        if(exchange_uv) then
           call mpp_complete_do_update(id_update, f_addrsx(1:l_size,1:ntile), f_addrsy(1:l_size,1:ntile), domain, &
                                       updatey, updatex, d_type, ke_max, ke_list(1:l_size,1:ntile),               &
                                       b_addrsx(1:l_size,1:ntile), b_addrsy(1:l_size,1:ntile), bsizex, bsizey,    &
                                       grid_offset_type, update_flags)
        else
           call mpp_complete_do_update(id_update, f_addrsx(1:l_size,1:ntile), f_addrsy(1:l_size,1:ntile), domain, &
                                       updatex, updatey, d_type, ke_max, ke_list(1:l_size,1:ntile),               &
                                       b_addrsx(1:l_size,1:ntile), b_addrsy(1:l_size,1:ntile), bsizex, bsizey,    &
                                       grid_offset_type, update_flags)
        endif
     endif
     l_size=0; f_addrsx=-9999; f_addrsy=-9999; ke_list=0
     bsizex=1; b_addrsx=-9999; bsizey=1; b_addrsy=-9999
     !--- For the last call of mpp_complete_update_domains
     !--- reset everything to init state
     if( num_update == 0) then
        recv_pos_list = 0
        send_pos_list = 0
        recv_msgsize_list = 0
        send_msgsize_list = 0
        update_flags_list = 0
        update_gridtype_list = 0
        update_whalo_list = 0
        update_ehalo_list = 0
        update_shalo_list = 0
        update_nhalo_list = 0
        current_id_update = 0
        nonblock_buffer_pos   = 0
     endif
  end if


end subroutine MPP_COMPLETE_UPDATE_DOMAINS_3D_V_

!####################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_4D_V_( id_update, fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery, complete )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete

  MPP_TYPE_,     intent(inout), optional :: bufferx(:), buffery(:)
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  call mpp_complete_update_domains(id_update, field3Dx, field3Dy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery, complete )

  return

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_4D_V_

!####################################################################################
subroutine MPP_COMPLETE_UPDATE_DOMAINS_5D_V_( id_update, fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery, complete )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:,:), fieldy(:,:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete

  MPP_TYPE_,     intent(inout), optional :: bufferx(:), buffery(:)
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4)*size(fieldx,5))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4)*size(fieldy,5))
  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx = LOC(fieldx)
  ptry = LOC(fieldy)

  call mpp_complete_update_domains(id_update, field3Dx, field3Dy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, bufferx, buffery, complete )

  return

end subroutine MPP_COMPLETE_UPDATE_DOMAINS_5D_V_

#endif
