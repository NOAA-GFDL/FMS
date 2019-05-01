! -*-f90-*- 

!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

! this routine is used to retrieve scalar boundary data for symmetric domain. 

subroutine MPP_GET_BOUNDARY_AD_2D_(field, domain, ebuffer, sbuffer, wbuffer, nbuffer, flags, &
                                position, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: field(:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffer(:), sbuffer(:), wbuffer(:), nbuffer(:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete
  
  MPP_TYPE_                             :: field3D(size(field,1),size(field,2),1)
  MPP_TYPE_, allocatable, dimension(:,:) :: ebuffer2D, sbuffer2D, wbuffer2D, nbuffer2D
  integer                               :: xcount, ycount


  integer                  :: ntile
  logical                  :: need_ebuffer, need_sbuffer, need_wbuffer, need_nbuffer
  integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),  save :: f_addrs=-9999
  integer(LONG_KIND),dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_addrs=-9999
  integer, save    :: bsize(4)=0, isize=0, jsize=0, ksize=0, pos, list=0, l_size=0, upflags
  integer          :: buffer_size(4)
  integer          :: max_ntile, tile, update_position, ishift, jshift
  logical          :: do_update, is_complete, set_mismatch
  character(len=3) :: text
  MPP_TYPE_        :: d_type
  type(overlapSpec), pointer :: bound => NULL()

  ntile = size(domain%x(:))

  if(present(flags)) then
     call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_: flags is a dummy optional argument")
  endif
  update_position = CENTER
  if(present(position)) update_position = position

  !--- check if the buffer are needed
  need_ebuffer=.false.; need_sbuffer=.false.; need_wbuffer=.false.; need_nbuffer=.false.
  if( domain%symmetry .AND. PRESENT(position) ) then
     select case(position)
     case(CORNER)
       need_ebuffer=.true.; need_sbuffer=.true.; need_wbuffer=.true.; need_nbuffer=.true.
     case(NORTH)
       need_sbuffer=.true.; need_nbuffer=.true.
     case(EAST)
       need_ebuffer=.true.; need_wbuffer=.true.
     end select
  end if
   
  tile = 1
  max_ntile = domain%max_ntile_pe
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete        
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  f_addrs(list, tile) = LOC(field)
  if(present(ebuffer)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_2D: ebuffer should not be present when north is folded')
     if(.not. need_ebuffer) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D: ebuffer should not be present')
     b_addrs(1, list, tile) = LOC(ebuffer)
     buffer_size(1) = size(ebuffer(:))
  else
     b_addrs(1, list, tile) = 0
     buffer_size(1) = 1
  end if
  if(present(sbuffer)) then
     if(.not. need_sbuffer) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D: sbuffer should not be present')
     b_addrs(2, list, tile) = LOC(sbuffer)
     buffer_size(2) = size(sbuffer(:))
  else
     b_addrs(2, list, tile) = 0
     buffer_size(2) = 1
  end if
  if(present(wbuffer)) then
     if(.not. need_wbuffer) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D: wbuffer should not be present')
     b_addrs(3, list, tile) = LOC(wbuffer)
     buffer_size(3) = size(wbuffer(:))
  else
     b_addrs(3, list, tile) = 0
     buffer_size(3) = 1
  end if
  if(present(nbuffer)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_2D: nbuffer should not be present when north is folded')
     if(.not. need_nbuffer) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D: nbuffer should be be present')
     b_addrs(4, list, tile) = LOC(nbuffer)
     buffer_size(4) = size(nbuffer(:))
  else
     b_addrs(4, list, tile) = 0
     buffer_size(4) = 1
  end if

  if(list == 1 .AND. tile == 1 )then
     isize=size(field,1); jsize=size(field,2); ksize = 1; pos = update_position
     bsize = buffer_size
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize .NE. size(field,1))
     set_mismatch = set_mismatch .OR. (jsize .NE. size(field,2))
     set_mismatch = set_mismatch .OR. ANY( bsize .NE. buffer_size )
     set_mismatch = set_mismatch .OR. (update_position .NE. pos)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  if(is_complete) then
     l_size = list
     list = 0
  end if

  if(do_update )then 
     !--- only non-center data in symmetry domain will be retrieved.
     if(position == CENTER .OR. (.NOT. domain%symmetry) ) return 
     bound => search_bound_overlap(domain, update_position)
     call mpp_get_domain_shift(domain, ishift, jshift, update_position)
     if(size(field,1) .NE. domain%x(1)%memory%size+ishift .OR. size(field,2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D: field is not on memory domain")
    if(ASSOCIATED(bound)) then
        call mpp_do_get_boundary(f_addrs(1:l_size,1:ntile), domain, bound, b_addrs(:,1:l_size,1:ntile), &
             bsize, ksize, d_type)
     endif
     l_size=0; f_addrs=-9999; bsize=0; b_addrs=-9999; isize=0;  jsize=0;  ksize=0
  end if

  return

end subroutine MPP_GET_BOUNDARY_AD_2D_


!###############################################################################################
subroutine MPP_GET_BOUNDARY_AD_3D_(field, domain, ebuffer, sbuffer, wbuffer, nbuffer, flags, &
                                position, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: field(:,:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffer(:,:), sbuffer(:,:), wbuffer(:,:), nbuffer(:,:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete

  integer                  :: ntile
  logical                  :: need_ebuffer, need_sbuffer, need_wbuffer, need_nbuffer
  integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),  save :: f_addrs=-9999
  integer(LONG_KIND),dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_addrs=-9999
  integer, save    :: bsize(4)=0, isize=0, jsize=0, ksize=0, pos, list=0, l_size=0, upflags
  integer          :: buffer_size(4)
  integer          :: max_ntile, tile, update_position, ishift, jshift
  logical          :: do_update, is_complete, set_mismatch
  character(len=3) :: text
  MPP_TYPE_        :: d_type
  type(overlapSpec), pointer :: bound => NULL()

  ntile = size(domain%x(:))

  if(present(flags)) then
     call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D_: flags is a dummy optional argument")
  endif
  update_position = CENTER
  if(present(position)) update_position = position

  !--- check if the suitable buffer are present
  need_ebuffer=.false.; need_sbuffer=.false.; need_wbuffer=.false.; need_nbuffer=.false.
  if( domain%symmetry .AND. PRESENT(position) ) then
     select case(position)
     case(CORNER)
       need_ebuffer=.true.; need_sbuffer=.true.; need_wbuffer=.true.; need_nbuffer=.true.
     case(NORTH)
       need_sbuffer=.true.; need_nbuffer=.true.
     case(EAST)
       need_ebuffer=.true.; need_wbuffer=.true.
     end select
  end if
   
  tile = 1
  max_ntile = domain%max_ntile_pe
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete        
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  f_addrs(list, tile) = LOC(field)
  if(present(ebuffer)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_3D: ebuffer should not be present when north is folded')
     if(.not. need_ebuffer) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: ebuffer should not be present')
     b_addrs(1, list, tile) = LOC(ebuffer)
     buffer_size(1) = size(ebuffer,1)
  else
     b_addrs(1, list, tile) = 0
     buffer_size(1) = 1
  end if
  if(present(sbuffer)) then
     if(.not. need_sbuffer) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: sbuffer should not be present')
     b_addrs(2, list, tile) = LOC(sbuffer)
     buffer_size(2) = size(sbuffer,1)
  else
     b_addrs(2, list, tile) = 0
     buffer_size(2) = 1
  end if
  if(present(wbuffer)) then
     if(.not. need_wbuffer) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: wbuffer should not be present')
     b_addrs(3, list, tile) = LOC(wbuffer)
     buffer_size(3) = size(wbuffer,1)
  else
     b_addrs(3, list, tile) = 0
     buffer_size(3) = 1
  end if
  if(present(nbuffer)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_3D: nbuffer should not be present when north is folded')
     if(.not. need_nbuffer) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: nbuffer should not be present')
     b_addrs(4, list, tile) = LOC(nbuffer)
     buffer_size(4) = size(nbuffer,1)
  else
     b_addrs(4, list, tile) = 0
     buffer_size(4) = 1
  end if


  if(list == 1 .AND. tile == 1 )then
     isize=size(field,1); jsize=size(field,2); ksize = size(field,3); pos = update_position
     bsize = buffer_size
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize .NE. size(field,1))
     set_mismatch = set_mismatch .OR. (jsize .NE. size(field,2))
     set_mismatch = set_mismatch .OR. (ksize .NE. size(field,3))
     set_mismatch = set_mismatch .OR. ANY( bsize .NE. buffer_size )
     set_mismatch = set_mismatch .OR. (update_position .NE. pos)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  if(is_complete) then
     l_size = list
     list = 0
  end if

  if(do_update )then 
     !--- only non-center data in symmetry domain will be retrieved.
     if(position == CENTER .OR. (.NOT. domain%symmetry) ) return 
     bound => search_bound_overlap(domain, update_position)
     call mpp_get_domain_shift(domain, ishift, jshift, update_position)
     if(size(field,1) .NE. domain%x(1)%memory%size+ishift .OR. size(field,2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D: field is not on memory domain")
    if(ASSOCIATED(bound)) then
        call mpp_do_get_boundary(f_addrs(1:l_size,1:ntile), domain, bound, b_addrs(:,1:l_size,1:ntile), &
             bsize, ksize, d_type)
     endif
     l_size=0; f_addrs=-9999; bsize=0; b_addrs=-9999; isize=0;  jsize=0;  ksize=0
  end if

end subroutine MPP_GET_BOUNDARY_AD_3D_


!####################################################################
! vector update
subroutine MPP_GET_BOUNDARY_AD_2D_V_(fieldx, fieldy, domain, ebufferx, sbufferx, wbufferx, nbufferx, &
                                  ebuffery, sbuffery, wbuffery, nbuffery, flags, gridtype, &
                                  complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: fieldx(:,:), fieldy(:,:)
  MPP_TYPE_, intent(inout), optional :: ebufferx(:), sbufferx(:), wbufferx(:), nbufferx(:)
  MPP_TYPE_, intent(inout), optional :: ebuffery(:), sbuffery(:), wbuffery(:), nbuffery(:)
  integer,      intent(in), optional :: flags, gridtype, tile_count
  logical,      intent(in), optional :: complete
  
  integer                 :: ntile, update_flags
  logical                 :: need_ebufferx, need_sbufferx, need_wbufferx, need_nbufferx
  logical                 :: need_ebuffery, need_sbuffery, need_wbuffery, need_nbuffery

  integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),  save :: f_addrsx=-9999
  integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),  save :: f_addrsy=-9999
  integer(LONG_KIND),dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_addrsx=-9999
  integer(LONG_KIND),dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_addrsy=-9999
  integer, save    :: bsizex(4)=0, bsizey(4)=0, isize(2)=0, jsize(2)=0, ksize=0, l_size=0, list=0
  integer, save    :: offset_type, upflags
  integer          :: bufferx_size(4), buffery_size(4)
  integer          :: max_ntile, tile, grid_offset_type
  logical          :: do_update, is_complete, set_mismatch
  character(len=3) :: text
  MPP_TYPE_        :: d_type
  type(overlapSpec), pointer :: boundx=>NULL()
  type(overlapSpec), pointer :: boundy=>NULL()
  integer                     :: position_x, position_y, ishift, jshift

  ntile = size(domain%x(:))
  update_flags = 0
  if( PRESENT(flags) ) then 
     update_flags = flags
  end if

  !--- check if the suitable buffer are present
  need_ebufferx=.FALSE.; need_sbufferx=.FALSE.
  need_wbufferx=.FALSE.; need_nbufferx=.FALSE.
  need_ebuffery=.FALSE.; need_sbuffery=.FALSE.
  need_wbuffery=.FALSE.; need_nbuffery=.FALSE.
  if( domain%symmetry .AND. PRESENT(gridtype) ) then
     select case(gridtype)
     case(BGRID_NE, BGRID_SW)
       need_ebufferx=.true.; need_sbufferx=.true.; need_wbufferx=.true.; need_nbufferx=.true.
       need_ebuffery=.true.; need_sbuffery=.true.; need_wbuffery=.true.; need_nbuffery=.true.
     case(CGRID_NE, CGRID_SW)
       need_ebufferx=.true.; need_wbufferx=.true.; need_sbuffery=.true.; need_nbuffery=.true.
     case(DGRID_NE, DGRID_SW)
       need_ebuffery=.true.; need_wbuffery=.true.; need_sbufferx=.true.; need_nbufferx=.true.
     end select
   end if

  tile = 1
  max_ntile = domain%max_ntile_pe
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete        
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  f_addrsx(list, tile) = LOC(fieldx)
  f_addrsy(list, tile) = LOC(fieldy)

  if(present(ebufferx)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_2D_V: ebufferx should not be present when north is folded')
     if(.not. need_ebufferx) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: ebufferx should not be present')
     b_addrsx(1, list, tile) = LOC(ebufferx)
     bufferx_size(1) = size(ebufferx,1)
  else
     b_addrsx(1, list, tile) = 0
     bufferx_size(1) = 1
  end if
  if(present(sbufferx)) then
     if(.not. need_sbufferx) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: sbufferx should not be present')
     b_addrsx(2, list, tile) = LOC(sbufferx)
     bufferx_size(2) = size(sbufferx,1)
  else
     b_addrsx(2, list, tile) = 0
     bufferx_size(2) = 1
  end if
  if(present(wbufferx)) then
     if(.not. need_wbufferx) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: wbufferx should not be present')
     b_addrsx(3, list, tile) = LOC(wbufferx)
     bufferx_size(3) = size(wbufferx,1)
  else
     b_addrsx(3, list, tile) = 0
     bufferx_size(3) = 1
  end if
  if(present(nbufferx)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_2D_V: nbufferx should not be present when north is folded')
     if(.not. need_nbufferx) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: nbufferx should not be present')
     b_addrsx(4, list, tile) = LOC(nbufferx)
     bufferx_size(4) = size(nbufferx,1)
  else
     b_addrsx(4, list, tile) = 0
     bufferx_size(4) = 1
  end if

  if(present(ebuffery)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_2D_V: ebuffery should not be present when north is folded')
     if(.not. need_ebuffery) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: ebuffery should not be present')
     b_addrsy(1, list, tile) = LOC(ebuffery)
     buffery_size(1) = size(ebuffery,1)
  else
     b_addrsy(1, list, tile) = 0
     buffery_size(1) = 1
  end if
  if(present(sbuffery)) then
     if(.not. need_sbuffery) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: sbuffery should not be present')
     b_addrsy(2, list, tile) = LOC(sbuffery)
     buffery_size(2) = size(sbuffery,1)
  else
     b_addrsy(2, list, tile) = 0
     buffery_size(2) = 1
  end if
  if(present(wbuffery)) then
     if(.not. need_wbuffery) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: wbuffery should not be present')
     b_addrsy(3, list, tile) = LOC(wbuffery)
     buffery_size(3) = size(wbuffery,1)
  else
     b_addrsy(3, list, tile) = 0
     buffery_size(3) = 1
  end if
  if(present(nbuffery)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_2D_V: nbuffery should not be present when north is folded')
     if(.not. need_nbuffery) call mpp_error(FATAL,'MPP_GET_BOUNDARY_2D_V: nbuffery should not be present')
     b_addrsy(4, list, tile) = LOC(nbuffery)
     buffery_size(4) = size(nbuffery,1)
  else
     b_addrsy(4, list, tile) = 0
     buffery_size(4) = 1
  end if

  grid_offset_type = AGRID
  if(present(gridtype)) grid_offset_type = gridtype
  if(list == 1 .AND. tile == 1 )then
     isize(1)=size(fieldx,1); jsize(1)=size(fieldx,2); isize(2)=size(fieldy,1); jsize(2)=size(fieldy,2)
     ksize = 1; offset_type = grid_offset_type
     bsizex = bufferx_size; bsizey = buffery_size; upflags = update_flags
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize(1) .NE. size(fieldx,1))
     set_mismatch = set_mismatch .OR. (jsize(1) .NE. size(fieldx,2))
     set_mismatch = set_mismatch .OR. (isize(2) .NE. size(fieldy,1))
     set_mismatch = set_mismatch .OR. (jsize(2) .NE. size(fieldy,2))
     set_mismatch = set_mismatch .OR. ANY( bsizex .NE. bufferx_size )
     set_mismatch = set_mismatch .OR. ANY( bsizey .NE. buffery_size )
     set_mismatch = set_mismatch .OR. (offset_type .NE. grid_offset_type)
     set_mismatch = set_mismatch .OR. (upflags .NE. update_flags)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  if(is_complete) then
     l_size = list
     list = 0
  end if

  if(do_update )then
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
     case (DGRID_NE, DGRID_SW)
        position_x = NORTH
        position_y = EAST
     case default
        call mpp_error(FATAL, "mpp_get_boundary.h: invalid value of grid_offset_type")
     end select

     boundx => search_bound_overlap(domain, position_x)
     boundy => search_bound_overlap(domain, position_y)  

     call mpp_get_domain_shift(domain, ishift, jshift, position_x)
     if(size(fieldx,1) .NE. domain%x(1)%memory%size+ishift .OR. size(fieldx,2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: fieldx is not on memory domain")
     call mpp_get_domain_shift(domain, ishift, jshift, position_y)
     if(size(fieldy,1) .NE. domain%x(1)%memory%size+ishift .OR. size(fieldy,2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: fieldy is not on memory domain")
     if(ASSOCIATED(boundx) ) then
        call mpp_do_get_boundary(f_addrsx(1:l_size,1:ntile), f_addrsy(1:l_size,1:ntile), domain, boundx, boundy, &
             b_addrsx(:,1:l_size,1:ntile), b_addrsy(:,1:l_size,1:ntile), bsizex, &
             bsizey, ksize, d_type, update_flags, grid_offset_type)
     endif
     l_size=0; f_addrsx=-9999; f_addrsy=-9999; bsizex=0; bsizey=0; 
     b_addrsx=-9999; b_addrsy=-9999; isize=0;  jsize=0;  ksize=0
  end if


  return

end subroutine MPP_GET_BOUNDARY_AD_2D_V_


!###############################################################################################
subroutine MPP_GET_BOUNDARY_AD_3D_V_(fieldx, fieldy, domain, ebufferx, sbufferx, wbufferx, nbufferx, &
                                  ebuffery, sbuffery, wbuffery, nbuffery, flags, gridtype, &
                                  complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: fieldx(domain%x(1)%memory%begin:,domain%y(1)%memory%begin:,:)
  MPP_TYPE_,            intent(in)   :: fieldy(domain%x(1)%memory%begin:,domain%y(1)%memory%begin:,:)
  MPP_TYPE_, intent(inout), optional :: ebufferx(:,:), sbufferx(:,:), wbufferx(:,:), nbufferx(:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffery(:,:), sbuffery(:,:), wbuffery(:,:), nbuffery(:,:)
  integer,      intent(in), optional :: flags, gridtype, tile_count
  logical,      intent(in), optional :: complete

  integer                 :: ntile, update_flags
  logical                 :: need_ebufferx, need_sbufferx, need_wbufferx, need_nbufferx
  logical                 :: need_ebuffery, need_sbuffery, need_wbuffery, need_nbuffery

  integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),  save :: f_addrsx=-9999
  integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),  save :: f_addrsy=-9999
  integer(LONG_KIND),dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_addrsx=-9999
  integer(LONG_KIND),dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_addrsy=-9999
  integer, save    :: bsizex(4)=0, bsizey(4)=0, isize(2)=0, jsize(2)=0, ksize=0, l_size=0, list=0
  integer, save    :: offset_type, upflags
  integer          :: bufferx_size(4), buffery_size(4)
  integer          :: max_ntile, tile, grid_offset_type
  logical          :: do_update, is_complete, set_mismatch
  character(len=3) :: text
  MPP_TYPE_        :: d_type
  type(overlapSpec), pointer :: boundx=>NULL()
  type(overlapSpec), pointer :: boundy=>NULL()
  integer                     :: position_x, position_y, ishift, jshift

  ntile = size(domain%x(:))
  update_flags = 0
  if( PRESENT(flags) ) then 
     update_flags = flags
  end if

  !--- check if the suitable buffer are present
  need_ebufferx=.FALSE.; need_sbufferx=.FALSE.
  need_wbufferx=.FALSE.; need_nbufferx=.FALSE.
  need_ebuffery=.FALSE.; need_sbuffery=.FALSE.
  need_wbuffery=.FALSE.; need_nbuffery=.FALSE.
  if( domain%symmetry .AND. PRESENT(gridtype) ) then
     select case(gridtype)
     case(BGRID_NE, BGRID_SW)
       need_ebufferx=.true.; need_sbufferx=.true.; need_wbufferx=.true.; need_nbufferx=.true.
       need_ebuffery=.true.; need_sbuffery=.true.; need_wbuffery=.true.; need_nbuffery=.true.
     case(CGRID_NE, CGRID_SW)
       need_ebufferx=.true.; need_wbufferx=.true.; need_sbuffery=.true.; need_nbuffery=.true.
     case(DGRID_NE, DGRID_SW)
       need_ebuffery=.true.; need_wbuffery=.true.; need_sbufferx=.true.; need_nbufferx=.true.
     end select
   end if

  tile = 1
  max_ntile = domain%max_ntile_pe
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D_V: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete        
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  f_addrsx(list, tile) = LOC(fieldx)
  f_addrsy(list, tile) = LOC(fieldy)

  if(present(ebufferx)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_3D_V: ebufferx should not be present when north is folded')
     if(.not. need_ebufferx) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: ebufferx should not be present')
     b_addrsx(1, list, tile) = LOC(ebufferx)
     bufferx_size(1) = size(ebufferx,1)
  else
     b_addrsx(1, list, tile) = 0
     bufferx_size(1) = 1
  end if
  if(present(sbufferx)) then
     if(.not. need_sbufferx) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: sbufferx should not be present')
     b_addrsx(2, list, tile) = LOC(sbufferx)
     bufferx_size(2) = size(sbufferx,1)
  else
     b_addrsx(2, list, tile) = 0
     bufferx_size(2) = 1
  end if
  if(present(wbufferx)) then
     if(.not. need_wbufferx) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: wbufferx should not be present')
     b_addrsx(3, list, tile) = LOC(wbufferx)
     bufferx_size(3) = size(wbufferx,1)
  else
     b_addrsx(3, list, tile) = 0
     bufferx_size(3) = 1
  end if
  if(present(nbufferx)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_3D_V: nbufferx should not be present when north is folded')
     if(.not. need_nbufferx) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: nbufferx should not be present')
     b_addrsx(4, list, tile) = LOC(nbufferx)
     bufferx_size(4) = size(nbufferx,1)
  else
     b_addrsx(4, list, tile) = 0
     bufferx_size(4) = 1
  end if

  if(present(ebuffery)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_3D_V: ebuffery should not be present when north is folded')
     if(.not. need_ebuffery) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: ebuffery should not be present')
     b_addrsy(1, list, tile) = LOC(ebuffery)
     buffery_size(1) = size(ebuffery,1)
  else
     b_addrsy(1, list, tile) = 0
     buffery_size(1) = 1
  end if
  if(present(sbuffery)) then
     if(.not. need_sbuffery) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: sbuffery should not be present')
     b_addrsy(2, list, tile) = LOC(sbuffery)
     buffery_size(2) = size(sbuffery,1)
  else
     b_addrsy(2, list, tile) = 0
     buffery_size(2) = 1
  end if
  if(present(wbuffery)) then
     if(.not. need_wbuffery) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: wbuffery should not be present')
     b_addrsy(3, list, tile) = LOC(wbuffery)
     buffery_size(3) = size(wbuffery,1)
  else
     b_addrsy(3, list, tile) = 0
     buffery_size(3) = 1
  end if
  if(present(nbuffery)) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'MPP_GET_BOUNDARY_3D_V: nbuffery should not be present when north is folded')
     if(.not. need_nbuffery) call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: nbuffery should not be present')
     b_addrsy(4, list, tile) = LOC(nbuffery)
     buffery_size(4) = size(nbuffery,1)
  else
     b_addrsy(4, list, tile) = 0
     buffery_size(4) = 1
  end if

  grid_offset_type = AGRID
  if(present(gridtype)) grid_offset_type = gridtype
  if(list == 1 .AND. tile == 1 )then
     isize(1)=size(fieldx,1); jsize(1)=size(fieldx,2); isize(2)=size(fieldy,1); jsize(2)=size(fieldy,2)
     ksize = size(fieldx,3); offset_type = grid_offset_type
     bsizex = bufferx_size; bsizey = buffery_size; upflags = update_flags
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize(1) .NE. size(fieldx,1))
     set_mismatch = set_mismatch .OR. (jsize(1) .NE. size(fieldx,2))
     set_mismatch = set_mismatch .OR. (ksize    .NE. size(fieldx,3))
     set_mismatch = set_mismatch .OR. (isize(2) .NE. size(fieldy,1))
     set_mismatch = set_mismatch .OR. (jsize(2) .NE. size(fieldy,2))
     set_mismatch = set_mismatch .OR. (ksize    .NE. size(fieldy,3))
     set_mismatch = set_mismatch .OR. ANY( bsizex .NE. bufferx_size )
     set_mismatch = set_mismatch .OR. ANY( bsizey .NE. buffery_size )
     set_mismatch = set_mismatch .OR. (offset_type .NE. grid_offset_type)
     set_mismatch = set_mismatch .OR. (upflags .NE. update_flags)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  if(is_complete) then
     l_size = list
     list = 0
  end if

  if(do_update )then
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
     case (DGRID_NE, DGRID_SW)
        position_x = NORTH
        position_y = EAST
     case default
        call mpp_error(FATAL, "mpp_get_boundary.h: invalid value of grid_offset_type")
     end select

     boundx => search_bound_overlap(domain, position_x)
     boundy => search_bound_overlap(domain, position_y)  

     call mpp_get_domain_shift(domain, ishift, jshift, position_x)
     if(size(fieldx,1) .NE. domain%x(1)%memory%size+ishift .OR. size(fieldx,2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D_V: fieldx is not on memory domain")
     call mpp_get_domain_shift(domain, ishift, jshift, position_y)
     if(size(fieldy,1) .NE. domain%x(1)%memory%size+ishift .OR. size(fieldy,2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D_V: fieldy is not on memory domain")
     if(ASSOCIATED(boundx) ) then
        call mpp_do_get_boundary_ad(f_addrsx(1:l_size,1:ntile), f_addrsy(1:l_size,1:ntile), domain, boundx, boundy, &
             b_addrsx(:,1:l_size,1:ntile), b_addrsy(:,1:l_size,1:ntile), bsizex, &
             bsizey, ksize, d_type, update_flags, grid_offset_type)
     endif
     l_size=0; f_addrsx=-9999; f_addrsy=-9999; bsizex=0; bsizey=0; 
     b_addrsx=-9999; b_addrsy=-9999; isize=0;  jsize=0;  ksize=0
  end if

end subroutine MPP_GET_BOUNDARY_AD_3D_V_

