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
subroutine MPP_UPDATE_NEST_FINE_2D_(field, nest_domain, wbuffer, ebuffer, sbuffer, nbuffer, &
                                    nest_level, flags, complete, position, extra_halo, name, tile_count)
      MPP_TYPE_,             intent(in)      :: field(:,:) !< field on the model grid
      type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                            !! between fine and coarse grid.
      MPP_TYPE_,             intent(inout)   :: wbuffer(:,:) !< west side buffer to be filled
                                                             !! with data on coarse grid.
      MPP_TYPE_,             intent(inout)   :: ebuffer(:,:) !< east side buffer to be filled
                                                             !! with data on coarse grid.
      MPP_TYPE_,             intent(inout)   :: sbuffer(:,:) !< south side buffer to be filled
                                                             !! with data on coarse grid.
      MPP_TYPE_,             intent(inout)   :: nbuffer(:,:) !< north side buffer to be filled
                                                             !! with data on coarse grid.
      integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
      integer,          intent(in), optional :: flags !< Specify the direction of fine grid halo buffer to be filled.
                                                      !! Default value is XUPDATE+YUPDATE.
      logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
                                                         !! Default value is .true.
      integer,          intent(in), optional :: position !< Cell position. It value should be
                                                         !! CENTER, EAST, CORNER, or NORTH. Default is CENTER.
      integer,          intent(in), optional :: extra_halo !< extra halo for passing data
                                                           !! from coarse grid to fine grid.
                                                           !! Default is 0 and currently only support extra_halo = 0.
      character(len=*), intent(in), optional :: name !< Name of the nest domain.
      integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                           !! default is 1 and currently only support tile_count = 1.

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
      MPP_TYPE_ :: wbuffer3D(size(wbuffer,1),size(wbuffer,2),1)
      MPP_TYPE_ :: ebuffer3D(size(ebuffer,1),size(ebuffer,2),1)
      MPP_TYPE_ :: sbuffer3D(size(sbuffer,1),size(sbuffer,2),1)
      MPP_TYPE_ :: nbuffer3D(size(nbuffer,1),size(nbuffer,2),1)
      pointer( ptr, field3D )
      pointer( ptr_w, wbuffer3D)
      pointer( ptr_e, ebuffer3D)
      pointer( ptr_s, sbuffer3D)
      pointer( ptr_n, nbuffer3D)
      ptr = LOC(field)
      ptr_w = LOC(wbuffer)
      ptr_e = LOC(ebuffer)
      ptr_s = LOC(sbuffer)
      ptr_n = LOC(nbuffer)
      call mpp_update_nest_fine( field3D, nest_domain, wbuffer3D, ebuffer3D, sbuffer3D, nbuffer3D, &
                                 nest_level, flags, complete, position, extra_halo, name, tile_count)

      return


end subroutine MPP_UPDATE_NEST_FINE_2D_

subroutine MPP_UPDATE_NEST_FINE_3D_(field, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, &
                                    nest_level, flags, complete, position, extra_halo, name, tile_count)
    MPP_TYPE_,             intent(in)      :: field(:,:,:) !< field on the model grid
    type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                          !! between fine and coarse grid.
    MPP_TYPE_,             intent(inout)   :: wbuffer(:,:,:) !< west side buffer to be filled
                                                             !! with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: ebuffer(:,:,:) !< east side buffer to be filled
                                                             !! with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: sbuffer(:,:,:) !< south side buffer to be filled
                                                             !! with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: nbuffer(:,:,:) !< north side buffer to be filled
                                                             !! with data on coarse grid.
    integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
    integer,          intent(in), optional :: flags !< Specify the direction of fine grid halo buffer to be filled.
                                                    !! Default value is XUPDATE+YUPDATE.
    logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
                                                       !! Default value is .true.
    integer,          intent(in), optional :: position !< Cell position. It value should be
                                                       !! CENTER, EAST, CORNER, or NORTH. Default is CENTER.
    integer,          intent(in), optional :: extra_halo !< extra halo for passing data
                                                         !! from coarse grid to fine grid.
                                                         !! Default is 0 and currently only support extra_halo = 0.
    character(len=*), intent(in), optional :: name !< Name of the nest domain.
    integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                         !! default is 1 and currently only support tile_count = 1.

   MPP_TYPE_        :: d_type
   type(nestSpec), pointer :: update=>NULL()
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: f_addrs=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: wb_addrs=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: eb_addrs=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: sb_addrs=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: nb_addrs=-9999
   character(len=3) :: text
   logical          :: is_complete, set_mismatch
   integer          :: tile
   integer          :: add_halo, update_flags, update_position
   integer          :: wbuffersz, ebuffersz, sbuffersz, nbuffersz
   integer          :: isize, jsize, ksize, l_size
   integer, save    :: isize_save=0, jsize_save=0, ksize_save=0
   integer, save    :: wbuffersz_save=0, ebuffersz_save=0, sbuffersz_save=0, nbuffersz_save=0
   integer, save    :: add_halo_save=0, update_flags_save=0, update_position_save=0
   integer, save    :: list=0
   integer          :: xbegin, xend, ybegin, yend

   add_halo = 0
   if(present(extra_halo)) add_halo = add_halo
   update_position = CENTER
   if(present(position)) update_position = position
   update_flags = XUPDATE+YUPDATE   !default
   if( PRESENT(flags) )update_flags = flags


   is_complete = .true.
   if(PRESENT(complete)) then
      is_complete = complete
   end if
   tile = 1
   if(present(tile_count)) tile = tile_count
   if( tile > 1 ) then
      call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D: currently do not support multiple tile per pe')
   endif

   list = list+1
   if(list > MAX_DOMAIN_FIELDS)then
      write( text,'(i2)' ) MAX_DOMAIN_FIELDS
      call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
   endif

   f_addrs(list) = LOC(field)
   wb_addrs(list) = LOC(wbuffer)
   eb_addrs(list) = LOC(ebuffer)
   sb_addrs(list) = LOC(sbuffer)
   nb_addrs(list) = LOC(nbuffer)

   wbuffersz = size(wbuffer); ebuffersz = size(ebuffer)
   sbuffersz = size(sbuffer); nbuffersz = size(nbuffer)
   isize=size(field,1); jsize=size(field,2); ksize = size(field,3)

   !---check data is on data domain or compute domain
   if( nest_domain%nest(nest_level)%is_coarse_pe ) then
      call mpp_get_data_domain(nest_domain%nest(nest_level)%domain_coarse, xbegin, xend, ybegin, yend, position=update_position)
      if(isize .NE. xend-xbegin+1 .OR. jsize .NE. yend-ybegin+1) then
         call mpp_get_compute_domain(nest_domain%nest(nest_level)%domain_coarse, xbegin, xend, ybegin, yend, &
                                     position=update_position)
         if(isize .NE. xend-xbegin+1 .OR. jsize .NE. yend-ybegin+1) then
            call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D_: field is neither on data domain nor on compute domain')
         endif
      endif
   else
      xbegin = 1; xend = 1
      ybegin = 1; yend = 1
   endif
   if(list == 1)then
      isize_save = isize; jsize_save = jsize; ksize_save = ksize
      update_position_save = update_position
      update_flags_save    = update_flags
      wbuffersz_save = wbuffersz; ebuffersz_save = ebuffersz
      sbuffersz_save = sbuffersz; nbuffersz_save = nbuffersz
      add_halo_save = add_halo
   else
      set_mismatch = .false.
      set_mismatch = set_mismatch .OR. (isize_save /= isize)
      set_mismatch = set_mismatch .OR. (jsize_save /= jsize)
      set_mismatch = set_mismatch .OR. (ksize_save /= ksize)
      set_mismatch = set_mismatch .OR. (update_position_save /= update_position)
      set_mismatch = set_mismatch .OR. (wbuffersz_save /= wbuffersz)
      set_mismatch = set_mismatch .OR. (ebuffersz_save /= ebuffersz)
      set_mismatch = set_mismatch .OR. (sbuffersz_save /= sbuffersz)
      set_mismatch = set_mismatch .OR. (nbuffersz_save /= nbuffersz)
      set_mismatch = set_mismatch .OR. (update_flags_save /= update_flags)
      set_mismatch = set_mismatch .OR. (add_halo_save /= add_halo)

      if(set_mismatch)then
         write( text,'(i2)' ) list
         call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D_: Incompatible field at count '//text//' for group update.' )
      endif
   endif

   if(is_complete) then
      l_size = list
      list = 0
   end if

   if(is_complete)then
      update => search_C2F_nest_overlap(nest_domain, nest_level, add_halo, update_position)
      call mpp_do_update_nest_fine(f_addrs(1:l_size), nest_domain%nest(nest_level), update, d_type, ksize, &
            wb_addrs(1:l_size), eb_addrs(1:l_size), sb_addrs(1:l_size), nb_addrs(1:l_size), update_flags, &
            xbegin, xend, ybegin, yend )
   endif


end subroutine MPP_UPDATE_NEST_FINE_3D_


!###############################################################################
subroutine MPP_UPDATE_NEST_FINE_4D_(field, nest_domain, wbuffer, ebuffer, sbuffer, nbuffer, &
                                    nest_level, flags, complete, position, extra_halo, name, tile_count)
      MPP_TYPE_,             intent(in)      :: field(:,:,:,:) !< field on the model grid
      type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                            !! between fine and coarse grid.
      MPP_TYPE_,             intent(inout)   :: wbuffer(:,:,:,:) !< west side buffer to be filled
                                                                 !! with data on coarse grid.
      MPP_TYPE_,             intent(inout)   :: ebuffer(:,:,:,:) !< east side buffer to be filled
                                                                 !! with data on coarse grid.
      MPP_TYPE_,             intent(inout)   :: sbuffer(:,:,:,:) !< south side buffer to be filled
                                                                 !! with data on coarse grid.
      MPP_TYPE_,             intent(inout)   :: nbuffer(:,:,:,:) !< north side buffer to be filled
                                                                 !! with data on coarse grid.
      integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
      integer,          intent(in), optional :: flags  !< Specify the direction of fine grid halo buffer to be filled.
                                                    !! Default value is XUPDATE+YUPDATE.
      logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
                                                         !! Default value is .true.
      integer,          intent(in), optional :: position !< Cell position. It value should be
                                                         !! CENTER, EAST, CORNER, or NORTH. Default is CENTER.
      integer,          intent(in), optional :: extra_halo !< extra halo for passing data
                                                           !! from coarse grid to fine grid.
                                                           !! Default is 0 and currently only support extra_halo = 0.
      character(len=*), intent(in), optional :: name !< Name of the nest domain.
      integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                           !! default is 1 and currently only support tile_count = 1.

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
      MPP_TYPE_ :: wbuffer3D(size(wbuffer,1),size(wbuffer,2),size(wbuffer,3)*size(wbuffer,4))
      MPP_TYPE_ :: ebuffer3D(size(ebuffer,1),size(ebuffer,2),size(ebuffer,3)*size(ebuffer,4))
      MPP_TYPE_ :: sbuffer3D(size(sbuffer,1),size(sbuffer,2),size(sbuffer,3)*size(sbuffer,4))
      MPP_TYPE_ :: nbuffer3D(size(nbuffer,1),size(nbuffer,2),size(nbuffer,3)*size(nbuffer,4))

      pointer( ptr, field3D )
      pointer( ptr_w, wbuffer3D)
      pointer( ptr_e, ebuffer3D)
      pointer( ptr_s, sbuffer3D)
      pointer( ptr_n, nbuffer3D)
      ptr = LOC(field)
      ptr_w = LOC(wbuffer)
      ptr_e = LOC(ebuffer)
      ptr_s = LOC(sbuffer)
      ptr_n = LOC(nbuffer)
      call mpp_update_nest_fine( field3D, nest_domain, wbuffer3D, ebuffer3D, sbuffer3D, nbuffer3D, &
                                 nest_level, flags, complete, position, extra_halo, name, tile_count)

      return


end subroutine MPP_UPDATE_NEST_FINE_4D_

#ifdef VECTOR_FIELD_
subroutine MPP_UPDATE_NEST_FINE_2D_V_(fieldx, fieldy, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, &
                                      ebufferx, ebuffery, nbufferx, nbuffery, nest_level, &
                                      flags, gridtype, complete, extra_halo, name, tile_count)
    MPP_TYPE_,             intent(in)      :: fieldx(:,:), fieldy(:,:) !< field x and y components on the model grid
    type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                          !! between fine and coarse grid.
    MPP_TYPE_,             intent(inout)   :: wbufferx(:,:), wbuffery(:,:) !< west side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: ebufferx(:,:), ebuffery(:,:) !< east side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: sbufferx(:,:), sbuffery(:,:) !< south side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: nbufferx(:,:), nbuffery(:,:) !< north side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
    integer,          intent(in), optional :: flags !< Specify the direction of fine grid halo buffer to be filled.
                                                    !! Default value is XUPDATE+YUPDATE.
    logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
                                                       !! Default value is .true.
    integer,          intent(in), optional :: gridtype
    integer,          intent(in), optional :: extra_halo !< extra halo for passing data
                                                         !! from coarse grid to fine grid.
                                                         !! Default is 0 and currently only support extra_halo = 0.
    character(len=*), intent(in), optional :: name !< Name of the nest domain.
    integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                         !! default is 1 and currently only support tile_count = 1.

    MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
    MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)
    MPP_TYPE_ :: wbuffer3Dx(size(wbufferx,1),size(wbufferx,2),1)
    MPP_TYPE_ :: ebuffer3Dx(size(ebufferx,1),size(ebufferx,2),1)
    MPP_TYPE_ :: sbuffer3Dx(size(sbufferx,1),size(sbufferx,2),1)
    MPP_TYPE_ :: nbuffer3Dx(size(nbufferx,1),size(nbufferx,2),1)
    MPP_TYPE_ :: wbuffer3Dy(size(wbuffery,1),size(wbuffery,2),1)
    MPP_TYPE_ :: ebuffer3Dy(size(ebuffery,1),size(ebuffery,2),1)
    MPP_TYPE_ :: sbuffer3Dy(size(sbuffery,1),size(sbuffery,2),1)
    MPP_TYPE_ :: nbuffer3Dy(size(nbuffery,1),size(nbuffery,2),1)
    pointer( ptrx, field3Dx )
    pointer( ptry, field3Dy )
    pointer( ptr_wx, wbuffer3Dx)
    pointer( ptr_ex, ebuffer3Dx)
    pointer( ptr_sx, sbuffer3Dx)
    pointer( ptr_nx, nbuffer3Dx)
    pointer( ptr_wy, wbuffer3Dy)
    pointer( ptr_ey, ebuffer3Dy)
    pointer( ptr_sy, sbuffer3Dy)
    pointer( ptr_ny, nbuffer3Dy)

    ptrx = LOC(fieldx)
    ptry = LOC(fieldy)
    ptr_wx = LOC(wbufferx)
    ptr_ex = LOC(ebufferx)
    ptr_sx = LOC(sbufferx)
    ptr_nx = LOC(nbufferx)
    ptr_wy = LOC(wbuffery)
    ptr_ey = LOC(ebuffery)
    ptr_sy = LOC(sbuffery)
    ptr_ny = LOC(nbuffery)

    call MPP_UPDATE_NEST_FINE_3D_V_(field3Dx, field3Dy, nest_domain, wbuffer3Dx, wbuffer3Dy, sbuffer3Dx, sbuffer3Dy, &
                                    ebuffer3Dx, ebuffer3Dy, nbuffer3Dx, nbuffer3Dy, nest_level, &
                                      flags, gridtype, complete, extra_halo, name, tile_count)

end subroutine MPP_UPDATE_NEST_FINE_2D_V_

subroutine MPP_UPDATE_NEST_FINE_3D_V_(fieldx, fieldy, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, &
                                      ebufferx, ebuffery, nbufferx, nbuffery, nest_level, &
                                      flags, gridtype, complete, extra_halo, name, tile_count)
    MPP_TYPE_,             intent(in)      :: fieldx(:,:,:), fieldy(:,:,:) !< field x and y components on the model grid
    type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                          !! between fine and coarse grid.
    MPP_TYPE_,             intent(inout)   :: wbufferx(:,:,:), wbuffery(:,:,:) !< west side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: ebufferx(:,:,:), ebuffery(:,:,:) !< east side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: sbufferx(:,:,:), sbuffery(:,:,:) !< south side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: nbufferx(:,:,:), nbuffery(:,:,:) !< north side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
    integer,          intent(in), optional :: flags !< Specify the direction of fine grid halo buffer to be filled.
                                                    !! Default value is XUPDATE+YUPDATE.
    logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
                                                       !! Default value is .true.
    integer,          intent(in), optional :: gridtype
    integer,          intent(in), optional :: extra_halo !< extra halo for passing data
                                                         !! from coarse grid to fine grid.
                                                         !! Default is 0 and currently only support extra_halo = 0.
    character(len=*), intent(in), optional :: name !< Name of the nest domain.
    integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                         !! default is 1 and currently only support tile_count = 1.

   MPP_TYPE_        :: d_type
   type(nestSpec), pointer :: updatex=>NULL()
   type(nestSpec), pointer :: updatey=>NULL()
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: f_addrsx=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: f_addrsy=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: wb_addrsx=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: eb_addrsx=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: sb_addrsx=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: nb_addrsx=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: wb_addrsy=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: eb_addrsy=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: sb_addrsy=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: nb_addrsy=-9999

   character(len=3) :: text
   logical          :: is_complete, set_mismatch
   integer          :: tile
   integer          :: add_halo, update_flags, grid_offset_type
   integer          :: wbufferszx, ebufferszx, sbufferszx, nbufferszx
   integer          :: wbufferszy, ebufferszy, sbufferszy, nbufferszy
   integer          :: isizex, jsizex, isizey, jsizey, ksize, l_size
   integer          :: position_x, position_y
   integer, save    :: isizex_save=0, jsizex_save=0, ksize_save=0
   integer, save    :: isizey_save=0, jsizey_save=0
   integer, save    :: wbufferszx_save=0, ebufferszx_save=0, sbufferszx_save=0, nbufferszx_save=0
   integer, save    :: wbufferszy_save=0, ebufferszy_save=0, sbufferszy_save=0, nbufferszy_save=0
   integer, save    :: add_halo_save=0, update_flags_save=0, grid_offset_type_save
   integer, save    :: list=0

   add_halo = 0
   if(present(extra_halo)) add_halo = add_halo

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
   grid_offset_type = AGRID
   if( PRESENT(gridtype) ) grid_offset_type = gridtype

   is_complete = .true.
   if(PRESENT(complete)) then
      is_complete = complete
   end if
   tile = 1
   if(present(tile_count)) tile = tile_count
   if( tile > 1 ) then
      call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D_V: currently do not support multiple tile per pe')
   endif

   list = list+1
   if(list > MAX_DOMAIN_FIELDS)then
      write( text,'(i2)' ) MAX_DOMAIN_FIELDS
      call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D_V: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
   endif

   isizex = 0; jsizex = 0
   isizey = 0; jsizey = 0
   ksize = 0
   wbufferszx = 0; wbufferszy = 0
   ebufferszx = 0; ebufferszy = 0
   sbufferszx = 0; sbufferszy = 0
   nbufferszx = 0; nbufferszy = 0

   if(size(fieldx,3) .NE. size(fieldy,3)) call mpp_error(FATAL, &
        'MPP_UPDATE_NEST_FINE_3D_V: size(fieldx,3) .NE. size(fieldy,3)')
   ksize = size(fieldx,3)

   if(nest_domain%nest(nest_level)%is_coarse_pe) then
      f_addrsx(list) = LOC(fieldx)
      f_addrsy(list) = LOC(fieldy)
      isizex=size(fieldx,1); jsizex=size(fieldx,2)
      isizey=size(fieldy,1); jsizey=size(fieldy,2)
      if(size(fieldx,3) .NE. size(fieldy,3)) call mpp_error(FATAL, &
           'MPP_UPDATE_NEST_FINE_3D_V: size(fieldx,3) .NE. size(fieldy,3)')
   endif

   if(nest_domain%nest(nest_level)%is_fine_pe) then
      wb_addrsx(list) = LOC(wbufferx)
      eb_addrsx(list) = LOC(ebufferx)
      sb_addrsx(list) = LOC(sbufferx)
      nb_addrsx(list) = LOC(nbufferx)
      wb_addrsy(list) = LOC(wbuffery)
      eb_addrsy(list) = LOC(ebuffery)
      sb_addrsy(list) = LOC(sbuffery)
      nb_addrsy(list) = LOC(nbuffery)

      wbufferszx = size(wbufferx); ebufferszx = size(ebufferx)
      sbufferszx = size(sbufferx); nbufferszx = size(nbufferx)
      wbufferszy = size(wbuffery); ebufferszy = size(ebuffery)
      sbufferszy = size(sbuffery); nbufferszy = size(nbuffery)
   endif

   if(list == 1)then
      isizex_save = isizex; jsizex_save = jsizex
      isizey_save = isizey; jsizey_save = jsizey
      ksize_save = ksize
      grid_offset_type_save = grid_offset_type
      update_flags_save    = update_flags
      wbufferszx_save = wbufferszx; ebufferszx_save = ebufferszx
      sbufferszx_save = sbufferszx; nbufferszx_save = nbufferszx
      wbufferszy_save = wbufferszy; ebufferszy_save = ebufferszy
      sbufferszy_save = sbufferszy; nbufferszy_save = nbufferszy
      add_halo_save = add_halo
   else
      set_mismatch = .false.
      set_mismatch = set_mismatch .OR. (isizex_save /= isizex)
      set_mismatch = set_mismatch .OR. (jsizex_save /= jsizex)
      set_mismatch = set_mismatch .OR. (isizey_save /= isizey)
      set_mismatch = set_mismatch .OR. (jsizey_save /= jsizey)
      set_mismatch = set_mismatch .OR. (ksize_save /= ksize)
      set_mismatch = set_mismatch .OR. (grid_offset_type_save /= grid_offset_type)
      set_mismatch = set_mismatch .OR. (wbufferszx_save /= wbufferszx)
      set_mismatch = set_mismatch .OR. (ebufferszx_save /= ebufferszx)
      set_mismatch = set_mismatch .OR. (sbufferszx_save /= sbufferszx)
      set_mismatch = set_mismatch .OR. (nbufferszx_save /= nbufferszx)
      set_mismatch = set_mismatch .OR. (wbufferszy_save /= wbufferszy)
      set_mismatch = set_mismatch .OR. (ebufferszy_save /= ebufferszy)
      set_mismatch = set_mismatch .OR. (sbufferszy_save /= sbufferszy)
      set_mismatch = set_mismatch .OR. (nbufferszy_save /= nbufferszy)
      set_mismatch = set_mismatch .OR. (update_flags_save /= update_flags)
      set_mismatch = set_mismatch .OR. (add_halo_save /= add_halo)

      if(set_mismatch)then
         write( text,'(i2)' ) list
         call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D_V_: Incompatible field at count '//text//' for group update.' )
      endif
   endif

   if(is_complete) then
      l_size = list
      list = 0
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
         position_y = EAST
         position_x = NORTH
      case default
         call mpp_error(FATAL, "MPP_UPDATE_NEST_FINE_3D_V: invalid value of grid_offset_type")
      end select

      updatex => search_C2F_nest_overlap(nest_domain, nest_level, add_halo, position_x)
      updatey => search_C2F_nest_overlap(nest_domain, nest_level, add_halo, position_y)
      !---make sure data size match the size specified in updatex and updatey
      if(nest_domain%nest(nest_level)%is_coarse_pe) then
         if(isizex .NE. updatex%xend-updatex%xbegin+1 .OR. jsizex .NE. updatex%yend-updatex%ybegin+1) &
            call mpp_error(FATAL, "MPP_UPDATE_NEST_FINE_3D_V: mismatch between size of fieldx and updatex")
         if(isizey .NE. updatey%xend-updatey%xbegin+1 .OR. jsizey .NE. updatey%yend-updatey%ybegin+1) &
            call mpp_error(FATAL, "MPP_UPDATE_NEST_FINE_3D_V: mismatch between size of fieldy and updatey")
      endif

      if(nest_domain%nest(nest_level)%is_fine_pe) then
         call check_data_size("MPP_UPDATE_NEST_FINE_3D_V", "wbufferx", wbufferszx, "updatex", &
                              (updatex%west%ie_you-updatex%west%is_you+1)*(updatex%west%je_you-updatex%west%js_you+1)*ksize )
         call check_data_size("MPP_UPDATE_NEST_FINE_3D_V", "wbuffery", wbufferszy, "updatey", &
                              (updatey%west%ie_you-updatey%west%is_you+1)*(updatey%west%je_you-updatey%west%js_you+1)*ksize )
         call check_data_size("MPP_UPDATE_NEST_FINE_3D_V", "ebufferx", ebufferszx, "updatex", &
                              (updatex%east%ie_you-updatex%east%is_you+1)*(updatex%east%je_you-updatex%east%js_you+1)*ksize )
         call check_data_size("MPP_UPDATE_NEST_FINE_3D_V", "ebuffery", ebufferszy, "updatey", &
                              (updatey%east%ie_you-updatey%east%is_you+1)*(updatey%east%je_you-updatey%east%js_you+1)*ksize )
         call check_data_size("MPP_UPDATE_NEST_FINE_3D_V", "sbufferx", sbufferszx, "updatex", &
                              (updatex%south%ie_you-updatex%south%is_you+1)*(updatex%south%je_you-updatex%south%js_you+1)*ksize )
         call check_data_size("MPP_UPDATE_NEST_FINE_3D_V", "sbuffery", sbufferszy, "updatey", &
                              (updatey%south%ie_you-updatey%south%is_you+1)*(updatey%south%je_you-updatey%south%js_you+1)*ksize )
         call check_data_size("MPP_UPDATE_NEST_FINE_3D_V", "nbufferx", nbufferszx, "updatex", &
                              (updatex%north%ie_you-updatex%north%is_you+1)*(updatex%north%je_you-updatex%north%js_you+1)*ksize )
         call check_data_size("MPP_UPDATE_NEST_FINE_3D_V", "nbuffery", nbufferszy, "updatey", &
                              (updatey%north%ie_you-updatey%north%is_you+1)*(updatey%north%je_you-updatey%north%js_you+1)*ksize )
      endif

      call mpp_do_update_nest_fine(f_addrsx(1:l_size), f_addrsy(1:l_size), nest_domain%nest(nest_level), updatex, updatey, &
            d_type, ksize, wb_addrsx(1:l_size), wb_addrsy(1:l_size), eb_addrsx(1:l_size), eb_addrsy(1:l_size),  &
            sb_addrsx(1:l_size), sb_addrsy(1:l_size), nb_addrsx(1:l_size), nb_addrsy(1:l_size), update_flags )

   endif


end subroutine MPP_UPDATE_NEST_FINE_3D_V_

subroutine MPP_UPDATE_NEST_FINE_4D_V_(fieldx, fieldy, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, &
                                      ebufferx, ebuffery, nbufferx, nbuffery, nest_level, &
                                      flags, gridtype, complete, extra_halo, name, tile_count)
    MPP_TYPE_,             intent(in)      :: fieldx(:,:,:,:), fieldy(:,:,:,:) !< field x and y
                                                                     !! components on the model grid
    type(nest_domain_type), intent(inout)  :: nest_domain
    MPP_TYPE_,             intent(inout)   :: wbufferx(:,:,:,:), wbuffery(:,:,:,:) !< west side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: ebufferx(:,:,:,:), ebuffery(:,:,:,:) !< east side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: sbufferx(:,:,:,:), sbuffery(:,:,:,:) !< south side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    MPP_TYPE_,             intent(inout)   :: nbufferx(:,:,:,:), nbuffery(:,:,:,:) !< north side buffer x and y  components
                                                                 !! to be filled with data on coarse grid.
    integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
    integer,          intent(in), optional :: flags !< Specify the direction of fine grid halo buffer to be filled.
                                                    !! Default value is XUPDATE+YUPDATE.
    logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
                                                       !! Default value is .true.
    integer,          intent(in), optional :: gridtype
    integer,          intent(in), optional :: extra_halo !< extra halo for passing data
                                                         !! from coarse grid to fine grid.
                                                         !! Default is 0 and currently only support extra_halo = 0.
    character(len=*), intent(in), optional :: name !< Name of the nest domain.
    integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                         !! default is 1 and currently only support tile_count = 1.

    MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
    MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
    MPP_TYPE_ :: wbuffer3Dx(size(wbufferx,1),size(wbufferx,2),size(wbufferx,3)*size(wbufferx,4))
    MPP_TYPE_ :: ebuffer3Dx(size(ebufferx,1),size(ebufferx,2),size(ebufferx,3)*size(ebufferx,4))
    MPP_TYPE_ :: sbuffer3Dx(size(sbufferx,1),size(sbufferx,2),size(sbufferx,3)*size(sbufferx,4))
    MPP_TYPE_ :: nbuffer3Dx(size(nbufferx,1),size(nbufferx,2),size(nbufferx,3)*size(nbufferx,4))
    MPP_TYPE_ :: wbuffer3Dy(size(wbuffery,1),size(wbuffery,2),size(wbuffery,3)*size(wbuffery,4))
    MPP_TYPE_ :: ebuffer3Dy(size(ebuffery,1),size(ebuffery,2),size(ebuffery,3)*size(ebuffery,4))
    MPP_TYPE_ :: sbuffer3Dy(size(sbuffery,1),size(sbuffery,2),size(sbuffery,3)*size(sbuffery,4))
    MPP_TYPE_ :: nbuffer3Dy(size(nbuffery,1),size(nbuffery,2),size(nbuffery,3)*size(nbuffery,4))
    pointer( ptrx, field3Dx )
    pointer( ptry, field3Dy )
    pointer( ptr_wx, wbuffer3Dx)
    pointer( ptr_ex, ebuffer3Dx)
    pointer( ptr_sx, sbuffer3Dx)
    pointer( ptr_nx, nbuffer3Dx)
    pointer( ptr_wy, wbuffer3Dy)
    pointer( ptr_ey, ebuffer3Dy)
    pointer( ptr_sy, sbuffer3Dy)
    pointer( ptr_ny, nbuffer3Dy)
    ptrx = LOC(fieldx)
    ptry = LOC(fieldy)
    ptr_wx = LOC(wbufferx)
    ptr_ex = LOC(ebufferx)
    ptr_sx = LOC(sbufferx)
    ptr_nx = LOC(nbufferx)
    ptr_wy = LOC(wbuffery)
    ptr_ey = LOC(ebuffery)
    ptr_sy = LOC(sbuffery)
    ptr_ny = LOC(nbuffery)

    call MPP_UPDATE_NEST_FINE_3D_V_(field3Dx, field3Dy, nest_domain, wbuffer3Dx, wbuffer3Dy, sbuffer3Dx, sbuffer3Dy, &
                                    ebuffer3Dx, ebuffer3Dy, nbuffer3Dx, nbuffer3Dy, nest_level, &
                                    flags, gridtype, complete, extra_halo, name, tile_count)

end subroutine MPP_UPDATE_NEST_FINE_4D_V_

#endif

subroutine MPP_UPDATE_NEST_COARSE_2D_(field_in, nest_domain, field_out, nest_level, complete, position, name, tile_count)
      MPP_TYPE_,             intent(in)      :: field_in(:,:) !< field on the model grid
      type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                            !! between fine and coarse grid.
      MPP_TYPE_,             intent(inout)   :: field_out(:,:) !< field_out to be filled with data on coarse grid
      integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
      logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
                                                         !! Default value is .true.
      integer,          intent(in), optional :: position !< Cell position. Its value should be CENTER, EAST, CORNER
                                                         !! or NORTH. Default is CENTER.
      character(len=*), intent(in), optional :: name !< Name of the nest domain optional argument
      integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                           !! default is 1 and currently only support tile_count = 1.

      MPP_TYPE_ :: field3D_in(size(field_in,1),size(field_in,2),1)
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),1)
      pointer( ptr_in, field3D_in )
      pointer( ptr_out, field3D_out)
      ptr_in = LOC(field_in)
      ptr_out = LOC(field_out)
      call mpp_update_nest_coarse( field3D_in, nest_domain, field3D_out, nest_level, complete, position, name, tile_count)

      return


end subroutine MPP_UPDATE_NEST_COARSE_2D_


!--- field_in is the data on fine grid pelist to be passed to coarse grid pelist.
!--- field_in and field_out are all on the coarse grid. field_in is remapped from fine grid to coarse grid.
subroutine MPP_UPDATE_NEST_COARSE_3D_(field_in, nest_domain, field_out, nest_level, complete, position, name, tile_count)
   MPP_TYPE_,             intent(in)      :: field_in(:,:,:) !< field on the model grid
   type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                            !! between fine and coarse grid.
   MPP_TYPE_,             intent(inout)   :: field_out(:,:,:) !< field_out to be filled with data on coarse grid
   integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
   logical,          intent(in), optional :: complete !< When .true., do the buffer filling. Default value is .true.
   integer,          intent(in), optional :: position !< Cell position. Its value should be CENTER, EAST, CORNER,
                                                      !! or NORTH. Default is CENTER.
   character(len=*), intent(in), optional :: name !< Name of the nest domain optional argument
   integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                           !! default is 1 and currently only support tile_count = 1.

   MPP_TYPE_        :: d_type
   type(nestSpec), pointer :: update=>NULL()
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: fin_addrs=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: fout_addrs=-9999
   character(len=3) :: text
   logical          :: is_complete, set_mismatch
   integer          :: tile
   integer          :: update_position
   integer          :: isize_in, jsize_in, l_size
   integer          :: isize_out, jsize_out, ksize
   integer, save    :: isize_in_save=0, jsize_in_save=0, ksize_save=0
   integer, save    :: isize_out_save=0, jsize_out_save=0
   integer, save    :: update_position_save=0
   integer, save    :: list=0

   update_position = CENTER
   if(present(position)) update_position = position

   is_complete = .true.
   if(PRESENT(complete)) then
      is_complete = complete
   end if
   tile = 1
   if(present(tile_count)) tile = tile_count
   if( tile > 1 ) then
      call mpp_error(FATAL,'MPP_UPDATE_NEST_COARSE_3D: currently do not support multiple tile per pe')
   endif

   list = list+1
   if(list > MAX_DOMAIN_FIELDS)then
      write( text,'(i2)' ) MAX_DOMAIN_FIELDS
      call mpp_error(FATAL,'MPP_UPDATE_NEST_COARSE_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
   endif

   isize_in = 0; jsize_in = 0
   isize_out = 0; jsize_out = 0

   if(nest_domain%nest(nest_level)%is_fine_pe) then
      fin_addrs(list) = LOC(field_in)
      isize_in=size(field_in,1); jsize_in=size(field_in,2)
      ksize = size(field_in,3)
   endif
   if(nest_domain%nest(nest_level)%is_coarse_pe) then
      fout_addrs(list) = LOC(field_out)
      isize_out=size(field_out,1); jsize_out=size(field_out,2)
      ksize = size(field_out,3)
   endif

   if(list == 1)then
      isize_in_save = isize_in; jsize_in_save = jsize_in; ksize_save = ksize
      isize_out_save = isize_out; jsize_out_save = jsize_out
      update_position_save = update_position
   else
      set_mismatch = .false.
      set_mismatch = set_mismatch .OR. (isize_in_save /= isize_in)
      set_mismatch = set_mismatch .OR. (jsize_in_save /= jsize_in)
      set_mismatch = set_mismatch .OR. (ksize_save /= ksize)
      set_mismatch = set_mismatch .OR. (isize_out_save /= isize_out)
      set_mismatch = set_mismatch .OR. (jsize_out_save /= jsize_out)
      set_mismatch = set_mismatch .OR. (update_position_save /= update_position)

      if(set_mismatch)then
         write( text,'(i2)' ) list
         call mpp_error(FATAL,'MPP_UPDATE_NEST_COARSE_3D_: Incompatible field at count '//text//' for group update.' )
      endif
   endif

   if(is_complete) then
      l_size = list
      list = 0
   end if

   if(is_complete)then
      update => search_F2C_nest_overlap(nest_domain, nest_level, update_position)
      if(nest_domain%nest(nest_level)%is_fine_pe) then
         if(isize_in .NE. update%xsize_c .OR. jsize_in .NE. update%ysize_c) then
           print*,"isize_in,jsize_in=", isize_in, jsize_in, update%xsize_c, update%ysize_c, mpp_pe()
           call mpp_error(FATAL, "MPP_UPDATE_NEST_COARSE_3D_: "// &
                "size of field_in does not match size specified in nest_domain_type")
         endif
      endif
      if(nest_domain%nest(nest_level)%is_coarse_pe) then
         if(isize_out .NE. update%xsize_c .OR. jsize_out .NE. update%ysize_c) then
           print*,"isize_out,jsize_out=", isize_out, jsize_out, update%xsize_c, update%ysize_c, mpp_pe()
           call mpp_error(FATAL, "MPP_UPDATE_NEST_COARSE_3D_: " // &
                "size of field_out does not match size specified in nest_domain_type")
         endif
      endif

      call mpp_do_update_nest_coarse(fin_addrs(1:l_size), fout_addrs(1:l_size), nest_domain, update, d_type, ksize)
   endif

end subroutine MPP_UPDATE_NEST_COARSE_3D_

!###############################################################################
subroutine MPP_UPDATE_NEST_COARSE_4D_(field_in, nest_domain, field_out, nest_level, complete, position, name, tile_count)
      MPP_TYPE_,             intent(in)      :: field_in(:,:,:,:) !< field on the model grid
      type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                            !! between fine and coarse grid.
      MPP_TYPE_,             intent(inout)   :: field_out(:,:,:,:) !< field_out to be filled with data on coarse grid
      integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
      logical,          intent(in), optional :: complete !< When .true., do the buffer filling. Default value is .true.
      integer,          intent(in), optional :: position !< Cell position. Its value should be CENTER, EAST, CORNER,
                                                      !! or NORTH. Default is CENTER.
      character(len=*), intent(in), optional :: name !< Name of the nest domain optional argument
      integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                           !! default is 1 and currently only support tile_count = 1.

      MPP_TYPE_ :: field3D_in(size(field_in,1),size(field_in,2),size(field_in,3)*size(field_in,4))
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),size(field_out,3)*size(field_out,4))
      pointer( ptr_in, field3D_in )
      pointer( ptr_out, field3D_out )
      ptr_in = LOC(field_in)
      ptr_out = LOC(field_out)
      call mpp_update_nest_coarse( field3D_in, nest_domain, field3D_out, nest_level, complete, position, name, tile_count)

      return


end subroutine MPP_UPDATE_NEST_COARSE_4D_


#ifdef VECTOR_FIELD_
!--- field_in is the data on fine grid pelist to be passed to coarse grid pelist.
!--- field_in and field_out are all on the coarse grid. field_in is remapped from fine grid to coarse grid.
subroutine MPP_UPDATE_NEST_COARSE_2D_V_(fieldx_in, fieldy_in, nest_domain, fieldx_out, fieldy_out, nest_level, &
                                        flags, gridtype, complete, name, tile_count)
   MPP_TYPE_,             intent(in)      :: fieldx_in(:,:) !< x component of field on the model grid
   MPP_TYPE_,             intent(in)      :: fieldy_in(:,:) !< y component of field on the model grid
   type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                         !! between fine and coarse grid.
   integer,          intent(in), optional :: flags, gridtype !< Specify the direction of fine grid halo buffer to be filled.
                                                    !! Default value is XUPDATE+YUPDATE.
   MPP_TYPE_,             intent(inout)   :: fieldx_out(:,:) !< x component of field_out to be
                                                             !! filled with data on coarse grid
   MPP_TYPE_,             intent(inout)   :: fieldy_out(:,:) !< y component of field_out to be
                                                             !! filled with data on coarse grid
   integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
   logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
   character(len=*), intent(in), optional :: name !< Name of the nest domain optional argument
   integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                        !! default is 1 and currently only support tile_count = 1.

   MPP_TYPE_ :: field3Dx_in(size(fieldx_in,1),size(fieldx_in,2),1)
   MPP_TYPE_ :: field3Dy_in(size(fieldy_in,1),size(fieldy_in,2),1)
   MPP_TYPE_ :: field3Dx_out(size(fieldx_out,1),size(fieldx_out,2),1)
   MPP_TYPE_ :: field3Dy_out(size(fieldy_out,1),size(fieldy_out,2),1)
   pointer( ptrx_in, field3Dx_in )
   pointer( ptry_in, field3Dy_in )
   pointer( ptrx_out, field3Dx_out )
   pointer( ptry_out, field3Dy_out )

   ptrx_in = LOC(fieldx_in)
   ptry_in = LOC(fieldy_in)
   ptrx_out = LOC(fieldx_out)
   ptry_out = LOC(fieldy_out)

   call MPP_UPDATE_NEST_COARSE_3D_V_(field3Dx_in, field3Dy_in, nest_domain, field3Dx_out, field3Dy_out, &
                                     nest_level, flags, gridtype, complete, name, tile_count)

end subroutine MPP_UPDATE_NEST_COARSE_2D_V_


!--- field_in is the data on fine grid pelist to be passed to coarse grid pelist.
!--- field_in and field_out are all on the coarse grid. field_in is remapped from fine grid to coarse grid.
subroutine MPP_UPDATE_NEST_COARSE_3D_V_(fieldx_in, fieldy_in, nest_domain, fieldx_out, fieldy_out, nest_level, &
                                        flags, gridtype, complete, name, tile_count)
   MPP_TYPE_,             intent(in)      :: fieldx_in(:,:,:) !< x component field on the model grid
   MPP_TYPE_,             intent(in)      :: fieldy_in(:,:,:) !< y component of field on the model grid
   type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                         !! between fine and coarse grid.
   integer,          intent(in), optional :: flags, gridtype !< Specify the direction of fine grid halo buffer to be filled.
                                                    !! Default value is XUPDATE+YUPDATE.
   MPP_TYPE_,             intent(inout)   :: fieldx_out(:,:,:) !< x component of field_out to be
                                                               !! filled with data on coarse grid
   MPP_TYPE_,             intent(inout)   :: fieldy_out(:,:,:) !< y component of field_out to be
                                                               !! filled with data on coarse grid
   integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
   logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
   character(len=*), intent(in), optional :: name !< Name of the nest domain optional argument
   integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                        !! default is 1 and currently only support tile_count = 1.

   MPP_TYPE_        :: d_type
   type(nestSpec), pointer :: updatex=>NULL()
   type(nestSpec), pointer :: updatey=>NULL()
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: fin_addrsx=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: fin_addrsy=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: fout_addrsx=-9999
   integer(i8_kind),dimension(MAX_DOMAIN_FIELDS),save :: fout_addrsy=-9999
   character(len=3) :: text
   logical          :: is_complete, set_mismatch
   integer          :: tile
   integer          :: grid_offset_type, update_flags
   integer          :: isizex_in, jsizex_in, l_size
   integer          :: isizex_out, jsizex_out, ksize
   integer          :: isizey_in, jsizey_in
   integer          :: isizey_out, jsizey_out
   integer          :: position_x, position_y
   integer, save    :: isizex_in_save=0, jsizex_in_save=0, ksize_save=0
   integer, save    :: isizex_out_save=0, jsizex_out_save=0
   integer, save    :: isizey_in_save=0, jsizey_in_save=0
   integer, save    :: isizey_out_save=0, jsizey_out_save=0
   integer, save    :: grid_offset_type_save=0
   integer, save    :: update_flags_save=0
   integer, save    :: list=0

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

   is_complete = .true.
   if(PRESENT(complete)) then
      is_complete = complete
   end if
   tile = 1
   if(present(tile_count)) tile = tile_count
   if( tile > 1 ) then
      call mpp_error(FATAL,'MPP_UPDATE_NEST_COARSE_3D_V: currently do not support multiple tile per pe')
   endif

   list = list+1
   if(list > MAX_DOMAIN_FIELDS)then
      write( text,'(i2)' ) MAX_DOMAIN_FIELDS
      call mpp_error(FATAL,'MPP_UPDATE_NEST_COARSE_3D_V: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
   endif

   isizex_in = 0; jsizex_in = 0
   isizex_out = 0; jsizex_out = 0
   isizey_in = 0; jsizey_in = 0
   isizey_out = 0; jsizey_out = 0
   ksize = 0

   if(nest_domain%nest(nest_level)%is_fine_pe) then
      fin_addrsx(list) = LOC(fieldx_in)
      fin_addrsy(list) = LOC(fieldy_in)
      isizex_in=size(fieldx_in,1); jsizex_in=size(fieldx_in,2)
      isizey_in=size(fieldy_in,1); jsizey_in=size(fieldy_in,2)
      ksize = size(fieldx_in,3)
      if(size(fieldx_in,3) .NE. size(fieldy_in,3)) call mpp_error(FATAL, &
           'MPP_UPDATE_NEST_COARSE_3D_V: size(fieldx_in,3) .NE. size(fieldy_in,3)')
   endif
   if(nest_domain%nest(nest_level)%is_coarse_pe) then
      fout_addrsx(list) = LOC(fieldx_out)
      fout_addrsy(list) = LOC(fieldy_out)
      isizex_out=size(fieldx_out,1); jsizex_out=size(fieldx_out,2)
      isizey_out=size(fieldy_out,1); jsizey_out=size(fieldy_out,2)
      ksize = size(fieldx_out,3)
      if(size(fieldx_out,3) .NE. size(fieldy_out,3)) call mpp_error(FATAL, &
           'MPP_UPDATE_NEST_COARSE_3D_V: size(fieldx_out,3) .NE. size(fieldy_out,3)')
   endif

   if(list == 1)then
      isizex_in_save = isizex_in; jsizex_in_save = jsizex_in; ksize_save = ksize
      isizex_out_save = isizex_out; jsizex_out_save = jsizex_out
      isizey_in_save = isizey_in; jsizey_in_save = jsizey_in
      isizey_out_save = isizey_out; jsizey_out_save = jsizey_out
      update_flags_save = update_flags
      grid_offset_type_save = grid_offset_type
   else
      set_mismatch = .false.
      set_mismatch = set_mismatch .OR. (isizex_in_save /= isizex_in)
      set_mismatch = set_mismatch .OR. (jsizex_in_save /= jsizex_in)
      set_mismatch = set_mismatch .OR. (isizey_in_save /= isizey_in)
      set_mismatch = set_mismatch .OR. (jsizey_in_save /= jsizey_in)
      set_mismatch = set_mismatch .OR. (ksize_save /= ksize)
      set_mismatch = set_mismatch .OR. (isizex_out_save /= isizex_out)
      set_mismatch = set_mismatch .OR. (jsizex_out_save /= jsizex_out)
      set_mismatch = set_mismatch .OR. (isizey_out_save /= isizey_out)
      set_mismatch = set_mismatch .OR. (jsizey_out_save /= jsizey_out)
      set_mismatch = set_mismatch .OR. (update_flags_save /= update_flags)
      set_mismatch = set_mismatch .OR. (grid_offset_type_save /= grid_offset_type)

      if(set_mismatch)then
         write( text,'(i2)' ) list
         call mpp_error(FATAL,'MPP_UPDATE_NEST_COARSE_3D_V: Incompatible field at count '//text//' for group update.' )
      endif
   endif

   if(is_complete) then
      l_size = list
      list = 0
   end if

   if(is_complete)then
      select case(grid_offset_type)
      case (AGRID)
         position_x = CENTER
         position_y = CENTER
      case (BGRID_NE, BGRID_SW)
         call mpp_error(FATAL, "MPP_UPDATE_NEST_COARSE_3D_V: currently does not support BGRID, contact developer")
         position_x = CORNER
         position_y = CORNER
      case (CGRID_NE, CGRID_SW)
         position_x = EAST
         position_y = NORTH
      case (DGRID_NE, DGRID_SW)
         position_y = EAST
         position_x = NORTH
      case default
         call mpp_error(FATAL, "MPP_UPDATE_NEST_COARSE_3D_V: invalid value of grid_offset_type")
      end select

      updatex => search_F2C_nest_overlap(nest_domain, nest_level, position_x)
      updatey => search_F2C_nest_overlap(nest_domain, nest_level, position_y)

      if(nest_domain%nest(nest_level)%is_fine_pe) then
         if(isizex_in .NE. updatex%xsize_c .OR. jsizex_in .NE. updatex%ysize_c) then
           print*,"isizex_in,jsizex_in=", isizex_in, jsizex_in, updatex%xsize_c, updatex%ysize_c, mpp_pe()
           call mpp_error(FATAL, "MPP_UPDATE_NEST_COARSE_3D_V: "// &
                "size of fieldx_in does not match size specified in nest_domain_type")
         endif
         if(isizey_in .NE. updatey%xsize_c .OR. jsizey_in .NE. updatey%ysize_c) then
           print*,"isizey_in,jsizey_in=", isizey_in, jsizey_in, updatey%xsize_c, updatey%ysize_c, mpp_pe()
           call mpp_error(FATAL, "MPP_UPDATE_NEST_COARSE_3D_V: "// &
                "size of fieldy_in does not match size specified in nest_domain_type")
         endif
      endif
      if(nest_domain%nest(nest_level)%is_coarse_pe) then
         if(isizex_out .NE. updatex%xsize_c .OR. jsizex_out .NE. updatex%ysize_c) then
           print*,"isizex_out,jsizex_out=", isizex_out, jsizex_out, updatex%xsize_c, updatex%ysize_c, mpp_pe()
           call mpp_error(FATAL, "MPP_UPDATE_NEST_COARSE_3D_V: " // &
                "size of fieldx_out does not match size specified in nest_domain_type")
         endif
         if(isizey_out .NE. updatey%xsize_c .OR. jsizey_out .NE. updatey%ysize_c) then
           print*,"isizey_out,jsizey_out=", isizey_out, jsizey_out, updatey%xsize_c, updatey%ysize_c, mpp_pe()
           call mpp_error(FATAL, "MPP_UPDATE_NEST_COARSE_3D_V: " // &
                "size of fieldy_out does not match size specified in nest_domain_type")
         endif
      endif

      call mpp_do_update_nest_coarse(fin_addrsx(1:l_size), fin_addrsy(1:l_size), fout_addrsx(1:l_size), fout_addrsy(1:l_size), &
                                     nest_domain, nest_domain%nest(nest_level), updatex, updatey, d_type, ksize, update_flags)
   endif

end subroutine MPP_UPDATE_NEST_COARSE_3D_V_

subroutine MPP_UPDATE_NEST_COARSE_4D_V_(fieldx_in, fieldy_in, nest_domain, fieldx_out, fieldy_out, nest_level, &
                                        flags, gridtype, complete, name, tile_count)
   MPP_TYPE_,             intent(in)      :: fieldx_in(:,:,:,:) !< x component field on the model grid
   MPP_TYPE_,             intent(in)      :: fieldy_in(:,:,:,:) !< y component field on the model grid
   type(nest_domain_type), intent(inout)  :: nest_domain !< Holds the information to pass data
                                                         !! between fine and coarse grid.
   integer,          intent(in), optional :: flags, gridtype !< Specify the direction of fine grid halo buffer to be filled.
                                                    !! Default value is XUPDATE+YUPDATE.
   MPP_TYPE_,             intent(inout)   :: fieldx_out(:,:,:,:) !< x component of field_out to be
                                                                 !! filled with data on coarse grid
   MPP_TYPE_,             intent(inout)   :: fieldy_out(:,:,:,:) !< y component of field_out to be
                                                                 !! filled with data on coarse grid
   integer,          intent(in)           :: nest_level !< level of the nest (> 1 implies a telescoping nest)
   logical,          intent(in), optional :: complete !< When .true., do the buffer filling.
   character(len=*), intent(in), optional :: name !< Name of the nest domain optional argument
   integer,          intent(in), optional :: tile_count !< Used to support multiple-tile-per-pe.
                                                        !! default is 1 and currently only support tile_count = 1.

   MPP_TYPE_ :: field3Dx_in(size(fieldx_in,1),size(fieldx_in,2),size(fieldx_in,3)*size(fieldx_in,4))
   MPP_TYPE_ :: field3Dy_in(size(fieldy_in,1),size(fieldy_in,2),size(fieldy_in,3)*size(fieldy_in,4))
   MPP_TYPE_ :: field3Dx_out(size(fieldx_out,1),size(fieldx_out,2),size(fieldx_out,3)*size(fieldx_out,4))
   MPP_TYPE_ :: field3Dy_out(size(fieldy_out,1),size(fieldy_out,2),size(fieldy_out,3)*size(fieldy_out,4))
   pointer( ptrx_in, field3Dx_in )
   pointer( ptry_in, field3Dy_in )
   pointer( ptrx_out, field3Dx_out )
   pointer( ptry_out, field3Dy_out )

   ptrx_in = LOC(fieldx_in)
   ptry_in = LOC(fieldy_in)
   ptrx_out = LOC(fieldx_out)
   ptry_out = LOC(fieldy_out)

   call MPP_UPDATE_NEST_COARSE_3D_V_(field3Dx_in, field3Dy_in, nest_domain, field3Dx_out, field3Dy_out, &
                                     nest_level, flags, gridtype, complete, name, tile_count)

end subroutine MPP_UPDATE_NEST_COARSE_4D_V_


#endif
