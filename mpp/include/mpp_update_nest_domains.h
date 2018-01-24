! -*-f90-*-
subroutine MPP_UPDATE_NEST_FINE_2D_(field, nest_domain, wbuffer, ebuffer, sbuffer, nbuffer, &
                                    flags, complete, position, extra_halo, name, tile_count) 
      MPP_TYPE_,             intent(in)      :: field(:,:)
      type(nest_domain_type), intent(inout)  :: nest_domain
      MPP_TYPE_,             intent(inout)   :: wbuffer(:,:)
      MPP_TYPE_,             intent(inout)   :: ebuffer(:,:)
      MPP_TYPE_,             intent(inout)   :: sbuffer(:,:)
      MPP_TYPE_,             intent(inout)   :: nbuffer(:,:)
      integer,          intent(in), optional :: flags
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      integer,          intent(in), optional :: extra_halo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

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
                                 flags, complete, position, extra_halo, name, tile_count) 

      return


end subroutine MPP_UPDATE_NEST_FINE_2D_

subroutine MPP_UPDATE_NEST_FINE_3D_(field, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, &
                                    flags, complete, position, extra_halo, name, tile_count) 
    MPP_TYPE_,             intent(in)      :: field(:,:,:)
    type(nest_domain_type), intent(inout)  :: nest_domain
    MPP_TYPE_,             intent(inout)   :: wbuffer(:,:,:)
    MPP_TYPE_,             intent(inout)   :: ebuffer(:,:,:)
    MPP_TYPE_,             intent(inout)   :: sbuffer(:,:,:)
    MPP_TYPE_,             intent(inout)   :: nbuffer(:,:,:)
    integer,          intent(in), optional :: flags
    logical,          intent(in), optional :: complete
    integer,          intent(in), optional :: position
    integer,          intent(in), optional :: extra_halo
    character(len=*), intent(in), optional :: name
    integer,          intent(in), optional :: tile_count

   MPP_TYPE_        :: d_type
   type(nestSpec), pointer :: update=>NULL()
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: f_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: wb_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: eb_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: sb_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: nb_addrs=-9999
   character(len=3) :: text
   logical          :: is_complete, set_mismatch
   integer          :: tile
   integer          :: add_halo, update_flags, update_position
   integer          :: wbuffersz, ebuffersz, sbuffersz, nbuffersz
   integer          :: isize, jsize, ksize, l_size
   integer, save    :: isize_save=0, jsize_save=0, ksize_save=0
   integer          :: wbuffersz_save=0, ebuffersz_save=0, sbuffersz_save=0, nbuffersz_save=0
   integer, save    :: add_halo_save=0, update_flags_save=0, update_position_save=0
   integer, save    :: list=0 

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
      update => search_C2F_nest_overlap(nest_domain, add_halo, update_position)
      call mpp_do_update_nest_fine(f_addrs(1:l_size), nest_domain, update, d_type, ksize, &
            wb_addrs(1:l_size), eb_addrs(1:l_size), sb_addrs(1:l_size), nb_addrs(1:l_size), update_flags )

   endif


end subroutine MPP_UPDATE_NEST_FINE_3D_


!###############################################################################
subroutine MPP_UPDATE_NEST_FINE_4D_(field, nest_domain, wbuffer, ebuffer, sbuffer, nbuffer, &
                                    flags, complete, position, extra_halo, name, tile_count) 
      MPP_TYPE_,             intent(in)      :: field(:,:,:,:)
      type(nest_domain_type), intent(inout)  :: nest_domain
      MPP_TYPE_,             intent(inout)   :: wbuffer(:,:,:,:)
      MPP_TYPE_,             intent(inout)   :: ebuffer(:,:,:,:)
      MPP_TYPE_,             intent(inout)   :: sbuffer(:,:,:,:)
      MPP_TYPE_,             intent(inout)   :: nbuffer(:,:,:,:)
      integer,          intent(in), optional :: flags
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      integer,          intent(in), optional :: extra_halo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

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
                                 flags, complete, position, extra_halo, name, tile_count) 

      return


end subroutine MPP_UPDATE_NEST_FINE_4D_



subroutine MPP_UPDATE_NEST_COARSE_2D_(field_in, nest_domain, field_out, complete, position, name, tile_count) 
      MPP_TYPE_,             intent(in)      :: field_in(:,:)
      type(nest_domain_type), intent(inout)  :: nest_domain
      MPP_TYPE_,             intent(inout)   :: field_out(:,:)
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3D_in(size(field_in,1),size(field_in,2),1)
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),1) 
      pointer( ptr_in, field3D_in )
      pointer( ptr_out, field3D_out)
      ptr_in = LOC(field_in)
      ptr_out = LOC(field_out)
      call mpp_update_nest_coarse( field3D_in, nest_domain, field3D_out, complete, position, name, tile_count) 

      return


end subroutine MPP_UPDATE_NEST_COARSE_2D_


!--- field_in is the data on fine grid pelist to be passed to coarse grid pelist.
!--- field_in and field_out are all on the coarse grid. field_in is remapped from fine grid to coarse grid.
subroutine MPP_UPDATE_NEST_COARSE_3D_(field_in, nest_domain, field_out, complete, position, name, tile_count) 
   MPP_TYPE_,             intent(in)      :: field_in(:,:,:)
   type(nest_domain_type), intent(inout)  :: nest_domain
   MPP_TYPE_,             intent(inout)   :: field_out(:,:,:)
   logical,          intent(in), optional :: complete
   integer,          intent(in), optional :: position
   character(len=*), intent(in), optional :: name
   integer,          intent(in), optional :: tile_count

   MPP_TYPE_        :: d_type
   type(nestSpec), pointer :: update=>NULL()
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: fin_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: fout_addrs=-9999
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

   if(nest_domain%is_fine_pe) then
      fin_addrs(list) = LOC(field_in)
      isize_in=size(field_in,1); jsize_in=size(field_in,2)
      ksize = size(field_in,3)
   endif
   if(nest_domain%is_coarse_pe) then
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
      update => search_F2C_nest_overlap(nest_domain, update_position)
      if(nest_domain%is_fine_pe) then
         if(isize_in .NE. update%xsize_c .OR. jsize_in .NE. update%ysize_c) then
           print*,"isize_in,jsize_in=", isize_in, jsize_in, update%xsize_c, update%ysize_c, mpp_pe()
           call mpp_error(FATAL, "MPP_UPDATE_NEST_COARSE_3D_: "// &
                "size of field_in does not match size specified in nest_domain_type")
         endif
      endif
      if(nest_domain%is_coarse_pe) then
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
subroutine MPP_UPDATE_NEST_COARSE_4D_(field_in, nest_domain, field_out, complete, position, name, tile_count) 
      MPP_TYPE_,             intent(in)      :: field_in(:,:,:,:)
      type(nest_domain_type), intent(inout)  :: nest_domain
      MPP_TYPE_,             intent(inout)   :: field_out(:,:,:,:)
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3D_in(size(field_in,1),size(field_in,2),size(field_in,3)*size(field_in,4))
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),size(field_out,3)*size(field_out,4))
      pointer( ptr_in, field3D_in )
      pointer( ptr_out, field3D_out )
      ptr_in = LOC(field_in)
      ptr_out = LOC(field_out)
      call mpp_update_nest_coarse( field3D_in, nest_domain, field3D_out, complete, position, name, tile_count) 

      return


end subroutine MPP_UPDATE_NEST_COARSE_4D_
