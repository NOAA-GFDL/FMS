    subroutine MPP_UPDATE_DOMAINS_2D_( field, domain, flags, complete, free, list_size, &
                                       dc_handle, position)
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:)
      type(domain2D), intent(inout) :: domain  ! Must be definable in mpp_update_init_comm
      integer, intent(in), optional :: flags
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      integer, intent(in), optional :: position

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( ptr, field3D )
      ptr = LOC(field)
      call mpp_update_domains( field3D, domain, flags, complete, free, list_size, &
                               dc_handle, position )
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      field3D = RESHAPE( field, SHAPE(field3D) )
      call mpp_update_domains( field3D, domain, flags, position=position)
      field = RESHAPE( field3D, SHAPE(field) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_2D_

    subroutine MPP_UPDATE_DOMAINS_3D_( field, domain, flags, complete, free, list_size, dc_handle, &
                                       position )
!updates data domain of 3D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:)
      type(domain2D), intent(inout) :: domain  ! Must be definable in mpp_update_init_comm
      integer, intent(in), optional :: flags
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      integer, intent(in), optional :: position

#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      type(DomainCommunicator2D),pointer,save :: d_comm =>NULL()
      logical                       :: do_update, free_comm
      integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: f_addrs=-9999
      integer, save :: isize(1)=0,jsize(1)=0,ke=0,l_size=0
      integer, save :: pos
      integer       :: update_position
      logical       :: set_mismatch
      character(len=2) :: text
      MPP_TYPE_ :: d_type

      if(PRESENT(complete) .or. PRESENT(free))then
         do_update=.true.; if(PRESENT(complete))do_update=complete
         free_comm=.false.; if(PRESENT(free))free_comm=free
         if(free_comm)then
            f_addrs(1) = LOC(field)
            ke = size(field,3)
            l_size=1; if(PRESENT(list_size))l_size=list_size
            call mpp_update_free_comm(domain,f_addrs(1),ke,l_size,flags=flags)
            l_size=0; f_addrs=-9999
         else
            l_size = l_size+1
            if(l_size > MAX_DOMAIN_FIELDS)then
               write( text,'(i2)' ) MAX_DOMAIN_FIELDS
               call mpp_error(FATAL,'MPP_UPDATE_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
            endif
            update_position = CENTER
            if(present(position)) update_position = position
            f_addrs(l_size) = LOC(field)
            if(l_size == 1)then
               isize(1)=size(field,1); jsize(1)=size(field,2); ke = size(field,3); pos = update_position
            else
               set_mismatch = .false.
               set_mismatch = set_mismatch .OR. (isize(1) /= size(field,1))
               set_mismatch = set_mismatch .OR. (jsize(1) /= size(field,2))
               set_mismatch = set_mismatch .OR. (ke /= size(field,3))
               set_mismatch = set_mismatch .OR. (update_position /= pos)
               if(set_mismatch)then
                  write( text,'(i2)' ) l_size
                  call mpp_error(FATAL,'MPP_UPDATE_3D: Incompatible field at count '//text//' for group update.' )
               endif
            endif
            if(do_update)then
               if(PRESENT(dc_handle))d_comm =>dc_handle  ! User has kept pointer to d_comm
               if(.not.ASSOCIATED(d_comm))then  ! d_comm needs initialization or lookup
                  d_comm =>mpp_update_init_comm(domain,f_addrs(1:l_size), &
                       isize,jsize,ke, flags=flags, position=update_position )
                  if(PRESENT(dc_handle))dc_handle =>d_comm  ! User wants to keep pointer to d_comm
               endif
               call mpp_do_update( f_addrs(1:l_size), d_comm, d_type )
               d_comm =>NULL(); l_size=0; f_addrs=-9999; isize=0;  jsize=0;  ke=0
            endif
         endif
      else
         call mpp_do_update( field, domain, flags, position=position )
      endif
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      call mpp_do_update( field, domain, flags, position = position )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_3D_

    subroutine MPP_UPDATE_DOMAINS_4D_( field, domain, flags, complete, free, list_size, &
                                       dc_handle, position )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:)
      type(domain2D), intent(inout) :: domain  ! Must be definable in mpp_update_init_comm
      integer, intent(in), optional :: flags
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      integer, intent(in), optional :: position

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( ptr, field3D )
      ptr = LOC(field)
      call mpp_update_domains( field3D, domain, flags, complete, free, list_size, dc_handle, position)
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      field3D = RESHAPE( field, SHAPE(field3D) )
      call mpp_update_domains( field3D, domain, flags, position=position)
      field = RESHAPE( field3D, SHAPE(field) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_

    subroutine MPP_UPDATE_DOMAINS_5D_( field, domain, flags, complete, free, list_size, &
                                       dc_handle, position )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:,:)
      type(domain2D), intent(inout) :: domain  ! Must be definable in mpp_update_init_comm
      integer, intent(in), optional :: flags
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      integer, intent(in), optional :: position

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( ptr, field3D )
      ptr = LOC(field)
      call mpp_update_domains( field3D, domain, flags, complete, free, list_size, dc_handle, position )
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      field3D = RESHAPE( field, SHAPE(field3D) )
      call mpp_update_domains( field3D, domain, flags, position = position)
      field = RESHAPE( field3D, SHAPE(field) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_

    subroutine MPP_REDISTRIBUTE_2D_( domain_in, field_in, domain_out, field_out, complete, free, list_size, dc_handle, position )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:)
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      integer, intent(in), optional :: position
      MPP_TYPE_ :: field3D_in (size(field_in, 1),size(field_in, 2),1)
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),1)
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( ptr_in,  field3D_in  )
      pointer( ptr_out, field3D_out )
      ptr_in  = LOC(field_in )
      ptr_out = LOC(field_out)
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out, complete, free, list_size, dc_handle, position )
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      field3D_in = RESHAPE( field_in, SHAPE(field3D_in) )
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out, position=position )
      field_out = RESHAPE( field3D_out, SHAPE(field_out) )
#endif
      return
    end subroutine MPP_REDISTRIBUTE_2D_


    subroutine MPP_REDISTRIBUTE_3D_( domain_in, field_in, domain_out, field_out, complete, free, list_size, dc_handle, position )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:,:)
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      integer, intent(in), optional :: position

#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      type(DomainCommunicator2D),pointer,save :: d_comm =>NULL()
      logical                       :: do_redist,free_comm
      integer                       :: lsize
      integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: l_addrs_in=-9999, l_addrs_out=-9999
      integer, save :: isize_in=0,jsize_in=0,ke_in=0,l_size=0
      integer, save :: isize_out=0,jsize_out=0,ke_out=0
      logical       :: set_mismatch
      integer       :: ke
      character(len=2) :: text
      MPP_TYPE_ :: d_type

      if(present(position)) then
         if(position .NE. CENTER) call mpp_error( FATAL,  &
             'MPP_REDISTRIBUTE_3Dold_: only position = CENTER is implemented, contact author')
      endif

      if(PRESENT(complete) .or. PRESENT(free))then
        do_redist=.true.; if(PRESENT(complete))do_redist=complete
        free_comm=.false.; if(PRESENT(free))free_comm=free
        if(free_comm)then
          l_addrs_in(1) = LOC(field_in); l_addrs_out(1) = LOC(field_out)
          if(l_addrs_out(1)>0)then
          ke = size(field_out,3)
          else
            ke = size(field_in,3)
      end if
          lsize=1; if(PRESENT(list_size))lsize=list_size
          call mpp_redistribute_free_comm(domain_in,l_addrs_in(1),domain_out,l_addrs_out(1),ke,lsize)
        else
          l_size = l_size+1
          if(l_size > MAX_DOMAIN_FIELDS)then
            write( text,'(i2)' ) MAX_DOMAIN_FIELDS
            call mpp_error(FATAL,'MPP_REDISTRIBUTE_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group redistribute.' )
      end if
          l_addrs_in(l_size) = LOC(field_in); l_addrs_out(l_size) = LOC(field_out)
          if(l_size == 1)then
            if(l_addrs_in(l_size) > 0)then
              isize_in=size(field_in,1); jsize_in=size(field_in,2); ke_in = size(field_in,3)
         end if
            if(l_addrs_out(l_size) > 0)then
              isize_out=size(field_out,1); jsize_out=size(field_out,2); ke_out = size(field_out,3)
            endif
          else   
            set_mismatch = .false.
            set_mismatch = l_addrs_in(l_size) == 0 .AND. l_addrs_in(l_size-1) /= 0
            set_mismatch = set_mismatch .OR. (l_addrs_in(l_size) > 0 .AND. l_addrs_in(l_size-1) == 0)
            set_mismatch = set_mismatch .OR. (l_addrs_out(l_size) == 0 .AND. l_addrs_out(l_size-1) /= 0)
            set_mismatch = set_mismatch .OR. (l_addrs_out(l_size) > 0 .AND. l_addrs_out(l_size-1) == 0)
           if(l_addrs_in(l_size) > 0)then
              set_mismatch = set_mismatch .OR. (isize_in /= size(field_in,1))
              set_mismatch = set_mismatch .OR. (jsize_in /= size(field_in,2))
              set_mismatch = set_mismatch .OR. (ke_in /= size(field_in,3))
            endif
            if(l_addrs_out(l_size) > 0)then
              set_mismatch = set_mismatch .OR. (isize_out /= size(field_out,1))
              set_mismatch = set_mismatch .OR. (jsize_out /= size(field_out,2))
              set_mismatch = set_mismatch .OR. (ke_out /= size(field_out,3))
            endif
            if(set_mismatch)then
              write( text,'(i2)' ) l_size
              call mpp_error(FATAL,'MPP_REDISTRIBUTE_3D: Incompatible field at count '//text//' for group redistribute.' )
            endif
          endif  
          if(do_redist)then
            if(PRESENT(dc_handle))d_comm =>dc_handle  ! User has kept pointer to d_comm
            if(.not.ASSOCIATED(d_comm))then  ! d_comm needs initialization or lookup
              d_comm =>mpp_redistribute_init_comm(domain_in,l_addrs_in(1:l_size),domain_out,l_addrs_out(1:l_size), &
                                                                  isize_in,jsize_in,ke_in,isize_out,jsize_out,ke_out)
              if(PRESENT(dc_handle))dc_handle =>d_comm  ! User wants to keep pointer to d_comm
            endif
            call mpp_do_redistribute( l_addrs_in(1:l_size), l_addrs_out(1:l_size), d_comm, d_type )
            l_size=0; l_addrs_in=-9999; l_addrs_out=-9999
            isize_in=0;  jsize_in=0;  ke_in=0
            isize_out=0; jsize_out=0; ke_out=0
            d_comm =>NULL()
          endif
        endif
      else
        call mpp_do_redistribute( domain_in, field_in, domain_out, field_out )
      endif
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      call mpp_do_redistribute( domain_in, field_in, domain_out, field_out )
#endif
    end subroutine MPP_REDISTRIBUTE_3D_


    subroutine MPP_REDISTRIBUTE_4D_( domain_in, field_in, domain_out, field_out, complete, free, list_size, dc_handle, position )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:,:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:,:,:)
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      integer, intent(in), optional :: position
      MPP_TYPE_ :: field3D_in (size(field_in, 1),size(field_in, 2),size(field_in ,3)*size(field_in ,4))
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),size(field_out,3)*size(field_out,4))
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( ptr_in,  field3D_in  )
      pointer( ptr_out, field3D_out )
      ptr_in  = LOC(field_in )
      ptr_out = LOC(field_out)
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out, complete, free, list_size, dc_handle, position  )
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      field3D_in = RESHAPE( field_in, SHAPE(field3D_in) )
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out, position=position )
      field_out = RESHAPE( field3D_out, SHAPE(field_out) )
#endif
      return
    end subroutine MPP_REDISTRIBUTE_4D_

    subroutine MPP_REDISTRIBUTE_5D_( domain_in, field_in, domain_out, field_out, complete, free, list_size, dc_handle, position )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:,:,:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:,:,:,:)
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      integer, intent(in), optional :: position
      MPP_TYPE_ :: field3D_in (size(field_in, 1),size(field_in, 2),size(field_in ,3)*size(field_in ,4)*size(field_in ,5))
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),size(field_out,3)*size(field_out,4)*size(field_out,5))
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( ptr_in,  field3D_in  )
      pointer( ptr_out, field3D_out )
      ptr_in  = LOC(field_in )
      ptr_out = LOC(field_out)
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out, complete, free, list_size, dc_handle, position  )
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      field3D_in = RESHAPE( field_in, SHAPE(field3D_in) )
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out, position=position )
      field_out = RESHAPE( field3D_out, SHAPE(field_out) )
#endif
      return
    end subroutine MPP_REDISTRIBUTE_5D_

#ifdef VECTOR_FIELD_

!VECTOR_FIELD_ is set to false for MPP_TYPE_ integer or logical.
!vector fields
    subroutine MPP_UPDATE_DOMAINS_2D_V_( fieldx, fieldy, domain, flags, gridtype, complete, free, list_size, dc_handle )
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_, intent(inout), dimension(:,:) :: fieldx, fieldy
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype, complete, free, list_size, dc_handle )
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      field3Dx = RESHAPE( fieldx, SHAPE(field3Dx) )
      field3Dy = RESHAPE( fieldy, SHAPE(field3Dy) )
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
      fieldx = RESHAPE( field3Dx, SHAPE(fieldx) )
      fieldy = RESHAPE( field3Dy, SHAPE(fieldy) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_2D_V_


    subroutine MPP_UPDATE_DOMAINS_3D_V_( fieldx, fieldy, domain, flags, gridtype, complete, free, list_size, dc_handle )
!updates data domain of 3D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: fieldx(:,:,:), fieldy(:,:,:)
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size

#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      type(DomainCommunicator2D),pointer,save :: d_comm =>NULL()
      logical                       :: do_update, free_comm
      integer                       :: grid_offset_type
      integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: f_addrsx=-9999, f_addrsy=-9999
      integer, save :: isize(2)=0,jsize(2)=0,ke=0,l_size=0, offset_type=0
      logical       :: set_mismatch
      character(len=2) :: text
      MPP_TYPE_ :: d_type

      if(PRESENT(complete) .or. PRESENT(free))then        
         do_update=.true.; if(PRESENT(complete))do_update=complete
         free_comm=.false.; if(PRESENT(free))free_comm=free
         if(free_comm)then
            f_addrsx(1) = LOC(fieldx)
            f_addrsy(1) = LOC(fieldy)
            ke = size(fieldx,3)
            l_size=1; if(PRESENT(list_size))l_size=list_size
            call mpp_update_free_comm(domain,f_addrsx(1),ke,l_size,f_addrsy(1),flags=flags,gridtype=gridtype)
            l_size=0; f_addrsx=-9999; f_addrsy=-9999
         else
            l_size = l_size+1
            if(l_size > MAX_DOMAIN_FIELDS)then
               write( text,'(i2)' ) MAX_DOMAIN_FIELDS
               call mpp_error(FATAL,'MPP_UPDATE_3D_V: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
            end if
            grid_offset_type = AGRID
            if(present(gridtype)) grid_offset_type = gridtype
            f_addrsx(l_size) = LOC(fieldx); f_addrsy(l_size) = LOC(fieldy)
            if(l_size == 1)then
               isize(1)=size(fieldx,1); jsize(1)=size(fieldx,2); ke = size(fieldx,3)
               isize(2)=size(fieldy,1); jsize(2)=size(fieldy,2)
               offset_type = grid_offset_type
            else
               set_mismatch = .false.
               set_mismatch = set_mismatch .OR. (isize(1) /= size(fieldx,1))
               set_mismatch = set_mismatch .OR. (jsize(1) /= size(fieldx,2))
               set_mismatch = set_mismatch .OR. (ke /= size(fieldx,3))
               set_mismatch = set_mismatch .OR. (isize(2) /= size(fieldy,1))
               set_mismatch = set_mismatch .OR. (jsize(2) /= size(fieldy,2))
               set_mismatch = set_mismatch .OR. (ke /= size(fieldy,3))
               set_mismatch = set_mismatch .OR. (grid_offset_type /= offset_type)
               if(set_mismatch)then
                  write( text,'(i2)' ) l_size
                  call mpp_error(FATAL,'MPP_UPDATE_3D_V: Incompatible field at count '//text//' for group vector update.' )
               end if
            end if
            if(do_update)then
               if(PRESENT(dc_handle))d_comm =>dc_handle   ! User has kept pointer to d_comm
               if(.not.ASSOCIATED(d_comm))then  ! d_comm needs initialization or lookup
                  d_comm =>mpp_update_init_comm(domain,f_addrsx(1:l_size), isize,jsize,ke, f_addrsy(1:l_size), &
                       flags=flags, gridtype=grid_offset_type )
                  if(PRESENT(dc_handle))dc_handle =>d_comm  ! User wants to keep pointer to d_comm
               endif

               call mpp_do_update(f_addrsx(1:l_size),f_addrsy(1:l_size), d_comm, isize,jsize, &
                    ke,d_type,flags,gridtype)
               d_comm=>NULL(); l_size=0; f_addrsx=-9999; f_addrsy=-9999; isize=0;  jsize=0;  ke=0
            end if
         end if
      else
         call mpp_do_update( fieldx, fieldy, domain, flags, gridtype )
      end if
#else   
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      call mpp_do_update( fieldx, fieldy, domain, flags, gridtype )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_3D_V_


    subroutine MPP_UPDATE_DOMAINS_4D_V_( fieldx, fieldy, domain, flags, gridtype, complete, free, list_size, dc_handle )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout), dimension(:,:,:,:) :: fieldx, fieldy
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype, complete, free, list_size, dc_handle )
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      field3Dx = RESHAPE( fieldx, SHAPE(field3Dx) )
      field3Dy = RESHAPE( fieldy, SHAPE(field3Dy) )
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
      fieldx = RESHAPE( field3Dx, SHAPE(fieldx) )
      fieldy = RESHAPE( field3Dy, SHAPE(fieldy) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_V_

    subroutine MPP_UPDATE_DOMAINS_5D_V_( fieldx, fieldy, domain, flags, gridtype, complete, free, list_size, dc_handle )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout), dimension(:,:,:,:,:) :: fieldx, fieldy
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      logical, intent(in), optional :: complete, free
      integer, intent(in), optional :: list_size
      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4)*size(fieldx,5))
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4)*size(fieldy,5))
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype, complete, free, list_size, dc_handle )
#else
      integer, optional :: dc_handle  ! Not used when there are no Cray pointers
      field3Dx = RESHAPE( fieldx, SHAPE(field3Dx) )
      field3Dy = RESHAPE( fieldy, SHAPE(field3Dy) )
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
      fieldx = RESHAPE( field3Dx, SHAPE(fieldx) )
      fieldy = RESHAPE( field3Dy, SHAPE(fieldy) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_V_
#endif VECTOR_FIELD_
