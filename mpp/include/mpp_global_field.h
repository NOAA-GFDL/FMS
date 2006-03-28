    subroutine MPP_GLOBAL_FIELD_2D_( domain, local, global, flags, new, dc_handle)
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:)
      MPP_TYPE_, intent(out) :: global(:,:)
      integer, intent(in), optional :: flags
      logical, intent(in), optional :: new

      MPP_TYPE_ :: local3D (size( local,1),size( local,2),1)
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),1)
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( lptr,  local3D )
      pointer( gptr, global3D )
      lptr = LOC( local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags, new, dc_handle )
#else
      integer, optional :: dc_handle  ! not used when there are no cray pointers
      local3D = RESHAPE( local, SHAPE(local3D) )
      call mpp_global_field( domain, local3D, global3D, flags )
      global = RESHAPE( global3D, SHAPE(global) )
#endif
    end subroutine MPP_GLOBAL_FIELD_2D_

    subroutine MPP_GLOBAL_FIELD_3D_( domain, local, global, flags, new, dc_handle)
!get a global field from a local field
!local field may be on compute OR data domain
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:)
      integer, intent(in), optional :: flags
      logical, intent(in), optional :: new

#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      type(DomainCommunicator2D),pointer,save :: d_comm =>NULL()
      integer(LONG_KIND)     :: l_addr
      logical                :: use_new
      integer :: isize_l,jsize_l,isize_g,jsize_g,ke

      use_new=.false.; if(PRESENT(new))use_new=new

      if(use_new)then
        ke = size(local,3)
        if(ke /= size(global,3))call mpp_error(FATAL, 'MPP_GLOBAL_FIELD: incoming arrays have mismatched K extent.')
        isize_l = size(local,1); jsize_l = size(local,2)
        isize_g = size(global,1); jsize_g = size(global,2)

        l_addr = LOC(local)
        if(PRESENT(dc_handle))d_comm =>dc_handle  ! User has kept pointer to d_comm
        if(.not.ASSOCIATED(d_comm))then  ! d_comm needs initialization or lookup
          d_comm =>mpp_global_field_init_comm(domain,l_addr,isize_g,jsize_g,isize_l,jsize_l,ke,flags=flags)
          if(PRESENT(dc_handle))dc_handle =>d_comm  ! User wants to keep pointer to d_comm
        endif
        call mpp_do_global_field( local, global, d_comm )
        d_comm =>NULL()
      else
        call mpp_do_global_field( domain, local, global, flags )
      end if
#else
      integer, optional :: dc_handle  ! not used when there are no cray pointers
      call mpp_do_global_field( domain, local, global, flags )
#endif
    end subroutine MPP_GLOBAL_FIELD_3D_

    subroutine MPP_GLOBAL_FIELD_4D_( domain, local, global, flags, new, dc_handle )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:,:)
      integer, intent(in), optional :: flags
      logical, intent(in), optional :: new
      MPP_TYPE_ :: local3D (size( local,1),size( local,2),size( local,3)*size(local,4))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),size(global,3)*size(local,4))
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags, new, dc_handle )
#else
      integer, optional :: dc_handle  ! not used when there are no cray pointers
      local3D = RESHAPE( local, SHAPE(local3D) )
      call mpp_global_field( domain, local3D, global3D, flags  )
      global = RESHAPE( global3D, SHAPE(global) )
#endif
    end subroutine MPP_GLOBAL_FIELD_4D_

    subroutine MPP_GLOBAL_FIELD_5D_( domain, local, global, flags, new, dc_handle )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:,:,:)
      integer, intent(in), optional :: flags
      logical, intent(in), optional :: new

      MPP_TYPE_ :: local3D (size( local,1),size( local,2),size( local,3)*size( local,4)*size(local,5))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),size(global,3)*size(global,4)*size(local,5))
#ifdef use_CRI_pointers
      type(DomainCommunicator2D),pointer,optional :: dc_handle
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags, new, dc_handle )
#else
      integer, optional :: dc_handle  ! not used when there are no cray pointers
      local3D = RESHAPE( local, SHAPE(local3D) )
      call mpp_global_field( domain, local3D, global3D, flags )
      global = RESHAPE( global3D, SHAPE(global) )
#endif
    end subroutine MPP_GLOBAL_FIELD_5D_
