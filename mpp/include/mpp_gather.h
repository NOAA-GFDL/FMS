subroutine MPP_GATHER_1D_(sbuf, rbuf,pelist)
! JWD: Did not create mpp_gather_2d because have no requirement for it
! JWD: See mpp_gather_2dv below
   MPP_TYPE_, dimension(:),    intent(in) :: sbuf
   MPP_TYPE_, dimension(:), intent(inout) :: rbuf
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: cnt, l, nproc, op_root
   integer, allocatable :: pelist2(:)


!  If pelist is provided, the first position must be
!  the operation root
   if(PRESENT(pelist))then
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   cnt = size(sbuf(:))
   if(size(rbuf(:)) < cnt*nproc) call mpp_error(FATAL, &
          "MPP_GATHER_1D_: size(rbuf) must be at least npes*size(sbuf) ")

   !--- pre-post receiving
   if(pe == op_root) then
      rbuf(1:cnt) = sbuf
      do l = 2, nproc
         call mpp_recv(rbuf((l-1)*cnt+1), glen=cnt, from_pe=pelist2(l), block=.FALSE., tag=COMM_TAG_1 )
      enddo
   else
      call mpp_send(sbuf(1), plen=cnt, to_pe=op_root, tag=COMM_TAG_1)
   endif

   call mpp_sync_self(check=EVENT_RECV)
   call mpp_sync_self()
   deallocate(pelist2)
end subroutine MPP_GATHER_1D_

subroutine MPP_GATHER_1DV_(sbuf, ssize, rbuf, rsize, pelist)
   MPP_TYPE_, dimension(:),    intent(in) :: sbuf
   MPP_TYPE_, dimension(:), intent(inout) :: rbuf
   integer,                    intent(in) :: ssize
   integer,   dimension(:),    intent(in) :: rsize
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: cnt, l, nproc, pos, op_root
   integer, allocatable :: pelist2(:)

!  If pelist is provided, the first position must be
!  the operation root
   if(PRESENT(pelist))then
      nproc = size(pelist)
      allocate(pelist2(nproc)) 
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))        
      pelist2 = (/ (l, l=0+root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)


   !--- pre-post receiving
   if (pe .eq. op_root) then
       pos = 1
       do l = 1,nproc   ! include op_root to simplify logic
           if (rsize(l) == 0) then
               cycle  ! avoid ranks with no data
           endif
           call mpp_recv(rbuf(pos),glen=rsize(l),from_pe=pelist2(l), &
                         block=.FALSE.,tag=COMM_TAG_2)
           pos = pos + rsize(l)
       enddo
   endif
   if (ssize .gt. 0) then
       call mpp_send(sbuf(1),plen=ssize,to_pe=op_root,tag=COMM_TAG_2) !avoid ranks with no data
   endif

   call mpp_sync_self(check=EVENT_RECV)
   call mpp_sync_self()
   deallocate(pelist2)
end subroutine MPP_GATHER_1DV_


subroutine MPP_GATHER_PELIST_2D_(is, ie, js, je, pelist, array_seg, data, is_root_pe, &
                                 ishift, jshift)
   integer,                           intent(in)    :: is, ie, js, je
   integer,   dimension(:),           intent(in)    :: pelist
   MPP_TYPE_, dimension(is:ie,js:je), intent(in)    :: array_seg
   MPP_TYPE_, dimension(:,:),         intent(inout) :: data
   logical,                           intent(in)    :: is_root_pe
   integer,   optional,               intent(in)    :: ishift, jshift

   MPP_TYPE_ ::  arr3D(size(array_seg,1),size(array_seg,2),1)
   MPP_TYPE_ :: data3D(size(     data,1),size(     data,2),1)
   pointer( aptr,  arr3D )
   pointer( dptr, data3D )
   aptr = LOC(array_seg)
   dptr = LOC(     data)

   call mpp_gather(is, ie, js, je, 1, pelist, arr3D, data3D, is_root_pe, &
                   ishift, jshift)
   return

end subroutine MPP_GATHER_PELIST_2D_


subroutine MPP_GATHER_PELIST_3D_(is, ie, js, je, nk, pelist, array_seg, data, is_root_pe, &
                                 ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   MPP_TYPE_, dimension(is:ie,js:je,1:nk), intent(in)    :: array_seg
   MPP_TYPE_, dimension(:,:,:),            intent(inout) :: data
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer :: i, msgsize, root_pe, root_pe_test
   integer :: i1, i2, j1, j2, ioff, joff
   integer :: my_ind(4), gind(4,size(pelist))
   type array3D
     MPP_TYPE_, dimension(:,:,:), allocatable :: data
   endtype array3D
   type(array3d), dimension(size(pelist)) :: temp

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif
! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): too many root_pes specified")


   ioff=0
   joff=0
   if (present(ishift)) ioff=ishift
   if (present(jshift)) joff=jshift

   my_ind(1) = is
   my_ind(2) = ie
   my_ind(3) = js
   my_ind(4) = je

! gather indices into global index on root_pe
   if (is_root_pe) then
     do i = 1, size(pelist)
! root_pe data copy - no send to self
       if (pelist(i).eq.root_pe) then
         gind(:,i) = my_ind(:)
       else
         call mpp_recv(gind(:,i:i), 4, pelist(i), .FALSE., COMM_TAG_1)
       endif
     enddo
     call mpp_sync_self(check=EVENT_RECV)
     gind(1,:)=gind(1,:)+ioff
     gind(2,:)=gind(2,:)+ioff
     gind(3,:)=gind(3,:)+joff
     gind(4,:)=gind(4,:)+joff
! check indices to make sure they are within the range of "data"
     if ((minval(gind).lt.1) .OR. (maxval(gind(1:2,:)).gt.size(data,1)) .OR. (maxval(gind(3:4,:)).gt.size(data,2))) &
         call mpp_error(FATAL,"fms_io(mpp_gather_pelist): specified indices (with shift) are outside of the & 
                        &range of the receiving array")
   else
! non root_pe's send indices to root_pe
     call mpp_send(my_ind(:), 4, root_pe, COMM_TAG_1)
     call mpp_sync_self(check=EVENT_SEND)
   endif

!  gather segments into data based on indices
   if (is_root_pe) then
     do i = 1, size(pelist)
       if (pelist(i).ne.root_pe) then    ! no send to self
         i1 = gind(1,i)
         i2 = gind(2,i)
         j1 = gind(3,i)
         j2 = gind(4,i)
         msgsize = (i2-i1+1)*(j2-j1+1)*nk
         allocate(temp(i)%data(i1:i2,j1:j2,1:nk))
         call mpp_recv(temp(i)%data(i1:i2,j1:j2,1:nk), msgsize, pelist(i), .FALSE., COMM_TAG_2)
       endif
     enddo
     call mpp_sync_self(check=EVENT_RECV)
!  unbuffer/copy the data into the return array
     do i = 1, size(pelist)
       if (pelist(i).eq.root_pe) then
!        data copy - no send to self
         data(is+ioff:ie+ioff,js+joff:je+joff,1:nk) = array_seg(is:ie,js:je,1:nk)
       else
         i1 = gind(1,i)
         i2 = gind(2,i)
         j1 = gind(3,i)
         j2 = gind(4,i)
         data(i1:i2,j1:j2,1:nk)=temp(i)%data(i1:i2,j1:j2,1:nk)
         deallocate(temp(i)%data)
       endif
     enddo
   else
!    non root_pe's send data to root_pe
     msgsize = (my_ind(2)-my_ind(1)+1) * (my_ind(4)-my_ind(3)+1) * nk
     call mpp_send(array_seg, msgsize, root_pe, COMM_TAG_2)
     call mpp_sync_self(check=EVENT_SEND)
   endif

   call mpp_sync_self()
   return

end subroutine MPP_GATHER_PELIST_3D_
