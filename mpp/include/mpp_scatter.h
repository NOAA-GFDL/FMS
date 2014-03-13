subroutine MPP_SCATTER_PELIST_2D_(is, ie, js, je, pelist, array_seg, data, is_root_pe, &
                                  ishift, jshift)
   integer,                           intent(in)    :: is, ie, js, je
   integer,   dimension(:),           intent(in)    :: pelist
   MPP_TYPE_, dimension(is:ie,js:je), intent(inout)    :: array_seg
   MPP_TYPE_, dimension(:,:),         intent(in) :: data
   logical,                           intent(in)    :: is_root_pe
   integer,   optional,               intent(in)    :: ishift, jshift

   MPP_TYPE_ ::  arr3D(size(array_seg,1),size(array_seg,2),1)
   MPP_TYPE_ :: data3D(size(     data,1),size(     data,2),1)
   pointer( aptr,  arr3D )
   pointer( dptr, data3D )
   aptr = LOC(array_seg)
   dptr = LOC(     data)

   call mpp_scatter(is, ie, js, je, 1, pelist, arr3D, data3D, is_root_pe, &
                    ishift, jshift)
   return

end subroutine MPP_SCATTER_PELIST_2D_


subroutine MPP_SCATTER_PELIST_3D_(is, ie, js, je, nk, pelist, array_seg, data, is_root_pe, &
                                  ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   MPP_TYPE_, dimension(is:ie,js:je,1:nk), intent(inout)    :: array_seg
   MPP_TYPE_, dimension(:,:,:),            intent(in) :: data
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
                "fms_io(mpp_scatter_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif
! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): too many root_pes specified")


   ioff=0
   joff=0
   if (present(ishift)) ioff=ishift
   if (present(jshift)) joff=jshift

   my_ind(1) = is
   my_ind(2) = ie
   my_ind(3) = js
   my_ind(4) = je

! scatter indices into global index on root_pe
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
         call mpp_error(FATAL,"fms_io(mpp_scatter_pelist): specified indices (with shift) are outside of the & 
                        &range of the receiving array")
   else
! non root_pe's send indices to root_pe
     call mpp_send(my_ind(:), 4, root_pe, COMM_TAG_1)
     call mpp_sync_self(check=EVENT_SEND)
   endif

!  scatter segments into data based on indices
   if (is_root_pe) then
     do i = 1, size(pelist)
       if (pelist(i).ne.root_pe) then    ! no send to self
         i1 = gind(1,i)
         i2 = gind(2,i)
         j1 = gind(3,i)
         j2 = gind(4,i)
         msgsize = (i2-i1+1)*(j2-j1+1)*nk
! allocate and copy data into a contiguous memory space
         allocate(temp(i)%data(i1:i2,j1:j2,1:nk))
         temp(i)%data(i1:i2,j1:j2,1:nk)=data(i1:i2,j1:j2,1:nk)
         call mpp_send(temp(i)%data, msgsize, pelist(i), COMM_TAG_2)
       else
!        data copy - no send to self
         array_seg(is:ie,js:je,1:nk) = data(is+ioff:ie+ioff,js+joff:je+joff,1:nk)
       endif
     enddo
     call mpp_sync_self(check=EVENT_SEND)
! deallocate the temporary array used for the send
     do i = 1, size(pelist)
       if (allocated(temp(i)%data)) deallocate(temp(i)%data)
     enddo
   else
!    non root_pe's recv data from root_pe
     msgsize = (my_ind(2)-my_ind(1)+1) * (my_ind(4)-my_ind(3)+1) * nk
     call mpp_recv(array_seg, msgsize, root_pe, .FALSE., COMM_TAG_2)
     call mpp_sync_self(check=EVENT_RECV)
   endif

   call mpp_sync_self()

   return

end subroutine MPP_SCATTER_PELIST_3D_
