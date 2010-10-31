subroutine MPP_GATHER_1D_(sbuf, rbuf)
   MPP_TYPE_, dimension(:),    intent(in) :: sbuf
   MPP_TYPE_, dimension(:), intent(inout) :: rbuf
   integer :: cnt, l, nproc

   nproc = mpp_npes()
   cnt = size(sbuf(:))
   if(size(rbuf(:)) .NE. cnt*nproc) call mpp_error(FATAL, &
          "MPP_GATHER_1D_: size(rbuf) should equal to npes*size(sbuf) ")

   !--- pre-post receiving
   if(pe == root_pe ) then
      rbuf(1:cnt) = sbuf
      do l = 1, nproc-1
         call mpp_recv(rbuf(l*cnt+1), glen=cnt, from_pe=root_pe+l, block=.FALSE. )
      enddo
   else
      call mpp_send(sbuf(1), plen=cnt, to_pe=root_pe)
   endif

   call mpp_sync_self(check=EVENT_RECV)
   call mpp_sync_self()

end subroutine MPP_GATHER_1D_
