    subroutine MPP_DO_GLOBAL_FIELD_3Dnew_( local, global, d_comm )
!get a global field from a local field
!local field may be on compute OR data domain
      type(DomainCommunicator2D), intent(in) :: d_comm
      MPP_TYPE_, intent(in)  ::  local(:,:,:)
      MPP_TYPE_, intent(out) :: global(d_comm%domain%x%global%begin:,d_comm%domain%y%global%begin:,:)
      integer :: i, j, k, m, n, nd, nwords, msgsize
      integer :: is, ie, js, je
      MPP_TYPE_ :: clocal (d_comm%domain%x%compute%size*    d_comm%domain%y%compute%size*    size(local,3))
      MPP_TYPE_ :: cremote(d_comm%domain%x%compute%max_size*d_comm%domain%y%compute%max_size*size(local,3))
      integer :: words_per_long, stackuse
      character(len=8) :: text
      pointer( ptr_local,  clocal  )
      pointer( ptr_remote, cremote )

      ptr_local  = LOC(mpp_domains_stack)
      ptr_remote = LOC(mpp_domains_stack(size(clocal(:))+1))

! make contiguous array from compute domain
      m = 0
      do k = 1,size(local,3)
         do j = d_comm%domain%y%compute%begin,d_comm%domain%y%compute%end
            do i = d_comm%domain%x%compute%begin,d_comm%domain%x%compute%end
               m = m + 1
               clocal(m) = local(i+d_comm%gf_ioff,j+d_comm%gf_joff,k)
               global(i,j,k) = clocal(m) !always fill local domain directly
            end do
         end do
      end do

      nd = d_comm%Rlist_size  ! same as size of send list
      msgsize = d_comm%S_msize(1,1)  ! constant for all sends
      do n = 1,nd-1
         call mpp_send( clocal(1), plen=msgsize, to_pe=d_comm%cto_pe(n) )
      end do
      do n = 1,nd-1
         nwords = d_comm%R_msize(1,n)
         call mpp_recv( cremote(1), glen=nwords, from_pe=d_comm%cfrom_pe(n) )
         is=d_comm%recvis(1,n); ie=d_comm%recvie(1,n)
         js=d_comm%recvjs(1,n); je=d_comm%recvje(1,n)
         m = 0
         do k = 1,size(global,3)
            do j = js,je
               do i = is,ie
                  m = m + 1
                  global(i,j,k) = cremote(m)
               end do
            end do
         end do
!     write(stdout(),*) 'new cremote chksum:', mpp_chksum(mpp_domains_stack(1:size(clocal(:))))
!     write(stdout(),*) 'new global chksum:', mpp_chksum(mpp_domains_stack(size(clocal(:))+1:size(global(:,:,:))))
      end do
      call mpp_sync_self()
    end subroutine MPP_DO_GLOBAL_FIELD_3Dnew_
