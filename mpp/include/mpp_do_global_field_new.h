    subroutine MPP_DO_GLOBAL_FIELD_3Dnew_( local, global, d_comm )
!get a global field from a local field
!local field may be on compute OR data domain
      type(DomainCommunicator2D), intent(in) :: d_comm
      MPP_TYPE_, intent(in)  ::  local(:,:,:)
      MPP_TYPE_, intent(out) :: global(d_comm%domain%x%global%begin:,d_comm%domain%y%global%begin:,:)
      integer :: i, j, k, m, n, nd, nwords, msgsize
      integer :: is, ie, js, je
      MPP_TYPE_ :: clocal ((d_comm%domain%x%compute%size+1)*    (d_comm%domain%y%compute%size+1)*    size(local,3))
      MPP_TYPE_ :: cremote((d_comm%domain%x%compute%max_size+1)*(d_comm%domain%y%compute%max_size+1)*size(local,3))
      integer :: words_per_long, stackuse
      character(len=8) :: text
      pointer( ptr_local,  clocal  )
      pointer( ptr_remote, cremote )

      ptr_local  = LOC(mpp_domains_stack)
      ptr_remote = LOC(mpp_domains_stack(size(clocal(:))+1))

#ifdef use_CRI_pointers
      words_per_long = size(transfer(local(1,1,1),mpp_domains_stack))
      stackuse = (size(clocal(:))+size(cremote(:)))*words_per_long
      if( stackuse.GT.mpp_domains_stack_size )then
          write( text, '(i8)' )stackuse
          call mpp_error( FATAL, &
               'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
      end if
      mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, stackuse )
#endif

! make contiguous array from compute domain
      m = 0
      is=d_comm%sendis(1,0); ie=d_comm%sendie(1,0)
      js=d_comm%sendjs(1,0); je=d_comm%sendje(1,0)
      do k = 1,size(local,3)
         do j = js, je
            do i = is, ie
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
