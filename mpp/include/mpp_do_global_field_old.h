    subroutine MPP_DO_GLOBAL_FIELD_3Dold_( domain, local, global, flags )
!get a global field from a local field
!local field may be on compute OR data domain
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:)
      MPP_TYPE_, intent(out) :: global(domain%x%global%begin:,domain%y%global%begin:,:)
      integer, intent(in), optional :: flags
      integer :: i, j, k, m, n, nd, nwords, lpos, rpos, ioff, joff
      logical :: xonly, yonly
      MPP_TYPE_ :: clocal (domain%x%compute%size*    domain%y%compute%size*    size(local,3))
      MPP_TYPE_ :: cremote(domain%x%compute%max_size*domain%y%compute%max_size*size(local,3))
      integer :: words_per_long, stackuse
      character(len=8) :: text
#ifdef use_CRI_pointers
      pointer( ptr_local,  clocal  )
      pointer( ptr_remote, cremote )

      ptr_local  = LOC(mpp_domains_stack)
      ptr_remote = LOC(mpp_domains_stack(size(clocal(:))+1))
#endif

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: must first call mpp_domains_init.' )
      xonly = .FALSE.
      yonly = .FALSE.
      if( PRESENT(flags) )then
          xonly = flags.EQ.XUPDATE
          yonly = flags.EQ.YUPDATE
          if( .NOT.xonly .AND. .NOT.yonly )call mpp_error( WARNING, 'MPP_GLOBAL_FIELD: you must have flags=XUPDATE or YUPDATE.' )
      end if

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
      if( size(global,1).NE.domain%x%global%size .OR. size(global,2).NE.domain%y%global%size .OR. &
          size(local,3).NE.size(global,3) ) &
           call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: incoming arrays do not match domain.' )
    if( size(local,1).EQ.domain%x%compute%size .AND. size(local,2).EQ.domain%y%compute%size )then
!local is on compute domain
        ioff = -domain%x%compute%begin + 1
        joff = -domain%y%compute%begin + 1
    else if( size(local,1).EQ.domain%x%data%size .AND. size(local,2).EQ.domain%y%data%size )then
!local is on data domain
        ioff = -domain%x%data%begin + 1
        joff = -domain%y%data%begin + 1
    else
        call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_: incoming field array must match either compute domain or data domain.' )
    end if

! make contiguous array from compute domain
      m = 0
      do k = 1,size(local,3)
         do j = domain%y%compute%begin,domain%y%compute%end
            do i = domain%x%compute%begin,domain%x%compute%end
               m = m + 1
               clocal(m) = local(i+ioff,j+joff,k)
               global(i,j,k) = clocal(m) !always fill local domain directly
            end do
         end do
      end do

!fill off-domains (note loops begin at an offset of 1)
      if( xonly )then
          nd = size(domain%x%list(:))
          do n = 1,nd-1
             lpos = mod(domain%x%pos+nd-n,nd)
             rpos = mod(domain%x%pos   +n,nd)
             nwords = domain%x%list(rpos)%compute%size * domain%x%list(rpos)%compute%size * size(local,3)
           ! Force use of scalar, integer ptr interface
             call mpp_transmit( put_data=clocal(1), plen=size(clocal(:)), to_pe=domain%x%list(lpos)%pe, &
                                get_data=cremote(1), glen=nwords, from_pe=domain%x%list(rpos)%pe )
             m = 0
             do k = 1,size(global,3)
                do j = domain%y%compute%begin,domain%y%compute%end
                   do i = domain%x%list(rpos)%compute%begin,domain%x%list(rpos)%compute%end
                      m = m + 1
                      global(i,j,k) = cremote(m)
                   end do
                end do
             end do
          end do
          call mpp_sync_self(domain%x%list(:)%pe)
      else if( yonly )then
          nd = size(domain%y%list(:))
          do n = 1,nd-1
             lpos = mod(domain%y%pos+nd-n,nd)
             rpos = mod(domain%y%pos   +n,nd)
             nwords = domain%y%list(rpos)%compute%size * domain%y%list(rpos)%compute%size * size(local,3)
           ! Force use of scalar, integer pointer interface
             call mpp_transmit( put_data=clocal(1), plen=size(clocal(:)), to_pe=domain%y%list(lpos)%pe, &
                                get_data=cremote(1), glen=nwords, from_pe=domain%y%list(rpos)%pe )
             m = 0
             do k = 1,size(global,3)
                do j = domain%y%list(rpos)%compute%begin,domain%y%list(rpos)%compute%end
                   do i = domain%x%compute%begin,domain%x%compute%end
                      m = m + 1
                      global(i,j,k) = cremote(m)
                   end do
                end do
             end do
          end do
          call mpp_sync_self(domain%y%list(:)%pe)
      else
          nd = size(domain%list(:))
          do n = 1,nd-1
             lpos = mod(domain%pos+nd-n,nd)
             call mpp_send( clocal(1), plen=size(clocal(:)), to_pe=domain%list(lpos)%pe )
          end do
          do n = 1,nd-1
             rpos = mod(domain%pos   +n,nd)
             nwords = domain%list(rpos)%x%compute%size * domain%list(rpos)%y%compute%size * size(local,3)
             call mpp_recv( cremote(1), glen=nwords, from_pe=domain%list(rpos)%pe )
             m = 0
             do k = 1,size(global,3)
                do j = domain%list(rpos)%y%compute%begin,domain%list(rpos)%y%compute%end
                   do i = domain%list(rpos)%x%compute%begin,domain%list(rpos)%x%compute%end
                      m = m + 1
                      global(i,j,k) = cremote(m)
                   end do
                end do
             end do
!     write(stdout(),*) 'old cremote chksum:', mpp_chksum(mpp_domains_stack(1:size(clocal(:))))
!     write(stdout(),*) 'old global chksum:', mpp_chksum(mpp_domains_stack(size(clocal(:))+1:size(global(:,:,:))))
          end do
!          call mpp_sync_self(domain%list(:)%pe)
          call mpp_sync_self()
      end if
          
      return
    end subroutine MPP_DO_GLOBAL_FIELD_3Dold_
