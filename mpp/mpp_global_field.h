    subroutine MPP_GLOBAL_FIELD_( domain, local, global, flags )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(domain%x%data%begin:  ,domain%y%data%begin:   MPP_EXTRA_INDICES )
      MPP_TYPE_, intent(out) :: global(domain%x%global%begin:,domain%y%global%begin: MPP_EXTRA_INDICES )
      integer, intent(in) :: flags
      integer :: i, n, nwords, lpos, rpos, size_extra_indices
      MPP_TYPE_, allocatable, dimension(:) :: clocal, cremote
      logical :: xonly, yonly

      xonly = .FALSE.
      yonly = .FALSE.
      if( PRESENT(flags) )then
          xonly = flags.EQ.XUPDATE
          yonly = flags.EQ.YUPDATE
          if( .NOT.xonly .AND. .NOT.yonly )call mpp_error( WARNING, 'MPP_GLOBAL_FIELD: you must have flags=XUPDATE or YUPDATE.' )
      end if

      size_extra_indices = size(local)/(size(local,1)*size(local,2))
      if( size( local,1).NE.domain%x%data%size   .OR. size( local,2).NE.domain%y%data%size   .OR. &
          size(global,1).NE.domain%x%global%size .OR. size(global,2).NE.domain%y%global%size .OR. &
          size_extra_indices.NE.size(global)/(size(global,1)*size(global,2)) ) &
          call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: incoming arrays do not match domain.' )

      allocate( clocal (domain%x%compute%size*    domain%y%compute%size*    size_extra_indices) )
      allocate( cremote(domain%x%compute%max_size*domain%y%compute%max_size*size_extra_indices) )
      clocal = TRANSFER( local(domain%x%compute%begin:domain%x%compute%end, &
                               domain%y%compute%begin:domain%y%compute%end,:), clocal )
      global = 0.
      global(domain%x%compute%begin:domain%x%compute%end,domain%y%compute%begin:domain%y%compute%end,:) = &
       local(domain%x%compute%begin:domain%x%compute%end,domain%y%compute%begin:domain%y%compute%end,:)

      if( xonly )then
          n = size(domain%x%list)
          do i = 1,n-1
             lpos = mod(domain%x%pos+n-i,n)
             rpos = mod(domain%x%pos  +i,n)
             nwords = domain%x%list(rpos)%compute%size * domain%x%list(rpos)%compute%size * size_extra_indices
             call mpp_transmit( clocal, size(clocal), domain%x%list(lpos)%pe, cremote, nwords, domain%x%list(rpos)%pe )
             global(domain%x%list(rpos)%compute%begin:domain%x%list(rpos)%compute%end, &
                    domain%y%compute%begin:domain%y%compute%end,:) = &
                  RESHAPE( cremote(1:nwords), (/domain%x%list(rpos)%compute%size,domain%y%compute%size,size_extra_indices/) )
          end do
          call mpp_sync_self(domain%x%list(:)%pe)
      else if( yonly )then
          n = size(domain%y%list)
          do i = 1,n-1
             lpos = mod(domain%y%pos+n-i,n)
             rpos = mod(domain%y%pos  +i,n)
             nwords = domain%y%list(rpos)%compute%size * domain%y%list(rpos)%compute%size * size_extra_indices
             call mpp_transmit( clocal, size(clocal), domain%y%list(lpos)%pe, cremote, nwords, domain%y%list(rpos)%pe )
             global(domain%x%compute%begin:domain%x%compute%end,&
                    domain%y%list(rpos)%compute%begin:domain%y%list(rpos)%compute%end,:) = &
                  RESHAPE( cremote(1:nwords), (/domain%x%compute%size,domain%y%list(rpos)%compute%size,size_extra_indices/) )
          end do
          call mpp_sync_self(domain%y%list(:)%pe)
      else
          n = size(domain%list)
          do i = 1,n-1
             lpos = mod(domain%pos+n-i,n)
             rpos = mod(domain%pos  +i,n)
             nwords = domain%list(rpos)%x%compute%size * domain%list(rpos)%y%compute%size * size_extra_indices
             call mpp_transmit( clocal, size(clocal), domain%list(lpos)%pe, cremote, nwords, domain%list(rpos)%pe )
             global(domain%list(rpos)%x%compute%begin:domain%list(rpos)%x%compute%end, &
                    domain%list(rpos)%y%compute%begin:domain%list(rpos)%y%compute%end,:) = &
                  RESHAPE( cremote(1:nwords), &
                  (/domain%list(rpos)%x%compute%size,domain%list(rpos)%y%compute%size,size_extra_indices/) )
          end do
          call mpp_sync_self(domain%list(:)%pe)
      end if
          
      return
    end subroutine MPP_GLOBAL_FIELD_
