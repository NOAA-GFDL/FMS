    subroutine MPP_DO_REDISTRIBUTE_3Dold_( domain_in, field_in, domain_out, field_out )
      type(domain2D), intent(in) :: domain_in, domain_out
!      MPP_TYPE_, intent(in)  :: field_in ( domain_in%x%data%begin:, domain_in%y%data%begin:,:)
!      MPP_TYPE_, intent(out) :: field_out(domain_out%x%data%begin:,domain_out%y%data%begin:,:)
      MPP_TYPE_, intent(in)  :: field_in (:,:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:,:)
      integer :: is, ie, js, je, ke, isc, iec, jsc, jec
      integer :: i, j, k
      integer :: list, m, n, pos, msgsize
      integer :: to_pe, from_pe
!      MPP_TYPE_, dimension(domain_in%x%compute%size*domain_in%y%compute%size*size(field_in,3)) :: send_buf, recv_buf
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif
      integer :: buffer_pos, wordlen
      character(len=8) :: text

!      ke = size(field_in,3)
!      if( ke.NE.size(field_out,3) )call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mismatch between field_in and field_out.' )
!      if( UBOUND(field_in,1).NE.domain_in%x%data%end .OR. &
!          UBOUND(field_in,2).NE.domain_in%y%data%end ) &
!          call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_in must be on data domain of domain_in.' )
!      if( UBOUND(field_out,1).NE.domain_out%x%data%end .OR. &
!          UBOUND(field_out,2).NE.domain_out%y%data%end ) &
!          call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_out must be on data domain of domain_out.' )

!fix ke
      ke = 0
      if( domain_in%pe.NE.NULL_PE )ke = size(field_in,3)
      if( domain_out%pe.NE.NULL_PE )then
          if( ke.NE.0 .AND. ke.NE.size(field_out,3) ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mismatch between field_in and field_out.' )
          ke = size(field_out,3)
      end if
      if( ke.EQ.0 )call mpp_error( FATAL, 'MPP_REDISTRIBUTE: either domain_in or domain_out must be native.' )
!check sizes
      if( domain_in%pe.NE.NULL_PE )then
          if( size(field_in,1).NE.domain_in%x%data%size .OR. size(field_in,2).NE.domain_in%y%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_in must be on data domain of domain_in.' )
      end if
      if( domain_out%pe.NE.NULL_PE )then
          if( size(field_out,1).NE.domain_out%x%data%size .OR. size(field_out,2).NE.domain_out%y%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_out must be on data domain of domain_out.' )
      end if

      buffer_pos = 0
#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
      wordlen = size(TRANSFER(buffer(1),mpp_domains_stack))
#endif
!send
      call mpp_get_compute_domain( domain_in, isc, iec, jsc, jec )
      n = size(domain_out%list(:))
      do list = 0,n-1
         m = mod( domain_out%pos+list+n, n )
         to_pe = domain_out%list(m)%pe
         call mpp_get_compute_domain( domain_out%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             msgsize = (ie-is+1)*(je-js+1)*ke
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             pos = buffer_pos
             do k = 1,ke
                do j = js-domain_in%y%data%begin+1,je-domain_in%y%data%begin+1
                   do i = is-domain_in%x%data%begin+1,ie-domain_in%x%data%begin+1
                      pos = pos+1
                      buffer(pos) = field_in(i,j,k)
                   end do
                end do
             end do
             if( debug )write( stderr(),* )'PE', pe, ' to PE ', to_pe, 'is,ie,js,je=', is, ie, js, je
             call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
             buffer_pos = pos
         end if
      end do
!recv
      call mpp_get_compute_domain( domain_out, isc, iec, jsc, jec )
      n = size(domain_in%list(:))
      do list = 0,n-1
         m = mod( domain_in%pos+n-list, n )
         from_pe = domain_in%list(m)%pe
         call mpp_get_compute_domain( domain_in%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             msgsize = (ie-is+1)*(je-js+1)*ke
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             if( debug )write( stderr(),* )'PE', pe, ' from PE ', from_pe, 'is,ie,js,je=', is, ie, js, je
             call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
             pos = buffer_pos
             do k = 1,ke
                do j = js-domain_out%y%data%begin+1,je-domain_out%y%data%begin+1
                   do i = is-domain_out%x%data%begin+1,ie-domain_out%x%data%begin+1
                      pos = pos+1
                      field_out(i,j,k) = buffer(pos)
                   end do
                end do
             end do
             buffer_pos = pos
         end if
      end do

!      call mpp_sync_self( domain_in%list(:)%pe )
      call mpp_sync_self()
      return
    end subroutine MPP_DO_REDISTRIBUTE_3Dold_
