    subroutine MPP_UPDATE_DOMAINS_2D_( field, domain )
!updates data domain of 2D field whose computational domains have been computed
      type(domain1D), intent(in), target :: domain
      MPP_TYPE_, intent(inout) :: field(domain%data%start_index:,:)

      type(domain1D), pointer :: put_domain, get_domain
!limits of computation, remote put and get domains
      integer :: isc, iec,  isp, iep,  isg, ieg,  to, from

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: You must first call mpp_domains_init.' )

      isc = domain%compute%start_index; iec = domain%compute%end_index
!update left ("prev") boundary
      put_domain => domain; get_domain => domain
      do while( put_domain.NE.NULL_DOMAIN1D .OR. get_domain.NE.NULL_DOMAIN1D )
         if( ASSOCIATED(put_domain%next) )then
             put_domain => put_domain%next
         else
             put_domain => NULL_DOMAIN1D
         end if
         if( ASSOCIATED(get_domain%prev) )then
             get_domain => get_domain%prev
         else
             get_domain => NULL_DOMAIN1D
         end if
         call get_halos_1D( domain, put_domain, get_domain, -1, to, from, isp, iep, isg, ieg )
         if( to.NE.NULL_PE .OR. from.NE.NULL_PE )call buffer_and_transmit
         if( isg.EQ.domain%data%start_index )get_domain => NULL_DOMAIN1D
         if( put_domain.EQ.domain )put_domain => NULL_DOMAIN1D
      end do
      call mpp_sync()
!update right ("next") boundary
      put_domain => domain; get_domain => domain
      do while( put_domain.NE.NULL_DOMAIN1D .OR. get_domain.NE.NULL_DOMAIN1D )
         if( ASSOCIATED(put_domain%prev) )then
             put_domain => put_domain%prev
         else
             put_domain => NULL_DOMAIN1D
         end if
         if( ASSOCIATED(get_domain%next) )then
             get_domain => get_domain%next
         else
             get_domain => NULL_DOMAIN1D
         end if
         call get_halos_1D( domain, put_domain, get_domain, +1, to, from, isp, iep, isg, ieg )
         if( to.NE.NULL_PE .OR. from.NE.NULL_PE )call buffer_and_transmit
         if( ieg.EQ.domain%data%end_index )get_domain => NULL_DOMAIN1D
         if( put_domain.EQ.domain )put_domain => NULL_DOMAIN1D
      end do
      call mpp_sync()

      return

      contains
        subroutine buffer_and_transmit
!buffer the halo region and send
!isp, iep, jsp, jep: limits of put domain
!to, from: pe of put and get domains
          integer :: i, k
          integer :: offset, put_len, get_len
          character(len=8) :: text
          MPP_TYPE_ :: work(mpp_domains_stack_size)
#ifdef use_CRI_pointers
          pointer( ptr, work )
          ptr = LOC(mpp_domains_stack)
#endif

          if( debug )then
              call SYSTEM_CLOCK(tk)
              write( stdout,'(a,i18,a,i5,a,2i2,4i4)' ) &
                   'T=',tk, ' PE=',pe, ' BUFFER_AND_TRANSMIT: ', to, from, isp, iep, isg, ieg
          end if
          if( to.NE.NULL_PE )then  !put is to be done: buffer input
              put_len = (iep-isp+1)*size(field,2)
              call mpp_sync_self()  !check if work is still in use
              work(1:put_len) = TRANSFER( field(isp:iep,:), work )
          else
              put_len = 1
          end if
          get_len = (ieg-isg+1)*size(field,2)
          if( from.EQ.NULL_PE )get_len = 1
#ifdef use_CRI_pointers
          i = (put_len+get_len)*size(transfer(work(1),mpp_domains_stack))
#else
          i = (put_len+get_len)
#endif
          if( i.GT.mpp_domains_stack_size )then
              write( text,'(i8)' )i
              call mpp_error( FATAL, 'BUFFER_AND_TRANSMIT user stack overflow: call mpp_domains_set_stack_size('//text// &
                   ') from all PEs.' )
          end if
          call mpp_transmit( work(1), put_len, to, work(put_len+1), get_len, from )
          if( from.NE.NULL_PE )then  !get was done: unbuffer output
              field(isg:ieg,:) = RESHAPE( work(put_len+1:put_len+get_len), (/ieg-isg+1,size(field,2)/) )
          end if

          return
        end subroutine buffer_and_transmit
    end subroutine MPP_UPDATE_DOMAINS_2D_

    subroutine MPP_UPDATE_DOMAINS_3D_( field, domain )
!updates data domain of 3D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:)
      type(domain1D), intent(in), target :: domain
      MPP_TYPE_ :: field2D(size(field,1),size(field,2)*size(field,3))
#ifdef use_CRI_pointers
      pointer( ptr, field2D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call MPP_UPDATE_DOMAINS_2D_( field2D, domain )
      return
    end subroutine MPP_UPDATE_DOMAINS_3D_

    subroutine MPP_UPDATE_DOMAINS_4D_( field, domain )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:)
      type(domain1D), intent(in), target :: domain
      MPP_TYPE_ :: field2D(size(field,1),size(field,2)*size(field,3)*size(field,4))
#ifdef use_CRI_pointers
      pointer( ptr, field2D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call MPP_UPDATE_DOMAINS_2D_( field2D, domain )
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_

    subroutine MPP_UPDATE_DOMAINS_5D_( field, domain )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:,:)
      type(domain1D), intent(in), target :: domain
      MPP_TYPE_ :: field2D(size(field,1),size(field,2)*size(field,3)*size(field,4)*size(field,5))
#ifdef use_CRI_pointers
      pointer( ptr, field2D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call MPP_UPDATE_DOMAINS_2D_( field2D, domain )
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_
