    subroutine MPP_UPDATE_DOMAINS_2D_( field, domain, flags )
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:)
      type(domain2D), intent(in), target :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call MPP_UPDATE_DOMAINS_3D_( field3D, domain, flags )
      return
    end subroutine MPP_UPDATE_DOMAINS_2D_

    subroutine MPP_UPDATE_DOMAINS_3D_( field, domain, flags_in )
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(in), target :: domain
      MPP_TYPE_, intent(inout) :: field(domain%x%data%start_index:,domain%y%data%start_index:,:)
      integer, intent(in), optional :: flags_in
      integer :: flags

      type(domain2D), pointer :: put_domain, get_domain
!limits of computation, remote put and get domains
      integer :: isc, iec, jsc, jec,  isp, iep, jsp, jep,  isg, ieg, jsg, jeg,  to, from

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: You must first call mpp_domains_init.' )

      flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags_in) )flags = flags_in
      if( flags.EQ.TRANSPOSE ) &
           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: if TRANSPOSE is set, an update direction must be set also.' )

      isc = domain%x%compute%start_index; iec = domain%x%compute%end_index
      jsc = domain%y%compute%start_index; jec = domain%y%compute%end_index
!      call mpp_sync()
      if( BTEST(flags,0) )then  !WUPDATE: update western halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%east) )then
                 put_domain => put_domain%east
                 if( BTEST(flags,4) )then
                     jsp = put_domain%y%compute%start_index; jep = put_domain%y%compute%end_index
                 else
                     jsp = jsc; jep = jec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%west) )then
                 get_domain => get_domain%west
                 jsg = jsc; jeg = jec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%x, put_domain%x, get_domain%x, -1, to, from, isp, iep, isg, ieg )
             if( to.NE.NULL_PE .OR. from.NE.NULL_PE )call buffer_and_transmit
             if( isg.EQ.domain%x%data%start_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
          isc = domain%x%data%start_index !reset compute domain left edge so that Y update will do corners
      end if
!      call mpp_sync()
      if( BTEST(flags,1) )then  !EUPDATE: update eastern halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%west) )then
                 put_domain => put_domain%west
                 if( BTEST(flags,4) )then
                     jsp = put_domain%y%compute%start_index; jep = put_domain%y%compute%end_index
                 else
                     jsp = jsc; jep = jec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%east) )then
                 get_domain => get_domain%east
                 jsg = jsc; jeg = jec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%x, put_domain%x, get_domain%x, +1, to, from, isp, iep, isg, ieg )
             if( to.NE.NULL_PE .OR. from.NE.NULL_PE )call buffer_and_transmit
             if( ieg.EQ.domain%x%data%end_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
          iec = domain%x%data%end_index !reset compute domain right edge so that Y update will do corners
      end if
!      call mpp_sync()
      if( BTEST(flags,2) )then  !SUPDATE: update southern halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%north) )then
                 put_domain => put_domain%north
                 if( BTEST(flags,4) )then
                     isp = put_domain%x%compute%start_index; iep = put_domain%x%compute%end_index
                 else
                     isp = isc; iep = iec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%south) )then
                 get_domain => get_domain%south
                 isg = isc; ieg = iec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%y, put_domain%y, get_domain%y, -1, to, from, jsp, jep, jsg, jeg )
             if( to.NE.NULL_PE .OR. from.NE.NULL_PE )call buffer_and_transmit
             if( jsg.EQ.domain%y%data%start_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
      end if
!      call mpp_sync()
      if( BTEST(flags,3) )then  !NUPDATE: update northern halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%south) )then
                 put_domain => put_domain%south
                 if( BTEST(flags,4) )then
                     isp = put_domain%x%compute%start_index; iep = put_domain%x%compute%end_index
                 else
                     isp = isc; iep = iec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%north) )then
                 get_domain => get_domain%north
                 isg = isc; ieg = iec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%y, put_domain%y, get_domain%y, +1, to, from, jsp, jep, jsg, jeg )
             if( to.NE.NULL_PE .OR. from.NE.NULL_PE )call buffer_and_transmit
             if( jeg.EQ.domain%y%data%end_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
      end if
      call mpp_sync()

      return

      contains
        subroutine buffer_and_transmit
!buffer the halo region and send
!isp, iep, jsp, jep: limits of put domain
!isg, ieg, jsg, jeg: limits of get domain
!to, from: pe of put and get domains
          integer :: i, j, k
          integer :: put_len, get_len
          character(len=8) :: text
          MPP_TYPE_ :: work(size(field))
#ifdef use_CRI_pointers
          pointer( ptr, work )
          ptr = LOC(mpp_domains_stack)
#else
          equivalence( work, mpp_domains_stack )
#endif

          if( debug )then
              call SYSTEM_CLOCK(tk)
              write( stdout,'(a,i18,a,i5,a,2i2,8i4)' ) &
                   'T=',tk, ' PE=',pe, ' BUFFER_AND_TRANSMIT: ', to, from, isp, iep, jsp, jep, isg, ieg, jsg, jeg
          end if
          if( to.NE.NULL_PE )then  !put is to be done: buffer input
              put_len = (iep-isp+1)*(jep-jsp+1)*size(field,3)
              call mpp_sync_self()  !check if put_r8 is still in use
              work(1:put_len) = TRANSFER( field(isp:iep,jsp:jep,:), work )
          else
              put_len = 1
          end if
          get_len = (ieg-isg+1)*(jeg-jsg+1)*size(field,3)
          if( from.EQ.NULL_PE )get_len = 1
          i = (put_len+get_len)*size(transfer(work(1),mpp_domains_stack))
          if( i.GT.mpp_domains_stack_size )then
              write( text,'(i8)' )i
              call mpp_error( FATAL, 'BUFFER_AND_TRANSMIT user stack overflow: call mpp_domains_set_stack_size('//text// &
                   ') from all PEs.' )
          end if
          
          call mpp_transmit( work(1), put_len, to, work(put_len+1), get_len, from )
          if( from.NE.NULL_PE )then  !get was done: unbuffer output
              field(isg:ieg,jsg:jeg,:) = &
                   RESHAPE( work(put_len+1:put_len+get_len), (/ieg-isg+1,jeg-jsg+1,size(field,3)/) )
          end if

          return
        end subroutine buffer_and_transmit
    end subroutine MPP_UPDATE_DOMAINS_3D_

    subroutine MPP_UPDATE_DOMAINS_4D_( field, domain, flags )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:)
      type(domain2D), intent(in), target :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call MPP_UPDATE_DOMAINS_3D_( field3D, domain, flags )
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_

    subroutine MPP_UPDATE_DOMAINS_5D_( field, domain, flags )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:,:)
      type(domain2D), intent(in), target :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call MPP_UPDATE_DOMAINS_3D_( field3D, domain, flags )
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_
