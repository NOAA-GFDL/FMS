    subroutine MPP_UPDATE_DOMAINS_2D_( field, domain, flags_in )
!updates data domain of 3D field whose computational domains have been computed
      type(domain1D), intent(inout), target :: domain
      MPP_TYPE_, intent(inout) :: field(domain%data%begin:,:)
      integer, intent(in), optional :: flags_in

      integer :: is, ie, i, j, k, l, flags
      integer :: stkpos, stkpos0, nwords
      logical :: mustput, mustget
      character(len=8) :: text
      MPP_TYPE_ :: work(size(field))
      integer :: words_per_long, mpp_domains_stack_pos
#ifdef use_CRI_pointers
      pointer( ptr, work )
      ptr = LOC(mpp_domains_stack)
#else
      equivalence( work, mpp_domains_stack )
#endif

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: You must first call mpp_domains_init.' )

      flags=XUPDATE
      if( PRESENT(flags_in) )flags = flags_in
!currently flags is unused in mpp_update_domains1D
      if( debug )write( stderr, '(a,3i4)' )'pe, lbound(field,1), size(field,1)=', &
                                            pe, lbound(field,1), size(field,1)

      call mpp_get_compute_domain( domain, is, ie )
      call mpp_set_active_domain ( domain, is, ie )

      words_per_long = size(transfer(work(1),mpp_domains_stack))
!put forward and backward halos
      stkpos = 0
      mpp_domains_stack_pos = 0
      do l = 0,size(domain%list)-1
         mustput = .FALSE.
         stkpos0 = stkpos
         if( BTEST(flags,0) .AND. domain%list(l)%mustputf )then
             mustput = .TRUE.
             nwords = domain%list(l)%putf%size * size(field,2)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,3i4,i5)')'->W: pe, (beg,end), is,ie,k,   to=', pe, stkpos+1, stkpos+nwords, &
                      domain%list(l)%putf%begin, domain%list(l)%putf%end, size(field,2), domain%list(l)%pe
             end if
             mpp_domains_stack_pos = (stkpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
!is reshape better than pack? 
!             work(stkpos+1:stkpos+nwords) = RESHAPE( field(domain%list(l)%putf%begin:domain%list(l)%putf%end,:), (/nwords/) )
!             stkpos = stkpos+nwords
             do j = 1,size(field,2)
                do i = domain%list(l)%putf%begin,domain%list(l)%putf%end
                   stkpos = stkpos + 1
                   work(stkpos) = field(i,j)
                end do
             end do
         end if
         if( BTEST(flags,1) .AND. domain%list(l)%mustputb )then
             mustput = .TRUE.
             nwords = domain%list(l)%putb%size * size(field,2)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,3i4,i5)')'->E: pe, (beg,end), is,ie,k,   to=', pe, stkpos+1, stkpos+nwords, &
                      domain%list(l)%putb%begin, domain%list(l)%putb%end, size(field,2), domain%list(l)%pe
             end if
             mpp_domains_stack_pos = (stkpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
!is reshape better than pack? 
!             work(stkpos+1:stkpos+nwords) = RESHAPE( field(domain%list(l)%putb%begin:domain%list(l)%putb%end,:), (/nwords/) )
!             stkpos = stkpos+nwords
             do j = 1,size(field,2)
                do i = domain%list(l)%putb%begin,domain%list(l)%putb%end
                   stkpos = stkpos + 1
                   work(stkpos) = field(i,j)
                end do
             end do
         end if
         if( mustput )call mpp_send( work(stkpos0+1:stkpos), stkpos-stkpos0, domain%list(l)%pe )
         if( debug )write( stderr, '(a,i4/(8f9.6))' )'pe, work=', pe, work(stkpos0+1:stkpos)
      end do
!get forward and backward halos
      do l = 0,size(domain%list)-1
         stkpos0 = stkpos
         mustget = .FALSE.
         if( BTEST(flags,0) .AND. domain%list(l)%mustgetb )then
             mustget = .TRUE.
             nwords = domain%list(l)%getb%size * size(field,2)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,3i4,i5)')'<-W: pe, (beg,end), is,ie,k, from=', pe, stkpos+1, stkpos+nwords, &
                      domain%list(l)%getb%begin, domain%list(l)%getb%end, size(field,2), domain%list(l)%pe
             end if
             mpp_domains_stack_pos = (stkpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             stkpos = stkpos+nwords
         end if
         if( BTEST(flags,1) .AND. domain%list(l)%mustgetf )then
             mustget = .TRUE.
             nwords = domain%list(l)%getf%size * size(field,2)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,3i4,i5)')'<-E: pe, (beg,end), is,ie,js,je,k, from=', pe, stkpos+1, stkpos+nwords, &
                      domain%list(l)%getf%begin, domain%list(l)%getf%end, size(field,2), domain%list(l)%pe
             end if
             mpp_domains_stack_pos = (stkpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             stkpos = stkpos+nwords
         end if
         if( mustget )then
             call mpp_recv( work(stkpos0+1:stkpos), stkpos-stkpos0, domain%list(l)%pe )
             if( debug )write( stderr, '(a,i4/(8f9.6))' )'pe, work=', pe, work(stkpos0+1:stkpos)
             stkpos = stkpos0
             if( BTEST(flags,0) .AND. domain%list(l)%mustgetb )then
!is reshape better than unpack? 
!                 field(domain%list(l)%getb%begin:domain%list(l)%getb%end,:) = &
!                      RESHAPE( work(stkpos+1:stkpos+nwords), (/domain%list(l)%getb%size,size(field,2)/) )
!                 stkpos = stkpos+nwords
                 do j = 1,size(field,2)
                    do i = domain%list(l)%getb%begin,domain%list(l)%getb%end
                       stkpos = stkpos + 1
                       field(i,j) = work(stkpos)
                    end do
                 end do
                 if( debug )write( stderr, '(a,3i4/(8f9.6))' ) 'getb: pe, is, ie=', &
                        pe, domain%list(l)%getb%begin,domain%list(l)%getb%end, &
                      field(domain%list(l)%getb%begin:domain%list(l)%getb%end,:) 
                 call mpp_set_active_domain( domain, begin=domain%list(l)%getb%begin )
             end if
             if( BTEST(flags,1) .AND. domain%list(l)%mustgetf )then
!is reshape better than unpack? 
!                 field(domain%list(l)%getf%begin:domain%list(l)%getf%end,:) = &
!                      RESHAPE( work(stkpos+1:stkpos+nwords), (/domain%list(l)%getf%size,size(field,2)/) )
!                 stkpos = stkpos+nwords
                 do j = 1,size(field,2)
                    do i = domain%list(l)%getf%begin,domain%list(l)%getf%end
                       stkpos = stkpos + 1
                       field(i,j) = work(stkpos)
                    end do
                 end do
                 if( debug )write( stderr, '(a,3i4/(8f9.6))' ) 'getf: pe, is, ie=', &
                        pe, domain%list(l)%getf%begin,domain%list(l)%getf%end, &
                      field(domain%list(l)%getf%begin:domain%list(l)%getf%end,:) 
                 call mpp_set_active_domain( domain, end=domain%list(l)%getf%end )
             end if
         end if
      end do
      call mpp_sync_self()

      return
    end subroutine MPP_UPDATE_DOMAINS_2D_

    subroutine MPP_UPDATE_DOMAINS_3D_( field, domain, flags )
!updates data domain of 3D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:)
      type(domain1D), intent(inout), target :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field2D(size(field,1),size(field,2)*size(field,3))
#ifdef use_CRI_pointers
      pointer( ptr, field2D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call MPP_UPDATE_DOMAINS_2D_( field2D, domain, flags )
      return
    end subroutine MPP_UPDATE_DOMAINS_3D_

    subroutine MPP_UPDATE_DOMAINS_4D_( field, domain, flags )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:)
      type(domain1D), intent(inout), target :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field2D(size(field,1),size(field,2)*size(field,3)*size(field,4))
#ifdef use_CRI_pointers
      pointer( ptr, field2D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call MPP_UPDATE_DOMAINS_2D_( field2D, domain, flags )
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_

    subroutine MPP_UPDATE_DOMAINS_5D_( field, domain, flags )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:,:)
      type(domain1D), intent(inout), target :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field2D(size(field,1),size(field,2)*size(field,3)*size(field,4)*size(field,5))
#ifdef use_CRI_pointers
      pointer( ptr, field2D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call MPP_UPDATE_DOMAINS_2D_( field2D, domain, flags )
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_
