    subroutine MPP_UPDATE_DOMAINS_2D_( field, domain, flags )
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:)
      type(domain2D), intent(inout), target :: domain
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
      type(domain2D), intent(inout), target :: domain
      MPP_TYPE_, intent(inout) :: field(domain%x%data%begin:,domain%y%data%begin:,:)
      integer, intent(in), optional :: flags_in

      integer :: is, ie, js, je, i, j, k, flags
      integer :: stkpos, stkpos0, nwords
      logical :: mustput, mustget
      character(len=8) :: text
      MPP_TYPE_ :: work(size(field)*2)
      integer :: words_per_long
#ifdef use_CRI_pointers
      pointer( ptr, work )
      ptr = LOC(mpp_domains_stack)
#else
      equivalence( work, mpp_domains_stack )
#endif

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: You must first call mpp_domains_init.' )

      flags=XUPDATE+YUPDATE
      if( PRESENT(flags_in) )flags = flags_in

      words_per_long = size(transfer(work(1),mpp_domains_stack))
!put west and east halos
      stkpos = 0
      do i = 0,size(domain%x%list)-1
         mustput = .FALSE.
         stkpos0 = stkpos
         if( BTEST(flags,0) .AND. domain%x%list(i)%mustputf )then
             mustput = .TRUE.
             nwords = domain%x%list(i)%putf%size * domain%y%active%size * size(field,3)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,5i4,i5)')'->W: pe, (beg,end), is,ie,js,je,k,   to=', pe, stkpos+1, stkpos+nwords, &
                      domain%x%list(i)%putf%begin, domain%x%list(i)%putf%end, domain%y%active%begin, domain%y%active%end, &
                      size(field,3), domain%x%list(i)%pe
             end if
             if( (stkpos+nwords)*words_per_long.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )stkpos+nwords
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
!is reshape better than pack? 
             work(stkpos+1:stkpos+nwords) = RESHAPE( field(domain%x%list(i)%putf%begin:domain%x%list(i)%putf%end, &
                                                        domain%y%active%begin:domain%y%active%end,:), (/nwords/) )
             stkpos = stkpos+nwords
         end if
         if( BTEST(flags,1) .AND. domain%x%list(i)%mustputb )then
             mustput = .TRUE.
             nwords = domain%x%list(i)%putb%size * domain%y%active%size * size(field,3)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,5i4,i5)')'->E: pe, (beg,end), is,ie,js,je,k,   to=', pe, stkpos+1, stkpos+nwords, &
                      domain%x%list(i)%putb%begin, domain%x%list(i)%putb%end, domain%y%active%begin, domain%y%active%end, &
                      size(field,3), domain%x%list(i)%pe
             end if
             if( (stkpos+nwords)*words_per_long.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )stkpos+nwords
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
!is reshape better than pack? 
             work(stkpos+1:stkpos+nwords) = RESHAPE( field(domain%x%list(i)%putb%begin:domain%x%list(i)%putb%end, &
                                                        domain%y%active%begin:domain%y%active%end,:), (/nwords/) )
             stkpos = stkpos+nwords
         end if
         if( mustput )call mpp_send( work(stkpos0+1:stkpos), stkpos-stkpos0, domain%x%list(i)%pe )
         if( debug )write( stderr, '(a,i4/(8f9.6))' )'pe, work=', pe, work(stkpos0+1:stkpos)
      end do
!get west and east halos
      do i = 0,size(domain%x%list)-1
         stkpos0 = stkpos
         mustget = .FALSE.
         if( BTEST(flags,0) .AND. domain%x%list(i)%mustgetb )then
             mustget = .TRUE.
             nwords = domain%x%list(i)%getb%size * domain%y%active%size * size(field,3)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,5i4,i5)')'<-W: pe, (beg,end), is,ie,js,je,k, from=', pe, stkpos+1, stkpos+nwords, &
                      domain%x%list(i)%getb%begin, domain%x%list(i)%getb%end, domain%y%active%begin, domain%y%active%end, &
                      size(field,3), domain%x%list(i)%pe
             end if
             if( (stkpos+nwords)*words_per_long.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )stkpos+nwords
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             stkpos = stkpos+nwords
         end if
         if( BTEST(flags,1) .AND. domain%x%list(i)%mustgetf )then
             mustget = .TRUE.
             nwords = domain%x%list(i)%getf%size * domain%y%active%size * size(field,3)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,5i4,i5)')'<-E: pe, (beg,end), is,ie,js,je,k, from=', pe, stkpos+1, stkpos+nwords, &
                      domain%x%list(i)%getf%begin, domain%x%list(i)%getf%end, domain%y%active%begin, domain%y%active%end, &
                      size(field,3), domain%x%list(i)%pe
             end if
             if( (stkpos+nwords)*words_per_long.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )stkpos+nwords
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             stkpos = stkpos+nwords
         end if
         if( mustget )then
             call mpp_recv( work(stkpos0+1:stkpos), stkpos-stkpos0, domain%x%list(i)%pe )
             if( debug )write( stderr, '(a,i4/(8f9.6))' )'pe, work=', pe, work(stkpos0+1:stkpos)
             stkpos = stkpos0
             if( BTEST(flags,0) .AND. domain%x%list(i)%mustgetb )then
!is reshape better than unpack? 
                 field(domain%x%list(i)%getb%begin:domain%x%list(i)%getb%end,domain%y%active%begin:domain%y%active%end,:) = &
                      RESHAPE( work(stkpos+1:stkpos+nwords), (/domain%x%list(i)%getb%size,domain%y%active%size,size(field,3)/) )
                 if( debug )write( stderr, '(a,5i4/(8f9.6))' ) 'getb: pe, is, ie, js, je=', &
                        pe, domain%x%list(i)%getb%begin,domain%x%list(i)%getb%end,domain%y%active%begin,domain%y%active%end, &
                      field(domain%x%list(i)%getb%begin:domain%x%list(i)%getb%end,domain%y%active%begin:domain%y%active%end,:) 
                 stkpos = stkpos+nwords
                 call mpp_set_active_domain( domain, xbegin=domain%x%list(i)%getb%begin )
             end if
             if( BTEST(flags,1) .AND. domain%x%list(i)%mustgetf )then
!is reshape better than unpack? 
                 field(domain%x%list(i)%getf%begin:domain%x%list(i)%getf%end,domain%y%active%begin:domain%y%active%end,:) = &
                      RESHAPE( work(stkpos+1:stkpos+nwords), (/domain%x%list(i)%getf%size,domain%y%active%size,size(field,3)/) )
                 if( debug )write( stderr, '(a,5i4/(8f9.6))' ) 'getf: pe, is, ie, js, je==', &
                        pe, domain%x%list(i)%getf%begin,domain%x%list(i)%getf%end,domain%y%active%begin,domain%y%active%end, &
                      field(domain%x%list(i)%getf%begin:domain%x%list(i)%getf%end,domain%y%active%begin:domain%y%active%end,:) 
                 stkpos = stkpos+nwords
                 call mpp_set_active_domain( domain, xend=domain%x%list(i)%getf%end )
             end if
         end if
      end do
      call mpp_sync_self( domain%x%list(:)%pe )

!put south and north halos
      stkpos = 0
      do i = 0,size(domain%y%list)-1
         stkpos0 = stkpos
         mustput = .FALSE.
         if( BTEST(flags,2) .AND. domain%y%list(i)%mustputf )then 
             mustput = .TRUE.
             nwords = domain%x%active%size * domain%y%list(i)%putf%size * size(field,3)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,5i4,i5)')'->S: pe, (beg,end), is,ie,js,je,k,   to=', pe, stkpos+1, stkpos+nwords, &
                      domain%x%active%begin, domain%x%active%end, domain%y%list(i)%putf%begin, domain%y%list(i)%putf%end, &
                      size(field,3), domain%y%list(i)%pe
             end if
             if( (stkpos+nwords)*words_per_long.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )stkpos+nwords
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
!is reshape better than pack? 
             work(stkpos+1:stkpos+nwords) = RESHAPE( field(domain%x%active%begin:domain%x%active%end, &
                  domain%y%list(i)%putf%begin:domain%y%list(i)%putf%end,:), (/nwords/) )
             stkpos = stkpos+nwords
         end if
         if( BTEST(flags,3) .AND. domain%y%list(i)%mustputb )then
             mustput = .TRUE.
             nwords = domain%x%active%size * domain%y%list(i)%putb%size * size(field,3)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,5i4,i5)')'->N: pe, (beg,end), is,ie,js,je,k,   to=', pe, stkpos+1, stkpos+nwords, &
                      domain%x%active%begin, domain%x%active%end, domain%y%list(i)%putb%begin, domain%y%list(i)%putb%end, &
                      size(field,3), domain%y%list(i)%pe
             end if
             if( (stkpos+nwords)*words_per_long.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )stkpos+nwords
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
!is reshape better than pack? 
             work(stkpos+1:stkpos+nwords) = RESHAPE( field(domain%x%active%begin:domain%x%active%end, &
                  domain%y%list(i)%putb%begin:domain%y%list(i)%putb%end,:), (/nwords/) )
             stkpos = stkpos+nwords
         end if
         if( mustput )call mpp_send( work(stkpos0+1:stkpos), stkpos-stkpos0, domain%y%list(i)%pe )
         if( debug )write( stderr, '(a,i4/(8f9.6))' )'pe, work=', pe, work(stkpos0+1:stkpos)
      end do
!get south and north halos
      do i = 0,size(domain%y%list)-1
         stkpos0 = stkpos
         mustget = .FALSE.
         if( BTEST(flags,2) .AND. domain%y%list(i)%mustgetb )then
             mustget = .TRUE.
             nwords = domain%x%active%size * domain%y%list(i)%getb%size * size(field,3)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,5i4,i5)')'<-S: pe, (beg,end), is,ie,js,je,k, from=', pe, stkpos+1, stkpos+nwords, &
                      domain%x%active%begin, domain%x%active%end, domain%y%list(i)%getb%begin, domain%y%list(i)%getb%end, &
                      size(field,3), domain%y%list(i)%pe
             end if
             if( (stkpos+nwords)*words_per_long.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )stkpos+nwords
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             stkpos = stkpos+nwords
         end if
         if( BTEST(flags,3) .AND. domain%y%list(i)%mustgetf )then
             mustget = .TRUE.
             nwords = domain%x%active%size * domain%y%list(i)%getf%size * size(field,3)
             if( debug )then
                 write(stderr,'(a,i3,x,2i4,x,5i4,i5)')'<-N: pe, (beg,end), is,ie,js,je,k, from=', pe, stkpos+1, stkpos+nwords, &
                      domain%x%active%begin, domain%x%active%end, domain%y%list(i)%getf%begin, domain%y%list(i)%getf%end, &
                      size(field,3), domain%y%list(i)%pe
             end if
             if( (stkpos+nwords)*words_per_long.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )stkpos+nwords
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             stkpos = stkpos+nwords
         end if
         if( mustget )then
             call mpp_recv( work(stkpos0+1:stkpos), stkpos-stkpos0, domain%y%list(i)%pe )
             if( debug )write( stderr, '(a,i4/(8f9.6))' )'pe, work=', pe, work(stkpos0+1:stkpos)
             stkpos = stkpos0
             if( BTEST(flags,2) .AND. domain%y%list(i)%mustgetb )then
!is reshape better than unpack? 
                 field(domain%x%active%begin:domain%x%active%end,domain%y%list(i)%getb%begin:domain%y%list(i)%getb%end,:) = &
                      RESHAPE( work(stkpos+1:stkpos+nwords), (/domain%x%active%size,domain%y%list(i)%getb%size,size(field,3)/) )
                 if( debug )write( stderr, '(a,5i4/(8f9.6))' ) 'getb: pe, is, ie, js, je==', &
                        pe, domain%x%active%begin,domain%x%active%end,domain%y%list(i)%getb%begin,domain%y%list(i)%getb%end, &
                      field(domain%x%active%begin:domain%x%active%end,domain%y%list(i)%getb%begin:domain%y%list(i)%getb%end,:) 
                 stkpos = stkpos+nwords
                 call mpp_set_active_domain( domain, ybegin=domain%y%list(i)%getb%begin )
             end if
             if( BTEST(flags,3) .AND. domain%y%list(i)%mustgetf )then
!is reshape better than unpack? 
                 field(domain%x%active%begin:domain%x%active%end,domain%y%list(i)%getf%begin:domain%y%list(i)%getf%end,:) = &
                      RESHAPE( work(stkpos+1:stkpos+nwords), (/domain%x%active%size,domain%y%list(i)%getf%size,size(field,3)/) )
                 if( debug )write( stderr, '(a,5i4/(8f9.6))' ) 'getf: pe, is, ie, js, je=', &
                        pe, domain%x%active%begin,domain%x%active%end,domain%y%list(i)%getf%begin,domain%y%list(i)%getf%end, &
                      field(domain%x%active%begin:domain%x%active%end,domain%y%list(i)%getf%begin:domain%y%list(i)%getf%end,:) 
                 stkpos = stkpos+nwords
                 call mpp_set_active_domain( domain, yend=domain%y%list(i)%getf%end )
             end if
         end if
      end do
      call mpp_sync_self( domain%y%list(:)%pe )

      return
    end subroutine MPP_UPDATE_DOMAINS_3D_

    subroutine MPP_UPDATE_DOMAINS_4D_( field, domain, flags )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:)
      type(domain2D), intent(inout), target :: domain
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
      type(domain2D), intent(inout), target :: domain
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
