    subroutine MPP_UPDATE_DOMAINS_2D_( field, domain, flags, type )
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:)
      type(domain2D), intent(inout), target :: domain
      integer, intent(in), optional :: flags, type
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call mpp_update_domains( field3D, domain, flags, type )
      return
    end subroutine MPP_UPDATE_DOMAINS_2D_

    subroutine MPP_UPDATE_DOMAINS_3D_( field, domain, flags_in, type )
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(inout), target :: domain
      MPP_TYPE_, intent(inout) :: field(domain%x%data%begin:,domain%y%data%begin:,:)
      integer, intent(in), optional :: flags_in, type

      integer :: is, ie, js, je, i, j, k, l, flags
      integer :: bufpos, lastpos, nwords
      logical :: mustput, mustget, vectorcomp
      character(len=8) :: text
      MPP_TYPE_ :: work(size(field)*2)
      integer :: words_per_long, mpp_domains_stack_pos
#ifdef use_CRI_pointers
      pointer( ptr, work )
      ptr = LOC(mpp_domains_stack)
#else
      equivalence( work, mpp_domains_stack )
#endif

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: You must first call mpp_domains_init.' )

      flags=XUPDATE+YUPDATE
      if( PRESENT(flags_in) )flags = flags_in
      vectorcomp = .FALSE.
      if( PRESENT(type) )then
          if( type.NE.VECTOR_COMPONENT ) &
               call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: the only value of type currently allowed is VECTOR_COMPONENT.' )
#ifdef MPP_TYPE_IS_LOGICAL_
          call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: a logical variable cannot have a VECTOR_COMPONENT type.' )
#endif
          vectorcomp = .TRUE.
      end if

      words_per_long = size(transfer(work(1),mpp_domains_stack))
!if update domains is being called, assume active domain is reduced to compute domain
      call mpp_get_compute_domain( domain, is, ie, js, je )
      call mpp_set_active_domain ( domain, is, ie, js, je )

!put south and north halos
      bufpos = 0
      mpp_domains_stack_pos = 0
      do l = 0,size(domain%y%list)-1
         lastpos = bufpos
         mustput = .FALSE.
         if( BTEST(flags,2) .AND. domain%y%list(l)%mustputf )then 
             mustput = .TRUE.
             nwords = domain%x%active%size * domain%y%list(l)%putf%size * size(field,3)
             mpp_domains_stack_pos = (bufpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             do k = 1,size(field,3)
                do j = domain%y%list(l)%putf%begin,domain%y%list(l)%putf%end
                   do i = domain%x%active%begin,domain%x%active%end
                      bufpos = bufpos + 1
                      work(bufpos) = field(i,j,k)
                   end do
                end do
             end do
         end if
         if( BTEST(flags,3) .AND. domain%y%list(l)%mustputb )then
             mustput = .TRUE.
             nwords = domain%x%active%size * domain%y%list(l)%putb%size * size(field,3)
             mpp_domains_stack_pos = (bufpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             do k = 1,size(field,3)
                do j = domain%y%list(l)%putb%begin,domain%y%list(l)%putb%end
                   do i = domain%x%active%begin,domain%x%active%end
                      bufpos = bufpos + 1
                      work(bufpos) = field(i,j,k)
                   end do
                end do
             end do
         end if
         if( mustput )then
             if( debug )write( stderr, '(a,3i4,i8)' )'YSEND from, to, l, length=', pe, domain%y%list(l)%pe, l, bufpos-lastpos
             call mpp_send( work(lastpos+1), bufpos-lastpos, domain%y%list(l)%pe )
         end if
      end do
!get south and north halos
      do l = 0,size(domain%y%list)-1
         lastpos = bufpos
         mustget = .FALSE.
         if( BTEST(flags,2) .AND. domain%y%list(l)%mustgetb )then
             mustget = .TRUE.
             nwords = domain%x%active%size * domain%y%list(l)%getb%size * size(field,3)
             mpp_domains_stack_pos = (bufpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             bufpos = bufpos+nwords
             call mpp_set_active_domain( domain, ybegin=domain%y%data%begin )
         end if
         if( BTEST(flags,3) .AND. domain%y%list(l)%mustgetf )then
             mustget = .TRUE.
             nwords = domain%x%active%size * domain%y%list(l)%getf%size * size(field,3)
             mpp_domains_stack_pos = (bufpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             bufpos = bufpos+nwords
             call mpp_set_active_domain( domain, yend=domain%y%data%end )
         end if
         if( mustget )then
             if( debug )write( stderr, '(a,3i4,i8)' )'YRECV to, from, l, length=', pe, domain%y%list(l)%pe, l, bufpos-lastpos
             call mpp_recv( work(lastpos+1), bufpos-lastpos, domain%y%list(l)%pe )
             bufpos = lastpos
             if( BTEST(flags,2) .AND. domain%y%list(l)%mustgetb )then
                 if( domain%y%list(l)%folded )then
                     if( vectorcomp )then
                         do k = 1,size(field,3)
                            do j = domain%y%list(l)%getb%end,domain%y%list(l)%getb%begin,-1
                               do i = domain%x%active%end,domain%x%active%begin,-1
                                  bufpos = bufpos + 1
#ifndef MPP_TYPE_IS_LOGICAL_
                                  field(i,j,k) = -work(bufpos)
#endif
                               end do
                            end do
                         end do
                     else
                         do k = 1,size(field,3)
                            do j = domain%y%list(l)%getb%end,domain%y%list(l)%getb%begin,-1
                               do i = domain%x%active%end,domain%x%active%begin,-1
                                  bufpos = bufpos + 1
                                  field(i,j,k) =  work(bufpos)
                               end do
                            end do
                         end do
                     end if
                 else
                     do k = 1,size(field,3)
                        do j = domain%y%list(l)%getb%begin,domain%y%list(l)%getb%end
                           do i = domain%x%active%begin,domain%x%active%end
                              bufpos = bufpos + 1
                              field(i,j,k) = work(bufpos)
                           end do
                        end do
                     end do
                 end if
                 call mpp_set_active_domain( domain, ybegin=domain%y%list(l)%getb%begin )
             end if
             if( BTEST(flags,3) .AND. domain%y%list(l)%mustgetf )then
                 if( domain%y%list(l)%folded )then
                     if( vectorcomp )then
                         do k = 1,size(field,3)
                            do j = domain%y%list(l)%getf%end,domain%y%list(l)%getf%begin,-1
                               do i = domain%x%active%end,domain%x%active%begin,-1
                                  bufpos = bufpos + 1
#ifndef MPP_TYPE_IS_LOGICAL_
                                  field(i,j,k) = -work(bufpos)
#endif
                               end do
                            end do
                         end do
                     else
                         do k = 1,size(field,3)
                            do j = domain%y%list(l)%getf%end,domain%y%list(l)%getf%begin,-1
                               do i = domain%x%active%end,domain%x%active%begin,-1
                                  bufpos = bufpos + 1
                                  field(i,j,k) =  work(bufpos)
                               end do
                            end do
                         end do
                     end if
                 else
                     do k = 1,size(field,3)
                        do j = domain%y%list(l)%getf%begin,domain%y%list(l)%getf%end
                           do i = domain%x%active%begin,domain%x%active%end
                              bufpos = bufpos + 1
                              field(i,j,k) = work(bufpos)
                           end do
                        end do
                     end do
                 end if
                 call mpp_set_active_domain( domain, yend=domain%y%list(l)%getf%end )
             end if
         end if
      end do
      call mpp_sync_self( domain%y%list(:)%pe )

!put west and east halos
      bufpos = 0
      mpp_domains_stack_pos = 0
      do l = 0,size(domain%x%list)-1
         mustput = .FALSE.
         lastpos = bufpos       !end of last send message buffer
         if( BTEST(flags,0) .AND. domain%x%list(l)%mustputf )then
             mustput = .TRUE.
             nwords = domain%x%list(l)%putf%size * domain%y%active%size * size(field,3)
             mpp_domains_stack_pos = (bufpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             do k = 1,size(field,3)
                do j = domain%y%active%begin,domain%y%active%end
                   do i = domain%x%list(l)%putf%begin,domain%x%list(l)%putf%end
                      bufpos = bufpos + 1
                      work(bufpos) = field(i,j,k)
                   end do
                end do
             end do
         end if
         if( BTEST(flags,1) .AND. domain%x%list(l)%mustputb )then
             mustput = .TRUE.
             nwords = domain%x%list(l)%putb%size * domain%y%active%size * size(field,3)
             mpp_domains_stack_pos = (bufpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             do k = 1,size(field,3)
                do j = domain%y%active%begin,domain%y%active%end
                   do i = domain%x%list(l)%putb%begin,domain%x%list(l)%putb%end
                      bufpos = bufpos + 1
                      work(bufpos) = field(i,j,k)
                   end do
                end do
             end do
         end if
         if( mustput )then
             if( debug )write( stderr, '(a,3i4,i8)' )'XSEND from, to, l, length=', pe, domain%x%list(l)%pe, l, bufpos-lastpos
             call mpp_send( work(lastpos+1), bufpos-lastpos, domain%x%list(l)%pe )
         end if
      end do
!get west and east halos
      do l = 0,size(domain%x%list)-1
         lastpos = bufpos
         mustget = .FALSE.
         if( BTEST(flags,0) .AND. domain%x%list(l)%mustgetb )then
             mustget = .TRUE.
             nwords = domain%x%list(l)%getb%size * domain%y%active%size * size(field,3)
             mpp_domains_stack_pos = (bufpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             bufpos = bufpos+nwords
             call mpp_set_active_domain( domain, xbegin=domain%x%data%begin )
         end if
         if( BTEST(flags,1) .AND. domain%x%list(l)%mustgetf )then
             mustget = .TRUE.
             nwords = domain%x%list(l)%getf%size * domain%y%active%size * size(field,3)
             mpp_domains_stack_pos = (bufpos+nwords)*words_per_long
             if( mpp_domains_stack_pos.GT.mpp_domains_stack_size )then
                 write( text, '(i8)' )mpp_domains_stack_pos
                 call mpp_error( FATAL, &
                      'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, mpp_domains_stack_pos )
             bufpos = bufpos+nwords
             call mpp_set_active_domain( domain, xend=domain%x%data%end )
         end if
         if( mustget )then
             if( debug )write( stderr, '(a,3i4,i8)' )'XRECV to, from, l, length=', pe, domain%x%list(l)%pe, l, bufpos-lastpos
             call mpp_recv( work(lastpos+1), bufpos-lastpos, domain%x%list(l)%pe )
             bufpos = lastpos
             if( BTEST(flags,0) .AND. domain%x%list(l)%mustgetb )then
                 if( domain%x%list(l)%folded )then
                     if( vectorcomp )then
                         do k = 1,size(field,3)
                            do j = domain%y%active%end,domain%y%active%begin,-1
                               do i = domain%x%list(l)%getb%end,domain%x%list(l)%getb%begin-1
                                  bufpos = bufpos + 1
#ifndef MPP_TYPE_IS_LOGICAL_
                                  field(i,j,k) = -work(bufpos)
#endif
                               end do
                            end do
                         end do
                     else
                         do k = 1,size(field,3)
                            do j = domain%y%active%end,domain%y%active%begin,-1
                               do i = domain%x%list(l)%getb%end,domain%x%list(l)%getb%begin-1
                                  bufpos = bufpos + 1
                                  field(i,j,k) =  work(bufpos)
                               end do
                            end do
                         end do
                     end if
                 else
                     do k = 1,size(field,3)
                        do j = domain%y%active%begin,domain%y%active%end
                           do i = domain%x%list(l)%getb%begin,domain%x%list(l)%getb%end
                              bufpos = bufpos + 1
                              field(i,j,k) = work(bufpos)
                           end do
                        end do
                     end do
                 end if
                 call mpp_set_active_domain( domain, xbegin=domain%x%list(l)%getb%begin )
             end if
             if( BTEST(flags,1) .AND. domain%x%list(l)%mustgetf )then
                 if( domain%x%list(l)%folded )then
                     if( vectorcomp )then
                         do k = 1,size(field,3)
                            do j = domain%y%active%end,domain%y%active%begin,-1
                               do i = domain%x%list(l)%getf%end,domain%x%list(l)%getf%begin-1
                                  bufpos = bufpos + 1
#ifndef MPP_TYPE_IS_LOGICAL_
                                  field(i,j,k) = -work(bufpos)
#endif
                               end do
                            end do
                         end do
                     else
                         do k = 1,size(field,3)
                            do j = domain%y%active%end,domain%y%active%begin,-1
                               do i = domain%x%list(l)%getf%end,domain%x%list(l)%getf%begin-1
                                  bufpos = bufpos + 1
                                  field(i,j,k) =  work(bufpos)
                               end do
                            end do
                         end do
                     end if
                 else
                     do k = 1,size(field,3)
                        do j = domain%y%active%begin,domain%y%active%end
                           do i = domain%x%list(l)%getf%begin,domain%x%list(l)%getf%end
                              bufpos = bufpos + 1
                              field(i,j,k) = work(bufpos)
                           end do
                        end do
                     end do
                 end if
                 call mpp_set_active_domain( domain, xend=domain%x%list(l)%getf%end )
             end if
         end if
      end do
      call mpp_sync_self( domain%x%list(:)%pe )

      return
    end subroutine MPP_UPDATE_DOMAINS_3D_

    subroutine MPP_UPDATE_DOMAINS_4D_( field, domain, flags, type )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:)
      type(domain2D), intent(inout), target :: domain
      integer, intent(in), optional :: flags, type
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call mpp_update_domains( field3D, domain, flags, type )
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_

    subroutine MPP_UPDATE_DOMAINS_5D_( field, domain, flags, type )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:,:)
      type(domain2D), intent(inout), target :: domain
      integer, intent(in), optional :: flags, type
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_2D_: requires Cray pointers.' )
#endif
      call mpp_update_domains( field3D, domain, flags, type )
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_
