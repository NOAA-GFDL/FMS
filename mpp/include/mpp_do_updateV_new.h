    subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_V_(f_addrsx,f_addrsy,isize,jsize,ke,d_type,domain,flags,gridtype,dc_handle)
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND), intent(in)         :: f_addrsx(:), f_addrsy(:)
      integer,            intent(in)         :: isize, jsize, ke
      MPP_TYPE_, intent(in)                  :: d_type  ! creates unique interface
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      type(DomainCommunicator2D),pointer,optional :: dc_handle

      type(DomainCommunicator2D),pointer,save :: dch_x =>NULL(), d_comm =>NULL()
      MPP_TYPE_ :: fieldx(domain%x%data%begin:domain%x%data%begin+isize-1,domain%y%data%begin:domain%y%data%begin+jsize-1,ke)
      pointer(ptr_fieldx, fieldx)
      MPP_TYPE_ :: fieldy(domain%x%data%begin:domain%x%data%begin+isize-1,domain%y%data%begin:domain%y%data%begin+jsize-1,ke)
      pointer(ptr_fieldy, fieldy)
      integer :: update_flags, gridtype_temp
      integer :: l_size, i,j,k,l,n, is, ie, js, je, pos
      integer :: ioff, joff
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer( ptr, buffer )
      integer :: buffer_pos, msgsize, wordlen
      logical :: complete(size(f_addrsy(:)))
      character(len=8) :: text


!for all gridtypes
      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) ) then
          update_flags = flags
          ! The following test is so that SCALAR_PAIR can be used alone with the
          ! same default update pattern as without.
          if (BTEST(update_flags,SCALAR_BIT)) then
            if (.NOT.(BTEST(update_flags,WEST) .OR. BTEST(update_flags,EAST) &
                 .OR. BTEST(update_flags,NORTH) .OR. BTEST(update_flags,SOUTH))) &
              update_flags = update_flags + XUPDATE+YUPDATE   !default with SCALAR_PAIR
          end if 
      end if  

!gridtype
      grid_offset_type = AGRID
      if( PRESENT(gridtype) )then
          if( gridtype.NE.AGRID .AND. &
              gridtype.NE.BGRID_NE .AND. gridtype.NE.BGRID_SW .AND. &
              gridtype.NE.CGRID_NE .AND. gridtype.NE.CGRID_SW ) &
               call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: gridtype must be one of AGRID|BGRID_NE|BGRID_SW|CGRID_NE|CGRID_SW.' )
!grid_offset_type used by update domains to determine shifts.
          grid_offset_type = gridtype
          call compute_overlaps(domain)
          if( grid_offset_type.NE.domain%gridtype ) &
               call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: gridtype cannot be changed during run.' )
      end if   
!need to add code for EWS boundaries
      if( BTEST(domain%fold,WEST) .AND. BTEST(update_flags,WEST) ) &
           call mpp_error( FATAL, 'velocity stencil not yet active for WEST fold, contact author.' )
      if( BTEST(domain%fold,EAST) .AND. BTEST(update_flags,EAST) ) &
           call mpp_error( FATAL, 'velocity stencil not yet active for EAST fold, contact author.' )
      if( BTEST(domain%fold,SOUTH) .AND. BTEST(update_flags,SOUTH) ) &
           call mpp_error( FATAL, 'velocity stencil not yet active for SOUTH fold, contact author.' )

      l_size = size(f_addrsx(:))
      complete(1:l_size-1)=.false.; complete(l_size)=.true.
      d_comm =>NULL(); if(PRESENT(dc_handle))d_comm =>dc_handle
      if (domain%gridtype.EQ.CGRID_NE) then
!On a CGRID, the x- and y-components are staggered differently.  Passing the
!x-component as though it were at the tracer point allows this routine to
!zonally shift it into place.
        gridtype_temp = grid_offset_type ; grid_offset_type = AGRID
!       dch_x =>NULL(); if(ASSOCIATED(d_comm%dch_x))dch_x =>d_comm%dch_x
        dch_x =>NULL()
        if(ASSOCIATED(d_comm))then
           if(ASSOCIATED(d_comm%dch_x))dch_x =>d_comm%dch_x
        endif
        do l=1,l_size
          ptr_fieldx = f_addrsx(l)
          call mpp_update_domains( fieldx, domain, flags, complete=complete(l), dc_handle=dch_x )
        end do
        grid_offset_type = gridtype_temp
        do l=1,l_size
          ptr_fieldy = f_addrsy(l)
          call mpp_update_domains( fieldy, domain, flags, complete=complete(l), dc_handle=d_comm )
        end do
        if(PRESENT(dc_handle))then
          d_comm%dch_x =>dch_x
          dc_handle =>d_comm 
        endif
      else
        do l=1,l_size
          ptr_fieldx = f_addrsx(l)
          ptr_fieldy = f_addrsy(l)
          call mpp_update_domains( fieldx, domain, flags, complete=.false. )
          call mpp_update_domains( fieldy, domain, flags, complete=complete(l), dc_handle=dc_handle )
        end do
      endif

!     do l=1,l_size
!       ptr_fieldx = f_addrsx(l)
!       ptr_fieldy = f_addrsy(l)
!       write(stdout(),*) 'fieldx chksum=',mpp_chksum(fieldx)
!       write(stdout(),*) 'fieldy chksum=',mpp_chksum(fieldy)
!     end do

      ptr = LOC(mpp_domains_stack)
      wordlen = size(TRANSFER(buffer(1),mpp_domains_stack))
      buffer_pos = 0
      call mpp_get_global_domain( domain, xsize=ioff, ysize=joff )
!northern boundary fold
      if( BTEST(domain%fold,NORTH) .AND. BTEST(update_flags,NORTH) )then
          js = domain%y%global%end + 1
          je = domain%y%data%end
          if( je.GE.js )then
!on offset grids, we need to move data leftward by one point
              pos = domain%x%pos - 1 !the one on your left
              if( pos.GE.0 )then
                  is = domain%x%list(pos)%data%end+1; ie=is
              else if( domain%x%cyclic )then
                  pos = pos + size(domain%x%list(:))
                  is = domain%x%list(pos)%data%end+1 - ioff; ie=is
              else
                  is=1; ie=0
              end if
              n = buffer_pos
              if( ie.EQ.is )then
                  msgsize = (je-js+1)*ke*2*l_size !only half this on CGRID actually
                  mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
                  if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                      write( text,'(i8)' )mpp_domains_stack_hwm
                      call mpp_error( FATAL, 'MPP_UPDATE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                           //trim(text)//') from all PEs.' )
                  end if   
                  select case(grid_offset_type)
                  case(BGRID_NE)
                    do l=1,l_size
                      ptr_fieldx = f_addrsx(l)
                      ptr_fieldy = f_addrsy(l)   
                      do k = 1,ke
                         do j = js,je
                            n = n + 2
                            buffer(n-1) = fieldx(is,j,k)
                            buffer(n  ) = fieldy(is,j,k)
                         end do
                      end do
                    end do
                    call mpp_send( buffer(buffer_pos+1), plen=n, to_pe=domain%x%list(pos)%pe )
                    buffer_pos = buffer_pos + n
                  case(CGRID_NE)
                    ! When the halo to the east has not been filled, copy at the end of the compute domain, not
                    ! the data domain.
                    if (.not.BTEST(update_flags,EAST)) is = is - domain%x%data%end + domain%x%compute%end
                    do l=1,l_size
                      ptr_fieldx = f_addrsx(l)
                      do k = 1,ke
                         do j = js,je
                            n = n + 1
                            buffer(n) = fieldx(is,j,k)
                         end do
                      end do
                    end do
                    call mpp_send( buffer(buffer_pos+1), plen=n, to_pe=domain%x%list(pos)%pe )
                    buffer_pos = buffer_pos + n
                  end select
!receive data at x%data%end
                  pos = domain%x%pos + 1 !the one on your right
                  if( pos.LT.size(domain%x%list(:)) )then
                      n = (je-js+1)*ke
                  else if( domain%x%cyclic )then
                      pos = pos - size(domain%x%list(:))
                      n = (je-js+1)*ke
                  else
                      n = 0
                  end if
                  n = n*l_size
                  if( n.GT.0 )then
                      mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+n)*wordlen*l_size )
                      if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                          write( text,'(i8)' )mpp_domains_stack_hwm
                          call mpp_error( FATAL, 'MPP_UPDATE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                               //trim(text)//') from all PEs.' )
                      end if   
                      select case(grid_offset_type)
                      case(BGRID_NE)
                        do l=1,l_size
                          ptr_fieldx = f_addrsx(l)
                          ptr_fieldy = f_addrsy(l)   
                          do k = 1,ke
                             do j = js,je
                                do i = domain%x%data%begin,domain%x%data%end-1
                                   fieldx(i,j,k) = fieldx(i+1,j,k)
                                   fieldy(i,j,k) = fieldy(i+1,j,k)
                                end do
                             end do
                          end do
                        end do
                        n = 2*n
                        call mpp_recv( buffer(buffer_pos+1), glen=n, from_pe=domain%x%list(pos)%pe )
                        i = domain%x%data%end
                        n = buffer_pos
                        do l=1,l_size
                          ptr_fieldx = f_addrsx(l)
                          ptr_fieldy = f_addrsy(l)   
                          do k = 1,ke
                             do j = js,je
                                n = n + 2
                                fieldx(i,j,k) = buffer(n-1)
                                fieldy(i,j,k) = buffer(n  )
                             end do
                          end do
                        end do
                      case(CGRID_NE)
                        do l=1,l_size
                          ptr_fieldx = f_addrsx(l)
                          ptr_fieldy = f_addrsy(l)   
                          do k = 1,ke
                             do j = js,je
                                do i = domain%x%data%begin,domain%x%data%end-1
                                   fieldx(i,j,k) = fieldx(i+1,j,k)
!                                   fieldy(i,j,k) = fieldy(i+1,j,k)
                                end do
                             end do
                          end do
                        end do
                        call mpp_recv( buffer(buffer_pos+1), glen=n, from_pe=domain%x%list(pos)%pe )
                        i = domain%x%data%end
! If the eastern halo has not been filled, it is the eastern end of the compute domain
! that needs to have the value shifted into it.
                        if (.not.BTEST(update_flags,EAST)) i = domain%x%compute%end
                        n = buffer_pos
                        do l=1,l_size
                          ptr_fieldx = f_addrsx(l)
                          do k = 1,ke
                             do j = js,je
                                n = n + 1
                                fieldx(i,j,k) = buffer(n)
                             end do
                          end do
                        end do
                      end select
                  end if  
              end if  

              if (.NOT.BTEST(update_flags,SCALAR_BIT)) then
!flip the sign if this is a vector
                is = domain%x%data%begin
                ie = domain%x%data%end
                do l=1,l_size
                  ptr_fieldx = f_addrsx(l)
                  ptr_fieldy = f_addrsy(l)   
                  do k = 1,ke
                     do j = js,je
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
                end do
              endif  ! if (.NOT.BTEST(update_flags,SCALAR_BIT))
          end if 

!         do l=1,l_size
!           ptr_fieldx = f_addrsx(l)
!           ptr_fieldy = f_addrsy(l)
!           write(stdout(),*) 'fieldx chksum=',mpp_chksum(fieldx)
!           write(stdout(),*) 'fieldy chksum=',mpp_chksum(fieldy)
!         end do

!eliminate redundant vector data at fold
          j = domain%y%global%end
          if( domain%y%data%begin.LE.j .AND. j.LE.domain%y%data%end )then !fold is within domain
!ship left-half data to right half: on BGRID_NE the x%data%end point is not in mirror domain and must be done separately.
              if( domain%x%pos.LT.(size(domain%x%list(:))+1)/2 )then
                  is = domain%x%data%begin
                  ie = min(domain%x%data%end,(domain%x%global%begin+domain%x%global%end)/2)
                  n = buffer_pos
                  select case(grid_offset_type)
                  case(BGRID_NE)
                    do l=1,l_size
                      ptr_fieldx = f_addrsx(l)
                      ptr_fieldy = f_addrsy(l)   
                      do k = 1,ke
                         do i = is,ie-1
                            n = n + 2
                            buffer(n-1) = fieldx(i,j,k)
                            buffer(n)   = fieldy(i,j,k)
                         end do
                      end do
                    end do
                    call mpp_send( buffer(buffer_pos+1), plen=n-buffer_pos, &
                                   to_pe=domain%x%list(size(domain%x%list(:))-domain%x%pos-1)%pe )
                    buffer_pos = n 
                  case(CGRID_NE)     
                    do l=1,l_size
                      ptr_fieldy = f_addrsy(l)   
                      do k = 1,ke
                         do i = is,ie
                            n = n + 1
                            buffer(n) = fieldy(i,j,k)
                         end do
                      end do
                    end do
                    call mpp_send( buffer(buffer_pos+1), plen=n-buffer_pos, &
                                   to_pe=domain%x%list(size(domain%x%list(:))-domain%x%pos-1)%pe )
                    buffer_pos = n 
                  end select         
              end if  
              if( domain%x%pos.GE.size(domain%x%list(:))/2 )then
                  is = max(domain%x%data%begin,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie = domain%x%data%end
                  select case(grid_offset_type)
                  case(BGRID_NE)
                    n = (ie-is+1)*ke*2*l_size
                    call mpp_recv( buffer(buffer_pos+1), glen=n, from_pe=domain%x%list(size(domain%x%list(:))-domain%x%pos-1)%pe )
                    n = buffer_pos
!get all values except at x%data%end
                    if (BTEST(update_flags,SCALAR_BIT)) then ! Do not change the signs if this is a pair of scalar fields.
                      do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)
                        do k = 1,ke
                           do i = ie-1,is,-1
                              n = n + 2
                              fieldx(i,j,k) = buffer(n-1)
                              fieldy(i,j,k) = buffer(n)
                           end do
                        end do
                      end do
                    else
                      do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                           do i = ie-1,is,-1
                              n = n + 2
                              fieldx(i,j,k) = -buffer(n-1)
                              fieldy(i,j,k) = -buffer(n)
                           end do
                        end do
                      end do
                    endif
!now get the value at domain%x%data%end
                    pos = domain%x%pos - 1
                    if( pos.GE.size(domain%x%list(:))/2 )then
                        i = domain%x%list(pos)%data%end
                        buffer_pos = n
                        do l=1,l_size
                          ptr_fieldx = f_addrsx(l)
                          ptr_fieldy = f_addrsy(l)   
                          do k = 1,ke
                             n = n + 2
                             buffer(n-1) = fieldx(i,j,k)
                             buffer(n  ) = fieldy(i,j,k)
                          end do
                        end do
                        call mpp_send( buffer(buffer_pos+1), plen=n-buffer_pos, to_pe=domain%x%list(pos)%pe )
                        buffer_pos = n
                    end if
                    pos = domain%x%pos + 1
                    if( pos.LT.size(domain%x%list(:)) )then
                        n = ke*2*l_size
                        call mpp_recv( buffer(buffer_pos+1), glen=n, from_pe=domain%x%list(pos)%pe )
                        n = buffer_pos
                        i = domain%x%data%end
                        do l=1,l_size
                          ptr_fieldx = f_addrsx(l)
                          ptr_fieldy = f_addrsy(l)   
                          do k = 1,ke
                             n = n + 2
                             fieldx(i,j,k) = buffer(n-1)
                             fieldy(i,j,k) = buffer(n  )
                          end do
                        end do
                    end if 
                  case(CGRID_NE)
                    n = (ie-is+1)*ke*l_size
                    call mpp_recv( buffer(buffer_pos+1), glen=n, from_pe=domain%x%list(size(domain%x%list(:))-domain%x%pos-1)%pe )
                    n = buffer_pos
                    if (BTEST(update_flags,SCALAR_BIT)) then ! Do not change the signs if this is a pair of scalar fields.
                      do l=1,l_size
                        ptr_fieldy = f_addrsy(l)
                        do k = 1,ke
                           do i = ie,is,-1
                              n = n + 1
                              fieldy(i,j,k) = buffer(n)
                           end do
                        end do
                      end do
                    else
                      do l=1,l_size
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                           do i = ie,is,-1
                              n = n + 1
                              fieldy(i,j,k) = -buffer(n)
                           end do
                        end do
                      end do
                    endif
                  end select
              end if  
!poles set to 0: BGRID only
              if( grid_offset_type.EQ.BGRID_NE.and. .not. BTEST(update_flags,SCALAR_BIT) )then
                  do i = domain%x%global%begin-1,domain%x%global%end,(domain%x%global%begin+domain%x%global%end)/2
                     if( domain%x%data%begin.LE.i .AND. i.LE.domain%x%data%end )then
                        do l=1,l_size
                          ptr_fieldx = f_addrsx(l)
                          ptr_fieldy = f_addrsy(l)   
                          do k = 1,ke
                            fieldx(i,j,k) = 0.
                            fieldy(i,j,k) = 0.
                          end do
                        end do
                     end if 
                  end do 
              end if 
!these last three code blocks correct an error where the data in your halo coming from other half may have the wrong sign
!off west edge
              select case(grid_offset_type)
              case(BGRID_NE)
                  is = domain%x%global%begin - 1
                  if( is.GT.domain%x%data%begin )then
                      if( 2*is-domain%x%data%begin.GT.domain%x%data%end ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: BGRID_NE west edge ubound error.' )
                      do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                          do i = domain%x%data%begin,is-1
                            fieldx(i,j,k) = fieldx(2*is-i,j,k)
                            fieldy(i,j,k) = fieldy(2*is-i,j,k)
                          end do
                        end do
                      end do
                  end if 
              case(CGRID_NE)
                  is = domain%x%global%begin
                  if( is.GT.domain%x%data%begin )then
                      if( 2*is-domain%x%data%begin-1.GT.domain%x%data%end ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: CGRID_NE west edge ubound error.' )
                      do l=1,l_size
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                          do i = domain%x%data%begin,is-1
                            fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                          end do
                        end do
                      end do
                  end if 
              end select
!right of midpoint
              is = (domain%x%global%begin+domain%x%global%end)/2
              if( domain%x%compute%begin.LE.is .AND. is.LT.domain%x%data%end )then
                  select case(grid_offset_type)
                  case(BGRID_NE)
                      ie = domain%x%data%end
                      if( 2*is-ie.LT.domain%x%data%begin )ie = ie - 1
                      if( 2*is-ie.LT.domain%x%data%begin ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: BGRID_NE midpoint lbound error.' )
                      if (BTEST(update_flags,SCALAR_BIT)) then
                        do l=1,l_size
                          ptr_fieldx = f_addrsx(l)
                          ptr_fieldy = f_addrsy(l)
                          do k = 1,ke
                            do i = is+1,ie
                              fieldx(i,j,k) = fieldx(2*is-i,j,k)
                              fieldy(i,j,k) = fieldy(2*is-i,j,k)
                            end do
                          end do
                        end do
                      else
                        do l=1,l_size
                          ptr_fieldx = f_addrsx(l)
                          ptr_fieldy = f_addrsy(l)   
                          do k = 1,ke
                            do i = is+1,ie
                              fieldx(i,j,k) = -fieldx(2*is-i,j,k)
                              fieldy(i,j,k) = -fieldy(2*is-i,j,k)
                            end do
                          end do
                        end do
                      endif
                  case(CGRID_NE)
                      if( 2*is-domain%x%data%end+1.LT.domain%x%data%begin ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: CGRID_NE midpoint lbound error.' )
                      if (BTEST(update_flags,SCALAR_BIT)) then
                        do l=1,l_size
                          ptr_fieldy = f_addrsy(l)
                          do k = 1,ke
                            do i = is+1,domain%x%data%end
                              fieldy(i,j,k) = fieldy(2*is-i+1,j,k)
                            end do
                          end do
                        end do
                      else
                        do l=1,l_size
                          ptr_fieldy = f_addrsy(l)   
                          do k = 1,ke
                            do i = is+1,domain%x%data%end
                              fieldy(i,j,k) = -fieldy(2*is-i+1,j,k)
                            end do
                          end do
                        end do
                      endif
                  end select
              end if  
!off east edge    
              is = domain%x%global%end
              if( is.LT.domain%x%data%end )then
                  select case(grid_offset_type)
                  case(BGRID_NE)
                      if( 2*is-domain%x%data%end.LT.domain%x%data%begin ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: BGRID_NE east edge lbound error.' )
                      do l=1,l_size
                        ptr_fieldx = f_addrsx(l)
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                          do i = is+1,domain%x%data%end
                            fieldx(i,j,k) = fieldx(2*is-i,j,k)
                            fieldy(i,j,k) = fieldy(2*is-i,j,k)
                          end do
                        end do
                      end do
                  case(CGRID_NE)
                      if( 2*is-domain%x%data%end+1.LT.domain%x%data%begin ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: CGRID_NE east edge lbound error.' )
                      do l=1,l_size
                        ptr_fieldy = f_addrsy(l)   
                        do k = 1,ke
                          do i = is+1,domain%x%data%end
                            fieldy(i,j,k) = fieldy(2*is-i+1,j,k)
                          end do
                        end do
                      end do
                  end select
              end if  
          end if  
!         do l=1,l_size
!           ptr_fieldx = f_addrsx(l)
!           ptr_fieldy = f_addrsy(l)
!           write(stdout(),*) 'fieldx chksum=',mpp_chksum(fieldx)
!           write(stdout(),*) 'fieldy chksum=',mpp_chksum(fieldy)
!         end do
      end if  
          
      grid_offset_type = AGRID  !reset
      call mpp_sync_self()
      return
    end subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_V_
