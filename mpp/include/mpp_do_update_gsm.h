    subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_( f_addrs, d_comm, d_type )
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND), intent(in)                 :: f_addrs(:)
      type(DomainCommunicator2D), intent(in)         :: d_comm
      MPP_TYPE_, intent(in)                          :: d_type  ! creates unique interface

      MPP_TYPE_ :: field(d_comm%domain%x%data%begin:d_comm%domain%x%data%begin+d_comm%isize-1, &
                         d_comm%domain%y%data%begin:d_comm%domain%y%data%begin+d_comm%jsize-1,d_comm%ke)
      pointer(ptr_field, field)
      MPP_TYPE_ :: field_in(d_comm%isize_max*d_comm%jsize_max*d_comm%ke)
      pointer( ptr_field_in,field_in )
      real(KIND(d_type)) :: r_field(d_comm%domain%x%data%size,d_comm%domain%y%data%size,d_comm%ke)
      pointer( r_ptr,r_field )
      integer :: i, j, k, l, m, n
      integer :: is, ie, js, je, ke, l_size, istat
      integer :: ii,jj,isizeR,jsizeR
      integer :: list, from_pe
!     integer,save :: done_cnt
!     logical, dimension(0:d_comm%Rlist_size-1) :: done


      ke    = d_comm%ke
      l_size = size(f_addrs(:))

      call mpp_clock_begin(unpk_clock)
!recv
!unpack halos in reverse order
      n = d_comm%Rlist_size
!     done = .false. .OR. .NOT.d_comm%R_do_buf(:)
!     done_cnt = n - count(done)
!     gsm_sync_start = Set_Sync_Logical(.true.)
      call mpp_sync(do_self=.false.)
!     do while(done_cnt /= 0); do list=n-1,0,-1
      do list=n-1,0,-1
         if( .NOT. d_comm%R_do_buf(list) )cycle
!        if(done(list))cycle
!        if( Is_Remote_False(d_comm%sync_start_list(list)) )cycle
!        done(list) = .true.; done_cnt = done_cnt - 1
         isizeR = d_comm%isizeR(list)
         jsizeR = d_comm%jsizeR(list)
         do l=1,l_size  ! loop over number of fields
           ptr_field = f_addrs(l)
           ptr_field_in = d_comm%rem_addrl(l,list)
           do m=8,1,-1
             if( d_comm%do_thisR(1,m,list) )then  ! ne
               is=d_comm%recvis(m,list); ie=d_comm%recvie(m,list)
               js=d_comm%recvjs(m,list); je=d_comm%recvje(m,list)
               if( d_comm%do_thisR(2,m,list) )then  ! folded?
                   do k = 1,ke
                       jj = d_comm%sendjsR(m,list)
                      do j = je,js,-1
                       ii = d_comm%sendisR(m,list)
                         do i = ie,is,-1
                            field(i,j,k) = field_in(ii+(jj-1)*isizeR+(k-1)*isizeR*jsizeR)
                            ii = ii + 1
                         end do
                        jj = jj + 1
                      end do
                   end do
               else
                  do k = 1,ke
                       jj = d_comm%sendjsR(m,list)
                      do j = js,je
                       ii = d_comm%sendisR(m,list)
                         do i = is,ie
                            field(i,j,k) = field_in(ii+(jj-1)*isizeR+(k-1)*isizeR*jsizeR)
                            ii = ii + 1
                         end do
                       jj = jj + 1
                      end do
                   end do
               end if
             end if
           end do
        end do
      end do
!     end do; end do

      if(debug_gsm)then
        do l=1,l_size
          r_ptr = f_addrs(l)
          write(stdout(),*) 'Domain checksum=',mpp_chksum(r_field)
        end do
      endif
 
      call mpp_sync(do_self=.false.)
!     gsm_sync_start = Set_Sync_Logical(.false.)

      call mpp_clock_end(unpk_clock)
      return
    end subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_
