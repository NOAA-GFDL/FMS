    subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_( f_addrs, d_comm, d_type )
      use  mpp_datatype_mod, only: CAFPNTR_TYPE_3D_
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND), intent(in)                 :: f_addrs(:)
      type(DomainCommunicator2D), intent(in), target :: d_comm
      MPP_TYPE_, intent(in)                          :: d_type  ! creates unique interface

      type d_comm_ptr_type
        type(DomainCommunicator2D), pointer :: dc_ptr
      end type
      type(d_comm_ptr_type), allocatable, save :: d_commG[:]
      type(CAFPNTR_TYPE_3D_), allocatable, save :: cafptr(:)[:]
      MPP_TYPE_ :: field(d_comm%domain%x%data%begin:d_comm%domain%x%data%begin+d_comm%isize-1, &
                         d_comm%domain%y%data%begin:d_comm%domain%y%data%begin+d_comm%jsize-1,d_comm%ke)
      pointer(ptr_field, field)
      real(KIND(d_type)) :: r_field(d_comm%domain%x%data%size,d_comm%domain%y%data%size,d_comm%ke)
      pointer( r_ptr,r_field )
      integer :: i, j, k, l, m, n
      integer :: is, ie, js, je, ke, l_size, istat
      integer :: ii,jj
      integer :: list, from_pe, r_idx

      logical, save :: first_time=.true.


      ke    = d_comm%ke
      l_size = size(f_addrs(:))

      if(first_time)then
        first_time = .false.
        ALLOCATE(cafptr(MAX_DOMAIN_FIELDS)[0:*],d_commG[0:*])
      endif

      do l=1,l_size
        ptr_field = f_addrs(l)
        call mpp_associate_caf_field(d_comm%domain%x%data%begin, d_comm%domain%x%data%end, &
                                     d_comm%domain%y%data%begin, d_comm%domain%y%data%end,d_comm%ke, &
                                     field,cafptr(l))
      end do
      d_commG%dc_ptr =>d_comm

      call mpp_clock_begin(unpk_clock)
!recv
!unpack halos in reverse order
      n = d_comm%Rlist_size
      call mpp_sync(do_self=.false.)
      do list=n-1,0,-1
         if( .NOT. d_comm%R_do_buf(list) )cycle
         from_pe = d_comm%cfrom_pe(list)
         r_idx = d_commG[from_pe]%dc_ptr%Rcaf_idx(pe)  ! determine this PEs position in remote list
         do l=1,l_size  ! loop over number of fields
           ptr_field = f_addrs(l)
           do m=8,1,-1
             if( d_comm%do_thisR(1,m,list) )then  ! ne
               is=d_comm%recvis(m,list); ie=d_comm%recvie(m,list)
               js=d_comm%recvjs(m,list); je=d_comm%recvje(m,list)
               if( d_comm%do_thisR(2,m,list) )then  ! folded?
                 do k = 1,ke
                    jj = d_commG[from_pe]%dc_ptr%sendjs(m,r_idx)
                    do j = je,js,-1
                      ii = d_commG[from_pe]%dc_ptr%sendis(m,r_idx)
                      do i = ie,is,-1
                         field(i,j,k) = cafptr(l)[from_pe]%pfield(ii,jj,k)
                         ii = ii + 1
                      end do
                      jj = jj + 1
                    end do
                 end do
               else
                 do k = 1,ke
                    jj = d_commG[from_pe]%dc_ptr%sendjs(m,r_idx)
                    do j = js,je
                      ii = d_commG[from_pe]%dc_ptr%sendis(m,r_idx)
                      do i = is,ie
                         field(i,j,k) = cafptr(l)[from_pe]%pfield(ii,jj,k)
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

      call mpp_sync(do_self=.false.)
      do l=1,l_size
        cafptr(l)%pfield =>NULL()
      end do
      d_commG%dc_ptr   =>NULL()

      if(debug_gsm)then
        do l=1,l_size
          r_ptr = f_addrs(l)
          write(stdout(),*) 'Domain checksum=',mpp_chksum(r_field)
        end do
      endif
 
      call mpp_clock_end(unpk_clock)
    end subroutine MPP_DOMAINS_DO_UPDATE_3Dnew_
