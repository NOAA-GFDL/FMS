    subroutine MPP_DO_REDISTRIBUTE_3Dnew_( f_in, f_out, d_comm, d_type )
      use  mpp_datatype_mod, only: CAFPNTR_TYPE_3D_
      integer(LONG_KIND), intent(in)                 :: f_in(:), f_out(:)
      type(DomainCommunicator2D), intent(in), target :: d_comm
      MPP_TYPE_, intent(in)                          :: d_type

      MPP_TYPE_ :: field_in(d_comm%domain_in%x%data%begin:d_comm%domain_in%x%data%end, &
                            d_comm%domain_in%y%data%begin:d_comm%domain_in%y%data%end,d_comm%ke)
      pointer( ptr_field_in, field_in)
      MPP_TYPE_ :: field_out(d_comm%domain_out%x%data%begin:d_comm%domain_out%x%data%end, &
                             d_comm%domain_out%y%data%begin:d_comm%domain_out%y%data%end,d_comm%ke)
      pointer( ptr_field_out, field_out)

      type d_comm_ptr_type
        type(DomainCommunicator2D), pointer :: dc_ptr
      end type
      type(d_comm_ptr_type), allocatable,save   :: d_commG[:]
      type(CAFPNTR_TYPE_3D_), allocatable, save :: cafptr(:)[:]
      real(KIND(d_type)) :: r_field(d_comm%domain_out%x%data%size,d_comm%domain_out%y%data%size,d_comm%ke)
      pointer( r_ptr,r_field )
      real(KIND(d_type)),save :: d_field(1,1,1)=0.d0  ! dummy for use with mpp_chksum
      integer :: i, j, k, l, m, n, l_size
      integer :: is, ie, js, je, ke, istat
      integer :: ii,jj,ii_start,jj_start
      integer :: list,from_pe,r_idx

      logical, save :: first_time=.true.


      ke    = d_comm%ke
      l_size = size(f_in(:))  ! same as size as f_out

      if(first_time)then
        first_time=.false.
        ALLOCATE(cafptr(MAX_DOMAIN_FIELDS)[0:*],d_commG[0:*])
      endif

      do l=1,l_size
        ptr_field_in = f_in(l)
        call mpp_associate_caf_field(d_comm%domain_in%x%data%begin, d_comm%domain_in%x%data%end, &
                                     d_comm%domain_in%y%data%begin, d_comm%domain_in%y%data%end,d_comm%ke, &
                                     field_in,cafptr(l))
      end do
      d_commG%dc_ptr =>d_comm
!recv    
      n = d_comm%Rlist_size
      call mpp_sync(do_self=.false.)
      do list=0,n-1
        if( .NOT. d_comm%R_do_buf(list) )cycle
        from_pe = d_comm%cfrom_pe(list)
        r_idx = d_commG[from_pe]%dc_ptr%Rcaf_idx(pe)  ! determine this PEs position in remote list
        is=d_comm%recvis(1,list); ie=d_comm%recvie(1,list)
        js=d_comm%recvjs(1,list); je=d_comm%recvje(1,list)
        ii_start = d_commG[from_pe]%dc_ptr%sendis(1,r_idx)
        jj_start = d_commG[from_pe]%dc_ptr%sendjs(1,r_idx)
        do l=1,l_size  ! loop over number of in/out fields
           ptr_field_out = f_out(l)
           do k = 1,ke
             jj = jj_start
             do j = js,je
               ii = ii_start
               do i = is,ie
                  field_out(i,j,k) = cafptr(l)[from_pe]%pfield(ii,jj,k)
                  ii = ii + 1
               end do
               jj = jj + 1
             end do
           end do
        end do
      end do
      
      call mpp_sync(do_self=.false.)
      do l=1,l_size
        cafptr(l)%pfield =>NULL()
      end do
      d_commG%dc_ptr   =>NULL()

      if(debug_gsm)then
        do l=1,l_size  ! loop over number of in/out fields
          if(ANY(d_comm%R_do_buf(:)))then
            r_ptr = f_out(l)
            write(stdout(),*) 'Domain checksum=',mpp_chksum(r_field)
          else
            write(stdout(),*) 'Domain checksum=',mpp_chksum(d_field)
          endif
        end do
      endif
    end subroutine MPP_DO_REDISTRIBUTE_3Dnew_
