    subroutine MPP_DO_GLOBAL_FIELD_3Dnew_( local, global, d_comm )
    !get a global field from a local field; local field may be on compute OR data domain
      use  mpp_datatype_mod, only: CAFPNTR_TYPE_3D_
      type(DomainCommunicator2D), intent(in), target :: d_comm
      MPP_TYPE_, intent(in),target  ::  local(:,:,:)
      MPP_TYPE_, intent(out) :: global(d_comm%domain%x%global%begin:,d_comm%domain%y%global%begin:,:)

      type d_comm_ptr_type
        type(DomainCommunicator2D), pointer :: dc_ptr
      end type
      type(d_comm_ptr_type), allocatable,save   :: d_commG[:]
      type(CAFPNTR_TYPE_3D_), allocatable, save :: cafptr[:]
      real(KIND(local(1,1,1))) :: r_field(d_comm%domain%x%global%size,d_comm%domain%y%global%size,d_comm%ke)
      pointer( r_ptr,r_field )
      integer :: i, j, k, m, n, nd
      integer :: is, ie, js, je, ke
      integer :: ii,jj,ioff,joff
      integer :: list,from_pe
      logical, save :: first_time=.true.


      ke = d_comm%ke

      if(first_time)then
        first_time = .false.
        ALLOCATE(cafptr[0:*],d_commG[0:*])
      endif

      cafptr%pfield =>local
      d_commG%dc_ptr =>d_comm

      n = d_comm%Rlist_size
      call mpp_sync(do_self=.false.)
      do list=0,n-1
         from_pe = d_comm%cfrom_pe(list)
         is=d_comm%recvis(1,list); ie=d_comm%recvie(1,list)
         js=d_comm%recvjs(1,list); je=d_comm%recvje(1,list)
         ioff = d_commG[from_pe]%dc_ptr%gf_ioff
         joff = d_commG[from_pe]%dc_ptr%gf_joff
         do k = 1,ke
           jj = d_commG[from_pe]%dc_ptr%domain%y%compute%begin
           do j = js,je
             ii = d_commG[from_pe]%dc_ptr%domain%x%compute%begin
             do i = is,ie
               global(i,j,k) = cafptr[from_pe]%pfield(ii+ioff,jj+joff,k)
               ii = ii + 1
             end do
             jj = jj + 1
           end do
         end do
      end do

      call mpp_sync(do_self=.false.)
      cafptr%pfield =>NULL()
      d_commG%dc_ptr =>NULL()

      if(debug_gsm)then
        r_ptr = LOC(global)
        write(stdout(),*) 'Global Field checksum=',mpp_chksum(r_field)
      endif
    end subroutine MPP_DO_GLOBAL_FIELD_3Dnew_
