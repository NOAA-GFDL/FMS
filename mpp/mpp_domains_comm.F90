module mpp_domains_comm_mod

  use mpp_mod,           only : mpp_root_pe, mpp_malloc, mpp_error, FATAL, WARNING, NOTE, ANY_PE, NULL_PE
  use mpp_mod,           only : mpp_send, mpp_recv, mpp_sync_self, mpp_pe, mpp_npes, ALL_PES, mpp_transmit, mpp_max
  use mpp_parameter_mod,    only : MPP_DOMAIN_TIME, XUPDATE, YUPDATE, BGRID_NE, BGRID_SW, SCALAR_BIT
  use mpp_parameter_mod,    only : CGRID_NE, CGRID_SW, WEST, EAST, SOUTH, NORTH, AGRID, DOMAIN_ID_BASE


  use mpp_datatype_mod,  only : domain1D, domain2D, DomainCommunicator2D
  use mpp_data_mod,      only : mpp_domains_stack, ptr_domains_stack, pe, mpp_domains_stack_size, mpp_domains_stack_hwm
  use mpp_data_mod,      only : debug_gsm,grid_offset_type, module_is_initialized=>mpp_domains_is_initialized

  use mpp_domains_util_mod,  only : mpp_get_compute_domain, compute_overlaps

  implicit none
  private

#ifdef use_GSM
#include "mpp/shmem.fh"
#endif
#include <os.h>

  public :: mpp_update_init_comm, mpp_update_free_comm
  public :: mpp_redistribute_init_comm, mpp_redistribute_free_comm
  public :: mpp_global_field_init_comm, mpp_global_field_free_comm
#ifdef use_CAF
  public :: mpp_associate_caf_field
#endif



  ! <INTERFACE NAME="mpp_domains_init_comm">
  !  <OVERVIEW>
  !    Instantiate mpp_domains_communicator derived type
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The call returns the allocated and filled mpp_communictor derived type
  !    This data object is designed to allows all data start and end points
  !    to be calculated once and re-used in an domains update, update_V and
  !    redistribute operation
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_domains_update_init_comm(d_comm,domain,flags )
  !  </TEMPLATE>
  !  <OUT NAME="d_comm"></OUT>
  !  <IN NAME="domain"></IN>
  !  <IN NAME="flags"></IN>
  ! </INTERFACE>

#ifdef use_CAF
  interface mpp_associate_caf_field
!     module procedure associate_caf_field_r8_1d
!     module procedure associate_caf_field_c8_1d
!     module procedure associate_caf_field_i8_1d
!     module procedure associate_caf_field_l8_1d
!     module procedure associate_caf_field_r4_1d
!     module procedure associate_caf_field_c4_1d
!     module procedure associate_caf_field_i4_1d
!     module procedure associate_caf_field_l4_1d
     module procedure associate_caf_field_r8_3d
     module procedure associate_caf_field_c8_3d
     module procedure associate_caf_field_i8_3d
     module procedure associate_caf_field_l8_3d
     module procedure associate_caf_field_r4_3d
     module procedure associate_caf_field_c4_3d
     module procedure associate_caf_field_i4_3d
     module procedure associate_caf_field_l4_3d
  end interface
#endif


      integer, parameter :: MAX_ADDRS=512
      integer(LONG_KIND),dimension(MAX_ADDRS),save :: addrs_sorted=-9999  ! list of sorted local addrs
      integer,           dimension(MAX_ADDRS),save :: addrs_idx=-9999     ! idx of addr2 assoicated w/ d_comm
      integer,           dimension(MAX_ADDRS),save :: a_salvage=-9999     ! freed idx list of addr
      integer,                                save :: a_sort_len=0        ! len sorted memory list
      integer,                                save :: n_addrs=0           ! num memory addresses used

      integer, parameter :: ADDR2_BASE=X'0000000000010000'
      integer, parameter :: MAX_ADDRS2=128
      integer(LONG_KIND),dimension(MAX_ADDRS2),save :: addrs2_sorted=-9999  ! list of sorted local addrs
      integer,           dimension(MAX_ADDRS2),save :: addrs2_idx=-9999     ! idx of addr2 assoicated w/ d_comm
      integer,           dimension(MAX_ADDRS2),save :: a2_salvage=-9999     ! freed indices of addr2
      integer,                                 save :: a2_sort_len=0        ! len sorted memory list
      integer,                                 save :: n_addrs2=0           ! num memory addresses used

      integer, parameter :: MAX_DOM_IDS=128
      integer(LONG_KIND),dimension(MAX_DOM_IDS),save :: ids_sorted=-9999 ! list of sorted domain identifiers
      integer,           dimension(MAX_DOM_IDS),save :: ids_idx=-9999    ! idx of d_comm associated w/ sorted addr
      integer,                                  save :: i_sort_len=0     ! len sorted domain ids list
      integer,                                  save :: n_ids=0          ! num domain ids used (=i_sort_len; dom ids never removed)

      integer, parameter :: MAX_FIELDS=1024
      integer(LONG_KIND),        dimension(MAX_FIELDS),save           :: dcKey_sorted=-9999  ! list of sorted local addrs
!     Not sure why static d_comm fails during deallocation of derived type members; allocatable works
!     type(DomainCommunicator2D),dimension(MAX_FIELDS),save,target    :: d_comm              ! domain communicators
      type(DomainCommunicator2D),dimension(:),allocatable,save,target :: d_comm              ! domain communicators
      integer,                   dimension(MAX_FIELDS),save           :: d_comm_idx=-9999    ! idx of d_comm associated w/ sorted addr
      integer,                   dimension(MAX_FIELDS),save           :: dc_salvage=-9999    ! freed indices of d_comm
      integer,                                         save           :: dc_sort_len=0       ! len sorted comm keys (=num active communicators)
      integer,                                         save           :: n_comm=0            ! num communicators used

!     integer(LONG_KIND), parameter :: GT_BASE=2**8
      integer(LONG_KIND), parameter :: GT_BASE=X'0000000000000100'  ! Workaround for 64bit int init problem

!     integer(LONG_KIND), parameter :: KE_BASE=2**48
      integer(LONG_KIND), parameter :: KE_BASE=X'0001000000000000'  ! Workaround for 64bit int init problem


contains


    function mpp_update_init_comm(domain,l_addr,isize,jsize,ksize,l_addr2,flags,gridtype) RESULT(d_comm)
      type(DomainCommunicator2D), pointer       :: d_comm
      type(domain2D),target,      intent(inout) :: domain  ! use of compute_overlaps requires intent(inout)
      integer(LONG_KIND),         intent(in)    :: l_addr(:)
      integer,                    intent(in)    :: isize
      integer,                    intent(in)    :: jsize
      integer,                    intent(in)    :: ksize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2(:)
      integer, optional,          intent(in)    :: flags
      integer, optional,          intent(in)    :: gridtype

      integer(LONG_KIND) :: domain_id
      integer :: m, n, list, buffer_pos
      integer :: update_flags
      integer :: is, ie, js, je, ke, ioff, joff, list_size
      integer :: i, lsize,rsize,msgsize,this_pe,my_pe,to_pe,from_pe
      integer,           allocatable,dimension(:,:)   :: isL, jsL
      logical :: recv_e, recv_se, recv_s, recv_sw, recv_w, recv_nw, recv_n, recv_ne
      logical :: send_e, send_se, send_s, send_sw, send_w, send_nw, send_n, send_ne
      logical :: folded
      real,dimension(*) :: f, f2
      pointer(ptr,f);  pointer(ptr2,f2)
      real,dimension(*) :: f_r, f2_r
      pointer(ptr_r,f_r);  pointer(ptr2_r,f2_r)
      logical :: use_shmem_ptr
      character(len=8) :: text


      use_shmem_ptr=.false.   ! False for any but GSM
#ifdef use_GSM
      call check_sma_env()
      ptr=l_addr(1); if(PRESENT(l_addr2))ptr2=l_addr2(1)  ! Set up target for shmem_ptr
      use_shmem_ptr=.false.; ptr_r=shmem_ptr(f,mpp_pe()); if(ptr_r>0)use_shmem_ptr=.true.
#endif

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

      if(PRESENT(gridtype))then
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
               call mpp_error( FATAL, 'MPP_UPDATE_INIT_COMM: gridtype cannot be changed during run.' )
        end if
      !need to add code for EWS boundaries
        if( BTEST(domain%fold,WEST) .AND. BTEST(update_flags,WEST) ) &
             call mpp_error( FATAL, 'velocity stencil not yet active for WEST fold, contact author.' )
        if( BTEST(domain%fold,EAST) .AND. BTEST(update_flags,EAST) ) &
             call mpp_error( FATAL, 'velocity stencil not yet active for EAST fold, contact author.' )
        if( BTEST(domain%fold,SOUTH) .AND. BTEST(update_flags,SOUTH) ) &
             call mpp_error( FATAL, 'velocity stencil not yet active for SOUTH fold, contact author.' )
      endif

      recv_w = BTEST(update_flags,WEST)
      recv_e = BTEST(update_flags,EAST)
      recv_s = BTEST(update_flags,SOUTH)
      recv_n = BTEST(update_flags,NORTH)
      recv_ne = recv_e .AND. recv_n
      recv_se = recv_e .AND. recv_s
      recv_sw = recv_w .AND. recv_s
      recv_nw = recv_w .AND. recv_n
      send_w = recv_e
      send_e = recv_w
      send_s = recv_n
      send_n = recv_s
      send_ne = send_e .AND. send_n
      send_se = send_e .AND. send_s
      send_sw = send_w .AND. send_s
      send_nw = send_w .AND. send_n

      if( recv_w .AND. BTEST(domain%fold,WEST)  .AND. BTEST(grid_offset_type,EAST)  ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
      if( recv_s .AND. BTEST(domain%fold,SOUTH) .AND. BTEST(grid_offset_type,NORTH) ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
      if( recv_e .AND. BTEST(domain%fold,EAST)  .AND. BTEST(grid_offset_type,WEST)  ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
      if( recv_n .AND. BTEST(domain%fold,NORTH) .AND. BTEST(grid_offset_type,SOUTH) ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )

    ! Create unique domain identifier
      list_size = size(l_addr(:))
      domain_id=set_domain_id(domain%id,ksize+list_size,update_flags,gridtype)
      if(PRESENT(l_addr2))then
        d_comm =>get_comm(domain_id,l_addr(1),l_addr2(1))
      else
        d_comm =>get_comm(domain_id,l_addr(1))
      endif

      if(d_comm%initialized)return   ! Found existing field/domain communicator

      d_comm%l_addr = l_addr(1)
      d_comm%domain =>domain
      d_comm%Slist_size = size(domain%list(:))
      d_comm%isize = isize
      d_comm%jsize = jsize
      d_comm%ke = ksize; ke = ksize

!send
      lsize = d_comm%Slist_size-1
      allocate(d_comm%sendis(8,0:lsize), d_comm%sendie(8,0:lsize), &
               d_comm%sendjs(8,0:lsize), d_comm%sendje(8,0:lsize), &
               d_comm%S_msize(8,0:lsize), d_comm%do_thisS(8,0:lsize), &
               isL(8,0:lsize),jsL(8,0:lsize))
      allocate(d_comm%cto_pe(0:lsize), d_comm%S_do_buf(0:lsize))
#ifdef use_CAF
      allocate(d_comm%Rcaf_idx(0:lsize))
      d_comm%Rcaf_idx = 0
#endif
      isL=0;jsL=0
      d_comm%cto_pe=-1
      d_comm%sendis=0; d_comm%sendie=0
      d_comm%sendjs=0; d_comm%sendje=0;
      d_comm%S_msize=0
      d_comm%S_do_buf=.false.
      d_comm%do_thisS=.false.

      ioff = domain%x%data%begin
      joff = domain%y%data%begin

      do list = 0,lsize
         m = mod( domain%pos+list, lsize+1 )
         if( .NOT.domain%list(m)%overlap )cycle
         d_comm%S_do_buf(list) = .true.
         d_comm%cto_pe(list) = domain%list(m)%pe
         to_pe = d_comm%cto_pe(list)
#ifdef use_CAF
         d_comm%Rcaf_idx(to_pe) = list   ! local CAF pe needs to know it's position in remote CAF list
#endif
         if( send_w .AND. domain%list(m)%send_w%overlap )then
             d_comm%do_thisS(1,list)=.true.
             is = domain%list(m)%send_w%is; ie = domain%list(m)%send_w%ie
             js = domain%list(m)%send_w%js; je = domain%list(m)%send_w%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_w_off%is; ie = domain%list(m)%send_w_off%ie
                 js = domain%list(m)%send_w_off%js; je = domain%list(m)%send_w_off%je
             end if
             d_comm%sendis(1,list)=is; d_comm%sendie(1,list)=ie
             d_comm%sendjs(1,list)=js; d_comm%sendje(1,list)=je
             d_comm%S_msize(1,list) = (ie-is+1)*(je-js+1)*ke
             isL(1,list) = is-ioff+1; jsL(1,list) = js-joff+1
         end if
         if( send_nw .AND. domain%list(m)%send_nw%overlap )then
             d_comm%do_thisS(2,list)=.true.
             is = domain%list(m)%send_nw%is; ie = domain%list(m)%send_nw%ie
             js = domain%list(m)%send_nw%js; je = domain%list(m)%send_nw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_nw_off%is; ie = domain%list(m)%send_nw_off%ie
                 js = domain%list(m)%send_nw_off%js; je = domain%list(m)%send_nw_off%je
             end if
             d_comm%sendis(2,list)=is; d_comm%sendie(2,list)=ie
             d_comm%sendjs(2,list)=js; d_comm%sendje(2,list)=je
             d_comm%S_msize(2,list) = (ie-is+1)*(je-js+1)*ke
             isL(2,list) = is-ioff+1; jsL(2,list) = js-joff+1
         end if
         if( send_n .AND. domain%list(m)%send_n%overlap )then
             d_comm%do_thisS(3,list)=.true.
             is = domain%list(m)%send_n%is; ie = domain%list(m)%send_n%ie
             js = domain%list(m)%send_n%js; je = domain%list(m)%send_n%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_n_off%is; ie = domain%list(m)%send_n_off%ie
                 js = domain%list(m)%send_n_off%js; je = domain%list(m)%send_n_off%je
             end if
             d_comm%sendis(3,list)=is; d_comm%sendie(3,list)=ie
             d_comm%sendjs(3,list)=js; d_comm%sendje(3,list)=je
             d_comm%S_msize(3,list) = (ie-is+1)*(je-js+1)*ke
             isL(3,list) = is-ioff+1; jsL(3,list) = js-joff+1
         end if
         if( send_ne .AND. domain%list(m)%send_ne%overlap )then
             d_comm%do_thisS(4,list)=.true.
             is = domain%list(m)%send_ne%is; ie = domain%list(m)%send_ne%ie
             js = domain%list(m)%send_ne%js; je = domain%list(m)%send_ne%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_ne_off%is; ie = domain%list(m)%send_ne_off%ie
                 js = domain%list(m)%send_ne_off%js; je = domain%list(m)%send_ne_off%je
             end if
             d_comm%sendis(4,list)=is; d_comm%sendie(4,list)=ie
             d_comm%sendjs(4,list)=js; d_comm%sendje(4,list)=je
             d_comm%S_msize(4,list) = (ie-is+1)*(je-js+1)*ke
             isL(4,list) = is-ioff+1; jsL(4,list) = js-joff+1
         end if
         if( send_e .AND. domain%list(m)%send_e%overlap )then
             d_comm%do_thisS(5,list)=.true.
             is = domain%list(m)%send_e%is; ie = domain%list(m)%send_e%ie
             js = domain%list(m)%send_e%js; je = domain%list(m)%send_e%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_e_off%is; ie = domain%list(m)%send_e_off%ie
                 js = domain%list(m)%send_e_off%js; je = domain%list(m)%send_e_off%je
             end if
             d_comm%sendis(5,list)=is; d_comm%sendie(5,list)=ie
             d_comm%sendjs(5,list)=js; d_comm%sendje(5,list)=je
             d_comm%S_msize(5,list) = (ie-is+1)*(je-js+1)*ke
             isL(5,list) = is-ioff+1; jsL(5,list) = js-joff+1
         end if
         if( send_se .AND. domain%list(m)%send_se%overlap )then
             d_comm%do_thisS(6,list)=.true.
             is = domain%list(m)%send_se%is; ie = domain%list(m)%send_se%ie
             js = domain%list(m)%send_se%js; je = domain%list(m)%send_se%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_se_off%is; ie = domain%list(m)%send_se_off%ie
                 js = domain%list(m)%send_se_off%js; je = domain%list(m)%send_se_off%je
             end if
             d_comm%sendis(6,list)=is; d_comm%sendie(6,list)=ie
             d_comm%sendjs(6,list)=js; d_comm%sendje(6,list)=je
             d_comm%S_msize(6,list) = (ie-is+1)*(je-js+1)*ke
             isL(6,list) = is-ioff+1; jsL(6,list) = js-joff+1
         end if
         if( send_s .AND. domain%list(m)%send_s%overlap )then
             d_comm%do_thisS(7,list)=.true.
             is = domain%list(m)%send_s%is; ie = domain%list(m)%send_s%ie
             js = domain%list(m)%send_s%js; je = domain%list(m)%send_s%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_s_off%is; ie = domain%list(m)%send_s_off%ie
                 js = domain%list(m)%send_s_off%js; je = domain%list(m)%send_s_off%je
             end if
             d_comm%sendis(7,list)=is; d_comm%sendie(7,list)=ie
             d_comm%sendjs(7,list)=js; d_comm%sendje(7,list)=je
             d_comm%S_msize(7,list) = (ie-is+1)*(je-js+1)*ke
             isL(7,list) = is-ioff+1; jsL(7,list) = js-joff+1
         end if
         if( send_sw .AND. domain%list(m)%send_sw%overlap )then
             d_comm%do_thisS(8,list)=.true.
             is = domain%list(m)%send_sw%is; ie = domain%list(m)%send_sw%ie
             js = domain%list(m)%send_sw%js; je = domain%list(m)%send_sw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_sw_off%is; ie = domain%list(m)%send_sw_off%ie
                 js = domain%list(m)%send_sw_off%js; je = domain%list(m)%send_sw_off%je
             end if
             d_comm%sendis(8,list)=is; d_comm%sendie(8,list)=ie
             d_comm%sendjs(8,list)=js; d_comm%sendje(8,list)=je
             d_comm%S_msize(8,list) = (ie-is+1)*(je-js+1)*ke
             isL(8,list) = is-ioff+1; jsL(8,list) = js-joff+1
         end if
         if(.not.ANY(d_comm%do_thisS(:,list)))d_comm%S_do_buf(list) = .false.
#ifdef use_GSM
         if(d_comm%S_do_buf(list))then
           if(.not.use_shmem_ptr)then  ! address is from GLOBAL_ALLOC
             if(.not.PRESENT(l_addr2))then
               call mpp_send(l_addr(1), plen=list_size, to_pe=to_pe )
             else
               call mpp_send(l_addr(1), plen=list_size, to_pe=to_pe )
               call mpp_send(l_addr2(1), plen=list_size, to_pe=to_pe )
             endif
           endif
           call mpp_send(isL(1,list), plen=8, to_pe=to_pe )
           call mpp_send(jsL(1,list), plen=8, to_pe=to_pe )
           call mpp_send(isize, plen=1, to_pe=to_pe )
           call mpp_send(jsize, plen=1, to_pe=to_pe )
         endif
#endif
      end do

      call mpp_sync_self()
!recv
      d_comm%Rlist_size = size(domain%list(:))

      rsize = d_comm%Rlist_size-1
      allocate(d_comm%recvis(8,0:rsize), d_comm%recvie(8,0:rsize), &
               d_comm%recvjs(8,0:rsize), d_comm%recvje(8,0:rsize), &
               d_comm%R_msize(8,0:rsize), d_comm%do_thisR(2,8,0:rsize))
      allocate(d_comm%cfrom_pe(0:rsize), d_comm%R_do_buf(0:rsize))
      allocate(d_comm%isizeR(0:rsize), d_comm%jsizeR(0:rsize))
      allocate(d_comm%sendisR(8,0:rsize), d_comm%sendjsR(8,0:rsize))
#ifdef use_GSM
      if(.not.PRESENT(l_addr2))then                                 
        allocate(d_comm%rem_addrl(list_size,0:rsize))
        d_comm%rem_addrl=-9999
      else
        allocate(d_comm%rem_addrlx(list_size,0:rsize),d_comm%rem_addrly(list_size,0:rsize))
        d_comm%rem_addrlx=-9999; d_comm%rem_addrly=-9999
      endif
#endif

      d_comm%cfrom_pe=-1
      d_comm%recvis=0; d_comm%recvie=0
      d_comm%recvjs=0; d_comm%recvje=0;
      d_comm%R_msize=0
      d_comm%R_do_buf=.false.
      d_comm%do_thisR=.false.
      d_comm%isizeR=0; d_comm%jsizeR=0
      d_comm%sendisR=0; d_comm%sendjsR=0
      do list = 0,rsize
         m = mod( domain%pos+rsize+1-list, rsize+1 )
         if( .NOT.domain%list(m)%overlap )cycle
         d_comm%cfrom_pe(list) = domain%list(m)%pe
         from_pe = d_comm%cfrom_pe(list)
         d_comm%R_do_buf(list) = .true.
         if( recv_e .AND. domain%list(m)%recv_e%overlap )then
             d_comm%do_thisR(1,1,list)=.true.
             is = domain%list(m)%recv_e%is; ie = domain%list(m)%recv_e%ie
             js = domain%list(m)%recv_e%js; je = domain%list(m)%recv_e%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_e_off%is; ie = domain%list(m)%recv_e_off%ie
                 js = domain%list(m)%recv_e_off%js; je = domain%list(m)%recv_e_off%je
             end if
             d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
             d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
             d_comm%R_msize(1,list) = (ie-is+1)*(je-js+1)*ke
             if(domain%list(m)%recv_e%folded)d_comm%do_thisR(2,1,list)=.true.
         end if
         if( recv_se .AND. domain%list(m)%recv_se%overlap )then
             d_comm%do_thisR(1,2,list)=.true.
             is = domain%list(m)%recv_se%is; ie = domain%list(m)%recv_se%ie
             js = domain%list(m)%recv_se%js; je = domain%list(m)%recv_se%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_se_off%is; ie = domain%list(m)%recv_se_off%ie
                 js = domain%list(m)%recv_se_off%js; je = domain%list(m)%recv_se_off%je
             end if
             d_comm%recvis(2,list)=is; d_comm%recvie(2,list)=ie
             d_comm%recvjs(2,list)=js; d_comm%recvje(2,list)=je
             d_comm%R_msize(2,list) = (ie-is+1)*(je-js+1)*ke
             if(domain%list(m)%recv_se%folded)d_comm%do_thisR(2,2,list)=.true.
         end if
         if( recv_s .AND. domain%list(m)%recv_s%overlap )then
             d_comm%do_thisR(1,3,list)=.true.
             is = domain%list(m)%recv_s%is; ie = domain%list(m)%recv_s%ie
             js = domain%list(m)%recv_s%js; je = domain%list(m)%recv_s%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_s_off%is; ie = domain%list(m)%recv_s_off%ie
                 js = domain%list(m)%recv_s_off%js; je = domain%list(m)%recv_s_off%je
             end if
             d_comm%recvis(3,list)=is; d_comm%recvie(3,list)=ie
             d_comm%recvjs(3,list)=js; d_comm%recvje(3,list)=je
             d_comm%R_msize(3,list) = (ie-is+1)*(je-js+1)*ke
             if(domain%list(m)%recv_s%folded)d_comm%do_thisR(2,3,list)=.true.
         end if
         if( recv_sw .AND. domain%list(m)%recv_sw%overlap )then
             d_comm%do_thisR(1,4,list)=.true.
             is = domain%list(m)%recv_sw%is; ie = domain%list(m)%recv_sw%ie
             js = domain%list(m)%recv_sw%js; je = domain%list(m)%recv_sw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_sw_off%is; ie = domain%list(m)%recv_sw_off%ie
                 js = domain%list(m)%recv_sw_off%js; je = domain%list(m)%recv_sw_off%je
             end if
             d_comm%recvis(4,list)=is; d_comm%recvie(4,list)=ie
             d_comm%recvjs(4,list)=js; d_comm%recvje(4,list)=je
             d_comm%R_msize(4,list) = (ie-is+1)*(je-js+1)*ke
             if(domain%list(m)%recv_sw%folded)d_comm%do_thisR(2,4,list)=.true.
         end if
         if( recv_w .AND. domain%list(m)%recv_w%overlap )then
             d_comm%do_thisR(1,5,list)=.true.
             is = domain%list(m)%recv_w%is; ie = domain%list(m)%recv_w%ie
             js = domain%list(m)%recv_w%js; je = domain%list(m)%recv_w%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_w_off%is; ie = domain%list(m)%recv_w_off%ie
                 js = domain%list(m)%recv_w_off%js; je = domain%list(m)%recv_w_off%je
             end if
             d_comm%recvis(5,list)=is; d_comm%recvie(5,list)=ie
             d_comm%recvjs(5,list)=js; d_comm%recvje(5,list)=je
             d_comm%R_msize(5,list) = (ie-is+1)*(je-js+1)*ke
             if(domain%list(m)%recv_w%folded)d_comm%do_thisR(2,5,list)=.true.
         end if
         if( recv_nw .AND. domain%list(m)%recv_nw%overlap )then
             d_comm%do_thisR(1,6,list)=.true.
             is = domain%list(m)%recv_nw%is; ie = domain%list(m)%recv_nw%ie
             js = domain%list(m)%recv_nw%js; je = domain%list(m)%recv_nw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_nw_off%is; ie = domain%list(m)%recv_nw_off%ie
                 js = domain%list(m)%recv_nw_off%js; je = domain%list(m)%recv_nw_off%je
             end if
             d_comm%recvis(6,list)=is; d_comm%recvie(6,list)=ie
             d_comm%recvjs(6,list)=js; d_comm%recvje(6,list)=je
             d_comm%R_msize(6,list) = (ie-is+1)*(je-js+1)*ke
             if(domain%list(m)%recv_nw%folded)d_comm%do_thisR(2,6,list)=.true.
         end if
         if( recv_n .AND. domain%list(m)%recv_n%overlap )then
             d_comm%do_thisR(1,7,list)=.true.
             is = domain%list(m)%recv_n%is; ie = domain%list(m)%recv_n%ie
             js = domain%list(m)%recv_n%js; je = domain%list(m)%recv_n%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_n_off%is; ie = domain%list(m)%recv_n_off%ie
                 js = domain%list(m)%recv_n_off%js; je = domain%list(m)%recv_n_off%je
             end if
             d_comm%recvis(7,list)=is; d_comm%recvie(7,list)=ie
             d_comm%recvjs(7,list)=js; d_comm%recvje(7,list)=je
             d_comm%R_msize(7,list) = (ie-is+1)*(je-js+1)*ke
             if(domain%list(m)%recv_n%folded)d_comm%do_thisR(2,7,list)=.true.
         end if
         if( recv_ne .AND. domain%list(m)%recv_ne%overlap )then
             d_comm%do_thisR(1,8,list)=.true.
             is = domain%list(m)%recv_ne%is; ie = domain%list(m)%recv_ne%ie
             js = domain%list(m)%recv_ne%js; je = domain%list(m)%recv_ne%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_ne_off%is; ie = domain%list(m)%recv_ne_off%ie
                 js = domain%list(m)%recv_ne_off%js; je = domain%list(m)%recv_ne_off%je
             end if
             d_comm%recvis(8,list)=is; d_comm%recvie(8,list)=ie
             d_comm%recvjs(8,list)=js; d_comm%recvje(8,list)=je
             d_comm%R_msize(8,list) = (ie-is+1)*(je-js+1)*ke
             if(domain%list(m)%recv_ne%folded)d_comm%do_thisR(2,8,list)=.true.
         end if
         if(.not.ANY(d_comm%do_thisR(1,:,list)))d_comm%R_do_buf(list) = .false.
#ifdef use_GSM
         if(d_comm%R_do_buf(list))then
           if(.not.PRESENT(l_addr2))then
             if(use_shmem_ptr)then                          
               do i=1,list_size   
                 ptr = l_addr(i)
                 d_comm%rem_addrl(i,list) = shmem_ptr(f,from_pe)
               end do
             else
               call mpp_recv(d_comm%rem_addrl(1,list), glen=list_size, from_pe=from_pe )
             endif
           else
             if(use_shmem_ptr)then
               do i=1,list_size
                 ptr  = l_addr(i)
                 ptr2 = l_addr2(i)
                 d_comm%rem_addrlx(i,list) = shmem_ptr(f,from_pe)
                 d_comm%rem_addrly(i,list) = shmem_ptr(f2,from_pe)
               end do
             else    
               call mpp_recv(d_comm%rem_addrlx(1,list), glen=list_size, from_pe=from_pe )
               call mpp_recv(d_comm%rem_addrly(1,list), glen=list_size, from_pe=from_pe )
             endif
           endif
           call mpp_recv(d_comm%sendisR(1,list), glen=8, from_pe=from_pe )
           call mpp_recv(d_comm%sendjsR(1,list), glen=8, from_pe=from_pe )
           call mpp_recv(d_comm%isizeR(list), glen=1, from_pe=from_pe )
           call mpp_recv(d_comm%jsizeR(list), glen=1, from_pe=from_pe )
         endif
#endif
      end do

      d_comm%isize_max = isize; call mpp_max(d_comm%isize_max)
      d_comm%jsize_max = jsize; call mpp_max(d_comm%jsize_max)

#if defined(use_libMPI) || defined(use_libSMA)
    ! Handles case where S_msize and/or R_msize are 0 size array
      msgsize = ( MAXVAL( (/0,sum(d_comm%S_msize(:,:))/) ) + MAXVAL( (/0,sum(d_comm%R_msize(:,:))/) ) ) * list_size
      if(msgsize>0)then
         mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, msgsize )
         if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
             write( text,'(i8)' )mpp_domains_stack_hwm
             call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                  //trim(text)//') from all PEs.' )
         end if
      end if
#endif

      DEALLOCATE(isL,jsL)

      d_comm%initialized = .true.

    end function mpp_update_init_comm


    function mpp_redistribute_init_comm(domain_in,l_addrs_in, domain_out,l_addrs_out, &
                                        isize_in,jsize_in,ksize_in,isize_out,jsize_out,ksize_out) RESULT(d_comm)
      type(DomainCommunicator2D), pointer       :: d_comm
      type(domain2D),target,      intent(in)    :: domain_in
      integer(LONG_KIND),         intent(in)    :: l_addrs_in(:)
      type(domain2D),target,      intent(in)    :: domain_out
      integer(LONG_KIND),         intent(in)    :: l_addrs_out(:)
      integer,                    intent(in)    :: isize_in
      integer,                    intent(in)    :: jsize_in
      integer,                    intent(in)    :: ksize_in
      integer,                    intent(in)    :: isize_out 
      integer,                    intent(in)    :: jsize_out
      integer,                    intent(in)    :: ksize_out 

      integer(LONG_KIND) :: domain_id
      integer :: m, n, list
      integer :: is, ie, js, je, ke, ioff, joff, list_size
      integer :: isc, iec, jsc, jec
      integer :: i, lsize,rsize,msgsize,this_pe,my_pe,to_pe,from_pe
      integer,           allocatable,dimension(:) :: isL, jsL
      integer(LONG_KIND),allocatable,dimension(:,:) :: slist_addr
      real,dimension(*) :: f
      pointer(ptr,f)
      real,dimension(*) :: f_r
      pointer(ptr_r,f_r)
      logical :: use_shmem_ptr
      character(len=8) :: text


    ! This test determines whether input fields are from allocated memory (LOC gets global
    ! address) or "static" memory (need shmem_ptr). This probably needs to be generalized
    ! to determine appropriate mechanism for each incoming address.

    ! "Concurrent" run mode may leave field_in or field_out unassociated if pe does not
    ! contain in/out data. Use of STATIC option for ocean complicates this as ocean component
    ! always defined. Field_out is always a boundary structure and so is always allocated or
    ! not depending on whether it's used. If field out is defined (>0), then it is used otherwise
    ! field in must be defined.

      use_shmem_ptr=.false.   ! False for any but use_GSM
#ifdef use_GSM
      call check_sma_env()
      if(l_addrs_in(1)>0)then  ! Set up target for shmem_ptr
        ptr=l_addrs_in(1)
        ptr_r=shmem_ptr(f,mpp_pe())
        if(ptr_r>0)use_shmem_ptr=.true.
      endif
#endif

!fix ke
      ke = 0
      if( domain_in%pe /= NULL_PE )ke = ksize_in
      if( domain_out%pe /= NULL_PE )then
          if( ke /= 0 .AND. ke /= ksize_out ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: mismatch between field_in and field_out.' )
          ke = ksize_out
      end if
      if( ke == 0 )call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: either domain_in or domain_out must be native.' )
!check sizes
      if( domain_in%pe /= NULL_PE )then
          if( isize_in /= domain_in%x%data%size .OR. jsize_in /= domain_in%y%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: field_in must be on data domain of domain_in.' )
      end if
      if( domain_out%pe /= NULL_PE )then
          if( isize_out /= domain_out%x%data%size .OR. jsize_out /= domain_out%y%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: field_out must be on data domain of domain_out.' )
      end if


    ! Create unique domain identifier
      list_size = size(l_addrs_in(:))
      if(l_addrs_out(1) > 0)then
        domain_id = set_domain_id(domain_out%id,ke+list_size)
      else
        domain_id = set_domain_id(domain_in%id,ke+list_size)
      endif

      d_comm =>get_comm(domain_id,l_addrs_in(1),l_addrs_out(1))

      if(d_comm%initialized)return   ! Found existing field/domain communicator

      d_comm%l_addr = l_addrs_in(1)
      d_comm%domain_in =>domain_in
      d_comm%Slist_size = size(domain_out%list(:))
      d_comm%isize_in = isize_in
      d_comm%jsize_in = jsize_in
      d_comm%ke = ke

!send 
      lsize = d_comm%Slist_size-1
      allocate(d_comm%sendis(1,0:lsize), d_comm%sendie(1,0:lsize), &
               d_comm%sendjs(1,0:lsize), d_comm%sendje(1,0:lsize), &
               d_comm%S_msize(1,0:lsize),isL(0:lsize),jsL(0:lsize))
      allocate(slist_addr(list_size,0:lsize))
      allocate(d_comm%cto_pe(0:lsize), d_comm%S_do_buf(0:lsize))
#ifdef use_CAF
      allocate(d_comm%Rcaf_idx(0:lsize))
      d_comm%Rcaf_idx = 0
#endif
      isL=0;jsL=0
      slist_addr = -9999
      d_comm%cto_pe=-1
      d_comm%sendis=0; d_comm%sendie=0
      d_comm%sendjs=0; d_comm%sendje=0;
      d_comm%S_msize=0
      d_comm%S_do_buf=.false.

      ioff = domain_in%x%data%begin
      joff = domain_in%y%data%begin

      call mpp_get_compute_domain( domain_in, isc, iec, jsc, jec )
      do list = 0,lsize
         m = mod( domain_out%pos+list+lsize+1, lsize+1 )
         d_comm%cto_pe(list) = domain_out%list(m)%pe
         to_pe =  d_comm%cto_pe(list)
#ifdef use_CAF
         d_comm%Rcaf_idx(to_pe) = list   ! local CAF pe needs to know it's position in remote CAF list
#endif
         call mpp_get_compute_domain( domain_out%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie >= is .AND. je >= js )then
           d_comm%S_do_buf(list) = .true.
           d_comm%sendis(1,list)=is; d_comm%sendie(1,list)=ie
           d_comm%sendjs(1,list)=js; d_comm%sendje(1,list)=je
           d_comm%S_msize(1,list) = (ie-is+1)*(je-js+1)*ke
           isL(list) = is-ioff+1; jsL(list) = js-joff+1
#ifdef use_GSM
           if(.not.use_shmem_ptr)then  ! address is from GLOBAL_ALLOC
             slist_addr(1:list_size,list) = l_addrs_in(1:list_size)
             call mpp_send(slist_addr(1,list), plen=list_size, to_pe=to_pe )
           endif
           call mpp_send(isL(list), plen=1, to_pe=to_pe )
           call mpp_send(jsL(list), plen=1, to_pe=to_pe )
           call mpp_send(isize_in, plen=1, to_pe=to_pe )   ! I extent of remote array
           call mpp_send(jsize_in, plen=1, to_pe=to_pe )   ! J extent of remote array
#endif
         end if
      end do 

      call mpp_sync_self()
!recv 
      d_comm%domain_out =>domain_out
      d_comm%Rlist_size = size(domain_in%list(:))
      d_comm%isize_out = isize_out
      d_comm%jsize_out = jsize_out

      rsize = d_comm%Rlist_size-1
      allocate(d_comm%recvis(1,0:rsize), d_comm%recvie(1,0:rsize), &
               d_comm%recvjs(1,0:rsize), d_comm%recvje(1,0:rsize), &
               d_comm%R_msize(1,0:rsize))
      allocate(d_comm%cfrom_pe(0:rsize), d_comm%R_do_buf(0:rsize))
      allocate(d_comm%isizeR(0:rsize), d_comm%jsizeR(0:rsize))
      allocate(d_comm%sendisR(1,0:rsize), d_comm%sendjsR(1,0:rsize))
      allocate(d_comm%rem_addrl(list_size,0:rsize))
      d_comm%rem_addrl=-9999
      d_comm%cfrom_pe=-1
      d_comm%recvis=0; d_comm%recvie=0
      d_comm%recvjs=0; d_comm%recvje=0;
      d_comm%R_msize=0
      d_comm%R_do_buf=.false.
      d_comm%isizeR=0; d_comm%jsizeR=0
      d_comm%sendisR=0; d_comm%sendjsR=0

      call mpp_get_compute_domain( domain_out, isc, iec, jsc, jec )
      do list = 0,rsize
         m = mod( domain_in%pos+rsize+1-list, rsize+1 )
         d_comm%cfrom_pe(list) = domain_in%list(m)%pe
         from_pe = d_comm%cfrom_pe(list)
         call mpp_get_compute_domain( domain_in%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie >= is .AND. je >= js )then
           d_comm%R_do_buf(list) = .true.
           d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
           d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
           d_comm%R_msize(1,list) = (ie-is+1)*(je-js+1)*ke
#ifdef use_GSM
           if(use_shmem_ptr)then
             do i=1,list_size
               ptr = l_addrs_in(i)
               d_comm%rem_addrl(i,list) = shmem_ptr(f,from_pe)
             end do
           else
             call mpp_recv(d_comm%rem_addrl(1,list), glen=list_size, from_pe=from_pe )
           endif
           call mpp_recv(d_comm%sendisR(1,list), glen=1, from_pe=from_pe )
           call mpp_recv(d_comm%sendjsR(1,list), glen=1, from_pe=from_pe )
           call mpp_recv(d_comm%isizeR(list), glen=1, from_pe=from_pe )
           call mpp_recv(d_comm%jsizeR(list), glen=1, from_pe=from_pe )
#endif
         end if
      end do 

#ifdef use_CAF
    ! just need address of local input array for CAF
      d_comm%rem_addrl(:,1) = l_addrs_in(1:list_size)
#endif

      d_comm%isize_max = isize_in; call mpp_max(d_comm%isize_max)
      d_comm%jsize_max = jsize_in; call mpp_max(d_comm%jsize_max)

#if defined(use_libMPI) || defined(use_libSMA)
    ! Handles case where S_msize and/or R_msize are 0 size array
      msgsize = ( MAXVAL( (/0,sum(d_comm%S_msize(:,:))/) ) + MAXVAL( (/0,sum(d_comm%R_msize(:,:))/) ) ) * list_size
      if(msgsize>0)then
         mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, msgsize )
         if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
             write( text,'(i8)' )mpp_domains_stack_hwm
             call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                  //trim(text)//') from all PEs.' )
         end if
      end if
#endif

      DEALLOCATE(slist_addr,isL,jsL)

      d_comm%initialized = .true.

    end function mpp_redistribute_init_comm


    function mpp_global_field_init_comm(domain,l_addr,isize_g,jsize_g,isize_l,jsize_l,ksize,l_addr2,flags) RESULT(d_comm)
      type(DomainCommunicator2D), pointer       :: d_comm
      type(domain2D),target,      intent(in)    :: domain
      integer(LONG_KIND),         intent(in)    :: l_addr
      integer,                    intent(in)    :: isize_g
      integer,                    intent(in)    :: jsize_g
      integer,                    intent(in)    :: isize_l
      integer,                    intent(in)    :: jsize_l 
      integer,                    intent(in)    :: ksize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2
      integer, optional,          intent(in)    :: flags

      integer(LONG_KIND) :: domain_id
      integer :: m, n, lpos, rpos, list
      integer :: update_flags
      logical :: xonly, yonly
      integer :: is, ie, js, je, ke, ioff, joff
      integer :: i, lsize,rsize,msgsize,this_pe,my_pe,to_pe,from_pe
      integer,           allocatable,dimension(:) :: isL, jsL
      integer(LONG_KIND),allocatable,dimension(:,:) :: slist_addr
      integer(LONG_KIND),save       ,dimension(2)   :: rem_addr
      real,dimension(*) :: f, f2
      pointer(ptr,f);  pointer(ptr2,f2)
      real,dimension(*) :: f_r, f2_r
      pointer(ptr_r,f_r);  pointer(ptr2_r,f2_r)
      logical :: use_shmem_ptr
      character(len=8) :: text


      use_shmem_ptr=.false.  ! False for any but use_GSM
#ifdef use_GSM
      call check_sma_env()
      ptr=l_addr; if(PRESENT(l_addr2))ptr2=l_addr2  ! Set up target for shmem_ptr
      use_shmem_ptr=.false.  ! False for CAF ; ptr_r=shmem_ptr(f,mpp_pe()); if(ptr_r>0)use_shmem_ptr=.true.
#endif

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: must first call mpp_domains_init.' )
      update_flags=0; xonly = .FALSE.;  yonly = .FALSE.
      if( PRESENT(flags) )then
          update_flags = flags
          xonly = flags == XUPDATE
          yonly = flags == YUPDATE
          if( .NOT.xonly .AND. .NOT.yonly )call mpp_error( WARNING, 'MPP_GLOBAL_FIELD: you must have flags=XUPDATE or YUPDATE.' )
      end if

      if( isize_g /= domain%x%global%size .OR. jsize_g /= domain%y%global%size ) &
           call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_INIT_COMM: incoming arrays do not match domain.' )

      if( isize_l == domain%x%compute%size .AND. jsize_l == domain%y%compute%size )then
!local is on compute domain                                                                      
        ioff = -domain%x%compute%begin + 1                                                       
        joff = -domain%y%compute%begin + 1
      elseif( isize_l == domain%x%data%size .AND. jsize_l == domain%y%data%size )then
!local is on data domain                                                                        
        ioff = -domain%x%data%begin + 1                                                         
        joff = -domain%y%data%begin + 1
      else
        call mpp_error(FATAL,'MPP_GLOBAL_FIELD_INIT_COMM: incoming field array must match either compute domain or data domain.')
      endif

    ! Create unique domain identifier
      domain_id=set_domain_id(domain%id,ksize,update_flags)
      d_comm =>get_comm(domain_id,l_addr,l_addr2)
      
      if(d_comm%initialized)return   ! Found existing field/domain communicator

      d_comm%domain =>domain
      d_comm%isize_in = isize_l; d_comm%isize_out = isize_g
      d_comm%jsize_in = jsize_l; d_comm%jsize_out = jsize_g
      d_comm%ke = ksize
      d_comm%gf_ioff=ioff; d_comm%gf_joff=joff

!fill off-domains (note loops begin at an offset of 1)
      if( xonly )then
!send
        d_comm%Slist_size = size(domain%x%list(:))
        lsize = d_comm%Slist_size-1
        allocate(d_comm%cto_pe(0:lsize))
        d_comm%cto_pe=-1
        do list = 0,lsize
          lpos = mod(domain%x%pos+lsize+1-list,lsize+1)
          d_comm%cto_pe(list) = domain%x%list(lpos)%pe
        end do
!recv
        d_comm%Rlist_size = d_comm%Slist_size
        rsize = d_comm%Rlist_size-1
        allocate(d_comm%cfrom_pe(0:rsize))
        allocate(d_comm%recvis(1,0:rsize), d_comm%recvie(1,0:rsize), &
                 d_comm%recvjs(1,0:rsize), d_comm%recvje(1,0:rsize), &
                 d_comm%R_msize(1,0:rsize))                          
        d_comm%cfrom_pe=-1
        d_comm%recvis=0; d_comm%recvie=0
        d_comm%recvjs=0; d_comm%recvje=0;
        d_comm%R_msize=0
        do list = 0,rsize
          rpos = mod(domain%x%pos+list,rsize+1)
          d_comm%cfrom_pe(list) = domain%x%list(rpos)%pe
          is = domain%x%list(rpos)%compute%begin; ie = domain%x%list(rpos)%compute%end
          js = domain%y%compute%begin; je = domain%y%compute%end
          d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
          d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
          d_comm%R_msize(1,list) = domain%x%list(rpos)%compute%size * domain%x%list(rpos)%compute%size * ksize
        end do

      elseif( yonly )then

!send  
        d_comm%Slist_size = size(domain%y%list(:))
        lsize = d_comm%Slist_size-1
        allocate(d_comm%cto_pe(0:lsize))
        d_comm%cto_pe=-1
        do list = 0,lsize
          lpos = mod(domain%y%pos+lsize+1-list,lsize+1)
          d_comm%cto_pe(list) = domain%y%list(lpos)%pe
        end do
!recv    
        d_comm%Rlist_size = d_comm%Slist_size
        rsize = d_comm%Rlist_size-1
        allocate(d_comm%cfrom_pe(0:rsize))
        allocate(d_comm%recvis(1,0:rsize), d_comm%recvie(1,0:rsize), &
                 d_comm%recvjs(1,0:rsize), d_comm%recvje(1,0:rsize), &
                 d_comm%R_msize(1,0:rsize))                          
        d_comm%cfrom_pe=-1
        d_comm%recvis=0; d_comm%recvie=0
        d_comm%recvjs=0; d_comm%recvje=0;
        d_comm%R_msize=0
        do list = 0,rsize
          rpos = mod(domain%y%pos+list,rsize+1)
          d_comm%cfrom_pe(list) = domain%y%list(rpos)%pe
          is = domain%x%compute%begin; ie = domain%x%compute%end
          js = domain%y%list(rpos)%compute%begin; je = domain%y%list(rpos)%compute%end
          d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
          d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
          d_comm%R_msize(1,list) = domain%y%list(rpos)%compute%size * domain%y%list(rpos)%compute%size * ksize
        end do

      else

!send
        d_comm%Slist_size =  size(domain%list(:))
        lsize = d_comm%Slist_size-1
        allocate(d_comm%cto_pe(0:lsize))
        d_comm%cto_pe=-1
        do list = 0,lsize
          lpos = mod(domain%pos+lsize+1-list,lsize+1)
          d_comm%cto_pe(list) = domain%list(lpos)%pe
        end do
!recv
        d_comm%Rlist_size = d_comm%Slist_size
        rsize = d_comm%Rlist_size-1
        allocate(d_comm%cfrom_pe(0:rsize))
        allocate(d_comm%recvis(1,0:rsize), d_comm%recvie(1,0:rsize), &
                 d_comm%recvjs(1,0:rsize), d_comm%recvje(1,0:rsize), &
                 d_comm%R_msize(1,0:rsize))
        d_comm%cfrom_pe=-1
        d_comm%recvis=0; d_comm%recvie=0
        d_comm%recvjs=0; d_comm%recvje=0;
        d_comm%R_msize=0
        do list = 0,rsize
          rpos = mod(domain%pos+list,rsize+1)
          d_comm%cfrom_pe(list) = domain%list(rpos)%pe
          is = domain%list(rpos)%x%compute%begin; ie = domain%list(rpos)%x%compute%end
          js = domain%list(rpos)%y%compute%begin; je = domain%list(rpos)%y%compute%end
          d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
          d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
          d_comm%R_msize(1,list) = domain%list(rpos)%x%compute%size * domain%list(rpos)%y%compute%size * ksize
        end do

      endif


!send 

      allocate(d_comm%sendis(1,0:lsize), d_comm%sendie(1,0:lsize), &
               d_comm%sendjs(1,0:lsize), d_comm%sendje(1,0:lsize), &
               d_comm%S_msize(1,0:lsize),isL(0:lsize),jsL(0:lsize))
      allocate(slist_addr(2,0:lsize))
      isL=0; jsL=0
      slist_addr = -9999
      d_comm%sendis=0; d_comm%sendie=0
      d_comm%sendjs=0; d_comm%sendje=0;
      d_comm%S_msize=0
      do list = 0,lsize
        to_pe = d_comm%cto_pe(list)
        is=domain%x%compute%begin; ie=domain%x%compute%end
        js=domain%y%compute%begin; je=domain%y%compute%end
        d_comm%sendis(1,list)=is; d_comm%sendie(1,list)=ie
        d_comm%sendjs(1,list)=js; d_comm%sendje(1,list)=je
        d_comm%S_msize(1,list) = domain%x%compute%size * domain%y%compute%size * ksize
        isL(list) = ioff+domain%x%compute%begin; jsL(list) = joff+domain%y%compute%begin
#ifdef use_GSM
        if(.not.use_shmem_ptr)then  ! address is from GLOBAL_ALLOC
          if(.not.PRESENT(l_addr2))then
            slist_addr(1,list) = d_comm%l_addr
            call mpp_send(slist_addr(1,list), plen=1, to_pe=to_pe )
          else
            slist_addr(1,list) = d_comm%l_addrx
            slist_addr(2,list) = d_comm%l_addry
            call mpp_send(slist_addr(1,list), plen=2, to_pe=to_pe )
          endif
        endif
        call mpp_send(isL(list), plen=1, to_pe=to_pe )
        call mpp_send(jsL(list), plen=1, to_pe=to_pe )
        call mpp_send(isize_l, plen=1, to_pe=to_pe )
        call mpp_send(jsize_l, plen=1, to_pe=to_pe )
#endif
      end do

      if( xonly ) then
        call mpp_sync_self(domain%x%list(:)%pe)
      elseif( yonly ) then
        call mpp_sync_self(domain%y%list(:)%pe)
      else
        call mpp_sync_self(domain%list(:)%pe)
      endif

!recv
      allocate(d_comm%isizeR(0:rsize), d_comm%jsizeR(0:rsize))
      allocate(d_comm%sendisR(1,0:rsize), d_comm%sendjsR(1,0:rsize))
      if(.not.PRESENT(l_addr2))then
        allocate(d_comm%rem_addr(0:rsize))
        d_comm%rem_addr=-9999
      else
        allocate(d_comm%rem_addrx(0:rsize),d_comm%rem_addry(0:rsize))
        d_comm%rem_addrx=-9999; d_comm%rem_addry=-9999
      endif
      d_comm%isizeR=0; d_comm%jsizeR=0
      d_comm%sendisR=0; d_comm%sendjsR=0
      rem_addr = -9999
      do list = 0,rsize
        from_pe = d_comm%cfrom_pe(list)
#ifdef use_GSM
        if(.not.PRESENT(l_addr2))then
          if(use_shmem_ptr)then
            d_comm%rem_addr(list) = shmem_ptr(f,from_pe)
          else
            call mpp_recv(rem_addr(1), glen=1, from_pe=from_pe )
            d_comm%rem_addr(list) = rem_addr(1)
          endif
        else   
          if(use_shmem_ptr)then
            d_comm%rem_addrx(list) = shmem_ptr(f,from_pe)
            d_comm%rem_addry(list) = shmem_ptr(f2,from_pe)
          else
            call mpp_recv(rem_addr(1), glen=2, from_pe=from_pe )
            d_comm%rem_addrx(list) = rem_addr(1)
            d_comm%rem_addry(list) = rem_addr(2)
          endif
        endif  
        call mpp_recv(d_comm%sendisR(1,list), glen=1, from_pe=from_pe )
        call mpp_recv(d_comm%sendjsR(1,list), glen=1, from_pe=from_pe )
        call mpp_recv(d_comm%isizeR(list), glen=1, from_pe=from_pe )
        call mpp_recv(d_comm%jsizeR(list), glen=1, from_pe=from_pe )
#endif
      end do

#ifdef use_CAF
      if(.not.PRESENT(l_addr2))then
        d_comm%rem_addr(0) = d_comm%l_addr
      else
        d_comm%rem_addrx(0) = d_comm%l_addrx
        d_comm%rem_addry(0) = d_comm%l_addry
      endif
#endif

      d_comm%isize_max = isize_l; call mpp_max(d_comm%isize_max)
      d_comm%jsize_max = jsize_l; call mpp_max(d_comm%jsize_max)

#if defined(use_libMPI) || defined(use_libSMA)
    ! Handles case where S_msize and/or R_msize are 0 size array
      msgsize = MAXVAL( (/0,sum(d_comm%S_msize(:,:))/) ) + MAXVAL( (/0,sum(d_comm%R_msize(:,:))/) )
      if(msgsize>0)then
         mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, msgsize )
         if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
             write( text,'(i8)' )mpp_domains_stack_hwm
             call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_INIT_COMM: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                  //trim(text)//') from all PEs.' )
         end if
      end if
#endif
            
      DEALLOCATE(slist_addr,isL,jsL)

      d_comm%initialized = .true.

    end function mpp_global_field_init_comm


    subroutine mpp_update_free_comm(domain,l_addr,ksize,lsize,l_addr2,flags,gridtype)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      type(domain2D),target,      intent(in)    :: domain
      integer(LONG_KIND),         intent(in)    :: l_addr
      integer,                    intent(in)    :: ksize
      integer,                    intent(in)    :: lsize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2
      integer, optional,          intent(in)    :: flags
      integer, optional,          intent(in)    :: gridtype

      integer(LONG_KIND) :: domain_id
      integer :: update_flags

      update_flags = XUPDATE+YUPDATE; if( PRESENT(flags) )update_flags = flags
      domain_id=set_domain_id(domain%id,ksize+lsize,update_flags,gridtype)
      call free_comm(domain_id,l_addr,l_addr2)
    end subroutine mpp_update_free_comm


    subroutine mpp_redistribute_free_comm(domain_in,l_addr,domain_out,l_addr2,ksize,lsize)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      type(domain2D),             intent(in)    :: domain_in
      integer(LONG_KIND),         intent(in)    :: l_addr
      type(domain2D),             intent(in)    :: domain_out  
      integer(LONG_KIND),         intent(in)    :: l_addr2
      integer,                    intent(in)    :: ksize,lsize

      integer(LONG_KIND) :: domain_id

      if(l_addr2 > 0)then
        domain_id = set_domain_id(domain_out%id,ksize+lsize)
      else
        domain_id = set_domain_id(domain_in%id,ksize+lsize)
      endif
      call free_comm(domain_id,l_addr,l_addr2)
    end subroutine mpp_redistribute_free_comm


    subroutine mpp_global_field_free_comm(domain,l_addr,ksize,l_addr2,flags)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      type(domain2D),             intent(in)    :: domain
      integer(LONG_KIND),         intent(in)    :: l_addr
      integer,                    intent(in)    :: ksize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2
      integer,           optional,intent(in)    :: flags

      integer :: update_flags
      integer(LONG_KIND) :: domain_id

      update_flags=0; if(PRESENT(flags))update_flags=flags
      domain_id=set_domain_id(domain%id,ksize,update_flags)
      call free_comm(domain_id,l_addr,l_addr2)
    end subroutine mpp_global_field_free_comm


    subroutine free_comm(domain_id,l_addr,l_addr2)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      integer(LONG_KIND),         intent(in) :: domain_id
      integer(LONG_KIND),         intent(in) :: l_addr
      integer(LONG_KIND),optional,intent(in) :: l_addr2

      integer(LONG_KIND) :: dc_key,a_key
      integer :: i,dc_idx,a_idx,i_idx,insert,insert_a,insert_i
      integer :: a2_idx,insert_a2


      i_idx = find_key(domain_id,ids_sorted(1:n_ids),insert_i)
      a_idx = find_key(l_addr,addrs_sorted(1:a_sort_len),insert_a)
      a_key = int(addrs_idx(a_idx),KIND(LONG_KIND))
      if(PRESENT(l_addr2))then
        a2_idx = find_key(l_addr2,addrs2_sorted(1:a2_sort_len),insert_a2)
        a_key = a_key + ADDR2_BASE*int(addrs2_idx(a2_idx),KIND(LONG_KIND))
      endif
      dc_key = DOMAIN_ID_BASE*int(ids_idx(i_idx),KIND(LONG_KIND)) + a_key
      dc_idx = find_key(dc_key,dcKey_sorted(1:dc_sort_len),insert)

      if(dc_idx < 0)then
        call mpp_error(FATAL,'FREE_COMM: attempt to remove nonexistent domains communicator key')
      endif
      call deallocate_comm(d_comm(dc_idx))
      call pop_key(dcKey_sorted,d_comm_idx,dc_sort_len,dc_idx)
      call pop_key(addrs_sorted,addrs_idx,a_sort_len,a_idx)
      if(PRESENT(l_addr2))call pop_key(addrs2_sorted,addrs2_idx,a2_sort_len,a2_idx)
    end subroutine free_comm


    function get_comm(domain_id,l_addr,l_addr2)
      integer(LONG_KIND),intent(in)       :: domain_id
      integer(LONG_KIND),intent(in)       :: l_addr
      integer(LONG_KIND),intent(in),optional :: l_addr2
      type(DomainCommunicator2D), pointer :: get_comm

      integer(LONG_KIND) :: dc_key,a_key
      integer :: i,dc_idx,a_idx,i_idx,insert,insert_a,insert_i
      integer :: a2_idx,insert_a2

      if(.not.ALLOCATED(d_comm))ALLOCATE(d_comm(MAX_FIELDS))
      i_idx   = find_key(domain_id,ids_sorted(1:n_ids),insert_i)
      a_idx = find_key(l_addr,addrs_sorted(1:a_sort_len),insert_a)
      a_key = int(addrs_idx(a_idx),KIND(LONG_KIND))
      if(PRESENT(l_addr2))then
        a2_idx = find_key(l_addr2,addrs2_sorted(1:a2_sort_len),insert_a2)
        a_key = a_key + ADDR2_BASE*int(addrs2_idx(a2_idx),KIND(LONG_KIND))
      endif
      dc_key   = DOMAIN_ID_BASE*int(ids_idx(i_idx),KIND(LONG_KIND)) + a_key
      dc_idx   = find_key(dc_key,dcKey_sorted(1:dc_sort_len),insert)
      if(dc_idx > 0)then
        get_comm =>d_comm(d_comm_idx(dc_idx))
      else
        if(i_idx<0)then
          if(n_ids == MAX_DOM_IDS)then
            call mpp_error(FATAL,'GET_COMM: Maximum number of domains exceeded')
          endif
          n_ids = n_ids+1
          i_idx = push_key(ids_sorted,ids_idx,i_sort_len,insert_i,domain_id,n_ids)
        endif
        if(a_idx<0)then
          if(n_addrs == MAX_ADDRS)then
             call mpp_error(FATAL,'GET_COMM: Maximum number of memory addresses exceeded')
          endif
          n_addrs = n_addrs + 1
          a_idx = push_key(addrs_sorted,addrs_idx,a_sort_len,insert_a,l_addr,n_addrs)
        endif
        if(PRESENT(l_addr2))then
          if(a2_idx<0)then
            if(n_addrs2 == MAX_ADDRS2)then
               call mpp_error(FATAL,'GET_COMM: Maximum number of 2nd memory addresses exceeded')
            endif
            n_addrs2 = n_addrs2 + 1
            a2_idx = push_key(addrs2_sorted,addrs2_idx,a2_sort_len,insert_a2,l_addr2,n_addrs2)
          endif
        endif
        if(n_comm == MAX_FIELDS)then
          call mpp_error(FATAL,'GET_COMM: Maximum number of fields exceeded')
        endif
        a_key = int(addrs_idx(a_idx),KIND(8))
        if(PRESENT(l_addr2))a_key = a_key + ADDR2_BASE*int(addrs2_idx(a2_idx),KIND(8))
        dc_key   = DOMAIN_ID_BASE*int(ids_idx(i_idx),KIND(LONG_KIND)) + a_key
        dc_idx   = find_key(dc_key,dcKey_sorted(1:dc_sort_len),insert)
        if(dc_idx /= -1)call mpp_error(FATAL,'GET_COMM: attempt to insert existing key')
        n_comm = n_comm + 1
        i = push_key(dcKey_sorted,d_comm_idx,dc_sort_len,insert,dc_key,n_comm)
        d_comm_idx(insert) = n_comm
        if(PRESENT(l_addr2))then
          d_comm(n_comm)%l_addrx = l_addr
          d_comm(n_comm)%l_addry = l_addr2
        else
          d_comm(n_comm)%l_addr = l_addr
        endif
        get_comm =>d_comm(n_comm)
      endif
    end function get_comm


    function push_key(sorted,idx,n_idx,insert,key,ival)
      integer(LONG_KIND),intent(inout),dimension(:) :: sorted
      integer,           intent(inout),dimension(:) :: idx
      integer,           intent(inout)              :: n_idx
      integer,           intent(in)                 :: insert
      integer(LONG_KIND),intent(in)                 :: key
      integer,           intent(in)                 :: ival

      integer                                       :: push_key,i

      do i=n_idx,insert,-1
        sorted(i+1) = sorted(i)
        idx(i+1)    = idx(i)
      end do
      sorted(insert) = key
      n_idx = n_idx + 1
      idx(insert) = ival
      push_key = insert
    end function push_key


    subroutine pop_key(sorted,idx,n_idx,key_idx)
      integer(LONG_KIND),intent(inout),dimension(:) :: sorted
      integer,           intent(inout),dimension(:) :: idx
      integer,           intent(inout)              :: n_idx
      integer,           intent(in)                 :: key_idx

      integer                                       :: i

      do i=key_idx,n_idx-1
        sorted(i) = sorted(i+1)
        idx(i)    = idx(i+1)  
      end do                  
      sorted(n_idx) = -9999 
      idx(n_idx) = -9999
      n_idx = n_idx - 1  
    end subroutine pop_key


    function find_key(key,sorted,insert) RESULT(n)
    ! The algorithm used here requires monotonic keys w/out repetition.
      integer(LONG_KIND),intent(in)              :: key        ! new address to be found in list
      integer(LONG_KIND),dimension(:),intent(in) :: sorted  ! list of sorted local addrs
      integer,                        intent(out) :: insert
      integer :: n, n_max, n_min, n_key
      logical :: not_found

      n_key = size(sorted(:))
      insert = 1
      n = -1  ! value not in list
      if(n_key == 0)return  ! first call
      
      if(key < sorted(1))then
        insert = 1; return
      elseif(key > sorted(n_key))then
        insert = n_key+1; return
      endif
        
      if(key == sorted(1))then
        n = 1; return
      elseif(key == sorted(n_key))then
        n = n_key; return
      endif
        
      not_found = .true.
      n = n_key/2 + 1
      n_min=1; n_max=n_key
      do while(not_found)
        if(key == sorted(n))then
          not_found = .false.
        elseif(key > sorted(n))then
          if(key < sorted(n+1))then
            insert = n+1; exit
          endif
          n_min = n
          n = (n+1+n_max)/2
        else
          if(key > sorted(n-1))then
            insert = n; exit
          endif
          n_max = n
          n = (n+n_min)/2
        endif
        if(n==1 .or. n==n_key)exit
      end do
      if(not_found)n = -1  ! value not in list
    end function find_key


    subroutine deallocate_comm(d_comm)
      type(DomainCommunicator2D), intent(inout) :: d_comm

      d_comm%domain  =>NULL()
      d_comm%domain_in  =>NULL()
      d_comm%domain_out =>NULL()

      d_comm%initialized=.false.
      d_comm%id=-9999
      d_comm%l_addr  =-9999
      d_comm%l_addrx =-9999
      d_comm%l_addry =-9999

      if( ASSOCIATED(d_comm%sendis) )   DEALLOCATE(d_comm%sendis);    d_comm%sendis =>NULL()
      if( ASSOCIATED(d_comm%sendie) )   DEALLOCATE(d_comm%sendie);    d_comm%sendie =>NULL()
      if( ASSOCIATED(d_comm%sendjs) )   DEALLOCATE(d_comm%sendjs);    d_comm%sendjs =>NULL()
      if( ASSOCIATED(d_comm%sendje) )   DEALLOCATE(d_comm%sendje);    d_comm%sendje =>NULL()
      if( ASSOCIATED(d_comm%S_msize) )  DEALLOCATE(d_comm%S_msize);   d_comm%S_msize =>NULL()
      if( ASSOCIATED(d_comm%do_thisS) ) DEALLOCATE(d_comm%do_thisS);  d_comm%do_thisS =>NULL()
      if( ASSOCIATED(d_comm%S_do_buf) ) DEALLOCATE(d_comm%S_do_buf);  d_comm%S_do_buf =>NULL()
      if( ASSOCIATED(d_comm%cto_pe) )   DEALLOCATE(d_comm%cto_pe);    d_comm%cto_pe  =>NULL()
      if( ASSOCIATED(d_comm%Rcaf_idx) ) DEALLOCATE(d_comm%Rcaf_idx);  d_comm%Rcaf_idx  =>NULL()
      if( ASSOCIATED(d_comm%recvis) )   DEALLOCATE(d_comm%recvis);    d_comm%recvis =>NULL()
      if( ASSOCIATED(d_comm%recvie) )   DEALLOCATE(d_comm%recvie);    d_comm%recvie =>NULL()
      if( ASSOCIATED(d_comm%recvjs) )   DEALLOCATE(d_comm%recvjs);    d_comm%recvjs =>NULL()
      if( ASSOCIATED(d_comm%recvje) )   DEALLOCATE(d_comm%recvje);    d_comm%recvje =>NULL()
      if( ASSOCIATED(d_comm%R_msize) )  DEALLOCATE(d_comm%R_msize);   d_comm%R_msize =>NULL()
      if( ASSOCIATED(d_comm%do_thisR) ) DEALLOCATE(d_comm%do_thisR);  d_comm%do_thisR =>NULL()
      if( ASSOCIATED(d_comm%R_do_buf) ) DEALLOCATE(d_comm%R_do_buf);  d_comm%R_do_buf =>NULL()
      if( ASSOCIATED(d_comm%cfrom_pe) ) DEALLOCATE(d_comm%cfrom_pe);  d_comm%cfrom_pe  =>NULL()
      d_comm%Slist_size=0; d_comm%Rlist_size=0
      d_comm%isize=0; d_comm%jsize=0; d_comm%ke=0
      d_comm%isize_in=0; d_comm%jsize_in=0
      d_comm%isize_out=0; d_comm%jsize_out=0
      d_comm%isize_max=0; d_comm%jsize_max=0
      d_comm%gf_ioff=0; d_comm%gf_joff=0
    ! Remote data
      if( ASSOCIATED(d_comm%isizeR) )   DEALLOCATE(d_comm%isizeR);    d_comm%isizeR =>NULL()
      if( ASSOCIATED(d_comm%jsizeR) )   DEALLOCATE(d_comm%jsizeR);    d_comm%jsizeR =>NULL()
      if( ASSOCIATED(d_comm%sendisR) )  DEALLOCATE(d_comm%sendisR);   d_comm%sendisR =>NULL()
      if( ASSOCIATED(d_comm%sendjsR) )  DEALLOCATE(d_comm%sendjsR);   d_comm%sendjsR =>NULL()
      if( ASSOCIATED(d_comm%rem_addr) ) DEALLOCATE(d_comm%rem_addr);  d_comm%rem_addr  =>NULL()
      if( ASSOCIATED(d_comm%rem_addrx) )DEALLOCATE(d_comm%rem_addrx); d_comm%rem_addrx =>NULL()
      if( ASSOCIATED(d_comm%rem_addry) )DEALLOCATE(d_comm%rem_addry); d_comm%rem_addry =>NULL()
      if( ASSOCIATED(d_comm%rem_addrl) )DEALLOCATE(d_comm%rem_addrl); d_comm%rem_addrl =>NULL()
    end subroutine deallocate_comm


    function set_domain_id(d_id,ksize,flags,gtype)
      integer(LONG_KIND), intent(in) :: d_id
      integer           , intent(in) :: ksize
      integer           , optional, intent(in) :: flags
      integer           , optional, intent(in) :: gtype

      integer(LONG_KIND)             :: set_domain_id
      
      set_domain_id=d_id + KE_BASE*int(ksize,KIND(d_id))
      if(PRESENT(flags))set_domain_id=set_domain_id+int(flags,KIND(d_id))
      if(PRESENT(gtype))set_domain_id=set_domain_id+GT_BASE*int(gtype,KIND(d_id))  ! Must be LONG_KIND arithmetic
    end function set_domain_id


#ifdef use_GSM
    subroutine check_sma_env()
      character(len=15) :: sma_env
      call getenv('SMA_GLOBAL_HEAP_SIZE',sma_env)
      if(sma_env == '               ') &
         call mpp_error( FATAL, 'Environment variables SMA_GLOBAL_ALLOC and SMA_GLOBAL_HEAP_SIZE must be set.' )
    end subroutine check_sma_env
#endif

!#######################################################################


#ifdef use_CAF
#define MPP_TYPE_ real(DOUBLE_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_r8_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_r8_1d
#define CAFPNTR_TYPE_3D_ cafptr_r8_3d_type
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_r8_3d
#include <mpp_domains_comm.h>

#define MPP_TYPE_ complex(DOUBLE_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_c8_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_c8_1d
#define CAFPNTR_TYPE_3D_ cafptr_c8_3d_type
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_c8_3d
#include <mpp_domains_comm.h>

#ifndef no_8byte_integers
#define MPP_TYPE_ integer(LONG_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_i8_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_i8_1d
#define CAFPNTR_TYPE_3D_ cafptr_i8_3d_type
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_i8_3d
#include <mpp_domains_comm.h>

#define MPP_TYPE_ logical(LONG_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_l8_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_l8_1d
#define CAFPNTR_TYPE_3D_ cafptr_l8_3d_type
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_l8_3d
#include <mpp_domains_comm.h>
#endif

#ifndef no_4byte_reals
#define MPP_TYPE_ real(FLOAT_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_r4_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_r4_1d
#define CAFPNTR_TYPE_3D_ cafptr_r4_3d_type
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_r4_3d
#include <mpp_domains_comm.h>

#define MPP_TYPE_ complex(FLOAT_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_c4_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_c4_1d
#define CAFPNTR_TYPE_3D_ cafptr_c4_3d_type
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_c4_3d
#include <mpp_domains_comm.h>
#endif

#define MPP_TYPE_ integer(INT_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_i4_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_i4_1d
#define CAFPNTR_TYPE_3D_ cafptr_i4_3d_type
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_i4_3d
#include <mpp_domains_comm.h>

#define MPP_TYPE_ logical(INT_KIND)
!#define CAFPNTR_TYPE_1D_ cafptr_l4_1d_type
!#define MPP_ASSOC_CAF_FIELD_1D_ associate_caf_field_l4_1d
#define CAFPNTR_TYPE_3D_ cafptr_l4_3d_type
#define MPP_ASSOC_CAF_FIELD_3D_ associate_caf_field_l4_3d
#include <mpp_domains_comm.h>                                                                      
#endif

end module mpp_domains_comm_mod
