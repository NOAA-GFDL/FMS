    subroutine MPP_REDUCE_( a, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      integer :: list,m,n
      MPP_TYPE_,volatile :: r_work(1)
      MPP_TYPE_,save :: l_work
      pointer( r_ptr,r_work )
      MPP_TYPE_, save :: work
      integer(LONG_KIND),allocatable,save :: work_addrs(:)
      integer :: done_cnt
      logical, dimension(0:npes-1) :: done


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE: You must first call mpp_init.' )
      if(.not.ALLOCATED(work_addrs))then
        allocate(work_addrs(0:npes-1))
        do list=0,npes-1
          work_addrs(list) = shmem_ptr(work,list)
        end do
      endif

      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)

      l_work = a; work = a
      call mpp_sync(pelist,do_self=.false.)
      if(n==world_peset_num)then
!       done = .false.; done_cnt = npes
!       gsm_sync_start = Set_Sync_Logical(.true.)
!       do while(done_cnt /= 0); do m=0,npes-1
        do m=0,npes-1
!         if(done(m))cycle
!         if( Is_Remote_False(gsm_sync_start_list(m)) )cycle
!         done(m) = .true.; done_cnt = done_cnt - 1
          r_ptr = work_addrs(m)
          if(r_work(1) GSM_REDUCE_ l_work)l_work = r_work(1)
        end do
!       end do; end do
      else
        list = size(peset(n)%list(:))
!       done = .false.; done_cnt = list
!       gsm_sync_start = Set_Sync_Logical(.true.)
!       do while(done_cnt /= 0); do m=0,list-1
        do m=0,list-1
!         if(done(m))cycle
!         if( Is_Remote_False(gsm_sync_start_list(peset(n)%list(m+1))) )cycle
!         done(m) = .true.; done_cnt = done_cnt - 1
          r_ptr = work_addrs(peset(n)%list(m+1))
          if(r_work(1) GSM_REDUCE_ l_work)l_work = r_work(1)
        end do
!       end do; end do
      endif
      a = l_work

      call mpp_sync(pelist,do_self=.false.)
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, MPP_TYPE_BYTELEN_ )
!     gsm_sync_start = Set_Sync_Logical(.false.)

    end subroutine MPP_REDUCE_
