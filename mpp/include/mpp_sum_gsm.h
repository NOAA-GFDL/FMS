    subroutine MPP_SUM_( a, length, pelist )
!sums array a over the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast: all PEs have the sum in a at the end
  !we are using f77-style call: array passed by address and not descriptor; further, 
  !the f90 conformance check is avoided.
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_, intent(inout) :: a(*)

      integer :: list,k,m,n
      integer, parameter :: MPP_SUM_BUFSIZE=1024
!     integer, parameter :: MPP_SUM_BUFSIZE=1048576
      MPP_TYPE_,volatile :: r_work(MPP_SUM_BUFSIZE)
      MPP_TYPE_ :: l_work(length)
      pointer( r_ptr,r_work )
      MPP_TYPE_,dimension(MPP_SUM_BUFSIZE), save :: work
      integer(LONG_KIND),allocatable,save :: work_addrs(:)


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      if(length > MPP_SUM_BUFSIZE)call mpp_error( FATAL, 'MPP_SUM: Sum lengths exceeds MPP_SUM_BUFSIZE.' )
      if(.not.ALLOCATED(work_addrs))then
        allocate(work_addrs(0:npes-1))
        do list=0,npes-1
          work_addrs(list) = shmem_ptr(work,list)
        end do
      endif

      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)

      l_work = 0.d0
      if(n==world_peset_num)then
         if(length == 1)then
           work(1) = a(1)  ! must be loaded prior to barrier
           call mpp_sync(pelist,do_self=.false.)
           do k=0,npes-1
             r_ptr = work_addrs(k)
             l_work(1) = l_work(1) + r_work(1)
           end do
         else
           work(1:length) = a(1:length)  ! must be loaded prior to barrier
           call mpp_sync(pelist,do_self=.false.)
           do k=0,npes-1
             r_ptr = work_addrs(k)
             do m=1,length
               l_work(m) = l_work(m) + r_work(m)
             end do
           end do
         endif

      else

        list = size(peset(n)%list(:))
        if(length == 1)then
           work(1) = a(1)  ! must be loaded prior to barrier
           call mpp_sync(pelist,do_self=.false.)
           do k=0,list-1
             r_ptr = work_addrs(peset(n)%list(k+1))
             l_work(1) = l_work(1) + r_work(1)
           end do
         else
           work(1:length) = a(1:length)  ! must be loaded prior to barrier
           call mpp_sync(pelist,do_self=.false.)
           do k=0,list-1
             r_ptr = work_addrs(peset(n)%list(k+1))
             do m=1,length
               l_work(m) = l_work(m) + r_work(m)
             end do
           end do
         endif
      endif
      a(1:length) = l_work(1:length)

      call mpp_sync(pelist,do_self=.false.)
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, length*MPP_TYPE_BYTELEN_ )
    end subroutine MPP_SUM_

!#######################################################################
#include <mpp_sum.inc>
