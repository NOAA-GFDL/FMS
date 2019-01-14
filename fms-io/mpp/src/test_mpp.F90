!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
#ifdef test_mpp
#ifdef SYSTEM_CLOCK
#undef SYSTEM_CLOCK
#endif

program test   !test various aspects of mpp_mod
#include <fms_platform.h>

#ifdef sgi_mipspro
  use shmem_interface
#endif

  use mpp_mod, only : mpp_init, mpp_exit, mpp_pe, mpp_npes, mpp_root_pe, stdout
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sync, mpp_malloc
  use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist, mpp_set_stack_size
  use mpp_mod, only : mpp_broadcast, mpp_transmit, mpp_sum, mpp_max, mpp_chksum, ALL_PES
  use mpp_mod, only : mpp_gather, mpp_error, FATAL, mpp_sync_self
  use mpp_io_mod, only: mpp_io_init, mpp_flush
#ifdef use_MPI_GSM
  use mpp_mod, only : mpp_gsm_malloc, mpp_gsm_free
#endif

  implicit none

  integer, parameter              :: n=1048576
  real, allocatable, dimension(:) :: a, b, c
#ifdef use_MPI_GSM
  real                            :: d(n)
  pointer (locd, d)
#else
  real, allocatable, dimension(:) :: d
  integer(LONG_KIND) :: locd
#endif
  integer                         :: tick, tick0, ticks_per_sec, id
  integer                         :: pe, npes, root, i, j, k, l, m, n2, istat
  integer                         :: out_unit
  real                            :: dt

  call mpp_init()
  call mpp_io_init()
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()

  out_unit = stdout()
  call test_gather(npes,pe,root,out_unit)
  call test_gatherV(npes,pe,root,out_unit)
  call test_gather2DV(npes,pe,root,out_unit)

  if(.false.) then

  ! first test broadcast
  call test_broadcast()

  call SYSTEM_CLOCK( count_rate=ticks_per_sec )

  allocate( a(n), b(n) )
  id = mpp_clock_id( 'Random number' )
  call mpp_clock_begin(id)
  call random_number(a)
  call mpp_clock_end  (id)
  !---------------------------------------------------------------------!
  !   time transmit, compare against shmem_put and get                  !
  !---------------------------------------------------------------------!
  if( pe.EQ.root )then
     print *, 'Time mpp_transmit for various lengths...'
#ifdef SGICRAY
     print *, 'For comparison, times for shmem_get and shmem_put are also provided.'
#endif
     print *
  end if
  id = mpp_clock_id( 'mpp_transmit' )
  call mpp_clock_begin(id)
  !timing is done for cyclical pass (more useful than ping-pong etc)
  l = n
  do while( l.GT.0 )
     !--- mpp_transmit -------------------------------------------------
     call mpp_sync()
     call SYSTEM_CLOCK(tick0)
     do i = 1,npes
        call mpp_transmit( put_data=a(1), plen=l, to_pe=modulo(pe+npes-i,npes), &
                           get_data=b(1), glen=l, from_pe=modulo(pe+i,npes) )
        call mpp_sync_self()
        !          call mpp_sync_self( (/modulo(pe+npes-i,npes)/) )
     end do
     call mpp_sync()
     call SYSTEM_CLOCK(tick)
     dt = real(tick-tick0)/(npes*ticks_per_sec)
     dt = max( dt, epsilon(dt) )
     if( pe.EQ.root )write( out_unit,'(/a,i8,f13.6,f8.2)' )'MPP_TRANSMIT length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!#ifdef SGICRAY
!     !--- shmem_put ----------------------------------------------------
!     call mpp_sync()
!     call SYSTEM_CLOCK(tick0)
!     do i = 1,npes
!       call shmem_real_put( b, a, l, modulo(pe+1,npes) )
!     end do
!     call mpp_sync()
!     call SYSTEM_CLOCK(tick)
!     dt = real(tick-tick0)/(npes*ticks_per_sec)
!     dt = max( dt, epsilon(dt) )
!     if( pe.EQ.root )write( stdout(),'( a,i8,f13.6,f8.2)' )'SHMEM_PUT    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!     !--- shmem_get ----------------------------------------------------
!     call mpp_sync()
!     call SYSTEM_CLOCK(tick0)
!     do i = 1,npes
!        call shmem_real_get( b, a, l, modulo(pe+1,npes) )
!     end do
!     call SYSTEM_CLOCK(tick)
!     dt = real(tick-tick0)/(npes*ticks_per_sec)
!     dt = max( dt, epsilon(dt) )
!     if( pe.EQ.root )write( stdout(),'( a,i8,f13.6,f8.2)' )'SHMEM_GET    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!#endif
     l = l/2
  end do
  !---------------------------------------------------------------------!
  !                   test mpp_sum                                      !
  !---------------------------------------------------------------------!
  if( pe.EQ.root )then
     print '(/a)', 'Time mpp_sum...'
  end if
  a = real(pe+1)
  call mpp_sync()
  call SYSTEM_CLOCK(tick0)
  call mpp_sum(a(1:1000),1000)
  call SYSTEM_CLOCK(tick)
  dt = real(tick-tick0)/ticks_per_sec
  dt = max( dt, epsilon(dt) )
  if( pe.EQ.root )write( out_unit,'(a,2i6,f9.1,i8,f13.6,f8.2/)' ) &
       'mpp_sum: pe, npes, sum(pe+1), length, time, bw(Mb/s)=', pe, npes, a(1), n, dt, n*8e-6/dt
  call mpp_clock_end(id)
  !---------------------------------------------------------------------!
  !                   test mpp_max                                      !
  !---------------------------------------------------------------------!
  if( pe.EQ.root )then
     print *
     print *, 'Test mpp_max...'
  end if
  a = real(pe+1)
  print *, 'pe,     pe+1 =', pe, a(1)
  call mpp_max( a(1) )
  print *, 'pe, max(pe+1)=', pe, a(1)
  !pelist check
  call mpp_sync()
  call flush(out_unit,istat)
  if( npes.GE.2 )then
     if( pe.EQ.root )print *, 'Test of pelists: bcast, sum and max using PEs 0...npes-2 (excluding last PE)'
     call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
     a = real(pe+1)
     if( pe.NE.npes-1 )call mpp_broadcast( a, n, npes-2, (/(i,i=0,npes-2)/) )
     print *, 'bcast(npes-1) from 0 to npes-2=', pe, a(1)
     a = real(pe+1)
     if( pe.NE.npes-1 )then
        call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
        id = mpp_clock_id( 'Partial mpp_sum' )
        call mpp_clock_begin(id)
        call mpp_sum( a(1:1000), 1000, (/(i,i=0,npes-2)/) )
        call mpp_clock_end  (id)
     end if
     if( pe.EQ.root )print *, 'sum(pe+1) from 0 to npes-2=', a(1)
     a = real(pe+1)
     if( pe.NE.npes-1 )call mpp_max( a(1), (/(i,i=0,npes-2)/) )
     if( pe.EQ.root )print *, 'max(pe+1) from 0 to npes-2=', a(1)
  end if
  call mpp_set_current_pelist()
  
#ifdef use_CRI_pointers
  !---------------------------------------------------------------------!
  !                   test mpp_chksum                                   !
  !---------------------------------------------------------------------!
  if( modulo(n,npes).EQ.0 )then  !only set up for even division
     n2 = 1024
     a = 0.d0
     if( pe.EQ.root )call random_number(a(1:n2))
!    if( pe.EQ.root )call random_number(a)
     call mpp_sync()
     call mpp_transmit( put_data=a(1), plen=n2, to_pe=ALL_PES, &
                        get_data=a(1), glen=n2, from_pe=root )
     call mpp_sync_self ()
!    call mpp_transmit( put_data=a(1), plen=n, to_pe=ALL_PES, &
!                       get_data=a(1), glen=n, from_pe=root )
     m= n2/npes
!    m= n/npes
     allocate( c(m) )
     c = a(pe*m+1:pe*m+m)
     
     if( pe.EQ.root )then
        print *
        print *, 'Test mpp_chksum...'
        print *, 'This test shows that a whole array and a distributed array give identical checksums.'
     end if
     print *, 'chksum(a(1:1024))=', mpp_chksum(a(1:n2),(/pe/))
     print *, 'chksum(c(1:1024))=', mpp_chksum(c)
!    print *, 'chksum(a)=', mpp_chksum(a,(/pe/))
!    print *, 'chksum(c)=', mpp_chksum(c)
  end if
!test of pointer sharing
#ifdef use_MPI_GSM
      call mpp_gsm_malloc( locd, sizeof(d) )
#else
  if( pe.EQ.root )then
      allocate( d(n) )
      locd = LOC(d)
  end if
  call mpp_broadcast(locd,root)
#endif
  if( pe.EQ.root )then
      call random_number(d)
  end if
  call mpp_sync()
!  call test_shared_pointers(locd,n)

#ifdef use_MPI_GSM
  call mpp_gsm_free( locd )
#else
  if( pe.EQ.root )then
      deallocate( d )
  end if
#endif
#endif
  endif  ! if(.false.)

  call mpp_exit()

contains

  !***********************************************
  !currently only test the mpp_broadcast_char
  subroutine test_broadcast()
     integer, parameter :: ARRAYSIZE = 3
     integer, parameter :: STRINGSIZE = 256
     character(len=STRINGSIZE), dimension(ARRAYSIZE) :: textA, textB
     integer :: n

     textA(1) = "This is line 1 "
     textA(2) = "Here comes the line 2 "
     textA(3) = "Finally is line 3 "  
     do n = 1, ARRAYSIZE  
        textB(n) = TextA(n)
     enddo

     if(mpp_pe() .NE. mpp_root_pe()) then
        do n =1, ARRAYSIZE
           textA(n) = ""
        enddo
     endif

     !--- comparing textA and textB. textA and textB are supposed to be different on pe other than root_pe
     if(mpp_pe() == mpp_root_pe()) then
        do n = 1, ARRAYSIZE         
           if(textA(n) .NE. textB(n)) call mpp_error(FATAL, "test_broadcast: on root_pe, textA should equal textB")
        enddo
     else
        do n = 1, ARRAYSIZE         
           if(textA(n) == textB(n)) call mpp_error(FATAL, "test_broadcast: on root_pe, textA should not equal textB")
        enddo 
     endif
     call mpp_broadcast(textA, STRINGSIZE, mpp_root_pe())
     !--- after broadcast, textA and textB should be the same
     do n = 1, ARRAYSIZE         
        if(textA(n) .NE. textB(n)) call mpp_error(FATAL, "test_broadcast: after broadcast, textA should equal textB")
     enddo

     write(out_unit,*) "==> NOTE from test_broadcast: The test is succesful"

  end subroutine test_broadcast

  subroutine test_gather(npes,pe,root,out_unit)
     integer, intent(in) :: npes,pe,root,out_unit

     integer :: pelist(npes)
     integer :: i
     real :: rdata(npes)
     real :: val

     if(npes < 3)then
       write(out_unit,*) "Minimum of 3 ranks required. Not testing gather; too few ranks."
       return
     endif
     write(out_unit,*)

     val = pe
     rdata = -1.0
     do i=1,npes
       pelist(i) = i-1
     enddo

     call mpp_gather((/val/),rdata)
     if(pe == root)then 
       do i=1,npes
        if(INT(rdata(i)) /= pelist(i))then
           write(6,*) "Gathered data ",INT(rdata(i)), " NE reference ",pelist(i), "at i=",i
           call mpp_error(FATAL, "Test gather uniform vector with global pelist failed")
        endif
       enddo
     endif

     call mpp_sync()
     write(out_unit,*) "Test gather uniform vector with global pelist successful"

     rdata = -1.0
     if(ANY(pe == pelist(2:npes)))call mpp_gather((/val/),rdata(2:npes),pelist(2:npes))
     if(pe == pelist(2))then
       do i=2,npes
        if(INT(rdata(i)) /= pelist(i))then
           write(6,*) "Gathered data ",INT(rdata(i)), " NE reference ",pelist(i), "at i=",i
           call mpp_error(FATAL, "Test gather uniform vector with reduced pelist failed")
        endif
       enddo
     endif
     call mpp_sync()
     write(out_unit,*) "Test gather uniform vector with reduced pelist successful"

  end subroutine test_gather


  subroutine test_gatherV(npes,pe,root,out_unit)
  implicit none
     integer, intent(in) :: npes,pe,root,out_unit

     integer :: pelist(npes),rsize(npes)
     integer :: i,j,k,dsize,ssize
     real,allocatable :: sdata(:), rdata(:), ref(:)

     if(npes < 3)then
       write(out_unit,*) "Minimum of 3 ranks required. Not testing gatherV; too few ranks."
       return
     elseif(npes > 9999)then
       write(out_unit,*) "Maximum of 9999 ranks supported. Not testing gatherV; too many ranks."
       return
     endif
     write(out_unit,*)

     ssize = pe+1
     allocate(sdata(ssize))
     do i=1,ssize
       sdata(i) = pe + 0.0001*i
     enddo
     do i=1,npes
       pelist(i) = i-1
       rsize(i) = i
     enddo

     dsize = sum(rsize)
     allocate(rdata(dsize),ref(dsize))
     rdata = -1.0
     k=1
     do j=1,npes
       do i=1,rsize(j)
          ref(k) = pelist(j) + 0.0001*i
          k = k+1
     enddo;enddo

     call mpp_gather(sdata,ssize,rdata,rsize)
 
     if(pe == root)then
       k = 1
       do j=1,npes
         do i=1,rsize(j)
           if(rdata(k) /= ref(k))then
              write(6,*) "Gathered data ",rdata(k), " NE reference ",ref(k), "at k=",k
              call mpp_error(FATAL, "Test gatherV global pelist failed")
           endif
           k = k+1
       enddo;enddo
     endif

     call mpp_sync()
     write(out_unit,*) "Test gatherV with global pelist successful"

     rdata = -1.0
     ref(1) = -1.0

     if(ANY(pe == pelist(2:npes)))call mpp_gather(sdata,ssize,rdata(2:),rsize(2:),pelist(2:npes))

     if(pe == pelist(2))then
       k = 1
       do j=1,npes
         do i=1,rsize(j)
           if(rdata(k) /= ref(k))then
              write(6,*) "Gathered data ",rdata(k), " NE reference ",ref(k), "at k=",k
              call mpp_error(FATAL, "Test gatherV with reduced pelist failed")
           endif
           k = k+1
       enddo;enddo
     endif
     call mpp_sync()

     write(out_unit,*) "Test gatherV with reduced pelist successful"
     deallocate(sdata,rdata,ref)
  end subroutine test_gatherV

subroutine test_gather2DV(npes,pe,root,out_unit)
  implicit none
     integer, intent(in) :: npes,pe,root,out_unit

     integer :: pelist(npes),rsize(npes)
     integer :: pelist2(npes),rsize2(npes)
     integer :: i,j,k,l,nz,ssize,nelems
     real,allocatable,dimension(:,:) :: data, cdata, sbuff,rbuff
     real,allocatable :: ref(:,:)
     integer, parameter :: KSIZE=10

     real :: sbuff1D(size(sbuff))
     real :: rbuff1D(size(rbuff))
     pointer(sptr,sbuff1D); pointer(rptr,rbuff1D)


     if(npes < 3)then
       write(out_unit,*) "Minimum of 3 ranks required. Not testing gather2DV; too few ranks."
       return
     elseif(npes > 9999)then
       write(out_unit,*) "Maximum of 9999 ranks supported. Not testing gather2DV; too many ranks."
       return
     endif
     write(out_unit,*)

     ssize = pe+1
     allocate(data(ssize,KSIZE))
     do k=1,KSIZE; do i=1,ssize
       data(i,k) = 10000.0*k + pe + 0.0001*i
     enddo; enddo
     do i=1,npes
       pelist(i) = i-1
       rsize(i) = i
     enddo

     nz = KSIZE
     nelems = sum(rsize(:))

     allocate(rbuff(nz,nelems)); rbuff = -1.0
     allocate(ref(nelems,nz),cdata(nelems,nz))
     ref = 0.0; cdata = 0.0
     if(pe == root)then
       do k=1,KSIZE
       l=1
       do j=1,npes
         do i=1,rsize(j)
            ref(l,k) = 10000.0*k + pelist(j) + 0.0001*i
            l = l+1
       enddo; enddo;enddo
     endif
     allocate(sbuff(nz,ssize))
     ! this matrix inversion makes for easy gather to the IO root
     ! and a clear, concise unpack
     do j=1,ssize
       do i=1,nz
         sbuff(i,j) = data(j,i)
     enddo; enddo

  !  Note that the gatherV implied here is asymmetric; only root needs to know the vector of recv size
     sptr = LOC(sbuff); rptr = LOC(rbuff)
     call mpp_gather(sbuff1D,size(sbuff),rbuff1D,nz*rsize(:))

     if(pe == root)then
        do j=1,nz
           do i=1,nelems
             cdata(i,j) = rbuff(j,i)
        enddo; enddo
        do j=1,nz
           do i=1,nelems
            if(cdata(i,j) /= ref(i,j))then
               write(6,*) "Gathered data ",cdata(i,j), " NE reference ",ref(i,j), "at i,j=",i,j
               call mpp_error(FATAL, "Test gather2DV global pelist failed")
            endif
       enddo;enddo
     endif

     call mpp_sync()
     write(out_unit,*) "Test gather2DV with global pelist successful"

     do i=1,npes
       pelist2(i) = pelist(npes-i+1) 
       rsize2(i) = rsize(npes-i+1)
     enddo

     rbuff = -1.0
     ref = 0.0; cdata = 0.0
     if(pe == pelist2(1))then
       do k=1,KSIZE
       l=1
       do j=1,npes
         do i=1,rsize2(j)
            ref(l,k) = 10000.0*k + pelist2(j) + 0.0001*i
            l = l+1
       enddo; enddo;enddo
     endif

     call mpp_gather(sbuff1D,size(sbuff),rbuff1D,nz*rsize2(:),pelist2)

     if(pe == pelist2(1))then
        do j=1,nz
           do i=1,nelems
             cdata(i,j) = rbuff(j,i)
        enddo; enddo
        do j=1,nz
           do i=1,nelems
            if(cdata(i,j) /= ref(i,j))then
               write(6,*) "Gathered data ",cdata(i,j), " NE reference ",ref(i,j), "at i,j=",i,j
               call mpp_error(FATAL, "Test gather2DV with reversed pelist failed")
            endif
       enddo;enddo
     endif
     call mpp_sync()
     write(out_unit,*) "Test gather2DV with reversed pelist successful"
     deallocate(data,sbuff,rbuff,cdata,ref)
  end subroutine test_gather2DV

  subroutine test_shared_pointers(locd,n)
    integer(LONG_KIND), intent(in) :: locd
    integer :: n
    real :: dd(n)
    pointer( p, dd )

    p = locd
    print *, 'TEST_SHARED_POINTERS: pe, locd=', pe, locd
!    print *, 'TEST_SHARED_POINTERS: pe, chksum(d)=', pe, mpp_chksum(dd,(/pe/))
    print *, 'TEST_SHARED_POINTERS: pe, sum(d)=', pe, sum(dd)
    return
  end subroutine test_shared_pointers
end program test

#else
module null_mpp_test
end module  

#endif /* test_mpp */
