!-----------------------------------------------------------------------
!                 Communication for message-passing codes
!
! AUTHOR: V. Balaji (V.Balaji@noaa.gov)
!         SGI/GFDL Princeton University
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------
module mpp_mod
!a generalized communication package for use with shmem and MPI
!will add: co_array_fortran, MPI2
!Balaji (V.Balaji@noaa.gov) 11 May 1998

! <CONTACT EMAIL="V.Balaji@noaa.gov">
!   V. Balaji
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <RCSLOG SRC="http://www.gfdl.noaa.gov/~vb/changes_mpp.html"/>

! <OVERVIEW>
!   <TT>mpp_mod</TT>, is a set of simple calls to provide a uniform interface
!   to different message-passing libraries. It currently can be
!   implemented either in the SGI/Cray native SHMEM library or in the MPI
!   standard. Other libraries (e.g MPI-2, Co-Array Fortran) can be
!   incorporated as the need arises.
! </OVERVIEW>

! <DESCRIPTION>
!   The data transfer between a processor and its own memory is based
!   on <TT>load</TT> and <TT>store</TT> operations upon
!   memory. Shared-memory systems (including distributed shared memory
!   systems) have a single address space and any processor can acquire any
!   data within the memory by <TT>load</TT> and
!   <TT>store</TT>. The situation is different for distributed
!   parallel systems. Specialized MPP systems such as the T3E can simulate
!   shared-memory by direct data acquisition from remote memory. But if
!   the parallel code is distributed across a cluster, or across the Net,
!   messages must be sent and received using the protocols for
!   long-distance communication, such as TCP/IP. This requires a
!   ``handshaking'' between nodes of the distributed system. One can think
!   of the two different methods as involving <TT>put</TT>s or
!   <TT>get</TT>s (e.g the SHMEM library), or in the case of
!   negotiated communication (e.g MPI), <TT>send</TT>s and
!   <TT>recv</TT>s.
!   
!   The difference between SHMEM and MPI is that SHMEM uses one-sided
!   communication, which can have very low-latency high-bandwidth
!   implementations on tightly coupled systems. MPI is a standard
!   developed for distributed computing across loosely-coupled systems,
!   and therefore incurs a software penalty for negotiating the
!   communication. It is however an open industry standard whereas SHMEM
!   is a proprietary interface. Besides, the <TT>put</TT>s or
!   <TT>get</TT>s on which it is based cannot currently be implemented in
!   a cluster environment (there are recent announcements from Compaq that
!   occasion hope).
!   
!   The message-passing requirements of climate and weather codes can be
!   reduced to a fairly simple minimal set, which is easily implemented in
!   any message-passing API. <TT>mpp_mod</TT> provides this API.
!
!    Features of <TT>mpp_mod</TT> include:
!   
!    1) Simple, minimal API, with free access to underlying API for
!       more complicated stuff.<BR/>
!    2) Design toward typical use in climate/weather CFD codes.<BR/>
!    3) Performance to be not significantly lower than any native API.
!   
!   This module is used to develop higher-level calls for <LINK 
!   SRC="mpp_domains.html">domain decomposition</LINK> and <LINK
!   SRC="mpp_io.html">parallel I/O</LINK>.
!   
!   Parallel computing is initially daunting, but it soon becomes
!   second nature, much the way many of us can now write vector code
!   without much effort. The key insight required while reading and
!   writing parallel code is in arriving at a mental grasp of several
!   independent parallel execution streams through the same code (the SPMD
!   model). Each variable you examine may have different values for each
!   stream, the processor ID being an obvious example. Subroutines and
!   function calls are particularly subtle, since it is not always obvious
!   from looking at a call what synchronization between execution streams
!   it implies. An example of erroneous code would be a global barrier
!   call (see <LINK SRC="#mpp_sync">mpp_sync</LINK> below) placed
!   within a code block that not all PEs will execute, e.g:
!   
!   <PRE>
!   if( pe.EQ.0 )call mpp_sync()
!   </PRE>
!   
!   Here only PE 0 reaches the barrier, where it will wait
!   indefinitely. While this is a particularly egregious example to
!   illustrate the coding flaw, more subtle versions of the same are
!   among the most common errors in parallel code.
!   
!   It is therefore important to be conscious of the context of a
!   subroutine or function call, and the implied synchronization. There
!   are certain calls here (e.g <TT>mpp_declare_pelist, mpp_init,
!   mpp_malloc, mpp_set_stack_size</TT>) which must be called by all
!   PEs. There are others which must be called by a subset of PEs (here
!   called a <TT>pelist</TT>) which must be called by all the PEs in the
!   <TT>pelist</TT> (e.g <TT>mpp_max, mpp_sum, mpp_sync</TT>). Still
!   others imply no synchronization at all. I will make every effort to
!   highlight the context of each call in the MPP modules, so that the
!   implicit synchronization is spelt out.  
!   
!   For performance it is necessary to keep synchronization as limited
!   as the algorithm being implemented will allow. For instance, a single
!   message between two PEs should only imply synchronization across the
!   PEs in question. A <I>global</I> synchronization (or <I>barrier</I>)
!   is likely to be slow, and is best avoided. But codes first
!   parallelized on a Cray T3E tend to have many global syncs, as very
!   fast barriers were implemented there in hardware.
!   
!   Another reason to use pelists is to run a single program in MPMD
!   mode, where different PE subsets work on different portions of the
!   code. A typical example is to assign an ocean model and atmosphere
!   model to different PE subsets, and couple them concurrently instead of
!   running them serially. The MPP module provides the notion of a
!   <I>current pelist</I>, which is set when a group of PEs branch off
!   into a subset. Subsequent calls that omit the <TT>pelist</TT> optional
!   argument (seen below in many of the individual calls) assume that the
!   implied synchronization is across the current pelist. The calls
!   <TT>mpp_root_pe</TT> and <TT>mpp_npes</TT> also return the values
!   appropriate to the current pelist. The <TT>mpp_set_current_pelist</TT>
!   call is provided to set the current pelist.

! </DESCRIPTION>
! <PUBLIC>
!  F90 is a strictly-typed language, and the syntax pass of the
!  compiler requires matching of type, kind and rank (TKR). Most calls
!  listed here use a generic type, shown here as <TT>MPP_TYPE_</TT>. This
!  is resolved in the pre-processor stage to any of a variety of
!  types. In general the MPP operations work on 4-byte and 8-byte
!  variants of <TT>integer, real, complex, logical</TT> variables, of
!  rank 0 to 5, leading to 48 specific module procedures under the same
!  generic interface. Any of the variables below shown as
!  <TT>MPP_TYPE_</TT> is treated in this way.
! </PUBLIC>

  use mpp_parameter_mod, only : MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE
  use mpp_parameter_mod, only : NOTE, WARNING, FATAL, MPP_CLOCK_DETAILED,MPP_CLOCK_SYNC
  use mpp_parameter_mod, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER
  use mpp_parameter_mod, only : CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA

  use mpp_data_mod,      only : request, mpp_record_timing_data

  use mpp_comm_mod,      only : mpp_chksum, mpp_max, mpp_min, mpp_sum, mpp_transmit
  use mpp_comm_mod,      only : mpp_send, mpp_recv, mpp_broadcast, mpp_malloc
  use mpp_comm_mod,      only : mpp_init, mpp_exit, mpp_set_stack_size
#ifdef use_MPI_GSM
  use mpp_comm_mod,      only : mpp_gsm_malloc, mpp_gsm_free
#endif

  use mpp_util_mod,      only : stdin, stdout, stderr, stdlog, lowercase, uppercase
  use mpp_util_mod,      only : mpp_error, mpp_error_state, mpp_set_warn_level, mpp_sync
  use mpp_util_mod,      only : mpp_sync_self, mpp_pe, mpp_node, mpp_npes, mpp_root_pe
  use mpp_util_mod,      only : mpp_set_root_pe, mpp_declare_pelist, mpp_get_current_pelist
  use mpp_util_mod,      only : mpp_set_current_pelist, mpp_clock_begin, mpp_clock_end
  use mpp_util_mod,      only : mpp_clock_id, mpp_clock_set_grain

  character(len=128), public :: version= &
       '$Id mpp.F90 $'
  character(len=128), public :: tagname= &
       '$Name: lima $'

  !--- public paramters  -----------------------------------------------
  public :: MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE, NOTE, WARNING, FATAL
  public :: MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
  public :: CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA

  !--- public data from mpp_data_mod ------------------------------
  public :: request

  !--- public interface from mpp_util_mod ------------------------------
  public :: stdin, stdout, stderr, stdlog, lowercase, uppercase, mpp_error, mpp_error_state
  public :: mpp_set_warn_level, mpp_sync, mpp_sync_self, mpp_set_stack_size, mpp_pe
  public :: mpp_node, mpp_npes, mpp_root_pe, mpp_set_root_pe, mpp_declare_pelist
  public :: mpp_get_current_pelist, mpp_set_current_pelist, mpp_clock_begin, mpp_clock_end
  public :: mpp_clock_id, mpp_clock_set_grain, mpp_record_timing_data

  !--- public interface from mpp_comm_mod ------------------------------
  public :: mpp_chksum, mpp_max, mpp_min, mpp_sum, mpp_transmit, mpp_send, mpp_recv
  public :: mpp_broadcast, mpp_malloc, mpp_init, mpp_exit
#ifdef use_MPI_GSM
  public :: mpp_gsm_malloc, mpp_gsm_free
#endif

  end module mpp_mod


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
  real                            :: dt

  call mpp_init()
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()
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
        !          call mpp_sync_self( (/modulo(pe+npes-i,npes)/) )
     end do
     call mpp_sync()
     call SYSTEM_CLOCK(tick)
     dt = real(tick-tick0)/(npes*ticks_per_sec)
     dt = max( dt, epsilon(dt) )
     if( pe.EQ.root )write( stdout(),'(/a,i8,f13.6,f8.2)' )'MPP_TRANSMIT length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
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
  if( pe.EQ.root )write( stdout(),'(a,2i4,f9.1,i8,f13.6,f8.2/)' ) &
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
  call flush(stdout(),istat)
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
  call test_shared_pointers(locd,n)

#ifdef use_MPI_GSM
  call mpp_gsm_free( locd )
#else
  if( pe.EQ.root )then
      deallocate( d )
  end if
#endif
#endif
  call mpp_exit()

contains

  subroutine test_shared_pointers(locd,n)
    integer(LONG_KIND), intent(in) :: locd
    integer :: n
    real :: dd(n)
    pointer( p, dd )

    p = locd
    print *, 'TEST_SHARED_POINTERS: pe, locd=', pe, locd
    print *, 'TEST_SHARED_POINTERS: pe, chksum(d)=', pe, mpp_chksum(dd,(/pe/))
    return
  end subroutine test_shared_pointers
end program test
  
#endif test_mpp

