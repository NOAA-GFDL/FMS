module mpp_comm_mod
#include <fms_platform.h>

#if defined(use_libSMA) && defined(sgi_mipspro)
  use shmem_interface
#endif

#if defined(use_libMPI) && defined(sgi_mipspro)
  use mpi
#endif

  use mpp_parameter_mod, only : MPP_CLOCK_SYNC, MPP_DEBUG, MPP_VERBOSE, NOTE, FATAL      
  use mpp_parameter_mod, only : MAX_EVENT_TYPES, EVENT_ALLREDUCE, NULL_PE, EVENT_SEND 
  use mpp_parameter_mod, only : ALL_PES, ANY_PE, NULL_PE, EVENT_RECV, EVENT_WAIT
  use mpp_parameter_mod, only : EVENT_BROADCAST, MPP_READY, MPP_WAIT
  use mpp_parameter_mod, only : mpp_parameter_version=>version, mpp_parameter_tagname=>tagname

  use mpp_data_mod,      only : peset, clocks, module_is_initialized=>mpp_is_initialized
  use mpp_data_mod,      only : error, pe, npes, root_pe, current_peset_num
  use mpp_data_mod,      only : world_peset_num, clock_num, current_clock, start_tick
  use mpp_data_mod,      only : tick,tick0, ticks_per_sec, max_ticks, tick_rate
  use mpp_data_mod,      only : debug=>debug_mpp, request, etc_unit, log_unit, configfile
  use mpp_data_mod,      only : stat, mpp_stack, ptr_stack, status, ptr_status, sync, ptr_sync
  use mpp_data_mod,      only : mpp_from_pe, ptr_from, remote_data_loc, ptr_remote, etcfile
  use mpp_data_mod,      only : first_call_system_clock_mpi, mpi_tick_rate, mpi_count0
  use mpp_data_mod,      only : mpp_data_version=>version, mpp_data_tagname=>tagname
  use mpp_data_mod,      only : mpp_comm_private

  use mpp_util_mod,      only : mpp_sync, mpp_error, mpp_npes, mpp_pe, stdlog, stdout, stderr
  use mpp_util_mod,      only : get_peset, get_unit, increment_current_clock, dump_clock_summary
  use mpp_util_mod,      only : sum_clock_data, mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_util_mod,      only : mpp_set_current_pelist, mpp_init_logfile
  use mpp_util_mod,      only : mpp_util_version=>version, mpp_util_tagname=>tagname

  implicit none
  private

#if defined(use_libSMA) || defined(use_GSM)
#include <mpp/shmem.fh>
#endif

#if defined(use_libMPI) && !defined(sgi_mipspro)
#include <mpif.h>   ! sgi_mipspro gets this from 'use mpi'
#endif

#ifdef use_libMPI
#ifdef _CRAYT3E
  !BWA: mpif.h on t3e currently does not contain MPI_INTEGER8 datatype
  !(O2k and t90 do)
  !(t3e: fixed on 3.3 I believe)
  integer, parameter :: MPI_INTEGER8=MPI_INTEGER
#endif
#endif use_libMPI

  integer            :: clock0    !measures total runtime from mpp_init to mpp_exit
  integer            :: mpp_stack_size=0, mpp_stack_hwm=0
  integer            :: tag=1
  logical            :: verbose=.FALSE.
#ifdef _CRAY
  integer(LONG_KIND) :: word(1)
#endif
#if defined(sgi_mipspro) || defined(__ia64)
  integer(INT_KIND)  :: word(1)
#endif

  character(len=128) :: version= &
       '$Id: mpp_comm.F90,v 12.0 2005/04/14 17:57:40 fms Exp $'
  character(len=128) :: tagname= &
       '$Name: lima $'

  public :: mpp_init, mpp_exit, mpp_min, mpp_max, mpp_sum, mpp_transmit, mpp_recv
  public :: mpp_send, mpp_broadcast, mpp_chksum, mpp_malloc, mpp_set_stack_size
#ifdef use_MPI_GSM
  public :: mpp_gsm_malloc, mpp_gsm_free
#endif

#ifdef use_libSMA
  !currently SMA contains no generic shmem_wait for different integer kinds:
  !I have inserted one here
  interface shmem_integer_wait
     module procedure shmem_int4_wait_local
     module procedure shmem_int8_wait_local
  end interface
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! <SUBROUTINE NAME="mpp_init">
  !  <OVERVIEW>
  !   Initialize <TT>mpp_mod</TT>.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Called to initialize the <TT>mpp_mod</TT> package. It is recommended
  !   that this call be the first executed line in your program. It sets the
  !   number of PEs assigned to this run (acquired from the command line, or
  !   through the environment variable <TT>NPES</TT>), and associates an ID
  !   number to each PE. These can be accessed by calling <LINK
  !   SRC="#mpp_npes"><TT>mpp_npes</TT></LINK> and <LINK
  !   SRC="#mpp_pe"><TT>mpp_pe</TT></LINK>.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call mpp_init( flags )
  !  </TEMPLATE>
  !  <IN NAME="flags" TYPE="integer">
  !   <TT>flags</TT> can be set to <TT>MPP_VERBOSE</TT> to
  !   have <TT>mpp_mod</TT> keep you informed of what it's up to.
  !  </IN>
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="mpp_exit">
  !  <OVERVIEW>
  !   Exit <TT>mpp_mod</TT>.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Called at the end of the run, or to re-initialize <TT>mpp_mod</TT>,
  !   should you require that for some odd reason.
  !
  !   This call implies synchronization across all PEs.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call mpp_exit()
  !  </TEMPLATE>
  ! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="mpp_malloc">
  !  <OVERVIEW>
  !    Symmetric memory allocation.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This routine is used on SGI systems when <TT>mpp_mod</TT> is
  !    invoked in the SHMEM library. It ensures that dynamically allocated
  !    memory can be used with <TT>shmem_get</TT> and
  !    <TT>shmem_put</TT>. This is called <I>symmetric
  !    allocation</I> and is described in the
  !    <TT>intro_shmem</TT> man page. <TT>ptr</TT> is a <I>Cray
  !    pointer</I> (see the section on <LINK
  !    SRC="#PORTABILITY">portability</LINK>).  The operation can be expensive
  !    (since it requires a global barrier). We therefore attempt to re-use
  !    existing allocation whenever possible. Therefore <TT>len</TT>
  !    and <TT>ptr</TT> must have the <TT>SAVE</TT> attribute
  !    in the calling routine, and retain the information about the last call
  !    to <TT>mpp_malloc</TT>. Additional memory is symmetrically
  !    allocated if and only if <TT>newlen</TT> exceeds
  !    <TT>len</TT>.
  !
  !    This is never required on Cray PVP or MPP systems. While the T3E
  !    manpages do talk about symmetric allocation, <TT>mpp_mod</TT>
  !    is coded to remove this restriction.
  !
  !    It is never required if <TT>mpp_mod</TT> is invoked in MPI.
  !
  !   This call implies synchronization across all PEs.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call mpp_malloc( ptr, newlen, len )
  !  </TEMPLATE>
  !  <IN NAME="ptr">
  !     a cray pointer, points to a dummy argument in this routine.
  !  </IN>
  !  <IN NAME="newlen" TYPE="integer">
  !     the required allocation length for the pointer ptr
  !  </IN>
  !  <IN NAME="len" TYPE="integer">
  !     the current allocation (0 if unallocated).
  !  </IN>
  ! </SUBROUTINE>

  !#####################################################################

  ! <SUBROUTINE NAME="mpp_set_stack_size">
  !  <OVERVIEW>
  !    Allocate module internal workspace.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    <TT>mpp_mod</TT> maintains a private internal array called
  !    <TT>mpp_stack</TT> for private workspace. This call sets the length,
  !    in words, of this array. 
  !
  !    The <TT>mpp_init</TT> call sets this
  !    workspace length to a default of 32768, and this call may be used if a
  !    longer workspace is needed.
  !    
  !    This call implies synchronization across all PEs.
  !    
  !    This workspace is symmetrically allocated, as required for
  !    efficient communication on SGI and Cray MPP systems. Since symmetric
  !    allocation must be performed by <I>all</I> PEs in a job, this call
  !    must also be called by all PEs, using the same value of
  !    <TT>n</TT>. Calling <TT>mpp_set_stack_size</TT> from a subset of PEs,
  !    or with unequal argument <TT>n</TT>, may cause the program to hang.
  !    
  !    If any MPP call using <TT>mpp_stack</TT> overflows the declared
  !    stack array, the program will abort with a message specifying the
  !    stack length that is required. Many users wonder why, if the required
  !    stack length can be computed, it cannot also be specified at that
  !    point. This cannot be automated because there is no way for the
  !    program to know if all PEs are present at that call, and with equal
  !    values of <TT>n</TT>. The program must be rerun by the user with the
  !    correct argument to <TT>mpp_set_stack_size</TT>, called at an
  !    appropriate point in the code where all PEs are known to be present.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_set_stack_size(n)
  !  </TEMPLATE>
  !  <IN NAME="n" TYPE="integer"></IN>
  ! </SUBROUTINE>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !            GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min             !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! <INTERFACE NAME="mpp_max">
  !  <OVERVIEW>
  !    Reduction operations.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Find the max of scalar a the PEs in pelist
  !    result is also automatically broadcast to all PEs
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call  mpp_max( a, pelist )
  !  </TEMPLATE>
  !  <IN NAME="a">
  !    <TT>real</TT> or <TT>integer</TT>, of 4-byte of 8-byte kind.
  !  </IN>
  !  <IN NAME="pelist">
  !    If <TT>pelist</TT> is omitted, the context is assumed to be the
  !    current pelist. This call implies synchronization across the PEs in
  !    <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !  </IN>
  ! </INTERFACE>

  interface mpp_max
     module procedure mpp_max_real8
#ifndef no_8byte_integers
     module procedure mpp_max_int8
#endif
#ifndef no_4byte_reals
     module procedure mpp_max_real4
#endif
     module procedure mpp_max_int4
  end interface

  interface mpp_min
     module procedure mpp_min_real8
#ifndef no_8byte_integers
     module procedure mpp_min_int8
#endif
#ifndef no_4byte_reals
     module procedure mpp_min_real4
#endif
     module procedure mpp_min_int4
  end interface


  ! <INTERFACE NAME="mpp_sum">
  !  <OVERVIEW>
  !    Reduction operation.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !    <TT>integer, real, complex</TT> variables, of rank 0 or 1. A
  !    contiguous block from a multi-dimensional array may be passed by its
  !    starting address and its length, as in <TT>f77</TT>.
  !
  !    Library reduction operators are not required or guaranteed to be
  !    bit-reproducible. In any case, changing the processor count changes
  !    the data layout, and thus very likely the order of operations. For
  !    bit-reproducible sums of distributed arrays, consider using the
  !    <TT>mpp_global_sum</TT> routine provided by the <LINK
  !    SRC="mpp_domains.html"><TT>mpp_domains</TT></LINK> module.
  !
  !    The <TT>bit_reproducible</TT> flag provided in earlier versions of
  !    this routine has been removed.
  !
  !
  !    If <TT>pelist</TT> is omitted, the context is assumed to be the
  !    current pelist. This call implies synchronization across the PEs in
  !    <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_sum( a, length, pelist )
  !  </TEMPLATE>
  !  <IN NAME="length"></IN>
  !  <IN NAME="pelist"></IN>
  !  <INOUT NAME="a"></INOUT>
  ! </INTERFACE>

  interface mpp_sum
#ifndef no_8byte_integers
     module procedure mpp_sum_int8
     module procedure mpp_sum_int8_scalar
     module procedure mpp_sum_int8_2d
     module procedure mpp_sum_int8_3d
     module procedure mpp_sum_int8_4d
     module procedure mpp_sum_int8_5d
#endif
     module procedure mpp_sum_real8
     module procedure mpp_sum_real8_scalar
     module procedure mpp_sum_real8_2d
     module procedure mpp_sum_real8_3d
     module procedure mpp_sum_real8_4d
     module procedure mpp_sum_real8_5d
     module procedure mpp_sum_cmplx8
     module procedure mpp_sum_cmplx8_scalar
     module procedure mpp_sum_cmplx8_2d
     module procedure mpp_sum_cmplx8_3d
     module procedure mpp_sum_cmplx8_4d
     module procedure mpp_sum_cmplx8_5d
     module procedure mpp_sum_int4
     module procedure mpp_sum_int4_scalar
     module procedure mpp_sum_int4_2d
     module procedure mpp_sum_int4_3d
     module procedure mpp_sum_int4_4d
     module procedure mpp_sum_int4_5d
#ifndef no_4byte_reals
     module procedure mpp_sum_real4
     module procedure mpp_sum_real4_scalar
     module procedure mpp_sum_real4_2d
     module procedure mpp_sum_real4_3d
     module procedure mpp_sum_real4_4d
     module procedure mpp_sum_real4_5d
     module procedure mpp_sum_cmplx4
     module procedure mpp_sum_cmplx4_scalar
     module procedure mpp_sum_cmplx4_2d
     module procedure mpp_sum_cmplx4_3d
     module procedure mpp_sum_cmplx4_4d
     module procedure mpp_sum_cmplx4_5d
#endif
  end interface

  !#####################################################################

  ! <INTERFACE NAME="mpp_transmit">
  !  <OVERVIEW>
  !    Basic message-passing call.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !    <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
  !    contiguous block from a multi-dimensional array may be passed by its
  !    starting address and its length, as in <TT>f77</TT>.
  !    
  !    <TT>mpp_transmit</TT> is currently implemented as asynchronous
  !    outward transmission and synchronous inward transmission. This follows
  !    the behaviour of <TT>shmem_put</TT> and <TT>shmem_get</TT>. In MPI, it
  !    is implemented as <TT>mpi_isend</TT> and <TT>mpi_recv</TT>. For most
  !    applications, transmissions occur in pairs, and are here accomplished
  !    in a single call.
  !    
  !    The special PE designations <TT>NULL_PE</TT>,
  !    <TT>ANY_PE</TT> and <TT>ALL_PES</TT> are provided by use
  !    association.
  !    
  !    <TT>NULL_PE</TT>: is used to disable one of the pair of
  !    transmissions.<BR/>
  !    <TT>ANY_PE</TT>: is used for unspecific remote
  !    destination. (Please note that <TT>put_pe=ANY_PE</TT> has no meaning
  !    in the MPI context, though it is available in the SHMEM invocation. If
  !    portability is a concern, it is best avoided).<BR/>
  !    <TT>ALL_PES</TT>: is used for broadcast operations.
  !    
  !    It is recommended that <LINK
  !    SRC="#mpp_broadcast"><TT>mpp_broadcast</TT></LINK> be used for
  !    broadcasts.
  !    
  !    The following example illustrates the use of
  !    <TT>NULL_PE</TT> and <TT>ALL_PES</TT>:
  !    
  !    <PRE>
  !    real, dimension(n) :: a
  !    if( pe.EQ.0 )then
  !        do p = 1,npes-1
  !           call mpp_transmit( a, n, p, a, n, NULL_PE )
  !        end do
  !    else
  !        call mpp_transmit( a, n, NULL_PE, a, n, 0 )
  !    end if
  !    
  !    call mpp_transmit( a, n, ALL_PES, a, n, 0 )
  !    </PRE>
  !    
  !    The do loop and the broadcast operation above are equivalent.
  !    
  !    Two overloaded calls <TT>mpp_send</TT> and
  !     <TT>mpp_recv</TT> have also been
  !    provided. <TT>mpp_send</TT> calls <TT>mpp_transmit</TT>
  !    with <TT>get_pe=NULL_PE</TT>. <TT>mpp_recv</TT> calls
  !    <TT>mpp_transmit</TT> with <TT>put_pe=NULL_PE</TT>. Thus
  !    the do loop above could be written more succinctly:
  !    
  !    <PRE>
  !    if( pe.EQ.0 )then
  !        do p = 1,npes-1
  !           call mpp_send( a, n, p )
  !        end do
  !    else
  !        call mpp_recv( a, n, 0 )
  !    end if
  !    </PRE>
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_transmit( put_data, put_len, put_pe, get_data, get_len, get_pe )
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_transmit
     module procedure mpp_transmit_real8
     module procedure mpp_transmit_real8_scalar
     module procedure mpp_transmit_real8_2d
     module procedure mpp_transmit_real8_3d
     module procedure mpp_transmit_real8_4d
     module procedure mpp_transmit_real8_5d
     module procedure mpp_transmit_cmplx8
     module procedure mpp_transmit_cmplx8_scalar
     module procedure mpp_transmit_cmplx8_2d
     module procedure mpp_transmit_cmplx8_3d
     module procedure mpp_transmit_cmplx8_4d
     module procedure mpp_transmit_cmplx8_5d
#ifndef no_8byte_integers
     module procedure mpp_transmit_int8
     module procedure mpp_transmit_int8_scalar
     module procedure mpp_transmit_int8_2d
     module procedure mpp_transmit_int8_3d
     module procedure mpp_transmit_int8_4d
     module procedure mpp_transmit_int8_5d
     module procedure mpp_transmit_logical8
     module procedure mpp_transmit_logical8_scalar
     module procedure mpp_transmit_logical8_2d
     module procedure mpp_transmit_logical8_3d
     module procedure mpp_transmit_logical8_4d
     module procedure mpp_transmit_logical8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_transmit_real4
     module procedure mpp_transmit_real4_scalar
     module procedure mpp_transmit_real4_2d
     module procedure mpp_transmit_real4_3d
     module procedure mpp_transmit_real4_4d
     module procedure mpp_transmit_real4_5d
     module procedure mpp_transmit_cmplx4
     module procedure mpp_transmit_cmplx4_scalar
     module procedure mpp_transmit_cmplx4_2d
     module procedure mpp_transmit_cmplx4_3d
     module procedure mpp_transmit_cmplx4_4d
     module procedure mpp_transmit_cmplx4_5d
#endif
     module procedure mpp_transmit_int4
     module procedure mpp_transmit_int4_scalar
     module procedure mpp_transmit_int4_2d
     module procedure mpp_transmit_int4_3d
     module procedure mpp_transmit_int4_4d
     module procedure mpp_transmit_int4_5d
     module procedure mpp_transmit_logical4
     module procedure mpp_transmit_logical4_scalar
     module procedure mpp_transmit_logical4_2d
     module procedure mpp_transmit_logical4_3d
     module procedure mpp_transmit_logical4_4d
     module procedure mpp_transmit_logical4_5d
  end interface
  interface mpp_recv
     module procedure mpp_recv_real8
     module procedure mpp_recv_real8_scalar
     module procedure mpp_recv_real8_2d
     module procedure mpp_recv_real8_3d
     module procedure mpp_recv_real8_4d
     module procedure mpp_recv_real8_5d
     module procedure mpp_recv_cmplx8
     module procedure mpp_recv_cmplx8_scalar
     module procedure mpp_recv_cmplx8_2d
     module procedure mpp_recv_cmplx8_3d
     module procedure mpp_recv_cmplx8_4d
     module procedure mpp_recv_cmplx8_5d
#ifndef no_8byte_integers
     module procedure mpp_recv_int8
     module procedure mpp_recv_int8_scalar
     module procedure mpp_recv_int8_2d
     module procedure mpp_recv_int8_3d
     module procedure mpp_recv_int8_4d
     module procedure mpp_recv_int8_5d
     module procedure mpp_recv_logical8
     module procedure mpp_recv_logical8_scalar
     module procedure mpp_recv_logical8_2d
     module procedure mpp_recv_logical8_3d
     module procedure mpp_recv_logical8_4d
     module procedure mpp_recv_logical8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_recv_real4
     module procedure mpp_recv_real4_scalar
     module procedure mpp_recv_real4_2d
     module procedure mpp_recv_real4_3d
     module procedure mpp_recv_real4_4d
     module procedure mpp_recv_real4_5d
     module procedure mpp_recv_cmplx4
     module procedure mpp_recv_cmplx4_scalar
     module procedure mpp_recv_cmplx4_2d
     module procedure mpp_recv_cmplx4_3d
     module procedure mpp_recv_cmplx4_4d
     module procedure mpp_recv_cmplx4_5d
#endif
     module procedure mpp_recv_int4
     module procedure mpp_recv_int4_scalar
     module procedure mpp_recv_int4_2d
     module procedure mpp_recv_int4_3d
     module procedure mpp_recv_int4_4d
     module procedure mpp_recv_int4_5d
     module procedure mpp_recv_logical4
     module procedure mpp_recv_logical4_scalar
     module procedure mpp_recv_logical4_2d
     module procedure mpp_recv_logical4_3d
     module procedure mpp_recv_logical4_4d
     module procedure mpp_recv_logical4_5d
  end interface
  interface mpp_send
     module procedure mpp_send_real8
     module procedure mpp_send_real8_scalar
     module procedure mpp_send_real8_2d
     module procedure mpp_send_real8_3d
     module procedure mpp_send_real8_4d
     module procedure mpp_send_real8_5d
     module procedure mpp_send_cmplx8
     module procedure mpp_send_cmplx8_scalar
     module procedure mpp_send_cmplx8_2d
     module procedure mpp_send_cmplx8_3d
     module procedure mpp_send_cmplx8_4d
     module procedure mpp_send_cmplx8_5d
#ifndef no_8byte_integers
     module procedure mpp_send_int8
     module procedure mpp_send_int8_scalar
     module procedure mpp_send_int8_2d
     module procedure mpp_send_int8_3d
     module procedure mpp_send_int8_4d
     module procedure mpp_send_int8_5d
     module procedure mpp_send_logical8
     module procedure mpp_send_logical8_scalar
     module procedure mpp_send_logical8_2d
     module procedure mpp_send_logical8_3d
     module procedure mpp_send_logical8_4d
     module procedure mpp_send_logical8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_send_real4
     module procedure mpp_send_real4_scalar
     module procedure mpp_send_real4_2d
     module procedure mpp_send_real4_3d
     module procedure mpp_send_real4_4d
     module procedure mpp_send_real4_5d
     module procedure mpp_send_cmplx4
     module procedure mpp_send_cmplx4_scalar
     module procedure mpp_send_cmplx4_2d
     module procedure mpp_send_cmplx4_3d
     module procedure mpp_send_cmplx4_4d
     module procedure mpp_send_cmplx4_5d
#endif
     module procedure mpp_send_int4
     module procedure mpp_send_int4_scalar
     module procedure mpp_send_int4_2d
     module procedure mpp_send_int4_3d
     module procedure mpp_send_int4_4d
     module procedure mpp_send_int4_5d
     module procedure mpp_send_logical4
     module procedure mpp_send_logical4_scalar
     module procedure mpp_send_logical4_2d
     module procedure mpp_send_logical4_3d
     module procedure mpp_send_logical4_4d
     module procedure mpp_send_logical4_5d
  end interface

  ! <INTERFACE NAME="mpp_broadcast">

  !   <OVERVIEW>
  !     Parallel broadcasts.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     The <TT>mpp_broadcast</TT> call has been added because the original
  !     syntax (using <TT>ALL_PES</TT> in <TT>mpp_transmit</TT>) did not
  !     support a broadcast across a pelist.
  !
  !     <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !     <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
  !     contiguous block from a multi-dimensional array may be passed by its
  !     starting address and its length, as in <TT>f77</TT>.
  !
  !     Global broadcasts through the <TT>ALL_PES</TT> argument to <LINK
  !     SRC="#mpp_transmit"><TT>mpp_transmit</TT></LINK> are still provided for
  !     backward-compatibility.
  !
  !     If <TT>pelist</TT> is omitted, the context is assumed to be the
  !     current pelist. <TT>from_pe</TT> must belong to the current
  !     pelist. This call implies synchronization across the PEs in
  !     <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call mpp_broadcast( data, length, from_pe, pelist )
  !   </TEMPLATE>
  !   <IN NAME="length"> </IN>
  !   <IN NAME="from_pe"> </IN>
  !   <IN NAME="pelist"> </IN>
  !   <INOUT NAME="data(*)"> </INOUT>
  ! </INTERFACE>
  interface mpp_broadcast
     module procedure mpp_broadcast_real8
     module procedure mpp_broadcast_real8_scalar
     module procedure mpp_broadcast_real8_2d
     module procedure mpp_broadcast_real8_3d
     module procedure mpp_broadcast_real8_4d
     module procedure mpp_broadcast_real8_5d
     module procedure mpp_broadcast_cmplx8
     module procedure mpp_broadcast_cmplx8_scalar
     module procedure mpp_broadcast_cmplx8_2d
     module procedure mpp_broadcast_cmplx8_3d
     module procedure mpp_broadcast_cmplx8_4d
     module procedure mpp_broadcast_cmplx8_5d
#ifndef no_8byte_integers
     module procedure mpp_broadcast_int8
     module procedure mpp_broadcast_int8_scalar
     module procedure mpp_broadcast_int8_2d
     module procedure mpp_broadcast_int8_3d
     module procedure mpp_broadcast_int8_4d
     module procedure mpp_broadcast_int8_5d
     module procedure mpp_broadcast_logical8
     module procedure mpp_broadcast_logical8_scalar
     module procedure mpp_broadcast_logical8_2d
     module procedure mpp_broadcast_logical8_3d
     module procedure mpp_broadcast_logical8_4d
     module procedure mpp_broadcast_logical8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_broadcast_real4
     module procedure mpp_broadcast_real4_scalar
     module procedure mpp_broadcast_real4_2d
     module procedure mpp_broadcast_real4_3d
     module procedure mpp_broadcast_real4_4d
     module procedure mpp_broadcast_real4_5d
     module procedure mpp_broadcast_cmplx4
     module procedure mpp_broadcast_cmplx4_scalar
     module procedure mpp_broadcast_cmplx4_2d
     module procedure mpp_broadcast_cmplx4_3d
     module procedure mpp_broadcast_cmplx4_4d
     module procedure mpp_broadcast_cmplx4_5d
#endif
     module procedure mpp_broadcast_int4
     module procedure mpp_broadcast_int4_scalar
     module procedure mpp_broadcast_int4_2d
     module procedure mpp_broadcast_int4_3d
     module procedure mpp_broadcast_int4_4d
     module procedure mpp_broadcast_int4_5d
     module procedure mpp_broadcast_logical4
     module procedure mpp_broadcast_logical4_scalar
     module procedure mpp_broadcast_logical4_2d
     module procedure mpp_broadcast_logical4_3d
     module procedure mpp_broadcast_logical4_4d
     module procedure mpp_broadcast_logical4_5d
  end interface

  !#####################################################################
  ! <INTERFACE NAME="mpp_chksum">

  !   <OVERVIEW>
  !     Parallel checksums.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     <TT>mpp_chksum</TT> is a parallel checksum routine that returns an
  !     identical answer for the same array irrespective of how it has been
  !     partitioned across processors. <TT>LONG_KIND</TT>is the <TT>KIND</TT>
  !     parameter corresponding to long integers (see discussion on
  !     OS-dependent preprocessor directives) defined in
  !     the header file <TT>fms_platform.h</TT>. <TT>MPP_TYPE_</TT> corresponds to any
  !     4-byte and 8-byte variant of <TT>integer, real, complex, logical</TT>
  !     variables, of rank 0 to 5.
  !
  !     Integer checksums on FP data use the F90 <TT>TRANSFER()</TT>
  !     intrinsic.
  !
  !     The <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/chksum/chksum.html">serial checksum module</LINK> is superseded
  !     by this function, and is no longer being actively maintained. This
  !     provides identical results on a single-processor job, and to perform
  !     serial checksums on a single processor of a parallel job, you only
  !     need to use the optional <TT>pelist</TT> argument.
  !     <PRE>
  !     use mpp_mod
  !     integer :: pe, chksum
  !     real :: a(:)
  !     pe = mpp_pe()
  !     chksum = mpp_chksum( a, (/pe/) )
  !     </PRE>
  !
  !     The additional functionality of <TT>mpp_chksum</TT> over
  !     serial checksums is to compute the checksum across the PEs in
  !     <TT>pelist</TT>. The answer is guaranteed to be the same for
  !     the same distributed array irrespective of how it has been
  !     partitioned.
  !
  !     If <TT>pelist</TT> is omitted, the context is assumed to be the
  !     current pelist. This call implies synchronization across the PEs in
  !     <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     mpp_chksum( var, pelist )
  !   </TEMPLATE>
  !   <IN NAME="pelist" TYPE="integer" DIM="(:)"> </IN>
  !   <IN NAME="var" TYPE="MPP_TYPE_"> </IN>
  ! </INTERFACE>
  interface mpp_chksum
#ifndef no_8byte_integers
     module procedure mpp_chksum_i8_1d
     module procedure mpp_chksum_i8_2d
     module procedure mpp_chksum_i8_3d
     module procedure mpp_chksum_i8_4d
#endif
     module procedure mpp_chksum_i4_1d
     module procedure mpp_chksum_i4_2d
     module procedure mpp_chksum_i4_3d
     module procedure mpp_chksum_i4_4d
     module procedure mpp_chksum_r8_0d
     module procedure mpp_chksum_r8_1d
     module procedure mpp_chksum_r8_2d
     module procedure mpp_chksum_r8_3d
     module procedure mpp_chksum_r8_4d
     module procedure mpp_chksum_r8_5d
     module procedure mpp_chksum_c8_0d
     module procedure mpp_chksum_c8_1d
     module procedure mpp_chksum_c8_2d
     module procedure mpp_chksum_c8_3d
     module procedure mpp_chksum_c8_4d
     module procedure mpp_chksum_c8_5d
#ifndef no_4byte_reals
     module procedure mpp_chksum_r4_0d
     module procedure mpp_chksum_r4_1d
     module procedure mpp_chksum_r4_2d
     module procedure mpp_chksum_r4_3d
     module procedure mpp_chksum_r4_4d
     module procedure mpp_chksum_r4_5d
     module procedure mpp_chksum_c4_0d
     module procedure mpp_chksum_c4_1d
     module procedure mpp_chksum_c4_2d
     module procedure mpp_chksum_c4_3d
     module procedure mpp_chksum_c4_4d
     module procedure mpp_chksum_c4_5d
#endif
  end interface

contains

#include <system_clock.h>

#ifdef use_libSMA
#include <mpp_comm_sma.inc>
#elif defined(use_libMPI)
#include <mpp_comm_mpi.inc>
#else
#include <mpp_comm_nocomm.inc>
#endif

#ifndef no_8byte_integers
#define MPP_CHKSUM_INT_ mpp_chksum_i8_1d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i8_2d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i8_3d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i8_4d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:,:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:,:,:,:,:)
#include <mpp_chksum_int.h>
#endif

#define MPP_CHKSUM_INT_ mpp_chksum_i4_1d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i4_2d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i4_3d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i4_4d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:,:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:,:,:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_ mpp_chksum_r8_0d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ !
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_1d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_2d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_3d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_4d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_0d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ !
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_1d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_2d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_3d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_4d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_5d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>

#ifndef no_4byte_reals
  !CAUTION: the r4 versions of these may produce
  !unpredictable results: I'm not sure what the result
  !of the TRANSFER() to integer(8) is from an odd number of real(4)s?
  !However the complex(4) will work, since it is guaranteed even.
#define MPP_CHKSUM_ mpp_chksum_r4_0d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ !
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_1d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_2d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_3d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_4d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_0d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ !
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_1d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_2d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_3d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_4d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_5d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>
#endif

end module mpp_comm_mod
