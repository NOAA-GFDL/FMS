!-----------------------------------------------------------------------
!                 Communication for message-passing codes
!
! AUTHOR: V. Balaji (vb@gfdl.gov)
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

!these are used to determine hardware/OS/compiler

#ifdef __sgi
#ifdef _COMPILER_VERSION
!the MIPSPro compiler defines _COMPILER_VERSION
#define sgi_mipspro
#else
#define sgi_generic
#endif
#endif

#if defined(_CRAY) || defined(sgi_mipspro)
#define SGICRAY
#endif

!compilers that support Cray pointers
#if defined(SGICRAY) || defined(__alpha)
#define use_CRI_pointers
#endif

!values of kind: double and long are 8-byte, float and int are 4-byte
#if defined(SGICRAY)
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#else
!these might be different on non-SGICRAY, I believe
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#endif

#ifdef sgi_generic
!this is for the Edinburgh n32/o32 compiler, which won't accept 8-byte ints
!at any price
#define LONG_KIND 4
#endif

!parallel machine types
#if defined(_CRAY) && !defined(_CRAYT3E) && !defined(_CRAYT3D)
#define CRAYPVP
#endif

#if defined(_CRAYT3E) || defined(_CRAYT3D) || defined(sgi_mipspro)
#define SGICRAY_MPP
#endif

!only one of SMA or MPI can be used (though mixing calls is allowed, this module will not)
#ifdef use_libSMA
#undef use_libMPI
#endif

!if using shmem calls on Origin, you may need to use shmalloc
#ifdef use_libSMA
#ifdef sgi_mipspro
#define use_shmalloc
#endif
#endif

!sums can be computed using library reduce calls (SHMEM_SUM, MPI_REDUCE) which do not reproduce the same bits. I have written a
! binary tree summation routine that reproduces bits and does sums over p processors in log(p) time. Specify -Dbit_reproducible to
! request the bit-reproducible version.
#ifndef bit_reproducible
#define use_library_summation_routines
#endif

!dimension of work array used by certain shmem collective operations, see mpp_max, mpp_sum
#ifndef WORK_ARRAY_SIZE
#define WORK_ARRAY_SIZE 100000
#endif

!dimension of arrays requiring at least NPES elements
#ifndef MAXPES
#define MAXPES 2048
#endif

!various lengths (see shpalloc) are estimated in "words" which are 32bit on SGI, 64bit on Cray
#ifdef _CRAY
#define WORDS_PER_REAL 1
#endif
#ifdef sgi_mipspro
#define WORDS_PER_REAL 2
#endif

module mpp_mod
!string BWA is used to tag lines that are bug workarounds and will disappear when offending compiler bug is fixed
!a generalized communication and domain decomposition package for use with shmem and MPI
!will add: co_array_fortran, MPI2
!Balaji (vb@gfdl.gov) 11 May 1998
#if defined(sgi_mipspro) && defined(use_libSMA)
  use shmem_interface
#endif
  implicit none
  private
  character(len=256), private :: version='$Id: mpp.F90,v 5.6 2000/03/01 15:57:19 vb Exp $'

#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
  integer, parameter, private :: stdin=100, stdout=101, stderr=102
#else
  integer, parameter, private :: stdin=5, stdout=6, stderr=0
#endif
  logical, private :: mpp_initialized=.FALSE.
  integer, private :: pe=0, npes=1
  integer, private :: error

!initialization flags
  integer, parameter, public :: MPP_VERBOSE=1
  logical, private :: verbose=.FALSE.

!flags to transmit routines
  integer, parameter, public :: ALL_PES=-1, ANY_PE=-2, NULL_PE=-3

!errortype flags
  integer, parameter, public :: NOTE=0, WARNING=1, FATAL=2
  logical, private :: warnings_are_fatal = .FALSE.

!timing
!since these are mainly timers associated with communication, we only do real times, not CPU times
  integer, private :: tick, ticks_per_sec

  integer(LONG_KIND), parameter, private :: MPP_WAIT_VALUE=-1
  integer, parameter, private :: MPP_PUT_WAIT=-1
  integer, private :: put_is_done !-1 if put is not done; get_pe returns its ID as tag after a shmem_get
#ifdef use_libSMA
#include <mpp/shmem.fh>
  integer :: rSync(SHMEM_REDUCE_SYNC_SIZE)
  integer :: bSync(SHMEM_BCAST_SYNC_SIZE)
  integer :: aSync(SHMEM_BARRIER_SYNC_SIZE)
  integer(LONG_KIND), private :: mpp_ready_to_recv(0:MAXPES-1)   !used to negotiate communication
  real(DOUBLE_KIND), private :: symm_work_array(WORK_ARRAY_SIZE) !symm_work_array may be used by SHMEM collective ops
  integer, private :: mpp_get_pe !used to announce from where data is coming from (only used by get_pe=ANY_PE)
  integer(LONG_KIND), private :: remote_data_ptr(0:MAXPES-1) !used to send pointer to remote data
#ifdef SGICRAY_MPP
!we call shpalloc in mpp_init() to ensure all these are remotely accessible
!on PVP where shpalloc doesn't exist, module variables are automatically guaranteed to be remotely accessible
  integer(INT_KIND) :: dum4     !BWA: this jogs all ptr addresses by 4-bytes and forces 8-byte alignment (see PV 696620)
  pointer( ptr$rSync, rSync )
  pointer( ptr$bSync, bSync )
  pointer( ptr$aSync, aSync )
  pointer( ptr$ready_to_recv, mpp_ready_to_recv )
  pointer( ptr$get_pe, mpp_get_pe )
  pointer( ptr$put_is_done, put_is_done )
  pointer( ptr$remote_data_ptr, remote_data_ptr )
  pointer( ptr$symm, symm_work_array )
#endif
#endif use_libSMA
#ifdef use_libMPI
#include <mpif.h>
!tag is never used, but must be initialized to non-negative integer
  integer, private :: tag=1, request(0:MAXPES-1)=MPI_REQUEST_NULL, stat(MPI_STATUS_SIZE), group_all
#ifdef _CRAYT3E
!BWA: mpif.h on t3e currently does not contain MPI_INTEGER8 datatype (O2k and t90 do)
!(t3e: fixed on 3.2 I believe)
  integer, parameter :: MPI_INTEGER8=MPI_INTEGER
#endif
#endif use_libMPI

!public interfaces
  interface mpp_max
     module procedure mpp_max_real8
     module procedure mpp_max_int8
  end interface
  interface mpp_min
     module procedure mpp_min_real8
     module procedure mpp_min_int8
  end interface
  interface mpp_sum
     module procedure mpp_sum_int8
     module procedure mpp_sum_real8
     module procedure mpp_sum_cmplx8
     module procedure mpp_sum_real8_scalar
     module procedure mpp_sum_cmplx8_scalar
  end interface
  interface mpp_transmit
     module procedure mpp_transmit_real8
     module procedure mpp_transmit_real8_scalar
     module procedure mpp_transmit_cmplx8
     module procedure mpp_transmit_cmplx8_scalar
     module procedure mpp_transmit_int8
     module procedure mpp_transmit_int8_scalar
  end interface
  interface mpp_recv
     module procedure mpp_recv_real8
     module procedure mpp_recv_real8_scalar
     module procedure mpp_recv_cmplx8
     module procedure mpp_recv_cmplx8_scalar
     module procedure mpp_recv_int8
     module procedure mpp_recv_int8_scalar
  end interface
  interface mpp_send
     module procedure mpp_send_real8
     module procedure mpp_send_real8_scalar
     module procedure mpp_send_cmplx8
     module procedure mpp_send_cmplx8_scalar
     module procedure mpp_send_int8
     module procedure mpp_send_int8_scalar
  end interface

  interface mpp_chksum
     module procedure mpp_chksum_int_1d
     module procedure mpp_chksum_int_2d
     module procedure mpp_chksum_int_3d
     module procedure mpp_chksum_int_4d
#ifdef use_CRI_pointers
     module procedure mpp_chksum_r8_0d
     module procedure mpp_chksum_r8_1d
     module procedure mpp_chksum_r8_2d
     module procedure mpp_chksum_r8_3d
     module procedure mpp_chksum_r8_4d
     module procedure mpp_chksum_c8_0d
     module procedure mpp_chksum_c8_1d
     module procedure mpp_chksum_c8_2d
     module procedure mpp_chksum_c8_3d
     module procedure mpp_chksum_c8_4d
#endif
  end interface

#ifdef use_libSMA
!currently SMA contains no generic shmem_wait for different integer kinds: I have inserted one here
  interface shmem_integer_wait
     module procedure shmem_int4_wait_local
     module procedure shmem_int8_wait_local
  end interface
#endif
  public :: mpp_chksum, mpp_error, mpp_exit, mpp_init, mpp_max, mpp_min, mpp_npes, mpp_pe, mpp_recv, mpp_send, &
       mpp_set_warn_level, mpp_sync, mpp_sync_self, mpp_sum, mpp_transmit
#ifdef use_CRI_pointers
  public :: mpp_malloc
#endif

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit         !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_init(flags)
      integer, optional, intent(in) :: flags
      integer :: my_pe, num_pes
#ifdef _CRAYT3E
      intrinsic my_pe
#endif

      if( mpp_initialized )return

#ifdef use_libSMA
      call START_PES(0)         !the argument 0 means extract from environment variable NPES on PVP/SGI, from mpprun -n on t3e
      pe = my_pe()
      npes = num_pes()
      if( npes.GT.MAXPES .AND. pe.EQ.0 )then
          print *, 'You have requested a processor count npes=', npes, ' larger than MAXPES=', MAXPES
          print *, 'Please recompile with the cpp flag -DMAXPES=',npes
          call abort()
      end if
#ifdef SGICRAY_MPP
!we use shpalloc to ensure all these are remotely accessible
      call SHMEM_BARRIER_ALL()
      call SHPALLOC( ptr$rsync, SHMEM_REDUCE_SYNC_SIZE,  error, -1 )
      call SHMEM_BARRIER_ALL()
      call SHPALLOC( ptr$bSync, SHMEM_BCAST_SYNC_SIZE,   error, -1 )
      call SHMEM_BARRIER_ALL()
      call SHPALLOC( ptr$aSync, SHMEM_BARRIER_SYNC_SIZE, error, -1 )
      call SHMEM_BARRIER_ALL()
      call SHPALLOC( ptr$symm, WORK_ARRAY_SIZE*WORDS_PER_REAL, error, -1 )
      call SHMEM_BARRIER_ALL()
      call SHPALLOC( ptr$ready_to_recv, MAXPES*WORDS_PER_REAL, error, -1 )
      call SHMEM_BARRIER_ALL()
      call SHPALLOC( ptr$remote_data_ptr, MAXPES*WORDS_PER_REAL, error, -1 )
      call SHMEM_BARRIER_ALL()
      call SHPALLOC( ptr$get_pe,1, error, -1 )
      call SHMEM_BARRIER_ALL()
      call SHPALLOC( ptr$put_is_done, 1, error, -1 )
      call SHMEM_BARRIER_ALL()
#endif
      rSync = SHMEM_SYNC_VALUE
      bSync = SHMEM_SYNC_VALUE
      aSync = SHMEM_SYNC_VALUE
#endif use_libSMA
#ifdef use_libMPI
      call MPI_INIT(error)
      call MPI_COMM_RANK ( MPI_COMM_WORLD, pe,        error )
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, npes,      error )
      call MPI_COMM_GROUP( MPI_COMM_WORLD, group_all, error )
#endif
      mpp_initialized = .TRUE.

      if( PRESENT(flags) )verbose = flags.EQ.MPP_VERBOSE
!messages
      if( verbose )call mpp_error( NOTE, 'MPP_INIT: initializing MPP module...' )
      if( pe.EQ.0 )then
          write( stdout,'(/a)' )'MPP module '//trim(version)
          write( stdout,* )'NPES=', npes
#ifdef use_libSMA
          write( stdout,* )'Using SMA (shmem) library for message passing...'
          if( verbose )write( stdout,* )'MPP_INIT: WORDS_PER_REAL=', WORDS_PER_REAL
#endif
#ifdef use_libMPI
          write( stdout,* )'Using MPI library for message passing...'
#endif
      end if
      call mpp_sync()

      return
    end subroutine mpp_init

    subroutine mpp_exit()
!to be called at the end of a run
!not strictly required for shmem_runs (but required for MPI)

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_EXIT: You must first call mpp_init.' )
      call mpp_sync()
#ifdef use_libMPI
      call MPI_FINALIZE(error)
#endif
      return
    end subroutine mpp_exit

    function mpp_pe()
      integer :: mpp_pe
      mpp_pe = pe
      return
    end function mpp_pe

    function mpp_npes()
      integer :: mpp_npes
      mpp_npes = npes
      return
    end function mpp_npes

    subroutine make_pe_set(pelist,peset)
!makes a PE set out of a PE list (list length must be .GE.2)
!a PE list is an ordered list of PEs
!a PE set is a triad (start,log2stride,size) for SHMEM, an a communicator for MPI (other two elements are unused)
!if stride is non-uniform or not a power of 2, SHMEM version will return error
      integer, intent(in) :: pelist(0:)
      integer, intent(out) :: peset(3)
#ifdef use_libSMA
      integer :: i, n, stride
      n = size(pelist)
      peset(1) = pelist(0)
      stride = pelist(1)-pelist(0)
      if( ANY(pelist(2:n-1)-pelist(1:n-2).NE.stride) ) &
           call mpp_error( FATAL, 'MAKE_PE_SET: pelist must have constant stride for SHMEM.' )
      peset(2) = log(float(stride))/log(2.)
      if( 2**peset(2).NE.stride )call mpp_error( FATAL, 'MAKE_PE_SET: pelist must have power-of-2 stride for SHMEM.' )
      peset(3) = n
#endif
#ifdef use_libMPI
      integer :: group
      call MPI_GROUP_INCL( group_all, size(pelist), pelist, group, error )
      call MPI_COMM_CREATE( MPI_COMM_WORLD, group, peset(1), error )
#endif

      return
    end subroutine make_pe_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!                BASIC MESSAGE PASSING ROUTINE: mpp_transmit                   !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_transmit_real8( put_data, put_len, put_pe, get_data, get_len, get_pe )
!a message-passing routine intended to be reminiscent equally of both MPI and SHMEM

!put_data and get_data are contiguous real*8 arrays

!at each call, your put_data array is put to   put_pe's get_data
!              your get_data array is got from get_pe's put_data
!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get

!special PE designations:
!      NULL_PE: to disable a put or a get (e.g at boundaries)
!      ANY_PE:  if remote PE for the put or get is to be unspecific
!      ALL_PES: broadcast and collect operations (collect not yet implemented)

!ideally we would not pass length, but this f77-style call performs better (arrays passed by address, not descriptor)
!further, this permits <length> contiguous words from an array of any rank to be passed (avoiding f90 rank conformance check)

!caller is responsible for completion checks before and after

!there are overloaded functions below for datatypes other than real*8; they fold into this call

      integer, intent(in) :: put_len, put_pe, get_len, get_pe
      real(DOUBLE_KIND), intent(in),  dimension(put_len) :: put_data
      real(DOUBLE_KIND), intent(out) :: get_data(get_len)
#ifdef use_libSMA
      integer :: np
      integer(LONG_KIND) :: data_loc
!pointer to remote data
      real(DOUBLE_KIND) :: remote_data(get_len)
      pointer( ptr$remote_data, remote_data )
      real(DOUBLE_KIND) :: broadcast_data(get_len)
      pointer( ptr$broadcast_data, broadcast_data )
      integer :: len$broadcast_data=0
      save :: ptr$broadcast_data, len$broadcast_data
#endif

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( put_pe.EQ.NULL_PE .AND. get_pe.EQ.NULL_PE )return
      
      if( verbose )then
          call SYSTEM_CLOCK(tick)
          write( stdout,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: put_pe, get_pe, put_len, get_len=', put_pe, get_pe, put_len, get_len
      end if

#ifdef use_libSMA
!tell get_pe we are expecting data
      if( get_pe.GE.0 .AND. get_pe.LT.npes )then
          data_loc = LOC(get_data)
          call SHMEM_PUT8( mpp_ready_to_recv(pe), data_loc, 1, get_pe )
      else if( get_pe.EQ.ANY_PE )then
          data_loc = LOC(get_data)
          do np=0,npes-1
             if( pe.NE.np )call SHMEM_PUT8( mpp_ready_to_recv(pe), data_loc, 1, np )
          end do
      else if( put_pe.EQ.pe )then
          mpp_ready_to_recv(pe) = MPP_WAIT_VALUE + 1 !some value not equal to MPP_WAIT_VALUE
      end if
#endif

!do put first and then get
      if( put_pe.GE.0 .AND. put_pe.LT.npes )then !announce to put_pe that data is available on pe <pe>
#ifdef use_libSMA
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          call SHMEM_INT8_WAIT( mpp_ready_to_recv(put_pe), MPP_WAIT_VALUE )
!at this point mpp_ready_to_recv(put_pe) contains location of put_pe's get_data
!we are not currently using this information, since we acquire data by get
!we just wait for it to be .NE.MPP_WAIT_VALUE
          mpp_ready_to_recv(pe) = MPP_WAIT_VALUE !reset
          data_loc = LOC(put_data)
          put_is_done = MPP_PUT_WAIT
          call SHMEM_INTEGER_PUT( mpp_get_pe, pe, 1, put_pe )
          call SHMEM_PUT8( remote_data_ptr(pe), data_loc, 1, put_pe )
#endif use_libSMA
#ifdef use_libMPI
!use non-blocking sends
          call MPI_ISEND( put_data, put_len, MPI_REAL8, put_pe, tag, MPI_COMM_WORLD, request(put_pe), error )
#endif

      else if( put_pe.EQ.ALL_PES )then !this is a broadcast from get_pe
          if( get_pe.LT.0 .OR. get_pe.GE.npes )call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len )call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
#ifdef use_libSMA
!since we're unable to detect (on Origin at any rate) whether array a is symmetric or not, we have assumed not
          if( get_len.GT.WORK_ARRAY_SIZE )then
#ifdef SGICRAY_MPP
              call mpp_malloc( ptr$broadcast_data, get_len, len$broadcast_data )
#else
              if( pe.EQ.0 )write( stderr,'(/2(a,i10))' ) &
                   'SHMEM collective operations are restricted to arrays of max length WORK_ARRAY_SIZE=', WORK_ARRAY_SIZE, &
                   'To accommodate the current request, recompile with -DWORK_ARRAY_SIZE=', get_len
              call ABORT()
#endif
          else
              ptr$broadcast_data = LOC(symm_work_array)
          end if
          if( npes.GT.1 )then
              broadcast_data(1:get_len) = put_data(1:get_len)
              call mpp_sync()
#ifdef _CRAYT90
              call SHMEM_UDCFLUSH !invalidate data cache
#endif
              call SHMEM_BROADCAST8( broadcast_data, broadcast_data, get_len, get_pe, 0,0,npes, bSync )
              call mpp_sync()
              get_data(1:get_len) = broadcast_data(1:get_len)
          end if
#endif
!SHMEM_BROADCAST does not broadcast to itself, need copies if get_data and put_data are not the same.
!likewise MPI_BCAST only copies itself, so need to make copy prior to calling MPI_BCAST.
!thus this copy operation is placed inbetween.
          if( pe.EQ.get_pe )then
#ifdef use_CRI_pointers
!dir$ IVDEP
              if( LOC(get_data).NE.LOC(put_data) ) &
#endif
                   get_data(1:put_len) = put_data(1:put_len)
          end if
#ifdef use_libMPI
          if( npes.GT.1 )call MPI_BCAST( get_data, put_len, MPI_REAL8, get_pe, MPI_COMM_WORLD, error )
#endif
          return

      else if( put_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets
#ifdef use_libSMA
          if( get_pe.LT.0 .OR. get_pe.GE.npes )call mpp_error( FATAL, 'MPP_TRANSMIT: invalid get_pe along with put_pe=ANY_PE.' )
          call SHMEM_GET8( get_data, put_data, get_len, get_pe )
          call SHMEM_INTEGER_PUT( put_is_done, pe, 1, get_pe )
          return
#endif
#ifdef use_libMPI
!...but you cannot have a pure get with MPI
          call mpp_error( FATAL, 'MPP_TRANSMIT: you cannot transmit to ANY_PE using MPI.' )
#endif

      else if( put_pe.NE.NULL_PE )then	!no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid put_pe.' )
      end if

!do the get: for libSMA, a get means do a wait to ensure put on remote PE is complete
      if( get_pe.GE.0 .AND. get_pe.LT.npes )then
#ifdef use_libSMA
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          call SHMEM_INT8_WAIT( remote_data_ptr(get_pe), MPP_WAIT_VALUE )
          ptr$remote_data = remote_data_ptr(get_pe)
          remote_data_ptr(get_pe) = MPP_WAIT_VALUE !reset
#if defined(CRAYPVP) || defined(sgi_mipspro)
!since we have the pointer to remote data, just retrieve it with a simple copy
!dir$ IVDEP
          if( LOC(get_data).NE.LOC(remote_data) )get_data(1:get_len) = remote_data(1:get_len)
#else
          call SHMEM_GET8( get_data, remote_data, get_len, get_pe )
#endif
          call SHMEM_INTEGER_PUT( put_is_done, pe, 1, get_pe )
#elif  use_libMPI
!receive from get_pe
          call MPI_RECV( get_data, get_len, MPI_REAL8, get_pe,         MPI_ANY_TAG, MPI_COMM_WORLD, stat, error )
#else !neither use_libSMA nor use_libMPI
          if( pe.EQ.get_pe )then
#ifdef use_CRI_pointers
!dir$ IVDEP
              if( LOC(get_data).NE.LOC(put_data) ) &
#endif
                   get_data(1:put_len) = put_data(1:put_len)
          end if
#endif

      else if( get_pe.EQ.ANY_PE )then
#ifdef use_libSMA
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
!since we don't know which PE is sending us data, we wait for remote PE to send us its ID
!this is only required for !CRAYPVP  && !sgi_mipspro, but is done there too, so that we can send put_is_done back.
          call shmem_integer_wait( mpp_get_pe, ANY_PE )
          call SHMEM_INT8_WAIT( remote_data_ptr(mpp_get_pe), MPP_WAIT_VALUE )
          ptr$remote_data = remote_data_ptr(mpp_get_pe)
          remote_data_ptr(mpp_get_pe) = MPP_WAIT_VALUE !reset
#if defined(CRAYPVP) || defined(sgi_mipspro)
!since we have the pointer to remote data, just retrieve it with a simple copy
!dir$ IVDEP
          if( LOC(get_data).NE.LOC(remote_data) )get_data(1:get_len) = remote_data(1:get_len)
#else
          call SHMEM_GET8( get_data, remote_data, get_len, mpp_get_pe )
#endif
          call SHMEM_INTEGER_PUT( put_is_done, pe, 1, mpp_get_pe )
          mpp_get_pe = ANY_PE   !reset
#endif use_libSMA
#ifdef use_libMPI
!receive from MPI_ANY_SOURCE
          call MPI_RECV( get_data, get_len, MPI_REAL8, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, stat, error )
#endif

      else if( get_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: get_pe=ALL_PES has ambiguous meaning, and hence is not implemented.' )

      else if( get_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid get_pe.' )
      end if

      if( verbose )then
          call SYSTEM_CLOCK(tick)
          write( stdout,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: put_pe, get_pe, put_len, get_len=', put_pe, get_pe, put_len, get_len
      end if
      return
    end subroutine mpp_transmit_real8

    subroutine mpp_recv_real8( get_data, get_len, get_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, get_pe
      real(DOUBLE_KIND), intent(out) :: get_data(get_len)
      real(DOUBLE_KIND) :: a(1)               !dummy real
      call mpp_transmit( a, 1, NULL_PE, get_data, get_len, get_pe )
    end subroutine mpp_recv_real8

    subroutine mpp_send_real8( put_data, put_len, put_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, put_pe
      real(DOUBLE_KIND), intent(in) :: put_data(put_len)
      real(DOUBLE_KIND) :: a(1)               !dummy real
      call mpp_transmit( put_data, put_len, put_pe, a, 1, NULL_PE )
    end subroutine mpp_send_real8

    subroutine mpp_recv_real8_scalar( get_data, get_len, get_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, get_pe
      real(DOUBLE_KIND), intent(out) :: get_data
      real(DOUBLE_KIND) :: get_data_r8(get_len)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, get_data_r8 )
      ptr = LOC(get_data)
      call mpp_transmit( a, 1, NULL_PE, get_data_r8, get_len, get_pe )
#else
      call mpp_error( FATAL, 'MPP_RECV_REAL8_SCALAR currently requires CRI pointers.' )
#endif
    end subroutine mpp_recv_real8_scalar

    subroutine mpp_send_real8_scalar( put_data, put_len, put_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, put_pe
      real(DOUBLE_KIND), intent(in) :: put_data
      real(DOUBLE_KIND) :: put_data_r8(put_len)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, put_data_r8 )
      ptr = LOC(put_data)
      call mpp_transmit( put_data_r8, put_len, put_pe, a, 1, NULL_PE )
#else
      call mpp_error( FATAL, 'MPP_SEND_REAL8_SCALAR currently requires CRI pointers.' )
#endif
    end subroutine mpp_send_real8_scalar

    subroutine mpp_recv_int8( get_data, get_len, get_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, get_pe
      integer(LONG_KIND), intent(out) :: get_data(get_len)
      real(DOUBLE_KIND) :: get_data_r8(get_len)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, get_data_r8 )
      ptr = LOC(get_data)
      call mpp_transmit( a, 1, NULL_PE, get_data_r8, get_len, get_pe )
#else
      call mpp_error( FATAL, 'MPP_RECV_INT8 currently requires CRI pointers.' )
#endif
    end subroutine mpp_recv_int8

    subroutine mpp_send_int8( put_data, put_len, put_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, put_pe
      integer(LONG_KIND), intent(in) :: put_data(put_len)
      real(DOUBLE_KIND) :: put_data_r8(put_len)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, put_data_r8 )
      ptr = LOC(put_data)
      call mpp_transmit( put_data_r8, put_len, put_pe, a, 1, NULL_PE )
#else
      call mpp_error( FATAL, 'MPP_SEND_INT8 currently requires CRI pointers.' )
#endif
    end subroutine mpp_send_int8

    subroutine mpp_recv_int8_scalar( get_data, get_len, get_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, get_pe
      integer(LONG_KIND), intent(out) :: get_data
      real(DOUBLE_KIND) :: get_data_r8(get_len)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, get_data_r8 )
      ptr = LOC(get_data)
      call mpp_transmit( a, 1, NULL_PE, get_data_r8, get_len, get_pe )
#else
      call mpp_error( FATAL, 'MPP_RECV_INT8_SCALAR currently requires CRI pointers.' )
#endif
    end subroutine mpp_recv_int8_scalar

    subroutine mpp_send_int8_scalar( put_data, put_len, put_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, put_pe
      integer(LONG_KIND), intent(in) :: put_data
      real(DOUBLE_KIND) :: put_data_r8(put_len)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, put_data_r8 )
      ptr = LOC(put_data)
      call mpp_transmit( put_data_r8, put_len, put_pe, a, 1, NULL_PE )
#else
      call mpp_error( FATAL, 'MPP_SEND_INT8_SCALAR currently requires CRI pointers.' )
#endif
    end subroutine mpp_send_int8_scalar

    subroutine mpp_recv_cmplx8( get_data, get_len, get_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, get_pe
      complex(DOUBLE_KIND), intent(out) :: get_data(get_len)
      real(DOUBLE_KIND) :: get_data_r8(get_len*2)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, get_data_r8 )
      ptr = LOC(get_data)
      call mpp_transmit( a, 1, NULL_PE, get_data_r8, get_len*2, get_pe )
#else
      call mpp_error( FATAL, 'MPP_RECV_CMPLX8 currently requires CRI pointers.' )
#endif
    end subroutine mpp_recv_cmplx8

    subroutine mpp_send_cmplx8( put_data, put_len, put_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, put_pe
      complex(DOUBLE_KIND), intent(in) :: put_data(put_len)
      real(DOUBLE_KIND) :: put_data_r8(put_len*2)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, put_data_r8 )
      ptr = LOC(put_data)
      call mpp_transmit( put_data_r8, put_len*2, put_pe, a, 1, NULL_PE )
#else
      call mpp_error( FATAL, 'MPP_SEND_CMPLX8 currently requires CRI pointers.' )
#endif
    end subroutine mpp_send_cmplx8

    subroutine mpp_recv_cmplx8_scalar( get_data, get_len, get_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, get_pe
      complex(DOUBLE_KIND), intent(out) :: get_data
      real(DOUBLE_KIND) :: get_data_r8(get_len*2)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, get_data_r8 )
      ptr = LOC(get_data)
      call mpp_transmit( a, 1, NULL_PE, get_data_r8, get_len*2, get_pe )
#else
      call mpp_error( FATAL, 'MPP_RECV_CMPLX8_SCALAR currently requires CRI pointers.' )
#endif
    end subroutine mpp_recv_cmplx8_scalar

    subroutine mpp_send_cmplx8_scalar( put_data, put_len, put_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, put_pe
      complex(DOUBLE_KIND), intent(in) :: put_data
      real(DOUBLE_KIND) :: put_data_r8(put_len*2)
      real(DOUBLE_KIND) :: a(1)               !dummy real
#ifdef use_CRI_pointers
      pointer( ptr, put_data_r8 )
      ptr = LOC(put_data)
      call mpp_transmit( put_data_r8, put_len*2, put_pe, a, 1, NULL_PE )
#else
      call mpp_error( FATAL, 'MPP_SEND_CMPLX8_SCALAR currently requires CRI pointers.' )
#endif
    end subroutine mpp_send_cmplx8_scalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!       GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min                   !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_max_real8( a, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      real(DOUBLE_KIND), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      real(DOUBLE_KIND) :: b
      integer :: peset(3)
#ifdef use_libSMA
      real :: pWrk(SHMEM_REDUCE_MIN_WRKDATA_SIZE)
#ifdef SGICRAY_MPP
      pointer( ptr$pWrk, pWrk )
      save :: ptr$pWrk
      integer, save :: len$pWrk=0
#endif
      pointer( ptr$b, b )
#endif use_libSMA

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_MAX: You must first call mpp_init.' )
      if( npes.EQ.1 )return

      if( PRESENT(pelist) )then
          if( size(pelist).EQ.1 )return
          call make_pe_set( pelist, peset )
      else
#ifdef use_libSMA
          peset(1) = 0
          peset(2) = 0
          peset(3) = npes
#endif
#ifdef use_libMPI
          peset(1) = MPI_COMM_WORLD
#endif
      end if
          
#ifdef use_libSMA
      if( verbose )call mpp_error( NOTE, 'MPP_MAX: using SHMEM_REAL8_MAX_TO_ALL...' )
#ifdef SGICRAY_MPP
      call mpp_malloc( ptr$pWrk, size(pWrk), len$pWrk )
#endif
!since we're unable to detect (on Origin at any rate) whether array a is symmetric or not, we have assumed not
      ptr$b = LOC(symm_work_array)
      b = a
      call SHMEM_REAL8_MAX_TO_ALL( b, b, 1, peset(1), peset(2), peset(3), pWrk, rSync )
      a = b
      call SHMEM_BARRIER_ALL()
#endif use_libSMA
#ifdef use_libMPI
!sum on pe 0 and broadcast
      if( verbose )call mpp_error( NOTE, 'MPP_MAX: using MPI_REDUCE and MPI_BCAST...' )
      call MPI_REDUCE( a, b, 1, MPI_REAL8, MPI_MAX, 0, peset(1), error )
      a = b
      call MPI_BCAST( a, 1, MPI_REAL8, 0, peset(1), error )
#endif
      return
    end subroutine mpp_max_real8

    subroutine mpp_max_int8( a, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      integer(LONG_KIND), intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      integer(LONG_KIND) :: b
      integer :: peset(3)
#ifdef use_libSMA
      integer(LONG_KIND) :: pWrk(SHMEM_REDUCE_MIN_WRKDATA_SIZE)
#ifdef SGICRAY_MPP
      pointer( ptr$pWrk, pWrk )
      save :: ptr$pWrk
      integer, save :: len$pWrk=0
#endif
      pointer( ptr$b, b )
#endif use_libSMA

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_MAX: You must first call mpp_init.' )
      if( npes.EQ.1 )return

      if( PRESENT(pelist) )then
          if( size(pelist).EQ.1 )return
          call make_pe_set( pelist, peset )
      else
#ifdef use_libSMA
          peset(1) = 0
          peset(2) = 0
          peset(3) = npes
#endif
#ifdef use_libMPI
          peset(1) = MPI_COMM_WORLD
#endif
      end if
          
#ifdef use_libSMA
      if( verbose )call mpp_error( NOTE, 'MPP_MAX: using SHMEM_INT8_MAX_TO_ALL...' )
#ifdef SGICRAY_MPP
      call mpp_malloc( ptr$pWrk, size(pWrk), len$pWrk )
#endif
!since we're unable to detect (on Origin at any rate) whether array a is symmetric or not, we have assumed not
      ptr$b = LOC(symm_work_array)
      b = a
      call SHMEM_INT8_MAX_TO_ALL( b, b, 1, peset(1), peset(2), peset(3), pWrk, rSync )
      a = b
      call SHMEM_BARRIER_ALL()
#endif use_libSMA
#ifdef use_libMPI
!sum on pe 0 and broadcast
      if( verbose )call mpp_error( NOTE, 'MPP_MAX: using MPI_REDUCE and MPI_BCAST...' )
      call MPI_REDUCE( a, b, 1, MPI_INTEGER8, MPI_MAX, 0, peset(1), error )
      a = b
      call MPI_BCAST( a, 1, MPI_INTEGER8, 0, peset(1), error )
#endif
      return
    end subroutine mpp_max_int8

    subroutine mpp_min_real8( a, pelist )
!find the min of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      real(DOUBLE_KIND), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      real(DOUBLE_KIND) :: b
      integer :: peset(3)
#ifdef use_libSMA
      real :: pWrk(SHMEM_REDUCE_MIN_WRKDATA_SIZE)
#ifdef SGICRAY_MPP
      pointer( ptr$pWrk, pWrk )
      save :: ptr$pWrk
      integer, save :: len$pWrk=0
#endif
      pointer( ptr$b, b )
#endif use_libSMA

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_MIN: You must first call mpp_init.' )
      if( npes.EQ.1 )return

      if( PRESENT(pelist) )then
          if( size(pelist).EQ.1 )return
          call make_pe_set( pelist, peset )
      else
#ifdef use_libSMA
          peset(1) = 0
          peset(2) = 0
          peset(3) = npes
#endif
#ifdef use_libMPI
          peset(1) = MPI_COMM_WORLD
#endif
      end if
          
#ifdef use_libSMA
      if( verbose )call mpp_error( NOTE, 'MPP_MIN: using SHMEM_REAL8_MIN_TO_ALL...' )
#ifdef SGICRAY_MPP
      call mpp_malloc( ptr$pWrk, size(pWrk), len$pWrk )
#endif
!since we're unable to detect (on Origin at any rate) whether array a is symmetric or not, we have assumed not
      ptr$b = LOC(symm_work_array)
      b = a
      call SHMEM_REAL8_MIN_TO_ALL( b, b, 1, peset(1), peset(2), peset(3), pWrk, rSync )
      a = b
      call SHMEM_BARRIER_ALL()
#endif use_libSMA
#ifdef use_libMPI
!sum on pe 0 and broadcast
      if( verbose )call mpp_error( NOTE, 'MPP_MIN: using MPI_REDUCE and MPI_BCAST...' )
      call MPI_REDUCE( a, b, 1, MPI_REAL8, MPI_MIN, 0, peset(1), error )
      a = b
      call MPI_BCAST( a, 1, MPI_REAL8, 0, peset(1), error )
#endif
      return
    end subroutine mpp_min_real8

    subroutine mpp_min_int8( a, pelist )
!find the min of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      integer(LONG_KIND), intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      integer(LONG_KIND) :: b
      integer :: peset(3)
#ifdef use_libSMA
      integer(LONG_KIND) :: pWrk(SHMEM_REDUCE_MIN_WRKDATA_SIZE)
#ifdef SGICRAY_MPP
      pointer( ptr$pWrk, pWrk )
      save :: ptr$pWrk
      integer, save :: len$pWrk=0
#endif
      pointer( ptr$b, b )
#endif use_libSMA

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_MIN: You must first call mpp_init.' )
      if( npes.EQ.1 )return

      if( PRESENT(pelist) )then
          if( size(pelist).EQ.1 )return
          call make_pe_set( pelist, peset )
      else
#ifdef use_libSMA
          peset(1) = 0
          peset(2) = 0
          peset(3) = npes
#endif
#ifdef use_libMPI
          peset(1) = MPI_COMM_WORLD
#endif
      end if
          
#ifdef use_libSMA
      if( verbose )call mpp_error( NOTE, 'MPP_MIN: using SHMEM_INT8_MIN_TO_ALL...' )
#ifdef SGICRAY_MPP
      call mpp_malloc( ptr$pWrk, size(pWrk), len$pWrk )
#endif
!since we're unable to detect (on Origin at any rate) whether array a is symmetric or not, we have assumed not
      ptr$b = LOC(symm_work_array)
      b = a
      call SHMEM_INT8_MIN_TO_ALL( b, b, 1, peset(1), peset(2), peset(3), pWrk, rSync )
      a = b
      call SHMEM_BARRIER_ALL()
#endif use_libSMA
#ifdef use_libMPI
!sum on pe 0 and broadcast
      if( verbose )call mpp_error( NOTE, 'MPP_MIN: using MPI_REDUCE and MPI_BCAST...' )
      call MPI_REDUCE( a, b, 1, MPI_INTEGER8, MPI_MIN, 0, peset(1), error )
      a = b
      call MPI_BCAST( a, 1, MPI_INTEGER8, 0, peset(1), error )
#endif
      return
    end subroutine mpp_min_int8

    subroutine mpp_sum_real8( a, length, pelist )
!sums array a over the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast: all PEs have the sum in a at the end
!we are using f77-style call: array passed by address and not descriptor; further, the f90 conformance check is avoided.
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      real(DOUBLE_KIND), intent(inout) :: a(length)
      integer :: level, lognpes, pedist, vpe, vpes
      real(DOUBLE_KIND), dimension(length) :: b, c !work arrays
      integer :: peset(3)
#ifdef use_libSMA
#ifdef bit_reproducible
#ifdef use_shmalloc
      pointer( ptr$b, b )
      save :: ptr$b
      integer, save :: len$b=0
      pointer( ptr$c, c )
      save :: ptr$c
      integer, save :: len$c=0
#endif
#else !bit_reproducible
      real :: pWrk(length/2+1+SHMEM_REDUCE_MIN_WRKDATA_SIZE)
#ifdef SGICRAY_MPP
      pointer( ptr$pWrk, pWrk )
      save :: ptr$pWrk
      integer, save :: len$pWrk=0
#endif
      pointer( ptr$b, b )
      save :: ptr$b
      integer, save :: len$b=0
#endif bit_reproducible
#endif use_libSMA

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      if( npes.EQ.1 )return

      if( PRESENT(pelist) )then
          if( size(pelist).EQ.1 )return
      end if

#ifdef bit_reproducible
      if( PRESENT(pelist) )call mpp_error( FATAL, 'MPP_SUM: currently only sums over all PEs, cannot use pelist.' )
#ifdef use_shmalloc
      call mpp_malloc( ptr$b, length, len$b )
      call mpp_malloc( ptr$c, length, len$c )
#endif
      if( verbose )call mpp_error( NOTE, 'MPP_SUM: using binary tree summation...' )
      lognpes = int( log(float(npes))/log(2.) - epsilon(1.) ) + 1
      if( verbose )then
          call SYSTEM_CLOCK(tick)
          write( stdout,'(a,i18,a,i5,a,i5)' )'T=',tick, ' PE=', pe, ' MPP_SUM: lognpes=', lognpes
      end if
      vpes = npes
      do level = 0,lognpes-1    !level on tree
         pedist = 2**level      !distance to sum over
         vpes = vpes + mod(vpes,pedist)
         if( verbose )then
             call SYSTEM_CLOCK(tick)
             write( stdout,'(a,i18,a,i5,a,3i5)' )'T=',tick, ' PE=', pe, ' MPP_SUM: level, pedist, vpes=', level, pedist, vpes
         end if
         b(1:length) = a(1:length)                  !initialize b for each level of the tree
         call mpp_sync()
         if( verbose )then
             call SYSTEM_CLOCK(tick)
             write( stdout,'(a,i18,a,i5,a,2i5)' )&
                  'T=',tick, ' PE=', pe, ' MPP_SUM: pedist, mod(pe,pedist*2)=', pedist, mod(pe,pedist*2)
         end if
         if( mod(pe,pedist*2).GE.pedist )then
             call mpp_transmit( b, length, pe-pedist, c, length, pe-pedist )
             a(1:length) = c(1:length) + b(1:length)          !if c came from the left,  sum on the left
             if( verbose )then
                 call SYSTEM_CLOCK(tick)
                 write( stdout,'(a,i18,a,i5,a,i5,3es23.15)' )&
                      'T=',tick, ' PE=', pe, ' MPP_SUM: Lsum: remote_pe=', pe-pedist, b(1), c(1), a(1)
             end if
         else if( pe+pedist.LT.npes )then
             call mpp_transmit( b, length, pe+pedist, c, length, pe+pedist )
             a(1:length) = a(1:length) + c(1:length)          !if c came from the right, sum on the right
             if( verbose )then
                 call SYSTEM_CLOCK(tick)
                 write( stdout,'(a,i18,a,i5,a,i5,3es23.15)' )&
                      'T=',tick, ' PE=', pe, ' MPP_SUM: Rsum: remote_pe=', pe+pedist, b(1), c(1), a(1)
             end if
         else if( pe+pedist.LT.vpes )then
             vpe = vpes - 1
             do while( vpe.GE.npes )
                vpe = (pe+vpe)/2
             end do
#ifdef use_libMPI
             call mpp_error( FATAL, 'The bit-reproducible MPI version of MPP_SUM _currently_ only works on power-of-2 NPES.' )
#endif
             call mpp_transmit( b, length, ANY_PE, c, length, vpe )
             a(1:length) = a(1:length) + c(1:length)          !if c came from the right, sum on the right
             if( verbose )then
                 call SYSTEM_CLOCK(tick)
                 write( stdout,'(a,i18,a,i5,a,i5,3es23.15)' )&
                      'T=',tick, ' PE=', pe, ' MPP_SUM: Vsum: remote_pe=', pe+pedist/2, b(1), c(1), a(1)
             end if
         end if
!still need to figure out put pattern for non-power-of-2 NPES
!         if( pe+pedist/2.LT.vpes .AND. pe+pedist.GT.vpes )then
!             call mpp_transmit( b, c, size(b), pe-pedist/2, NULL_PE )
!             if( verbose )&
!                  write( stdout,* )'MPP_SUM: Vsum: pe, remote_pe=', pe, pe-pedist/2
!         end if
      end do
      call mpp_sync()
#else! bit_reproducible
      if( PRESENT(pelist) )then
          call make_pe_set(pelist,peset)
      else
#ifdef use_libSMA
          peset(1) = 0
          peset(2) = 0
          peset(3) = npes
#endif
#ifdef use_libMPI
          peset(1) = MPI_COMM_WORLD
#endif
      end if
#ifdef use_libSMA
      if( verbose )call mpp_error( NOTE, 'MPP_SUM: using SHMEM_REAL8_SUM_TO_ALL...' )
#ifdef SGICRAY_MPP
      call mpp_malloc( ptr$pWrk, size(pWrk), len$pWrk )
#endif
!since we're unable to detect (on Origin at any rate) whether array a is symmetric or not, we have assumed not
      if( length.GT.WORK_ARRAY_SIZE )then
#ifdef SGICRAY_MPP
          call mpp_malloc( ptr$b, length, len$b )
#else
          if( pe.EQ.0 )write( stderr,'(/2(a,i10))' ) &
               'SHMEM collective operations are restricted to arrays of max length WORK_ARRAY_SIZE=', WORK_ARRAY_SIZE, &
               'To accommodate the current request, recompile with -DWORK_ARRAY_SIZE=', length
          call ABORT()
#endif
      else
          ptr$b = LOC(symm_work_array)
      end if
      b(1:length) = a(1:length)
      call SHMEM_REAL8_SUM_TO_ALL( b, b, length, peset(1), peset(2), peset(3), pWrk, rSync )
      a(1:length) = b(1:length)
      call SHMEM_BARRIER_ALL()
#endif use_libSMA
#ifdef use_libMPI
!sum on pe 0 and broadcast
      if( verbose )call mpp_error( NOTE, 'MPP_SUM: using MPI_REDUCE and MPI_BCAST...' )
      call MPI_REDUCE( a, b, length, MPI_REAL8, MPI_SUM, 0, peset(1), error )
      a = b
      call MPI_BCAST( a, length, MPI_REAL8, 0, peset(1), error )
#endif
#endif bit_reproducible
      return
    end subroutine mpp_sum_real8

    subroutine mpp_sum_int8( a, pelist )
!sums long integer a over the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast: all PEs have the sum in a at the end
      integer, intent(in), optional :: pelist(:)
      integer(LONG_KIND), intent(inout) :: a
      integer :: peset(3)
      integer(LONG_KIND) :: b
#ifdef use_libSMA
      integer(LONG_KIND) :: pWrk(SHMEM_REDUCE_MIN_WRKDATA_SIZE)
#ifdef SGICRAY_MPP
      pointer( ptr$pWrk, pWrk )
      save :: ptr$pWrk
      integer, save :: len$pWrk=0
#endif
      pointer( ptr$b, b )
      save :: ptr$b
      integer, save :: len$b=0
#endif use_libSMA

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      if( npes.EQ.1 )return

      if( PRESENT(pelist) )then
          if( size(pelist).EQ.1 )return
          call make_pe_set(pelist,peset)
      else
#ifdef use_libSMA
          peset(1) = 0
          peset(2) = 0
          peset(3) = npes
#endif
#ifdef use_libMPI
          peset(1) = MPI_COMM_WORLD
#endif
      end if
#ifdef use_libSMA
      if( verbose )call mpp_error( NOTE, 'MPP_SUM: using SHMEM_REAL8_SUM_TO_ALL...' )
#ifdef SGICRAY_MPP
      call mpp_malloc( ptr$pWrk, size(pWrk), len$pWrk )
#endif
      ptr$b = LOC(symm_work_array)
      b = a
      call SHMEM_INT8_SUM_TO_ALL( b, b, 1, peset(1), peset(2), peset(3), pWrk, rSync )
      a = b
      call SHMEM_BARRIER_ALL()
#endif use_libSMA
#ifdef use_libMPI
!sum on pe 0 and broadcast
      if( verbose )call mpp_error( NOTE, 'MPP_SUM: using MPI_REDUCE and MPI_BCAST...' )
      call MPI_REDUCE( a, b, 1, MPI_INTEGER8, MPI_SUM, 0, peset(1), error )
      a = b
      call MPI_BCAST( a, 1, MPI_INTEGER8, 0, peset(1), error )
#endif
      return
    end subroutine mpp_sum_int8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!           SYNCHRONIZATION ROUTINES: mpp_sync, mpp_sync_self                  !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_sync( pelist )
!synchronize PEs in list
      integer, intent(in), optional :: pelist(:)
#ifdef use_libSMA
      mpp_ready_to_recv(0:npes-1) = MPP_WAIT_VALUE
      remote_data_ptr(0:npes-1) = MPP_WAIT_VALUE
      put_is_done = pe
      mpp_get_pe = ANY_PE
      call SHMEM_BARRIER_ALL()
#endif
#ifdef use_libMPI
      call mpp_sync_self()
#endif
      return
    end subroutine mpp_sync

    subroutine mpp_sync_self( remote_pe )
!this is to check if current PE's outstanding puts are complete
!but we can't use shmem_fence because we are actually waiting for
!a remote PE to complete its get
!so, I'm using the put_is_done variable:
!this has the value MPP_PUT_WAIT if there is an outstanding put
!otherwise contains the value of remote pe that acquired data
!typically we check is put_is_done.NE.MPP_PUT_WAIT
!can also compare with value of remote_pe if supplied
      integer, intent(in), optional :: remote_pe
      integer :: i
#ifdef use_libSMA
#ifdef _CRAYT90
      call SHMEM_UDCFLUSH !invalidate data cache
#endif
      call shmem_integer_wait( put_is_done, MPP_PUT_WAIT )
      if( PRESENT(remote_pe) )then
          if( put_is_done.NE.remote_pe )call mpp_error( FATAL, 'MPP_SYNC_SELF: communication mismatch.' )
      end if
#endif use_libSMA
#ifdef use_libMPI
      if( PRESENT(remote_pe) )then
!check if outstanding non-blocking messages to remote_pe are complete
          if( request(remote_pe).NE.MPI_REQUEST_NULL )call MPI_WAIT( request(remote_pe), stat, error )
      else
!check if outstanding non-blocking messages to all PEs are complete
          do i = 0,npes-1
             if( request(i).NE.MPI_REQUEST_NULL )call MPI_WAIT( request(i), stat, error )
          end do
      end if
#endif
      put_is_done = pe
      return
    end subroutine mpp_sync_self

#ifdef use_libSMA
!these local versions are written for grouping into shmem_integer_wait
    subroutine shmem_int4_wait_local( ivar, cmp_value )
!dir$ INLINEALWAYS shmem_int4_wait_local
      integer(INT_KIND), intent(in) :: cmp_value
      integer(INT_KIND), intent(inout) :: ivar
      call SHMEM_INT4_WAIT( ivar, cmp_value )
      return
    end subroutine shmem_int4_wait_local
    subroutine shmem_int8_wait_local( ivar, cmp_value )
!dir$ INLINEALWAYS shmem_int8_wait_local
      integer(LONG_KIND), intent(in) :: cmp_value
      integer(LONG_KIND), intent(inout) :: ivar
      call SHMEM_INT8_WAIT( ivar, cmp_value )
      return
    end subroutine shmem_int8_wait_local
#endif
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!MISCELLANEOUS UTILITIES: mpp_error, mpp_chksum, mpp_malloc, mpp_transmissible !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_error( errortype, errormsg )
!a very basic error handler
      integer, intent(in) :: errortype
      character(len=*), intent(in), optional :: errormsg
      character(len=256) :: text
      logical :: opened
      
      if( .NOT.mpp_initialized )call abort()
      if( errortype.EQ.NOTE    )text = 'NOTE'	!just FYI
      if( errortype.EQ.WARNING )text = 'WARNING'	!probable error
      if( errortype.EQ.FATAL   )text = 'FATAL'	!fatal error

      if( errortype.EQ.FATAL )call FLUSH(stdout)

      if( npes.GT.1 )write( text,'(a,i5)' )trim(text)//' from PE', pe	!this is the mpp part
      if( PRESENT(errormsg) )text = trim(text)//': '//trim(errormsg)

      if( errortype.NE.NOTE )then
          write( stderr,'(/a/)' )trim(text)
          if( errortype.EQ.FATAL .OR. warnings_are_fatal )then
#ifdef sgi_mipspro
              call TRACE_BACK_STACK_AND_PRINT()
#endif
              call ABORT()	!automatically calls traceback on Cray systems
          end if
      end if
      return
    end subroutine mpp_error

    subroutine mpp_set_warn_level(flag)
      integer, intent(in) :: flag

      if( flag.EQ.WARNING )then
          warnings_are_fatal = .FALSE.
      else if( flag.EQ.FATAL )then
          warnings_are_fatal = .TRUE.
      else
          call mpp_error( FATAL, 'MPP_SET_WARN_LEVEL: warning flag must be set to WARNING or FATAL.' )
      end if
      return
    end subroutine mpp_set_warn_level

#ifdef use_CRI_pointers
    subroutine mpp_malloc( ptr, newlen, len )
!routine to perform symmetric allocation:
!this is required on the t3e/O2k for variables that will be non-local arguments
!to a shmem call (see man intro_shmem(3F)).
!newlen is the required allocation length for the pointer ptr
!   len is the current allocation (0 if unallocated)
      integer, intent(in) :: newlen
      integer, intent(inout) :: len
      real :: dummy
!argument ptr is a cray pointer, points to a dummy argument in this routine
      pointer( ptr, dummy )
      integer(LONG_KIND) :: error_8

#if defined(SGICRAY_MPP) && defined(use_libSMA)
      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_MALLOC: You must first call mpp_init.' )
!use existing allocation if it is enough
      if( newlen.LE.len )return

      call SHMEM_BARRIER_ALL()
!if the pointer is already allocated, deallocate
      if( len.NE.0 )call SHPDEALLC( ptr, error_8, -1 ) !BWA: error_8 instead of error, see PV 682618 (fixed in mpt.1.3.0.1)
!allocate new length
      call SHPALLOC( ptr, newlen*WORDS_PER_REAL, error, -1 )
      len = newlen
      call SHMEM_BARRIER_ALL()

      if( verbose )then
          call SYSTEM_CLOCK(tick)
          write( stdout,'(a,i18,a,i5,a,2i8,i16)' )'T=', tick, ' PE=', pe, ' MPP_MALLOC: len, newlen, ptr=', len, newlen, ptr
      end if
#endif      SGICRAY_MPP  &&         use_libSMA
      return
    end subroutine mpp_malloc
#endif use_CRI_pointers

    function mpp_chksum_int_1d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_int_1d
      integer(LONG_KIND), intent(in) :: var(:)
      integer, optional :: pelist(:)
      mpp_chksum_int_1d = sum(var)
      call mpp_sum( mpp_chksum_int_1d, pelist )
      return
    end function mpp_chksum_int_1d

    function mpp_chksum_int_2d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_int_2d
      integer(LONG_KIND), intent(in) :: var(:,:)
      integer, optional :: pelist(:)
      mpp_chksum_int_2d = sum(var)
      call mpp_sum( mpp_chksum_int_2d, pelist )
      return
    end function mpp_chksum_int_2d

    function mpp_chksum_int_3d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_int_3d
      integer(LONG_KIND), intent(in) :: var(:,:,:)
      integer, optional :: pelist(:)
      mpp_chksum_int_3d = sum(var)
      call mpp_sum( mpp_chksum_int_3d, pelist )
      return
    end function mpp_chksum_int_3d

    function mpp_chksum_int_4d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_int_4d
      integer(LONG_KIND), intent(in) :: var(:,:,:,:)
      integer, optional :: pelist(:)
      mpp_chksum_int_4d = sum(var)
      call mpp_sum( mpp_chksum_int_4d, pelist )
      return
    end function mpp_chksum_int_4d

#ifdef use_CRI_pointers
    function mpp_chksum_r8_0d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_r8_0d
      real(DOUBLE_KIND), intent(in) :: var
      integer, optional :: pelist(:)
      integer(LONG_KIND) :: i(2)
      pointer( ptr, i )
      ptr = LOC(var)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_r8_0d = mpp_chksum_int_1d( i(1:2), pelist )
      else
          mpp_chksum_r8_0d = mpp_chksum_int_1d( i(1:1), pelist )
      endif
      return
    end function mpp_chksum_r8_0d

    function mpp_chksum_r8_1d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_r8_1d
      real(DOUBLE_KIND), intent(in) :: var(:)
      integer, optional :: pelist(:)
      real(DOUBLE_KIND) :: g(size(var))
      integer(LONG_KIND) :: i(size(var)*2)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_r8_1d = mpp_chksum_int_1d( i(1:size(var)*2), pelist )
      else
          mpp_chksum_r8_1d = mpp_chksum_int_1d( i(1:size(var)), pelist )
      endif
      return
    end function mpp_chksum_r8_1d

    function mpp_chksum_r8_2d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_r8_2d
      real(DOUBLE_KIND), intent(in) :: var(:,:)
      integer, optional :: pelist(:)
      real(DOUBLE_KIND) :: g(size(var,1),size(var,2))
      integer(LONG_KIND) :: i(size(var)*2)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_r8_2d = mpp_chksum_int_1d( i(1:size(var)*2), pelist )
      else
          mpp_chksum_r8_2d = mpp_chksum_int_1d( i(1:size(var)), pelist )
      endif
      return
    end function mpp_chksum_r8_2d

    function mpp_chksum_r8_3d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_r8_3d
      real(DOUBLE_KIND), intent(in) :: var(:,:,:)
      integer, optional :: pelist(:)
      real(DOUBLE_KIND) :: g(size(var,1),size(var,2),size(var,3))
      integer(LONG_KIND) :: i(size(var)*2)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_r8_3d = mpp_chksum_int_1d( i(1:size(var)*2), pelist )
      else
          mpp_chksum_r8_3d = mpp_chksum_int_1d( i(1:size(var)), pelist )
      endif
      return
    end function mpp_chksum_r8_3d

    function mpp_chksum_r8_4d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_r8_4d
      real(DOUBLE_KIND), intent(in) :: var(:,:,:,:)
      integer, optional :: pelist(:)
      real(DOUBLE_KIND) :: g(size(var,1),size(var,2),size(var,3),size(var,4))
     integer(LONG_KIND) :: i(size(var)*2)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_r8_4d = mpp_chksum_int_1d( i(1:size(var)*2), pelist )
      else
          mpp_chksum_r8_4d = mpp_chksum_int_1d( i(1:size(var)), pelist )
      endif
      return
    end function mpp_chksum_r8_4d

    function mpp_chksum_c8_0d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_c8_0d
      complex(DOUBLE_KIND), intent(in) :: var
      integer, optional :: pelist(:)
      integer(LONG_KIND) :: i(4)
      pointer( ptr, i )
      ptr = LOC(var)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_c8_0d = mpp_chksum_int_1d( i(1:4), pelist )
      else
          mpp_chksum_c8_0d = mpp_chksum_int_1d( i(1:2), pelist )
      endif
      return
    end function mpp_chksum_c8_0d

    function mpp_chksum_c8_1d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_c8_1d
      complex(DOUBLE_KIND), intent(in) :: var(:)
      integer, optional :: pelist(:)
      complex(DOUBLE_KIND) :: g(size(var))
      integer(LONG_KIND) :: i(size(var)*4)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_c8_1d = mpp_chksum_int_1d( i(1:size(var)*4), pelist )
      else
          mpp_chksum_c8_1d = mpp_chksum_int_1d( i(1:size(var)*2), pelist )
      endif
      return
    end function mpp_chksum_c8_1d

    function mpp_chksum_c8_2d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_c8_2d
      complex(DOUBLE_KIND), intent(in) :: var(:,:)
      integer, optional :: pelist(:)
      complex(DOUBLE_KIND) :: g(size(var,1),size(var,2))
      integer(LONG_KIND) :: i(size(var)*4)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_c8_2d = mpp_chksum_int_1d( i(1:size(var)*4), pelist )
      else
          mpp_chksum_c8_2d = mpp_chksum_int_1d( i(1:size(var)*2), pelist )
      endif
      return
    end function mpp_chksum_c8_2d

    function mpp_chksum_c8_3d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_c8_3d
      complex(DOUBLE_KIND), intent(in) :: var(:,:,:)
      integer, optional :: pelist(:)
      complex(DOUBLE_KIND) :: g(size(var,1),size(var,2),size(var,3))
      integer(LONG_KIND) :: i(size(var)*4)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_c8_3d = mpp_chksum_int_1d( i(1:size(var)*4), pelist )
      else
          mpp_chksum_c8_3d = mpp_chksum_int_1d( i(1:size(var)*2), pelist )
      endif
      return
    end function mpp_chksum_c8_3d

    function mpp_chksum_c8_4d( var, pelist )
      integer(LONG_KIND) :: mpp_chksum_c8_4d
      complex(DOUBLE_KIND), intent(in) :: var(:,:,:,:)
      integer, optional :: pelist(:)
      complex(DOUBLE_KIND) :: g(size(var,1),size(var,2),size(var,3),size(var,4))
      integer(LONG_KIND) :: i(size(var)*4)
      pointer( ptr, i )
      g = var                   !g is guaranteed contiguous
      ptr = LOC(g)
      if( INT_KIND.EQ.LONG_KIND )then
!8-byte ints are not permitted
          mpp_chksum_c8_4d = mpp_chksum_int_1d( i(1:size(var)*4), pelist )
      else
          mpp_chksum_c8_4d = mpp_chksum_int_1d( i(1:size(var)*2), pelist )
      endif
      return
    end function mpp_chksum_c8_4d
#endif use_CRI_pointers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!       OVERLOADED FUNCTIONS CONGRUENT TO VARIOUS ROUTINES ABOVE               !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_transmit_real8_scalar( put_data, put_len, put_pe, get_data, get_len, get_pe )
!overloaded routine for scalar arguments: convert to real 1D array and call mpp_transmit_real8
      integer, intent(in) :: put_len, put_pe, get_len, get_pe
      real(DOUBLE_KIND), intent(in)  :: put_data
      real(DOUBLE_KIND), intent(out) :: get_data
      real(DOUBLE_KIND) :: put_r8(put_len)
      real(DOUBLE_KIND) :: get_r8(get_len)
#ifdef use_CRI_pointers
      pointer( ptr$put, put_r8 )
      pointer( ptr$get, get_r8 )
      ptr$put = LOC(put_data)
      ptr$get = LOC(get_data)
      if( verbose )call mpp_error( NOTE, 'MPP_TRANSMIT_REAL8_SCALAR: calling mpp_transmit_real8...' )
      call mpp_transmit_real8( put_r8, put_len, put_pe, get_r8, get_len, get_pe )
#else
      call mpp_error( FATAL, 'MPP_TRANSMIT_REAL8_SCALAR: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_transmit_real8_scalar

    subroutine mpp_transmit_cmplx8( put_data, put_len, put_pe, get_data, get_len, get_pe )
!overloaded routine for complex arguments: convert to real 1D array and call mpp_transmit_real8
      integer, intent(in) :: put_len, put_pe, get_len, get_pe
      complex(DOUBLE_KIND), intent(in),  dimension(put_len) :: put_data
      complex(DOUBLE_KIND), intent(out) :: get_data(get_len)
      real(DOUBLE_KIND) :: put_r8(put_len*2)
      real(DOUBLE_KIND) :: get_r8(get_len*2)
#ifdef use_CRI_pointers
      pointer( ptr$put, put_r8 )
      pointer( ptr$get, get_r8 )
      ptr$put = LOC(put_data)
      ptr$get = LOC(get_data)
      if( verbose )call mpp_error( NOTE, 'MPP_TRANSMIT_CMPLX8: calling mpp_transmit_real8...' )
      call mpp_transmit_real8( put_r8, put_len*2, put_pe, get_r8, get_len*2, get_pe )
#else
      call mpp_error( FATAL, 'MPP_TRANSMIT_CMPLX8: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_transmit_cmplx8

    subroutine mpp_transmit_cmplx8_scalar( put_data, put_len, put_pe, get_data, get_len, get_pe )
!overloaded routine for scalar arguments: convert to real 1D array and call mpp_transmit_real8
      integer, intent(in) :: put_len, put_pe, get_len, get_pe
      complex(DOUBLE_KIND), intent(in)  :: put_data
      complex(DOUBLE_KIND), intent(out) :: get_data
      real(DOUBLE_KIND) :: put_r8(put_len*2)
      real(DOUBLE_KIND) :: get_r8(get_len*2)
#ifdef use_CRI_pointers
      pointer( ptr$put, put_r8 )
      pointer( ptr$get, get_r8 )
      ptr$put = LOC(put_data)
      ptr$get = LOC(get_data)
      if( verbose )call mpp_error( NOTE, 'MPP_TRANSMIT_CMPLX8_SCALAR: calling mpp_transmit_real8...' )
      call mpp_transmit_real8( put_r8, put_len*2, put_pe, get_r8, get_len*2, get_pe )
#else
      call mpp_error( FATAL, 'MPP_TRANSMIT_CMPLX8_SCALAR: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_transmit_cmplx8_scalar

    subroutine mpp_transmit_int8( put_data, put_len, put_pe, get_data, get_len, get_pe )
!overloaded routine for complex arguments: convert to real 1D array and call mpp_transmit_real8
      integer, intent(in) :: put_len, put_pe, get_len, get_pe
      integer(LONG_KIND), intent(in),  dimension(put_len) :: put_data
      integer(LONG_KIND), intent(out) :: get_data(get_len)
      real(DOUBLE_KIND) :: put_r8(put_len)
      real(DOUBLE_KIND) :: get_r8(get_len)
#ifdef use_CRI_pointers
      pointer( ptr$put, put_r8 )
      pointer( ptr$get, get_r8 )
      ptr$put = LOC(put_data)
      ptr$get = LOC(get_data)
      if( verbose )call mpp_error( NOTE, 'MPP_TRANSMIT_INT8: calling mpp_transmit_real8...' )
      call mpp_transmit_real8( put_r8, put_len, put_pe, get_r8, get_len, get_pe )
#else
      call mpp_error( FATAL, 'MPP_TRANSMIT_INT8: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_transmit_int8

    subroutine mpp_transmit_int8_scalar( put_data, put_len, put_pe, get_data, get_len, get_pe )
!overloaded routine for scalar arguments: convert to real 1D array and call mpp_transmit_real8
      integer, intent(in) :: put_len, put_pe, get_len, get_pe
      integer(LONG_KIND), intent(in)  :: put_data
      integer(LONG_KIND), intent(out) :: get_data
      real(DOUBLE_KIND) :: put_r8(put_len)
      real(DOUBLE_KIND) :: get_r8(get_len)
#ifdef use_CRI_pointers
      pointer( ptr$put, put_r8 )
      pointer( ptr$get, get_r8 )
      ptr$put = LOC(put_data)
      ptr$get = LOC(get_data)
      if( verbose )call mpp_error( NOTE, 'MPP_TRANSMIT_INT8_SCALAR: calling mpp_transmit_real8...' )
      call mpp_transmit_real8( put_r8, put_len, put_pe, get_r8, get_len, get_pe )
#else
      call mpp_error( FATAL, 'MPP_TRANSMIT_INT8_SCALAR: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_transmit_int8_scalar

    subroutine mpp_sum_cmplx8( a, length, pelist )
!sums complex array a: this routine just converts to a call to mpp_sum_real8
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      complex(DOUBLE_KIND), intent(inout) :: a(length)
      real(DOUBLE_KIND) :: b(length*2)
#ifdef use_CRI_pointers
      pointer( ptr, b )
      ptr = loc(a)
      if( verbose )call mpp_error( NOTE, 'MPP_SUM_CMPLX8: calling mpp_sum_real8...' )
      call mpp_sum_real8( b, length*2, pelist )
#else
      call mpp_error( FATAL, 'MPP_SUM_CMPLX8: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_sum_cmplx8

    subroutine mpp_sum_real8_scalar( a, length, pelist )
!sums array a when only first element is passed: this routine just converts to a call to mpp_sum_real8
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      real(DOUBLE_KIND), intent(inout) :: a
      real(DOUBLE_KIND) :: b(length)
#ifdef use_CRI_pointers
      pointer( ptr, b )
      ptr = loc(a)
      if( verbose )call mpp_error( NOTE, 'MPP_SUM_REAL8_SCALAR: calling mpp_sum_real8...' )
      call mpp_sum_real8( b, length, pelist )
#else
      call mpp_error( FATAL, 'MPP_SUM_REAL8_SCALAR: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_sum_real8_scalar

    subroutine mpp_sum_cmplx8_scalar( a, length, pelist )
!sums complex array a when only address of first element is passed: this routine just converts to a call to mpp_sum_real8
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      complex(DOUBLE_KIND), intent(inout) :: a
      real(DOUBLE_KIND) :: b(length*2)
#ifdef use_CRI_pointers
      pointer( ptr, b )
      ptr = loc(a)
      if( verbose )call mpp_error( NOTE, 'MPP_SUM_CMPLX8_SCALAR: calling mpp_sum_real8...' )
      call mpp_sum_real8( b, length*2, pelist )
#else
      call mpp_error( FATAL, 'MPP_SUM_CMPLX8_SCALAR: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_sum_cmplx8_scalar

  end module mpp_mod

#ifdef test_mpp
  program test
!test various aspects of mpp_mod
    use mpp_mod
    integer :: pe, npes
#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
    integer, parameter :: stdin=100, stdout=101, stderr=102
#else
    integer, parameter :: stdin=5, stdout=6, stderr=0
#endif
    integer, parameter :: n=1048576
    real :: a(n), b(n)
    common /junk/ a, b
    real, allocatable :: c(:)
    integer :: tick, tick0, ticks_per_sec

    call mpp_init()
    pe = mpp_pe()
    npes = mpp_npes()

    call system_clock( count_rate=ticks_per_sec )
    call random_number(a)
!time transmit, compare against shmem_put and get
    if( pe.EQ.0 )then
        print *, 'Time mpp_transmit for various lengths...'
#ifdef SGICRAY
        print *, 'For comparison, times for shmem_get and shmem_put are also provided.'
#endif
        print *
    end if
!timing is done for cyclical pass (more useful than ping-pong etc)
    l = n
    do while( l.GT.0 )
!mpp_transmit
       call system_clock(tick0)
       do nn = 1,npes
          call mpp_transmit( a, l, mod(pe+npes-1,npes), b, l, mod(pe+1,npes) )
          call mpp_sync()
       end do
       call system_clock(tick)
       dt = float(tick-tick0)/(npes*ticks_per_sec)
       if( pe.EQ.0 )write( stdout,'(/a,i8,f13.6,f8.2)' )'MPP_TRANSMIT length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
#ifdef SGICRAY
!shmem_put
       call system_clock(tick0)
       do nn = 1,npes
          call shmem_put8( b, a, l, mod(pe+1,npes) )
          call mpp_sync()
       end do
       call system_clock(tick)
       dt = float(tick-tick0)/(npes*ticks_per_sec)
       if( pe.EQ.0 )write( stdout,'( a,i8,f13.6,f8.2)' )'SHMEM_PUT    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!shmem_get
       call system_clock(tick0)
       do nn = 1,npes
          call shmem_get8( b, a, l, mod(pe+1,npes) )
          call mpp_sync()
       end do
       call system_clock(tick)
       dt = float(tick-tick0)/(npes*ticks_per_sec)
       if( pe.EQ.0 )write( stdout,'( a,i8,f13.6,f8.2)' )'SHMEM_GET    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
#endif
       l = l/2
    end do

!test mpp_sum
    if( pe.EQ.0 )then
        print *
        print *, 'Time mpp_sum...'
#ifdef bit_reproducible
        print *, 'The bit reproducibility flag for mpp_sum is turned on...'
#endif
    end if
    a = float(pe+1)
    call system_clock(tick0)
    call mpp_sum(a,n)
    call system_clock(tick)
    dt = float(tick-tick0)/ticks_per_sec
    if( pe.EQ.0 )write( stdout,'(a,2i4,f9.1,i8,f13.6,f8.2/)' ) &
         'mpp_sum: pe, npes, sum(pe+1), length, time, bw(Mb/s)=', pe, npes, a(1), n, dt, n*8e-6/dt

!test mpp_max
    if( pe.EQ.0 )then
        print *
        print *, 'Test mpp_max...'
    end if
    a = float(pe+1)
    print *, 'pe,     pe+1 =', pe, a(1)
    call mpp_max( a(1) )
    print *, 'pe, max(pe+1)=', pe, a(1)

#ifdef use_CRI_pointers
!test mpp_chksum
    if( mod(n,npes).EQ.0 )then  !only set up for even division
        call random_number(a)
        m= n/npes
        allocate( c(m) )
        c = a(pe*m+1:pe*m+m)

        if( pe.EQ.0 )then
            print *
            print *, 'Test mpp_chksum...'
            print *, 'This test shows that a whole array and a distributed array give identical checksums.'
        end if
        print *, 'chksum(a)=', mpp_chksum(a,(/pe/))
        print *, 'chksum(c)=', mpp_chksum(c)
    end if
#endif
  end program test
#endif test_mpp
