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
#include <os.h>

!only one of SMA or MPI can be used
!(though mixing calls is allowed, this module will not)
#ifdef use_libSMA
#undef use_libMPI
#endif

!shmalloc is used on MPP SGI/Cray systems for shmem
#if defined(use_libSMA) && defined(SGICRAY_MPP)
#define use_shmalloc
#endif

#ifdef __sgi
#define SYSTEM_CLOCK system_clock_sgi
#endif

module mpp_mod
!string BWA is used to tag lines that are bug workarounds and will disappear
!when offending compiler bug is fixed
!a generalized communication package for use with shmem and MPI
!will add: co_array_fortran, MPI2
!Balaji (vb@gfdl.gov) 11 May 1998
#ifdef sgi_mipspro
#ifdef use_libSMA
  use shmem_interface
#endif
#ifdef use_libMPI
  use mpi
#endif
#endif
  implicit none
  private
  character(len=128), private :: version= &
       '$Id: mpp.F90,v 6.4 2002/02/22 19:09:16 fms Exp $'
  character(len=128), private :: name= &
       '$Name: galway $'

!various lengths (see shpalloc) are estimated in "words" which are 32bit on SGI, 64bit on Cray
!these are also the expected sizeof of args to MPI/shmem libraries
#ifdef _CRAY
  integer(LONG_KIND), private :: word(1)
#endif
#ifdef sgi_mipspro
  integer(INT_KIND), private :: word(1)
#endif

#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
  integer, private :: in_unit=100, out_unit=101, err_unit=102, log_unit=103
#else
  integer, private :: in_unit=5, out_unit=6, err_unit=0, log_unit=1
#endif
  logical, private :: mpp_initialized=.FALSE.
  integer, private :: pe=0, node=0, npes=1, root_pe=0
  integer, private :: error
  integer, parameter, private :: MAXPES=2048 !used for dimensioning stuff that might be indexed by pe

!initialization flags
  integer, parameter, public :: MPP_VERBOSE=1, MPP_DEBUG=2
  logical, private :: verbose=.FALSE., debug=.FALSE.

!flags to transmit routines
  integer, parameter, public :: ALL_PES=-1, ANY_PE=-2, NULL_PE=-3

!errortype flags
  integer, parameter, public :: NOTE=0, WARNING=1, FATAL=2
  logical, private :: warnings_are_fatal = .FALSE.
  integer, private :: error_state=0

  integer(LONG_KIND), parameter, private :: MPP_WAIT=-1, MPP_READY=-2
#ifdef use_libSMA
#include <mpp/shmem.fh>
  integer :: sync(SHMEM_REDUCE_SYNC_SIZE+SHMEM_BCAST_SYNC_SIZE+SHMEM_BARRIER_SYNC_SIZE)
!status and remote_data_loc are used to synchronize communication is MPP_TRANSMIT
#ifdef use_shmalloc
  integer(LONG_KIND), private, dimension(0:MAXPES) :: status, remote_data_loc
#else
  integer(LONG_KIND), private, allocatable, dimension(:) :: status, remote_data_loc
#endif
  integer, private :: mpp_from_pe !used to announce from where data is coming from
#ifdef use_shmalloc
!we call shpalloc in mpp_init() to ensure all these are remotely accessible
!on PVP where shpalloc doesn't exist, module variables are automatically
!guaranteed to be remotely accessible
  pointer( ptr_sync, sync )
  pointer( ptr_status, status )
  pointer( ptr_from, mpp_from_pe )
  pointer( ptr_remote, remote_data_loc )
#endif
#endif use_libSMA
#ifdef use_libMPI
#ifndef sgi_mipspro
!sgi_mipspro gets this from 'use mpi'
#include <mpif.h>
#endif
!tag is never used, but must be initialized to non-negative integer
  integer, private :: tag=1, stat(MPI_STATUS_SIZE)
!  integer, private, allocatable :: request(:)
  integer, public, allocatable :: request(:)
#ifdef _CRAYT3E
!BWA: mpif.h on t3e currently does not contain MPI_INTEGER8 datatype
!(O2k and t90 do)
!(t3e: fixed on 3.3 I believe)
  integer, parameter :: MPI_INTEGER8=MPI_INTEGER
#endif
#endif use_libMPI

!mpp_stack is used by SHMEM collective ops
!must be SHPALLOC'd on SGICRAY_MPP, but is allocatable on PVP
#ifdef use_shmalloc
  real(DOUBLE_KIND), private :: mpp_stack(1)
  pointer( ptr_stack, mpp_stack )
#else
  real(DOUBLE_KIND), private, allocatable :: mpp_stack(:)
#endif
  integer, private :: mpp_stack_size=0, mpp_stack_hwm=0

!peset hold communicators as SHMEM-compatible triads (start, log2(stride), num)
  type, private :: communicator
     integer, pointer :: list(:)
     integer :: count
#ifdef use_libSMA
     integer :: start, log2stride
#elif use_libMPI
     integer :: id, group    !MPI communicator and group id for this PE set
#endif
  end type
  integer, parameter :: PESET_MAX=32 !should be .LE. max num of MPI communicators
  type(communicator) :: peset(0:PESET_MAX) !0 is a dummy used to hold single-PE "self" communicator
  integer :: peset_num=0, current_peset_num=0, peset_work_array_ptr=0
  integer, target :: peset_work_array(MAXPES*PESET_MAX) !used to hold peset%list pointers
  integer :: world_peset_num !the world communicator

!performance profiling
!  This profiles every type of MPI/SHMEM call within
!    a specified region of high-level code
!  Initialize or retrieve a clock with
!  id = mpp_clock_id( 'Region identifier name' )
!  Then set caliper points around the region using:
!  call mpp_clock_begin(id)
!  ...
!  call mpp_clock_end(id)
!  mpp_exit will print out the results.
#ifdef __sgi
  integer(LONG_KIND), private :: tick, ticks_per_sec, max_ticks, start_tick, end_tick, tick0=0
#else
  integer, private :: tick, ticks_per_sec, max_ticks, start_tick, end_tick, tick0=0
#endif
  real, private :: tick_rate
  integer, private, parameter :: MAX_CLOCKS=100, MAX_EVENT_TYPES=5, MAX_EVENTS=40000
!event types
  integer, private, parameter :: EVENT_ALLREDUCE=1, EVENT_BROADCAST=2, EVENT_RECV=3, EVENT_SEND=4, EVENT_WAIT=5
  integer, private :: clock_num=0, current_clock=0
  integer, private :: clock0    !measures total runtime from mpp_init to mpp_exit
!the event contains information for each type of event (e.g SHMEM_PUT)
  type, private :: event
     character(len=16)  :: name
     integer(LONG_KIND) :: ticks(MAX_EVENTS), bytes(MAX_EVENTS)
     integer            :: calls
  end type event
!a clock contains an array of event profiles for a region
  integer, parameter, public :: MPP_CLOCK_SYNC=1, MPP_CLOCK_DETAILED=2
  type, private :: clock
     character(len=32) :: name
     integer(LONG_KIND) :: tick, total_ticks
     logical :: sync_on_begin, detailed
     type(event), pointer :: events(:) !if needed, allocate to MAX_EVENT_TYPES
  end type
  type(clock) :: clocks(MAX_CLOCKS)

  integer,parameter :: MAX_BINS=20
  TYPE :: Clock_Data_Summary
    character(len=16) :: name
    real(DOUBLE_KIND) :: msg_size_sums(MAX_BINS)
    real(DOUBLE_KIND) :: msg_time_sums(MAX_BINS)
    real(DOUBLE_KIND) :: total_data
    real(DOUBLE_KIND) :: total_time
    integer(LONG_KIND) :: msg_size_cnts(MAX_BINS)
    integer(LONG_KIND) :: total_cnts
  END TYPE Clock_Data_Summary

  TYPE :: Summary_Struct
    character(len=16)         :: name
    type (Clock_Data_Summary) :: event(MAX_EVENT_TYPES)
  END TYPE Summary_Struct
  type(Summary_Struct) :: clock_summary(MAX_CLOCKS)
  
!public interfaces
  interface mpp_max
     module procedure mpp_max_real8
#ifndef no_8byte_integers
     module procedure mpp_max_int8
#endif
     module procedure mpp_max_real4
     module procedure mpp_max_int4
  end interface
  interface mpp_min
     module procedure mpp_min_real8
#ifndef no_8byte_integers
     module procedure mpp_min_int8
#endif
     module procedure mpp_min_real4
     module procedure mpp_min_int4
  end interface
  interface mpp_sum
#ifndef no_8byte_integers
     module procedure mpp_sum_int8
     module procedure mpp_sum_int8_scalar
#endif
     module procedure mpp_sum_real8
     module procedure mpp_sum_real8_scalar
     module procedure mpp_sum_cmplx8
     module procedure mpp_sum_cmplx8_scalar
     module procedure mpp_sum_int4
     module procedure mpp_sum_int4_scalar
     module procedure mpp_sum_real4
     module procedure mpp_sum_real4_scalar
     module procedure mpp_sum_cmplx4
     module procedure mpp_sum_cmplx4_scalar
  end interface
  interface mpp_transmit
     module procedure mpp_transmit_real8
     module procedure mpp_transmit_real8_scalar
     module procedure mpp_transmit_cmplx8
     module procedure mpp_transmit_cmplx8_scalar
#ifndef no_8byte_integers
     module procedure mpp_transmit_int8
     module procedure mpp_transmit_int8_scalar
     module procedure mpp_transmit_logical8
     module procedure mpp_transmit_logical8_scalar
#endif
     module procedure mpp_transmit_real4
     module procedure mpp_transmit_real4_scalar
     module procedure mpp_transmit_cmplx4
     module procedure mpp_transmit_cmplx4_scalar
     module procedure mpp_transmit_int4
     module procedure mpp_transmit_int4_scalar
     module procedure mpp_transmit_logical4
     module procedure mpp_transmit_logical4_scalar
  end interface
  interface mpp_recv
     module procedure mpp_recv_real8
     module procedure mpp_recv_real8_scalar
     module procedure mpp_recv_cmplx8
     module procedure mpp_recv_cmplx8_scalar
#ifndef no_8byte_integers
     module procedure mpp_recv_int8
     module procedure mpp_recv_int8_scalar
     module procedure mpp_recv_logical8
     module procedure mpp_recv_logical8_scalar
#endif
     module procedure mpp_recv_real4
     module procedure mpp_recv_real4_scalar
     module procedure mpp_recv_cmplx4
     module procedure mpp_recv_cmplx4_scalar
     module procedure mpp_recv_int4
     module procedure mpp_recv_int4_scalar
     module procedure mpp_recv_logical4
     module procedure mpp_recv_logical4_scalar
  end interface
  interface mpp_send
     module procedure mpp_send_real8
     module procedure mpp_send_real8_scalar
     module procedure mpp_send_cmplx8
     module procedure mpp_send_cmplx8_scalar
#ifndef no_8byte_integers
     module procedure mpp_send_int8
     module procedure mpp_send_int8_scalar
     module procedure mpp_send_logical8
     module procedure mpp_send_logical8_scalar
#endif
     module procedure mpp_send_real4
     module procedure mpp_send_real4_scalar
     module procedure mpp_send_cmplx4
     module procedure mpp_send_cmplx4_scalar
     module procedure mpp_send_int4
     module procedure mpp_send_int4_scalar
     module procedure mpp_send_logical4
     module procedure mpp_send_logical4_scalar
  end interface

  interface mpp_broadcast
     module procedure mpp_broadcast_real8
     module procedure mpp_broadcast_real8_scalar
     module procedure mpp_broadcast_cmplx8
     module procedure mpp_broadcast_cmplx8_scalar
#ifndef no_8byte_integers
     module procedure mpp_broadcast_int8
     module procedure mpp_broadcast_int8_scalar
     module procedure mpp_broadcast_logical8
     module procedure mpp_broadcast_logical8_scalar
#endif
     module procedure mpp_broadcast_real4
     module procedure mpp_broadcast_real4_scalar
     module procedure mpp_broadcast_cmplx4
     module procedure mpp_broadcast_cmplx4_scalar
     module procedure mpp_broadcast_int4
     module procedure mpp_broadcast_int4_scalar
     module procedure mpp_broadcast_logical4
     module procedure mpp_broadcast_logical4_scalar
  end interface

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
  end interface

  interface mpp_error
     module procedure mpp_error_basic
     module procedure mpp_error_mesg
     module procedure mpp_error_noargs
  end interface

#ifdef use_libSMA
!currently SMA contains no generic shmem_wait for different integer kinds:
!I have inserted one here
  interface shmem_integer_wait
     module procedure shmem_int4_wait_local
     module procedure shmem_int8_wait_local
  end interface
#endif
  public :: mpp_broadcast, mpp_chksum, mpp_clock_begin, mpp_clock_end, mpp_clock_id, mpp_declare_pelist, mpp_error, &
            mpp_error_state, mpp_exit, mpp_init, mpp_max, mpp_min, mpp_node, mpp_npes, mpp_pe, mpp_recv, mpp_root_pe, mpp_send, &
            mpp_set_root_pe, mpp_set_warn_level, mpp_sum, mpp_set_stack_size, mpp_sync, mpp_sync_self, mpp_transmit
  public :: stdin, stdout, stderr, stdlog
#ifdef use_shmalloc
  public :: mpp_malloc
#endif

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_init( flags, in, out, err, log )
      integer, optional, intent(in) :: flags, in, out, err, log
      integer :: my_pe, num_pes, len
      integer :: i
      logical :: opened
#ifdef _CRAYT3E
      intrinsic my_pe
#endif

      if( mpp_initialized )return

#ifdef use_libSMA
      call START_PES(0)         !the argument 0 means extract from environment variable NPES on PVP/SGI, from mpprun -n on t3e
      pe = my_pe()
      node = pe                 !on an SMP this should return node ID rather than PE number.
      npes = num_pes()
#elif use_libMPI
      call MPI_INITIALIZED( opened, error ) !in case called from another MPI package
      if( .NOT.opened )call MPI_INIT(error)
      call MPI_COMM_RANK( MPI_COMM_WORLD, pe,   error )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, npes, error )
      allocate( request(0:npes-1) )
      request(:) = MPI_REQUEST_NULL
#endif
      mpp_initialized = .TRUE.

!PEsets: make defaults illegal
      peset(:)%count = -1
#ifdef use_libSMA
      peset(:)%start = -1
      peset(:)%log2stride = -1
#elif use_libMPI
      peset(:)%id = -1
      peset(:)%group = -1
#endif
!0=single-PE, initialized so that count returns 1
      peset(0)%count = 1
#ifdef use_libMPI
      current_peset_num = 0
      peset(0)%id = MPI_COMM_WORLD
      call MPI_COMM_GROUP( MPI_COMM_WORLD, peset(0)%group, error )
#endif
      world_peset_num = get_peset( (/(i,i=0,npes-1)/) )
      current_peset_num = world_peset_num !initialize current PEset to world

!initialize clocks
      call SYSTEM_CLOCK( count=tick0, count_rate=ticks_per_sec, count_max=max_ticks )
      tick_rate = 1./ticks_per_sec
      clock0 = mpp_clock_id( 'Total runtime', flags=MPP_CLOCK_SYNC )

      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
      end if
!logunit: log messages are written to logfile.out by default
!if optional argument logunit=stdout, write messages to stdout instead.
!if specifying non-defaults, you must specify units not yet in use.
      if( PRESENT(in) )then
          inquire( unit=in, opened=opened )
          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stdin.' )
          in_unit=in
      end if
      if( PRESENT(out) )then
          inquire( unit=out, opened=opened )
          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stdout.' )
          out_unit=out
      end if
      if( PRESENT(err) )then
          inquire( unit=err, opened=opened )
          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stderr.' )
          err_unit=err
      end if
      if( PRESENT(log) )then
          inquire( unit=log, opened=opened )
          if( opened .AND. log.NE.out_unit )call mpp_error( FATAL, 'MPP_INIT: unable to open stdlog.' )
          log_unit=log
      end if
!log_unit can be written to only from root_pe, all others write to stdout
      if( pe.EQ.root_pe .AND. log_unit.NE.out_unit )then
          inquire( unit=log_unit, opened=opened )
          if( opened )call mpp_error( FATAL, 'MPP_INIT: specified unit for stdlog already in use.' )
          open( unit=log_unit, status='REPLACE', file='logfile.out' )
      end if
!messages
      if( verbose )call mpp_error( NOTE, 'MPP_INIT: initializing MPP module...' )
      if( pe.EQ.root_pe )then
          write( stdlog(),'(/a)' )'MPP module '//trim(version)
          write( stdlog(),'(a,i4)' )'MPP started with NPES=', npes
#ifdef use_libSMA
          write( stdlog(),'(a)' )'Using SMA (shmem) library for message passing...'
#endif
#ifdef use_libMPI
          write( stdlog(),'(a)' )'Using MPI library for message passing...'
#endif
          write( stdlog(), '(a,es12.4,a,i10,a)' ) &
               'Realtime clock resolution=', tick_rate, ' sec (', ticks_per_sec, ' ticks/sec)'
          write( stdlog(), '(a,es12.4,a,i20,a)' ) &
               'Clock rolls over after ', max_ticks*tick_rate, ' sec (', max_ticks, ' ticks)'
      end if

#ifdef use_libSMA
#ifdef use_shmalloc
!we use shpalloc to ensure all these are remotely accessible
      len=0; ptr_sync = LOC(pe)   !null initialization
      call mpp_malloc( ptr_sync,        size(TRANSFER(sync,word)),            len )
      len=0; ptr_status = LOC(pe)  !null initialization
      call mpp_malloc( ptr_status, npes*size(TRANSFER(status(0),word)),   len )
      len=0; ptr_remote = LOC(pe) !null initialization
      call mpp_malloc( ptr_remote, npes*size(TRANSFER(remote_data_loc(0),word)), len )
      len=0; ptr_from = LOC(pe)   !null initialization
      call mpp_malloc( ptr_from,        size(TRANSFER(mpp_from_pe,word)),     len )
#else
      allocate( status(0:npes-1) )
      allocate( remote_data_loc(0:npes-1) )
#endif
      sync(:) = SHMEM_SYNC_VALUE
      status(0:npes-1) = MPP_READY
      remote_data_loc(0:npes-1) = MPP_WAIT
      call mpp_set_stack_size(32768) !default initial value
#endif
      call mpp_clock_begin(clock0)

      return
    end subroutine mpp_init

    function stdin()
      integer :: stdin
      stdin = in_unit
      return
    end function stdin

    function stdout()
      integer :: stdout
      stdout = out_unit
      return
    end function stdout

    function stderr()
      integer :: stderr
      stderr = err_unit
      return
    end function stderr

    function stdlog()
      integer :: stdlog
      if( pe.EQ.root_pe )then
          stdlog = log_unit
      else
          stdlog = out_unit
      end if
      return
    end function stdlog

    subroutine mpp_exit()
!to be called at the end of a run
      integer :: i, j, k, n, nmax
      real :: t, tmin, tmax, tavg, tstd
      real :: m, mmin, mmax, mavg, mstd
      real :: t_total

      if( .NOT.mpp_initialized )return
      call mpp_clock_end(clock0)
      t_total = clocks(clock0)%total_ticks*tick_rate
      if( clock_num.GT.0 )then
          if( ANY(clocks(1:clock_num)%detailed) )then
              call sum_clock_data; call dump_clock_summary
          end if
          if( pe.EQ.root_pe )then
              write( stdout(),'(/a)' ) 'Tabulating mpp_clock statistics across PEs...'
              if( ANY(clocks(1:clock_num)%detailed) ) &
                   write( stdout(),'(a)' )'   ... see mpp_clock.out.#### for details on individual PEs.'
              write( stdout(),'(/32x,a)' ) '          tmin          tmax          tavg          tstd  tfrac'
          end if
          do i = 1,clock_num
!times between mpp_clock ticks
             t = clocks(i)%total_ticks*tick_rate
             tmin = t; call mpp_min(tmin)
             tmax = t; call mpp_max(tmax)
             tavg = t; call mpp_sum(tavg,1); tavg = tavg/mpp_npes()
             tstd = (t-tavg)**2; call mpp_sum(tstd,1); tstd = sqrt( tstd/mpp_npes() )
             if( pe.EQ.root_pe )write( stdout(),'(a32,4f14.6,f7.3)' ) &
                  clocks(i)%name, tmin, tmax, tavg, tstd, tavg/t_total
          end do
          if( ANY(clocks(1:clock_num)%detailed) .AND. pe.EQ.root_pe )write( stdout(),'(/32x,a)' ) &
               '       tmin       tmax       tavg       tstd       mmin       mmax       mavg       mstd  mavg/tavg'
          do i = 1,clock_num
!messages: bytelengths and times
	     if( .NOT.clocks(i)%detailed )cycle
             do j = 1,MAX_EVENT_TYPES
                n = clocks(i)%events(j)%calls; nmax = n
                call mpp_max(nmax)
                if( nmax.NE.0 )then
!don't divide by n because n might be 0
                    m = 0
                    if( n.GT.0 )m = sum(clocks(i)%events(j)%bytes(1:n))
                    mmin = m; call mpp_min(mmin)
                    mmax = m; call mpp_max(mmax)
                    mavg = m; call mpp_sum(mavg,1); mavg = mavg/mpp_npes()
                    mstd = (m-mavg)**2; call mpp_sum(mstd,1); mstd = sqrt( mstd/mpp_npes() )
                    t = 0
                    if( n.GT.0 )t = sum(clocks(i)%events(j)%ticks(1:n))*tick_rate
                    tmin = t; call mpp_min(tmin)
                    tmax = t; call mpp_max(tmax)
                    tavg = t; call mpp_sum(tavg,1); tavg = tavg/mpp_npes()
                    tstd = (t-tavg)**2; call mpp_sum(tstd,1); tstd = sqrt( tstd/mpp_npes() )
                    if( pe.EQ.root_pe )write( stdout(),'(a32,4f11.3,5es11.3)' ) &
                         trim(clocks(i)%name)//' '//trim(clocks(i)%events(j)%name), &
                         tmin, tmax, tavg, tstd, mmin, mmax, mavg, mstd, mavg/tavg
                end if
             end do
          end do
      end if
      call mpp_set_current_pelist()
      call mpp_sync()
      call mpp_max(mpp_stack_hwm)
      if( pe.EQ.root_pe )write( stdout(),* )'MPP_STACK high water mark=', mpp_stack_hwm
#ifdef use_libMPI
      call MPI_FINALIZE(error)
#endif

      return
    end subroutine mpp_exit

    function mpp_pe()
      integer :: mpp_pe

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_PE: You must first call mpp_init.' )
      mpp_pe = pe
      return
    end function mpp_pe

    function mpp_node()
      integer :: mpp_node

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_NODE: You must first call mpp_init.' )
      mpp_node = node
      return
    end function mpp_node

    function mpp_npes()
      integer :: mpp_npes

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_NPES: You must first call mpp_init.' )
!      mpp_npes = npes
      mpp_npes = size(peset(current_peset_num)%list)
      return
    end function mpp_npes

    function mpp_root_pe()
      integer :: mpp_root_pe

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_ROOT_PE: You must first call mpp_init.' )
      mpp_root_pe = root_pe
      return
    end function mpp_root_pe

    subroutine mpp_set_root_pe(num)
      integer, intent(in) :: num
      logical :: opened

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_SET_ROOT_PE: You must first call mpp_init.' )
      if( .NOT.(ANY(num.EQ.peset(current_peset_num)%list)) ) &
           call mpp_error( FATAL, 'MPP_SET_ROOT_PE: you cannot set a root PE outside the current pelist.' )
!actions to take if root_pe has changed:
! open log_unit on new root_pe, close it on old root_pe and point its log_unit to stdout.
      if( num.NE.root_pe )then  !root_pe has changed
          if( pe.EQ.num )then
!on the new root_pe
              if( log_unit.NE.out_unit )then
                  inquire( unit=log_unit, opened=opened )
                  if( .NOT.opened )open( unit=log_unit, status='OLD', file='logfile.out', position='APPEND' )
              end if
          else if( pe.EQ.root_pe )then
!on the old root_pe
              if( log_unit.NE.out_unit )then
                  inquire( unit=log_unit, opened=opened )
                  if( opened )close(log_unit)
                  log_unit = out_unit
              end if
          end if
      end if
      root_pe = num
      return
    end subroutine mpp_set_root_pe

    subroutine mpp_declare_pelist( pelist )
!this call is written specifically to accommodate a brain-dead MPI restriction
!that requires a parent communicator to create a child communicator:
!in other words: a pelist cannot go off and declare a communicator, but every PE
!in the parent, including those not in pelist(:), must get together for the
!MPI_COMM_CREATE call. The parent is typically MPI_COMM_WORLD, though it could also
!be a subset that includes all PEs in pelist.
!This restriction does not apply to SMA but to have uniform code,
!you may as well call it. It must be placed in a context where all PEs call it.
!Subsequent calls that use the pelist should be called PEs in the pelist only.
      integer, intent(in) :: pelist(:)
      integer :: i

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_DECLARE_PELIST: You must first call mpp_init.' )
      i = get_peset(pelist)
      return
    end subroutine mpp_declare_pelist

    subroutine mpp_set_current_pelist( pelist )
!Once we branch off into a PE subset, we want subsequent "global" calls to
!sync only across this subset. This is declared as the current pelist (peset(current_peset_num)%list)
!when current_peset all pelist ops with no pelist should apply the current pelist.
!also, we set the start PE in this pelist to be the root_pe.
!unlike mpp_declare_pelist, this is called by the PEs in the pelist only
!so if the PEset has not been previously declared, this will hang in MPI.
!if pelist is omitted, we reset pelist to the world pelist.
      integer, intent(in), optional :: pelist(:)
      integer :: i

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_SET_CURRENT_PELIST: You must first call mpp_init.' )
      if( PRESENT(pelist) )then
          if( .NOT.ANY(pe.EQ.pelist) )call mpp_error( FATAL, 'MPP_SET_CURRENT_PELIST: pe must be in pelist.' )
          current_peset_num = get_peset(pelist)
      else
          current_peset_num = world_peset_num
      end if
      call mpp_set_root_pe(peset(current_peset_num)%list(1))
      call mpp_sync()           !this is called to make sure everyone in the current pelist is here.
      return
    end subroutine mpp_set_current_pelist

    function get_peset(pelist)
      integer :: get_peset
!makes a PE set out of a PE list
!a PE list is an ordered list of PEs
!a PE set is a triad (start,log2stride,size) for SHMEM, an a communicator for MPI
!if stride is non-uniform or not a power of 2, will return error (not required for MPI but enforced for uniformity)
      integer, intent(in), optional :: pelist(:)
      integer :: group
      integer :: i, n, stride
      integer, allocatable :: sorted(:)

      if( .NOT.PRESENT(pelist) )then !set it to current_peset_num
          get_peset = current_peset_num; return
      end if
      if( size(pelist).EQ.1 .AND. npes.GT.1 )then    !collective ops on single PEs should return
          get_peset = 0; return
      end if
!make a sorted list
      n = 1
      if( ascend_sort(pelist).NE.1 )call mpp_error( FATAL, 'GET_PESET: sort error.' )   !result is the array sorted(:)
      if( debug )write( stderr(),* )'pelist=', pelist, ' sorted=', sorted
!find if this array matches any existing peset
      do i = 1,peset_num
         if( debug )write( stderr(),'(a,3i4)' )'pe, i, peset_num=', pe, i, peset_num
         if( size(sorted).EQ.size(peset(i)%list) )then
             if( ALL(sorted.EQ.peset(i)%list) )then
                 get_peset = i; return
             end if
         end if
      end do
!not found, so create new peset
      peset_num = peset_num + 1
      if( peset_num.GE.PESET_MAX )call mpp_error( FATAL, 'GET_PESET: number of PE sets exceeds PESET_MAX.' )
      i = peset_num             !shorthand
!create list
      peset(i)%list => peset_work_array(peset_work_array_ptr+1:peset_work_array_ptr+size(sorted))
      peset_work_array_ptr = peset_work_array_ptr + size(sorted)
      peset(i)%list(:) = sorted(:)
      peset(i)%count = size(sorted)
#ifdef use_libSMA
      peset(i)%start = sorted(1)
      if( size(sorted).GT.1 )then
          stride = sorted(2)-sorted(1)
          if( ANY(sorted(2:n)-sorted(1:n-1).NE.stride) ) &
               call mpp_error( WARNING, 'GET_PESET: pelist must have constant stride.' )
          peset(i)%log2stride = nint( log(real(stride))/log(2.) )
          if( 2**peset(i)%log2stride.NE.stride )call mpp_error( WARNING, 'GET_PESET: pelist must have power-of-2 stride.' )
      else
          peset(i)%log2stride = 0
      end if
#elif use_libMPI
      call MPI_GROUP_INCL( peset(current_peset_num)%group, size(sorted), sorted, peset(i)%group, error )
      call MPI_COMM_CREATE( peset(current_peset_num)%id, peset(i)%group, peset(i)%id, error )
#endif
      deallocate(sorted)
      get_peset = i

      return

      contains
        
        recursive function ascend_sort(a)
          integer :: ascend_sort
          integer, intent(in) :: a(:)
          integer :: b, i
          if( size(a).EQ.1 .OR. ALL(a.EQ.a(1)) )then
              allocate( sorted(n) )
              sorted(n) = a(1)
              ascend_sort = n
              return
          end if
          b = minval(a)
          n = n + 1
          i = ascend_sort( pack(a,mask=a.NE.b) )
          ascend_sort = i - 1
          sorted(i-1) = b
          return
        end function ascend_sort

    end function get_peset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                        PERFORMANCE PROFILING CALLS                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpp_clock_id( name, flags )
!return an ID for a new or existing clock
      integer :: mpp_clock_id
      character(len=*), intent(in) :: name
      integer, intent(in), optional :: flags
      integer :: i

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_CLOCK_ID: You must first call mpp_init.' )
      mpp_clock_id = 1
      if( clock_num.EQ.0 )then  !first
!         allocate( clocks(MAX_CLOCKS) )
          clock_num = mpp_clock_id
          call clock_init(mpp_clock_id,name,flags)
      else
          FIND_CLOCK: do while( trim(name).NE.trim(clocks(mpp_clock_id)%name) )
             mpp_clock_id = mpp_clock_id + 1
             if( mpp_clock_id.GT.clock_num )then
                 if( mpp_clock_id.GT.MAX_CLOCKS )then
                     call mpp_error( WARNING, 'MPP_CLOCK_ID: too many clock requests, this one is ignored.' )
                 else               !new clock: initialize
                     clock_num = mpp_clock_id
                     call clock_init(mpp_clock_id,name,flags)
                     exit FIND_CLOCK
                 end if
             end if
          end do FIND_CLOCK
      endif
      return
    end function mpp_clock_id

    subroutine mpp_clock_begin(id)
      integer, intent(in) :: id

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: You must first call mpp_init.' )
      if( id.EQ.0 )return
      if( id.LT.0 .OR. id.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid id.' )

      if( clocks(id)%sync_on_begin )then
!do an untimed sync at the beginning of the clock
!this puts all PEs in the current pelist on par, so that measurements begin together
!ending time will be different, thus measuring load imbalance for this clock.
          current_clock = 0; call mpp_sync()
      end if
      current_clock = id
      call SYSTEM_CLOCK( clocks(id)%tick )
      return
    end subroutine mpp_clock_begin

    subroutine mpp_clock_end(id)
!the id argument is not used for anything at present
      integer, intent(in), optional :: id
      integer(LONG_KIND) :: delta

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_CLOCK_END: You must first call mpp_init.' )
      if( id.EQ.0 )return
      if( id.LT.0 .OR. id.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid id.' )
      call SYSTEM_CLOCK(end_tick)
      delta = end_tick - clocks(id)%tick
      if( delta.LT.0 )then
          write( stderr(),* )'pe, id, start_tick, end_tick, delta, max_ticks=', pe, id, clocks(id)%tick, end_tick, delta, max_ticks
          delta = delta + max_ticks + 1
          call mpp_error( WARNING, 'MPP_CLOCK_END: Clock rollover, assumed single roll.' )
      end if
      clocks(id)%total_ticks = clocks(id)%total_ticks + delta
      current_clock = 0
      return
    end subroutine mpp_clock_end

    subroutine increment_current_clock( event_id, bytes )
      integer, intent(in) :: event_id
      integer, intent(in), optional :: bytes
      integer :: n
      integer(LONG_KIND) :: delta

      if( current_clock.EQ.0 )return
      if( current_clock.LT.0 .OR. current_clock.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid current_clock.' )
      if( .NOT.clocks(current_clock)%detailed )return
      call SYSTEM_CLOCK(end_tick)
      n = clocks(current_clock)%events(event_id)%calls + 1

      if( n.EQ.MAX_EVENTS )call mpp_error( WARNING, &
           'MPP_CLOCK: events exceed MAX_EVENTS, ignore detailed profiling data for clock '//trim(clocks(current_clock)%name) )
      if( n.GT.MAX_EVENTS )return

      clocks(current_clock)%events(event_id)%calls = n
      delta = end_tick - start_tick
      if( delta.LT.0 )then
          delta = delta + max_ticks + 1
          call mpp_error( WARNING, 'MPP_CLOCK_END: Clock rollover, assumed single roll.' )
      end if
      clocks(current_clock)%events(event_id)%ticks(n) = delta
      if( PRESENT(bytes) )clocks(current_clock)%events(event_id)%bytes(n) = bytes
      return
    end subroutine increment_current_clock

  subroutine dump_clock_summary()
    implicit none

    real :: total_time,total_time_all,total_data
    real :: msg_size,eff_BW,s
    integer :: SD_UNIT
    integer :: total_calls
    integer :: i,j,k,ct
    integer :: msg_cnt
    character(len=2)  :: u
    character(len=18) :: filename
    character(len=20),dimension(MAX_BINS),save :: bin

    data bin( 1)  /'  0   -    8    B:  '/
    data bin( 2)  /'  8   -   16    B:  '/
    data bin( 3)  /' 16   -   32    B:  '/
    data bin( 4)  /' 32   -   64    B:  '/
    data bin( 5)  /' 64   -  128    B:  '/
    data bin( 6)  /'128   -  256    B:  '/
    data bin( 7)  /'256   -  512    B:  '/
    data bin( 8)  /'512   - 1024    B:  '/
    data bin( 9)  /'  1.0 -    2.1 KB:  '/
    data bin(10)  /'  2.1 -    4.1 KB:  '/
    data bin(11)  /'  4.1 -    8.2 KB:  '/
    data bin(12)  /'  8.2 -   16.4 KB:  '/
    data bin(13)  /' 16.4 -   32.8 KB:  '/
    data bin(14)  /' 32.8 -   65.5 KB:  '/
    data bin(15)  /' 65.5 -  131.1 KB:  '/
    data bin(16)  /'131.1 -  262.1 KB:  '/
    data bin(17)  /'262.1 -  524.3 KB:  '/
    data bin(18)  /'524.3 - 1048.6 KB:  '/
    data bin(19)  /'  1.0 -    2.1 MB:  '/
    data bin(20)  /' >2.1          MB:  '/

    if( .NOT.ANY(clocks(1:clock_num)%detailed) )return
    write( filename,'(a,i4.4)' )'mpp_clock.out.', pe

    SD_UNIT = get_unit()
    open(SD_UNIT,file=trim(filename),form='formatted')

    COMM_TYPE: do ct = 1,clock_num

      if( .NOT.clocks(ct)%detailed )cycle
      write(SD_UNIT,*) &
          clock_summary(ct)%name(1:15),' Communication Data for PE ',pe

      write(SD_UNIT,*) ' '
      write(SD_UNIT,*) ' '

      total_time_all = 0.0
      EVENT_TYPE: do k = 1,MAX_EVENT_TYPES-1

        if(clock_summary(ct)%event(k)%total_time == 0.0)cycle

        total_time = clock_summary(ct)%event(k)%total_time
        total_time_all = total_time_all + total_time
        total_data = clock_summary(ct)%event(k)%total_data
        total_calls = clock_summary(ct)%event(k)%total_cnts

        write(SD_UNIT,1000) clock_summary(ct)%event(k)%name(1:9) // ':'

        write(SD_UNIT,1001) 'Total Data: ',total_data*1.0e-6, &
                            'MB; Total Time: ', total_time, &
                            'secs; Total Calls: ',total_calls

        write(SD_UNIT,*) ' '
        write(SD_UNIT,1002) '     Bin            Counts      Avg Size        Eff B/W'
        write(SD_UNIT,*) ' '

        BIN_LOOP: do j=1,MAX_BINS

          if(clock_summary(ct)%event(k)%msg_size_cnts(j)==0)cycle

          if(j<=8)then
            s = 1.0
            u = ' B'
          elseif(j<=18)then
            s = 1.0e-3
            u = 'KB'
          else
            s = 1.0e-6
            u = 'MB'
          endif

          msg_cnt = clock_summary(ct)%event(k)%msg_size_cnts(j)
          msg_size = &
            s*(clock_summary(ct)%event(k)%msg_size_sums(j)/real(msg_cnt))
          eff_BW = (1.0e-6)*( clock_summary(ct)%event(k)%msg_size_sums(j) / &
                                  clock_summary(ct)%event(k)%msg_time_sums(j) )

          write(SD_UNIT,1003) bin(j),msg_cnt,msg_size,u,eff_BW

        end do BIN_LOOP

        write(SD_UNIT,*) ' '
        write(SD_UNIT,*) ' '
      end do EVENT_TYPE

   ! "Data-less" WAIT

      if(clock_summary(ct)%event(MAX_EVENT_TYPES)%total_time>0.0)then

        total_time = clock_summary(ct)%event(MAX_EVENT_TYPES)%total_time
        total_time_all = total_time_all + total_time
        total_calls = clock_summary(ct)%event(MAX_EVENT_TYPES)%total_cnts

        write(SD_UNIT,1000) clock_summary(ct)%event(MAX_EVENT_TYPES)%name(1:9) // ':'

        write(SD_UNIT,1004) 'Total Calls: ',total_calls,'; Total Time: ', &
                             total_time,'secs'

      endif

      write(SD_UNIT,*) ' '
      write(SD_UNIT,1005) 'Total communication time spent for ' // &
                      clock_summary(ct)%name(1:9) // ': ',total_time_all,'secs'
      write(SD_UNIT,*) ' '
      write(SD_UNIT,*) ' '
      write(SD_UNIT,*) ' '

    end do COMM_TYPE

    close(SD_UNIT)

1000  format(a)
1001  format(a,f8.2,a,f8.2,a,i6)
1002  format(a)
1003  format(a,i6,'    ','  ',f6.1,a,'    ',f7.3,'MB/sec')
1004  format(a,i8,a,f9.2,a)
1005  format(a,f9.2,a)
    return
    contains
      integer function get_unit()
        implicit none

        integer,save :: i
        logical :: l_open

        i = 10
        do i=10,99
           inquire(unit=i,opened=l_open)
           if(.not.l_open)exit
        end do

        if(i==100)then
            call mpp_error(FATAL,'Unable to get I/O unit')
        else
            get_unit = i
        endif

        return
      end function get_unit

  end subroutine dump_clock_summary


  subroutine sum_clock_data()
    implicit none

    integer :: i,j,k,ct,event_size,event_cnt
    real    :: msg_time

    CLOCK_TYPE: do ct=1,clock_num
      if( .NOT.clocks(ct)%detailed )cycle
      EVENT_TYPE: do j=1,MAX_EVENT_TYPES-1
        event_cnt = clocks(ct)%events(j)%calls
        EVENT_SUMMARY: do i=1,event_cnt

        clock_summary(ct)%event(j)%total_cnts = &
              clock_summary(ct)%event(j)%total_cnts + 1

        event_size = clocks(ct)%events(j)%bytes(i)

        k = find_bin(event_size)

        clock_summary(ct)%event(j)%msg_size_cnts(k) = &
              clock_summary(ct)%event(j)%msg_size_cnts(k) + 1

        clock_summary(ct)%event(j)%msg_size_sums(k) = &
              clock_summary(ct)%event(j)%msg_size_sums(k) &
            + clocks(ct)%events(j)%bytes(i)

        clock_summary(ct)%event(j)%total_data = &
              clock_summary(ct)%event(j)%total_data &
            + clocks(ct)%events(j)%bytes(i)

        msg_time = clocks(ct)%events(j)%ticks(i)
        msg_time = tick_rate * real( clocks(ct)%events(j)%ticks(i) )

        clock_summary(ct)%event(j)%msg_time_sums(k) = &
              clock_summary(ct)%event(j)%msg_time_sums(k) + msg_time

        clock_summary(ct)%event(j)%total_time = &
              clock_summary(ct)%event(j)%total_time + msg_time

        end do EVENT_SUMMARY
      end do EVENT_TYPE

      j = MAX_EVENT_TYPES ! WAITs
           ! "msg_size_cnts" doesn't really mean anything for WAIT
           ! but position will be used to store number of counts for now.

      event_cnt = clocks(ct)%events(j)%calls
      clock_summary(ct)%event(j)%msg_size_cnts(1) = event_cnt
      clock_summary(ct)%event(j)%total_cnts       = event_cnt

      msg_time = tick_rate * real( sum ( clocks(ct)%events(j)%ticks(1:event_cnt) ) )
      clock_summary(ct)%event(j)%msg_time_sums(1) = &
              clock_summary(ct)%event(j)%msg_time_sums(1) + msg_time

      clock_summary(ct)%event(j)%total_time = clock_summary(ct)%event(j)%msg_time_sums(1)

    end do CLOCK_TYPE

    return
    contains
      integer function find_bin(event_size)
        implicit none

        integer,intent(in) :: event_size
        integer :: k,msg_size

        msg_size = 8
        k = 1
        do while(event_size>msg_size .and. k<MAX_BINS)
           k = k+1
           msg_size = msg_size*2
        end do
        find_bin = k
        return
      end function find_bin

  end subroutine sum_clock_data


  subroutine clock_init(id,name,flags)
    integer, intent(in) :: id
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: flags
    integer :: i

    clocks(id)%name = name
    clocks(id)%tick = 0
    clocks(id)%total_ticks = 0
    clocks(id)%sync_on_begin = .FALSE.
    clocks(id)%detailed      = .FALSE.
    if( PRESENT(flags) )then
        if( BTEST(flags,0) )clocks(id)%sync_on_begin = .TRUE.
        if( BTEST(flags,1) )clocks(id)%detailed      = .TRUE.
    end if
    if( clocks(id)%detailed )then
        allocate( clocks(id)%events(MAX_EVENT_TYPES) )
        clocks(id)%events(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
        clocks(id)%events(EVENT_BROADCAST)%name = 'BROADCAST'
        clocks(id)%events(EVENT_RECV)%name = 'RECV'
        clocks(id)%events(EVENT_SEND)%name = 'SEND'
        clocks(id)%events(EVENT_WAIT)%name = 'WAIT'
        do i=1,MAX_EVENT_TYPES
           clocks(id)%events(i)%ticks(:) = 0
           clocks(id)%events(i)%bytes(:) = 0
           clocks(id)%events(i)%calls = 0
        end do
        clock_summary(id)%name = name
        clock_summary(id)%event(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
        clock_summary(id)%event(EVENT_BROADCAST)%name = 'BROADCAST'
        clock_summary(id)%event(EVENT_RECV)%name = 'RECV'
        clock_summary(id)%event(EVENT_SEND)%name = 'SEND'
        clock_summary(id)%event(EVENT_WAIT)%name = 'WAIT'
        do i=1,MAX_EVENT_TYPES
           clock_summary(id)%event(i)%msg_size_sums(:) = 0.0
           clock_summary(id)%event(i)%msg_time_sums(:) = 0.0
           clock_summary(id)%event(i)%total_data = 0.0
           clock_summary(id)%event(i)%total_time = 0.0
           clock_summary(id)%event(i)%msg_size_cnts(:) = 0
           clock_summary(id)%event(i)%total_cnts = 0
        end do
    end if
    return
  end subroutine clock_init

#ifdef __sgi
    subroutine system_clock_sgi( count, count_rate, count_max )
!mimics F90 SYSTEM_CLOCK intrinsic
      integer(LONG_KIND), intent(out), optional :: count, count_rate, count_max
      integer(LONG_KIND) :: sgi_tick, sgi_ticks_per_sec, sgi_max_tick
!sgi_max_tick currently returns 64
!count must return a number between 0 and count_max
      integer(LONG_KIND), save :: maxtick=0
      if( maxtick.EQ.0 )then
          maxtick = sgi_max_tick() !actually reports #bits in maxtick
          if( maxtick.LT.BIT_SIZE(maxtick) )then
              maxtick = 2**maxtick
          else
              maxtick = huge(maxtick)
          end if
      end if
      if( PRESENT(count) )then
          count = modulo( sgi_tick()-tick0, maxtick )
!          count = sgi_tick()
      end if
      if( PRESENT(count_rate) )then
          count_rate = sgi_ticks_per_sec()
      end if
      if( PRESENT(count_max) )then
          count_max = maxtick-1
      end if
      return
    end subroutine system_clock_sgi
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                BASIC MESSAGE PASSING ROUTINE: mpp_transmit                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_TRANSMIT_ mpp_transmit_real8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_real8_scalar
#define MPP_RECV_ mpp_recv_real8
#define MPP_RECV_SCALAR_ mpp_recv_real8_scalar
#define MPP_SEND_ mpp_send_real8
#define MPP_SEND_SCALAR_ mpp_send_real8_scalar
#define MPP_BROADCAST_ mpp_broadcast_real8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_real8_scalar
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_REAL8
#define SHMEM_BROADCAST_ SHMEM_BROADCAST8
#define SHMEM_GET_ SHMEM_GET8
#include <mpp_transmit.h>

#define MPP_TRANSMIT_ mpp_transmit_real4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_real4_scalar
#define MPP_RECV_ mpp_recv_real4
#define MPP_RECV_SCALAR_ mpp_recv_real4_scalar
#define MPP_SEND_ mpp_send_real4
#define MPP_SEND_SCALAR_ mpp_send_real4_scalar
#define MPP_BROADCAST_ mpp_broadcast_real4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_real4_scalar
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_REAL4
#define SHMEM_BROADCAST_ SHMEM_BROADCAST4
#define SHMEM_GET_ SHMEM_GET4
#include <mpp_transmit.h>

#define MPP_TRANSMIT_ mpp_transmit_cmplx8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_cmplx8_scalar
#define MPP_RECV_ mpp_recv_cmplx8
#define MPP_RECV_SCALAR_ mpp_recv_cmplx8_scalar
#define MPP_SEND_ mpp_send_cmplx8
#define MPP_SEND_SCALAR_ mpp_send_cmplx8_scalar
#define MPP_BROADCAST_ mpp_broadcast_cmplx8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_cmplx8_scalar
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_TYPE_BYTELEN_ 16
#define MPI_TYPE_ MPI_DOUBLE_COMPLEX
#define SHMEM_BROADCAST_ SHMEM_BROADCAST8
#define SHMEM_GET_ SHMEM_GET128
#include <mpp_transmit.h>

#define MPP_TRANSMIT_ mpp_transmit_cmplx4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_cmplx4_scalar
#define MPP_RECV_ mpp_recv_cmplx4
#define MPP_RECV_SCALAR_ mpp_recv_cmplx4_scalar
#define MPP_SEND_ mpp_send_cmplx4
#define MPP_SEND_SCALAR_ mpp_send_cmplx4_scalar
#define MPP_BROADCAST_ mpp_broadcast_cmplx4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_cmplx4_scalar
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_COMPLEX
#define SHMEM_BROADCAST_ SHMEM_BROADCAST4
#define SHMEM_GET_ SHMEM_GET64
#include <mpp_transmit.h>

#ifndef no_8byte_integers
#define MPP_TRANSMIT_ mpp_transmit_int8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_int8_scalar
#define MPP_RECV_ mpp_recv_int8
#define MPP_RECV_SCALAR_ mpp_recv_int8_scalar
#define MPP_SEND_ mpp_send_int8
#define MPP_SEND_SCALAR_ mpp_send_int8_scalar
#define MPP_BROADCAST_ mpp_broadcast_int8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_int8_scalar
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_INTEGER8
#define SHMEM_BROADCAST_ SHMEM_BROADCAST8
#define SHMEM_GET_ SHMEM_GET8
#include <mpp_transmit.h>
#endif

#define MPP_TRANSMIT_ mpp_transmit_int4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_int4_scalar
#define MPP_RECV_ mpp_recv_int4
#define MPP_RECV_SCALAR_ mpp_recv_int4_scalar
#define MPP_SEND_ mpp_send_int4
#define MPP_SEND_SCALAR_ mpp_send_int4_scalar
#define MPP_BROADCAST_ mpp_broadcast_int4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_int4_scalar
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_INTEGER4
#define SHMEM_BROADCAST_ SHMEM_BROADCAST4
#define SHMEM_GET_ SHMEM_GET4
#include <mpp_transmit.h>

#ifndef no_8byte_integers
#define MPP_TRANSMIT_ mpp_transmit_logical8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_logical8_scalar
#define MPP_RECV_ mpp_recv_logical8
#define MPP_RECV_SCALAR_ mpp_recv_logical8_scalar
#define MPP_SEND_ mpp_send_logical8
#define MPP_SEND_SCALAR_ mpp_send_logical8_scalar
#define MPP_BROADCAST_ mpp_broadcast_logical8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_logical8_scalar
#define MPP_TYPE_ logical(LONG_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_INTEGER8
#define SHMEM_BROADCAST_ SHMEM_BROADCAST8
#define SHMEM_GET_ SHMEM_GET8
#include <mpp_transmit.h>
#endif

#define MPP_TRANSMIT_ mpp_transmit_logical4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_logical4_scalar
#define MPP_RECV_ mpp_recv_logical4
#define MPP_RECV_SCALAR_ mpp_recv_logical4_scalar
#define MPP_SEND_ mpp_send_logical4
#define MPP_SEND_SCALAR_ mpp_send_logical4_scalar
#define MPP_BROADCAST_ mpp_broadcast_logical4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_logical4_scalar
#define MPP_TYPE_ logical(INT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_INTEGER4
#define SHMEM_BROADCAST_ SHMEM_BROADCAST4
#define SHMEM_GET_ SHMEM_GET4
#include <mpp_transmit.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!            GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_REDUCE_ mpp_max_real8
#define MPP_TYPE_ real(DOUBLE_KIND)
#define SHMEM_REDUCE_ SHMEM_REAL8_MAX_TO_ALL
#define MPI_TYPE_ MPI_REAL8
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce.h>

#define MPP_REDUCE_ mpp_max_real4
#define MPP_TYPE_ real(FLOAT_KIND)
#define SHMEM_REDUCE_ SHMEM_REAL4_MAX_TO_ALL
#define MPI_TYPE_ MPI_REAL4
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce.h>

#ifndef no_8byte_integers   
#define MPP_REDUCE_ mpp_max_int8
#define MPP_TYPE_ integer(LONG_KIND)
#define SHMEM_REDUCE_ SHMEM_INT8_MAX_TO_ALL
#define MPI_TYPE_ MPI_INTEGER8
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce.h>
#endif

#define MPP_REDUCE_ mpp_max_int4
#define MPP_TYPE_ integer(INT_KIND)
#define SHMEM_REDUCE_ SHMEM_INT4_MAX_TO_ALL
#define MPI_TYPE_ MPI_INTEGER4
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce.h>

#define MPP_REDUCE_ mpp_min_real8
#define MPP_TYPE_ real(DOUBLE_KIND)
#define SHMEM_REDUCE_ SHMEM_REAL8_MIN_TO_ALL
#define MPI_TYPE_ MPI_REAL8
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce.h>

#define MPP_REDUCE_ mpp_min_real4
#define MPP_TYPE_ real(FLOAT_KIND)
#define SHMEM_REDUCE_ SHMEM_REAL4_MIN_TO_ALL
#define MPI_TYPE_ MPI_REAL4
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce.h>

#ifndef no_8byte_integers   
#define MPP_REDUCE_ mpp_min_int8
#define MPP_TYPE_ integer(LONG_KIND)
#define SHMEM_REDUCE_ SHMEM_INT8_MIN_TO_ALL
#define MPI_TYPE_ MPI_INTEGER8
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce.h>
#endif

#define MPP_REDUCE_ mpp_min_int4
#define MPP_TYPE_ integer(INT_KIND)
#define SHMEM_REDUCE_ SHMEM_INT4_MIN_TO_ALL
#define MPI_TYPE_ MPI_INTEGER4
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce.h>

#define MPP_SUM_ mpp_sum_real8
#define MPP_SUM_SCALAR_ mpp_sum_real8_scalar
#define MPP_TYPE_ real(DOUBLE_KIND)
#define SHMEM_SUM_ SHMEM_REAL8_SUM_TO_ALL
#define MPI_TYPE_ MPI_REAL8
#define MPP_TYPE_BYTELEN_ 8
#include <mpp_sum.h>

#define MPP_SUM_ mpp_sum_real4
#define MPP_SUM_SCALAR_ mpp_sum_real4_scalar
#define MPP_TYPE_ real(FLOAT_KIND)
#define SHMEM_SUM_ SHMEM_REAL4_SUM_TO_ALL
#define MPI_TYPE_ MPI_REAL4
#define MPP_TYPE_BYTELEN_ 4
#include <mpp_sum.h>

#define MPP_SUM_ mpp_sum_cmplx8
#define MPP_SUM_SCALAR_ mpp_sum_cmplx8_scalar
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define SHMEM_SUM_ SHMEM_COMP8_SUM_TO_ALL
#define MPI_TYPE_ MPI_DOUBLE_COMPLEX
#define MPP_TYPE_BYTELEN_ 16
#include <mpp_sum.h>

#define MPP_SUM_ mpp_sum_cmplx4
#define MPP_SUM_SCALAR_ mpp_sum_cmplx4_scalar
#define MPP_TYPE_ complex(FLOAT_KIND)
#define SHMEM_SUM_ SHMEM_COMP4_SUM_TO_ALL
#define MPI_TYPE_ MPI_COMPLEX
#define MPP_TYPE_BYTELEN_ 8
#include <mpp_sum.h>

#ifndef no_8byte_integers
#define MPP_SUM_ mpp_sum_int8
#define MPP_SUM_SCALAR_ mpp_sum_int8_scalar
#define MPP_TYPE_ integer(LONG_KIND)
#define SHMEM_SUM_ SHMEM_INT8_SUM_TO_ALL
#define MPI_TYPE_ MPI_INTEGER8
#define MPP_TYPE_BYTELEN_ 8
#include <mpp_sum.h>
#endif

#define MPP_SUM_ mpp_sum_int4
#define MPP_SUM_SCALAR_ mpp_sum_int4_scalar
#define MPP_TYPE_ integer(INT_KIND)
#define SHMEM_SUM_ SHMEM_INT4_SUM_TO_ALL
#define MPI_TYPE_ MPI_INTEGER4
#define MPP_TYPE_BYTELEN_ 4
#include <mpp_sum.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!           SYNCHRONIZATION ROUTINES: mpp_sync, mpp_sync_self                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_sync( pelist )
!synchronize PEs in list
      integer, intent(in), optional :: pelist(:)
      integer :: n

      call mpp_sync_self(pelist)

      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef use_libSMA
      if( n.EQ.world_peset_num )then
          call SHMEM_BARRIER_ALL() !special call is faster
      else
          call SHMEM_BARRIER( peset(n)%start, peset(n)%log2stride, peset(n)%count, sync )
      end if
#endif
#ifdef use_libMPI
      call MPI_BARRIER( peset(n)%id, error )
#endif
      if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)

      return
    end subroutine mpp_sync

    subroutine mpp_sync_self( pelist )
!this is to check if current PE's outstanding puts are complete
!but we can't use shmem_fence because we are actually waiting for
!a remote PE to complete its get
      integer, intent(in), optional :: pelist(:)
      integer :: i, m, n, stride

      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef use_libSMA
#ifdef _CRAYT90
      call SHMEM_UDCFLUSH !invalidate data cache
#endif
#endif
      do m = 1,peset(n)%count
         i = peset(n)%list(m)
#ifdef use_libSMA
         call SHMEM_INT8_WAIT( status(i), MPP_WAIT ) !wait for status.NE.MPP_WAIT
#endif
#ifdef use_libMPI
         if( request(i).NE.MPI_REQUEST_NULL )call MPI_WAIT( request(i), stat, error )
#endif
      end do
      if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
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
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!         MISCELLANEOUS UTILITIES: mpp_error, mpp_chksum, mpp_malloc          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_error_basic( errortype, errormsg )
!a very basic error handler
      integer, intent(in) :: errortype
      character(len=*), intent(in), optional :: errormsg
      character(len=128) :: text
      logical :: opened
      
      if( .NOT.mpp_initialized )call abort()
      if( errortype.EQ.NOTE    )text = 'NOTE'	!just FYI
      if( errortype.EQ.WARNING )text = 'WARNING'	!probable error
      if( errortype.EQ.FATAL   )text = 'FATAL'	!fatal error

      if( npes.GT.1 )write( text,'(a,i5)' )trim(text)//' from PE', pe	!this is the mpp part
      if( PRESENT(errormsg) )text = trim(text)//': '//trim(errormsg)

      if( errortype.EQ.NOTE )then
          write( stdout(),'(a)' )trim(text)
      else
          write( stderr(),'(/a/)' )trim(text)
          if( errortype.EQ.FATAL .OR. warnings_are_fatal )then
              call FLUSH(stdout())
#ifdef sgi_mipspro
              call TRACE_BACK_STACK_AND_PRINT()
#endif
#ifdef use_libMPI
              call MPI_ABORT( MPI_COMM_WORLD, 1, error )
#endif
              call ABORT()	!automatically calls traceback on Cray systems
          end if
      end if
      error_state = errortype
      return
    end subroutine mpp_error_basic
!overloads to mpp_error_basic
    subroutine mpp_error_mesg( routine, errormsg, errortype )
!support for error_mesg routine in FMS
      character(len=*), intent(in) :: routine, errormsg
      integer, intent(in) :: errortype
      call mpp_error( errortype, trim(routine)//': '//trim(errormsg) )
      return
    end subroutine mpp_error_mesg
    subroutine mpp_error_noargs()
      call mpp_error(FATAL)
    end subroutine mpp_error_noargs
      
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

    function mpp_error_state()
      integer :: mpp_error_state
      mpp_error_state = error_state
      return
    end function mpp_error_state

#ifdef use_shmalloc
    subroutine mpp_malloc( ptr, newlen, len )
!routine to perform symmetric allocation:
!this is required on the t3e/O2k for variables that will be non-local arguments
!to a shmem call (see man intro_shmem(3F)).
!newlen is the required allocation length for the pointer ptr
!   len is the current allocation (0 if unallocated)
      integer, intent(in) :: newlen
      integer, intent(inout) :: len
      real :: dummy
      integer :: words_per_long
      integer(LONG_KIND) :: long
!argument ptr is a cray pointer, points to a dummy argument in this routine
      pointer( ptr, dummy )
!      integer(LONG_KIND) :: error_8

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_MALLOC: You must first call mpp_init.' )
!use existing allocation if it is enough
      if( newlen.LE.len )return

      call SHMEM_BARRIER_ALL()
!if the pointer is already allocated, deallocate
!      if( len.NE.0 )call SHPDEALLC( ptr, error_8, -1 ) !BWA: error_8 instead of error, see PV 682618 (fixed in mpt.1.3.0.1)
      if( len.NE.0 )call SHPDEALLC( ptr, error, -1 )
!allocate new length: assume that the array is KIND=8
      words_per_long = size(transfer(long,word))
      call SHPALLOC( ptr, newlen*words_per_long, error, -1 )
      len = newlen
      call SHMEM_BARRIER_ALL()

      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( stdout(),'(a,i18,a,i5,a,2i8,i16)' )'T=', tick, ' PE=', pe, ' MPP_MALLOC: len, newlen, ptr=', len, newlen, ptr
      end if
      return
    end subroutine mpp_malloc
#endif use_shmalloc

    subroutine mpp_set_stack_size(n)
!set the mpp_stack variable to be at least n LONG words long
      integer, intent(in) :: n
      character(len=8) :: text
#ifdef use_shmalloc
      call mpp_malloc( ptr_stack, n, mpp_stack_size )
#else
      if( n.GT.mpp_stack_size .AND. allocated(mpp_stack) )deallocate(mpp_stack)
      if( .NOT.allocated(mpp_stack) )then
          allocate( mpp_stack(n) )
          mpp_stack_size = n
      end if
#endif
      write( text,'(i8)' )n
      if( pe.EQ.root_pe )call mpp_error( NOTE, 'MPP_SET_STACK_SIZE: stack size set to '//text//'.' )

      return
    end subroutine mpp_set_stack_size

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

  end module mpp_mod

#ifdef test_mpp
#ifdef SYSTEM_CLOCK
#undef SYSTEM_CLOCK
#endif
  program test
!test various aspects of mpp_mod
#ifdef sgi_mipspro
    use shmem_interface
#endif
    use mpp_mod
    implicit none
    integer :: pe, npes, root
    integer, parameter :: n=1048576
    real, allocatable, dimension(:) :: a, b, c
    integer :: tick, tick0, ticks_per_sec, id
    integer :: i, j, k, l, m
    real :: dt

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
!time transmit, compare against shmem_put and get
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
!mpp_transmit
       call mpp_sync()
       call SYSTEM_CLOCK(tick0)
       do i = 1,npes
          call mpp_transmit( a, l, modulo(pe+npes-i,npes), b, l, modulo(pe+i,npes) )
!          call mpp_sync_self( (/modulo(pe+npes-i,npes)/) )
       end do
       call mpp_sync()
       call SYSTEM_CLOCK(tick)
       dt = real(tick-tick0)/(npes*ticks_per_sec)
       dt = max( dt, epsilon(dt) )
       if( pe.EQ.root )write( stdout(),'(/a,i8,f13.6,f8.2)' )'MPP_TRANSMIT length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
#ifdef SGICRAY
!shmem_put
       call mpp_sync()
       call SYSTEM_CLOCK(tick0)
       do i = 1,npes
          call shmem_put8( b, a, l, modulo(pe+1,npes) )
       end do
       call mpp_sync()
       call SYSTEM_CLOCK(tick)
       dt = real(tick-tick0)/(npes*ticks_per_sec)
       dt = max( dt, epsilon(dt) )
       if( pe.EQ.root )write( stdout(),'( a,i8,f13.6,f8.2)' )'SHMEM_PUT    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!shmem_get
       call mpp_sync()
       call SYSTEM_CLOCK(tick0)
       do i = 1,npes
          call shmem_get8( b, a, l, modulo(pe+1,npes) )
       end do
       call SYSTEM_CLOCK(tick)
       dt = real(tick-tick0)/(npes*ticks_per_sec)
       dt = max( dt, epsilon(dt) )
       if( pe.EQ.root )write( stdout(),'( a,i8,f13.6,f8.2)' )'SHMEM_GET    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
#endif
       l = l/2
    end do

!test mpp_sum
    if( pe.EQ.root )then
        print '(/a)', 'Time mpp_sum...'
    end if
    a = real(pe+1)
    call mpp_sync()
    call SYSTEM_CLOCK(tick0)
    call mpp_sum(a,n)
    call SYSTEM_CLOCK(tick)
    dt = real(tick-tick0)/ticks_per_sec
    dt = max( dt, epsilon(dt) )
    if( pe.EQ.root )write( stdout(),'(a,2i4,f9.1,i8,f13.6,f8.2/)' ) &
         'mpp_sum: pe, npes, sum(pe+1), length, time, bw(Mb/s)=', pe, npes, a(1), n, dt, n*8e-6/dt
    call mpp_clock_end(id)

!test mpp_max
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
    call flush(stdout())
    if( npes.GE.2 )then
        if( pe.EQ.root )print *, 'Test of pelists: bcast, sum and max using PEs 0...npes-2 (excluding last PE)'
        call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
            
        a = real(pe+1)
        if( pe.NE.npes-1 )call mpp_broadcast( a, n, npes-2, (/(i,i=0,npes-2)/) )
        print *, 'bcast(npes-1) from 0 to npes-2=', pe, a(1)
        a = real(pe+1)
        if( pe.NE.npes-1 )call mpp_sum( a, n, (/(i,i=0,npes-2)/) )
        if( pe.EQ.root )print *, 'sum(pe+1) from 0 to npes-2=', a(1)
        a = real(pe+1)
        if( pe.NE.npes-1 )call mpp_max( a(1), (/(i,i=0,npes-2)/) )
        if( pe.EQ.root )print *, 'max(pe+1) from 0 to npes-2=', a(1)
    end if
#ifdef use_CRI_pointers
!test mpp_chksum
    if( modulo(n,npes).EQ.0 )then  !only set up for even division
        if( pe.EQ.root )call random_number(a)
        call mpp_sync()
        call mpp_transmit( a, n, ALL_PES, a, n, root )
        m= n/npes
        allocate( c(m) )
        c = a(pe*m+1:pe*m+m)

        if( pe.EQ.root )then
            print *
            print *, 'Test mpp_chksum...'
            print *, 'This test shows that a whole array and a distributed array give identical checksums.'
        end if
        print *, 'chksum(a)=', mpp_chksum(a,(/pe/))
        print *, 'chksum(c)=', mpp_chksum(c)
    end if
#endif
    call mpp_exit()
  end program test
#endif test_mpp
