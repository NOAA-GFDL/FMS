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
       '$Id: mpp.F90,v 6.0 2001/03/06 20:26:41 fms Exp $'
  character(len=128), private :: name= &
       '$Name: damascus $'

!various lengths (see shpalloc) are estimated in "words" which are 32bit on SGI, 64bit on Cray
!these are also the expected sizeof of args to MPI/shmem libraries
#ifdef _CRAY
  integer(LONG_KIND), private :: word(1)
#endif
#ifdef sgi_mipspro
  integer(INT_KIND), private :: word(1)
#endif
  integer, parameter :: word_kind=KIND(word)

#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
  integer, parameter, private :: stdin=100, stdout=101, stderr=102
#else
  integer, parameter, private :: stdin=5, stdout=6, stderr=0
#endif
  logical, private :: mpp_initialized=.FALSE.
  integer, private :: pe=0, node=0, npes=1, root_pe=0
  integer, private :: error

!initialization flags
  integer, parameter, public :: MPP_VERBOSE=1, MPP_DEBUG=2
  logical, private :: verbose=.FALSE., debug=.FALSE.

!flags to transmit routines
  integer, parameter, public :: ALL_PES=-1, ANY_PE=-2, NULL_PE=-3

!errortype flags
  integer, parameter, public :: NOTE=0, WARNING=1, FATAL=2
  logical, private :: warnings_are_fatal = .FALSE.

!timing
!since these are mainly timers associated with communication
!we only do real times, not CPU times
  integer, private :: tick, ticks_per_sec

  integer(LONG_KIND), parameter, private :: MPP_WAIT=-1, MPP_READY=-2
#ifdef use_libSMA
#include <mpp/shmem.fh>
  integer :: sync(SHMEM_REDUCE_SYNC_SIZE+SHMEM_BCAST_SYNC_SIZE+SHMEM_BARRIER_SYNC_SIZE)
!status and remote_data_loc are used to synchronize communication is MPP_TRANSMIT
#ifdef use_shmalloc
!9999 has no particular significance: the actual length is set by mpp_init
!but this value will prevent out-of-bound warnings until NPES>9999.
  integer(LONG_KIND), private, dimension(0:9999) :: status, remote_data_loc
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
  integer, private :: tag=1, stat(MPI_STATUS_SIZE), group_all
  integer, private, allocatable :: request(:)
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
  integer, private :: mpp_stack_size=0

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
  integer, private :: ticks_per_sec, max_ticks, start_tick, end_tick
  real, private :: tick_rate
  integer, private, parameter :: MAX_CLOCKS=10, MAX_EVENTS=6, MAX_EVENT_CNT=40000
  integer, private, parameter :: EVENT_ALLREDUCE=1, EVENT_BROADCAST=2, EVENT_RECV=3, EVENT_REDUCE=4, EVENT_SEND=5, EVENT_WAIT=6
  integer, private :: clock_num=0, current_clock=0
!the event contains information for each type of event (e.g SHMEM_PUT)
  type, private :: event
     character(len=16)                :: name
     integer,dimension(MAX_EVENT_CNT) :: ticks, bytes
     integer                          :: calls
  end type event
!a clock contains an array of event profiles for a region
  type, private :: clock
     character(len=32) :: name
     type(event)       :: events(MAX_EVENTS)
  end type
  type(clock), allocatable :: clocks(:)

  integer,parameter :: MAX_BINS=20
  TYPE :: Clock_Data_Summary
    character(len=16)                      :: name
    real(DOUBLE_KIND), dimension(MAX_BINS) :: msg_size_sums
    real(DOUBLE_KIND), dimension(MAX_BINS) :: msg_time_sums
    real(DOUBLE_KIND)                      :: total_data
    real(DOUBLE_KIND)                      :: total_time
    integer, dimension(MAX_BINS)           :: msg_size_cnts
    integer                                :: total_cnts
  END TYPE Clock_Data_Summary

  TYPE :: Summary_Struct
    character(len=16)                               :: name
    type (Clock_Data_Summary),dimension(MAX_EVENTS) :: event
  END TYPE Summary_Struct
  type(Summary_Struct), dimension(MAX_CLOCKS) :: clock_summary
  
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

#ifdef use_libSMA
!currently SMA contains no generic shmem_wait for different integer kinds:
!I have inserted one here
  interface shmem_integer_wait
     module procedure shmem_int4_wait_local
     module procedure shmem_int8_wait_local
  end interface
#endif
  public :: mpp_chksum, mpp_clock_begin, mpp_clock_end, mpp_clock_id, mpp_error, mpp_exit, mpp_init, mpp_max, mpp_min, &
            mpp_node, mpp_npes, mpp_pe, mpp_recv, mpp_root_pe, mpp_send, mpp_set_root_pe, mpp_set_warn_level, mpp_sum, &
            mpp_set_stack_size, mpp_sync, mpp_sync_self, mpp_transmit
#ifdef use_shmalloc
  public :: mpp_malloc
#endif

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_init(flags)
      integer, optional, intent(in) :: flags
      integer :: my_pe, num_pes, len
#ifdef _CRAYT3E
      intrinsic my_pe
#endif

      if( mpp_initialized )return

#ifdef use_libSMA
      call START_PES(0)         !the argument 0 means extract from environment variable NPES on PVP/SGI, from mpprun -n on t3e
      pe = my_pe()
      node = pe                 !on an SMP this should return node ID rather than PE number.
      npes = num_pes()
#endif use_libSMA
#ifdef use_libMPI
      call MPI_INIT(error)
      call MPI_COMM_RANK ( MPI_COMM_WORLD, pe,        error )
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, npes,      error )
      call MPI_COMM_GROUP( MPI_COMM_WORLD, group_all, error )
      allocate( request(0:npes-1) )
      request(:) = MPI_REQUEST_NULL
#endif
      mpp_initialized = .TRUE.
!initialize clocks
      call system_clock( count_rate=ticks_per_sec, count_max=max_ticks )
      tick_rate = 1./ticks_per_sec

      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
      end if
!messages
      if( verbose )call mpp_error( NOTE, 'MPP_INIT: initializing MPP module...' )
      if( pe.EQ.root_pe )then
          write( stdout,'(/a)' )'MPP module '//trim(version)
          write( stdout,* )'NPES=', npes
#ifdef use_libSMA
          write( stdout,* )'Using SMA (shmem) library for message passing...'
#endif
#ifdef use_libMPI
          write( stdout,* )'Using MPI library for message passing...'
#endif
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
      
      call mpp_sync()

      return
    end subroutine mpp_init

    subroutine mpp_exit()
!to be called at the end of a run
!not strictly required for shmem_runs (but required for MPI)
      real :: total_time, time_per_call, bandwidth
      integer :: i, j

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_EXIT: You must first call mpp_init.' )
      call mpp_sync()
      if( clock_num.GT.0 )then
          call sum_clock_data()
          call dump_clock_summary()
      end if
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

    function mpp_node()
      integer :: mpp_node
      mpp_node = node
      return
    end function mpp_node

    function mpp_npes()
      integer :: mpp_npes
      mpp_npes = npes
      return
    end function mpp_npes

    function mpp_root_pe()
      integer :: mpp_root_pe
      mpp_root_pe = root_pe
      return
    end function mpp_root_pe

    subroutine mpp_set_root_pe(num)
      integer, intent(in) :: num
      root_pe = num
      return
    end subroutine mpp_set_root_pe

    subroutine make_pe_set(pelist,peset)
!makes a PE set out of a PE list (list length must be .GE.2)
!a PE list is an ordered list of PEs
!a PE set is a triad (start,log2stride,size) for SHMEM, an a communicator for MPI (other two elements are unused)
!if stride is non-uniform or not a power of 2, SHMEM version will return error
      integer, intent(in) :: pelist(0:)
      integer, intent(out) :: peset(3)
      integer :: group
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
      call MPI_GROUP_INCL( group_all, size(pelist), pelist, group, error )
      call MPI_COMM_CREATE( MPI_COMM_WORLD, group, peset(1), error )
#endif

      return
    end subroutine make_pe_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                        PERFORMANCE PROFILING CALLS                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpp_clock_id(name)
!return an ID for a new or existing clock
      integer :: mpp_clock_id
      character(len=32), intent(in) :: name
      integer :: i

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_CLOCK_ID: You must first call mpp_init.' )
      mpp_clock_id = 1
      if( clock_num.EQ.0 )then  !first
          allocate( clocks(MAX_CLOCKS) )
          clock_num = mpp_clock_id
          call clock_init(mpp_clock_id,name)
      else
          FIND_CLOCK: do while( name.NE.clocks(mpp_clock_id)%name )
             mpp_clock_id = mpp_clock_id + 1
             if( mpp_clock_id.GT.clock_num )then
                 if( clock_num.GT.max_clocks )then
                     call mpp_error( WARNING, 'MPP_CLOCK_ID: too many clock requests, this one is ignored.' )
                 else               !new clock: initialize
                     clock_num = mpp_clock_id
                     call clock_init(mpp_clock_id,name)
                     exit FIND_CLOCK
                 end if
             end if
          end do FIND_CLOCK
      endif
      return
    end function mpp_clock_id

    subroutine mpp_clock_begin(id)
      integer, intent(in) :: id

      if( id.LE.0 .OR. id.GT.clock_num )return
      current_clock = id
      return
    end subroutine mpp_clock_begin

    subroutine mpp_clock_end(id)
!the id argument is not used for anything at present
      integer, intent(in), optional :: id
      current_clock = 0
    end subroutine mpp_clock_end

    subroutine increment_current_clock( event_id, bytes )
      integer, intent(in) :: event_id
      integer, intent(in), optional :: bytes
      integer :: n

      call SYSTEM_CLOCK(end_tick)
      n = clocks(current_clock)%events(event_id)%calls + 1

      if(n>MAX_EVENT_CNT)call mpp_error(FATAL,'Profiling event count exceeds MAX_EVENT_CNT')

      clocks(current_clock)%events(event_id)%calls = n
      clocks(current_clock)%events(event_id)%ticks(n) = &
                      clocks(current_clock)%events(event_id)%ticks(n) &
                    + end_tick - start_tick + 1
      if( PRESENT(bytes) )clocks(current_clock)%events(event_id)%bytes(n) = &
                    clocks(current_clock)%events(event_id)%bytes(n) + bytes
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
    character(len=3)  :: f_ext
    character(len=15) :: filename
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

    write(f_ext,'(i3)') pe
    i = index(f_ext,' ',BACK=.true.)
    filename = 'comm_summary.' // f_ext(i+1:3)

    SD_UNIT = get_unit()
    open(SD_UNIT,file=trim(filename),form='formatted')

    COMM_TYPE: do ct = 1,clock_num

      write(SD_UNIT,*) &
          clock_summary(ct)%name(1:15),' Communication Data for PE ',pe

      write(SD_UNIT,*) ' '
      write(SD_UNIT,*) ' '

      total_time_all = 0.0
      EVENT_TYPE: do k = 1,MAX_EVENTS-1

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
            s*(clock_summary(ct)%event(k)%msg_size_sums(j)/float(msg_cnt))
          eff_BW = (1.0e-6)*( clock_summary(ct)%event(k)%msg_size_sums(j) / &
                                  clock_summary(ct)%event(k)%msg_time_sums(j) )

          write(SD_UNIT,1003) bin(j),msg_cnt,msg_size,u,eff_BW

        end do BIN_LOOP

        write(SD_UNIT,*) ' '
        write(SD_UNIT,*) ' '
      end do EVENT_TYPE

   ! "Data-less" WAIT

      if(clock_summary(ct)%event(MAX_EVENTS)%total_time>0.0)then

        total_time = clock_summary(ct)%event(MAX_EVENTS)%total_time
        total_time_all = total_time_all + total_time
        total_calls = clock_summary(ct)%event(MAX_EVENTS)%total_cnts

        write(SD_UNIT,1000) clock_summary(ct)%event(MAX_EVENTS)%name(1:9) // ':'

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
  end subroutine dump_clock_summary


  subroutine sum_clock_data()
    implicit none

    integer :: i,j,k,ct,event_size,event_cnt
    real    :: msg_time

    CLOCK_TYPE: do ct=1,clock_num
      EVENT_TYPE: do j=1,MAX_EVENTS-1
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
        msg_time = tick_rate * float( clocks(ct)%events(j)%ticks(i) )

        clock_summary(ct)%event(j)%msg_time_sums(k) = &
              clock_summary(ct)%event(j)%msg_time_sums(k) + msg_time

        clock_summary(ct)%event(j)%total_time = &
              clock_summary(ct)%event(j)%total_time + msg_time

        end do EVENT_SUMMARY
      end do EVENT_TYPE

      j = MAX_EVENTS ! WAITs
           ! "msg_size_cnts" doesn't really mean anything for WAIT
           ! but position will be used to store number of counts for now.

      event_cnt = clocks(ct)%events(j)%calls
      clock_summary(ct)%event(j)%msg_size_cnts(1) = event_cnt
      clock_summary(ct)%event(j)%total_cnts       = event_cnt

      msg_time = tick_rate * float( sum ( clocks(ct)%events(j)%ticks(1:event_cnt) ) )
      clock_summary(ct)%event(j)%msg_time_sums(1) = &
              clock_summary(ct)%event(j)%msg_time_sums(1) + msg_time

      clock_summary(ct)%event(j)%total_time = clock_summary(ct)%event(j)%msg_time_sums(1)

    end do CLOCK_TYPE

    return
  end subroutine sum_clock_data


  subroutine clock_init(id,name)
    integer, intent(in) :: id
    character(len=24), intent(in) :: name
    integer :: i

    clocks(id)%name = name
    clocks(id)%events(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
    clocks(id)%events(EVENT_BROADCAST)%name = 'BROADCAST'
    clocks(id)%events(EVENT_RECV)%name = 'RECV'
    clocks(id)%events(EVENT_REDUCE)%name = 'REDUCE'
    clocks(id)%events(EVENT_SEND)%name = 'SEND'
    clocks(id)%events(EVENT_WAIT)%name = 'WAIT'
    do i=1,MAX_EVENTS
      clocks(id)%events(i)%ticks(:) = 0
      clocks(id)%events(i)%bytes(:) = 0
      clocks(id)%events(i)%calls = 0
    end do
    clock_summary(id)%name = name
    clock_summary(id)%event(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
    clock_summary(id)%event(EVENT_BROADCAST)%name = 'BROADCAST'
    clock_summary(id)%event(EVENT_RECV)%name = 'RECV'
    clock_summary(id)%event(EVENT_REDUCE)%name = 'REDUCE'
    clock_summary(id)%event(EVENT_SEND)%name = 'SEND'
    clock_summary(id)%event(EVENT_WAIT)%name = 'WAIT'
    do i=1,MAX_EVENTS
      clock_summary(id)%event(i)%msg_size_sums(:) = 0.0
      clock_summary(id)%event(i)%msg_time_sums(:) = 0.0
      clock_summary(id)%event(i)%total_data = 0.0
      clock_summary(id)%event(i)%total_time = 0.0
      clock_summary(id)%event(i)%msg_size_cnts(:) = 0
      clock_summary(id)%event(i)%total_cnts = 0
    end do
    return
  end subroutine clock_init


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
      integer :: peset(3)

      call mpp_sync_self(pelist)
      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
      if( PRESENT(pelist) )then
          call make_pe_set(pelist,peset)
#ifdef use_libSMA
          call SHMEM_BARRIER( peset(1), peset(2), peset(3), sync )
#endif
#ifdef use_libMPI
          call MPI_BARRIER( peset(1), error )
#endif
      else
#ifdef use_libSMA
          call SHMEM_BARRIER_ALL()
#endif
#ifdef use_libMPI
          call MPI_BARRIER( MPI_COMM_WORLD, error )
#endif
      end if
      if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)

      return
    end subroutine mpp_sync

    subroutine mpp_sync_self( pelist )
!this is to check if current PE's outstanding puts are complete
!but we can't use shmem_fence because we are actually waiting for
!a remote PE to complete its get
      integer, intent(in), optional :: pelist(:)
      integer :: i
      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef use_libSMA
#ifdef _CRAYT90
      call SHMEM_UDCFLUSH !invalidate data cache
#endif
      if( PRESENT(pelist) )then
          do i = pelist(1),pelist(size(pelist))
             call SHMEM_INT8_WAIT( status(i), MPP_WAIT ) !wait for status.NE.MPP_WAIT
          end do
      else
          do i = 0,npes-1
             call SHMEM_INT8_WAIT( status(i), MPP_WAIT ) !wait for status.NE.MPP_WAIT
          end do
      end if
#endif use_libSMA
#ifdef use_libMPI
      if( PRESENT(pelist) )then
!check if outstanding non-blocking messages to pelist are complete
          do i = pelist(1),pelist(size(pelist))
             if( request(i).NE.MPI_REQUEST_NULL )call MPI_WAIT( request(i), stat, error )
          end do
      else
!check if outstanding non-blocking messages to all PEs are complete
          do i = 0,npes-1
             if( request(i).NE.MPI_REQUEST_NULL )call MPI_WAIT( request(i), stat, error )
          end do
      end if
#endif
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

      if( npes.GT.1 )write( text,'(a,i5)' )trim(text)//' from PE', pe	!this is the mpp part
      if( PRESENT(errormsg) )text = trim(text)//': '//trim(errormsg)

      if( errortype.EQ.NOTE )then
          write( stdout,'(/a/)' )trim(text)
      else
          write( stderr,'(/a/)' )trim(text)
          if( errortype.EQ.FATAL .OR. warnings_are_fatal )then
              call FLUSH(stdout)
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
          write( stdout,'(a,i18,a,i5,a,2i8,i16)' )'T=', tick, ' PE=', pe, ' MPP_MALLOC: len, newlen, ptr=', len, newlen, ptr
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
  program test
!test various aspects of mpp_mod
#ifdef sgi_mipspro
    use shmem_interface
#endif
    use mpp_mod
    integer :: pe, npes
#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
    integer, parameter :: stdin=100, stdout=101, stderr=102
#else
    integer, parameter :: stdin=5, stdout=6, stderr=0
#endif
    integer, parameter :: n=1048576
    real, allocatable, dimension(:) :: a, b, c
    integer :: tick, tick0, ticks_per_sec, id

    call mpp_init()
    call mpp_set_stack_size(3145746)
    pe = mpp_pe()
    npes = mpp_npes()

    call system_clock( count_rate=ticks_per_sec )
    allocate( a(n), b(n) )
    id = mpp_clock_id( 'Random number' )
    call mpp_clock_begin(id)
    call random_number(a)
    call mpp_clock_end  (id)
!time transmit, compare against shmem_put and get
    if( pe.EQ.mpp_root_pe() )then
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
       call mpp_sync()
       call system_clock(tick0)
       do nn = 1,npes
          call mpp_transmit( a, l, mod(pe+npes-nn,npes), b, l, mod(pe+nn,npes) )
!          call mpp_sync_self( (/mod(pe+npes-nn,npes)/) )
       end do
       call mpp_sync()
       call system_clock(tick)
       dt = float(tick-tick0)/(npes*ticks_per_sec)
       if( pe.EQ.mpp_root_pe() )write( stdout,'(/a,i8,f13.6,f8.2)' )'MPP_TRANSMIT length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
#ifdef SGICRAY
!shmem_put
       call mpp_sync()
       call system_clock(tick0)
       do nn = 1,npes
          call shmem_put8( b, a, l, mod(pe+1,npes) )
       end do
       call mpp_sync()
       call system_clock(tick)
       dt = float(tick-tick0)/(npes*ticks_per_sec)
       if( pe.EQ.mpp_root_pe() )write( stdout,'( a,i8,f13.6,f8.2)' )'SHMEM_PUT    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!shmem_get
       call mpp_sync()
       call system_clock(tick0)
       do nn = 1,npes
          call shmem_get8( b, a, l, mod(pe+1,npes) )
       end do
       call system_clock(tick)
       dt = float(tick-tick0)/(npes*ticks_per_sec)
       if( pe.EQ.mpp_root_pe() )write( stdout,'( a,i8,f13.6,f8.2)' )'SHMEM_GET    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
#endif
       l = l/2
    end do

!test mpp_sum
    if( pe.EQ.mpp_root_pe() )then
        print *
        print *, 'Time mpp_sum...'
#ifdef bit_reproducible
        print *, 'The bit reproducibility flag for mpp_sum is turned on...'
#endif
    end if
    a = float(pe+1)
    call mpp_sync()
    call system_clock(tick0)
    call mpp_sum(a,n)
    call system_clock(tick)
    dt = float(tick-tick0)/ticks_per_sec
    if( pe.EQ.mpp_root_pe() )write( stdout,'(a,2i4,f9.1,i8,f13.6,f8.2/)' ) &
         'mpp_sum: pe, npes, sum(pe+1), length, time, bw(Mb/s)=', pe, npes, a(1), n, dt, n*8e-6/dt

!test mpp_max
    if( pe.EQ.mpp_root_pe() )then
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
        if( pe.EQ.mpp_root_pe() )call random_number(a)
        call mpp_sync()
        call mpp_transmit( a, n, ALL_PES, a, n, mpp_root_pe() )
        m= n/npes
        allocate( c(m) )
        c = a(pe*m+1:pe*m+m)

        if( pe.EQ.mpp_root_pe() )then
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
