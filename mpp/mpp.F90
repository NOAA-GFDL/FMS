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
!-----------------------------------------------------------------------
!                 Communication for message-passing codes
!
! AUTHOR: V. Balaji (V.Balaji@noaa.gov)
!         SGI/GFDL Princeton University
!
!-----------------------------------------------------------------------

!> @defgroup mpp_mod mpp_mod
!> @ingroup mpp
!> @brief This module defines interfaces for common operations using message-passing libraries.
!!
!> @author V. Balaji <"V.Balaji@noaa.gov">
!!
!!   A set of simple calls to provide a uniform interface
!!   to different message-passing libraries. It currently can be
!!   implemented either in the SGI/Cray native SHMEM library or in the MPI
!!   standard. Other libraries (e.g MPI-2, Co-Array Fortran) can be
!!   incorporated as the need arises.
!!
!!   The data transfer between a processor and its own memory is based
!!   on <TT>load</TT> and <TT>store</TT> operations upon
!!   memory. Shared-memory systems (including distributed shared memory
!!   systems) have a single address space and any processor can acquire any
!!   data within the memory by <TT>load</TT> and
!!   <TT>store</TT>. The situation is different for distributed
!!   parallel systems. Specialized MPP systems such as the T3E can simulate
!!   shared-memory by direct data acquisition from remote memory. But if
!!   the parallel code is distributed across a cluster, or across the Net,
!!   messages must be sent and received using the protocols for
!!   long-distance communication, such as TCP/IP. This requires a
!!   ``handshaking'' between nodes of the distributed system. One can think
!!   of the two different methods as involving <TT>put</TT>s or
!!   <TT>get</TT>s (e.g the SHMEM library), or in the case of
!!   negotiated communication (e.g MPI), <TT>send</TT>s and
!!   <TT>recv</TT>s.
!!
!!   The difference between SHMEM and MPI is that SHMEM uses one-sided
!!   communication, which can have very low-latency high-bandwidth
!!   implementations on tightly coupled systems. MPI is a standard
!!   developed for distributed computing across loosely-coupled systems,
!!   and therefore incurs a software penalty for negotiating the
!!   communication. It is however an open industry standard whereas SHMEM
!!   is a proprietary interface. Besides, the <TT>put</TT>s or
!!   <TT>get</TT>s on which it is based cannot currently be implemented in
!!   a cluster environment (there are recent announcements from Compaq that
!!   occasion hope).
!!
!!   The message-passing requirements of climate and weather codes can be
!!   reduced to a fairly simple minimal set, which is easily implemented in
!!   any message-passing API. <TT>mpp_mod</TT> provides this API.
!!
!!    Features of <TT>mpp_mod</TT> include:
!!    <ol>
!!     <li> Simple, minimal API, with free access to underlying API for </li>
!!       more complicated stuff.<BR/>
!!     <li> Design toward typical use in climate/weather CFD codes. </li>
!!     <li> Performance to be not significantly lower than any native API. </li>
!!    </ol>
!!
!!   This module is used to develop higher-level calls for <LINK
!!   SRC="mpp_domains.html">domain decomposition</LINK> and <LINK
!!   SRC="mpp_io.html">parallel I/O</LINK>.
!! <br/>
!!   Parallel computing is initially daunting, but it soon becomes
!!   second nature, much the way many of us can now write vector code
!!   without much effort. The key insight required while reading and
!!   writing parallel code is in arriving at a mental grasp of several
!!   independent parallel execution streams through the same code (the SPMD
!!   model). Each variable you examine may have different values for each
!!   stream, the processor ID being an obvious example. Subroutines and
!!   function calls are particularly subtle, since it is not always obvious
!!   from looking at a call what synchronization between execution streams
!!   it implies. An example of erroneous code would be a global barrier
!!   call (see <LINK SRC="#mpp_sync">mpp_sync</LINK> below) placed
!!   within a code block that not all PEs will execute, e.g:
!!
!!   <PRE>
!!   if( pe.EQ.0 )call mpp_sync()
!!   </PRE>
!!
!!   Here only PE 0 reaches the barrier, where it will wait
!!   indefinitely. While this is a particularly egregious example to
!!   illustrate the coding flaw, more subtle versions of the same are
!!   among the most common errors in parallel code.
!!  <br/>
!!   It is therefore important to be conscious of the context of a
!!   subroutine or function call, and the implied synchronization. There
!!   are certain calls here (e.g <TT>mpp_declare_pelist, mpp_init,
!!   mpp_set_stack_size</TT>) which must be called by all
!!   PEs. There are others which must be called by a subset of PEs (here
!!   called a <TT>pelist</TT>) which must be called by all the PEs in the
!!   <TT>pelist</TT> (e.g <TT>mpp_max, mpp_sum, mpp_sync</TT>). Still
!!   others imply no synchronization at all. I will make every effort to
!!   highlight the context of each call in the MPP modules, so that the
!!   implicit synchronization is spelt out.
!! <br/>
!!   For performance it is necessary to keep synchronization as limited
!!   as the algorithm being implemented will allow. For instance, a single
!!   message between two PEs should only imply synchronization across the
!!   PEs in question. A <I>global</I> synchronization (or <I>barrier</I>)
!!   is likely to be slow, and is best avoided. But codes first
!!   parallelized on a Cray T3E tend to have many global syncs, as very
!!   fast barriers were implemented there in hardware.
!! <br/>
!!   Another reason to use pelists is to run a single program in MPMD
!!   mode, where different PE subsets work on different portions of the
!!   code. A typical example is to assign an ocean model and atmosphere
!!   model to different PE subsets, and couple them concurrently instead of
!!   running them serially. The MPP module provides the notion of a
!!   <I>current pelist</I>, which is set when a group of PEs branch off
!!   into a subset. Subsequent calls that omit the <TT>pelist</TT> optional
!!   argument (seen below in many of the individual calls) assume that the
!!   implied synchronization is across the current pelist. The calls
!!   <TT>mpp_root_pe</TT> and <TT>mpp_npes</TT> also return the values
!!   appropriate to the current pelist. The <TT>mpp_set_current_pelist</TT>
!!   call is provided to set the current pelist.
!! </DESCRIPTION>
!! <br/>
!!
!!  F90 is a strictly-typed language, and the syntax pass of the
!!  compiler requires matching of type, kind and rank (TKR). Most calls
!!  listed here use a generic type, shown here as <TT>MPP_TYPE_</TT>. This
!!  is resolved in the pre-processor stage to any of a variety of
!!  types. In general the MPP operations work on 4-byte and 8-byte
!!  variants of <TT>integer, real, complex, logical</TT> variables, of
!!  rank 0 to 5, leading to 48 specific module procedures under the same
!!  generic interface. Any of the variables below shown as
!!  <TT>MPP_TYPE_</TT> is treated in this way.

!> @file
!> @brief File for @ref mpp_mod

module mpp_mod

! Define rank(X) for PGI compiler
#if defined( __PGI) || defined (__FLANG)
#define rank(X) size(shape(X))
#endif


#if defined(use_libMPI)
  use mpi
#endif

  use iso_fortran_env,   only : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
  use mpp_parameter_mod, only : MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE
  use mpp_parameter_mod, only : NOTE, WARNING, FATAL, MPP_CLOCK_DETAILED,MPP_CLOCK_SYNC
  use mpp_parameter_mod, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER
  use mpp_parameter_mod, only : CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
  use mpp_parameter_mod, only : MAX_EVENTS, MAX_BINS, MAX_EVENT_TYPES, MAX_CLOCKS
  use mpp_parameter_mod, only : MAXPES, EVENT_WAIT, EVENT_ALLREDUCE, EVENT_BROADCAST
  use mpp_parameter_mod, only : EVENT_ALLTOALL
  use mpp_parameter_mod, only : EVENT_TYPE_CREATE, EVENT_TYPE_FREE
  use mpp_parameter_mod, only : EVENT_RECV, EVENT_SEND, MPP_READY, MPP_WAIT
  use mpp_parameter_mod, only : mpp_parameter_version=>version
  use mpp_parameter_mod, only : DEFAULT_TAG
  use mpp_parameter_mod, only : COMM_TAG_1,  COMM_TAG_2,  COMM_TAG_3,  COMM_TAG_4
  use mpp_parameter_mod, only : COMM_TAG_5,  COMM_TAG_6,  COMM_TAG_7,  COMM_TAG_8
  use mpp_parameter_mod, only : COMM_TAG_9,  COMM_TAG_10, COMM_TAG_11, COMM_TAG_12
  use mpp_parameter_mod, only : COMM_TAG_13, COMM_TAG_14, COMM_TAG_15, COMM_TAG_16
  use mpp_parameter_mod, only : COMM_TAG_17, COMM_TAG_18, COMM_TAG_19, COMM_TAG_20
  use mpp_parameter_mod, only : MPP_FILL_INT,MPP_FILL_DOUBLE
  use mpp_data_mod,      only : stat, mpp_stack, ptr_stack, status, ptr_status, sync, ptr_sync
  use mpp_data_mod,      only : mpp_from_pe, ptr_from, remote_data_loc, ptr_remote
  use mpp_data_mod,      only : mpp_data_version=>version
  use platform_mod

implicit none
private

  !--- public paramters  -----------------------------------------------
  public :: MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE, NOTE, WARNING, FATAL
  public :: MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
  public :: CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
  public :: MAXPES, EVENT_RECV, EVENT_SEND
  public :: COMM_TAG_1,  COMM_TAG_2,  COMM_TAG_3,  COMM_TAG_4
  public :: COMM_TAG_5,  COMM_TAG_6,  COMM_TAG_7,  COMM_TAG_8
  public :: COMM_TAG_9,  COMM_TAG_10, COMM_TAG_11, COMM_TAG_12
  public :: COMM_TAG_13, COMM_TAG_14, COMM_TAG_15, COMM_TAG_16
  public :: COMM_TAG_17, COMM_TAG_18, COMM_TAG_19, COMM_TAG_20
  public :: MPP_FILL_INT,MPP_FILL_DOUBLE
  public :: mpp_init_test_full_init, mpp_init_test_init_true_only, mpp_init_test_peset_allocated
  public :: mpp_init_test_clocks_init, mpp_init_test_datatype_list_init, mpp_init_test_logfile_init
  public :: mpp_init_test_read_namelist, mpp_init_test_etc_unit, mpp_init_test_requests_allocated

  !--- public data from mpp_data_mod ------------------------------
!  public :: request

  !--- public interface from mpp_util.h ------------------------------
  public :: stdin, stdout, stderr, stdlog, lowercase, uppercase, mpp_error, mpp_error_state
  public :: mpp_set_warn_level, mpp_sync, mpp_sync_self, mpp_pe
  public :: mpp_npes, mpp_root_pe, mpp_set_root_pe, mpp_declare_pelist
  public :: mpp_get_current_pelist, mpp_set_current_pelist, mpp_get_current_pelist_name
  public :: mpp_clock_id, mpp_clock_set_grain, mpp_record_timing_data, get_unit
  public :: read_ascii_file, read_input_nml, mpp_clock_begin, mpp_clock_end
  public :: get_ascii_file_num_lines, get_ascii_file_num_lines_and_length
  public :: mpp_record_time_start, mpp_record_time_end

  !--- public interface from mpp_comm.h ------------------------------
  public :: mpp_chksum, mpp_max, mpp_min, mpp_sum, mpp_transmit, mpp_send, mpp_recv
  public :: mpp_sum_ad
  public :: mpp_broadcast, mpp_init, mpp_exit
  public :: mpp_gather, mpp_scatter, mpp_alltoall
  public :: mpp_type, mpp_byte, mpp_type_create, mpp_type_free

  !*********************************************************************
  !
  !    public data type
  !
  !*********************************************************************
  !> Communication information for message passing libraries
  !!
  !> peset hold communicators as SHMEM-compatible triads (start, log2(stride), num)
  !> @ingroup mpp_mod
  type :: communicator
     private
     character(len=32) :: name
     integer, pointer  :: list(:) =>NULL()
     integer           :: count
     integer           :: start, log2stride !< dummy variables when libMPI is defined.
     integer           :: id, group         !< MPI communicator and group id for this PE set.
  end type communicator

  !> Communication event profile
  !> @ingroup mpp_mod
  type :: event
     private
     character(len=16)                         :: name
     integer(i8_kind), dimension(MAX_EVENTS)   :: ticks, bytes
     integer                                   :: calls
  end type event

  !> a clock contains an array of event profiles for a region
  !> @ingroup mpp_mod
  type :: clock
     private
     character(len=32)    :: name
     integer(i8_kind)     :: hits
     integer(i8_kind)     :: tick
     integer(i8_kind)     :: total_ticks
     integer              :: peset_num
     logical              :: sync_on_begin, detailed
     integer              :: grain
     type(event), pointer :: events(:) =>NULL() !> if needed, allocate to MAX_EVENT_TYPES
     logical              :: is_on              !> initialize to false. set true when calling mpp_clock_begin
                                                !! set false when calling mpp_clock_end
  end type clock

  !> Summary of information from a clock run
  !> @ingroup mpp_mod
  type :: Clock_Data_Summary
     private
     character(len=16)  :: name
     real(r8_kind)      :: msg_size_sums(MAX_BINS)
     real(r8_kind)      :: msg_time_sums(MAX_BINS)
     real(r8_kind)      :: total_data
     real(r8_kind)      :: total_time
     integer(i8_kind)   :: msg_size_cnts(MAX_BINS)
     integer(i8_kind)   :: total_cnts
  end type Clock_Data_Summary

  !> holds name and clock data for use in @ref mpp_util.h
  !> @ingroup mpp_mod
  type :: Summary_Struct
     private
     character(len=16)         :: name
     type (Clock_Data_Summary) :: event(MAX_EVENT_TYPES)
  end type Summary_Struct

  !> Data types for generalized data transfer (e.g. MPI_Type)
  !> @ingroup mpp_mod
  type :: mpp_type
     private
     integer :: counter !> Number of instances of this type
     integer :: ndims
     integer, allocatable :: sizes(:)
     integer, allocatable :: subsizes(:)
     integer, allocatable :: starts(:)
     integer :: etype   !> Elementary data type (e.g. MPI_BYTE)
     integer :: id      !> Identifier within message passing library (e.g. MPI)

     type(mpp_type), pointer :: prev => null()
     type(mpp_type), pointer :: next => null()
  end type mpp_type

  !> Persisent elements for linked list interaction
  !> @ingroup mpp_mod
  type :: mpp_type_list
      private
      type(mpp_type), pointer :: head => null()
      type(mpp_type), pointer :: tail => null()
      integer :: length
  end type mpp_type_list

!***********************************************************************
!
!     public interface from mpp_util.h
!
!***********************************************************************
  !> @brief Error handler.
  !!
  !>    It is strongly recommended that all error exits pass through
  !!    <TT>mpp_error</TT> to assure the program fails cleanly. An individual
  !!    PE encountering a <TT>STOP</TT> statement, for instance, can cause the
  !!    program to hang. The use of the <TT>STOP</TT> statement is strongly
  !!    discouraged.
  !!
  !!    Calling mpp_error with no arguments produces an immediate error
  !!    exit, i.e:
  !!    <PRE>
  !!                    call mpp_error
  !!                    call mpp_error()
  !!    </PRE>
  !!    are equivalent.
  !!
  !!    The argument order
  !!    <PRE>
  !!                    call mpp_error( routine, errormsg, errortype )
  !!    </PRE>
  !!    is also provided to support legacy code. In this version of the
  !!    call, none of the arguments may be omitted.
  !!
  !!    The behaviour of <TT>mpp_error</TT> for a <TT>WARNING</TT> can be
  !!    controlled with an additional call <TT>mpp_set_warn_level</TT>.
  !!    <PRE>
  !!                    call mpp_set_warn_level(ERROR)
  !!    </PRE>
  !!    causes <TT>mpp_error</TT> to treat <TT>WARNING</TT>
  !!    exactly like <TT>FATAL</TT>.
  !!    <PRE>
  !!                    call mpp_set_warn_level(WARNING)
  !!    </PRE>
  !!    resets to the default behaviour described above.
  !!
  !!    <TT>mpp_error</TT> also has an internal error state which
  !!    maintains knowledge of whether a warning has been issued. This can be
  !!    used at startup in a subroutine that checks if the model has been
  !!    properly configured. You can generate a series of warnings using
  !!    <TT>mpp_error</TT>, and then check at the end if any warnings has been
  !!    issued using the function <TT>mpp_error_state()</TT>. If the value of
  !!    this is <TT>WARNING</TT>, at least one warning has been issued, and
  !!    the user can take appropriate action:
  !!
  !!    <PRE>
  !!                    if( ... )call mpp_error( WARNING, '...' )
  !!                    if( ... )call mpp_error( WARNING, '...' )
  !!                    if( ... )call mpp_error( WARNING, '...' )
  !!                    ...
  !!                    if( mpp_error_state().EQ.WARNING )call mpp_error( FATAL, '...' )
  !!    </PRE>
  !!  </DESCRIPTION>
  !! <br> Example usage:
  !! @code{.F90}
  !! call mpp_error( errortype, routine, errormsg )
  !! @endcode
  !! @param errortype
  !!    One of <TT>NOTE</TT>, <TT>WARNING</TT> or <TT>FATAL</TT>
  !!    (these definitions are acquired by use association).
  !!    <TT>NOTE</TT> writes <TT>errormsg</TT> to <TT>STDOUT</TT>.
  !!    <TT>WARNING</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>.
  !!    <TT>FATAL</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>,
  !!    and induces a clean error exit with a call stack traceback.
  !! @param routine Calling routine name
  !! @param errmsg Message to output
  !!  </IN>
  !> @ingroup mpp_mod
  interface mpp_error
     module procedure mpp_error_basic
     module procedure mpp_error_mesg
     module procedure mpp_error_noargs
     module procedure mpp_error_is
     module procedure mpp_error_rs
     module procedure mpp_error_ia
     module procedure mpp_error_ra
     module procedure mpp_error_ia_ia
     module procedure mpp_error_ia_ra
     module procedure mpp_error_ra_ia
     module procedure mpp_error_ra_ra
     module procedure mpp_error_ia_is
     module procedure mpp_error_ia_rs
     module procedure mpp_error_ra_is
     module procedure mpp_error_ra_rs
     module procedure mpp_error_is_ia
     module procedure mpp_error_is_ra
     module procedure mpp_error_rs_ia
     module procedure mpp_error_rs_ra
     module procedure mpp_error_is_is
     module procedure mpp_error_is_rs
     module procedure mpp_error_rs_is
     module procedure mpp_error_rs_rs
  end interface
  !> Takes a given integer or real array and returns it as a string
  !> @ingroup mpp_mod
  interface array_to_char
     module procedure iarray_to_char
     module procedure rarray_to_char
  end interface

!***********************************************************************
!
!    public interface from mpp_comm.h
!
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @fn subroutine mpp_init( flags, localcomm, test_level)
  !> @brief Initialize @ref mpp_mod
  !!
  !> Called to initialize the <TT>mpp_mod</TT> package. It is recommended
  !! that this call be the first executed line in your program. It sets the
  !! number of PEs assigned to this run (acquired from the command line, or
  !! through the environment variable <TT>NPES</TT>), and associates an ID
  !! number to each PE. These can be accessed by calling @ref mpp_npes and
  !! @ref mpp_pe.
  !! <br> Example usage:
  !!
  !!            call mpp_init( flags )
  !!
  !! @param flags
  !!   <TT>flags</TT> can be set to <TT>MPP_VERBOSE</TT> to
  !!   have <TT>mpp_mod</TT> keep you informed of what it's up to.
  !! @param test_level
  !!   Debugging flag to set amount of initialization tasks performed
  !> @ingroup mpp_mod

  !> @fn subroutine mpp_exit()
  !> @brief Exit <TT>@ref mpp_mod</TT>.
  !!
  !> Called at the end of the run, or to re-initialize <TT>mpp_mod</TT>,
  !! should you require that for some odd reason.
  !!
  !! This call implies synchronization across all PEs.
  !!
  !! <br>Example usage:
  !!
  !!            call mpp_exit()
  !> @ingroup mpp_mod

  !#####################################################################

  !> @fn subroutine mpp_set_stack_size(n)
  !> @brief Allocate module internal workspace.
  !> @param Integer to set stack size to(in words)
  !> <TT>mpp_mod</TT> maintains a private internal array called
  !! <TT>mpp_stack</TT> for private workspace. This call sets the length,
  !! in words, of this array.
  !!
  !! The <TT>mpp_init</TT> call sets this
  !! workspace length to a default of 32768, and this call may be used if a
  !! longer workspace is needed.
  !!
  !! This call implies synchronization across all PEs.
  !!
  !! This workspace is symmetrically allocated, as required for
  !! efficient communication on SGI and Cray MPP systems. Since symmetric
  !! allocation must be performed by <I>all</I> PEs in a job, this call
  !! must also be called by all PEs, using the same value of
  !! <TT>n</TT>. Calling <TT>mpp_set_stack_size</TT> from a subset of PEs,
  !! or with unequal argument <TT>n</TT>, may cause the program to hang.
  !!
  !! If any MPP call using <TT>mpp_stack</TT> overflows the declared
  !! stack array, the program will abort with a message specifying the
  !! stack length that is required. Many users wonder why, if the required
  !! stack length can be computed, it cannot also be specified at that
  !! point. This cannot be automated because there is no way for the
  !! program to know if all PEs are present at that call, and with equal
  !! values of <TT>n</TT>. The program must be rerun by the user with the
  !! correct argument to <TT>mpp_set_stack_size</TT>, called at an
  !! appropriate point in the code where all PEs are known to be present.
  !!        @verbose call mpp_set_stack_size(n)
  !!
  !> @ingroup mpp_mod
  public :: mpp_set_stack_size
  ! from mpp_util.h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              DATA TRANSFER TYPES: mpp_type_create                           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Create a mpp_type variable
  !> @param[in] field A field of any numerical or logical type
  !> @param[in] array_of_subsizes Integer array of subsizes
  !> @param[in] array_of_starts Integer array of starts
  !> @param[out] dtype_out Output variable for created @ref mpp_type
  !> @ingroup mpp_mod
  interface mpp_type_create
      module procedure mpp_type_create_int4
      module procedure mpp_type_create_int8
      module procedure mpp_type_create_real4
      module procedure mpp_type_create_real8
      module procedure mpp_type_create_cmplx4
      module procedure mpp_type_create_cmplx8
      module procedure mpp_type_create_logical4
      module procedure mpp_type_create_logical8
  end interface mpp_type_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !            GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min             !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Reduction operations.
  !>    Find the max of scalar a the PEs in pelist
  !!    result is also automatically broadcast to all PEs
  !!    @code{.F90}
  !!            call  mpp_max( a, pelist )
  !!    @endcode
  !> @param a <TT>real</TT> or <TT>integer</TT>, of 4-byte of 8-byte kind.
  !> @param pelist If <TT>pelist</TT> is omitted, the context is assumed to be the
  !!    current pelist. This call implies synchronization across the PEs in
  !!    <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !> @ingroup mpp_mod
  interface mpp_max
     module procedure mpp_max_real8_0d
     module procedure mpp_max_real8_1d
     module procedure mpp_max_int8_0d
     module procedure mpp_max_int8_1d
     module procedure mpp_max_real4_0d
     module procedure mpp_max_real4_1d
     module procedure mpp_max_int4_0d
     module procedure mpp_max_int4_1d
  end interface

  !> @brief Get minimum value out of the PEs in pelist
  !> Result is also broadcast to all PEs
  !> @ingroup mpp_mod
  interface mpp_min
     module procedure mpp_min_real8_0d
     module procedure mpp_min_real8_1d
     module procedure mpp_min_int8_0d
     module procedure mpp_min_int8_1d
     module procedure mpp_min_real4_0d
     module procedure mpp_min_real4_1d
     module procedure mpp_min_int4_0d
     module procedure mpp_min_int4_1d
  end interface


  !> @brief Reduction operation.
  !!
  !> <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !! <TT>integer, real, complex</TT> variables, of rank 0 or 1. A
  !! contiguous block from a multi-dimensional array may be passed by its
  !! starting address and its length, as in <TT>f77</TT>.
  !!
  !! Library reduction operators are not required or guaranteed to be
  !! bit-reproducible. In any case, changing the processor count changes
  !! the data layout, and thus very likely the order of operations. For
  !! bit-reproducible sums of distributed arrays, consider using the
  !! <TT>mpp_global_sum</TT> routine provided by the <LINK
  !! SRC="mpp_domains.html"><TT>mpp_domains</TT></LINK> module.
  !!
  !! The <TT>bit_reproducible</TT> flag provided in earlier versions of
  !! this routine has been removed.
  !!
  !!
  !! If <TT>pelist</TT> is omitted, the context is assumed to be the
  !! current pelist. This call implies synchronization across the PEs in
  !! <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !! Example usage:
  !!            call mpp_sum( a, length, pelist )
  !!
  !> @ingroup mpp_mod
  interface mpp_sum
     module procedure mpp_sum_int8
     module procedure mpp_sum_int8_scalar
     module procedure mpp_sum_int8_2d
     module procedure mpp_sum_int8_3d
     module procedure mpp_sum_int8_4d
     module procedure mpp_sum_int8_5d
     module procedure mpp_sum_real8
     module procedure mpp_sum_real8_scalar
     module procedure mpp_sum_real8_2d
     module procedure mpp_sum_real8_3d
     module procedure mpp_sum_real8_4d
     module procedure mpp_sum_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_sum_cmplx8
     module procedure mpp_sum_cmplx8_scalar
     module procedure mpp_sum_cmplx8_2d
     module procedure mpp_sum_cmplx8_3d
     module procedure mpp_sum_cmplx8_4d
     module procedure mpp_sum_cmplx8_5d
#endif
     module procedure mpp_sum_int4
     module procedure mpp_sum_int4_scalar
     module procedure mpp_sum_int4_2d
     module procedure mpp_sum_int4_3d
     module procedure mpp_sum_int4_4d
     module procedure mpp_sum_int4_5d
     module procedure mpp_sum_real4
     module procedure mpp_sum_real4_scalar
     module procedure mpp_sum_real4_2d
     module procedure mpp_sum_real4_3d
     module procedure mpp_sum_real4_4d
     module procedure mpp_sum_real4_5d
#ifdef OVERLOAD_C4
     module procedure mpp_sum_cmplx4
     module procedure mpp_sum_cmplx4_scalar
     module procedure mpp_sum_cmplx4_2d
     module procedure mpp_sum_cmplx4_3d
     module procedure mpp_sum_cmplx4_4d
     module procedure mpp_sum_cmplx4_5d
#endif
  end interface

  !> Calculates sum of a given numerical array across pe's for adjoint domains
  !> @ingroup mpp_mod
  interface mpp_sum_ad
     module procedure mpp_sum_int8_ad
     module procedure mpp_sum_int8_scalar_ad
     module procedure mpp_sum_int8_2d_ad
     module procedure mpp_sum_int8_3d_ad
     module procedure mpp_sum_int8_4d_ad
     module procedure mpp_sum_int8_5d_ad
     module procedure mpp_sum_real8_ad
     module procedure mpp_sum_real8_scalar_ad
     module procedure mpp_sum_real8_2d_ad
     module procedure mpp_sum_real8_3d_ad
     module procedure mpp_sum_real8_4d_ad
     module procedure mpp_sum_real8_5d_ad
#ifdef OVERLOAD_C8
     module procedure mpp_sum_cmplx8_ad
     module procedure mpp_sum_cmplx8_scalar_ad
     module procedure mpp_sum_cmplx8_2d_ad
     module procedure mpp_sum_cmplx8_3d_ad
     module procedure mpp_sum_cmplx8_4d_ad
     module procedure mpp_sum_cmplx8_5d_ad
#endif
     module procedure mpp_sum_int4_ad
     module procedure mpp_sum_int4_scalar_ad
     module procedure mpp_sum_int4_2d_ad
     module procedure mpp_sum_int4_3d_ad
     module procedure mpp_sum_int4_4d_ad
     module procedure mpp_sum_int4_5d_ad
     module procedure mpp_sum_real4_ad
     module procedure mpp_sum_real4_scalar_ad
     module procedure mpp_sum_real4_2d_ad
     module procedure mpp_sum_real4_3d_ad
     module procedure mpp_sum_real4_4d_ad
     module procedure mpp_sum_real4_5d_ad
#ifdef OVERLOAD_C4
     module procedure mpp_sum_cmplx4_ad
     module procedure mpp_sum_cmplx4_scalar_ad
     module procedure mpp_sum_cmplx4_2d_ad
     module procedure mpp_sum_cmplx4_3d_ad
     module procedure mpp_sum_cmplx4_4d_ad
     module procedure mpp_sum_cmplx4_5d_ad
#endif
  end interface

  !> @brief Gather information onto root pe
  !> @ingroup mpp_mod
  interface mpp_gather
     module procedure mpp_gather_logical_1d
     module procedure mpp_gather_int4_1d
     module procedure mpp_gather_real4_1d
     module procedure mpp_gather_real8_1d
     module procedure mpp_gather_logical_1dv
     module procedure mpp_gather_int4_1dv
     module procedure mpp_gather_real4_1dv
     module procedure mpp_gather_real8_1dv
     module procedure mpp_gather_pelist_logical_2d
     module procedure mpp_gather_pelist_logical_3d
     module procedure mpp_gather_pelist_int4_2d
     module procedure mpp_gather_pelist_int4_3d
     module procedure mpp_gather_pelist_real4_2d
     module procedure mpp_gather_pelist_real4_3d
     module procedure mpp_gather_pelist_real8_2d
     module procedure mpp_gather_pelist_real8_3d
  end interface

  !> @brief Scatter information to the given pelist
  !> @ingroup mpp_mod
  interface mpp_scatter
     module procedure mpp_scatter_pelist_int4_2d
     module procedure mpp_scatter_pelist_int4_3d
     module procedure mpp_scatter_pelist_real4_2d
     module procedure mpp_scatter_pelist_real4_3d
     module procedure mpp_scatter_pelist_real8_2d
     module procedure mpp_scatter_pelist_real8_3d
  end interface

  !#####################################################################
  !> @brief Scatter a vector across all PEs
  !!
  !> Transpose the vector and PE index
  !> @ingroup mpp_mod
  interface mpp_alltoall
     module procedure mpp_alltoall_int4
     module procedure mpp_alltoall_int8
     module procedure mpp_alltoall_real4
     module procedure mpp_alltoall_real8
#ifdef OVERLOAD_C4
     module procedure mpp_alltoall_cmplx4
#endif
#ifdef OVERLOAD_C8
     module procedure mpp_alltoall_cmplx8
#endif
     module procedure mpp_alltoall_logical4
     module procedure mpp_alltoall_logical8
     module procedure mpp_alltoall_int4_v
     module procedure mpp_alltoall_int8_v
     module procedure mpp_alltoall_real4_v
     module procedure mpp_alltoall_real8_v
#ifdef OVERLOAD_C4
     module procedure mpp_alltoall_cmplx4_v
#endif
#ifdef OVERLOAD_C8
     module procedure mpp_alltoall_cmplx8_v
#endif
     module procedure mpp_alltoall_logical4_v
     module procedure mpp_alltoall_logical8_v
     module procedure mpp_alltoall_int4_w
     module procedure mpp_alltoall_int8_w
     module procedure mpp_alltoall_real4_w
     module procedure mpp_alltoall_real8_w
#ifdef OVERLOAD_C4
     module procedure mpp_alltoall_cmplx4_w
#endif
#ifdef OVERLOAD_C8
     module procedure mpp_alltoall_cmplx8_w
#endif
     module procedure mpp_alltoall_logical4_w
     module procedure mpp_alltoall_logical8_w
  end interface


  !#####################################################################
  !> @brief Basic message-passing call.
  !!
  !>    <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !!    <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
  !!    contiguous block from a multi-dimensional array may be passed by its
  !!    starting address and its length, as in <TT>f77</TT>.
  !!
  !!    <TT>mpp_transmit</TT> is currently implemented as asynchronous
  !!    outward transmission and synchronous inward transmission. This follows
  !!    the behaviour of <TT>shmem_put</TT> and <TT>shmem_get</TT>. In MPI, it
  !!    is implemented as <TT>mpi_isend</TT> and <TT>mpi_recv</TT>. For most
  !!    applications, transmissions occur in pairs, and are here accomplished
  !!    in a single call.
  !!
  !!    The special PE designations <TT>NULL_PE</TT>,
  !!    <TT>ANY_PE</TT> and <TT>ALL_PES</TT> are provided by use
  !!    association.
  !!
  !!    <TT>NULL_PE</TT>: is used to disable one of the pair of
  !!    transmissions.<BR/>
  !!    <TT>ANY_PE</TT>: is used for unspecific remote
  !!    destination. (Please note that <TT>put_pe=ANY_PE</TT> has no meaning
  !!    in the MPI context, though it is available in the SHMEM invocation. If
  !!    portability is a concern, it is best avoided).<BR/>
  !!    <TT>ALL_PES</TT>: is used for broadcast operations.
  !!
  !!    It is recommended that <LINK
  !!    SRC="#mpp_broadcast"><TT>mpp_broadcast</TT></LINK> be used for
  !!    broadcasts.
  !!
  !!    The following example illustrates the use of
  !!    <TT>NULL_PE</TT> and <TT>ALL_PES</TT>:
  !!
  !!    <PRE>
  !!    real, dimension(n) :: a
  !!    if( pe.EQ.0 )then
  !!        do p = 1,npes-1
  !!           call mpp_transmit( a, n, p, a, n, NULL_PE )
  !!        end do
  !!    else
  !!        call mpp_transmit( a, n, NULL_PE, a, n, 0 )
  !!    end if
  !!
  !!    call mpp_transmit( a, n, ALL_PES, a, n, 0 )
  !!    </PRE>
  !!
  !!    The do loop and the broadcast operation above are equivalent.
  !!
  !!    Two overloaded calls <TT>mpp_send</TT> and
  !!     <TT>mpp_recv</TT> have also been
  !!    provided. <TT>mpp_send</TT> calls <TT>mpp_transmit</TT>
  !!    with <TT>get_pe=NULL_PE</TT>. <TT>mpp_recv</TT> calls
  !!    <TT>mpp_transmit</TT> with <TT>put_pe=NULL_PE</TT>. Thus
  !!    the do loop above could be written more succinctly:
  !!
  !!    <PRE>
  !!    if( pe.EQ.0 )then
  !!        do p = 1,npes-1
  !!           call mpp_send( a, n, p )
  !!        end do
  !!    else
  !!        call mpp_recv( a, n, 0 )
  !!    end if
  !!    </PRE>
  !! <br>Example call:
  !! @code{.F90}
  !!    call mpp_transmit( put_data, put_len, put_pe, get_data, get_len, get_pe )
  !! @endcode
  !> @ingroup mpp_mod
  interface mpp_transmit
     module procedure mpp_transmit_real8
     module procedure mpp_transmit_real8_scalar
     module procedure mpp_transmit_real8_2d
     module procedure mpp_transmit_real8_3d
     module procedure mpp_transmit_real8_4d
     module procedure mpp_transmit_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_transmit_cmplx8
     module procedure mpp_transmit_cmplx8_scalar
     module procedure mpp_transmit_cmplx8_2d
     module procedure mpp_transmit_cmplx8_3d
     module procedure mpp_transmit_cmplx8_4d
     module procedure mpp_transmit_cmplx8_5d
#endif
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

     module procedure mpp_transmit_real4
     module procedure mpp_transmit_real4_scalar
     module procedure mpp_transmit_real4_2d
     module procedure mpp_transmit_real4_3d
     module procedure mpp_transmit_real4_4d
     module procedure mpp_transmit_real4_5d

#ifdef OVERLOAD_C4
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
  !> @brief Recieve data to another PE
  !!
  !> @param[out] get_data scalar or array to get written with received data
  !> @param get_len size of array to recv from get_data
  !> @param from_pe PE number to receive from
  !> @param block true for blocking, false for non-blocking. Defaults to true
  !> @param tag communication tag
  !> @param[out] request MPI request handle
  !> @ingroup mpp_mod
  interface mpp_recv
     module procedure mpp_recv_real8
     module procedure mpp_recv_real8_scalar
     module procedure mpp_recv_real8_2d
     module procedure mpp_recv_real8_3d
     module procedure mpp_recv_real8_4d
     module procedure mpp_recv_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_recv_cmplx8
     module procedure mpp_recv_cmplx8_scalar
     module procedure mpp_recv_cmplx8_2d
     module procedure mpp_recv_cmplx8_3d
     module procedure mpp_recv_cmplx8_4d
     module procedure mpp_recv_cmplx8_5d
#endif
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

     module procedure mpp_recv_real4
     module procedure mpp_recv_real4_scalar
     module procedure mpp_recv_real4_2d
     module procedure mpp_recv_real4_3d
     module procedure mpp_recv_real4_4d
     module procedure mpp_recv_real4_5d

#ifdef OVERLOAD_C4
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
  !> Send data to a receiving PE.
  !!
  !> @param put_data scalar or array to get sent to a receiving PE
  !> @param put_len size of data to send from put_data
  !> @param to_pe PE number to send to
  !> @param block true for blocking, false for non-blocking. Defaults to true
  !> @param tag communication tag
  !> @param[out] request MPI request handle
  !! <br> Example usage:
  !! @code{.F90} call mpp_send(data, ie, pe) @endcode
  !> @ingroup mpp_mod
  interface mpp_send
     module procedure mpp_send_real8
     module procedure mpp_send_real8_scalar
     module procedure mpp_send_real8_2d
     module procedure mpp_send_real8_3d
     module procedure mpp_send_real8_4d
     module procedure mpp_send_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_send_cmplx8
     module procedure mpp_send_cmplx8_scalar
     module procedure mpp_send_cmplx8_2d
     module procedure mpp_send_cmplx8_3d
     module procedure mpp_send_cmplx8_4d
     module procedure mpp_send_cmplx8_5d
#endif
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

     module procedure mpp_send_real4
     module procedure mpp_send_real4_scalar
     module procedure mpp_send_real4_2d
     module procedure mpp_send_real4_3d
     module procedure mpp_send_real4_4d
     module procedure mpp_send_real4_5d

#ifdef OVERLOAD_C4
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


  !> @brief Perform parallel broadcasts
  !!
  !> The <TT>mpp_broadcast</TT> call has been added because the original
  !! syntax (using <TT>ALL_PES</TT> in <TT>mpp_transmit</TT>) did not
  !! support a broadcast across a pelist.
  !!
  !! <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !! <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
  !! contiguous block from a multi-dimensional array may be passed by its
  !! starting address and its length, as in <TT>f77</TT>.
  !!
  !! Global broadcasts through the <TT>ALL_PES</TT> argument to <LINK
  !! SRC="#mpp_transmit"><TT>mpp_transmit</TT></LINK> are still provided for
  !! backward-compatibility.
  !!
  !! If <TT>pelist</TT> is omitted, the context is assumed to be the
  !! current pelist. <TT>from_pe</TT> must belong to the current
  !! pelist. This call implies synchronization across the PEs in
  !! <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !!
  !! <br>Example usage:
  !!
  !!            call mpp_broadcast( data, length, from_pe, pelist )
  !!
  !> @param[inout] data Data to broadcast
  !> @param length Length of data to broadcast
  !> @param from_pe PE to send the data from
  !> @param pelist List of PE's to broadcast across, if not provided uses current list
  !> @ingroup mpp_mod
  interface mpp_broadcast
     module procedure mpp_broadcast_char
     module procedure mpp_broadcast_real8
     module procedure mpp_broadcast_real8_scalar
     module procedure mpp_broadcast_real8_2d
     module procedure mpp_broadcast_real8_3d
     module procedure mpp_broadcast_real8_4d
     module procedure mpp_broadcast_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_broadcast_cmplx8
     module procedure mpp_broadcast_cmplx8_scalar
     module procedure mpp_broadcast_cmplx8_2d
     module procedure mpp_broadcast_cmplx8_3d
     module procedure mpp_broadcast_cmplx8_4d
     module procedure mpp_broadcast_cmplx8_5d
#endif
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

     module procedure mpp_broadcast_real4
     module procedure mpp_broadcast_real4_scalar
     module procedure mpp_broadcast_real4_2d
     module procedure mpp_broadcast_real4_3d
     module procedure mpp_broadcast_real4_4d
     module procedure mpp_broadcast_real4_5d

#ifdef OVERLOAD_C4
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

  !> @brief Calculate parallel checksums
  !!
  !> \e mpp_chksum is a parallel checksum routine that returns an
  !! identical answer for the same array irrespective of how it has been
  !! partitioned across processors. \eint_kind is the KIND
  !! parameter corresponding to long integers (see discussion on
  !! OS-dependent preprocessor directives) defined in
  !! the file platform.F90. \eMPP_TYPE_ corresponds to any
  !! 4-byte and 8-byte variant of \einteger, \ereal, \ecomplex, \elogical
  !! variables, of rank 0 to 5.
  !!
  !! Integer checksums on FP data use the F90 <TT>TRANSFER()</TT>
  !! intrinsic.
  !!
  !! The <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/chksum/chksum.html">serial checksum module</LINK> is superseded
  !! by this function, and is no longer being actively maintained. This
  !! provides identical results on a single-processor job, and to perform
  !! serial checksums on a single processor of a parallel job, you only
  !! need to use the optional <TT>pelist</TT> argument.
  !! <PRE>
  !! use mpp_mod
  !! integer :: pe, chksum
  !! real :: a(:)
  !! pe = mpp_pe()
  !! chksum = mpp_chksum( a, (/pe/) )
  !! </PRE>
  !!
  !! The additional functionality of <TT>mpp_chksum</TT> over
  !! serial checksums is to compute the checksum across the PEs in
  !! <TT>pelist</TT>. The answer is guaranteed to be the same for
  !! the same distributed array irrespective of how it has been
  !! partitioned.
  !!
  !! If <TT>pelist</TT> is omitted, the context is assumed to be the
  !! current pelist. This call implies synchronization across the PEs in
  !! <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !! <br> Example usage:
  !!
  !!            mpp_chksum( var, pelist )
  !!
  !! @param var Data to calculate checksum of
  !! @param pelist Optional list of PE's to include in checksum calculation if not using
  !! current pelist
  !! @return Parallel checksum of var across given or implicit pelist
  !> @ingroup mpp_mod
  interface mpp_chksum
     module procedure mpp_chksum_i8_1d
     module procedure mpp_chksum_i8_2d
     module procedure mpp_chksum_i8_3d
     module procedure mpp_chksum_i8_4d
     module procedure mpp_chksum_i8_5d
     module procedure mpp_chksum_i8_1d_rmask
     module procedure mpp_chksum_i8_2d_rmask
     module procedure mpp_chksum_i8_3d_rmask
     module procedure mpp_chksum_i8_4d_rmask
     module procedure mpp_chksum_i8_5d_rmask

     module procedure mpp_chksum_i4_1d
     module procedure mpp_chksum_i4_2d
     module procedure mpp_chksum_i4_3d
     module procedure mpp_chksum_i4_4d
     module procedure mpp_chksum_i4_5d
     module procedure mpp_chksum_i4_1d_rmask
     module procedure mpp_chksum_i4_2d_rmask
     module procedure mpp_chksum_i4_3d_rmask
     module procedure mpp_chksum_i4_4d_rmask
     module procedure mpp_chksum_i4_5d_rmask

     module procedure mpp_chksum_r8_0d
     module procedure mpp_chksum_r8_1d
     module procedure mpp_chksum_r8_2d
     module procedure mpp_chksum_r8_3d
     module procedure mpp_chksum_r8_4d
     module procedure mpp_chksum_r8_5d

     module procedure mpp_chksum_r4_0d
     module procedure mpp_chksum_r4_1d
     module procedure mpp_chksum_r4_2d
     module procedure mpp_chksum_r4_3d
     module procedure mpp_chksum_r4_4d
     module procedure mpp_chksum_r4_5d
#ifdef OVERLOAD_C8
     module procedure mpp_chksum_c8_0d
     module procedure mpp_chksum_c8_1d
     module procedure mpp_chksum_c8_2d
     module procedure mpp_chksum_c8_3d
     module procedure mpp_chksum_c8_4d
     module procedure mpp_chksum_c8_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_chksum_c4_0d
     module procedure mpp_chksum_c4_1d
     module procedure mpp_chksum_c4_2d
     module procedure mpp_chksum_c4_3d
     module procedure mpp_chksum_c4_4d
     module procedure mpp_chksum_c4_5d
#endif
  end interface

!> @addtogroup mpp_mod
!> @{
!***********************************************************************
!
!            module variables
!
!***********************************************************************
  integer, parameter   :: PESET_MAX = 10000
  integer              :: current_peset_max = 32
  type(communicator), allocatable :: peset(:) !< Will be allocated starting from 0, 0 is a dummy used to hold single-PE "self" communicator
  logical              :: module_is_initialized = .false.
  logical              :: debug = .false.
  integer              :: npes=1, root_pe=0, pe=0
  integer(i8_kind)     :: tick, ticks_per_sec, max_ticks, start_tick, end_tick, tick0=0
  integer              :: mpp_comm_private
  logical              :: first_call_system_clock_mpi=.TRUE.
  real(r8_kind)        :: mpi_count0=0  !< use to prevent integer overflow
  real(r8_kind)        :: mpi_tick_rate=0.d0  !< clock rate for mpi_wtick()
  logical              :: mpp_record_timing_data=.TRUE.
  type(clock),save     :: clocks(MAX_CLOCKS)
  integer              :: log_unit, etc_unit
  character(len=32)    :: configfile='logfile'
  integer              :: peset_num=0, current_peset_num=0
  integer              :: world_peset_num                  !<the world communicator
  integer              :: error
  integer              :: clock_num=0, num_clock_ids=0,current_clock=0, previous_clock(MAX_CLOCKS)=0
  real                 :: tick_rate

  type(mpp_type_list)    :: datatypes
  type(mpp_type), target :: mpp_byte

  integer              :: cur_send_request = 0
  integer              :: cur_recv_request = 0
  integer, allocatable :: request_send(:)
  integer, allocatable :: request_recv(:)
  integer, allocatable :: size_recv(:)
  integer, allocatable :: type_recv(:)
! if you want to save the non-root PE information uncomment out the following line
! and comment out the assigment of etcfile to '/dev/null'
#ifdef NO_DEV_NULL
  character(len=32)    :: etcfile='._mpp.nonrootpe.msgs'
#else
  character(len=32)    :: etcfile='/dev/null'
#endif

!> Use the intrinsics in iso_fortran_env
  integer :: in_unit=INPUT_UNIT, out_unit=OUTPUT_UNIT, err_unit=ERROR_UNIT
  integer :: stdout_unit

  !--- variables used in mpp_util.h
  type(Summary_Struct) :: clock_summary(MAX_CLOCKS)
  logical              :: warnings_are_fatal = .FALSE.
  integer              :: error_state=0
  integer              :: clock_grain=CLOCK_LOOP-1

  !--- variables used in mpp_comm.h
  integer            :: clock0    !<measures total runtime from mpp_init to mpp_exit
  integer            :: mpp_stack_size=0, mpp_stack_hwm=0
  logical            :: verbose=.FALSE.

  integer :: get_len_nocomm = 0 !< needed for mpp_transmit_nocomm.h

  !--- variables used in mpp_comm_mpi.inc
  integer, parameter :: mpp_init_test_full_init = -1
  integer, parameter :: mpp_init_test_init_true_only = 0
  integer, parameter :: mpp_init_test_peset_allocated = 1
  integer, parameter :: mpp_init_test_clocks_init = 2
  integer, parameter :: mpp_init_test_datatype_list_init = 3
  integer, parameter :: mpp_init_test_logfile_init = 4
  integer, parameter :: mpp_init_test_read_namelist = 5
  integer, parameter :: mpp_init_test_etc_unit = 6
  integer, parameter :: mpp_init_test_requests_allocated = 7


!***********************************************************************
!  variables needed for subroutine read_input_nml (include/mpp_util.inc)
!
! public variable needed for reading input nml file from an internal file
  character(len=:), dimension(:), allocatable, target, public :: input_nml_file
  logical :: read_ascii_file_on = .FALSE.
!***********************************************************************

! Include variable "version" to be written to log file.
#include<file_version.h>
  public version

  integer, parameter :: MAX_REQUEST_MIN  = 10000
  integer            :: request_multiply = 20

  logical :: etc_unit_is_stderr = .false.
  integer :: max_request = 0
  logical :: sync_all_clocks = .false.
  namelist /mpp_nml/ etc_unit_is_stderr, request_multiply, mpp_record_timing_data, sync_all_clocks

  contains
#include <system_clock.h>
#include <mpp_util.inc>
#include <mpp_comm.inc>

  end module mpp_mod
!> @}
! close documentation grouping
