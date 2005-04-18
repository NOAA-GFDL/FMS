module mpp_data_mod
#include <fms_platform.h>

#if defined(use_libMPI) && defined(sgi_mipspro)
  use mpi
#endif

  use mpp_parameter_mod, only : PESET_MAX,  MAX_CLOCKS, MAXPES, AGRID
  use mpp_datatype_mod,  only : communicator, clock, domain1D, domain2D
  use mpp_datatype_mod,  only : filetype, axistype, fieldtype, atttype

  implicit none
  private

  character(len=128), public :: version= &
       '$Id mpp_data.F90 $'
  character(len=128), public :: tagname= &
       '$Name: lima $'

#if defined(use_libSMA) || defined(use_libGSM)
#include <mpp/shmem.fh>
#endif

#if defined(use_libMPI) && !defined(sgi_mipspro)
#include <mpif.h>  !sgi_mipspro gets this from 'use mpi'
#endif

  !--- public data which is shared by all the modules in mpp package.
  public :: pe, npes

  !--- public data  which is used by mpp_mod and its components. 
  !--- All othere modules should import these parameters from mpp_mod.  
  public :: peset, clocks, mpp_is_initialized, error, root_pe, debug_mpp
  public :: configfile, etcfile, current_peset_num, peset_num, world_peset_num               
  public :: clock_num, previous_clock, current_clock, num_clock_ids
  public :: tick_rate, tick0, ticks_per_sec, max_ticks, mpi_tick_rate, mpi_count0
  public :: first_call_system_clock_mpi, mpp_record_timing_data, start_tick, tick, end_tick
  public :: stat, request, mpp_stack, ptr_stack, status
  public :: ptr_status, sync, ptr_sync,  mpp_from_pe, ptr_from, remote_data_loc 
  public :: ptr_remote, in_unit, out_unit, err_unit, log_unit, etc_unit
  public :: mpp_comm_private

  !--- public data which is used by mpp_domains_mod and its components. 
  !--- All othere modules should import these parameters from mpp_domains_mod. 
  public :: mpp_domains_stack_size, mpp_domains_stack, ptr_domains_stack
  public :: mpp_domains_stack_hwm, domain_info_buf, ptr_info, grid_offset_type
  public :: NULL_DOMAIN1D, NULL_DOMAIN2D, debug_mpp_domains, verbose
  public :: debug_gsm,mpp_domains_is_initialized

  !--- public data which is used by mpp_io_mod and its components. All othere modules should import
  !--- these parameters from mpp_io_mod. 
  public :: mpp_io_is_initialized, mpp_file, default_axis, default_field, default_att
  public :: maxunits, unit_begin, unit_end, mpp_io_stack_size, mpp_io_stack_hwm
  public :: mpp_io_stack, debug_mpp_io, verbose_mpp_io, text, varnum, records_per_pe

  !--- The following data is shared by all the modules in mpp package
  integer              :: pe=0
  !----------------------------------------------------------------------!
  !  The following data are used by mpp_mod and its components.          !
  !----------------------------------------------------------------------!
  type(communicator),save :: peset(0:PESET_MAX) !0 is a dummy used to hold single-PE "self" communicator
  integer(LONG_KIND)   :: tick, ticks_per_sec, max_ticks, start_tick, end_tick, tick0=0
  integer              :: mpp_comm_private
  logical              :: first_call_system_clock_mpi=.TRUE.
  real(DOUBLE_KIND)    :: mpi_count0=0  ! use to prevent integer overflow
  real(DOUBLE_KIND)    :: mpi_tick_rate=0.d0  ! clock rate for mpi_wtick()
  logical              :: mpp_record_timing_data=.TRUE.
  type(clock),save     :: clocks(MAX_CLOCKS)
  logical              :: mpp_is_initialized=.FALSE.
  integer              :: npes=1, root_pe=0
  integer              :: log_unit, etc_unit
  character(len=32)    :: configfile='logfile'
  integer              :: peset_num=0, current_peset_num=0
  integer              :: world_peset_num                  !the world communicator
  logical              :: debug_mpp=.FALSE.
  integer              :: error
  integer              :: clock_num=0, num_clock_ids=0,current_clock=0, previous_clock(MAX_CLOCKS)=0
  real                 :: tick_rate
  integer, allocatable :: request(:)
  character(len=32)    :: etcfile='._mpp.nonrootpe.stdout'
#ifdef SGICRAY
  integer :: in_unit=100, out_unit=101, err_unit=102 !see intro_io(3F): to see why these values are used rather than 5,6,0
#else
  integer :: in_unit=5, out_unit=6, err_unit=0
#endif
  !----------------------------------------------------------------------------!
  ! The following data types are used by mpp_domains_mod and its components.   !
  !----------------------------------------------------------------------------!
  logical        :: mpp_domains_is_initialized = .false.
  logical        :: verbose=.FALSE., debug_mpp_domains =.FALSE.
  logical        :: debug_gsm=.false.
  integer        :: grid_offset_type = AGRID
  integer        :: mpp_domains_stack_size=0
  integer        :: mpp_domains_stack_hwm=0
  integer        :: domain_info_buf(16)
  type(domain1D),save :: NULL_DOMAIN1D
  type(domain2D),save :: NULL_DOMAIN2D

  !----------------------------------------------------------------------------!
  ! The following data types are used by mpp_io_mod and its components.   !
  !----------------------------------------------------------------------------!
  logical            :: verbose_mpp_io =.FALSE.
  logical            :: debug_mpp_io = .FALSE.
  logical            :: mpp_io_is_initialized = .FALSE.
  type(axistype),save     :: default_axis      !provided to users with default components
  type(fieldtype),save    :: default_field     !provided to users with default components
  type(atttype),save      :: default_att       !provided to users with default components
  integer            :: maxunits, unit_begin, unit_end
  integer            :: mpp_io_stack_size=0, mpp_io_stack_hwm=0
  integer            :: varnum=0
  character(len=256) :: text
  integer            :: records_per_pe
  type(filetype),    allocatable :: mpp_file(:)
  real(DOUBLE_KIND), allocatable :: mpp_io_stack(:)
  !-------------------------------------------------------------------------------!
  ! The following data included in the .inc file are diffrent for sma or mpi case !
  !-------------------------------------------------------------------------------!

#ifdef use_libSMA
#include <mpp_data_sma.inc>
#elif use_libMPI
#include <mpp_data_mpi.inc>
#else
#include <mpp_data_nocomm.inc>
#endif

end module mpp_data_mod
