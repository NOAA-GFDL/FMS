module mpp_datatype_mod
#include <fms_platform.h>

  use mpp_parameter_mod, only : MAX_EVENTS, MAX_EVENT_TYPES, MAX_BINS, PESET_MAX, MAX_CLOCKS
#ifdef use_CAF
  use mpp_parameter_mod, only : MAX_DOMAIN_FIELDS
#endif

  implicit none
  private

  character(len=128), public :: version= &
       '$Id mpp_datatype.F90 $'
  character(len=128), public :: tagname= &
       '$Name: lima $'

  !--- public data type which is used by mpp_mod and its components. 
  !--- All othere modules should import these parameters from mpp_mod.  
  public :: communicator, event, clock, Clock_Data_Summary, Summary_Struct

  !--- public data type which is used by mpp_domains_mod and its components. All othere modules should import
  !--- these parameters from mpp_domains_mod. 
  public :: domain_axis_spec, domain1D, domain2D, DomainCommunicator2D, rectangle

#ifdef use_CAF
  public :: cafptr_r8_3d_type
  public :: cafptr_r8_3d
  public :: cafptr_r8_1d_type
  public :: cafptr_r8_1d
  public :: cafptr_c8_3d_type
  public :: cafptr_c8_3d
  public :: cafptr_c8_1d_type
  public :: cafptr_c8_1d
#ifndef no_8byte_integers
  public :: cafptr_i8_3d_type
  public :: cafptr_i8_3d
  public :: cafptr_i8_1d_type
  public :: cafptr_i8_1d
  public :: cafptr_l8_3d_type
  public :: cafptr_l8_3d
  public :: cafptr_l8_1d_type
  public :: cafptr_l8_1d

#endif
#ifndef no_4byte_reals
  public :: cafptr_r4_3d_type
  public :: cafptr_r4_3d
  public :: cafptr_r4_1d_type
  public :: cafptr_r4_1d
  public :: cafptr_c4_3d_type
  public :: cafptr_c4_3d
  public :: cafptr_c4_1d_type
  public :: cafptr_c4_1d
#endif
  public :: cafptr_i4_3d_type
  public :: cafptr_i4_3d
  public :: cafptr_i4_1d_type
  public :: cafptr_i4_1d
  public :: cafptr_l4_3d_type
  public :: cafptr_l4_3d
  public :: cafptr_l4_1d_type
  public :: cafptr_l4_1d

#endif


  !--- public data type which is used by mpp_io_mod and its components. All othere modules should import
  !--- these parameters from mpp_io_mod. 
  public :: axistype, atttype, validtype, fieldtype, filetype 

  !----------------------------------------------------------------------
  !--- The following data types are used by mpp_mod and its components.
  !----------------------------------------------------------------------
  !peset hold communicators as SHMEM-compatible triads (start, log2(stride), num)
  type :: communicator
     character(len=32) :: name
     integer, pointer  :: list(:) =>NULL()
     integer           :: count
     integer           :: start, log2stride ! dummy variables when libMPI is defined.
     integer           :: id, group         ! MPI communicator and group id for this PE set.
                                            ! dummy variables when libSMA is defined.
  end type communicator

  type :: event
     character(len=16)                         :: name
     integer(LONG_KIND), dimension(MAX_EVENTS) :: ticks, bytes
     integer                                   :: calls
  end type event

  !a clock contains an array of event profiles for a region
  type :: clock
     character(len=32)    :: name
     integer(LONG_KIND)   :: tick
     integer(LONG_KIND)   :: total_ticks
     integer              :: peset_num
     logical              :: sync_on_begin, detailed
     integer              :: grain
     type(event), pointer :: events(:) =>NULL() !if needed, allocate to MAX_EVENT_TYPES
  end type clock

  type :: Clock_Data_Summary
     character(len=16)  :: name
     real(DOUBLE_KIND)  :: msg_size_sums(MAX_BINS)
     real(DOUBLE_KIND)  :: msg_time_sums(MAX_BINS)
     real(DOUBLE_KIND)  :: total_data
     real(DOUBLE_KIND)  :: total_time
     integer(LONG_KIND) :: msg_size_cnts(MAX_BINS)
     integer(LONG_KIND) :: total_cnts
  end type Clock_Data_Summary

  type :: Summary_Struct
     character(len=16)         :: name
     type (Clock_Data_Summary) :: event(MAX_EVENT_TYPES)
  end type Summary_Struct


  !-----------------------------------------------------------------------------
  !--- The following data types are used by mpp_domains_mod and its components.
  !-----------------------------------------------------------------------------
  type domain_axis_spec        !type used to specify index limits along an axis of a domain
     integer :: begin, end, size, max_size      !start, end of domain axis, size, max size in set
     logical :: is_global       !TRUE if domain axis extent covers global domain
  end type domain_axis_spec
  type domain1D
     type(domain_axis_spec) :: compute, data, global
     logical :: cyclic
     type(domain1D), pointer :: list(:) =>NULL()
     integer :: pe              !PE to which this domain is assigned
     integer :: pos             !position of this PE within link list, i.e domain%list(pos)%pe = pe
  end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain, not once per axis
  type rectangle
     integer :: is, ie, js, je
     logical :: overlap, folded
  end type rectangle
  type domain2D
     integer(LONG_KIND) :: id 
     type(domain1D) :: x
     type(domain1D) :: y
     type(domain2D), pointer :: list(:) =>NULL()
     integer :: pe              !PE to which this domain is assigned
     integer :: pos             !position of this PE within link list, i.e domain%list(pos)%pe = pe
     integer :: fold, gridtype
     logical :: overlap
     type(rectangle) :: recv_e, recv_se, recv_s, recv_sw, &
                        recv_w, recv_nw, recv_n, recv_ne
     type(rectangle) :: send_e, send_se, send_s, send_sw, &
                        send_w, send_nw, send_n, send_ne
     logical :: remote_domains_initialized
     type(rectangle) :: recv_e_off, recv_se_off, recv_s_off, recv_sw_off, &
                        recv_w_off, recv_nw_off, recv_n_off, recv_ne_off
     type(rectangle) :: send_e_off, send_se_off, send_s_off, send_sw_off, &
                        send_w_off, send_nw_off, send_n_off, send_ne_off
     logical :: remote_off_domains_initialized
  end type domain2D     

  type DomainCommunicator2D
     logical            :: initialized=.false.
     integer(LONG_KIND) :: id=-9999
     integer(LONG_KIND) :: l_addr  =-9999
     integer(LONG_KIND) :: l_addrx =-9999
     integer(LONG_KIND) :: l_addry =-9999
     type(domain2D), pointer :: domain     =>NULL()
     type(domain2D), pointer :: domain_in  =>NULL()
     type(domain2D), pointer :: domain_out =>NULL()
     integer, dimension(:,:), _ALLOCATABLE :: sendis _NULL
     integer, dimension(:,:), _ALLOCATABLE :: sendie _NULL
     integer, dimension(:,:), _ALLOCATABLE :: sendjs _NULL
     integer, dimension(:,:), _ALLOCATABLE :: sendje _NULL
     integer, dimension(:,:), _ALLOCATABLE :: S_msize _NULL
     logical, dimension(:,:), _ALLOCATABLE :: do_thisS _NULL
     logical, dimension(:), _ALLOCATABLE :: S_do_buf _NULL
     integer, dimension(:), _ALLOCATABLE :: cto_pe  _NULL
     integer, dimension(:), _ALLOCATABLE :: Rcaf_idx  _NULL
     integer, dimension(:,:), _ALLOCATABLE :: recvis _NULL
     integer, dimension(:,:), _ALLOCATABLE :: recvie _NULL
     integer, dimension(:,:), _ALLOCATABLE :: recvjs _NULL
     integer, dimension(:,:), _ALLOCATABLE :: recvje _NULL
     integer, dimension(:,:), _ALLOCATABLE :: R_msize _NULL
     logical, dimension(:,:,:), _ALLOCATABLE :: do_thisR _NULL
     logical, dimension(:), _ALLOCATABLE :: R_do_buf _NULL
     integer, dimension(:), _ALLOCATABLE :: cfrom_pe  _NULL
     integer :: Slist_size=0, Rlist_size=0
     integer :: isize=0, jsize=0, ke=0
     integer :: isize_in=0, jsize_in=0
     integer :: isize_out=0, jsize_out=0
     integer :: isize_max=0, jsize_max=0
     integer :: gf_ioff=0, gf_joff=0
  ! Remote data
     integer, dimension(:)  , _ALLOCATABLE :: isizeR _NULL
     integer, dimension(:)  , _ALLOCATABLE :: jsizeR _NULL
     integer, dimension(:,:), _ALLOCATABLE :: sendisR _NULL
     integer, dimension(:,:), _ALLOCATABLE :: sendjsR _NULL
     integer(LONG_KIND), dimension(:), _ALLOCATABLE :: rem_addr  _NULL
     integer(LONG_KIND), dimension(:), _ALLOCATABLE :: rem_addrx _NULL
     integer(LONG_KIND), dimension(:), _ALLOCATABLE :: rem_addry _NULL
     integer(LONG_KIND), dimension(:,:), _ALLOCATABLE :: rem_addrl  _NULL
     integer(LONG_KIND), dimension(:,:), _ALLOCATABLE :: rem_addrlx  _NULL
     integer(LONG_KIND), dimension(:,:), _ALLOCATABLE :: rem_addrly  _NULL
     integer(LONG_KIND), dimension(:), _ALLOCATABLE :: sync_start_list _NULL
     integer(LONG_KIND), dimension(:), _ALLOCATABLE :: sync_end_list _NULL
     type(DomainCommunicator2D), pointer :: dch_x =>NULL()
  end type DomainCommunicator2D


  !-----------------------------------------------------------------------------
  !--- The following data types are used by mpp_io_mod and its components.
  !-----------------------------------------------------------------------------
  type :: atttype
     integer             :: type, len
     character(len=128)  :: name
     character(len=1280) :: catt
     real, pointer       :: fatt(:) =>NULL() ! just use type conversion for integers
  end type atttype

  type :: axistype
     character(len=128) :: name
     character(len=128) :: units
     character(len=256) :: longname
     character(len=8)   :: cartesian
     character(len=24)  :: calendar
     integer            :: sense, len          !+/-1, depth or height?
     type(domain1D)     :: domain              !if pointer is associated, it is a distributed data axis
     real, pointer      :: data(:) =>NULL()    !axis values (not used if time axis)
     integer            :: id, did, type, natt !id is the "variable ID", did is the "dimension ID": 
                                               !netCDF requires 2 IDs for axes
     type(atttype), pointer :: Att(:) =>NULL()
  end type axistype

  type :: validtype
     logical :: is_range ! if true, then the data represent the valid range
     real    :: min,max  ! boundaries of the valid range or missing value
  end type validtype

  type :: fieldtype
     character(len=128)      :: name
     character(len=128)      :: units
     character(len=256)      :: longname
     character(len=128)      :: standard_name   ! CF standard name
     real                    :: min, max, missing, fill, scale, add
     integer                 :: pack
     type(axistype), pointer :: axes(:) =>NULL() !axes associated with field size, time_axis_index redundantly 
                                        !hold info already contained in axes. it's clunky and inelegant, 
                                        !but required so that axes can be shared among multiple files
     integer, pointer        :: size(:) =>NULL()
     integer                 :: time_axis_index
     integer                 :: id, type, natt, ndim
     type(atttype), pointer  :: Att(:) =>NULL()
  end type fieldtype

  type :: filetype
     character(len=256) :: name
     integer            :: action, format, access, threading, fileset, record, ncid
     logical            :: opened, initialized, nohdrs
     integer            :: time_level
     real(DOUBLE_KIND)  :: time
     integer            :: id       !variable ID of time axis associated with file (only one time axis per file)
     integer            :: recdimid !dim ID of time axis associated with file (only one time axis per file)
     real(DOUBLE_KIND), pointer :: time_values(:) =>NULL() ! time axis values are stored here instead of axis%data 
                                                  ! since mpp_write assumes these values are not time values. 
                                                  ! Not used in mpp_write
     ! additional elements of filetype for mpp_read (ignored for mpp_write)
     integer :: ndim, nvar, natt  ! number of dimensions, non-dimension variables and global attributes
                                  ! redundant axis types stored here and in associated fieldtype
                                  ! some axes are not used by any fields, i.e. "edges"
     type(axistype), pointer  :: axis(:) =>NULL()
     type(fieldtype), pointer :: var(:) =>NULL()
     type(atttype), pointer   :: att(:) =>NULL()
  end type filetype


!#######################################################################

#ifdef use_CAF
#define MPP_TYPE_ real(DOUBLE_KIND)
#define CAFPNTR_TYPE_3D_ cafptr_r8_3d_type
#define CAFPNTR_3D_ cafptr_r8_3d
#define CAFPNTR_TYPE_1D_ cafptr_r8_1d_type
#define CAFPNTR_1D_ cafptr_r8_1d
#include <mpp_datatype.h>

#define MPP_TYPE_ complex(DOUBLE_KIND)
#define CAFPNTR_TYPE_3D_ cafptr_c8_3d_type
#define CAFPNTR_3D_ cafptr_c8_3d
#define CAFPNTR_TYPE_1D_ cafptr_c8_1d_type
#define CAFPNTR_1D_ cafptr_c8_1d
#include <mpp_datatype.h>


#ifndef no_8byte_integers
#define MPP_TYPE_ integer(LONG_KIND)
#define CAFPNTR_TYPE_3D_ cafptr_i8_3d_type
#define CAFPNTR_3D_ cafptr_i8_3d
#define CAFPNTR_TYPE_1D_ cafptr_i8_1d_type
#define CAFPNTR_1D_ cafptr_i8_1d
#include <mpp_datatype.h>

#define MPP_TYPE_ logical(LONG_KIND)
#define CAFPNTR_TYPE_3D_ cafptr_l8_3d_type
#define CAFPNTR_3D_ cafptr_l8_3d
#define CAFPNTR_TYPE_1D_ cafptr_l8_1d_type
#define CAFPNTR_1D_ cafptr_l8_1d
#include <mpp_datatype.h>
#endif

#ifndef no_4byte_reals
#define MPP_TYPE_ real(FLOAT_KIND)
#define CAFPNTR_TYPE_3D_ cafptr_r4_3d_type
#define CAFPNTR_3D_ cafptr_r4_3d
#define CAFPNTR_TYPE_1D_ cafptr_r4_1d_type
#define CAFPNTR_1D_ cafptr_r4_1d
#include <mpp_datatype.h>

#define MPP_TYPE_ complex(FLOAT_KIND)
#define CAFPNTR_TYPE_3D_ cafptr_c4_3d_type
#define CAFPNTR_3D_ cafptr_c4_3d
#define CAFPNTR_TYPE_1D_ cafptr_c4_1d_type
#define CAFPNTR_1D_ cafptr_c4_1d
#include <mpp_datatype.h>
#endif

#define MPP_TYPE_ integer(INT_KIND)
#define CAFPNTR_TYPE_3D_ cafptr_i4_3d_type
#define CAFPNTR_3D_ cafptr_i4_3d
#define CAFPNTR_TYPE_1D_ cafptr_i4_1d_type
#define CAFPNTR_1D_ cafptr_i4_1d
#include <mpp_datatype.h>

#define MPP_TYPE_ logical(INT_KIND)
#define CAFPNTR_TYPE_3D_ cafptr_l4_3d_type
#define CAFPNTR_3D_ cafptr_l4_3d
#define CAFPNTR_TYPE_1D_ cafptr_l4_1d_type
#define CAFPNTR_1D_ cafptr_l4_1d
#include <mpp_datatype.h>
#endif

end module mpp_datatype_mod
