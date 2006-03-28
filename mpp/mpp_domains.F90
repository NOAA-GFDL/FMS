!-----------------------------------------------------------------------
!   Domain decomposition and domain update for message-passing codes
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

! <CONTACT EMAIL="vb@gfdl.noaa.gov">
!   V. Balaji
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <RCSLOG SRC="http://www.gfdl.noaa.gov/~vb/changes_mpp_domains.html"/>

! <OVERVIEW>
!   <TT>mpp_domains_mod</TT> is a set of simple calls for domain
!   decomposition and domain updates on rectilinear grids. It requires the
!   module <LINK SRC="mpp.html">mpp_mod</LINK>, upon which it is built.
! </OVERVIEW>

! <DESCRIPTION>
!   Scalable implementations of finite-difference codes are generally
!   based on decomposing the model domain into subdomains that are
!   distributed among processors. These domains will then be obliged to
!   exchange data at their boundaries if data dependencies are merely
!   nearest-neighbour, or may need to acquire information from the global
!   domain if there are extended data dependencies, as in the spectral
!   transform. The domain decomposition is a key operation in the
!   development of parallel codes.
!   
!   <TT>mpp_domains_mod</TT> provides a domain decomposition and domain
!   update API for <I>rectilinear</I> grids, built on top of the <LINK
!   SRC="mpp.html">mpp_mod</LINK> API for message passing. Features
!   of <TT>mpp_domains_mod</TT> include:
! 
!   Simple, minimal API, with free access to underlying API for more complicated stuff.
!
!   Design toward typical use in climate/weather CFD codes.
!  
!   <H4>Domains</H4>
! 
!   I have assumed that domain decomposition will mainly be in 2
!   horizontal dimensions, which will in general be the two
!   fastest-varying indices. There is a separate implementation of 1D
!   decomposition on the fastest-varying index, and 1D decomposition on
!   the second index, treated as a special case of 2D decomposition, is
!   also possible. We define <I>domain</I> as the grid associated with a <I>task</I>.
!   We define the <I>compute domain</I> as the set of gridpoints that are
!   computed by a task, and the <I>data domain</I> as the set of points
!   that are required by the task for the calculation. There can in
!   general be more than 1 task per PE, though often
!   the number of domains is the same as the processor count. We define
!   the <I>global domain</I> as the global computational domain of the
!   entire model (i.e, the same as the computational domain if run on a
!   single processor). 2D domains are defined using a derived type <TT>domain2D</TT>,
!   constructed as follows (see comments in code for more details):
!   
!   <PRE>
!     type, public :: domain_axis_spec
!        private
!        integer :: begin, end, size, max_size
!        logical :: is_global
!     end type domain_axis_spec
!     type, public :: domain1D
!        private
!        type(domain_axis_spec) :: compute, data, global, active
!        logical :: mustputb, mustgetb, mustputf, mustgetf, folded
!        type(domain1D), pointer, dimension(:) :: list
!        integer :: pe              !PE to which this domain is assigned
!        integer :: pos
!     end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain
!     type, public :: domain2D
!        private
!        type(domain1D) :: x
!        type(domain1D) :: y
!        type(domain2D), pointer, dimension(:) :: list
!        integer :: pe              !PE to which this domain is assigned
!        integer :: pos
!     end type domain2D
!     type(domain1D), public :: NULL_DOMAIN1D
!     type(domain2D), public :: NULL_DOMAIN2D
!   </PRE>

!   The <TT>domain2D</TT> type contains all the necessary information to
!   define the global, compute and data domains of each task, as well as the PE
!   associated with the task. The PEs from which remote data may be
!   acquired to update the data domain are also contained in a linked list
!   of neighbours.
! </DESCRIPTION>

module mpp_domains_mod
!a generalized domain decomposition package for use with mpp_mod
!Balaji (vb@gfdl.gov) 15 March 1999
  use mpp_parameter_mod,      only : MPP_DEBUG, MPP_VERBOSE, MPP_DOMAIN_TIME
  use mpp_parameter_mod,      only : GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, GLOBAL,CYCLIC 
  use mpp_parameter_mod,      only : AGRID, BGRID_SW, BGRID_NE, CGRID_NE, CGRID_SW, FOLD_WEST_EDGE
  use mpp_parameter_mod,      only : FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE
  use mpp_parameter_mod,      only : WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE
  use mpp_parameter_mod,      only : NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, MPP_DOMAIN_TIME
  use mpp_parameter_mod,      only : CENTER, CORNER, SCALAR_PAIR, SCALAR_BIT
  use mpp_parameter_mod,      only : NORTH, NORTH_EAST, EAST, SOUTH_EAST
  use mpp_parameter_mod,      only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST
  use mpp_parameter_mod,      only : MAX_DOMAIN_FIELDS, NULL_PE, CUBIC_GRID, REGULAR, DOMAIN_ID_BASE
  use mpp_parameter_mod,      only : ZERO, NINETY, MINUS_NINETY, ONE_HUNDRED_EIGHTY 
  use mpp_data_mod,           only : mpp_domains_stack, ptr_domains_stack
  use mpp_data_mod,           only : domain_info_buf, ptr_info
  use mpp_mod,                only : mpp_pe, mpp_root_pe, mpp_npes, mpp_error, FATAL, WARNING, NOTE
  use mpp_mod,                only : stdout, stderr, stdlog, mpp_send, mpp_recv, mpp_transmit, mpp_sync_self
  use mpp_mod,                only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_mod,                only : mpp_max, mpp_min, mpp_sum, mpp_get_current_pelist, mpp_broadcast
  use mpp_mod,                only : mpp_sync, mpp_init, mpp_malloc
  use mpp_pset_mod, only: mpp_pset_init
  implicit none
  private

#include <fms_platform.h>

  !--- public paramters imported from mpp_domains_parameter_mod
  public :: GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, BGRID_NE, BGRID_SW, CGRID_NE
  public :: CGRID_SW, FOLD_WEST_EDGE, FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE
  public :: WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE
  public :: NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, MPP_DOMAIN_TIME
  public :: CENTER, CORNER, SCALAR_PAIR
  public :: NORTH, NORTH_EAST, EAST, SOUTH_EAST
  public :: SOUTH, SOUTH_WEST, WEST, NORTH_WEST

  !--- public data imported from mpp_data_mod
  public :: NULL_DOMAIN1D, NULL_DOMAIN2D

  public :: domain_axis_spec, domain1D, domain2D, DomainCommunicator2D

  !--- public interface from mpp_domains_util.h
  public :: mpp_domains_set_stack_size, mpp_get_compute_domain, mpp_get_compute_domains
  public :: mpp_get_data_domain, mpp_get_global_domain, mpp_get_domain_components
  public :: mpp_get_layout, mpp_get_pelist, operator(.EQ.), operator(.NE.) 
  public :: mpp_get_global_shift, mpp_domain_is_symmetry
  public :: mpp_get_neighbor_pe, mpp_nullify_domain_list
  public :: mpp_set_compute_domain, mpp_set_data_domain, mpp_set_global_domain
  !--- public interface from mpp_domains_reduce.h
  public :: mpp_global_field, mpp_global_max, mpp_global_min, mpp_global_sum
  !--- public interface from mpp_domains_misc.h
  public :: mpp_broadcast_domain, mpp_domains_init, mpp_domains_exit, mpp_redistribute
  public ::  mpp_update_domains, mpp_check_field
  !--- public interface from mpp_domains_define.h
  public :: mpp_define_layout, mpp_define_domains, mpp_modify_domain, mpp_define_mosaic


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

  !--- data types used mpp_domains_mod.
  type domain_axis_spec        !type used to specify index limits along an axis of a domain
     private
     integer :: begin, end, size, max_size      !start, end of domain axis, size, max size in set
     logical :: is_global       !TRUE if domain axis extent covers global domain
  end type domain_axis_spec

  type domain1D
     private
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
     private
     integer :: is, ie, js, je
     logical :: overlap, folded
  end type rectangle

  type overlapList
     private
     integer          :: n              ! number of points to send/recv
     integer, pointer :: i(:) => NULL() ! i-index of each point to send/recv
     integer, pointer :: j(:) => NULL() ! j-index of each point to send/recv
     integer          :: is, ie, js, je ! overlapping region.
     logical          :: overlap(4)     ! 1 for rect overlap within tile, 2 for list overlap within tile
                                        ! 3 for overlap between tile, 4 for corner fix between tiles.
     logical          :: folded         ! indicate if the overlap is folded.
     integer          :: rotation       ! rotate angle between tiles.
     integer          :: i2, j2         ! for corner point update ( cubic grid )
   end type overlapList

  type domain2D
     private
     integer(LONG_KIND)      :: id 
     type(domain1D)          :: x
     type(domain1D)          :: y
     type(domain2D), pointer :: list(:) =>NULL()
     integer,        pointer :: pearray(:,:) =>NULL()
     integer                 :: pe            !PE to which this domain is assigned
     integer                 :: pos           !position of this PE within link list, i.e domain%list(pos)%pe = pe
     integer                 :: fold
     logical                 :: overlap(4)
     logical                 :: symmetry      !indicate the domain is symmetric or non-symmetric.
     type(overlapList)       :: send(8)
     type(overlapList)       :: recv(8)
     type(domain2D), pointer :: T =>NULL()    ! T-cell domain
     type(domain2D), pointer :: E =>NULL()    ! E-cell domain
     type(domain2D), pointer :: N =>NULL()    ! N-cell domain
     type(domain2D), pointer :: C =>NULL()    ! C-cell domain
     integer                 :: position      ! position of the domain, CENTER, EAST, NORTH, CORNER
     integer                 :: tile_id       ! tile number on current pe
     integer                 :: ntiles        ! number of tiles on the mosaic
     integer                 :: ncontacts     ! number of contacts in the mosaic
     integer                 :: topology_type ! its value can be REGULAR or CUBIC_GRID
  end type domain2D     

  type DomainCommunicator2D
     private
     logical            :: initialized=.false.
     integer(LONG_KIND) :: id=-9999
     integer(LONG_KIND) :: l_addr  =-9999
     integer(LONG_KIND) :: l_addrx =-9999
     integer(LONG_KIND) :: l_addry =-9999
     type(domain2D), pointer :: domain     =>NULL()
     type(domain2D), pointer :: domain_in  =>NULL()
     type(domain2D), pointer :: domain_out =>NULL()
     type(overlapList), pointer :: send(:,:) => NULL()
     type(overlapList), pointer :: recv(:,:) => NULL()
     integer, dimension(:,:),   _ALLOCATABLE :: sendis _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: sendie _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: sendjs _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: sendje _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: S_msize _NULL
     logical, dimension(:,:),   _ALLOCATABLE :: do_thisS _NULL
     logical, dimension(:),     _ALLOCATABLE :: S_do_buf _NULL
     logical, dimension(:,:),   _ALLOCATABLE :: do_thisS2 _NULL
     logical, dimension(:),     _ALLOCATABLE :: S_do_buf2 _NULL
     logical, dimension(:,:),   _ALLOCATABLE :: do_thisS3 _NULL
     logical, dimension(:),     _ALLOCATABLE :: S_do_buf3 _NULL
     logical, dimension(:,:),   _ALLOCATABLE :: do_thisS4 _NULL
     logical, dimension(:),     _ALLOCATABLE :: S_do_buf4 _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: S_msize2 _NULL
     integer, dimension(:),     _ALLOCATABLE :: cto_pe  _NULL
     integer, dimension(:),     _ALLOCATABLE :: Rcaf_idx  _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: recvis _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: recvie _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: recvjs _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: recvje _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: R_msize _NULL
     logical, dimension(:,:),   _ALLOCATABLE :: do_thisR _NULL
     logical, dimension(:),     _ALLOCATABLE :: R_do_buf _NULL
     logical, dimension(:,:),   _ALLOCATABLE :: do_thisR2 _NULL
     logical, dimension(:),     _ALLOCATABLE :: R_do_buf2 _NULL
     logical, dimension(:,:),   _ALLOCATABLE :: do_thisR3 _NULL
     logical, dimension(:),     _ALLOCATABLE :: R_do_buf3 _NULL
     logical, dimension(:,:),   _ALLOCATABLE :: do_thisR4 _NULL
     logical, dimension(:),     _ALLOCATABLE :: R_do_buf4 _NULL
     integer, dimension(:,:),   _ALLOCATABLE :: R_msize2 _NULL
     integer, dimension(:),     _ALLOCATABLE :: cfrom_pe  _NULL
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
     type(DomainCommunicator2D), pointer :: y_comm =>NULL() ! domain update information for y-component
     logical                             :: staggered ! true means x and y-component are staggered grid.
     integer                             :: position        ! data location. T, E, C, or N.
  end type DomainCommunicator2D

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

!***********************************************************************
!
!     module variables 
!
!***********************************************************************
  integer             :: pe
  logical             :: module_is_initialized = .false.
  logical             :: debug                 = .FALSE.
  logical             :: verbose=.FALSE.
  logical             :: debug_gsm=.false.
  integer             :: mpp_domains_stack_size=0
  integer             :: mpp_domains_stack_hwm=0
  type(domain1D),save :: NULL_DOMAIN1D
  type(domain2D),save :: NULL_DOMAIN2D

  !-------- The following variables are used in mpp_domains_comm.h
  
  integer, parameter :: MAX_ADDRS=512
  integer(LONG_KIND),dimension(MAX_ADDRS),save :: addrs_sorted=-9999  ! list of sorted local addrs
  integer,           dimension(-1:MAX_ADDRS),save :: addrs_idx=-9999  ! idx of addr assoicated w/ d_comm
  integer,           dimension(MAX_ADDRS),save :: a_salvage=-9999     ! freed idx list of addr
  integer,                                save :: a_sort_len=0        ! len sorted memory list
  integer,                                save :: n_addrs=0           ! num memory addresses used

  integer(LONG_KIND), parameter :: ADDR2_BASE=Z'0000000000010000'
  integer, parameter :: MAX_ADDRS2=128
  integer(LONG_KIND),dimension(MAX_ADDRS2),save :: addrs2_sorted=-9999  ! list of sorted local addrs
  integer,           dimension(-1:MAX_ADDRS2),save :: addrs2_idx=-9999  ! idx of addr2 assoicated w/ d_comm
  integer,           dimension(MAX_ADDRS2),save :: a2_salvage=-9999     ! freed indices of addr2
  integer,                                 save :: a2_sort_len=0        ! len sorted memory list
  integer,                                 save :: n_addrs2=0           ! num memory addresses used

  integer, parameter :: MAX_DOM_IDS=128
  integer(LONG_KIND),dimension(MAX_DOM_IDS),save :: ids_sorted=-9999 ! list of sorted domain identifiers
  integer,           dimension(-1:MAX_DOM_IDS),save :: ids_idx=-9999 ! idx of d_comm associated w/ sorted addr
  integer,                                  save :: i_sort_len=0     ! len sorted domain ids list
  integer,                                  save :: n_ids=0          ! num domain ids used (=i_sort_len; dom ids never removed)

  integer, parameter :: MAX_FIELDS=1024
  integer(LONG_KIND),        dimension(MAX_FIELDS),save           :: dcKey_sorted=-9999  ! list of sorted local addrs
  !     Not sure why static d_comm fails during deallocation of derived type members; allocatable works
  !     type(DomainCommunicator2D),dimension(MAX_FIELDS),save,target    :: d_comm              ! domain communicators
  type(DomainCommunicator2D),dimension(:),allocatable,save,target :: d_comm              ! domain communicators
  integer,                   dimension(-1:MAX_FIELDS),save           :: d_comm_idx=-9999 ! idx of d_comm associated w/ sorted addr
  integer,                   dimension(MAX_FIELDS),save           :: dc_salvage=-9999    ! freed indices of d_comm
  integer,                                         save           :: dc_sort_len=0       ! len sorted comm keys (=num active communicators)
  integer,                                         save           :: n_comm=0            ! num communicators used

  !     integer(LONG_KIND), parameter :: GT_BASE=2**8
  integer(LONG_KIND), parameter :: GT_BASE=Z'0000000000000100'  ! Workaround for 64bit int init problem

  !     integer(LONG_KIND), parameter :: KE_BASE=2**48
  integer(LONG_KIND), parameter :: KE_BASE=Z'0001000000000000'  ! Workaround for 64bit int init problem


  !--- the following variables are used in mpp_domains_misc.h
  logical :: domain_clocks_on=.FALSE.
  integer :: send_clock=0, recv_clock=0, unpk_clock=0
  integer :: wait_clock=0, pack_clock=0, pack_loop_clock=0

  !***********************************************************************
  !
  !           public interface for mpp_domains_comm.h
  !
  !***********************************************************************

#ifdef use_CAF
  interface mpp_associate_caf_field
     !     module procedure associate_caf_field_r8_1d
     !     module procedure associate_caf_field_c8_1d
     !     module procedure associate_caf_field_i8_1d
     !     module procedure associate_caf_field_l8_1d
     !     module procedure associate_caf_field_r4_1d
     !     module procedure associate_caf_field_c4_1d
     !     module procedure associate_caf_field_i4_1d
     !     module procedure associate_caf_field_l4_1d
     module procedure associate_caf_field_r8_3d
     module procedure associate_caf_field_c8_3d
     module procedure associate_caf_field_i8_3d
     module procedure associate_caf_field_l8_3d
     module procedure associate_caf_field_r4_3d
     module procedure associate_caf_field_c4_3d
     module procedure associate_caf_field_i4_3d
     module procedure associate_caf_field_l4_3d
  end interface
#endif

!***********************************************************************
!
!         public interface from mpp_domains_define.h
!
!***********************************************************************

  ! <INTERFACE NAME="mpp_define_layout">
  !  <OVERVIEW>
  !    Retrieve layout associated with a domain decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Given a global 2D domain and the number of divisions in the
  !    decomposition (<TT>ndivs</TT>: usually the PE count unless some
  !    domains are masked) this calls returns a 2D domain layout.
  !    
  !    By default, <TT>mpp_define_layout</TT> will attempt to divide the
  !    2D index space into domains that maintain the aspect ratio of the
  !    global domain. If this cannot be done, the algorithm favours domains
  !    that are longer in <TT>x</TT> than <TT>y</TT>, a preference that could
  !    improve vector performance.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_define_layout( global_indices, ndivs, layout )
  !  </TEMPLATE>
  !  <IN NAME="global_indices"></IN>
  !  <IN NAME="ndivs"></IN>
  !  <OUT NAME="layout"></OUT>
  ! </INTERFACE>

  interface mpp_define_layout
     module procedure mpp_define_layout2D
  end interface


  ! <INTERFACE NAME="mpp_define_domains">

  !   <OVERVIEW>
  !     Set up a domain decomposition.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     There are two forms for the <TT>mpp_define_domains</TT> call. The 2D
  !     version is generally to be used but is built by repeated calls to the
  !     1D version, also provided.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call mpp_define_domains( global_indices, ndivs, domain, &
  !                                   pelist, flags, halo, extent, maskmap )
  !   </TEMPLATE>
  !  <TEMPLATE>
  !    call mpp_define_domains( global_indices, layout, domain, pelist, &
  !                                   xflags, yflags, xhalo, yhalo,           &
  !                                   xextent, yextent, maskmap, name )
  !  </TEMPLATE>
  !   <IN NAME="global_indices" >
  !     Defines the global domain.
  !   </IN>
  !   <IN NAME="ndivs">
  !     Is the number of domain divisions required.
  !   </IN>
  !   <INOUT NAME="domain">
  !     Holds the resulting domain decomposition.
  !   </INOUT>
  !   <IN NAME="pelist">
  !     List of PEs to which the domains are to be assigned.
  !   </IN>
  !   <IN NAME="flags">
  !      An optional flag to pass additional information
  !      about the desired domain topology. Useful flags in a 1D decomposition
  !      include <TT>GLOBAL_DATA_DOMAIN</TT> and
  !      <TT>CYCLIC_GLOBAL_DOMAIN</TT>. Flags are integers: multiple flags may
  !      be added together. The flag values are public parameters available by
  !      use association.
  !   </IN>
  !   <IN NAME="halo">
  !     Width of the halo.
  !   </IN>
  !   <IN NAME="extent">
  !      Normally <TT>mpp_define_domains</TT> attempts
  !      an even division of the global domain across <TT>ndivs</TT>
  !      domains. The <TT>extent</TT> array can be used by the user to pass a
  !      custom domain division. The <TT>extent</TT> array has <TT>ndivs</TT>
  !      elements and holds the compute domain widths, which should add up to
  !      cover the global domain exactly.
  !   </IN>
  !   <IN NAME="maskmap">
  !     Some divisions may be masked
  !     (<TT>maskmap=.FALSE.</TT>) to exclude them from the computation (e.g
  !     for ocean model domains that are all land). The <TT>maskmap</TT> array
  !     is dimensioned <TT>ndivs</TT> and contains <TT>.TRUE.</TT> values for
  !     any domain that must be <I>included</I> in the computation (default
  !     all). The <TT>pelist</TT> array length should match the number of
  !     domains included in the computation.
  !    </IN>   

  !  <IN NAME="layout"></IN>
  !  <IN NAME="xflags, yflags"></IN>
  !  <IN NAME="xhalo, yhalo"></IN>
  !  <IN NAME="xextent, yextent"></IN>
  !  <IN NAME="name" ></IN>

  !  <NOTE>    
  !    For example:
  !    
  !    <PRE>
  !    call mpp_define_domains( (/1,100/), 10, domain, &
  !         flags=GLOBAL_DATA_DOMAIN+CYCLIC_GLOBAL_DOMAIN, halo=2 )
  !    </PRE>
  !    
  !    defines 10 compute domains spanning the range [1,100] of the global
  !    domain. The compute domains are non-overlapping blocks of 10. All the data
  !    domains are global, and with a halo of 2 span the range [-1:102]. And
  !    since the global domain has been declared to be cyclic,
  !    <TT>domain(9)%next => domain(0)</TT> and <TT>domain(0)%prev =>
  !    domain(9)</TT>. A field is allocated on the data domain, and computations proceed on
  !    the compute domain. A call to <LINK
  !    SRC="#mpp_update_domains"><TT>mpp_update_domains</TT></LINK> would fill in
  !    the values in the halo region:

  !    <PRE>
  !    call mpp_get_data_domain( domain, isd, ied ) !returns -1 and 102
  !    call mpp_get_compute_domain( domain, is, ie ) !returns (1,10) on PE 0 ...
  !    allocate( a(isd:ied) )
  !    do i = is,ie
  !       a(i) = &lt;perform computations&gt;
  !    end do
  !    call mpp_update_domains( a, domain )
  !    </PRE>

  !    The call to <TT>mpp_update_domains</TT> fills in the regions outside
  !    the compute domain. Since the global domain is cyclic, the values at
  !    <TT>i=(-1,0)</TT> are the same as at <TT>i=(99,100)</TT>; and
  !    <TT>i=(101,102)</TT> are the same as <TT>i=(1,2)</TT>.
  !    
  !    The 2D version is just an extension of this syntax to two
  !    dimensions.
  !
  !    The 2D version of the above should generally be used in
  !    codes, including 1D-decomposed ones, if there is a possibility of
  !    future evolution toward 2D decomposition. The arguments are similar to
  !    the 1D case, except that now we have optional arguments
  !    <TT>flags</TT>, <TT>halo</TT>, <TT>extent</TT> and <TT>maskmap</TT>
  !    along two axes.
  !    
  !    <TT>flags</TT> can now take an additional possible value to fold
  !    one or more edges. This is done by using flags
  !    <TT>FOLD_WEST_EDGE</TT>, <TT>FOLD_EAST_EDGE</TT>,
  !    <TT>FOLD_SOUTH_EDGE</TT> or <TT>FOLD_NORTH_EDGE</TT>. When a fold
  !    exists (e.g cylindrical domain), vector fields reverse sign upon
  !    crossing the fold. This parity reversal is performed only in the
  !    vector version of <LINK
  !    SRC="#mpp_update_domains"><TT>mpp_update_domains</TT></LINK>. In
  !    addition, shift operations may need to be applied to vector fields on
  !    staggered grids, also described in the vector interface to
  !    <TT>mpp_update_domains</TT>.
  !    
  !    <TT>name</TT> is the name associated with the decomposition,
  !    e.g <TT>'Ocean model'</TT>. If this argument is present,
  !    <TT>mpp_define_domains</TT> will print the domain decomposition
  !    generated to <TT>stdlog</TT>.
  !    
  !    Examples:
  !    
  !    <PRE>
  !    call mpp_define_domains( (/1,100,1,100/), (/2,2/), domain, xhalo=1 )
  !    </PRE>
  !    
  !    will create the following domain layout:
  !    <PRE>
  !                   |---------|-----------|-----------|-------------|
  !                   |domain(1)|domain(2)  |domain(3)  |domain(4)    |
  !    |--------------|---------|-----------|-----------|-------------|
  !    |Compute domain|1,50,1,50|51,100,1,50|1,50,51,100|51,100,51,100|
  !    |--------------|---------|-----------|-----------|-------------|
  !    |Data domain   |0,51,1,50|50,101,1,50|0,51,51,100|50,101,51,100|
  !    |--------------|---------|-----------|-----------|-------------|
  !    </PRE>
  !    
  !    Again, we allocate arrays on the data domain, perform computations
  !    on the compute domain, and call <TT>mpp_update_domains</TT> to update
  !    the halo region.
  !    
  !    If we wished to perfom a 1D decomposition along <TT>Y</TT>
  !    on the same global domain, we could use:

  !    <PRE>
  !    call mpp_define_domains( (/1,100,1,100/), layout=(/4,1/), domain, xhalo=1 )
  !    </PRE>

  !    This will create the following domain layout:
  !    <PRE>
  !                   |----------|-----------|-----------|------------|
  !                   |domain(1) |domain(2)  |domain(3)  |domain(4)   |
  !    |--------------|----------|-----------|-----------|------------|
  !    |Compute domain|1,100,1,25|1,100,26,50|1,100,51,75|1,100,76,100|
  !    |--------------|----------|-----------|-----------|------------|
  !    |Data domain   |0,101,1,25|0,101,26,50|0,101,51,75|1,101,76,100|
  !    |--------------|----------|-----------|-----------|------------|
  !    </PRE>
  !   </NOTE>
  ! </INTERFACE>
  interface mpp_define_domains
     module procedure mpp_define_domains1D
     module procedure mpp_define_domains2D
  end interface

! <INTERFACE NAME="mpp_modify_domain">
!   <OVERVIEW>
!     modifies the extents (compute, data and global) of domain
!   </OVERVIEW>
!   <IN NAME="domain_in">
!     The source domain.
!   </IN>
!   <IN NAME="halo">
!     Halo size of the returned 1D doamin. Default value is 0.
!   </IN>
!   <IN NAME="cbegin,cend">
!    Axis specifications associated with the compute domain of the returned 1D domain.
!   </IN>
!   <IN NAME="gbegin,gend">
!    Axis specifications associated with the global domain of the returned 1D domain.
!   </IN>
!   <IN NAME="isc,iec">
!    Zonal axis specifications associated with the compute domain of the returned 2D domain.
!   </IN>
!   <IN NAME="jsc,jec">
!    Meridinal axis specifications associated with the compute domain of the returned 2D domain.
!   </IN>
!   <IN NAME="isg,ieg">
!    Zonal axis specifications associated with the global domain of the returned 2D domain.
!   </IN>
!   <IN NAME="jsg,jeg">
!    Meridinal axis specifications associated with the global domain of the returned 2D domain.
!   </IN>
!   <IN NAME="xhalo,yhalo">
!     Halo size of the returned 2D doamin. Default value is 0.
!   </IN>
!   <INOUT NAME="domain_out">
!     The returned domain.
!   </INOUT>

! </INTERFACE>

  interface mpp_modify_domain
     module procedure mpp_modify_domain1D
     module procedure mpp_modify_domain2D
  end interface


!***********************************************************************
!
!        public interface from mpp_domains_misc.h
!
!***********************************************************************

! <INTERFACE NAME="mpp_update_domains">
!  <OVERVIEW>
!     Halo updates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_update_domains</TT> is used to perform a halo update of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>complex</TT>, <TT>integer</TT>, <TT>logical</TT> or <TT>real</TT>;
!    of 4-byte or 8-byte kind; of rank up to 5. The vector version (with
!    two input data fields) is only present for <TT>real</TT> types.
!    
!    For 2D domain updates, if there are halos present along both
!    <TT>x</TT> and <TT>y</TT>, we can choose to update one only, by
!    specifying <TT>flags=XUPDATE</TT> or <TT>flags=YUPDATE</TT>. In
!    addition, one-sided updates can be performed by setting <TT>flags</TT>
!    to any combination of <TT>WUPDATE</TT>, <TT>EUPDATE</TT>,
!    <TT>SUPDATE</TT> and <TT>NUPDATE</TT>, to update the west, east, north
!    and south halos respectively. Any combination of halos may be used by
!    adding the requisite flags, e.g: <TT>flags=XUPDATE+SUPDATE</TT> or
!    <TT>flags=EUPDATE+WUPDATE+SUPDATE</TT> will update the east, west and
!    south halos.
!    
!    If a call to <TT>mpp_update_domains</TT> involves at least one E-W
!    halo and one N-S halo, the corners involved will also be updated, i.e,
!    in the example above, the SE and SW corners will be updated.
!    
!    If <TT>flags</TT> is not supplied, that is
!    equivalent to <TT>flags=XUPDATE+YUPDATE</TT>.
!    
!    The vector version is passed the <TT>x</TT> and <TT>y</TT>
!    components of a vector field in tandem, and both are updated upon
!    return. They are passed together to treat parity issues on various
!    grids. For example, on a cubic sphere projection, the <TT>x</TT> and
!    <TT>y</TT> components may be interchanged when passing from an
!    equatorial cube face to a polar face. For grids with folds, vector
!    components change sign on crossing the fold.  Paired scalar quantities
!    can also be passed with the vector version if flags=SCALAR_PAIR, in which
!    case components are appropriately interchanged, but signs are not.
!    
!    Special treatment at boundaries such as folds is also required for
!    staggered grids. The following types of staggered grids are
!    recognized:
!    
!    1) <TT>AGRID</TT>: values are at grid centers.<BR/>
!    2) <TT>BGRID_NE</TT>: vector fields are at the NE vertex of a grid
!    cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT> are
!    actually at (i+&#189;,j+&#189;) with respect to the grid centers.<BR/>
!    3) <TT>BGRID_SW</TT>: vector fields are at the SW vertex of a grid
!    cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT> are
!    actually at (i-&#189;,j-&#189;) with respect to the grid centers.<BR/>
!    4) <TT>CGRID_NE</TT>: vector fields are at the N and E faces of a
!    grid cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT>
!    are actually at (i+&#189;,j) and (i,j+&#189;) with respect to the
!    grid centers.<BR/>
!    5) <TT>CGRID_SW</TT>: vector fields are at the S and W faces of a
!    grid cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT>
!    are actually at (i-&#189;,j) and (i,j-&#189;) with respect to the
!    grid centers.
!
!    The gridtypes listed above are all available by use association as
!    integer parameters. The scalar version of <TT>mpp_update_domains</TT>
!    assumes that the values of a scalar field are always at <TT>AGRID</TT>
!    locations, and no special boundary treatment is required. If vector
!    fields are at staggered locations, the optional argument
!    <TT>gridtype</TT> must be appropriately set for correct treatment at
!    boundaries.
!    
!    It is safe to apply vector field updates to the appropriate arrays
!    irrespective of the domain topology: if the topology requires no
!    special treatment of vector fields, specifying <TT>gridtype</TT> will
!    do no harm.
!
!    <TT>mpp_update_domains</TT> internally buffers the date being sent
!    and received into single messages for efficiency. A turnable internal
!    buffer area in memory is provided for this purpose by
!    <TT>mpp_domains_mod</TT>. The size of this buffer area can be set by
!    the user by calling <LINK SRC="mpp_domains.html#mpp_domains_set_stack_size">
!    <TT>mpp_domains_set_stack_size</TT></LINK>.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_update_domains( field, domain, flags )
!  </TEMPLATE>
!  <TEMPLATE>
!    call mpp_update_domains( fieldx, fieldy, domain, flags, gridtype )
!  </TEMPLATE>
! </INTERFACE>
  interface mpp_update_domains
     module procedure mpp_update_domain2D_r8_2d
     module procedure mpp_update_domain2D_r8_3d
     module procedure mpp_update_domain2D_r8_4d
     module procedure mpp_update_domain2D_r8_5d
     module procedure mpp_update_domain2D_r8_2dv
     module procedure mpp_update_domain2D_r8_3dv
     module procedure mpp_update_domain2D_r8_4dv
     module procedure mpp_update_domain2D_r8_5dv
     module procedure mpp_update_domain2D_c8_2d
     module procedure mpp_update_domain2D_c8_3d
     module procedure mpp_update_domain2D_c8_4d
     module procedure mpp_update_domain2D_c8_5d
#ifndef no_8byte_integers
     module procedure mpp_update_domain2D_i8_2d
     module procedure mpp_update_domain2D_i8_3d
     module procedure mpp_update_domain2D_i8_4d
     module procedure mpp_update_domain2D_i8_5d
     module procedure mpp_update_domain2D_l8_2d
     module procedure mpp_update_domain2D_l8_3d
     module procedure mpp_update_domain2D_l8_4d
     module procedure mpp_update_domain2D_l8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_update_domain2D_r4_2d
     module procedure mpp_update_domain2D_r4_3d
     module procedure mpp_update_domain2D_r4_4d
     module procedure mpp_update_domain2D_r4_5d
     module procedure mpp_update_domain2D_c4_2d
     module procedure mpp_update_domain2D_c4_3d
     module procedure mpp_update_domain2D_c4_4d
     module procedure mpp_update_domain2D_c4_5d
     module procedure mpp_update_domain2D_r4_2dv
     module procedure mpp_update_domain2D_r4_3dv
     module procedure mpp_update_domain2D_r4_4dv
     module procedure mpp_update_domain2D_r4_5dv
#endif
     module procedure mpp_update_domain2D_i4_2d
     module procedure mpp_update_domain2D_i4_3d
     module procedure mpp_update_domain2D_i4_4d
     module procedure mpp_update_domain2D_i4_5d
     module procedure mpp_update_domain2D_l4_2d
     module procedure mpp_update_domain2D_l4_3d
     module procedure mpp_update_domain2D_l4_4d
     module procedure mpp_update_domain2D_l4_5d
  end interface


  interface mpp_do_update
     module procedure mpp_do_update_new_r8_3d
     module procedure mpp_do_update_old_r8_3d
     module procedure mpp_do_update_new_r8_3dv
     module procedure mpp_do_update_old_r8_3dv
     module procedure mpp_do_update_new_c8_3d
     module procedure mpp_do_update_old_c8_3d
#ifndef no_8byte_integers
     module procedure mpp_do_update_new_i8_3d
     module procedure mpp_do_update_old_i8_3d
     module procedure mpp_do_update_new_l8_3d
     module procedure mpp_do_update_old_l8_3d
#endif
#ifndef no_4byte_reals
     module procedure mpp_do_update_new_r4_3d
     module procedure mpp_do_update_old_r4_3d
     module procedure mpp_do_update_new_r4_3dv
     module procedure mpp_do_update_old_r4_3dv
     module procedure mpp_do_update_new_c4_3d
     module procedure mpp_do_update_old_c4_3d
#endif
     module procedure mpp_do_update_new_i4_3d
     module procedure mpp_do_update_old_i4_3d
     module procedure mpp_do_update_new_l4_3d
     module procedure mpp_do_update_old_l4_3d
  end interface

! <INTERFACE NAME="mpp_redistribute">
!  <OVERVIEW>
!    Reorganization of distributed global arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_redistribute</TT> is used to reorganize a distributed
!    array.  <TT>MPP_TYPE_</TT> can be of type <TT>integer</TT>,
!    <TT>complex</TT>, or <TT>real</TT>; of 4-byte or 8-byte kind; of rank
!    up to 5.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_redistribute( domain_in, field_in, domain_out, field_out )
!  </TEMPLATE>
!  <IN NAME="field_in" TYPE="MPP_TYPE_">
!    <TT>field_in</TT> is dimensioned on the data domain of <TT>domain_in</TT>.
!  </IN>
!  <OUT NAME="field_out" TYPE="MPP_TYPE_">
!    <TT>field_out</TT> on the data domain of <TT>domain_out</TT>.
!  </OUT>
! </INTERFACE>
  interface mpp_redistribute
     module procedure mpp_redistribute_r8_2D
     module procedure mpp_redistribute_r8_3D
     module procedure mpp_redistribute_r8_4D
     module procedure mpp_redistribute_r8_5D
     module procedure mpp_redistribute_c8_2D
     module procedure mpp_redistribute_c8_3D
     module procedure mpp_redistribute_c8_4D
     module procedure mpp_redistribute_c8_5D
#ifndef no_8byte_integers
     module procedure mpp_redistribute_i8_2D
     module procedure mpp_redistribute_i8_3D
     module procedure mpp_redistribute_i8_4D
     module procedure mpp_redistribute_i8_5D
     module procedure mpp_redistribute_l8_2D
     module procedure mpp_redistribute_l8_3D
     module procedure mpp_redistribute_l8_4D
     module procedure mpp_redistribute_l8_5D
#endif
#ifndef no_4byte_reals
     module procedure mpp_redistribute_r4_2D
     module procedure mpp_redistribute_r4_3D
     module procedure mpp_redistribute_r4_4D
     module procedure mpp_redistribute_r4_5D
     module procedure mpp_redistribute_c4_2D
     module procedure mpp_redistribute_c4_3D
     module procedure mpp_redistribute_c4_4D
     module procedure mpp_redistribute_c4_5D
#endif
     module procedure mpp_redistribute_i4_2D
     module procedure mpp_redistribute_i4_3D
     module procedure mpp_redistribute_i4_4D
     module procedure mpp_redistribute_i4_5D
     module procedure mpp_redistribute_l4_2D
     module procedure mpp_redistribute_l4_3D
     module procedure mpp_redistribute_l4_4D
     module procedure mpp_redistribute_l4_5D
  end interface

  interface mpp_do_redistribute
     module procedure mpp_do_redistribute_new_r8_3D
     module procedure mpp_do_redistribute_old_r8_3D
     module procedure mpp_do_redistribute_new_c8_3D
     module procedure mpp_do_redistribute_old_c8_3D
#ifndef no_8byte_integers
     module procedure mpp_do_redistribute_new_i8_3D
     module procedure mpp_do_redistribute_old_i8_3D
     module procedure mpp_do_redistribute_new_l8_3D
     module procedure mpp_do_redistribute_old_l8_3D
#endif
#ifndef no_4byte_reals
     module procedure mpp_do_redistribute_new_r4_3D
     module procedure mpp_do_redistribute_old_r4_3D
     module procedure mpp_do_redistribute_new_c4_3D
     module procedure mpp_do_redistribute_old_c4_3D
#endif
     module procedure mpp_do_redistribute_new_i4_3D
     module procedure mpp_do_redistribute_old_i4_3D
     module procedure mpp_do_redistribute_new_l4_3D
     module procedure mpp_do_redistribute_old_l4_3D
  end interface


! <INTERFACE NAME="mpp_check_field">
!   <OVERVIEW>
!     Parallel checking between two ensembles which run
!     on different set pes at the same time.
!   </OVERVIEW>
!   <DESCRIPTION>
!     There are two forms for the <TT>mpp_check_field</TT> call. The 2D
!     version is generally to be used and 3D version is  built by repeated calls to the
!     2D version.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_check_field(field_in, pelist1, pelist2, domain, mesg, &
!                                w_halo, s_halo, e_halo, n_halo, force_abort  )
!   </TEMPLATE>
!   <IN NAME="field_in" >
!     Field to be checked
!   </IN>
!   <IN NAME="pelist1, pelist2">
!     Pelist of the two ensembles to be compared
!   </IN>
!   <IN NAME="domain">
!     Domain of current pe
!   </IN>
!   <IN NAME="mesg" >
!     Message to be printed out
!   </IN>
!   <IN NAME="w_halo, s_halo, e_halo, n_halo">
!     Halo size to be checked. Default value is 0.
!   </IN>
!   <IN NAME="force_abort">
!     When true, abort program when any difference found. Default value is false.
!   </IN>
! </INTERFACE>

  interface mpp_check_field
     module procedure mpp_check_field_2D
     module procedure mpp_check_field_3D
  end interface

!***********************************************************************
!
!         public interface from mpp_domains_reduce.h
!
!***********************************************************************

! <INTERFACE NAME="mpp_global_field">
!  <OVERVIEW>
!    Fill in a global array from domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_field</TT> is used to get an entire
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>complex</TT>, <TT>integer</TT>, <TT>logical</TT> or <TT>real</TT>;
!    of 4-byte or 8-byte kind; of rank up to 5.
!    
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_field</TT>, and each will have a complete global field
!    at the end. Please note that a global array of rank 3 or higher could
!    occupy a lot of memory.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_global_field( domain, local, global, flags )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="local" TYPE="MPP_TYPE_">
!    <TT>local</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <OUT NAME="global" TYPE="MPP_TYPE_">
!    <TT>global</TT> is dimensioned on the corresponding global domain.
!  </OUT>
!  <IN NAME="flags" TYPE="integer">
!    <TT>flags</TT> can be given the value <TT>XONLY</TT> or
!    <TT>YONLY</TT>, to specify a globalization on one axis only.
!  </IN>
! </INTERFACE>
  interface mpp_global_field
     module procedure mpp_global_field2D_r8_2d
     module procedure mpp_global_field2D_r8_3d
     module procedure mpp_global_field2D_r8_4d
     module procedure mpp_global_field2D_r8_5d
     module procedure mpp_global_field2D_c8_2d
     module procedure mpp_global_field2D_c8_3d
     module procedure mpp_global_field2D_c8_4d
     module procedure mpp_global_field2D_c8_5d
#ifndef no_8byte_integers
     module procedure mpp_global_field2D_i8_2d
     module procedure mpp_global_field2D_i8_3d
     module procedure mpp_global_field2D_i8_4d
     module procedure mpp_global_field2D_i8_5d
     module procedure mpp_global_field2D_l8_2d
     module procedure mpp_global_field2D_l8_3d
     module procedure mpp_global_field2D_l8_4d
     module procedure mpp_global_field2D_l8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_global_field2D_r4_2d
     module procedure mpp_global_field2D_r4_3d
     module procedure mpp_global_field2D_r4_4d
     module procedure mpp_global_field2D_r4_5d
     module procedure mpp_global_field2D_c4_2d
     module procedure mpp_global_field2D_c4_3d
     module procedure mpp_global_field2D_c4_4d
     module procedure mpp_global_field2D_c4_5d
#endif
     module procedure mpp_global_field2D_i4_2d
     module procedure mpp_global_field2D_i4_3d
     module procedure mpp_global_field2D_i4_4d
     module procedure mpp_global_field2D_i4_5d
     module procedure mpp_global_field2D_l4_2d
     module procedure mpp_global_field2D_l4_3d
     module procedure mpp_global_field2D_l4_4d
     module procedure mpp_global_field2D_l4_5d
  end interface

  interface mpp_do_global_field
     module procedure mpp_do_global_field2Dnew_r8_3d
     module procedure mpp_do_global_field2Dold_r8_3d
     module procedure mpp_do_global_field2Dnew_c8_3d
     module procedure mpp_do_global_field2Dold_c8_3d
#ifndef no_8byte_integers
     module procedure mpp_do_global_field2Dnew_i8_3d
     module procedure mpp_do_global_field2Dold_i8_3d
     module procedure mpp_do_global_field2Dnew_l8_3d
     module procedure mpp_do_global_field2Dold_l8_3d
#endif
#ifndef no_4byte_reals
     module procedure mpp_do_global_field2Dnew_r4_3d
     module procedure mpp_do_global_field2Dold_r4_3d
     module procedure mpp_do_global_field2Dnew_c4_3d
     module procedure mpp_do_global_field2Dold_c4_3d
#endif
     module procedure mpp_do_global_field2Dnew_i4_3d
     module procedure mpp_do_global_field2Dold_i4_3d
     module procedure mpp_do_global_field2Dnew_l4_3d
     module procedure mpp_do_global_field2Dold_l4_3d
  end interface

! <INTERFACE NAME="mpp_global_max">
!  <OVERVIEW>
!    Global max/min of domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_max</TT> is used to get the maximum value of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>integer</TT> or <TT>real</TT>; of 4-byte or 8-byte kind; of rank
!    up to 5. The dimension of <TT>locus</TT> must equal the rank of
!    <TT>field</TT>.
!    
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_max</TT>, and each will have the result upon exit.
!    
!    The function <TT>mpp_global_min</TT>, with an identical syntax. is
!    also available.
!  </DESCRIPTION>
!  <TEMPLATE>
!    mpp_global_max( domain, field, locus )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="field" TYPE="MPP_TYPE_">  
!    <TT>field</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <OUT NAME="locus" TYPE="integer" DIM="(:)">
!    <TT>locus</TT>, if present, can be used to retrieve the location of
!    the maximum (as in the <TT>MAXLOC</TT> intrinsic of f90).
!  </OUT>
! </INTERFACE>

  interface mpp_global_max
     module procedure mpp_global_max_r8_2d
     module procedure mpp_global_max_r8_3d
     module procedure mpp_global_max_r8_4d
     module procedure mpp_global_max_r8_5d
#ifndef no_4byte_reals
     module procedure mpp_global_max_r4_2d
     module procedure mpp_global_max_r4_3d
     module procedure mpp_global_max_r4_4d
     module procedure mpp_global_max_r4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_max_i8_2d
     module procedure mpp_global_max_i8_3d
     module procedure mpp_global_max_i8_4d
     module procedure mpp_global_max_i8_5d
#endif
     module procedure mpp_global_max_i4_2d
     module procedure mpp_global_max_i4_3d
     module procedure mpp_global_max_i4_4d
     module procedure mpp_global_max_i4_5d
  end interface

  interface mpp_global_min
     module procedure mpp_global_min_r8_2d
     module procedure mpp_global_min_r8_3d
     module procedure mpp_global_min_r8_4d
     module procedure mpp_global_min_r8_5d
#ifndef no_4byte_reals
     module procedure mpp_global_min_r4_2d
     module procedure mpp_global_min_r4_3d
     module procedure mpp_global_min_r4_4d
     module procedure mpp_global_min_r4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_min_i8_2d
     module procedure mpp_global_min_i8_3d
     module procedure mpp_global_min_i8_4d
     module procedure mpp_global_min_i8_5d
#endif
     module procedure mpp_global_min_i4_2d
     module procedure mpp_global_min_i4_3d
     module procedure mpp_global_min_i4_4d
     module procedure mpp_global_min_i4_5d
  end interface

! <INTERFACE NAME="mpp_global_sum">
!  <OVERVIEW>
!    Global sum of domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_sum</TT> is used to get the sum of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>integer</TT>, <TT>complex</TT>, or <TT>real</TT>; of 4-byte or
!    8-byte kind; of rank up to 5.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_global_sum( domain, field, flags )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="field" TYPE="MPP_TYPE_">
!    <TT>field</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <IN NAME="flags" TYPE="integer">
!    <TT>flags</TT>, if present, must have the value
!    <TT>BITWISE_EXACT_SUM</TT>. This produces a sum that is guaranteed to
!    produce the identical result irrespective of how the domain is
!    decomposed. This method does the sum first along the ranks beyond 2,
!    and then calls <LINK
!    SRC="#mpp_global_field"><TT>mpp_global_field</TT></LINK> to produce a
!    global 2D array which is then summed. The default method, which is
!    considerably faster, does a local sum followed by <LINK
!    SRC="mpp.html#mpp_sum"><TT>mpp_sum</TT></LINK> across the domain
!    decomposition.
!  </IN>
!  <NOTE>
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_sum</TT>, and each will have the result upon exit.
!  </NOTE>
! </INTERFACE>

  interface mpp_global_sum
     module procedure mpp_global_sum_r8_2d
     module procedure mpp_global_sum_r8_3d
     module procedure mpp_global_sum_r8_4d
     module procedure mpp_global_sum_r8_5d
     module procedure mpp_global_sum_c8_2d
     module procedure mpp_global_sum_c8_3d
     module procedure mpp_global_sum_c8_4d
     module procedure mpp_global_sum_c8_5d
#ifndef no_4byte_reals
     module procedure mpp_global_sum_r4_2d
     module procedure mpp_global_sum_r4_3d
     module procedure mpp_global_sum_r4_4d
     module procedure mpp_global_sum_r4_5d
     module procedure mpp_global_sum_c4_2d
     module procedure mpp_global_sum_c4_3d
     module procedure mpp_global_sum_c4_4d
     module procedure mpp_global_sum_c4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_sum_i8_2d
     module procedure mpp_global_sum_i8_3d
     module procedure mpp_global_sum_i8_4d
     module procedure mpp_global_sum_i8_5d
#endif
     module procedure mpp_global_sum_i4_2d
     module procedure mpp_global_sum_i4_3d
     module procedure mpp_global_sum_i4_4d
     module procedure mpp_global_sum_i4_5d
  end interface


!***********************************************************************
!
!            public interface from mpp_domain_util.h
!
!***********************************************************************

  ! <INTERFACE NAME="mpp_get_neighbor_pe">
  !  <OVERVIEW>
  !    Retrieve PE number of a neighboring domain.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Given a 1-D or 2-D domain decomposition, this call allows users to retrieve 
  !    the PE number of an adjacent PE-domain while taking into account that the 
  !    domain may have holes (masked) and/or have cyclic boundary conditions and/or a 
  !    folded edge. Which PE-domain will be retrived will depend on "direction": 
  !    +1 (right) or -1 (left) for a 1-D domain decomposition and either NORTH, SOUTH, 
  !    EAST, WEST, NORTH_EAST, SOUTH_EAST, SOUTH_WEST, or NORTH_WEST for a 2-D 
  !    decomposition. If no neighboring domain exists (masked domain), then the 
  !    returned "pe" value will be set to NULL_PE.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_neighbor_pe( domain1d, direction=+1   , pe)
  !    call mpp_get_neighbor_pe( domain2d, direction=NORTH, pe)
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_neighbor_pe
     module procedure mpp_get_neighbor_pe_1d
     module procedure mpp_get_neighbor_pe_2d
  end interface 
  ! <INTERFACE NAME="operator">
  !  <OVERVIEW>
  !    Equality/inequality operators for domaintypes.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The module provides public operators to check for
  !    equality/inequality of domaintypes, e.g:
  !    
  !    <PRE>
  !    type(domain1D) :: a, b
  !    type(domain2D) :: c, d
  !    ...
  !    if( a.NE.b )then
  !        ...
  !    end if
  !    if( c==d )then
  !        ...
  !    end if
  !    </PRE>
  !    
  !    Domains are considered equal if and only if the start and end
  !    indices of each of their component global, data and compute domains
  !    are equal.
  !  </DESCRIPTION>
  ! </INTERFACE>
  interface operator(.EQ.)
     module procedure mpp_domain1D_eq
     module procedure mpp_domain2D_eq
  end interface

  interface operator(.NE.)
     module procedure mpp_domain1D_ne
     module procedure mpp_domain2D_ne
  end interface

  ! <INTERFACE NAME="mpp_get_compute_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the compute domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the compute domains
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_compute_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_compute_domain
     module procedure mpp_get_compute_domain1D
     module procedure mpp_get_compute_domain2D
  end interface

  ! <INTERFACE NAME="mpp_get_compute_domains">
  !  <OVERVIEW>
  !    Retrieve the entire array of compute domain extents associated with a decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Retrieve the entire array of compute domain extents associated with a decomposition.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_compute_domains( domain, xbegin, xend, xsize, &
  !                                                ybegin, yend, ysize )
  !  </TEMPLATE>
  !  <IN NAME="domain" TYPE="type(domain2D)"></IN>
  !  <OUT NAME="xbegin,ybegin" TYPE="integer" DIM="(:)"></OUT>
  !  <OUT NAME="xend,yend" TYPE="integer" DIM="(:)"></OUT>
  !  <OUT NAME="xsize,ysize" TYPE="integer" DIM="(:)"></OUT>
  ! </INTERFACE>
  interface mpp_get_compute_domains
     module procedure mpp_get_compute_domains1D
     module procedure mpp_get_compute_domains2D
  end interface

  ! <INTERFACE NAME="mpp_get_data_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the data domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the data domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_data_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_data_domain
     module procedure mpp_get_data_domain1D
     module procedure mpp_get_data_domain2D
  end interface

  ! <INTERFACE NAME="mpp_get_global_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the global domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the global domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_global_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_global_domain
     module procedure mpp_get_global_domain1D
     module procedure mpp_get_global_domain2D
  end interface

  ! <INTERFACE NAME="mpp_set_compute_domain">
  !  <OVERVIEW>
  !    These routines set the axis specifications associated with the compute domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    set the axis specifications associated with the compute domains
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_set_compute_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_set_compute_domain
     module procedure mpp_set_compute_domain1D
     module procedure mpp_set_compute_domain2D
  end interface

  ! <INTERFACE NAME="mpp_set_data_domain">
  !  <OVERVIEW>
  !    These routines set the axis specifications associated with the data domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    set the axis specifications associated with the data domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_set_data_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_set_data_domain
     module procedure mpp_set_data_domain1D
     module procedure mpp_set_data_domain2D
  end interface

  ! <INTERFACE NAME="mpp_set_global_domain">
  !  <OVERVIEW>
  !    These routines set the axis specifications associated with the global domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    set the axis specifications associated with the global domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_set_global_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_set_global_domain
     module procedure mpp_set_global_domain1D
     module procedure mpp_set_global_domain2D
  end interface


  ! <INTERFACE NAME="mpp_get_pelist">
  !  <OVERVIEW>
  !    Retrieve list of PEs associated with a domain decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The 1D version of this call returns an array of the PEs assigned to this 1D domain
  !    decomposition. In addition the optional argument <TT>pos</TT> may be
  !    used to retrieve the 0-based position of the domain local to the
  !    calling PE, i.e <TT>domain%list(pos)%pe</TT> is the local PE,
  !    as returned by <LINK SRC="mpp.html#mpp_pe"><TT>mpp_pe()</TT></LINK>.
  !    The 2D version of this call is identical to 1D version.
  !  </DESCRIPTION>
  !  <IN NAME="domain"></IN>
  !  <OUT NAME="pelist"></OUT>
  !  <OUT NAME="pos"></OUT>
  ! </INTERFACE>
  interface mpp_get_pelist
     module procedure mpp_get_pelist1D
     module procedure mpp_get_pelist2D
  end interface

  ! <INTERFACE NAME="mpp_get_layout">
  !  <OVERVIEW>
  !    Retrieve layout associated with a domain decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The 1D version of this call returns the number of divisions that was assigned to this
  !    decomposition axis. The 2D version of this call returns an array of
  !    dimension 2 holding the results on two axes.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_layout( domain, layout )
  !  </TEMPLATE>
  !  <IN NAME="domain"></IN>
  !  <OUT NAME="layout"></OUT>
  ! </INTERFACE>
  interface mpp_get_layout
     module procedure mpp_get_layout1D
     module procedure mpp_get_layout2D
  end interface

  ! <INTERFACE NAME="mpp_nullify_domain_list">
  !  <OVERVIEW>
  !    nullify domain list.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Nullify domain list. This interface is needed in mpp_domains_test.
  !    1-D case can be added in if needed.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_nullify_domain_list( domain)
  !  </TEMPLATE>
  !  <INOUT NAME="domain"></INOUT>
  ! </INTERFACE>
  interface mpp_nullify_domain_list
     module procedure nullify_domain2d_list
  end interface  

  !--- version information variables
  character(len=128), public :: version= &
       '$Id: mpp_domains.F90,v 13.0 2006/03/28 21:42:16 fms Exp $'
  character(len=128), public :: tagname= &
       '$Name: memphis $'


contains

#include <mpp_domains_util.inc>
#include <mpp_domains_comm.inc>
#include <mpp_domains_define.inc>
#include <mpp_domains_misc.inc>
#include <mpp_domains_reduce.inc>


end module mpp_domains_mod

#ifdef test_mpp_domains
program mpp_domains_test
  use mpp_mod,         only : FATAL, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_node, mpp_root_pe, mpp_error, mpp_set_warn_level
  use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync
  use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod,         only : mpp_init, mpp_exit, mpp_chksum, stdout, stderr
  use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, FOLD_NORTH_EDGE, CGRID_NE
  use mpp_domains_mod, only : MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
  use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod, only : mpp_update_domains, mpp_check_field, mpp_redistribute
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains, mpp_modify_domain
  use mpp_domains_mod, only : mpp_get_neighbor_pe, mpp_define_mosaic, mpp_nullify_domain_list
  use mpp_domains_mod, only : NORTH, NORTH_EAST, EAST, SOUTH_EAST
  use mpp_domains_mod, only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST


  implicit none
#include <fms_platform.h>
  integer :: pe, npes
  integer :: nx=128, ny=128, nz=40, halo=2, stackmax=4000000
  real, dimension(:,:,:), allocatable :: global
  integer :: unit=7
  logical :: debug=.FALSE., opened
  logical :: check_parallel = .FALSE.  ! when check_parallel set to false,
                                       ! mpes should be equal to npes     
  integer :: mpes = 0
  integer :: xhalo =0, yhalo =0
  namelist / mpp_domains_nml / nx, ny, nz, halo, stackmax, debug, mpes, check_parallel, xhalo, yhalo
  integer :: i, j, k
  integer :: layout(2)
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: id

  call mpp_init()

  call mpp_set_warn_level(FATAL)
!possibly open a file called mpp_domains.nml
  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, status='OLD', file='mpp_domains.nml', err=10 )
  read( unit,mpp_domains_nml )
  close(unit)
10 continue
  
  pe = mpp_pe()
  npes = mpp_npes()
  
  if( debug )then
      call mpp_domains_init(MPP_DEBUG)
  else
      call mpp_domains_init(MPP_DOMAIN_TIME)
  end if
  call mpp_domains_set_stack_size(stackmax)
  
  if( pe.EQ.mpp_root_pe() )print '(a,6i4)', 'npes, mpes, nx, ny, nz, halo=', npes, mpes, nx, ny, nz, halo
  
  allocate( global(1-halo:nx+halo,1-halo:ny+halo,nz) )
  
!fill in global array: with k.iiijjj
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx
           global(i,j,k) = k + i*1e-3 + j*1e-6
        end do
     end do
  end do
     
  if( .not. check_parallel) then
      call test_modify_domain()
      call test_mosaic('Single tile')
      call test_mosaic('Cyclic mosaic')
      call test_mosaic('Cubic-grid')

      call test_halo_update( 'Simple' ) !includes global field, global sum tests
      call test_halo_update( 'Cyclic' )
      call test_halo_update( 'Folded' ) !includes vector field test
      call test_halo_update( 'Masked' ) !includes vector field test

      call test_halo_update( 'Simple symmetry' ) !includes global field, global sum tests
      call test_halo_update( 'Cyclic symmetry' )
      call test_halo_update( 'Folded symmetry' ) !includes vector field test
      !--- z1l: The following will not work due to symmetry and domain%x is cyclic.
      !--- Will solve this problem in the future if needed.
      ! call test_halo_update( 'Masked symmetry' ) !includes vector field test

      call test_global_field( 'Non-symmetry' )
      call test_global_field( 'Symmetry center' )
      call test_global_field( 'Symmetry corner' )
      call test_global_field( 'Symmetry east' )
      call test_global_field( 'Symmetry north' )

      call test_global_reduce( 'Simple')
      call test_global_reduce( 'Simple symmetry center')
      call test_global_reduce( 'Simple symmetry corner')
      call test_global_reduce( 'Simple symmetry east')
      call test_global_reduce( 'Simple symmetry north')
      call test_global_reduce( 'Cyclic symmetry center')
      call test_global_reduce( 'Cyclic symmetry corner')
      call test_global_reduce( 'Cyclic symmetry east')
      call test_global_reduce( 'Cyclic symmetry north')

      call test_redistribute( 'Complete pelist' )
      call test_redistribute( 'Overlap  pelist' )
      call test_redistribute( 'Disjoint pelist' )
  else
      call test_parallel( )
  endif

!Balaji adding openMP tests
  call test_openmp()
 
! Alewxander.Pletzer get_neighbor tests
  call test_get_neighbor_1d
  call test_get_neighbor_non_cyclic
  call test_get_neighbor_cyclic
  call test_get_neighbor_folded_north
  call test_get_neighbor_mask

  call mpp_domains_exit()
  call mpp_exit()
  
contains
  subroutine test_openmp()
#ifdef _OPENMP
    integer :: omp_get_num_thread, omp_get_max_threads, omp_get_thread_num
    real, allocatable :: a(:,:,:)
    type(domain2D) :: domain
    integer :: layout(2)
    integer :: i,j,k, jthr
    integer :: thrnum, maxthr
    integer(LONG_KIND) :: sum1, sum2
    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain )
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( a(isd:ied,jsd:jed,nz) )
    maxthr = omp_get_max_threads()
    write( stdout(),'(a,4i4)' )'pe,js,je,maxthr=', pe, js, je, maxthr
!    write( stderr(),'(a,2i4)' )'pe,mldid=', pe, mld_id()
    if( mod(je-js+1,maxthr).NE.0 ) &
         call mpp_error( FATAL, 'maxthr must divide domain (TEMPORARY).' )
    jthr = (je-js+1)/maxthr
!$OMP PARALLEL PRIVATE(i,j,k,thrnum)
    thrnum = omp_get_thread_num()
    write( stdout(),'(a,4i4)' )'pe,thrnum,js,je=', &
         pe, thrnum, js+thrnum*jthr,js+(thrnum+1)*jthr-1
    write( stdout(),'(a,3i4)' )'pe,thrnum,node=', pe, thrnum, mpp_node()
!!$OMP DO
    do k = 1,nz
!when omp DO is commented out, user must compute j loop limits
!with omp DO, let OMP figure it out
       do j = js+thrnum*jthr,js+(thrnum+1)*jthr-1
!       do j = js,je
          do i = is,ie
             a(i,j,k) = global(i,j,k)
          end do
       end do
    end do
!!$OMP END DO
!$OMP END PARALLEL
    sum1 = mpp_chksum( a(is:ie,js:je,:) )
    sum2 = mpp_chksum( global(is:ie,js:je,:) )
    if( sum1.EQ.sum2 )then
        call mpp_error( NOTE, 'OMP parallel test OK.' )
    else
        if( mpp_pe().EQ.mpp_root_pe() )write( stderr(),'(a,2z18)' )'OMP checksums: ', sum1, sum2
        call mpp_error( FATAL, 'OMP parallel test failed.' )
    end if
#endif
    return
  end subroutine test_openmp

  subroutine test_redistribute( type )
!test redistribute between two domains
    character(len=*), intent(in) :: type
    type(domain2D) :: domainx, domainy
    type(DomainCommunicator2D), pointer, save :: dch =>NULL()
    real, allocatable, dimension(:,:,:)       :: gcheck
    real, allocatable, dimension(:,:,:), save :: x, y
    real, allocatable, dimension(:,:,:), save :: x2, y2
    real, allocatable, dimension(:,:,:), save :: x3, y3
    real, allocatable, dimension(:,:,:), save :: x4, y4
    real, allocatable, dimension(:,:,:), save :: x5, y5
    real, allocatable, dimension(:,:,:), save :: x6, y6
    integer, allocatable :: pelist(:)
    integer :: pemax
    
    pemax = npes/2              !the partial pelist will run from 0...pemax
    !--- nullify domain list otherwise it retains memory between calls.
    call mpp_nullify_domain_list(domainx)
    call mpp_nullify_domain_list(domainy)

    allocate( gcheck(nx,ny,nz) )
    
!select pelists
    select case(type)
    case( 'Complete pelist' )
!both pelists run from 0...npes-1
        allocate( pelist(0:npes-1) )
        pelist = (/ (i,i=0,npes-1) /)
        call mpp_declare_pelist( pelist )
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)
        call mpp_declare_pelist( pelist )
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1
        if( pemax+1.GE.npes )return
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)

        call mpp_declare_pelist( pelist )
        ! z1l: the follwing will cause deadlock will happen
        ! for npes = 6, x- mpp_global_field will call mpp_sync
        call mpp_declare_pelist( (/ (i,i=pemax+1,npes-1) /))
    case default
        call mpp_error( FATAL, 'TEST_REDISTRIBUTE: no such test: '//type )
    end select
        
!set up x and y arrays
    select case(type)
    case( 'Complete pelist' )
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        allocate( x2(isd:ied,jsd:jed,nz) )
        allocate( x3(isd:ied,jsd:jed,nz) )
        allocate( x4(isd:ied,jsd:jed,nz) )
        allocate( x5(isd:ied,jsd:jed,nz) )
        allocate( x6(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
        x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
!set up y array
        call mpp_define_domains( (/1,nx,1,ny/), (/npes,1/), domainy, name=type )
        call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
        allocate( y(isd:ied,jsd:jed,nz) )
        allocate( y2(isd:ied,jsd:jed,nz) )
        allocate( y3(isd:ied,jsd:jed,nz) )
        allocate( y4(isd:ied,jsd:jed,nz) )
        allocate( y5(isd:ied,jsd:jed,nz) )
        allocate( y6(isd:ied,jsd:jed,nz) )
        y = 0.
        y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        allocate( x2(isd:ied,jsd:jed,nz) )
        allocate( x3(isd:ied,jsd:jed,nz) )
        allocate( x4(isd:ied,jsd:jed,nz) )
        allocate( x5(isd:ied,jsd:jed,nz) )
        allocate( x6(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
        x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            allocate( y2(isd:ied,jsd:jed,nz) )
            allocate( y3(isd:ied,jsd:jed,nz) )
            allocate( y4(isd:ied,jsd:jed,nz) )
            allocate( y5(isd:ied,jsd:jed,nz) )
            allocate( y6(isd:ied,jsd:jed,nz) )
            y = 0.
            y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
        end if
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1
    
!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            allocate( y2(isd:ied,jsd:jed,nz) )
            allocate( y3(isd:ied,jsd:jed,nz) )
            allocate( y4(isd:ied,jsd:jed,nz) )
            allocate( y5(isd:ied,jsd:jed,nz) )
            allocate( y6(isd:ied,jsd:jed,nz) )
            y = 0.
            y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
        else
!set up x array
            call mpp_set_current_pelist( (/ (i,i=pemax+1,npes-1) /) )
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
            call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
            allocate( x(isd:ied,jsd:jed,nz) )
            allocate( x2(isd:ied,jsd:jed,nz) )
            allocate( x3(isd:ied,jsd:jed,nz) )
            allocate( x4(isd:ied,jsd:jed,nz) )
            allocate( x5(isd:ied,jsd:jed,nz) )
            allocate( x6(isd:ied,jsd:jed,nz) )
            x = 0.
            x(is:ie,js:je,:) = global(is:ie,js:je,:)
            x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
         end if
    end select
         
!go global and redistribute
    call mpp_set_current_pelist()
    call mpp_broadcast_domain(domainx)
    call mpp_broadcast_domain(domainy)
    
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_redistribute( domainx, x, domainy, y )
    call mpp_clock_end  (id)
    
!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if
        
    call mpp_set_current_pelist()

    call mpp_clock_begin(id)
    if(ALLOCATED(y))y=0.
    call mpp_redistribute( domainx, x, domainy, y, complete=.false. )
    call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
    call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
    call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
    call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
    call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
    call mpp_clock_end  (id)
    
!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y2, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y3, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y4, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y5, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y6, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if

    call mpp_set_current_pelist()

    if(type == 'Complete pelist')then
      write(stdout(),*) 'Use domain communicator handle'
      call mpp_clock_begin(id)
      if(ALLOCATED(y))then
         y=0.; y2=0.; y3=0.; y4=0.; y5=0.; y6=0.
      endif
      call mpp_redistribute( domainx, x, domainy, y, complete=.false. )
      call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
      call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
      call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
      call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
      call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
      call mpp_clock_end  (id)
    
!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y2, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y3, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y4, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y5, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y6, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if
    endif
    dch =>NULL()
        
    call mpp_set_current_pelist()

    call mpp_sync()
    
    deallocate(gcheck)

    if(ALLOCATED(x))then
      call mpp_redistribute( domainx, x, domainy, y, free=.true.,list_size=6 )
      deallocate(x,x2,x3,x4,x5,x6)
    endif
    if(ALLOCATED(y))deallocate(y,y2,y3,y4,y5,y6)
  end subroutine test_redistribute


  subroutine test_mosaic( type )
    character(len=*), intent(in) :: type

    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile, tile
    integer        :: i, j, k, l, n, shift, xtile, ytile, ctile

    type(DomainCommunicator2D), pointer,save :: dch =>NULL()
    integer, allocatable, dimension(:)       :: npes_tile, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:,:)   :: global1, global2, gcheck
    real,    allocatable, dimension(:,:,:)   :: x, y, x1, x2, x3, x4, y1, y2, y3, y4
    real,    allocatable, dimension(:,:,:,:) :: global1_all, global2_all, global_all
    character(len=64) :: type2

    !--- check the type
    select case(type)
    case ( 'Single tile' )   !--- single with cyclic along x- and y-direction
       ntiles = 1
       num_contact = 2
       if( mod(npes,ntiles) .NE. 0 ) then 
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Single tile mosaic, npes should be multiple of ntiles(4). ' // &
                              'No test is done for Cyclic mosaic. ' )
          return
       end if
    case ( 'Cyclic mosaic' ) !--- cyclic along both x- and y-direction. 
       ntiles = 4
       num_contact = 8
       if( mod(npes,ntiles) .NE. 0 ) then 
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cyclic mosaic, npes should be multiple of ntiles(4). ' // &
                              'No test is done for Cyclic mosaic. ' )
          return
       end if
    case ( 'Cubic-grid' )
       ntiles = 6
       num_contact = 12
       !--- cubic grid always have six tiles, so npes should be multiple of 6
       if( mod(npes,ntiles) .NE. 0 .OR. nx .NE. ny) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cubic_grid mosaic, npes should be multiple of ntiles(6) ' // &
                              'and nx should equal ny, No test is done for Cubic-grid mosaic. ' )
          return
       end if
    case default
       call mpp_error(FATAL, 'TEST_MPP_DOMAINS: no such test: '//type)
    end select
       
    npes_per_tile = npes/ntiles
    call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
    allocate(layout2D(2,ntiles), global_indices(4,ntiles), npes_tile(ntiles) )
    npes_tile = npes_per_tile
    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do
    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) ) 
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) ) 

    !--- define domain
    select case(type)
    case( 'Single tile' )
       !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
       istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              npes_tile, xhalo = halo, yhalo = halo, name = type  )
    case( 'Cyclic mosaic' )
       !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
       tile1(1) = 1; tile2(1) = 2
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
       tile1(2) = 1; tile2(2) = 3
       istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
       istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
       !--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
       tile1(3) = 1; tile2(3) = 2
       istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
       istart2(3) = nx; iend2(3) = nx; jstart2(3) = 1;  jend2(3) = ny
       !--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH) 
       tile1(4) = 1; tile2(4) = 3
       istart1(4) = 1;  iend1(4) = nx; jstart1(4) = ny;  jend1(4) = ny
       istart2(4) = 1;  iend2(4) = nx; jstart2(4) = 1;   jend2(4) = 1
       !--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
       tile1(5) = 2; tile2(5) = 4
       istart1(5) = 1;  iend1(5) = nx; jstart1(5) = 1;  jend1(5) = 1
       istart2(5) = 1;  iend2(5) = nx; jstart2(5) = ny; jend2(5) = ny
       !--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
       tile1(6) = 2; tile2(6) = 4
       istart1(6) = 1;  iend1(6) = nx; jstart1(6) = ny;  jend1(6) = ny
       istart2(6) = 1;  iend2(6) = nx; jstart2(6) = 1;   jend2(6) = 1
       !--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST) 
       tile1(7) = 3; tile2(7) = 4
       istart1(7) = nx; iend1(7) = nx; jstart1(7) = 1;  jend1(7) = ny
       istart2(7) = 1;  iend2(7) = 1;  jstart2(7) = 1;  jend2(7) = ny
       !--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
       tile1(8) = 3; tile2(8) = 4
       istart1(8) = 1;  iend1(8) = 1;  jstart1(8) = 1;  jend1(8) = ny
       istart2(8) = nx; iend2(8) = nx; jstart2(8) = 1;  jend2(8) = ny
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              npes_tile, xhalo = halo, yhalo = halo, name = type  )
    case( 'Cubic-grid' )
       !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
       tile1(1) = 1; tile2(1) = 2
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
       tile1(2) = 1; tile2(2) = 3
       istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
       istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
       !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
       tile1(3) = 1; tile2(3) = 5
       istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
       istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
       !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
       tile1(4) = 1; tile2(4) = 6
       istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
       istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny       
       !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
       tile1(5) = 2; tile2(5) = 3
       istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
       istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
       !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
       tile1(6) = 2; tile2(6) = 4
       istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
       istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
       !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
       tile1(7) = 2; tile2(7) = 6
       istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
       istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
       !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
       tile1(8) = 3; tile2(8) = 4
       istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
       istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
       !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
       tile1(9) = 3; tile2(9) = 5
       istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
       istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
       !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
       tile1(10) = 4; tile2(10) = 5
       istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny; jend1(10) = ny
       istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;  jend2(10) = 1
       !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
       tile1(11) = 4; tile2(11) = 6
       istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
       istart2(11) = nx; iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = 1
       !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
       tile1(12) = 5; tile2(12) = 6
       istart1(12) = nx; iend1(12) = nx; jstart1(12) = 1;  jend1(12) = ny
       istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;  jend2(12) = ny

       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              npes_tile, topology_type = "cubic-grid", xhalo = halo, yhalo = halo, name = type )
    end select

    !--- find the tile number
    tile = mpp_pe()/npes_per_tile+1

    !--- setup data
    allocate(global2(1-halo:nx+halo,1-halo:ny+halo,nz) ) 
    allocate(global_all(1:nx,1:ny,nz, ntiles) )    
    global2 = 0
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                global_all(i,j,k,l) = l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    global2(1:nx,1:ny,:) = global_all(:,:,:,tile)

    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( gcheck(nx, ny, nz) )
    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( x1(isd:ied,jsd:jed,nz) )
    allocate( x2(isd:ied,jsd:jed,nz) )
    allocate( x3(isd:ied,jsd:jed,nz) )
    allocate( x4(isd:ied,jsd:jed,nz) )
    x = 0.
    x(is:ie,js:je,:) = global2(is:ie,js:je,:)

    x1 = x; x2 = x; x3 = x; x4 = x;

    !--- test mpp_global_field 

    gcheck = 0.    
    id = mpp_clock_id( type//' global field ', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global2(1:nx,1:ny,:), gcheck, type//' mpp_global_field ' )

    gcheck=0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, new=.true. )
    call mpp_clock_end  (id)                 
    !compare checksums between global and x arrays
    call compare_checksums( global2(1:nx,1:ny,:), gcheck, type//' mpp_global_field_new ' )

!!$    !xupdate
!!$    gcheck = 0.
!!$    call mpp_clock_begin(id)
!!$    call mpp_global_field( domain, x, gcheck, flags = XUPDATE )
!!$    call mpp_clock_end  (id)
!!$    !compare checksums between global and x arrays
!!$    call compare_checksums( global2(1:nx,js:je,:), gcheck(1:nx,js:je,:), type//' mpp_global_field xupdate only')
!!$
!!$    !yupdate
!!$    gcheck = 0.
!!$    call mpp_clock_begin(id)
!!$    call mpp_global_field( domain, x, gcheck, flags = YUPDATE )
!!$    call mpp_clock_end  (id)
!!$    !compare checksums between global and x arrays
!!$    call compare_checksums( global2(is:ie,1:ny,:), gcheck(is:ie,1:ny,:), type//' mpp_global_field yupdate only' )

!full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain )
    call mpp_clock_end  (id)

!partial update
    id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x1, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x2, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x3, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x4, domain, NUPDATE+EUPDATE, complete=.true. )
    call mpp_clock_end  (id)

    !--- fill up the value at halo points.
    select case ( type )
    case ( 'Single tile')
       global2(nx+1:nx+halo, 1:ny, :)       = global_all(1:halo, 1:ny, :, 1)            ! east
       global2(1:nx, 1-halo:0, :)           = global_all(1:nx, ny-halo+1:ny, :, 1)      ! south 
       global2(1-halo:0, 1:ny, :)           = global_all(nx-halo+1:nx, 1:ny, :, 1)      ! west
       global2(1:nx, ny+1:ny+halo, :)       = global_all(1:nx, 1:halo, :, 1)            ! north  
       global2(nx+1:nx+halo,1-halo:0,:)     = global_all(1:halo,ny-halo+1:ny,:,1)       ! southeast
       global2(1-halo:0,1-halo:0,:)         = global_all(nx-halo+1:nx,ny-halo+1:ny,:,1) ! southwest
       global2(nx+1:nx+halo,ny+1:ny+halo,:) = global_all(1:halo,1:halo,:,1)             ! northeast
       global2(1-halo:0,ny+1:ny+halo,:)     = global_all(nx-halo+1:nx,1:halo,:,1)       ! northwest       
    case ( 'Cyclic mosaic')
       select case ( tile )
       case (1)
          xtile = 2; ytile = 3; ctile = 4
       case (2)
          xtile = 1; ytile = 4; ctile = 3
       case (3)
          xtile = 4; ytile = 1; ctile = 2
       case (4)
          xtile = 3; ytile = 2; ctile = 1
       end select
       global2(nx+1:nx+halo, 1:ny, :)       = global_all(1:halo, 1:ny, :, xtile)            ! east
       global2(1:nx, 1-halo:0, :)           = global_all(1:nx, ny-halo+1:ny, :, ytile)      ! south 
       global2(1-halo:0, 1:ny, :)           = global_all(nx-halo+1:nx, 1:ny, :, xtile)      ! west
       global2(1:nx, ny+1:ny+halo, :)       = global_all(1:nx, 1:halo, :, ytile)            ! north  
       global2(nx+1:nx+halo,1-halo:0,:)     = global_all(1:halo,ny-halo+1:ny,:,ctile)       ! southeast
       global2(1-halo:0,1-halo:0,:)         = global_all(nx-halo+1:nx,ny-halo+1:ny,:,ctile) ! southwest
       global2(nx+1:nx+halo,ny+1:ny+halo,:) = global_all(1:halo,1:halo,:,ctile)             ! northeast
       global2(1-halo:0,ny+1:ny+halo,:)     = global_all(nx-halo+1:nx,1:halo,:,ctile)       ! northwest
    case ( 'Cubic-grid' )
       select case ( tile )
       case (1)
          global2(nx+1:nx+halo, 1:ny, :) = global_all(1:halo, 1:ny, :, 2) ! east 
          global2(1:nx, 1-halo:0, :)     = global_all(1:nx, ny-halo+1:ny, :, 6) ! south 
          do i = 1, halo 
             global2(1:nx, ny+i, :)    = global_all(i, ny:1:-1, :, 3) ! north 
             global2(1-i, 1:ny, :)     = global_all(nx:1:-1, ny-i+1, :, 5) ! west 
          end do
       case (2)
          global2(1:nx, ny+1:ny+halo, :) = global_all(1:nx, 1:halo, :, 3) ! north 
          global2(1-halo:0, 1:ny, :)     = global_all(nx-halo+1:nx, 1:ny, :, 1) ! west 
          do i = 1, halo 
             global2(nx+i, 1:ny, :)    = global_all(nx:1:-1, i, :, 4) ! east 
             global2(1:nx, 1-i, :)     = global_all(nx-i+1, ny:1:-1, :, 6) ! south 
          end do
       case (3)
          global2(nx+1:nx+halo, 1:ny, :) = global_all(1:halo, 1:ny, :, 4) ! east 
          global2(1:nx, 1-halo:0, :)     = global_all(1:nx, ny-halo+1:ny, :, 2) ! south 
          do i = 1, halo 
             global2(1:nx, ny+i, :)    = global_all(i, ny:1:-1, :, 5) ! north 
             global2(1-i, 1:ny, :)     = global_all(nx:1:-1, ny-i+1, :, 1) ! west 
          end do
       case (4)
          global2(1:nx, ny+1:ny+halo, :) = global_all(1:nx, 1:halo, :, 5) ! north 
          global2(1-halo:0, 1:ny, :)     = global_all(nx-halo+1:nx, 1:ny, :, 3) ! west 
          do i = 1, halo 
             global2(nx+i, 1:ny, :)    = global_all(nx:1:-1, i, :, 6) ! east 
             global2(1:nx, 1-i, :)     = global_all(nx-i+1, ny:1:-1, :, 2) ! south 
          end do
       case (5)
          global2(nx+1:nx+halo, 1:ny, :) = global_all(1:halo, 1:ny, :, 6) ! east 
          global2(1:nx, 1-halo:0, :)     = global_all(1:nx, ny-halo+1:ny, :, 4) ! south 
          do i = 1, halo 
             global2(1:nx, ny+i, :)    = global_all(i, ny:1:-1, :, 1) ! north 
             global2(1-i, 1:ny, :)     = global_all(nx:1:-1, ny-i+1, :, 3) ! west 
          end do
       case (6)
          global2(1:nx, ny+1:ny+halo, :) = global_all(1:nx, 1:halo, :, 1) ! north 
          global2(1-halo:0, 1:ny, :)     = global_all(nx-halo+1:nx, 1:ny, :, 5) ! west 
          do i = 1, halo 
             global2(nx+i, 1:ny, :)    = global_all(nx:1:-1, i, :, 2) ! east 
             global2(1:nx, 1-i, :)     = global_all(nx-i+1, ny:1:-1, :, 4) ! south 
          end do
       end select
    end select 

    call compare_checksums( x, global2(isd:ied,jsd:jed,:), type )
    call compare_checksums( x1(is:ied,js:jed,:), global2(is:ied,js:jed,:), type//' partial x1' )
    call compare_checksums( x2(is:ied,js:jed,:), global2(is:ied,js:jed,:), type//' partial x2' )
    call compare_checksums( x3(is:ied,js:jed,:), global2(is:ied,js:jed,:), type//' partial x3' )
    call compare_checksums( x4(is:ied,js:jed,:), global2(is:ied,js:jed,:), type//' partial x4' )

    deallocate(global2, global_all, x, x1, x2, x3, x4 )
    !------------------------------------------------------------------
    !              vector update : BGRID_NE, one extra point in each direction for cubic-grid
    !------------------------------------------------------------------
    !--- setup data
    shift = 0
    select case ( type ) 
    case ( 'Cyclic mosaic', 'Single tile' )
       shift = 0
    case ( 'Cubic-grid' ) 
       shift = 1
    end select

    allocate(global1(1-halo:nx+shift+halo,1-halo:ny+shift+halo,nz) ) 
    allocate(global2(1-halo:nx+shift+halo,1-halo:ny+shift+halo,nz) ) 
    allocate(global1_all(nx+shift,ny+shift,nz, ntiles),  global2_all(nx+shift,ny+shift,nz, ntiles))    
    global1 = 0; global2 = 0
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny+shift
             do i = 1, nx+shift
                global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    global1(1:nx+shift,1:ny+shift,:) = global1_all(:,:,:,tile)
    global2(1:nx+shift,1:ny+shift,:) = global2_all(:,:,:,tile)

    allocate( x (isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y (isd:ied+shift,jsd:jed+shift,nz) )
    allocate( x1(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( x2(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( x3(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( x4(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y1(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y2(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y3(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y4(isd:ied+shift,jsd:jed+shift,nz) )

    x = 0.; y = 0
    x (is:ie+shift,js:je+shift,:) = global1(is:ie+shift,js:je+shift,:)
    y (is:ie+shift,js:je+shift,:) = global2(is:ie+shift,js:je+shift,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    select case ( type ) 
    case ( 'Cyclic mosaic', 'Single tile' )
       id = mpp_clock_id( type//' vector BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
       call mpp_clock_begin(id)
       call mpp_update_domains( x,  y,  domain, gridtype=BGRID_NE)
       call mpp_update_domains( x1, y1, domain, gridtype=BGRID_NE, complete=.false. )
       call mpp_update_domains( x2, y2, domain, gridtype=BGRID_NE, complete=.true., dc_handle=dch )
       call mpp_update_domains( x3, y3, domain, gridtype=BGRID_NE, complete=.false. )
       call mpp_update_domains( x4, y4, domain, gridtype=BGRID_NE, complete=.true., dc_handle=dch )
       call mpp_clock_end  (id)
    case ( 'Cubic-grid' )
       id = mpp_clock_id( type//' paired-scalar BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
       call mpp_clock_begin(id)
       call mpp_update_domains( x,  y,  domain, flags=SCALAR_PAIR, gridtype=BGRID_NE)
       call mpp_update_domains( x1, y1, domain, flags=SCALAR_PAIR, gridtype=BGRID_NE, complete=.false. )
       call mpp_update_domains( x2, y2, domain, flags=SCALAR_PAIR, gridtype=BGRID_NE, complete=.true., dc_handle=dch )
       call mpp_update_domains( x3, y3, domain, flags=SCALAR_PAIR, gridtype=BGRID_NE, complete=.false. )
       call mpp_update_domains( x4, y4, domain, flags=SCALAR_PAIR, gridtype=BGRID_NE, complete=.true., dc_handle=dch )
       call mpp_clock_end  (id)
    end select
    dch => NULL()

    !-----------------------------------------------------------------------
    !                   fill up the value at halo points. For the cubic-grid, 
    !   On the contact line, the following relation will be used to 
    !   --- fill the value on contact line ( balance send and recv).
    !       2W --> 1E, 1S --> 6N, 3W --> 1N, 4S --> 2E
    !       4W --> 3E, 3S --> 2N, 1W --> 5N, 2S --> 6E
    !       6W --> 5E, 5S --> 4N, 5W --> 3N, 6S --> 4E
    !   --- the 8 corner points will be decided by the following.
    !   --- 1,2,3: take value at 3, need to send southwest point from 3 to northeast point at 1.
    !   --- 2,3,4: take value at 4, need to send southwest point from 4 to northeast point at 2.
    !   --- 3,4,5: take value at 5, need to send southwest point from 5 to northeast point at 3.
    !   --- 4,5,6: take value at 6, need to send southwest point from 6 to northeast point at 4.
    !   --- 5,6,1: take value at 1, need to send southwest point from 1 to northeast point at 5.
    !   --- 6,1,2: take value at 2, need to send southwest point from 2 to northeast point at 6.
    !   --- 1,3,5: take value at 3
    !   --- 2,4,6: take value at 4
    select case ( type )
    case ( 'Single tile' )
       global1(nx+1:nx+halo, 1:ny, :)       = global1_all(1:halo, 1:ny, :, 1)            ! east
       global1(1:nx, 1-halo:0, :)           = global1_all(1:nx, ny-halo+1:ny, :, 1)      ! south 
       global1(1-halo:0, 1:ny, :)           = global1_all(nx-halo+1:nx, 1:ny, :, 1)      ! west
       global1(1:nx, ny+1:ny+halo, :)       = global1_all(1:nx, 1:halo, :, 1)            ! north  
       global1(nx+1:nx+halo,1-halo:0,:)     = global1_all(1:halo,ny-halo+1:ny,:,1)       ! southeast
       global1(1-halo:0,1-halo:0,:)         = global1_all(nx-halo+1:nx,ny-halo+1:ny,:,1) ! southwest
       global1(nx+1:nx+halo,ny+1:ny+halo,:) = global1_all(1:halo,1:halo,:,1)             ! northeast
       global1(1-halo:0,ny+1:ny+halo,:)     = global1_all(nx-halo+1:nx,1:halo,:,1)       ! northwest
       global2(nx+1:nx+halo, 1:ny, :)       = global2_all(1:halo, 1:ny, :, 1)            ! east
       global2(1:nx, 1-halo:0, :)           = global2_all(1:nx, ny-halo+1:ny, :, 1)      ! south 
       global2(1-halo:0, 1:ny, :)           = global2_all(nx-halo+1:nx, 1:ny, :, 1)      ! west
       global2(1:nx, ny+1:ny+halo, :)       = global2_all(1:nx, 1:halo, :, 1)            ! north 
       global2(nx+1:nx+halo,1-halo:0,:)     = global2_all(1:halo,ny-halo+1:ny,:,1)       ! southeast
       global2(1-halo:0,1-halo:0,:)         = global2_all(nx-halo+1:nx,ny-halo+1:ny,:,1) ! southwest
       global2(nx+1:nx+halo,ny+1:ny+halo,:) = global2_all(1:halo,1:halo,:,1)             ! northeast
       global2(1-halo:0,ny+1:ny+halo,:)     = global2_all(nx-halo+1:nx,1:halo,:,1)       ! northwest      
    case ( 'Cyclic mosaic' )
       select case ( tile )
       case (1)
          xtile = 2; ytile = 3; ctile = 4
       case (2)
          xtile = 1; ytile = 4; ctile = 3
       case (3)
          xtile = 4; ytile = 1; ctile = 2
       case (4)
          xtile = 3; ytile = 2; ctile = 1
       end select
       global1(nx+1:nx+halo, 1:ny, :)       = global1_all(1:halo, 1:ny, :, xtile)            ! east
       global1(1:nx, 1-halo:0, :)           = global1_all(1:nx, ny-halo+1:ny, :, ytile)      ! south 
       global1(1-halo:0, 1:ny, :)           = global1_all(nx-halo+1:nx, 1:ny, :, xtile)      ! west
       global1(1:nx, ny+1:ny+halo, :)       = global1_all(1:nx, 1:halo, :, ytile)            ! north  
       global1(nx+1:nx+halo,1-halo:0,:)     = global1_all(1:halo,ny-halo+1:ny,:,ctile)       ! southeast
       global1(1-halo:0,1-halo:0,:)         = global1_all(nx-halo+1:nx,ny-halo+1:ny,:,ctile) ! southwest
       global1(nx+1:nx+halo,ny+1:ny+halo,:) = global1_all(1:halo,1:halo,:,ctile)             ! northeast
       global1(1-halo:0,ny+1:ny+halo,:)     = global1_all(nx-halo+1:nx,1:halo,:,ctile)       ! northwest
       global2(nx+1:nx+halo, 1:ny, :)       = global2_all(1:halo, 1:ny, :, xtile)            ! east
       global2(1:nx, 1-halo:0, :)           = global2_all(1:nx, ny-halo+1:ny, :, ytile)      ! south 
       global2(1-halo:0, 1:ny, :)           = global2_all(nx-halo+1:nx, 1:ny, :, xtile)      ! west
       global2(1:nx, ny+1:ny+halo, :)       = global2_all(1:nx, 1:halo, :, ytile)            ! north 
       global2(nx+1:nx+halo,1-halo:0,:)     = global2_all(1:halo,ny-halo+1:ny,:,ctile)       ! southeast
       global2(1-halo:0,1-halo:0,:)         = global2_all(nx-halo+1:nx,ny-halo+1:ny,:,ctile) ! southwest
       global2(nx+1:nx+halo,ny+1:ny+halo,:) = global2_all(1:halo,1:halo,:,ctile)             ! northeast
       global2(1-halo:0,ny+1:ny+halo,:)     = global2_all(nx-halo+1:nx,1:halo,:,ctile)       ! northwest          
    case ( 'Cubic-grid' )
       select case ( tile )
       case (1)
          global1(nx+1:nx+1+halo, 1:ny+1, :) = global1_all(1:halo+1, 1:ny+1, :, 2) ! east 
          global1(1:nx+1, 1-halo:0, :)       = global1_all(1:nx+1, ny-halo+1:ny, :, 6) ! south 
          do i = 1, halo 
             global1(1-i, 1:ny+1, :)     = global1_all(nx+1:1:-1, ny-i+1, :, 5) ! west 
          end do
          do i = 1, halo+1 
             global1(1:nx+1, ny+i, :)    = global1_all(i, ny+1:1:-1, :, 3) ! north 
          end do
          global2(nx+1:nx+1+halo, 1:ny+1, :) = global2_all(1:halo+1, 1:ny+1, :, 2) ! east 
          global2(1:nx+1, 1-halo:0, :)       = global2_all(1:nx+1, ny-halo+1:ny, :, 6) ! south 
          do i = 1, halo 
             global2(1-i, 1:ny+1, :)     = global2_all(nx+1:1:-1, ny-i+1, :, 5) ! west 
          end do
          do i = 1, halo+1 
             global2(1:nx+1, ny+i, :)    = global2_all(i, ny+1:1:-1, :, 3) ! north 
          end do

          global1(nx+1,ny+1,:)   = global1_all(1,1,:,3)
          global2(nx+1,ny+1,:)   = global2_all(1,1,:,3)
          global1(1,   ny+1,:)   = global1_all(1,ny+1,:,3)
          global2(1,   ny+1,:)   = global2_all(1,ny+1,:,3)
       case (2)  
          global1(1:nx+1, ny+1:ny+1+halo, :) = global1_all(1:nx+1, 1:halo+1, :, 3) ! north 
          global1(1-halo:0, 1:ny+1, :)       = global1_all(nx-halo+1:nx, 1:ny+1, :, 1) ! west
          do i = 1, halo 
             global1(1:nx+1, 1-i, :)     = global1_all(nx-i+1, ny+1:1:-1, :, 6) ! south 
          end do
          do i = 1, halo+1 
             global1(nx+i, 1:ny+1, :)    = global1_all(nx+1:1:-1, i, :, 4) ! east 
          end do
          global2(1:nx+1, ny+1:ny+1+halo, :) = global2_all(1:nx+1, 1:halo+1, :, 3) ! north 
          global2(1-halo:0, 1:ny+1, :)       = global2_all(nx-halo+1:nx, 1:ny+1, :, 1) ! west
          do i = 1, halo 
             global2(1:nx+1, 1-i, :)     = global2_all(nx-i+1, ny+1:1:-1, :, 6) ! south 
          end do
          do i = 1, halo+1 
             global2(nx+i, 1:ny+1, :)    = global2_all(nx+1:1:-1, i, :, 4) ! east 
          end do
          global1(nx+1,ny+1,:)   = global1_all(1,1,:,4)
          global2(nx+1,ny+1,:)   = global2_all(1,1,:,4)
          global1(nx+1,   1,:)   = global1_all(nx+1,1,:,4)
          global2(nx+1,   1,:)   = global2_all(nx+1,1,:,4)
       case (3)
          global1(nx+1:nx+1+halo, 1:ny+1, :) = global1_all(1:halo+1, 1:ny+1, :, 4) ! east 
          global1(1:nx+1, 1-halo:0, :)       = global1_all(1:nx+1, ny-halo+1:ny, :, 2) ! south 
          do i = 1, halo 
             global1(1-i, 1:ny+1, :)     = global1_all(nx+1:1:-1, ny-i+1, :, 1) ! west 
          end do
          do i = 1, halo+1
             global1(1:nx+1, ny+i, :)    = global1_all(i, ny+1:1:-1, :, 5) ! north 
          end do
          global2(nx+1:nx+1+halo, 1:ny+1, :) = global2_all(1:halo+1, 1:ny+1, :, 4) ! east 
          global2(1:nx+1, 1-halo:0, :)       = global2_all(1:nx+1, ny-halo+1:ny, :, 2) ! south 
          do i = 1, halo 
             global2(1-i, 1:ny+1, :)     = global2_all(nx+1:1:-1, ny-i+1, :, 1) ! west 
          end do
          do i = 1, halo+1
             global2(1:nx+1, ny+i, :)    = global2_all(i, ny+1:1:-1, :, 5) ! north 
          end do

          global1(nx+1,ny+1,:)   = global1_all(1,1,:,5)
          global2(nx+1,ny+1,:)   = global2_all(1,1,:,5)
          global1(1,   ny+1,:)   = global1_all(1,ny+1,:,3)
          global2(1,   ny+1,:)   = global2_all(1,ny+1,:,3)
       case (4)
          global1(1:nx+1, ny+1:ny+1+halo, :) = global1_all(1:nx+1, 1:halo+1, :, 5) ! north 
          global1(1-halo:0, 1:ny+1, :)       = global1_all(nx-halo+1:nx, 1:ny+1, :, 3) ! west 
          do i = 1, halo 
             global1(1:nx+1, 1-i, :)     = global1_all(nx-i+1, ny+1:1:-1, :, 2) ! south 
          end do
          do i = 1, halo+1
             global1(nx+i, 1:ny+1, :)    = global1_all(nx+1:1:-1, i, :, 6) ! east 
          end do
          global2(1:nx+1, ny+1:ny+1+halo, :) = global2_all(1:nx+1, 1:halo+1, :, 5) ! north 
          global2(1-halo:0, 1:ny+1, :)       = global2_all(nx-halo+1:nx, 1:ny+1, :, 3) ! west 
          do i = 1, halo 
             global2(1:nx+1, 1-i, :)     = global2_all(nx-i+1, ny+1:1:-1, :, 2) ! south 
          end do
          do i = 1, halo+1
             global2(nx+i, 1:ny+1, :)    = global2_all(nx+1:1:-1, i, :, 6) ! east 
          end do
          global1(nx+1,ny+1,:)   = global1_all(1,1,:,6)
          global2(nx+1,ny+1,:)   = global2_all(1,1,:,6)
          global1(nx+1,   1,:)   = global1_all(nx+1,1,:,4)
          global2(nx+1,   1,:)   = global2_all(nx+1,1,:,4)
       case (5)
          global1(nx+1:nx+1+halo, 1:ny+1, :) = global1_all(1:halo+1, 1:ny+1, :, 6) ! east 
          global1(1:nx+1, 1-halo:0, :)       = global1_all(1:nx+1, ny-halo+1:ny, :, 4) ! south 
          do i = 1, halo 
             global1(1-i, 1:ny+1, :)     = global1_all(nx+1:1:-1, ny-i+1, :, 3) ! west 
          end do
          do i = 1, halo+1
             global1(1:nx+1, ny+i, :)    = global1_all(i, ny+1:1:-1, :, 1) ! north 
          end do
          global2(nx+1:nx+1+halo, 1:ny+1, :) = global2_all(1:halo+1, 1:ny+1, :, 6) ! east 
          global2(1:nx+1, 1-halo:0, :)       = global2_all(1:nx+1, ny-halo+1:ny, :, 4) ! south 
          do i = 1, halo 
             global2(1-i, 1:ny+1, :)     = global2_all(nx+1:1:-1, ny-i+1, :, 3) ! west 
          end do
          do i = 1, halo+1
             global2(1:nx+1, ny+i, :)    = global2_all(i, ny+1:1:-1, :, 1) ! north 
          end do
          global1(nx+1,ny+1,:)   = global1_all(1,1,:,1)
          global2(nx+1,ny+1,:)   = global2_all(1,1,:,1)
          global1(1,   ny+1,:)   = global1_all(1,ny+1,:,3)
          global2(1,   ny+1,:)   = global2_all(1,ny+1,:,3)
       case (6)
          global1(1:nx+1, ny+1:ny+1+halo, :) = global1_all(1:nx+1, 1:halo+1, :, 1) ! north 
          global1(1-halo:0, 1:ny+1, :)       = global1_all(nx-halo+1:nx, 1:ny+1, :, 5) ! west 
          do i = 1, halo 
             global1(1:nx+1, 1-i, :)     = global1_all(nx-i+1, ny+1:1:-1, :, 4) ! south 
          end do
          do i = 1, halo+1 
             global1(nx+i, 1:ny+1, :)    = global1_all(nx+1:1:-1, i, :, 2) ! east 
          end do
          global2(1:nx+1, ny+1:ny+1+halo, :) = global2_all(1:nx+1, 1:halo+1, :, 1) ! north 
          global2(1-halo:0, 1:ny+1, :)       = global2_all(nx-halo+1:nx, 1:ny+1, :, 5) ! west 
          do i = 1, halo 
             global2(1:nx+1, 1-i, :)     = global2_all(nx-i+1, ny+1:1:-1, :, 4) ! south 
          end do
          do i = 1, halo+1 
             global2(nx+i, 1:ny+1, :)    = global2_all(nx+1:1:-1, i, :, 2) ! east 
          end do
          global1(nx+1,ny+1,:)   = global1_all(1,1,:,2)
          global2(nx+1,ny+1,:)   = global2_all(1,1,:,2)
          global1(nx+1,   1,:)   = global1_all(nx+1,1,:,4)
          global2(nx+1,   1,:)   = global2_all(nx+1,1,:,4)
       end select
    end select

    select case ( type )
    case ( 'Cyclic mosaic', 'Single tile' )
       type2 = type//' vector BGRID_NE'
    case ( 'Cubic-grid' )
       type2 = type//' paired-scalar BGRID_NE'
    end select
    call compare_checksums( x,  global1(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' X' )
    call compare_checksums( y,  global2(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' Y' )
    call compare_checksums( x1, global1(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' X1' )
    call compare_checksums( x2, global1(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' X2' )
    call compare_checksums( x3, global1(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' X3' )
    call compare_checksums( x4, global1(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' X4' )
    call compare_checksums( y1, global2(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' Y1' )
    call compare_checksums( y2, global2(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' Y2' )
    call compare_checksums( y3, global2(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' Y3' )
    call compare_checksums( y4, global2(isd:ied+shift,jsd:jed+shift,:), trim(type2)//' Y4' )

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- setup data
    if( type == 'Cubic-grid' ) then
       deallocate(global1, global2, x, y, x1, x2, x3, x4, y1, y2, y3, y4)
       allocate(global1(1-halo:nx+shift+halo,1-halo:ny  +halo,nz) ) 
       allocate(global2(1-halo:nx  +halo,1-halo:ny+shift+halo,nz) ) 
       global1 = 0; global2 = 0

       global1(1:nx+shift,1:ny  ,:) = global1_all(1:nx+shift,1:ny,  :,tile)
       global2(1:nx  ,1:ny+shift,:) = global2_all(1:nx  ,1:ny+shift,:,tile)

       allocate( x (isd:ied+shift,jsd:jed  ,nz) )
       allocate( y (isd:ied  ,jsd:jed+shift,nz) )
       allocate( x1(isd:ied+shift,jsd:jed  ,nz) )
       allocate( x2(isd:ied+shift,jsd:jed  ,nz) )
       allocate( x3(isd:ied+shift,jsd:jed  ,nz) )
       allocate( x4(isd:ied+shift,jsd:jed  ,nz) )
       allocate( y1(isd:ied  ,jsd:jed+shift,nz) )
       allocate( y2(isd:ied  ,jsd:jed+shift,nz) )
       allocate( y3(isd:ied  ,jsd:jed+shift,nz) )
       allocate( y4(isd:ied  ,jsd:jed+shift,nz) )
    end if
    x = 0.; y = 0.
    x (is:ie+shift,js:je  ,:) = global1(is:ie+shift,js:je  ,:)
    y (is:ie  ,js:je+shift,:) = global2(is:ie  ,js:je+shift,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    id = mpp_clock_id( type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=CGRID_NE, complete=.true., dc_handle=dch )
    call mpp_update_domains( x3, y3, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=CGRID_NE, complete=.true., dc_handle=dch )
    call mpp_clock_end  (id)
    dch => NULL()
    !-----------------------------------------------------------------------
    !                   fill up the value at halo points for cubic-grid.
    !   On the contact line, the following relation will be used to 
    !   --- fill the value on contact line ( balance send and recv).
    !       2W --> 1E, 1S --> 6N, 3W --> 1N, 4S --> 2E
    !       4W --> 3E, 3S --> 2N, 1W --> 5N, 2S --> 6E
    !       6W --> 5E, 5S --> 4N, 5W --> 3N, 6S --> 4E
    !---------------------------------------------------------------------------
    if( type == 'Cubic-grid' ) then    
       select case ( tile )
       case (1)
          global1(nx+1:nx+1+halo, 1:ny, :) = global1_all(1:halo+1, 1:ny, :, 2) ! east 
          global1(1:nx+1, 1-halo:0, :)     = global1_all(1:nx+1, ny-halo+1:ny, :, 6) ! south 
          do i = 1, halo 
             global1(1-i, 1:ny, :)     = global2_all(nx:1:-1, ny-i+1, :, 5) ! west 
          end do
          do i = 1, halo 
             global1(1:nx+1, ny+i, :)    = -global2_all(i, ny+1:1:-1, :, 3) ! north 
          end do
          global2(nx+1:nx+halo, 1:ny+1, :) = global2_all(1:halo, 1:ny+1, :, 2) ! east 
          global2(1:nx, 1-halo:0, :)       = global2_all(1:nx, ny-halo+1:ny, :, 6) ! south 
          do i = 1, halo 
             global2(1-i, 1:ny+1, :)     = -global1_all(nx+1:1:-1, ny-i+1, :, 5) ! west 
          end do
          do i = 1, halo+1
             global2(1:nx, ny+i, :)    = global1_all(i, ny:1:-1, :, 3) ! north 
          end do
       case (2)  
          global1(1:nx+1, ny+1:ny+halo, :) = global1_all(1:nx+1, 1:halo, :, 3) ! north 
          global1(1-halo:0, 1:ny, :)       = global1_all(nx-halo+1:nx, 1:ny, :, 1) ! west
          do i = 1, halo 
             global1(1:nx+1, 1-i, :)     = -global2_all(nx-i+1, ny+1:1:-1, :, 6) ! south 
          end do
          do i = 1, halo+1 
             global1(nx+i, 1:ny, :)    = global2_all(nx:1:-1, i, :, 4) ! east 
          end do
          global2(1:nx, ny+1:ny+1+halo, :) = global2_all(1:nx, 1:halo+1, :, 3) ! north 
          global2(1-halo:0, 1:ny+1, :)       = global2_all(nx-halo+1:nx, 1:ny+1, :, 1) ! west
          do i = 1, halo 
             global2(1:nx, 1-i, :)     = global1_all(nx-i+1, ny:1:-1, :, 6) ! south 
          end do
          do i = 1, halo 
             global2(nx+i, 1:ny+1, :)    = -global1_all(nx+1:1:-1, i, :, 4) ! east 
          end do
       case (3)
          global1(nx+1:nx+1+halo, 1:ny, :) = global1_all(1:halo+1, 1:ny, :, 4) ! east 
          global1(1:nx+1, 1-halo:0, :)       = global1_all(1:nx+1, ny-halo+1:ny, :, 2) ! south 
          do i = 1, halo 
             global1(1-i, 1:ny, :)     = global2_all(nx:1:-1, ny-i+1, :, 1) ! west 
          end do
          do i = 1, halo
             global1(1:nx+1, ny+i, :)    = -global2_all(i, ny+1:1:-1, :, 5) ! north 
          end do
          global2(nx+1:nx+halo, 1:ny+1, :) = global2_all(1:halo, 1:ny+1, :, 4) ! east 
          global2(1:nx, 1-halo:0, :)       = global2_all(1:nx, ny-halo+1:ny, :, 2) ! south 
          do i = 1, halo 
             global2(1-i, 1:ny+1, :)     = -global1_all(nx+1:1:-1, ny-i+1, :, 1) ! west 
          end do
          do i = 1, halo+1
             global2(1:nx, ny+i, :)    = global1_all(i, ny:1:-1, :, 5) ! north 
          end do
       case (4)
          global1(1:nx+1, ny+1:ny+halo, :) = global1_all(1:nx+1, 1:halo, :, 5) ! north 
          global1(1-halo:0, 1:ny, :)       = global1_all(nx-halo+1:nx, 1:ny, :, 3) ! west 
          do i = 1, halo 
             global1(1:nx+1, 1-i, :)     = -global2_all(nx-i+1, ny+1:1:-1, :, 2) ! south 
          end do
          do i = 1, halo+1
             global1(nx+i, 1:ny, :)    = global2_all(nx:1:-1, i, :, 6) ! east 
          end do
          global2(1:nx, ny+1:ny+1+halo, :) = global2_all(1:nx, 1:halo+1, :, 5) ! north 
          global2(1-halo:0, 1:ny+1, :)       = global2_all(nx-halo+1:nx, 1:ny+1, :, 3) ! west 
          do i = 1, halo 
             global2(1:nx, 1-i, :)     = global1_all(nx-i+1, ny:1:-1, :, 2) ! south 
          end do
          do i = 1, halo
             global2(nx+i, 1:ny+1, :)    = -global1_all(nx+1:1:-1, i, :, 6) ! east 
          end do
       case (5)
          global1(nx+1:nx+1+halo, 1:ny, :) = global1_all(1:halo+1, 1:ny, :, 6) ! east 
          global1(1:nx+1, 1-halo:0, :)       = global1_all(1:nx+1, ny-halo+1:ny, :, 4) ! south 
          do i = 1, halo 
             global1(1-i, 1:ny, :)     = global2_all(nx:1:-1, ny-i+1, :, 3) ! west 
          end do
          do i = 1, halo
             global1(1:nx+1, ny+i, :)    = -global2_all(i, ny+1:1:-1, :, 1) ! north 
          end do
          global2(nx+1:nx+halo, 1:ny+1, :) = global2_all(1:halo, 1:ny+1, :, 6) ! east 
          global2(1:nx, 1-halo:0, :)       = global2_all(1:nx, ny-halo+1:ny, :, 4) ! south 
          do i = 1, halo 
             global2(1-i, 1:ny+1, :)     = -global1_all(nx+1:1:-1, ny-i+1, :, 3) ! west 
          end do
          do i = 1, halo+1
             global2(1:nx, ny+i, :)    = global1_all(i, ny:1:-1, :, 1) ! north 
          end do
       case (6)
          global1(1:nx+1, ny+1:ny+halo, :) = global1_all(1:nx+1, 1:halo, :, 1) ! north 
          global1(1-halo:0, 1:ny, :)       = global1_all(nx-halo+1:nx, 1:ny, :, 5) ! west 
          do i = 1, halo 
             global1(1:nx+1, 1-i, :)     = -global2_all(nx-i+1, ny+1:1:-1, :, 4) ! south 
          end do
          do i = 1, halo+1 
             global1(nx+i, 1:ny, :)    = global2_all(nx:1:-1, i, :, 2) ! east 
          end do
          global2(1:nx, ny+1:ny+1+halo, :) = global2_all(1:nx, 1:halo+1, :, 1) ! north 
          global2(1-halo:0, 1:ny+1, :)       = global2_all(nx-halo+1:nx, 1:ny+1, :, 5) ! west 
          do i = 1, halo 
             global2(1:nx, 1-i, :)     = global1_all(nx-i+1, ny:1:-1, :, 4) ! south 
          end do
          do i = 1, halo
             global2(nx+i, 1:ny+1, :)    = -global1_all(nx+1:1:-1, i, :, 2) ! east 
          end do
       end select
    end if

    call compare_checksums( x,  global1(isd:ied+shift,jsd:jed,  :), type//' CGRID_NE X' )
    call compare_checksums( y,  global2(isd:ied,  jsd:jed+shift,:), type//' CGRID_NE Y' )
    call compare_checksums( x1, global1(isd:ied+shift,jsd:jed,  :), type//' CGRID_NE X1' )
    call compare_checksums( x2, global1(isd:ied+shift,jsd:jed,  :), type//' CGRID_NE X2' )
    call compare_checksums( x3, global1(isd:ied+shift,jsd:jed,  :), type//' CGRID_NE X3' )
    call compare_checksums( x4, global1(isd:ied+shift,jsd:jed,  :), type//' CGRID_NE X4' )
    call compare_checksums( y1, global2(isd:ied,  jsd:jed+shift,:), type//' CGRID_NE Y1' )
    call compare_checksums( y2, global2(isd:ied,  jsd:jed+shift,:), type//' CGRID_NE Y2' )
    call compare_checksums( y3, global2(isd:ied,  jsd:jed+shift,:), type//' CGRID_NE Y3' )
    call compare_checksums( y4, global2(isd:ied,  jsd:jed+shift,:), type//' CGRID_NE Y4' )

    deallocate(global1, global2, x, y, x1, x2, x3, x4, y1, y2, y3, y4, global1_all, global2_all)
    deallocate(layout2D, global_indices, npes_tile, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2 ) 

  end subroutine test_mosaic
    
  subroutine test_halo_update( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, x1, x2, x3, x4
    real, allocatable, dimension(:,:,:) :: y, y1, y2, y3, y4
    type(domain2D) :: domain
    type(DomainCommunicator2D), pointer, save :: dch =>NULL()
    real,    allocatable :: global1(:,:,:), global2(:,:,:)
    logical, allocatable :: maskmap(:,:)
    integer :: shift, i

    ! when testing maskmap option, nx*ny should be able to be divided by both npes and npes+1
    if(type == 'Masked' .or. type == 'Masked symmetry') then
       if(mod(nx*ny, npes) .NE. 0 .OR. mod(nx*ny, npes+1) .NE. 0 ) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: nx*ny can not be divided by both npes and npes+1, '//&
               'Masked test_halo_update will not be tested')
          return
       end if
    end if

    allocate(global2(1-halo:nx+halo,1-halo:ny+halo,nz) )
    global2 = global
    
    select case(type)
    case( 'Simple', 'Simple symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        if(type == 'Simple') then
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name=type )
        else
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name=type, symmetry = .true. )
        endif
        global2(1-halo:0,    :,:) = 0
        global2(nx+1:nx+halo,:,:) = 0
        global2(:,    1-halo:0,:) = 0
        global2(:,ny+1:ny+halo,:) = 0
    case( 'Cyclic', 'Cyclic symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        if(type == 'Cyclic') then
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
                xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, name=type )
        else
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
                xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, name=type, symmetry = .true. )
        endif
        global2(1-halo:0,    1:ny,:) = global2(nx-halo+1:nx,1:ny,:)
        global2(nx+1:nx+halo,1:ny,:) = global2(1:halo,      1:ny,:)
        global2(1-halo:nx+halo,    1-halo:0,:) = global2(1-halo:nx+halo,ny-halo+1:ny,:)
        global2(1-halo:nx+halo,ny+1:ny+halo,:) = global2(1-halo:nx+halo,1:halo,      :)
    case( 'Folded', 'Folded symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        if(type == 'Folded') then
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
                xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, name=type )
        else
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
                xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, name=type, symmetry = .true. )
        endif
        global2(1-halo:0,    1:ny,:) = global2(nx-halo+1:nx,1:ny,:)
        global2(nx+1:nx+halo,1:ny,:) = global2(1:halo,      1:ny,:)
        global2(1-halo:nx+halo,ny+1:ny+halo,:) = global2(nx+halo:1-halo:-1,ny:ny-halo+1:-1,:)
        global2(1-halo:nx+halo,1-halo:0,:) = 0
    case( 'Masked', 'Masked symmetry' )
!with fold and cyclic, assign to npes+1 and mask out the top-rightdomain
        call mpp_define_layout( (/1,nx,1,ny/), npes+1, layout )
        allocate( maskmap(layout(1),layout(2)) )
        maskmap(:,:) = .TRUE.; maskmap(layout(1),layout(2)) = .FALSE.
        if(type == 'Masked') then
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
                xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, maskmap=maskmap, name=type )
        else 
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
                xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, maskmap=maskmap, name=type, symmetry = .true. )
        endif
        deallocate(maskmap)
!we need to zero out the global data on the missing domain.
!this logic assumes top-right, in an even division
        if( mod(nx,layout(1)).NE.0 .OR. mod(ny,layout(2)).NE.0 )call mpp_error( FATAL, &
             'TEST_MPP_DOMAINS: test for masked domains needs (nx,ny) to divide evenly on npes+1 PEs.' )
        global2(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny,:) = 0
!then apply same folded logic as above
        global2(1-halo:0,    1:ny,:) = global2(nx-halo+1:nx,1:ny,:)
        global2(nx+1:nx+halo,1:ny,:) = global2(1:halo,      1:ny,:)
        global2(1-halo:nx+halo,ny+1:ny+halo,:) = global2(nx+halo:1-halo:-1,ny:ny-halo+1:-1,:)
        global2(1-halo:nx+halo,1-halo:0,:) = 0
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select
        
!set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( x1(isd:ied,jsd:jed,nz) )
    allocate( x2(isd:ied,jsd:jed,nz) )
    allocate( x3(isd:ied,jsd:jed,nz) )
    allocate( x4(isd:ied,jsd:jed,nz) )
    x = 0.; x1 = 0.; x2 = 0.; x3 = 0.; x4 = 0.
    x (is:ie,js:je,:) = global2(is:ie,js:je,:)
    x1(is:ie,js:je,:) = global2(is:ie,js:je,:)
    x2(is:ie,js:je,:) = global2(is:ie,js:je,:)
    x3(is:ie,js:je,:) = global2(is:ie,js:je,:)
    x4(is:ie,js:je,:) = global2(is:ie,js:je,:)

!full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain )
    call mpp_clock_end  (id)
    call compare_checksums( x, global2(isd:ied,jsd:jed,:), type )

!partial update
    id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x1, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x2, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x3, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x4, domain, NUPDATE+EUPDATE, complete=.true. )
    call mpp_clock_end  (id)
    call compare_checksums( x1(is:ied,js:jed,:), global2(is:ied,js:jed,:), type//' partial x1' )
    call compare_checksums( x2(is:ied,js:jed,:), global2(is:ied,js:jed,:), type//' partial x2' )
    call compare_checksums( x3(is:ied,js:jed,:), global2(is:ied,js:jed,:), type//' partial x3' )
    call compare_checksums( x4(is:ied,js:jed,:), global2(is:ied,js:jed,:), type//' partial x4' )
    
    !--- test vector update for FOLDED and MASKED case.
    if(type == 'Simple' .or. type == 'Simple symmetry' .or. type == 'Cyclic' .or. type == 'Cyclic symmetry') then
       deallocate(x,x1,x2,x3,x4)
       return       
    end if

    !------------------------------------------------------------------
    !              vector update : BGRID_NE
    !------------------------------------------------------------------
    if(type == 'Masked' .OR. type =='Folded') then
       shift = 0
    else
       shift = 1
       deallocate(global2)
       allocate(global2(1-halo:nx+halo+shift,1-halo:ny+halo+shift,nz) )
       global2 = 0.0
       do k = 1,nz
          do j = 1,ny+1
             do i = 1,nx+1
                global2(i,j,k) = k + i*1e-3 + j*1e-6
             end do
          end do
       end do
       if(type == 'Masked symmetry') then
           global2(nx-nx/layout(1)+1:nx+1,ny-ny/layout(2)+1:ny+1,:) = 0
       endif
       deallocate(x, x1, x2, x3, x4)
       allocate( x (isd:ied+1,jsd:jed+1,nz) )
       allocate( x1(isd:ied+1,jsd:jed+1,nz) )
       allocate( x2(isd:ied+1,jsd:jed+1,nz) )
       allocate( x3(isd:ied+1,jsd:jed+1,nz) )
       allocate( x4(isd:ied+1,jsd:jed+1,nz) )
    endif

    select case (type)
    case ('Folded', 'Masked')
       !fill in folded north edge, cyclic east and west edge
       global2(1-halo:0,    1:ny,:) = global2(nx-halo+1:nx,1:ny,:)
       global2(nx+1:nx+halo,1:ny,:) = global2(1:halo,      1:ny,:)
       global2(1-halo:nx+halo-1,ny+1:ny+halo,:) = -global2(nx+halo-1:1-halo:-1,ny-1:ny-halo:-1,:)
       global2(nx+halo,         ny+1:ny+halo,:) = -global2(nx-halo,ny-1:ny-halo:-1,:)
       global2(1-halo:nx+halo,      1-halo:0,:) = 0
    case ('Folded symmetry', 'Masked symmetry' )
       global2(1-halo:0,    1:ny+1,:)   = global2(nx-halo+1:nx,1:ny+1,:)
       global2(nx+1:nx+halo+1,1:ny+1,:) = global2(1:halo+1,    1:ny+1,:)
       global2(1-halo:nx+halo+1,ny+2:ny+halo+1,:) = -global2(nx+halo+1:1-halo:-1,ny:ny-halo+1:-1,:)
       global2(1-halo:nx+halo+1,1-halo:0,:) = 0
    end select

    x = 0.
    x(is:ie+shift,js:je+shift,:) = global2(is:ie+shift,js:je+shift,:)
    !set up y array
    allocate( y (isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y1(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y2(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y3(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y4(isd:ied+shift,jsd:jed+shift,nz) )
    y = x; x1 = x; x2 = x; x3 = x; x4 = x
    y = x; y1 = x; y2 = x; y3 = x; y4 = x

    id = mpp_clock_id( type//' vector BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=BGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=BGRID_NE, complete=.true., dc_handle=dch )
    call mpp_update_domains( x3, y3, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=BGRID_NE, complete=.true., dc_handle=dch )
    call mpp_clock_end  (id)
    dch => NULL()

    !redundant points must be equal and opposite
    global2(nx/2+shift,ny+shift,:) = 0.  !pole points must have 0 velocity
    global2(nx+shift  ,ny+shift,:) = 0.  !pole points must have 0 velocity
    global2(nx/2+1+shift:nx-1+shift,  ny+shift,:) = -global2(nx/2-1+shift:1+shift:-1, ny+shift,:)
    global2(1-halo:shift,             ny+shift,:) = -global2(nx-halo+1:nx+shift,      ny+shift,:)
    global2(nx+1+shift:nx+halo+shift, ny+shift,:) = -global2(1+shift:halo+shift,      ny+shift,:)
    !--- the following will fix the +0/-0 problem on altix
    if(halo >0) global2(shift,ny+shift,:) = 0.  !pole points must have 0 velocity

    call compare_checksums( x,  global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X' )
    call compare_checksums( y,  global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y' )
    call compare_checksums( x1, global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X1' )
    call compare_checksums( x2, global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X2' )
    call compare_checksums( x3, global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X3' )
    call compare_checksums( x4, global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X4' )
    call compare_checksums( y1, global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y1' )
    call compare_checksums( y2, global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y2' )
    call compare_checksums( y3, global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y3' )
    call compare_checksums( y4, global2(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y4' )
    !--- delete the communicator
!!$    call mpp_update_domains( x1, domain, free=.true., list_size=4 )
    deallocate(global2, x, x1, x2, x3, x4, y, y1, y2, y3, y4)

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- global1 is x-component and global2 is y-component
    allocate(global1(1-halo:nx+halo+shift, 1-halo:ny+halo,nz), global2(1-halo:nx+halo, 1-halo:ny+halo+shift,nz))
    allocate(x (isd:ied+shift,jsd:jed,nz), y (isd:ied,jsd:jed+shift,nz) )
    allocate(x1(isd:ied+shift,jsd:jed,nz), y1(isd:ied,jsd:jed+shift,nz) )
    allocate(x2(isd:ied+shift,jsd:jed,nz), y2(isd:ied,jsd:jed+shift,nz) )
    allocate(x3(isd:ied+shift,jsd:jed,nz), y3(isd:ied,jsd:jed+shift,nz) )
    allocate(x4(isd:ied+shift,jsd:jed,nz), y4(isd:ied,jsd:jed+shift,nz) )
    if( type == 'Masked' .OR. type =='Folded') then
       global1 = global
       global2 = global
    else
       global1 = 0.0
       global2 = 0.0
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx+1
                global1(i,j,k) = k + i*1e-3 + j*1e-6
             end do
          end do
          do j = 1,ny+1
             do i = 1,nx
                global2(i,j,k) = k + i*1e-3 + j*1e-6
             end do
          end do
       end do
    endif

    if(type == 'Masked' .or. type == 'Masked symmetry') then
       global1(nx-nx/layout(1)+1:nx+shift,ny-ny/layout(2)+1:ny,:) = 0
       global2(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny+shift,:) = 0
    end if

    select case (type)
    case ('Folded', 'Masked')
       !fill in folded north edge, cyclic east and west edge
       global1(1-halo:0,    1:ny,:) = global1(nx-halo+1:nx,1:ny,:)
       global1(nx+1:nx+halo,1:ny,:) = global1(1:halo,      1:ny,:)
       global1(1-halo:nx+halo-1,ny+1:ny+halo,:) = -global1(nx+halo-1:1-halo:-1,ny:ny-halo+1:-1,:)
       global1(nx+halo,         ny+1:ny+halo,:) = -global1(nx-halo,ny:ny-halo+1:-1,:)
       global1(1-halo:nx+halo,      1-halo:0,:) = 0
       global2(1-halo:0,    1:ny,:) = global2(nx-halo+1:nx,1:ny,:)
       global2(nx+1:nx+halo,1:ny,:) = global2(1:halo,      1:ny,:)
       global2(1-halo:nx+halo,ny+1:ny+halo,:) = -global2(nx+halo:1-halo:-1,ny-1:ny-halo:-1,:)
       global2(1-halo:nx+halo,      1-halo:0,:) = 0
    case ('Folded symmetry')
       global1(1-halo:0,                1:ny,:) = global1(nx-halo+1:nx,1:ny,:)
       global1(nx+1:nx+halo+1,          1:ny,:) = global1(1:halo+1,    1:ny,:)
       global1(1-halo:nx+halo+1,ny+1:ny+halo,:) = -global1(nx+halo+1:1-halo:-1,ny:ny-halo+1:-1,:)
       global1(1-halo:nx+halo+1,    1-halo:0,:) = 0
       global2(1-halo:0,              1:ny+1,:) = global2(nx-halo+1:nx,                 1:ny+1,:)
       global2(nx+1:nx+halo,          1:ny+1,:) = global2(1:halo,                       1:ny+1,:)
       global2(1-halo:nx+halo,ny+2:ny+halo+1,:) = -global2(nx+halo:1-halo:-1,  ny:ny-halo+1:-1,:)
       global2(1-halo:nx+halo,      1-halo:0,:) = 0
    end select

    x = 0.; y = 0.
    x(is:ie+shift,js:je,      :) = global1(is:ie+shift,js:je,      :)
    y(is:ie      ,js:je+shift,:) = global2(is:ie,      js:je+shift,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    id = mpp_clock_id( type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=CGRID_NE, complete=.true., dc_handle=dch )
    call mpp_update_domains( x3, y3, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=CGRID_NE, complete=.true., dc_handle=dch )
    call mpp_clock_end  (id)
    dch => NULL()

    !redundant points must be equal and opposite
    global2(nx/2+1:nx,    ny+shift,:) = -global2(nx/2:1:-1, ny+shift,:)
    global2(1-halo:0,     ny+shift,:) = -global2(nx-halo+1:nx, ny+shift,:)
    global2(nx+1:nx+halo, ny+shift,:) = -global2(1:halo,      ny+shift,:)

    call compare_checksums( x,  global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X' )
    call compare_checksums( y,  global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y' )
    call compare_checksums( x1, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X1' )
    call compare_checksums( x2, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X2' )
    call compare_checksums( x3, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X3' )
    call compare_checksums( x4, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X4' )
    call compare_checksums( y1, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y1' )
    call compare_checksums( y2, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y2' )
    call compare_checksums( y3, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y3' )
    call compare_checksums( y4, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y4' )

!!$    !--- delete the communicator
!!$    if(type == 'Simple' .or. type == 'Cyclic' .or. type == 'Folded') then
!!$       call mpp_update_domains( x1, domain, free=.true., list_size=4 )
!!$    else
!!$       call mpp_update_domains( x1, domain, free=.true., list_size=2 )
!!$    endif
    deallocate(global1, global2, x, x1, x2, x3, x4, y, y1, y2, y3, y4)


  end subroutine test_halo_update


  subroutine test_global_field( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, gcheck
    type(domain2D) :: domain
    type(DomainCommunicator2D), pointer, save :: dch =>NULL()
    real, allocatable    :: global1(:,:,:)
    integer              :: ishift, jshift, ni, nj, i, j
    integer, allocatable :: pelist(:)

    !--- set up domain    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name=type, symmetry = .true. )
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
        
    !--- determine if an extra point is needed
    ishift = 0; jshift = 0
    select case(type)
    case ('Symmetry corner')
       ishift = 1; jshift = 1
    case ('Symmetry east')
       ishift = 1; jshift = 0
    case ('Symmetry north')
       ishift = 0; jshift = 1
    end select

    ie  = ie+ishift;  je  = je+jshift
    ied = ied+ishift; jed = jed+jshift
    ni  = nx+ishift;  nj  = ny+jshift
    allocate(global1(1-halo:ni+halo, 1-halo:nj+halo, nz))
    global1 = 0.0
    do k = 1,nz
       do j = 1,nj
          do i = 1,ni
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    enddo

    allocate( gcheck(ni, nj, nz) )
    allocate( x (isd:ied,jsd:jed,nz) )

    x(:,:,:) = global1(isd:ied,jsd:jed,:)

    !--- test the data on data domain
    gcheck = 0.    
    id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on data domain' )

    !--- Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !--- will be declared. But for the x-direction global field, mpp_sync_self will
    !--- be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !--- in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !--- deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !--- will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !--- on all pe is needed for those partial pelist. But for y-update, it is ok. 
    !--- because the pelist in y-update is not continous.
    allocate(pelist(0:layout(1)-1))    
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !xupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags = XUPDATE )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
                            type//' mpp_global_field xupdate only on data domain' )

    !yupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags = YUPDATE )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
                            type//' mpp_global_field yupdate only on data domain' )

    gcheck=0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, new=.true., dc_handle=dch )
    call mpp_clock_end  (id)                 
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field_new on data domain' )

    gcheck=0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, new=.true., dc_handle=dch )
    dch =>NULL()
    call mpp_clock_end  (id)                                          
    !compare checksums between global and x arrays  
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, &
                            type//' mpp_global_fld_new w/ dom comm handle on data domain' )

    !--- test the data on compute domain
    gcheck = 0.    
    id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je, :), gcheck )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on compute domain' )

    !xupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, flags = XUPDATE )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
                            type//' mpp_global_field xupdate only on compute domain' )

    !yupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, flags = YUPDATE )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
                            type//' mpp_global_field yupdate only on compute domain' )

    gcheck=0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, new=.true., dc_handle=dch )
    call mpp_clock_end  (id)                 
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field_new on compute domain' )

    gcheck=0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, new=.true., dc_handle=dch )
    dch =>NULL()
    call mpp_clock_end  (id)                                          
    !compare checksums between global and x arrays  
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, &
                            type//' mpp_global_fld_new w/ dom comm handle on compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field

    !--- test mpp_global_sum, mpp_global_min and mpp_global_max
  subroutine test_global_reduce (type)
    character(len=*), intent(in) :: type
    real    :: lsum, gsum, lmax, gmax, lmin, gmin
    integer :: ni, nj, ishift, jshift
    type(domain2D) :: domain
    real, allocatable, dimension(:,:,:) :: global1, x, gcheck
    real, allocatable, dimension(:,:)   :: global2D
    !--- set up domain    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Simple' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name=type )
    case( 'Simple symmetry center', 'Simple symmetry corner', 'Simple symmetry east', 'Simple symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name=type, symmetry = .true. )
    case( 'Cyclic symmetry center', 'Cyclic symmetry corner', 'Cyclic symmetry east', 'Cyclic symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name=type, symmetry = .true., &
                                    xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN )
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
        
    !--- determine if an extra point is needed
    ishift = 0; jshift = 0
    select case(type)
    case ('Simple symmetry corner', 'Cyclic symmetry corner')
       ishift = 1; jshift = 1
    case ('Simple symmetry east', 'Cyclic symmetry east' )
       ishift = 1; jshift = 0
    case ('Simple symmetry north', 'Cyclic symmetry north')
       ishift = 0; jshift = 1
    end select

    ie  = ie+ishift;  je  = je+jshift
    ied = ied+ishift; jed = jed+jshift
    ni  = nx+ishift;  nj  = ny+jshift
    allocate(global1(1-halo:ni+halo, 1-halo:nj+halo, nz))
    global1 = 0.0
    do k = 1,nz
       do j = 1,nj
          do i = 1,ni
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    enddo

    !--- NOTE: even though the domain is cyclic, no need to apply cyclic condition on the global data

    allocate( gcheck(ni, nj, nz) )
    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( global2D(ni,nj))

    x(:,:,:) = global1(isd:ied,jsd:jed,:)
    do j = 1, nj
       do i = 1, ni
          global2D(i,j) = sum(global1(i,j,:))
       enddo 
    enddo
    !test mpp_global_sum
   
    if(type(1:6) == 'Simple') then
       gsum = sum( global2D(1:ni,1:nj) )
    else
       gsum = sum( global2D(1+ishift:ni, 1+jshift:nj) )
    endif
    id = mpp_clock_id( type//' sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lsum = mpp_global_sum( domain, x )
    call mpp_clock_end  (id)
    if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', type//' Fast sum=', lsum, gsum

    !test exact mpp_global_sum
    id = mpp_clock_id( type//' exact sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lsum = mpp_global_sum( domain, x, BITWISE_EXACT_SUM )
    call mpp_clock_end  (id)
    !--- The following check will fail on altix in normal mode, but it is ok
    !--- in debugging mode. It is ok on irix.
!    call compare_data_scalar(lsum, gsum, FATAL, type//' mpp_global_exact_sum')
    if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', type//' Bitwise-exact sum=', lsum, gsum

    !test mpp_global_min
    gmin = minval(global1(1:ni, 1:nj, :))
    id = mpp_clock_id( type//' min', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lmin = mpp_global_min( domain, x )
    call mpp_clock_end  (id)
    call compare_data_scalar(lmin, gmin, FATAL, type//' mpp_global_min')

    !test mpp_global_max
    gmax = maxval(global1(1:ni, 1:nj, :))
    id = mpp_clock_id( type//' max', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lmax = mpp_global_max( domain, x )
    call mpp_clock_end  (id)
    call compare_data_scalar(lmax, gmax, FATAL, type//' mpp_global_max' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_reduce


  subroutine test_parallel ( )
  
    integer :: npes, layout(2), i, j, k,is, ie, js, je, isd, ied, jsd, jed
    real, dimension(:,:), allocatable :: field, lfield
    real, dimension(:,:,:), allocatable :: field3d, lfield3d
    type(domain2d) :: domain
    integer, dimension(:), allocatable :: pelist1 , pelist2
    logical :: group1, group2
    character(len=128)  :: mesg
    
    npes = mpp_npes()
    allocate(pelist1(npes-mpes), pelist2(mpes))
    pelist1 = (/(i, i = 0, npes-mpes -1)/)
    pelist2 = (/(i, i = npes-mpes, npes - 1)/)
    call mpp_declare_pelist(pelist1)
    call mpp_declare_pelist(pelist2)
    group1 = .FALSE. ; group2 = .FALSE.
    if(any(pelist1==pe)) group1 = .TRUE.
    if(any(pelist2==pe)) group2 = .TRUE.
    mesg = 'parallel checking'
    
  if(group1) then
     call mpp_set_current_pelist(pelist1)
     call mpp_define_layout( (/1,nx,1,ny/), npes-mpes, layout )
  else if(group2) then
     call mpp_set_current_pelist(pelist2)
     call mpp_define_layout( (/1,nx,1,ny/), mpes, layout )
  endif
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain,&
                                   xhalo=xhalo, yhalo=yhalo)
                                   
     call mpp_set_current_pelist() 
     
     call mpp_get_compute_domain(domain, is, ie, js, je)
     call mpp_get_data_domain(domain, isd, ied, jsd, jed)
     allocate(lfield(is:ie,js:je),field(isd:ied,jsd:jed))
     allocate(lfield3d(is:ie,js:je,nz),field3d(isd:ied,jsd:jed,nz))
     
     do i = is, ie
     do j = js, je
        lfield(i,j) = real(i)+real(j)*0.001
     enddo
     enddo
     do i = is, ie
     do j = js, je
     do k = 1, nz
        lfield3d(i,j,k) = real(i)+real(j)*0.001+real(k)*0.00001
     enddo
     enddo
     enddo
     field = 0.0
     field3d = 0.0
     field(is:ie,js:je)= lfield(is:ie,js:je)
     field3d(is:ie,js:je,:) = lfield3d(is:ie,js:je,:)
     call mpp_update_domains(field,domain)
     call mpp_update_domains(field3d,domain)
     
    call mpp_check_field(field, pelist1, pelist2,domain, '2D '//mesg, w_halo = xhalo, &
                            s_halo = yhalo, e_halo = xhalo, n_halo = yhalo)
    call mpp_check_field(field3d, pelist1, pelist2,domain, '3D '//mesg, w_halo = xhalo, &
                            s_halo = yhalo, e_halo = xhalo, n_halo = yhalo)
                            
  end subroutine test_parallel
  
  subroutine test_modify_domain( )
  
    type(domain2D) :: domain2d_no_halo, domain2d_with_halo
    integer :: is,ie,js,je,isd,ied,jsd,jed
    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain2d_no_halo,   &
                            yflags=CYCLIC_GLOBAL_DOMAIN, xhalo=0, yhalo=0)
    call mpp_get_compute_domain(domain2d_no_halo, is, ie, js, je)
    call mpp_get_data_domain(domain2d_no_halo, isd, ied, jsd, jed)
    print*, "at pe ", mpp_pe(), "the compute domain decomposition of domain without halo", is, ie, js, je
    print*, "at pe ", mpp_pe(), "the data domain decomposition of domain without halo", isd, ied, jsd, jed
    call mpp_modify_domain(domain2d_no_halo, domain2d_with_halo, xhalo=xhalo,yhalo=yhalo)
    call mpp_get_compute_domain(domain2d_with_halo, is, ie, js, je)
    call mpp_get_data_domain(domain2d_with_halo, isd, ied, jsd, jed)
    print*, "at pe ", mpp_pe(), "the compute domain decomposition of domain with halo", is, ie, js, je
    print*, "at pe ", mpp_pe(), "the data domain decomposition of domain with halo", isd, ied, jsd, jed
    
    return
    
end subroutine test_modify_domain

  subroutine compare_checksums( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(LONG_KIND) :: sum1, sum2
    integer :: i, j, k
    call mpp_sync()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
         call mpp_error(FATAL,'compare_chksum: size of a and b does not match')

    do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,k) .ne. b(i,j,k)) then
                print*," at pe", mpp_pe(), "at point (",i,", ", j, ", ", k, ")", ", a = ", &
                     a(i,j,k), ", b = ", b(i,j,k)
                call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
             endif
          enddo
       enddo
    enddo

    sum1 = mpp_chksum( a, (/pe/) )
    sum2 = mpp_chksum( b, (/pe/) )

    if( sum1.EQ.sum2 )then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
        !--- in some case, even though checksum agree, the two arrays 
        !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
        !--- hence we need to check the value point by point.
    else
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums

  subroutine compare_data_scalar( a, b, action, string )
    real,             intent(in) :: a, b
    integer,          intent(in) :: action
    character(len=*), intent(in) :: string
    if( a .EQ. b)then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': data comparison are OK.' )
    else
        print*,' on pe ', mpp_pe(),' a = ', a, ', b = ', b, ', a - b =', a-b
        call mpp_error( action, trim(string)//': data comparison are not OK.' )
    end if

  end subroutine compare_data_scalar

  subroutine test_get_neighbor_1d
    type(domain1d) :: dmn1d
    integer npes, peN, peS
    npes = mpp_npes()
    call mpp_define_domains((/1,npes/), npes, dmn1d)
    call mpp_get_neighbor_pe(dmn1d, direction=+1, pe=peN)
    call mpp_get_neighbor_pe(dmn1d, direction=-1, pe=peS)
    print '(a,i2,a,2i3)', 'PE: ', mpp_pe(), ' R/L pes: ', peN, peS
  end subroutine test_get_neighbor_1d

  subroutine test_get_neighbor_non_cyclic
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    if( npes .NE. 8 ) then 
       call mpp_error(NOTE, 'mpp_domains_test: test_get_neighbor_non_cyclic '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i2,a,2i2,a,8i3)','PE: ', mpp_pe(), ' layout (non-cyclic): ', layout,  &
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_non_cyclic

  subroutine test_get_neighbor_cyclic
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    if( npes .NE. 8 ) then 
       call mpp_error(NOTE, 'mpp_domains_test: test_get_neighbor_cyclic '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
         xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i2,a,2i2,a,8i3)','PE: ', mpp_pe(), ' layout (cyclic)    : ', layout, & 
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_cyclic

  subroutine test_get_neighbor_folded_north
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    if( npes .NE. 8 ) then 
       call mpp_error(NOTE, 'mpp_domains_test: test_get_neighbor_folded_north '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
         xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i2,a,2i2,a,8i3)','PE: ', mpp_pe(), ' layout (folded N)  : ', layout, & 
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_folded_north

  subroutine test_get_neighbor_mask
    logical, allocatable ::  mask(:,:)
    integer :: im, jm, n_remove
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    
    n_remove = 2
    if( npes .NE. 8 ) then 
       call mpp_error(NOTE, 'mpp_domains_test: test_get_neighbor_mask '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes+n_remove, layout )
    allocate(mask(layout(1), layout(2)))
    mask = .TRUE.  ! activate domains
    im = min(layout(1), ceiling(layout(1)/2.0))
    jm = min(layout(2), ceiling(layout(2)/2.0))
    mask(im  ,jm  ) = .FALSE. ! deactivate domain
    mask(im  ,jm-1) = .FALSE. ! deactivate domain
    print '(a,2i3,a,2i3)', 'Masked out domains ', im, jm, ' and ', im,jm-1
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
         maskmap=mask)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i3,a,2i3,a,8i3)','PE: ', mpp_pe(), ' layout (mask   )  : ', layout, & 
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_mask


end program mpp_domains_test
#endif


! <INFO>

!   <COMPILER NAME="">     
!     Any module or program unit using <TT>mpp_domains_mod</TT>
!     must contain the line

!     <PRE>
!     use mpp_domains_mod
!     </PRE>

!     <TT>mpp_domains_mod</TT> <TT>use</TT>s <LINK
!     SRC="mpp.html">mpp_mod</LINK>, and therefore is subject to the <LINK
!     SRC="mpp.html#COMPILING AND LINKING SOURCE">compiling and linking requirements of that module.</LINK>
!   </COMPILER>
!   <PRECOMP FLAG="">      
!     <TT>mpp_domains_mod</TT> uses standard f90, and has no special
!     requirements. There are some OS-dependent
!     pre-processor directives that you might need to modify on
!     non-SGI/Cray systems and compilers. The <LINK
!     SRC="mpp.html#PORTABILITY">portability of mpp_mod</LINK>
!     obviously is a constraint, since this module is built on top of
!     it. Contact me, Balaji, SGI/GFDL, with questions.
!   </PRECOMP> 
!   <LOADER FLAG="">       
!     The <TT>mpp_domains</TT> source consists of the main source file
!     <TT>mpp_domains.F90</TT> and also requires the following include files:
!    <PRE>
!     <TT>fms_platform.h</TT>
!     <TT>mpp_update_domains2D.h</TT>
!     <TT>mpp_global_reduce.h</TT>
!     <TT>mpp_global_sum.h</TT>
!     <TT>mpp_global_field.h</TT>
!    </PRE>
!    GFDL users can check it out of the main CVS repository as part of
!    the <TT>mpp</TT> CVS module. The current public tag is <TT>galway</TT>.
!    External users can download the latest <TT>mpp</TT> package <LINK SRC=
!    "ftp://ftp.gfdl.gov/pub/vb/mpp/mpp.tar.Z">here</LINK>. Public access
!    to the GFDL CVS repository will soon be made available.

!   </LOADER>

! </INFO>
