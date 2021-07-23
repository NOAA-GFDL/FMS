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
!> @defgroup mpp_domains_mod mpp_domains_mod
!> @ingroup mpp
!> @brief Domain decomposition and domain update for message-passing codes
!> @author V. Balaji SGI/GFDL Princeton University
!!
!>  A set of simple calls for domain
!!  decomposition and domain updates on rectilinear grids. It requires the
!!  module mpp.F90, upon which it is built.\n
!!  Scalable implementations of finite-difference codes are generally
!!  based on decomposing the model domain into subdomains that are
!!  distributed among processors. These domains will then be obliged to
!!  exchange data at their boundaries if data dependencies are merely
!!  nearest-neighbour, or may need to acquire information from the global
!!  domain if there are extended data dependencies, as in the spectral
!!  transform. The domain decomposition is a key operation in the
!!  development of parallel codes.\n
!!\n
!! mpp_domains_mod provides a domain decomposition and domain
!! update API for rectilinear grids, built on top of the mpp_mod API for message passing.
!! Features of mpp_domains_mod include:\n
!!\n
!! Simple, minimal API, with free access to underlying API for more complicated stuff.\n
!!\n
!! Design toward typical use in climate/weather CFD codes.\n
!!
!> @par[Domains]
!! It is assumed that domain decomposition will mainly be in 2
!! horizontal dimensions, which will in general be the two
!! fastest-varying indices. There is a separate implementation of 1D
!! decomposition on the fastest-varying index, and 1D decomposition on
!! the second index, treated as a special case of 2D decomposition, is
!! also possible. We define domain as the grid associated with a <I>task</I>.
!! We define the compute domain as the set of gridpoints that are
!! computed by a task, and the data domain as the set of points
!! that are required by the task for the calculation. There can in
!! general be more than 1 task per PE, though often
!! the number of domains is the same as the processor count. We define
!! the global domain as the global computational domain of the
!! entire model (i.e, the same as the computational domain if run on a
!! single processor). 2D domains are defined using a derived type domain2D,
!! constructed as follows (see comments in code for more details).
!!
!!      type, public :: domain_axis_spec
!!        private
!!        integer :: begin, end, size, max_size
!!        logical :: is_global
!!      end type domain_axis_spec
!!
!!      type, public :: domain1D
!!        private
!!        type(domain_axis_spec) :: compute, data, global, active
!!        logical :: mustputb, mustgetb, mustputf, mustgetf, folded
!!        type(domain1D), pointer, dimension(:) :: list
!!        integer :: pe  ! pe to which the domain is assigned
!!        integer :: pos
!!      end type domain1D
!!
!!      type, public :: domain2D
!!        private
!!        type(domain1D) :: x
!!        type(domain1D) :: y
!!        type(domain2D), pointer, dimension(:) :: list
!!        integer :: pe ! PE to which this domain is assigned
!!        integer :: pos
!!      end type domain2D
!!
!!      type(domain1D), public :: NULL_DOMAIN1D
!!      type(domain2D), public :: NULL_DOMAIN2D

!> @file
!> @brief File for @ref mpp_domains_mod

!> @addtogroup mpp_domains_mod
!> @{

module mpp_domains_mod

#if defined(use_libMPI)
  use mpi
#endif

  use mpp_parameter_mod,      only : MPP_DEBUG, MPP_VERBOSE, MPP_DOMAIN_TIME
  use mpp_parameter_mod,      only : GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, GLOBAL,CYCLIC
  use mpp_parameter_mod,      only : AGRID, BGRID_SW, BGRID_NE, CGRID_NE, CGRID_SW, DGRID_NE, DGRID_SW
  use mpp_parameter_mod,      only : FOLD_WEST_EDGE, FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE
  use mpp_parameter_mod,      only : WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE
  use mpp_parameter_mod,      only : NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, MPP_DOMAIN_TIME
  use mpp_parameter_mod,      only : CENTER, CORNER, SCALAR_PAIR, SCALAR_BIT, BITWISE_EFP_SUM
  use mpp_parameter_mod,      only : NORTH, NORTH_EAST, EAST, SOUTH_EAST
  use mpp_parameter_mod,      only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST
  use mpp_parameter_mod,      only : MAX_DOMAIN_FIELDS, NULL_PE, DOMAIN_ID_BASE
  use mpp_parameter_mod,      only : ZERO, NINETY, MINUS_NINETY, ONE_HUNDRED_EIGHTY, MAX_TILES
  use mpp_parameter_mod,      only : EVENT_SEND, EVENT_RECV, ROOT_GLOBAL
  use mpp_parameter_mod,      only : NONBLOCK_UPDATE_TAG, EDGEONLY, EDGEUPDATE
  use mpp_parameter_mod,      only : NONSYMEDGE, NONSYMEDGEUPDATE
  use mpp_data_mod,           only : mpp_domains_stack, ptr_domains_stack
  use mpp_data_mod,           only : mpp_domains_stack_nonblock, ptr_domains_stack_nonblock
  use mpp_mod,                only : mpp_pe, mpp_root_pe, mpp_npes, mpp_error, FATAL, WARNING, NOTE
  use mpp_mod,                only : stdout, stderr, stdlog, mpp_send, mpp_recv, mpp_transmit, mpp_sync_self
  use mpp_mod,                only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_mod,                only : mpp_max, mpp_min, mpp_sum, mpp_get_current_pelist, mpp_broadcast
  use mpp_mod,                only : mpp_sum_ad
  use mpp_mod,                only : mpp_sync, mpp_init, lowercase
  use mpp_mod,                only : input_nml_file, mpp_alltoall
  use mpp_mod,                only : mpp_type, mpp_byte
  use mpp_mod,                only : mpp_type_create, mpp_type_free
  use mpp_mod,                only : COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4
  use mpp_mod,                only : mpp_declare_pelist, mpp_set_current_pelist
  use mpp_memutils_mod,       only : mpp_memuse_begin, mpp_memuse_end
  use mpp_efp_mod,            only : mpp_reproducing_sum
  use platform_mod
  implicit none
  private

  !--- public paramters imported from mpp_domains_parameter_mod
  public :: GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, BGRID_NE, BGRID_SW, CGRID_NE, CGRID_SW, AGRID
  public :: DGRID_NE, DGRID_SW, FOLD_WEST_EDGE, FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE
  public :: WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE
  public :: NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, MPP_DOMAIN_TIME, BITWISE_EFP_SUM
  public :: CENTER, CORNER, SCALAR_PAIR
  public :: NORTH, NORTH_EAST, EAST, SOUTH_EAST
  public :: SOUTH, SOUTH_WEST, WEST, NORTH_WEST
  public :: ZERO, NINETY, MINUS_NINETY, ONE_HUNDRED_EIGHTY
  public :: EDGEUPDATE, NONSYMEDGEUPDATE

  !--- public data imported from mpp_data_mod
  public :: NULL_DOMAIN1D, NULL_DOMAIN2D

  public :: domain_axis_spec, domain1D, domain2D, DomainCommunicator2D
  public :: nest_domain_type, mpp_group_update_type

  !--- public interface from mpp_domains_util.h
  public :: mpp_domains_set_stack_size, mpp_get_compute_domain, mpp_get_compute_domains
  public :: mpp_get_data_domain, mpp_get_global_domain, mpp_get_domain_components
  public :: mpp_get_layout, mpp_get_pelist, operator(.EQ.), operator(.NE.)
  public :: mpp_domain_is_symmetry, mpp_domain_is_initialized
  public :: mpp_get_neighbor_pe, mpp_nullify_domain_list
  public :: mpp_set_compute_domain, mpp_set_data_domain, mpp_set_global_domain
  public :: mpp_get_memory_domain, mpp_get_domain_shift, mpp_domain_is_tile_root_pe
  public :: mpp_get_tile_id, mpp_get_domain_extents, mpp_get_current_ntile, mpp_get_ntile_count
  public :: mpp_get_tile_list
  public :: mpp_get_tile_npes, mpp_get_domain_root_pe, mpp_get_tile_pelist, mpp_get_tile_compute_domains
  public :: mpp_get_num_overlap, mpp_get_overlap
  public :: mpp_get_io_domain, mpp_get_domain_pe, mpp_get_domain_tile_root_pe
  public :: mpp_get_domain_name, mpp_get_io_domain_layout
  public :: mpp_copy_domain, mpp_set_domain_symmetry
  public :: mpp_get_update_pelist, mpp_get_update_size
  public :: mpp_get_domain_npes, mpp_get_domain_pelist
  public :: mpp_clear_group_update
  public :: mpp_group_update_initialized, mpp_group_update_is_set
  public :: mpp_get_global_domains
  public :: mpp_create_super_grid_domain

  !--- public interface from mpp_domains_reduce.h
  public :: mpp_global_field, mpp_global_max, mpp_global_min, mpp_global_sum
  public :: mpp_global_sum_tl, mpp_global_sum_ad
  !--- public interface from mpp_domains_misc.h
  public :: mpp_broadcast_domain, mpp_domains_init, mpp_domains_exit, mpp_redistribute
  public :: mpp_update_domains, mpp_check_field
  public :: mpp_start_update_domains, mpp_complete_update_domains
  public :: mpp_create_group_update, mpp_do_group_update
  public :: mpp_start_group_update, mpp_complete_group_update
  public :: mpp_reset_group_update_field
  public :: mpp_update_nest_fine, mpp_update_nest_coarse
  public :: mpp_get_boundary
  public :: mpp_update_domains_ad
  public :: mpp_get_boundary_ad
  public :: mpp_pass_SG_to_UG, mpp_pass_UG_to_SG
  !--- public interface from mpp_domains_define.h
  public :: mpp_define_layout, mpp_define_domains, mpp_modify_domain, mpp_define_mosaic
  public :: mpp_define_mosaic_pelist, mpp_define_null_domain, mpp_mosaic_defined
  public :: mpp_define_io_domain, mpp_deallocate_domain
  public :: mpp_compute_extent, mpp_compute_block_extent

  !--- public interface for unstruct domain
  public :: mpp_define_unstruct_domain, domainUG, mpp_get_UG_io_domain
  public :: mpp_get_UG_domain_npes, mpp_get_UG_compute_domain, mpp_get_UG_domain_tile_id
  public :: mpp_get_UG_domain_pelist, mpp_get_ug_domain_grid_index
  public :: mpp_get_UG_domain_ntiles, mpp_get_UG_global_domain
  public :: mpp_global_field_ug, mpp_get_ug_domain_tile_list, mpp_get_UG_compute_domains
  public :: mpp_define_null_UG_domain, NULL_DOMAINUG, mpp_get_UG_domains_index
  public :: mpp_get_UG_SG_domain, mpp_get_UG_domain_tile_pe_inf

  !--- public interface from mpp_define_domains.inc
  public :: mpp_define_nest_domains, mpp_get_C2F_index, mpp_get_F2C_index
  public :: mpp_get_nest_coarse_domain, mpp_get_nest_fine_domain
  public :: mpp_is_nest_coarse, mpp_is_nest_fine
  public :: mpp_get_nest_pelist, mpp_get_nest_npes
  public :: mpp_get_nest_fine_pelist, mpp_get_nest_fine_npes

!----------
!ug support
  public :: mpp_domain_UG_is_tile_root_pe
  public :: mpp_deallocate_domainUG
  public :: mpp_get_io_domain_UG_layout
!----------

  integer, parameter :: NAME_LENGTH = 64
  integer, parameter :: MAXLIST = 100
  integer, parameter :: MAXOVERLAP = 200
  integer, parameter :: FIELD_S = 0
  integer, parameter :: FIELD_X = 1
  integer, parameter :: FIELD_Y = 2

  !> @}

  ! data types used by mpp_domains_mod

  !> @brief Private type for axis specification data for an unstructured grid
  !> @ingroup mpp_domains_mod
  type :: unstruct_axis_spec
     private
     integer :: begin, end, size, max_size
     integer :: begin_index, end_index
  end type unstruct_axis_spec

  !> Private type for axis specification data for an unstructured domain
  !> @ingroup mpp_domains_mod
  type :: unstruct_domain_spec
     private
     type(unstruct_axis_spec) :: compute
     integer :: pe
     integer :: pos
     integer :: tile_id
  end type unstruct_domain_spec

  !> Private type
  !> @ingroup mpp_domains_mod
  type :: unstruct_overlap_type
     private
     integer :: count = 0
     integer :: pe
     integer, pointer :: i(:)=>NULL()
     integer, pointer :: j(:)=>NULL()
  end type unstruct_overlap_type

  !> Private type
  !> @ingroup mpp_domains_mod
  type :: unstruct_pass_type
     private
     integer :: nsend, nrecv
     type(unstruct_overlap_type), pointer :: recv(:)=>NULL()
     type(unstruct_overlap_type), pointer :: send(:)=>NULL()
  end type unstruct_pass_type

  !> Domain information for managing data on unstructured grids
  !> @ingroup mpp_domains_mod
  type :: domainUG
     private
     type(unstruct_axis_spec) :: compute, global !< axis specifications
     type(unstruct_domain_spec), pointer :: list(:)=>NULL() !<
     type(domainUG), pointer :: io_domain=>NULL() !<
     type(unstruct_pass_type) :: SG2UG
     type(unstruct_pass_type) :: UG2SG
     integer, pointer :: grid_index(:) => NULL() !< index of grid on current pe
     type(domain2d), pointer :: SG_domain => NULL()
     integer :: pe
     integer :: pos
     integer :: ntiles
     integer :: tile_id
     integer :: tile_root_pe
     integer :: tile_npes
     integer :: npes_io_group
     integer(i4_kind) :: io_layout
  end type domainUG

  !> Used to specify index limits along an axis of a domain
  !> @ingroup mpp_domains_mod
  type :: domain_axis_spec
     private
     integer :: begin !< start of domain axis
     integer :: end !< end of domain axis
     integer :: size !< size of domain axis
     integer :: max_size !< max size in set
     logical :: is_global !< .true. if domain axis extent covers global domain
  end type domain_axis_spec

  !> One dimensional domain used to manage shared data access between pes
  !> @ingroup mpp_domains_mod
  type :: domain1D
     private
     type(domain_axis_spec) :: compute, data, global, memory !> index limits for different domains
     logical :: cyclic
     type(domain1D), pointer :: list(:) =>NULL() !> list of each pe's domains
     integer :: pe !<PE to which this domain is assigned
     integer :: pos !< position of this PE within link list, i.e domain%list(pos)%pe = pe
     integer :: goffset, loffset !< needed for global sum
  end type domain1D

  !> Private type used to specify index limits for a domain decomposition
  !> @ingroup mpp_domains_mod
  type :: domain1D_spec
     private
     type(domain_axis_spec) :: compute
     type(domain_axis_spec) :: global
     integer                :: pos
  end type domain1D_spec

  !> @brief Private type to specify multiple index limits and pe information for a 2D domain
  !> @ingroup mpp_domains_mod
  type :: domain2D_spec
     private
     type(domain1D_spec), pointer :: x(:)       => NULL() !< x-direction domain decomposition
     type(domain1D_spec), pointer :: y(:)       => NULL() !< y-direction domain decomposition
     integer,        pointer :: tile_id(:) => NULL() !< tile id of each tile
     integer                 :: pe                   !< PE to which this domain is assigned
     integer                 :: pos                  !< position of this PE within link list
     integer                 :: tile_root_pe         !< root pe of tile.
  end type domain2D_spec

  !> Type for overlapping data
  !> @ingroup mpp_domains_mod
  type :: overlap_type
     private
     integer                  :: count = 0                 !< number of overlapping
     integer                  :: pe
     integer                  :: start_pos                 !< start position in the buffer
     integer                  :: totsize                   !< all message size
     integer ,        pointer :: msgsize(:)      => NULL() !< overlapping msgsize to be sent or received
     integer,         pointer :: tileMe(:)       => NULL() !< my tile id for this overlap
     integer,         pointer :: tileNbr(:)      => NULL() !< neighbor tile id for this overlap
     integer,         pointer :: is(:)           => NULL() !< starting i-index
     integer,         pointer :: ie(:)           => NULL() !< ending   i-index
     integer,         pointer :: js(:)           => NULL() !< starting j-index
     integer,         pointer :: je(:)           => NULL() !< ending   j-index
     integer,         pointer :: dir(:)          => NULL() !< direction ( value 1,2,3,4 = E,S,W,N)
     integer,         pointer :: rotation(:)     => NULL() !< rotation angle.
     integer,         pointer :: index(:)        => NULL() !< for refinement
     logical,         pointer :: from_contact(:) => NULL() !< indicate if the overlap is computed from define_contact_overlap
  end type overlap_type

  !> Private type for overlap specifications
  !> @ingroup mpp_domains_mod
  type :: overlapSpec
     private
     integer                     :: whalo, ehalo, shalo, nhalo !< halo size
     integer                     :: xbegin, xend, ybegin, yend
     integer                     :: nsend, nrecv
     integer                     :: sendsize, recvsize
     type(overlap_type), pointer :: send(:) => NULL()
     type(overlap_type), pointer :: recv(:) => NULL()
     type(overlapSpec),  pointer :: next => NULL()
  end type overlapSpec

  !> @brief Upper and lower x and y bounds for a tile
  !> @ingroup mpp_domains_mod
  type :: tile_type
     integer :: xbegin, xend, ybegin, yend
  end type tile_type

  !> The domain2D type contains all the necessary information to
  !! define the global, compute and data domains of each task, as well as the PE
  !! associated with the task. The PEs from which remote data may be
  !! acquired to update the data domain are also contained in a linked list of neighbours.
  !!
  !! Domain types of higher rank can be constructed from type domain1D
  !! typically we only need 1 and 2D, but could need higher (e.g 3D LES)
  !! some elements are repeated below if they are needed once per domain, not once per axis
  !> @ingroup mpp_domains_mod
  TYPE :: domain2D
     private
     character(len=NAME_LENGTH)  :: name='unnamed' !< name of the domain, default is "unspecified"
     integer(i8_kind)            :: id
     integer                     :: pe             !< PE to which this domain is assigned
     integer                     :: fold
     integer                     :: pos            !< position of this PE within link list
     logical                     :: symmetry       !< indicate the domain is symmetric or non-symmetric.
     integer                     :: whalo, ehalo   !< halo size in x-direction
     integer                     :: shalo, nhalo   !< halo size in y-direction
     integer                     :: ntiles         !< number of tiles within mosaic
     integer                     :: max_ntile_pe   !< maximum value in the pelist of number of tiles on each pe.
     integer                     :: ncontacts     !< number of contact region within mosaic.
     logical                     :: rotated_ninety  !< indicate if any contact rotate NINETY or MINUS_NINETY
     logical                     :: initialized=.FALSE. !< indicate if the overlapping is computed or not.
     integer                     :: tile_root_pe  !< root pe of current tile.
     integer                     :: io_layout(2)  !< io_layout, will be set through mpp_define_io_domain
                                                  !! default = domain layout
     integer,            pointer :: pearray(:,:)  => NULL() !< pe of each layout position
     integer,            pointer :: tile_id(:)    => NULL() !< tile id of each tile on current processor
     integer,            pointer :: tile_id_all(:)=> NULL() !< tile id of all the tiles of domain
     type(domain1D),     pointer :: x(:)          => NULL() !< x-direction domain decomposition
     type(domain1D),     pointer :: y(:)          => NULL() !< y-direction domain decomposition
     type(domain2D_spec),pointer :: list(:)       => NULL() !< domain decomposition on pe list
     type(tile_type),    pointer :: tileList(:)   => NULL() !< store tile information
     type(overlapSpec),  pointer :: check_C       => NULL() !< send and recv information for boundary consistency check of C-cell
     type(overlapSpec),  pointer :: check_E       => NULL() !< send and recv information for boundary consistency check of E-cell
     type(overlapSpec),  pointer :: check_N       => NULL() !< send and recv information for boundary consistency check of N-cell
     type(overlapSpec),  pointer :: bound_C       => NULL() !< send information for getting boundary value for symmetry domain.
     type(overlapSpec),  pointer :: bound_E       => NULL() !< send information for getting boundary value for symmetry domain.
     type(overlapSpec),  pointer :: bound_N       => NULL() !< send information for getting boundary value for symmetry domain.
     type(overlapSpec),  pointer :: update_T      => NULL() !< send and recv information for halo update of T-cell.
     type(overlapSpec),  pointer :: update_E      => NULL() !< send and recv information for halo update of E-cell.
     type(overlapSpec),  pointer :: update_C      => NULL() !< send and recv information for halo update of C-cell.
     type(overlapSpec),  pointer :: update_N      => NULL() !< send and recv information for halo update of N-cell.
     type(domain2d),     pointer :: io_domain     => NULL() !< domain for IO, will be set through calling mpp_set_io_domain ( this will be changed).
  END TYPE domain2D

  !> Type used to represent the contact between tiles.
  !> @note This type will only be used in mpp_domains_define.inc
  !> @ingroup mpp_domains_mod
  type, private :: contact_type
     private
     integer          :: ncontact                               !< number of neighbor tile.
     integer, pointer :: tile(:) =>NULL()                       !< neighbor tile
     integer, pointer :: align1(:)=>NULL(), align2(:)=>NULL()   !< alignment of me and neighbor
     real,    pointer :: refine1(:)=>NULL(), refine2(:)=>NULL() !
     integer, pointer :: is1(:)=>NULL(), ie1(:)=>NULL()         !< i-index of current tile repsenting contact
     integer, pointer :: js1(:)=>NULL(), je1(:)=>NULL()         !< j-index of current tile repsenting contact
     integer, pointer :: is2(:)=>NULL(), ie2(:)=>NULL()         !< i-index of neighbor tile repsenting contact
     integer, pointer :: js2(:)=>NULL(), je2(:)=>NULL()         !< j-index of neighbor tile repsenting contact
  end type contact_type

  !> index bounds for use in @ref nestSpec
  !> @ingroup mpp_domains_mod
  type :: index_type
     integer :: is_me, ie_me, js_me, je_me
     integer :: is_you, ie_you, js_you, je_you
  end type index_type

  !> Used to specify bounds and index information for nested tiles as a linked list
  !> @ingroup mpp_domains_mod
  type :: nestSpec
     private
     integer                     :: xbegin, xend, ybegin, yend
     integer                     :: xbegin_c, xend_c, ybegin_c, yend_c
     integer                     :: xbegin_f, xend_f, ybegin_f, yend_f
     integer                     :: xsize_c, ysize_c
     type(index_type)            :: west, east, south, north, center
     integer                     :: nsend, nrecv
     integer                     :: extra_halo
     type(overlap_type), pointer :: send(:) => NULL()
     type(overlap_type), pointer :: recv(:) => NULL()
     type(nestSpec),     pointer :: next => NULL()

  end type nestSpec

  !> @brief domain with nested fine and course tiles
  !> @ingroup mpp_domains_mod
  type :: nest_domain_type
     character(len=NAME_LENGTH)     :: name
     integer                        :: num_level
     type(nest_level_type), pointer :: nest(:) => NULL()
     integer                        :: num_nest
     integer,               pointer :: tile_fine(:), tile_coarse(:)
     integer,               pointer :: istart_fine(:), iend_fine(:), jstart_fine(:), jend_fine(:)
     integer,               pointer :: istart_coarse(:), iend_coarse(:), jstart_coarse(:), jend_coarse(:)
  end type nest_domain_type

  !> Private type to hold data for each level of nesting
  !> @ingroup mpp_domains_mod
  type :: nest_level_type
     private
     logical                    :: on_level
     logical                    :: is_fine, is_coarse
     integer                    :: num_nest
     integer                    :: my_num_nest
     integer,           pointer :: my_nest_id(:)
     integer,           pointer :: tile_fine(:), tile_coarse(:)
     integer,           pointer :: istart_fine(:), iend_fine(:), jstart_fine(:), jend_fine(:)
     integer,           pointer :: istart_coarse(:), iend_coarse(:), jstart_coarse(:), jend_coarse(:)
     integer                    :: x_refine, y_refine
     logical                    :: is_fine_pe, is_coarse_pe
     integer,           pointer :: pelist(:) => NULL()
     integer,           pointer :: pelist_fine(:) => NULL()
     integer,           pointer :: pelist_coarse(:) => NULL()
     type(nestSpec), pointer :: C2F_T => NULL()
     type(nestSpec), pointer :: C2F_C => NULL()
     type(nestSpec), pointer :: C2F_E => NULL()
     type(nestSpec), pointer :: C2F_N => NULL()
     type(nestSpec), pointer :: F2C_T => NULL()
     type(nestSpec), pointer :: F2C_C => NULL()
     type(nestSpec), pointer :: F2C_E => NULL()
     type(nestSpec), pointer :: F2C_N => NULL()
     type(domain2d), pointer :: domain_fine   => NULL()
     type(domain2d), pointer :: domain_coarse => NULL()
  end type nest_level_type


  !> Used for sending domain data between pe's
  !> @ingroup mpp_domains_mod
  type :: domaincommunicator2D
     private
     logical            :: initialized=.false.
     integer(i8_kind) :: id=-9999
     integer(i8_kind) :: l_addr  =-9999
     integer(i8_kind) :: l_addrx =-9999
     integer(i8_kind) :: l_addry =-9999
     type(domain2D), pointer :: domain     =>NULL()
     type(domain2D), pointer :: domain_in  =>NULL()
     type(domain2D), pointer :: domain_out =>NULL()
     type(overlapSpec), pointer :: send(:,:,:,:) => NULL()
     type(overlapSpec), pointer :: recv(:,:,:,:) => NULL()
     integer, dimension(:,:),       allocatable :: sendis
     integer, dimension(:,:),       allocatable :: sendie
     integer, dimension(:,:),       allocatable :: sendjs
     integer, dimension(:,:),       allocatable :: sendje
     integer, dimension(:,:),       allocatable :: recvis
     integer, dimension(:,:),       allocatable :: recvie
     integer, dimension(:,:),       allocatable :: recvjs
     integer, dimension(:,:),       allocatable :: recvje
     logical, dimension(:),         allocatable :: S_do_buf
     logical, dimension(:),         allocatable :: R_do_buf
     integer, dimension(:),         allocatable :: cto_pe
     integer, dimension(:),         allocatable :: cfrom_pe
     integer, dimension(:),         allocatable :: S_msize
     integer, dimension(:),         allocatable :: R_msize
     integer :: Slist_size=0, Rlist_size=0
     integer :: isize=0, jsize=0, ke=0
     integer :: isize_in=0, jsize_in=0
     integer :: isize_out=0, jsize_out=0
     integer :: isize_max=0, jsize_max=0
     integer :: gf_ioff=0, gf_joff=0
  ! Remote data
     integer, dimension(:)  , allocatable :: isizeR
     integer, dimension(:)  , allocatable :: jsizeR
     integer, dimension(:,:), allocatable :: sendisR
     integer, dimension(:,:), allocatable :: sendjsR
     integer(i8_kind), dimension(:), allocatable :: rem_addr
     integer(i8_kind), dimension(:), allocatable :: rem_addrx
     integer(i8_kind), dimension(:), allocatable :: rem_addry
     integer(i8_kind), dimension(:,:), allocatable :: rem_addrl
     integer(i8_kind), dimension(:,:), allocatable :: rem_addrlx
     integer(i8_kind), dimension(:,:), allocatable :: rem_addrly
     integer                             :: position !< data location. T, E, C, or N.
  end type DomainCommunicator2D

  integer, parameter :: MAX_REQUEST = 100

  !> Used for nonblocking data transfer
  !> @ingroup mpp_domains_mod
  type :: nonblock_type
     integer                         :: recv_pos
     integer                         :: send_pos
     integer                         :: recv_msgsize
     integer                         :: send_msgsize
     integer                         :: update_flags
     integer                         :: update_position
     integer                         :: update_gridtype
     integer                         :: update_whalo
     integer                         :: update_ehalo
     integer                         :: update_shalo
     integer                         :: update_nhalo
     integer                         :: request_send_count
     integer                         :: request_recv_count
     integer, dimension(MAX_REQUEST) :: request_send
     integer, dimension(MAX_REQUEST) :: request_recv
     integer, dimension(MAX_REQUEST) :: size_recv
     integer, dimension(MAX_REQUEST) :: type_recv
     integer, dimension(MAX_REQUEST) :: buffer_pos_send
     integer, dimension(MAX_REQUEST) :: buffer_pos_recv
     integer(i8_kind)              :: field_addrs(MAX_DOMAIN_FIELDS)
     integer(i8_kind)              :: field_addrs2(MAX_DOMAIN_FIELDS)
     integer                         :: nfields
  end type nonblock_type

  !> used for updates on a group
  !> @ingroup mpp_domains_mod
  type :: mpp_group_update_type
     private
     logical            :: initialized = .FALSE.
     logical            :: k_loop_inside = .TRUE.
     logical            :: nonsym_edge = .FALSE.
     integer            :: nscalar = 0
     integer            :: nvector = 0
     integer            :: flags_s=0, flags_v=0
     integer            :: whalo_s=0, ehalo_s=0, shalo_s=0, nhalo_s=0
     integer            :: isize_s=0, jsize_s=0, ksize_s=1
     integer            :: whalo_v=0, ehalo_v=0, shalo_v=0, nhalo_v=0
     integer            :: isize_x=0, jsize_x=0, ksize_v=1
     integer            :: isize_y=0, jsize_y=0
     integer            :: position=0, gridtype=0
     logical            :: recv_s(8), recv_x(8), recv_y(8)
     integer            :: is_s=0, ie_s=0, js_s=0, je_s=0
     integer            :: is_x=0, ie_x=0, js_x=0, je_x=0
     integer            :: is_y=0, ie_y=0, js_y=0, je_y=0
     integer            :: nrecv=0, nsend=0
     integer            :: npack=0, nunpack=0
     integer            :: reset_index_s = 0
     integer            :: reset_index_v = 0
     integer            :: tot_msgsize = 0
     integer            :: from_pe(MAXOVERLAP)
     integer            :: to_pe(MAXOVERLAP)
     integer            :: recv_size(MAXOVERLAP)
     integer            :: send_size(MAXOVERLAP)
     integer            :: buffer_pos_recv(MAXOVERLAP)
     integer            :: buffer_pos_send(MAXOVERLAP)
     integer            :: pack_type(MAXOVERLAP)
     integer            :: pack_buffer_pos(MAXOVERLAP)
     integer            :: pack_rotation(MAXOVERLAP)
     integer            :: pack_size(MAXOVERLAP)
     integer            :: pack_is(MAXOVERLAP)
     integer            :: pack_ie(MAXOVERLAP)
     integer            :: pack_js(MAXOVERLAP)
     integer            :: pack_je(MAXOVERLAP)
     integer            :: unpack_type(MAXOVERLAP)
     integer            :: unpack_buffer_pos(MAXOVERLAP)
     integer            :: unpack_rotation(MAXOVERLAP)
     integer            :: unpack_size(MAXOVERLAP)
     integer            :: unpack_is(MAXOVERLAP)
     integer            :: unpack_ie(MAXOVERLAP)
     integer            :: unpack_js(MAXOVERLAP)
     integer            :: unpack_je(MAXOVERLAP)
     integer(i8_kind) :: addrs_s(MAX_DOMAIN_FIELDS)
     integer(i8_kind) :: addrs_x(MAX_DOMAIN_FIELDS)
     integer(i8_kind) :: addrs_y(MAX_DOMAIN_FIELDS)
     integer            :: buffer_start_pos = -1
     integer            :: request_send(MAX_REQUEST)
     integer            :: request_recv(MAX_REQUEST)
     integer            :: type_recv(MAX_REQUEST)
  end type mpp_group_update_type

!#######################################################################

!> @addtogroup mpp_domains_mod
!> @{
!***********************************************************************
!
!     module variables
!
!***********************************************************************
  integer              :: pe
  logical              :: module_is_initialized = .false.
  logical              :: debug                 = .FALSE.
  logical              :: verbose=.FALSE.
  logical              :: mosaic_defined = .false.
  integer              :: mpp_domains_stack_size=0
  integer              :: mpp_domains_stack_hwm=0
  type(domain1D),save  :: NULL_DOMAIN1D
  type(domain2D),save  :: NULL_DOMAIN2D
  type(domainUG),save  :: NULL_DOMAINUG
  integer              :: current_id_update = 0
  integer                         :: num_update = 0
  integer                         :: num_nonblock_group_update = 0
  integer                         :: nonblock_buffer_pos = 0
  integer                         :: nonblock_group_buffer_pos = 0
  logical                         :: start_update = .true.
  logical                         :: complete_update = .false.
  type(nonblock_type), allocatable :: nonblock_data(:)
  integer, parameter              :: MAX_NONBLOCK_UPDATE = 100

  integer                         :: group_update_buffer_pos = 0
  logical                         :: complete_group_update_on = .false.
  !-------- The following variables are used in mpp_domains_comm.h

  integer, parameter :: MAX_ADDRS=512
  integer(i8_kind),dimension(MAX_ADDRS),save :: addrs_sorted=-9999 !< list of sorted local addresses
  integer,           dimension(-1:MAX_ADDRS),save :: addrs_idx=-9999 !< index of address associated with d_comm
  integer,           dimension(MAX_ADDRS),save :: a_salvage=-9999 !< freed index list of addresses
  integer,                                save :: a_sort_len=0 !< length sorted memory list
  integer,                                save :: n_addrs=0   !< number of memory addresses used

  integer(i8_kind), parameter :: ADDR2_BASE = int(Z'0000000000010000', kind=i8_kind)
  integer, parameter :: MAX_ADDRS2=128
  integer(i8_kind),dimension(MAX_ADDRS2),save :: addrs2_sorted=-9999 !< list of sorted local addresses
  integer,           dimension(-1:MAX_ADDRS2),save :: addrs2_idx=-9999 !< index of addr2 associated with d_comm
  integer,           dimension(MAX_ADDRS2),save :: a2_salvage=-9999 !< freed indices of addr2
  integer,                                 save :: a2_sort_len=0   !< length sorted memory list
  integer,                                 save :: n_addrs2=0  !< number of memory addresses used

  integer, parameter :: MAX_DOM_IDS=128
  integer(i8_kind),dimension(MAX_DOM_IDS),save :: ids_sorted=-9999 !< list of sorted domain identifiers
  integer,           dimension(-1:MAX_DOM_IDS),save :: ids_idx=-9999 !< index of d_comm associated with sorted addesses
  integer,                                  save :: i_sort_len=0 !< length sorted domain ids list
  integer,                                  save :: n_ids=0 !< number of domain ids used
                                                            !!(=i_sort_len; domain ids are never removed)

  integer, parameter :: MAX_FIELDS=1024
  integer(i8_kind),        dimension(MAX_FIELDS),save           :: dcKey_sorted=-9999 !< list of sorted local addresses
  !     Not sure why static d_comm fails during deallocation of derived type members; allocatable works
  !     type(DomainCommunicator2D),dimension(MAX_FIELDS),save,target    :: d_comm !< domain communicators
  type(DomainCommunicator2D),dimension(:),allocatable,save,target :: d_comm  !< domain communicators
  integer,                   dimension(-1:MAX_FIELDS),save           :: d_comm_idx=-9999 !< index of d_comm associated with sorted addresses
  integer,                   dimension(MAX_FIELDS),save           :: dc_salvage=-9999 !< freed indices of d_comm
  integer,                                         save           :: dc_sort_len=0 !< length sorted comm keys
!! (=num active communicators)
  integer,                                         save           :: n_comm=0  !< number of communicators used

  !     integer(i8_kind), parameter :: GT_BASE=2**8
  integer(i8_kind), parameter :: GT_BASE = int(Z'0000000000000100', kind=i8_kind)

  !     integer(i8_kind), parameter :: KE_BASE=2**48
  integer(i8_kind), parameter :: KE_BASE = int(Z'0001000000000000', kind=i8_kind)

  integer(i8_kind) :: domain_cnt=0

  !--- the following variables are used in mpp_domains_misc.h
  logical :: domain_clocks_on=.FALSE.
  integer :: send_clock=0, recv_clock=0, unpk_clock=0
  integer :: wait_clock=0, pack_clock=0
  integer :: send_pack_clock_nonblock=0, recv_clock_nonblock=0, unpk_clock_nonblock=0
  integer :: wait_clock_nonblock=0
  integer :: nest_send_clock=0, nest_recv_clock=0, nest_unpk_clock=0
  integer :: nest_wait_clock=0, nest_pack_clock=0
  integer :: group_recv_clock=0, group_send_clock=0, group_pack_clock=0, group_unpk_clock=0, group_wait_clock=0
  integer :: nonblock_group_recv_clock=0, nonblock_group_send_clock=0, nonblock_group_pack_clock=0
  integer :: nonblock_group_unpk_clock=0, nonblock_group_wait_clock=0

!> namelist interface
  character(len=32) :: debug_update_domain = "none" !< when debug_update_domain = none, no debug will be done.
                                                    !! When debug_update_domain is set to fatal,
                                                    !! the run will be exited with fatal error message
                                                    !! When debug_update_domain is set to warning,
                                                    !! the run will output warning message.
                                                    !! When debug update_domain is set to note,
                                                    !! the run will output some note message.
  logical           :: debug_message_passing = .false. !<  Will check the consistency on the boundary between
                                                       !! processor/tile when updating domain for symmetric domain and
                                                       !! check the consistency on the north folded edge.
  integer           :: nthread_control_loop = 8 !< Determine the loop order for packing and unpacking.
                                                !! When number of threads is greater than nthread_control_loop,
                                                !! the k-loop will be moved outside and combined with number
                                                !! of pack and unpack. When the number of threads is
                                                !! less than or equal to nthread_control_loop, the k-loop
                                                !! is moved inside, but still outside, of j,i loop.
  logical           :: efp_sum_overflow_check = .false. !< If .true., always do overflow_check
                                                        !! when doing EFP bitwise mpp_global_sum.
  logical           :: use_alltoallw = .false.
  namelist /mpp_domains_nml/ debug_update_domain, domain_clocks_on, debug_message_passing, nthread_control_loop, &
                             efp_sum_overflow_check, use_alltoallw

  !***********************************************************************

  integer, parameter :: NO_CHECK = -1
  integer            :: debug_update_level = NO_CHECK
!> @}
!***********************************************************************
!
!         public interface from mpp_domains_define.h
!
!***********************************************************************
  !> @brief Retrieve layout associated with a domain decomposition.
  !> Given a global 2D domain and the number of divisions in the
  !! decomposition ndivs (usually the PE count unless some
  !! domains are \e masked) this calls returns a 2D domain layout.
  !! By default, mpp_define_layout will attempt to divide the
  !! 2D index space into domains that maintain the aspect ratio of the
  !! global domain. If this cannot be done, the algorithm favours domains
  !! that are longer in \e x than \e y, a preference that could improve vector performance.
  !! <br>Example usage:
  !! @code{.F90}call mpp_define_layout( global_indices, ndivs, layout )@endcode
  !> @ingroup mpp_domains_mod
  interface mpp_define_layout
     module procedure mpp_define_layout2D
  end interface

  !> @brief Set up a domain decomposition.
  !!
  !> There are two forms for the \e mpp_define_domains call. The 2D version is generally
  !! to be used but is built by repeated calls to the 1D version, also provided.
  !!
  !! <br>Example usage:
  !!
  !!                    call mpp_define_domains( global_indices, ndivs, domain, &
  !!                                   pelist, flags, halo, extent, maskmap )
  !!                    call mpp_define_domains( global_indices, layout, domain, pelist, &
  !!                                   xflags, yflags, xhalo, yhalo,           &
  !!                                   xextent, yextent, maskmap, name )
  !!
  !! @param global_indices Defines the global domain.
  !! @param ndivs The number of domain divisions required.
  !! @param [inout] domain Holds the resulting domain decomposition.
  !! @param pelist List of PEs to which the domains are to be assigned.
  !! @param flags An optional flag to pass additional information
  !! about the desired domain topology. Useful flags in a 1D decomposition
  !! include <TT>GLOBAL_DATA_DOMAIN</TT> and
  !! <TT>CYCLIC_GLOBAL_DOMAIN</TT>. Flags are integers: multiple flags may
  !! be added together. The flag values are public parameters available by
  !! use association.
  !! @param halo Width of the halo.
  !! @param extent Normally <TT>mpp_define_domains</TT> attempts
  !! an even division of the global domain across <TT>ndivs</TT>
  !! domains. The <TT>extent</TT> array can be used by the user to pass a
  !! custom domain division. The <TT>extent</TT> array has <TT>ndivs</TT>
  !! elements and holds the compute domain widths, which should add up to
  !! cover the global domain exactly.
  !! @param maskmap Some divisions may be masked
  !! (<TT>maskmap=.FALSE.</TT>) to exclude them from the computation (e.g
  !! for ocean model domains that are all land). The <TT>maskmap</TT> array
  !! is dimensioned <TT>ndivs</TT> and contains <TT>.TRUE.</TT> values for
  !! any domain that must be <I>included</I> in the computation (default
  !! all). The <TT>pelist</TT> array length should match the number of
  !! domains included in the computation.
  !!
  !! <br>Example usage:
  !! @code{.F90}
  !!    call mpp_define_domains( (/1,100/), 10, domain, &
  !!                             flags=GLOBAL_DATA_DOMAIN+CYCLIC_GLOBAL_DOMAIN, halo=2 )
  !! @endcode
  !!
  !! defines 10 compute domains spanning the range [1,100] of the global
  !! domain. The compute domains are non-overlapping blocks of 10. All the data
  !! domains are global, and with a halo of 2 span the range [-1:102]. And
  !! since the global domain has been declared to be cyclic,
  !! domain(9)%next => domain(0) and domain(0)%prev =>
  !! domain(9). A field is allocated on the data domain, and computations proceed on
  !! the compute domain. A call to mpp_update_domains would fill in the values
  !! in the halo region:<br>
  !!<br>
  !! @code{.F90}
  !!            call mpp_get_data_domain( domain, isd, ied ) !returns -1 and 102
  !!            call mpp_get_compute_domain( domain, is, ie ) !returns (1,10) on PE 0 ...
  !!            allocate( a(isd:ied) )
  !!            do i = is,ie
  !!              a(i) = &lt;perform computations&gt;
  !!            end do
  !!            call mpp_update_domains( a, domain )
  !! @endcode
  !!<br>
  !! The call to mpp_update_domainsfills in the regions outside
  !! the compute domain. Since the global domain is cyclic, the values at
  !! \e i=(-1,0) are the same as at \e i=(99,100); and \e i=(101,102)
  !! are the same as \e i=(1,2).
  !!
  !! The 2D version is just an extension of this syntax to two dimensions.
  !!
  !! The 2D version of the above should generally be used in
  !! codes, including 1D-decomposed ones, if there is a possibility of
  !! future evolution toward 2D decomposition. The arguments are similar to
  !! the 1D case, except that now we have optional arguments
  !! flags, halo, extent and maskmap along two axes.
  !!
  !! flags can now take an additional possible value to fold one or more edges.
  !! This is done by using flags \e FOLD_WEST_EDGE, \e FOLD_EAST_EDGE, \e FOLD_SOUTH_EDGE or
  !! \e FOLD_NORTH_EDGE. When a fold exists (e.g cylindrical domain),
  !! vector fields reverse sign upon
  !! crossing the fold. This parity reversal is performed only in the vector version of
  !! mpp_update_domains. In addition, shift operations may need to be applied to vector fields on
  !! staggered grids, also described in the vector interface to mpp_update_domains.
  !!
  !! <TT>name</TT> is the name associated with the decomposition,
  !! e.g <TT>'Ocean model'</TT>. If this argument is present,
  !! <TT>mpp_define_domains</TT> will print the domain decomposition
  !! generated to <TT>stdlog</TT>.
  !!
  !! <br>Examples:
  !!                    call mpp_define_domains( (/1,100,1,100/), (/2,2/), domain, xhalo=1 )
  !! will create the following domain layout:<br>
  !!<br>
  !! |    domain    |domain(1)|domain(2)  |domain(3)  |domain(4)    |
  !! |--------------|---------|-----------|-----------|-------------|
  !! |Compute domain|1,50,1,50|51,100,1,50|1,50,51,100|51,100,51,100|
  !! |Data domain   |0,51,1,50|50,101,1,50|0,51,51,100|50,101,51,100|
  !!
  !!
  !! Again, we allocate arrays on the data domain, perform computations
  !! on the compute domain, and call mpp_update_domains to update the halo region.
  !!
  !! If we wished to perfom a 1D decomposition along Y on the same global domain,
  !! we could use:
  !!
  !!                    call mpp_define_domains( (/1,100,1,100/), layout=(/4,1/), domain, xhalo=1 )
  !!
  !! This will create the following domain layout:<br>
  !!<br>
  !!    |    domain    |domain(1) |domain(2)  |domain(3)  |domain(4)   |
  !!    |--------------|----------|-----------|-----------|------------|
  !!    |Compute domain|1,100,1,25|1,100,26,50|1,100,51,75|1,100,76,100|
  !!    |Data domain   |0,101,1,25|0,101,26,50|0,101,51,75|1,101,76,100|
  !> @ingroup mpp_domains_mod
  interface mpp_define_domains
     module procedure mpp_define_domains1D
     module procedure mpp_define_domains2D
  end interface

  !> Defines a nullified 1D or 2D domain
  !!
  !> <br> Example usage:
  !! @code{.F90}
  !! call mpp_define_null_domain(domain)
  !! @endcode
  !> @ingroup mpp_domains_mod
  interface mpp_define_null_domain
     module procedure mpp_define_null_domain1D
     module procedure mpp_define_null_domain2D
  end interface

  !> Copy 1D or 2D domain
  !> @param domain_in Input domain to get read
  !> @param domain_out Output domain to get written to
  !> @ingroup mpp_domains_mod
  interface mpp_copy_domain
     module procedure mpp_copy_domain1D
     module procedure mpp_copy_domain2D
  end interface mpp_copy_domain
  !> Deallocate given 1D or 2D domain
  !> @ingroup mpp_domains_mod
  interface mpp_deallocate_domain
     module procedure mpp_deallocate_domain1D
     module procedure mpp_deallocate_domain2D
  end interface

  !> Modifies the extents (compute, data and global) of a given domain
  !> @ingroup mpp_domains_mod
  interface mpp_modify_domain
     module procedure mpp_modify_domain1D
     module procedure mpp_modify_domain2D
  end interface


!***********************************************************************
!
!        public interface from mpp_domains_misc.h
!
!***********************************************************************

!> Performs halo updates for a given domain.<br>
!!
!! Used to perform a halo update of a
!! domain-decomposed array on each PE. \e MPP_TYPE can be of type
!! complex, integer, logical or real of 4-byte or 8-byte kind; of rank up to 5.
!! The vector version (with two input data fields) is only present for real types.
!! For 2D domain updates, if there are halos present along both
!! x and y, we can choose to update one only, by specifying \e flags=XUPDATE or \e flags=YUPDATE.
!! In addition, one-sided updates can be performed by setting flags
!! to any combination of WUPDATE, EUPDATE, SUPDATE and NUPDATE
!! to update the west, east, north and south halos respectively.
!! Any combination of halos may be used by adding the requisite flags, e.g:
!! \e flags=XUPDATE+SUPDATE or \e flags=EUPDATE+WUPDATE+SUPDATE will update the east,
!! west and south halos.<br>
!!<br>
!! If a call to \e mpp_update_domains involves at least one E-W
!! halo and one N-S halo, the corners involved will also be updated, i.e,
!! in the example above, the SE and SW corners will be updated.<br>
!! If \e flags is not supplied, that is equivalent to \e flags=XUPDATE+YUPDATE.<br>
!!<br>
!! The vector version is passed the \e x and \e y components of a vector field in tandem,
!! and both are updated upon return. They are passed together to treat parity issues on various
!! grids. For example, on a cubic sphere projection, the \e x \e y components may be
!! interchanged when passing from an equatorial cube face to a polar face.
!! For grids with folds, vector components change sign on crossing the fold. Paired scalar
!! quantities can also be passed with the vector version if \e flags=SCALAR_PAIR, in which
!! case components are appropriately interchanged, but signs are not.<br>
!!<br>
!!    Special treatment at boundaries such as folds is also required for
!!    staggered grids. The following types of staggered grids are
!!    recognized:<br>
!!<br>
!!    1) AGRID: values are at grid centers.<br>
!!    2) BGRID_NE: vector fields are at the NE vertex of a grid
!!    cell, i.e: the array elements \e u(i,j)and \e v(i,j)are
!!    actually at (i,j;) with respect to the grid centers.<br>
!!    3) BGRID_SW: vector fields are at the SW vertex of a grid
!!    cell, i.e: the array elements \e u(i,j) and \e v(i,j) are
!!    actually at (i;,j;) with respect to the grid centers<br>
!!    4) CGRID_NE: vector fields are at the N and E faces of a
!!    grid cell, i.e: the array elements \e u(i,j) and \e v(i,j)
!!    are actually at (i;,j) and (i,j+&#189;) with respect to the
!!    grid centers.<br>
!!    5) CGRID_SW: vector fields are at the S and W faces of a
!!    grid cell, i.e: the array elements \e u(i,j)and \e v(i,j)
!!    are actually at (i;,j) and (i,j;) with respect to the
!!    grid centers.<br>
!!<br>
!!    The gridtypes listed above are all available by use association as
!!   integer parameters. The scalar version of \e mpp_update_domains
!!    assumes that the values of a scalar field are always at \e AGRID
!!    locations, and no special boundary treatment is required. If vector
!!    fields are at staggered locations, the optional argument
!!    \e gridtype must be appropriately set for correct treatment at
!!    boundaries.
!!<br>
!!    It is safe to apply vector field updates to the appropriate arrays
!!    irrespective of the domain topology: if the topology requires no
!!    special treatment of vector fields, specifying \e gridtype will
!!    do no harm.<br>
!!<br>
!!    \e mpp_update_domains internally buffers the date being sent
!!    and received into single messages for efficiency. A turnable internal
!!    buffer area in memory is provided for this purpose by
!!    \e mpp_domains_mod. The size of this buffer area can be set by
!!   the user by calling mpp_domains
!!   \e mpp_domains_set_stack_size.
!!
!!   Example usage:
!!              call mpp_update_domains( field, domain, flags )
!!   Update a 1D domain for the given field.
!!              call mpp_update_domains( fieldx, fieldy, domain, flags, gridtype )
!!   Update a 2D domain for the given fields.
!> @ingroup mpp_domains_mod
  interface mpp_update_domains
     module procedure mpp_update_domain2D_r8_2d
     module procedure mpp_update_domain2D_r8_3d
     module procedure mpp_update_domain2D_r8_4d
     module procedure mpp_update_domain2D_r8_5d
     module procedure mpp_update_domain2D_r8_2dv
     module procedure mpp_update_domain2D_r8_3dv
     module procedure mpp_update_domain2D_r8_4dv
     module procedure mpp_update_domain2D_r8_5dv
#ifdef OVERLOAD_C8
     module procedure mpp_update_domain2D_c8_2d
     module procedure mpp_update_domain2D_c8_3d
     module procedure mpp_update_domain2D_c8_4d
     module procedure mpp_update_domain2D_c8_5d
#endif
     module procedure mpp_update_domain2D_i8_2d
     module procedure mpp_update_domain2D_i8_3d
     module procedure mpp_update_domain2D_i8_4d
     module procedure mpp_update_domain2D_i8_5d
     module procedure mpp_update_domain2D_r4_2d
     module procedure mpp_update_domain2D_r4_3d
     module procedure mpp_update_domain2D_r4_4d
     module procedure mpp_update_domain2D_r4_5d
     module procedure mpp_update_domain2D_r4_2dv
     module procedure mpp_update_domain2D_r4_3dv
     module procedure mpp_update_domain2D_r4_4dv
     module procedure mpp_update_domain2D_r4_5dv
#ifdef OVERLOAD_C4
     module procedure mpp_update_domain2D_c4_2d
     module procedure mpp_update_domain2D_c4_3d
     module procedure mpp_update_domain2D_c4_4d
     module procedure mpp_update_domain2D_c4_5d
#endif
     module procedure mpp_update_domain2D_i4_2d
     module procedure mpp_update_domain2D_i4_3d
     module procedure mpp_update_domain2D_i4_4d
     module procedure mpp_update_domain2D_i4_5d
  end interface

!> Interface to start halo updates
!! \e mpp_start_update_domains is used to start a halo update of a
!!    domain-decomposed array on each PE. \e MPP_TYPE_ can be of type
!!    \e complex, \e integer, \e logical or \e real;
!!    of 4-byte or 8-byte kind; of rank up to 5. The vector version (with
!!    two input data fields) is only present for \ereal types.<br>
!!<br>
!!    \empp_start_update_domains must be paired together with
!!    \empp_complete_update_domains. In \e mpp_start_update_domains,
!!    a buffer will be pre-post to receive (non-blocking) the
!!    data and data on computational domain will be packed and sent (non-blocking send)
!!    to other processor. In \e mpp_complete_update_domains, buffer will
!!    be unpacked to fill the halo and mpp_sync_self will be called to
!!    to ensure communication safe at the last call of mpp_complete_update_domains.<br>
!!<br>
!!    Each mpp_update_domains can be replaced by the combination of mpp_start_update_domains
!!    and mpp_complete_update_domains. The arguments in mpp_start_update_domains
!!    and mpp_complete_update_domains should be the exact the same as in
!!    mpp_update_domains to be replaced except no optional argument "complete".
!!    The following are examples on how to replace mpp_update_domains with
!!    mpp_start_update_domains/mpp_complete_update_domains
!!
!> @par Example 1: Replace one scalar mpp_update_domains.<br>
!!<br>
!!    Replace<br>
!!<br>
!!              call mpp_update_domains(data, domain, flags=update_flags)<br>
!!
!!    with<br>
!!<br>
!!              id_update = mpp_start_update_domains(data, domain, flags=update_flags)<br>
!!              ...( doing some computation )<br>
!!              call mpp_complete_update_domains(id_update, data, domain, flags=update_flags)<br>
!!
!> @par Example 2: Replace group scalar mpp_update_domains<br>
!!<br>
!!    Replace<br>
!!<br>
!!              call mpp_update_domains(data_1, domain, flags=update_flags, complete=.false.)<br>
!!              .... ( other n-2 call mpp_update_domains with complete = .false. )<br>
!!              call mpp_update_domains(data_n, domain, flags=update_flags, complete=.true. )<br>
!!<br>
!!    With<br>
!!<br>
!!              id_up_1 = mpp_start_update_domains(data_1, domain, flags=update_flags)<br>
!!              .... ( other n-2 call mpp_start_update_domains )<br>
!!              id_up_n = mpp_start_update_domains(data_n, domain, flags=update_flags)<br>
!!<br>
!!              ..... ( doing some computation )<br>
!!<br>
!!              call mpp_complete_update_domains(id_up_1, data_1, domain, flags=update_flags)<br>
!!              .... ( other n-2 call mpp_complete_update_domains  )<br>
!!              call mpp_complete_update_domains(id_up_n, data_n, domain, flags=update_flags)<br>
!!
!> @par Example 3: Replace group CGRID_NE vector, mpp_update_domains<br>
!!<br>
!!    Replace<br>
!!<br>
!!              call mpp_update_domains(u_1, v_1, domain, flags=update_flgs, gridtype=CGRID_NE, complete=.false.)<br>
!!              .... ( other n-2 call mpp_update_domains with complete = .false. )<br>
!!              call mpp_update_domains(u_1, v_1, domain, flags=update_flags, gridtype=CGRID_NE, complete=.true. )<br>
!!<br>
!!    with<br>
!!<br>
!!              id_up_1 = mpp_start_update_domains(u_1, v_1, domain, flags=update_flags, gridtype=CGRID_NE)<br>
!!              .... ( other n-2 call mpp_start_update_domains )<br>
!!              id_up_n = mpp_start_update_domains(u_n, v_n, domain, flags=update_flags, gridtype=CGRID_NE)<br>
!!<br>
!!              ..... ( doing some computation )<br>
!!<br>
!!              call mpp_complete_update_domains(id_up_1, u_1, v_1, domain, flags=update_flags, gridtype=CGRID_NE)<br>
!!              .... ( other n-2 call mpp_complete_update_domains  )<br>
!!              call mpp_complete_update_domains(id_up_n, u_n, v_n, domain, flags=update_flags, gridtype=CGRID_NE)<br>
!!<br>
!!    For 2D domain updates, if there are halos present along both
!!   \e x and \e y, we can choose to update one only, by
!!    specifying \e flags=XUPDATE or \e flags=YUPDATE. In
!!    addition, one-sided updates can be performed by setting \e flags
!!    to any combination of \e WUPDATE, \e EUPDATE,
!!    \e SUPDATE and \e NUPDATE, to update the west, east, north
!!    and south halos respectively. Any combination of halos may be used by
!!   adding the requisite flags, e.g: \e flags=XUPDATE+SUPDATE or
!!    \e flags=EUPDATE+WUPDATE+SUPDATE will update the east, west and
!!    south halos.<br>
!!<br>
!!   If a call to \e mpp_start_update_domains/mpp_complete_update_domains involves at least one E-W
!!    halo and one N-S halo, the corners involved will also be updated, i.e,
!!    in the example above, the SE and SW corners will be updated.<br>
!!<br>
!!    If \e flags is not supplied, that is
!!    equivalent to \e flags=XUPDATE+YUPDATE.<br>
!!<br>
!!   The vector version is passed the \e x and \e y
!!   components of a vector field in tandem, and both are updated upon
!!    return. They are passed together to treat parity issues on various
!!    grids. For example, on a cubic sphere projection, the \e x and
!!    \e y components may be interchanged when passing from an
!!    equatorial cube face to a polar face. For grids with folds, vector
!!    components change sign on crossing the fold.  Paired scalar quantities
!!    can also be passed with the vector version if flags=SCALAR_PAIR, in which
!!    case components are appropriately interchanged, but signs are not.<br>
!!<br>
!!    Special treatment at boundaries such as folds is also required for
!!    staggered grids. The following types of staggered grids are
!!    recognized:
!!<br>
!!    1) \e AGRID: values are at grid centers.<br>
!!    2) \e BGRID_NE: vector fields are at the NE vertex of a grid
!!    cell, i.e: the array elements \e u(i,j) and \e v(i,j) are
!!    actually at (i+&#189;,j+&#189;) with respect to the grid centers.<br>
!!    3) \e BGRID_SW: vector fields are at the SW vertex of a grid
!!    cell, i.e., the array elements \e u(i,j) and \e v(i,j) are
!!    actually at (i-&#189;,j-&#189;) with respect to the grid centers.<br>
!!    4) \e CGRID_NE: vector fields are at the N and E faces of a
!!    grid cell, i.e: the array elements \e u(i,j) and \e v(i,j)
!!    are actually at (i+&#189;,j) and (i,j+&#189;) with respect to the
!!    grid centers.<br>
!!    5) \e CGRID_SW: vector fields are at the S and W faces of a
!!    grid cell, i.e: the array elements \e u(i,j) and \e v(i,j)
!!    are actually at (i-&#189;,j) and (i,j-&#189;) with respect to the
!!    grid centers.<br>
!!<br>
!!    The gridtypes listed above are all available by use association as
!!    integer parameters. If vector fields are at staggered locations, the
!!    optional argument \e gridtype must be appropriately set for
!!    correct treatment at boundaries.
!!<br>
!!    It is safe to apply vector field updates to the appropriate arrays
!!    irrespective of the domain topology: if the topology requires no
!!    special treatment of vector fields, specifying \e gridtype will
!!    do no harm.<br>
!!<br>
!!    \e mpp_start_update_domains/mpp_complete_update_domains internally
!!    buffers the data being sent and received into single messages for efficiency.
!!    A turnable internal buffer area in memory is provided for this purpose by
!!    \e mpp_domains_mod. The size of this buffer area can be set by
!!    the user by calling \e mpp_domains_set_stack_size.<br>
!!    Example usage:
!!
!!                    call mpp_start_update_domains( field, domain, flags )
!!                    call mpp_complete_update_domains( field, domain, flags )
!> @ingroup mpp_domains_mod
  interface mpp_start_update_domains
     module procedure mpp_start_update_domain2D_r8_2d
     module procedure mpp_start_update_domain2D_r8_3d
     module procedure mpp_start_update_domain2D_r8_4d
     module procedure mpp_start_update_domain2D_r8_5d
     module procedure mpp_start_update_domain2D_r8_2dv
     module procedure mpp_start_update_domain2D_r8_3dv
     module procedure mpp_start_update_domain2D_r8_4dv
     module procedure mpp_start_update_domain2D_r8_5dv
#ifdef OVERLOAD_C8
     module procedure mpp_start_update_domain2D_c8_2d
     module procedure mpp_start_update_domain2D_c8_3d
     module procedure mpp_start_update_domain2D_c8_4d
     module procedure mpp_start_update_domain2D_c8_5d
#endif
     module procedure mpp_start_update_domain2D_i8_2d
     module procedure mpp_start_update_domain2D_i8_3d
     module procedure mpp_start_update_domain2D_i8_4d
     module procedure mpp_start_update_domain2D_i8_5d
     module procedure mpp_start_update_domain2D_r4_2d
     module procedure mpp_start_update_domain2D_r4_3d
     module procedure mpp_start_update_domain2D_r4_4d
     module procedure mpp_start_update_domain2D_r4_5d
     module procedure mpp_start_update_domain2D_r4_2dv
     module procedure mpp_start_update_domain2D_r4_3dv
     module procedure mpp_start_update_domain2D_r4_4dv
     module procedure mpp_start_update_domain2D_r4_5dv
#ifdef OVERLOAD_C4
     module procedure mpp_start_update_domain2D_c4_2d
     module procedure mpp_start_update_domain2D_c4_3d
     module procedure mpp_start_update_domain2D_c4_4d
     module procedure mpp_start_update_domain2D_c4_5d
#endif
     module procedure mpp_start_update_domain2D_i4_2d
     module procedure mpp_start_update_domain2D_i4_3d
     module procedure mpp_start_update_domain2D_i4_4d
     module procedure mpp_start_update_domain2D_i4_5d
  end interface

  !> Must be used after a call to @ref mpp_start_update_domains
  !! in order to complete a nonblocking domain update. See @ref mpp_start_update_domains
  !! for more info.
  !> @ingroup mpp_domains_mod
  interface mpp_complete_update_domains
     module procedure mpp_complete_update_domain2D_r8_2d
     module procedure mpp_complete_update_domain2D_r8_3d
     module procedure mpp_complete_update_domain2D_r8_4d
     module procedure mpp_complete_update_domain2D_r8_5d
     module procedure mpp_complete_update_domain2D_r8_2dv
     module procedure mpp_complete_update_domain2D_r8_3dv
     module procedure mpp_complete_update_domain2D_r8_4dv
     module procedure mpp_complete_update_domain2D_r8_5dv
#ifdef OVERLOAD_C8
     module procedure mpp_complete_update_domain2D_c8_2d
     module procedure mpp_complete_update_domain2D_c8_3d
     module procedure mpp_complete_update_domain2D_c8_4d
     module procedure mpp_complete_update_domain2D_c8_5d
#endif
     module procedure mpp_complete_update_domain2D_i8_2d
     module procedure mpp_complete_update_domain2D_i8_3d
     module procedure mpp_complete_update_domain2D_i8_4d
     module procedure mpp_complete_update_domain2D_i8_5d
     module procedure mpp_complete_update_domain2D_r4_2d
     module procedure mpp_complete_update_domain2D_r4_3d
     module procedure mpp_complete_update_domain2D_r4_4d
     module procedure mpp_complete_update_domain2D_r4_5d
     module procedure mpp_complete_update_domain2D_r4_2dv
     module procedure mpp_complete_update_domain2D_r4_3dv
     module procedure mpp_complete_update_domain2D_r4_4dv
     module procedure mpp_complete_update_domain2D_r4_5dv
#ifdef OVERLOAD_C4
     module procedure mpp_complete_update_domain2D_c4_2d
     module procedure mpp_complete_update_domain2D_c4_3d
     module procedure mpp_complete_update_domain2D_c4_4d
     module procedure mpp_complete_update_domain2D_c4_5d
#endif
     module procedure mpp_complete_update_domain2D_i4_2d
     module procedure mpp_complete_update_domain2D_i4_3d
     module procedure mpp_complete_update_domain2D_i4_4d
     module procedure mpp_complete_update_domain2D_i4_5d
  end interface

  !> Private interface used for non blocking updates
  !> @ingroup mpp_domains_mod
  interface mpp_start_do_update
     module procedure mpp_start_do_update_r8_3d
     module procedure mpp_start_do_update_r8_3dv
#ifdef OVERLOAD_C8
     module procedure mpp_start_do_update_c8_3d
#endif
     module procedure mpp_start_do_update_i8_3d
     module procedure mpp_start_do_update_r4_3d
     module procedure mpp_start_do_update_r4_3dv
#ifdef OVERLOAD_C4
     module procedure mpp_start_do_update_c4_3d
#endif
     module procedure mpp_start_do_update_i4_3d
  end interface

  !> Private interface used for non blocking updates
  !> @ingroup mpp_domains_mod
  interface mpp_complete_do_update
     module procedure mpp_complete_do_update_r8_3d
     module procedure mpp_complete_do_update_r8_3dv
#ifdef OVERLOAD_C8
     module procedure mpp_complete_do_update_c8_3d
#endif
     module procedure mpp_complete_do_update_i8_3d
     module procedure mpp_complete_do_update_r4_3d
     module procedure mpp_complete_do_update_r4_3dv
#ifdef OVERLOAD_C4
     module procedure mpp_complete_do_update_c4_3d
#endif
     module procedure mpp_complete_do_update_i4_3d
  end interface


  !> @ingroup mpp_domains_mod
  interface mpp_create_group_update
     module procedure mpp_create_group_update_r4_2d
     module procedure mpp_create_group_update_r4_3d
     module procedure mpp_create_group_update_r4_4d
     module procedure mpp_create_group_update_r4_2dv
     module procedure mpp_create_group_update_r4_3dv
     module procedure mpp_create_group_update_r4_4dv
     module procedure mpp_create_group_update_r8_2d
     module procedure mpp_create_group_update_r8_3d
     module procedure mpp_create_group_update_r8_4d
     module procedure mpp_create_group_update_r8_2dv
     module procedure mpp_create_group_update_r8_3dv
     module procedure mpp_create_group_update_r8_4dv
  end interface mpp_create_group_update

  !> @ingroup mpp_domains_mod
  interface mpp_do_group_update
     module procedure mpp_do_group_update_r4
     module procedure mpp_do_group_update_r8
  end interface mpp_do_group_update

  !> @ingroup mpp_domains_mod
  interface mpp_start_group_update
     module procedure mpp_start_group_update_r4
     module procedure mpp_start_group_update_r8
  end interface mpp_start_group_update

  !> @ingroup mpp_domains_mod
  interface mpp_complete_group_update
     module procedure mpp_complete_group_update_r4
     module procedure mpp_complete_group_update_r8
  end interface mpp_complete_group_update

  !> @ingroup mpp_domains_mod
  interface mpp_reset_group_update_field
     module procedure mpp_reset_group_update_field_r4_2d
     module procedure mpp_reset_group_update_field_r4_3d
     module procedure mpp_reset_group_update_field_r4_4d
     module procedure mpp_reset_group_update_field_r4_2dv
     module procedure mpp_reset_group_update_field_r4_3dv
     module procedure mpp_reset_group_update_field_r4_4dv
     module procedure mpp_reset_group_update_field_r8_2d
     module procedure mpp_reset_group_update_field_r8_3d
     module procedure mpp_reset_group_update_field_r8_4d
     module procedure mpp_reset_group_update_field_r8_2dv
     module procedure mpp_reset_group_update_field_r8_3dv
     module procedure mpp_reset_group_update_field_r8_4dv
  end interface mpp_reset_group_update_field

  !> Pass the data from coarse grid to fill the buffer to be ready to be interpolated
  !! nto fine grid.
  !! <br>Example usage:
  !!
  !!                call mpp_update_nest_fine(field, nest_domain, wbuffer, ebuffer, sbuffer,
  !!                            nbuffer, nest_level, flags, complete, position, extra_halo, name,
  !!                            tile_count)
  !> @ingroup mpp_domains_mod
  interface mpp_update_nest_fine
     module procedure mpp_update_nest_fine_r8_2d
     module procedure mpp_update_nest_fine_r8_3d
     module procedure mpp_update_nest_fine_r8_4d
     module procedure mpp_update_nest_fine_r8_2dv
     module procedure mpp_update_nest_fine_r8_3dv
     module procedure mpp_update_nest_fine_r8_4dv
#ifdef OVERLOAD_C8
     module procedure mpp_update_nest_fine_c8_2d
     module procedure mpp_update_nest_fine_c8_3d
     module procedure mpp_update_nest_fine_c8_4d
#endif
     module procedure mpp_update_nest_fine_i8_2d
     module procedure mpp_update_nest_fine_i8_3d
     module procedure mpp_update_nest_fine_i8_4d
     module procedure mpp_update_nest_fine_r4_2d
     module procedure mpp_update_nest_fine_r4_3d
     module procedure mpp_update_nest_fine_r4_4d
     module procedure mpp_update_nest_fine_r4_2dv
     module procedure mpp_update_nest_fine_r4_3dv
     module procedure mpp_update_nest_fine_r4_4dv
#ifdef OVERLOAD_C4
     module procedure mpp_update_nest_fine_c4_2d
     module procedure mpp_update_nest_fine_c4_3d
     module procedure mpp_update_nest_fine_c4_4d
#endif
     module procedure mpp_update_nest_fine_i4_2d
     module procedure mpp_update_nest_fine_i4_3d
     module procedure mpp_update_nest_fine_i4_4d
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_do_update_nest_fine
     module procedure mpp_do_update_nest_fine_r8_3d
     module procedure mpp_do_update_nest_fine_r8_3dv
#ifdef OVERLOAD_C8
     module procedure mpp_do_update_nest_fine_c8_3d
#endif
     module procedure mpp_do_update_nest_fine_i8_3d
     module procedure mpp_do_update_nest_fine_r4_3d
     module procedure mpp_do_update_nest_fine_r4_3dv
#ifdef OVERLOAD_C4
     module procedure mpp_do_update_nest_fine_c4_3d
#endif
     module procedure mpp_do_update_nest_fine_i4_3d
  end interface

  !> Pass the data from fine grid to fill the buffer to be ready to be interpolated
  !! onto coarse grid.
  !! <br>Example usage:
  !!
  !!               call mpp_update_nest_coarse(field, nest_domain, field_out, nest_level, complete,
  !!                                 position, name, tile_count)
  !> @ingroup mpp_domains_mod
  interface mpp_update_nest_coarse
     module procedure mpp_update_nest_coarse_r8_2d
     module procedure mpp_update_nest_coarse_r8_3d
     module procedure mpp_update_nest_coarse_r8_4d
     module procedure mpp_update_nest_coarse_r8_2dv
     module procedure mpp_update_nest_coarse_r8_3dv
     module procedure mpp_update_nest_coarse_r8_4dv
#ifdef OVERLOAD_C8
     module procedure mpp_update_nest_coarse_c8_2d
     module procedure mpp_update_nest_coarse_c8_3d
     module procedure mpp_update_nest_coarse_c8_4d
#endif
     module procedure mpp_update_nest_coarse_i8_2d
     module procedure mpp_update_nest_coarse_i8_3d
     module procedure mpp_update_nest_coarse_i8_4d
     module procedure mpp_update_nest_coarse_r4_2d
     module procedure mpp_update_nest_coarse_r4_3d
     module procedure mpp_update_nest_coarse_r4_4d
     module procedure mpp_update_nest_coarse_r4_2dv
     module procedure mpp_update_nest_coarse_r4_3dv
     module procedure mpp_update_nest_coarse_r4_4dv
#ifdef OVERLOAD_C4
     module procedure mpp_update_nest_coarse_c4_2d
     module procedure mpp_update_nest_coarse_c4_3d
     module procedure mpp_update_nest_coarse_c4_4d
#endif
     module procedure mpp_update_nest_coarse_i4_2d
     module procedure mpp_update_nest_coarse_i4_3d
     module procedure mpp_update_nest_coarse_i4_4d
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_do_update_nest_coarse
     module procedure mpp_do_update_nest_coarse_r8_3d
     module procedure mpp_do_update_nest_coarse_r8_3dv
#ifdef OVERLOAD_C8
     module procedure mpp_do_update_nest_coarse_c8_3d
#endif
     module procedure mpp_do_update_nest_coarse_i8_3d
     module procedure mpp_do_update_nest_coarse_r4_3d
     module procedure mpp_do_update_nest_coarse_r4_3dv
#ifdef OVERLOAD_C4
     module procedure mpp_do_update_nest_coarse_c4_3d
#endif
     module procedure mpp_do_update_nest_coarse_i4_3d
  end interface

  !> Get the index of the data passed from fine grid to coarse grid.
  !! <br>Example usage:
  !!
  !!            call mpp_get_F2C_index(nest_domain, is_coarse, ie_coarse, js_coarse, je_coarse,
  !!                            is_fine, ie_fine, js_fine, je_fine, nest_level, position)
  !> @ingroup mpp_domains_mod
  interface mpp_get_F2C_index
    module procedure mpp_get_F2C_index_fine
    module procedure mpp_get_F2C_index_coarse
  end interface

  !> Send domain to every pe
  !> @ingroup mpp_domains_mod
  interface mpp_broadcast_domain
    module procedure mpp_broadcast_domain_1
    module procedure mpp_broadcast_domain_2
    module procedure mpp_broadcast_domain_ug
    module procedure mpp_broadcast_domain_nest_fine
    module procedure mpp_broadcast_domain_nest_coarse
  end interface

!--------------------------------------------------------------
! for adjoint update
!--------------------------------------------------------------
  !> Similar to @ref mpp_update_domains , updates adjoint domains
  !> @ingroup mpp_domains_mod
  interface mpp_update_domains_ad
     module procedure mpp_update_domains_ad_2D_r8_2d
     module procedure mpp_update_domains_ad_2D_r8_3d
     module procedure mpp_update_domains_ad_2D_r8_4d
     module procedure mpp_update_domains_ad_2D_r8_5d
     module procedure mpp_update_domains_ad_2D_r8_2dv
     module procedure mpp_update_domains_ad_2D_r8_3dv
     module procedure mpp_update_domains_ad_2D_r8_4dv
     module procedure mpp_update_domains_ad_2D_r8_5dv
     module procedure mpp_update_domains_ad_2D_r4_2d
     module procedure mpp_update_domains_ad_2D_r4_3d
     module procedure mpp_update_domains_ad_2D_r4_4d
     module procedure mpp_update_domains_ad_2D_r4_5d
     module procedure mpp_update_domains_ad_2D_r4_2dv
     module procedure mpp_update_domains_ad_2D_r4_3dv
     module procedure mpp_update_domains_ad_2D_r4_4dv
     module procedure mpp_update_domains_ad_2D_r4_5dv
  end interface
!
  !> Private interface used for @ref mpp_update_domains
  !> @ingroup mpp_domains_mod
  interface mpp_do_update
     module procedure mpp_do_update_r8_3d
     module procedure mpp_do_update_r8_3dv
#ifdef OVERLOAD_C8
     module procedure mpp_do_update_c8_3d
#endif
     module procedure mpp_do_update_i8_3d
     module procedure mpp_do_update_r4_3d
     module procedure mpp_do_update_r4_3dv
#ifdef OVERLOAD_C4
     module procedure mpp_do_update_c4_3d
#endif
     module procedure mpp_do_update_i4_3d
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_do_check
     module procedure mpp_do_check_r8_3d
     module procedure mpp_do_check_r8_3dv
#ifdef OVERLOAD_C8
     module procedure mpp_do_check_c8_3d
#endif
     module procedure mpp_do_check_i8_3d
     module procedure mpp_do_check_r4_3d
     module procedure mpp_do_check_r4_3dv
#ifdef OVERLOAD_C4
     module procedure mpp_do_check_c4_3d
#endif
     module procedure mpp_do_check_i4_3d
  end interface

  !> Passes data from a structured grid to an unstructured grid
  !! <br>Example usage:
  !!
  !!            call mpp_pass_SG_to_UG(domain, sg_data, ug_data)
  !> @ingroup mpp_domains_mod
  interface mpp_pass_SG_to_UG
     module procedure mpp_pass_SG_to_UG_r8_2d
     module procedure mpp_pass_SG_to_UG_r8_3d
     module procedure mpp_pass_SG_to_UG_r4_2d
     module procedure mpp_pass_SG_to_UG_r4_3d
     module procedure mpp_pass_SG_to_UG_i4_2d
     module procedure mpp_pass_SG_to_UG_i4_3d
     module procedure mpp_pass_SG_to_UG_l4_2d
     module procedure mpp_pass_SG_to_UG_l4_3d
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_pass_UG_to_SG
     module procedure mpp_pass_UG_to_SG_r8_2d
     module procedure mpp_pass_UG_to_SG_r8_3d
     module procedure mpp_pass_UG_to_SG_r4_2d
     module procedure mpp_pass_UG_to_SG_r4_3d
     module procedure mpp_pass_UG_to_SG_i4_2d
     module procedure mpp_pass_UG_to_SG_i4_3d
     module procedure mpp_pass_UG_to_SG_l4_2d
     module procedure mpp_pass_UG_to_SG_l4_3d
  end interface


!!$     module procedure mpp_do_update_ad_i4_3d
!!$  end interface
!
  !> @ingroup mpp_domains_mod
  interface mpp_do_update_ad
     module procedure mpp_do_update_ad_r8_3d
     module procedure mpp_do_update_ad_r8_3dv
     module procedure mpp_do_update_ad_r4_3d
     module procedure mpp_do_update_ad_r4_3dv
  end interface
!
!> Get the boundary data for symmetric domain when the data is at C, E, or N-cell center.<br>
!! \e mpp_get_boundary is used to get the boundary data for symmetric domain
!! when the data is at C, E, or N-cell center. For cubic grid, the data should always
!! at C-cell center.
!! <br>Example usage:
!!
!!                    call mpp_get_boundary(domain, field, ebuffer, sbuffer, wbuffer, nbuffer)
!! Get boundary information from domain and field and store in buffers
!> @ingroup mpp_domains_mod
  interface mpp_get_boundary
     module procedure mpp_get_boundary_r8_2d
     module procedure mpp_get_boundary_r8_3d
!     module procedure mpp_get_boundary_r8_4d
!     module procedure mpp_get_boundary_r8_5d
     module procedure mpp_get_boundary_r8_2dv
     module procedure mpp_get_boundary_r8_3dv
!     module procedure mpp_get_boundary_r8_4dv
!     module procedure mpp_get_boundary_r8_5dv
     module procedure mpp_get_boundary_r4_2d
     module procedure mpp_get_boundary_r4_3d
!     module procedure mpp_get_boundary_r4_4d
!     module procedure mpp_get_boundary_r4_5d
     module procedure mpp_get_boundary_r4_2dv
     module procedure mpp_get_boundary_r4_3dv
!     module procedure mpp_get_boundary_r4_4dv
!     module procedure mpp_get_boundary_r4_5dv
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_get_boundary_ad
     module procedure mpp_get_boundary_ad_r8_2d
     module procedure mpp_get_boundary_ad_r8_3d
     module procedure mpp_get_boundary_ad_r8_2dv
     module procedure mpp_get_boundary_ad_r8_3dv
     module procedure mpp_get_boundary_ad_r4_2d
     module procedure mpp_get_boundary_ad_r4_3d
     module procedure mpp_get_boundary_ad_r4_2dv
     module procedure mpp_get_boundary_ad_r4_3dv
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_do_get_boundary
     module procedure mpp_do_get_boundary_r8_3d
     module procedure mpp_do_get_boundary_r8_3dv
     module procedure mpp_do_get_boundary_r4_3d
     module procedure mpp_do_get_boundary_r4_3dv
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_do_get_boundary_ad
     module procedure mpp_do_get_boundary_ad_r8_3d
     module procedure mpp_do_get_boundary_ad_r8_3dv
     module procedure mpp_do_get_boundary_ad_r4_3d
     module procedure mpp_do_get_boundary_ad_r4_3dv
  end interface

!> Reorganization of distributed global arrays.<br>
!! \e mpp_redistribute is used to reorganize a distributed array.
!! \e MPP_TYPE_can be of type \e integer, \e complex, or \e real;
!! of 4-byte or 8-byte kind; of rank up to 5.
!! <br>Example usage:
!!              call mpp_redistribute( domain_in, field_in, domain_out, field_out )
!> @ingroup mpp_domains_mod
  interface mpp_redistribute
     module procedure mpp_redistribute_r8_2D
     module procedure mpp_redistribute_r8_3D
     module procedure mpp_redistribute_r8_4D
     module procedure mpp_redistribute_r8_5D
#ifdef OVERLOAD_C8
     module procedure mpp_redistribute_c8_2D
     module procedure mpp_redistribute_c8_3D
     module procedure mpp_redistribute_c8_4D
     module procedure mpp_redistribute_c8_5D
#endif
     module procedure mpp_redistribute_i8_2D
     module procedure mpp_redistribute_i8_3D
     module procedure mpp_redistribute_i8_4D
     module procedure mpp_redistribute_i8_5D
!!$     module procedure mpp_redistribute_l8_2D
!!$     module procedure mpp_redistribute_l8_3D
!!$     module procedure mpp_redistribute_l8_4D
!!$     module procedure mpp_redistribute_l8_5D
     module procedure mpp_redistribute_r4_2D
     module procedure mpp_redistribute_r4_3D
     module procedure mpp_redistribute_r4_4D
     module procedure mpp_redistribute_r4_5D
#ifdef OVERLOAD_C4
     module procedure mpp_redistribute_c4_2D
     module procedure mpp_redistribute_c4_3D
     module procedure mpp_redistribute_c4_4D
     module procedure mpp_redistribute_c4_5D
#endif
     module procedure mpp_redistribute_i4_2D
     module procedure mpp_redistribute_i4_3D
     module procedure mpp_redistribute_i4_4D
     module procedure mpp_redistribute_i4_5D
!!$     module procedure mpp_redistribute_l4_2D
!!$     module procedure mpp_redistribute_l4_3D
!!$     module procedure mpp_redistribute_l4_4D
!!$     module procedure mpp_redistribute_l4_5D
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_do_redistribute
     module procedure mpp_do_redistribute_r8_3D
#ifdef OVERLOAD_C8
     module procedure mpp_do_redistribute_c8_3D
#endif
     module procedure mpp_do_redistribute_i8_3D
     module procedure mpp_do_redistribute_l8_3D
     module procedure mpp_do_redistribute_r4_3D
#ifdef OVERLOAD_C4
     module procedure mpp_do_redistribute_c4_3D
#endif
     module procedure mpp_do_redistribute_i4_3D
     module procedure mpp_do_redistribute_l4_3D
  end interface

!> Parallel checking between two ensembles which run on different set pes at the same time<br>
!! There are two forms for the <TT>mpp_check_field</TT> call. The 2D
!! version is generally to be used and 3D version is  built by repeated calls to the
!! 2D version.<br>
!! <br>Example usage:
!! @code{.F90}
!!     call mpp_check_field(field_in, pelist1, pelist2, domain, mesg, &
!!                                w_halo, s_halo, e_halo, n_halo, force_abort  )
!! @endcode
!! @param field_in Field to be checked
!! @param domain Domain of current pe
!! @param mesg Message to be printed out
!! @param w_halo Halo size to be checked, default is 0
!! @param s_halo Halo size to be checked, default is 0
!! @param e_halo Halo size to be checked, default is 0
!! @param n_halo Halo size to be checked, default is 0
!! @param force_abort When true, abort program when any difference found. Default is false.
!> @ingroup mpp_domains_mod
  interface mpp_check_field
     module procedure mpp_check_field_2D
     module procedure mpp_check_field_3D
  end interface

!***********************************************************************
!
!         public interface from mpp_domains_reduce.h
!
!***********************************************************************

!> Fill in a global array from domain-decomposed arrays.<br>
!!
!> <TT>mpp_global_field</TT> is used to get an entire
!! domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!! <TT>complex</TT>, <TT>integer</TT>, <TT>logical</TT> or <TT>real</TT>;
!! of 4-byte or 8-byte kind; of rank up to 5.<br>
!!
!! All PEs in a domain decomposition must call
!! <TT>mpp_global_field</TT>, and each will have a complete global field
!! at the end. Please note that a global array of rank 3 or higher could
!! occupy a lot of memory.
!!
!! @param domain 2D domain
!! @param local Data dimensioned on either the compute or data domains of 'domain'
!! @param[out] global output data dimensioned on the corresponding global domain
!! @param flags can be either XONLY or YONLY parameters to specify a globalization on one axis only
!!
!! <br> Example usage:
!! @code{.F90}
!! call mpp_global_field( domain, local, global, flags )
!! @endcode
!> @ingroup mpp_domains_mod
  interface mpp_global_field
     module procedure mpp_global_field2D_r8_2d
     module procedure mpp_global_field2D_r8_3d
     module procedure mpp_global_field2D_r8_4d
     module procedure mpp_global_field2D_r8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_global_field2D_c8_2d
     module procedure mpp_global_field2D_c8_3d
     module procedure mpp_global_field2D_c8_4d
     module procedure mpp_global_field2D_c8_5d
#endif
     module procedure mpp_global_field2D_i8_2d
     module procedure mpp_global_field2D_i8_3d
     module procedure mpp_global_field2D_i8_4d
     module procedure mpp_global_field2D_i8_5d
     module procedure mpp_global_field2D_l8_2d
     module procedure mpp_global_field2D_l8_3d
     module procedure mpp_global_field2D_l8_4d
     module procedure mpp_global_field2D_l8_5d
     module procedure mpp_global_field2D_r4_2d
     module procedure mpp_global_field2D_r4_3d
     module procedure mpp_global_field2D_r4_4d
     module procedure mpp_global_field2D_r4_5d
#ifdef OVERLOAD_C4
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

!> @ingroup mpp_domains_mod
  interface mpp_global_field_ad
     module procedure mpp_global_field2D_r8_2d_ad
     module procedure mpp_global_field2D_r8_3d_ad
     module procedure mpp_global_field2D_r8_4d_ad
     module procedure mpp_global_field2D_r8_5d_ad
#ifdef OVERLOAD_C8
     module procedure mpp_global_field2D_c8_2d_ad
     module procedure mpp_global_field2D_c8_3d_ad
     module procedure mpp_global_field2D_c8_4d_ad
     module procedure mpp_global_field2D_c8_5d_ad
#endif
     module procedure mpp_global_field2D_i8_2d_ad
     module procedure mpp_global_field2D_i8_3d_ad
     module procedure mpp_global_field2D_i8_4d_ad
     module procedure mpp_global_field2D_i8_5d_ad
     module procedure mpp_global_field2D_l8_2d_ad
     module procedure mpp_global_field2D_l8_3d_ad
     module procedure mpp_global_field2D_l8_4d_ad
     module procedure mpp_global_field2D_l8_5d_ad
     module procedure mpp_global_field2D_r4_2d_ad
     module procedure mpp_global_field2D_r4_3d_ad
     module procedure mpp_global_field2D_r4_4d_ad
     module procedure mpp_global_field2D_r4_5d_ad
#ifdef OVERLOAD_C4
     module procedure mpp_global_field2D_c4_2d_ad
     module procedure mpp_global_field2D_c4_3d_ad
     module procedure mpp_global_field2D_c4_4d_ad
     module procedure mpp_global_field2D_c4_5d_ad
#endif
     module procedure mpp_global_field2D_i4_2d_ad
     module procedure mpp_global_field2D_i4_3d_ad
     module procedure mpp_global_field2D_i4_4d_ad
     module procedure mpp_global_field2D_i4_5d_ad
     module procedure mpp_global_field2D_l4_2d_ad
     module procedure mpp_global_field2D_l4_3d_ad
     module procedure mpp_global_field2D_l4_4d_ad
     module procedure mpp_global_field2D_l4_5d_ad
  end interface

!> @ingroup mpp_domains_mod
  interface mpp_do_global_field
     module procedure mpp_do_global_field2D_r8_3d
#ifdef OVERLOAD_C8
     module procedure mpp_do_global_field2D_c8_3d
#endif
     module procedure mpp_do_global_field2D_i8_3d
     module procedure mpp_do_global_field2D_l8_3d
     module procedure mpp_do_global_field2D_r4_3d
#ifdef OVERLOAD_C4
     module procedure mpp_do_global_field2D_c4_3d
#endif
     module procedure mpp_do_global_field2D_i4_3d
     module procedure mpp_do_global_field2D_l4_3d
  end interface

  interface mpp_do_global_field_a2a
     module procedure mpp_do_global_field2D_a2a_r8_3d
#ifdef OVERLOAD_C8
     module procedure mpp_do_global_field2D_a2a_c8_3d
#endif
     module procedure mpp_do_global_field2D_a2a_i8_3d
     module procedure mpp_do_global_field2D_a2a_l8_3d
     module procedure mpp_do_global_field2D_a2a_r4_3d
#ifdef OVERLOAD_C4
     module procedure mpp_do_global_field2D_a2a_c4_3d
#endif
     module procedure mpp_do_global_field2D_a2a_i4_3d
     module procedure mpp_do_global_field2D_a2a_l4_3d
  end interface

!> @ingroup mpp_domains_mod
  interface mpp_global_field_ug
     module procedure mpp_global_field2D_ug_r8_2d
     module procedure mpp_global_field2D_ug_r8_3d
     module procedure mpp_global_field2D_ug_r8_4d
     module procedure mpp_global_field2D_ug_r8_5d
     module procedure mpp_global_field2D_ug_i8_2d
     module procedure mpp_global_field2D_ug_i8_3d
     module procedure mpp_global_field2D_ug_i8_4d
     module procedure mpp_global_field2D_ug_i8_5d
     module procedure mpp_global_field2D_ug_r4_2d
     module procedure mpp_global_field2D_ug_r4_3d
     module procedure mpp_global_field2D_ug_r4_4d
     module procedure mpp_global_field2D_ug_r4_5d
     module procedure mpp_global_field2D_ug_i4_2d
     module procedure mpp_global_field2D_ug_i4_3d
     module procedure mpp_global_field2D_ug_i4_4d
     module procedure mpp_global_field2D_ug_i4_5d
  end interface

!> @ingroup mpp_domains_mod
  interface mpp_do_global_field_ad
     module procedure mpp_do_global_field2D_r8_3d_ad
#ifdef OVERLOAD_C8
     module procedure mpp_do_global_field2D_c8_3d_ad
#endif
     module procedure mpp_do_global_field2D_i8_3d_ad
     module procedure mpp_do_global_field2D_l8_3d_ad
     module procedure mpp_do_global_field2D_r4_3d_ad
#ifdef OVERLOAD_C4
     module procedure mpp_do_global_field2D_c4_3d_ad
#endif
     module procedure mpp_do_global_field2D_i4_3d_ad
     module procedure mpp_do_global_field2D_l4_3d_ad
  end interface

!> Global max/min of domain-decomposed arrays.<br>
!! \e mpp_global_max is used to get the maximum value of a
!! domain-decomposed array on each PE. \e MPP_TYPE_can be of type
!! \e integer or \e real; of 4-byte or 8-byte kind; of rank
!! up to 5. The dimension of \e locus must equal the rank of \e field.<br>
!!<br>
!! All PEs in a domain decomposition must call \e mpp_global_max,
!! and each will have the result upon exit.
!! The function \e mpp_global_min, with an identical syntax. is also available.
!!
!! @param domain 2D domain
!! @param field field data dimensioned on either the compute or data domains of 'domain'
!! @param locus If present, van be used to retrieve the location of the maximum
!!
!! <br>Example usage:
!!              mpp_global_max( domain, field, locus )
!> @ingroup mpp_domains_mod
  interface mpp_global_max
     module procedure mpp_global_max_r8_2d
     module procedure mpp_global_max_r8_3d
     module procedure mpp_global_max_r8_4d
     module procedure mpp_global_max_r8_5d
     module procedure mpp_global_max_r4_2d
     module procedure mpp_global_max_r4_3d
     module procedure mpp_global_max_r4_4d
     module procedure mpp_global_max_r4_5d
     module procedure mpp_global_max_i8_2d
     module procedure mpp_global_max_i8_3d
     module procedure mpp_global_max_i8_4d
     module procedure mpp_global_max_i8_5d
     module procedure mpp_global_max_i4_2d
     module procedure mpp_global_max_i4_3d
     module procedure mpp_global_max_i4_4d
     module procedure mpp_global_max_i4_5d
  end interface

!> @ingroup mpp_domains_mod
  interface mpp_global_min
     module procedure mpp_global_min_r8_2d
     module procedure mpp_global_min_r8_3d
     module procedure mpp_global_min_r8_4d
     module procedure mpp_global_min_r8_5d
     module procedure mpp_global_min_r4_2d
     module procedure mpp_global_min_r4_3d
     module procedure mpp_global_min_r4_4d
     module procedure mpp_global_min_r4_5d
     module procedure mpp_global_min_i8_2d
     module procedure mpp_global_min_i8_3d
     module procedure mpp_global_min_i8_4d
     module procedure mpp_global_min_i8_5d
     module procedure mpp_global_min_i4_2d
     module procedure mpp_global_min_i4_3d
     module procedure mpp_global_min_i4_4d
     module procedure mpp_global_min_i4_5d
  end interface

!> Global sum of domain-decomposed arrays.<br>
!! \e mpp_global_sum is used to get the sum of a domain-decomposed array
!! on each PE. \e MPP_TYPE_ can be of type \e integer, \e complex, or \e real; of 4-byte or
!! 8-byte kind; of rank up to 5.
!!
!! @param domain 2D domain
!! @param field field data dimensioned on either the compute or data domain of 'domain'
!! @param flags If present must have the value BITWISE_EXACT_SUM. This produces a sum that
!! is guaranteed to produce the identical result irrespective of how the domain is decomposed.
!! This method does the sum first along the ranks beyond 2, and then calls mpp_global_field
!! to produce a global 2D array which is then summed. The default method, which is
!! considerably faster, does a local sum followed by mpp_sum across the domain
!! decomposition.
!!
!! <br>Example usage:
!!              call mpp_global_sum( domain, field, flags )
!! @note All PEs in a domain decomposition must call \e mpp_global_sum,
!! and each will have the result upon exit.
!> @ingroup mpp_domains_mod
  interface mpp_global_sum
     module procedure mpp_global_sum_r8_2d
     module procedure mpp_global_sum_r8_3d
     module procedure mpp_global_sum_r8_4d
     module procedure mpp_global_sum_r8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_global_sum_c8_2d
     module procedure mpp_global_sum_c8_3d
     module procedure mpp_global_sum_c8_4d
     module procedure mpp_global_sum_c8_5d
#endif
     module procedure mpp_global_sum_r4_2d
     module procedure mpp_global_sum_r4_3d
     module procedure mpp_global_sum_r4_4d
     module procedure mpp_global_sum_r4_5d
#ifdef OVERLOAD_C4
     module procedure mpp_global_sum_c4_2d
     module procedure mpp_global_sum_c4_3d
     module procedure mpp_global_sum_c4_4d
     module procedure mpp_global_sum_c4_5d
#endif
     module procedure mpp_global_sum_i8_2d
     module procedure mpp_global_sum_i8_3d
     module procedure mpp_global_sum_i8_4d
     module procedure mpp_global_sum_i8_5d
     module procedure mpp_global_sum_i4_2d
     module procedure mpp_global_sum_i4_3d
     module procedure mpp_global_sum_i4_4d
     module procedure mpp_global_sum_i4_5d
  end interface

!gag
!> @ingroup mpp_domains_mod
  interface mpp_global_sum_tl
     module procedure mpp_global_sum_tl_r8_2d
     module procedure mpp_global_sum_tl_r8_3d
     module procedure mpp_global_sum_tl_r8_4d
     module procedure mpp_global_sum_tl_r8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_global_sum_tl_c8_2d
     module procedure mpp_global_sum_tl_c8_3d
     module procedure mpp_global_sum_tl_c8_4d
     module procedure mpp_global_sum_tl_c8_5d
#endif
     module procedure mpp_global_sum_tl_r4_2d
     module procedure mpp_global_sum_tl_r4_3d
     module procedure mpp_global_sum_tl_r4_4d
     module procedure mpp_global_sum_tl_r4_5d
#ifdef OVERLOAD_C4
     module procedure mpp_global_sum_tl_c4_2d
     module procedure mpp_global_sum_tl_c4_3d
     module procedure mpp_global_sum_tl_c4_4d
     module procedure mpp_global_sum_tl_c4_5d
#endif
     module procedure mpp_global_sum_tl_i8_2d
     module procedure mpp_global_sum_tl_i8_3d
     module procedure mpp_global_sum_tl_i8_4d
     module procedure mpp_global_sum_tl_i8_5d
     module procedure mpp_global_sum_tl_i4_2d
     module procedure mpp_global_sum_tl_i4_3d
     module procedure mpp_global_sum_tl_i4_4d
     module procedure mpp_global_sum_tl_i4_5d
  end interface
!gag

!bnc
!> @ingroup mpp_domains_mod
  interface mpp_global_sum_ad
     module procedure mpp_global_sum_ad_r8_2d
     module procedure mpp_global_sum_ad_r8_3d
     module procedure mpp_global_sum_ad_r8_4d
     module procedure mpp_global_sum_ad_r8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_global_sum_ad_c8_2d
     module procedure mpp_global_sum_ad_c8_3d
     module procedure mpp_global_sum_ad_c8_4d
     module procedure mpp_global_sum_ad_c8_5d
#endif
     module procedure mpp_global_sum_ad_r4_2d
     module procedure mpp_global_sum_ad_r4_3d
     module procedure mpp_global_sum_ad_r4_4d
     module procedure mpp_global_sum_ad_r4_5d
#ifdef OVERLOAD_C4
     module procedure mpp_global_sum_ad_c4_2d
     module procedure mpp_global_sum_ad_c4_3d
     module procedure mpp_global_sum_ad_c4_4d
     module procedure mpp_global_sum_ad_c4_5d
#endif
     module procedure mpp_global_sum_ad_i8_2d
     module procedure mpp_global_sum_ad_i8_3d
     module procedure mpp_global_sum_ad_i8_4d
     module procedure mpp_global_sum_ad_i8_5d
     module procedure mpp_global_sum_ad_i4_2d
     module procedure mpp_global_sum_ad_i4_3d
     module procedure mpp_global_sum_ad_i4_4d
     module procedure mpp_global_sum_ad_i4_5d
  end interface
!bnc

!***********************************************************************
!
!            public interface from mpp_domain_util.h
!
!***********************************************************************
  !> @brief Retrieve PE number of a neighboring domain.
  !!
  !> Given a 1-D or 2-D domain decomposition, this call allows users to retrieve
  !! the PE number of an adjacent PE-domain while taking into account that the
  !! domain may have holes (masked) and/or have cyclic boundary conditions and/or a
  !! folded edge. Which PE-domain will be retrived will depend on "direction":
  !! +1 (right) or -1 (left) for a 1-D domain decomposition and either NORTH, SOUTH,
  !! EAST, WEST, NORTH_EAST, SOUTH_EAST, SOUTH_WEST, or NORTH_WEST for a 2-D
  !! decomposition. If no neighboring domain exists (masked domain), then the
  !! returned "pe" value will be set to NULL_PE.<br>
  !! <br>Example usage:
  !!
  !!                    call mpp_get_neighbor_pe( domain1d, direction=+1   , pe)
  !!
  !! Set pe to the neighbor pe number that is to the right of the current pe
  !!
  !!                    call mpp_get_neighbor_pe( domain2d, direction=NORTH, pe)
  !!
  !! Get neighbor pe number that's above/north of the current pe
  !> @ingroup mpp_domains_mod
  interface mpp_get_neighbor_pe
     module procedure mpp_get_neighbor_pe_1d
     module procedure mpp_get_neighbor_pe_2d
  end interface

  !> Equality/inequality operators for domaintypes. <br>
  !!
  !! <br>The module provides public operators to check for
  !! equality/inequality of domaintypes, e.g:<br>
  !!
  !!            type(domain1D) :: a, b
  !!            type(domain2D) :: c, d
  !!            ...
  !!            if( a.NE.b )then
  !!            ...
  !!            end if
  !!            if( c==d )then
  !!            ...
  !!            end if
  !!<br>
  !! Domains are considered equal if and only if the start and end
  !! indices of each of their component global, data and compute domains are equal.
  !> @ingroup mpp_domains_mod
  interface operator(.EQ.)
     module procedure mpp_domain1D_eq
     module procedure mpp_domain2D_eq
     module procedure mpp_domainUG_eq
  end interface

  !> @ingroup mpp_domains_mod
  interface operator(.NE.)
     module procedure mpp_domain1D_ne
     module procedure mpp_domain2D_ne
     module procedure mpp_domainUG_ne
  end interface

  !> These routines retrieve the axis specifications associated with the compute domains.
  !! The domain is a derived type with private elements. These routines
  !! retrieve the axis specifications associated with the compute domains
  !! The 2D version of these is a simple extension of 1D.
  !! <br>Example usage:
  !!
  !!            call mpp_get_compute_domain(domain_1D, is, ie)
  !!            call mpp_get_compute_domain(domain_2D, is, ie, js, je)
  !> @ingroup mpp_domains_mod
  interface mpp_get_compute_domain
     module procedure mpp_get_compute_domain1D
     module procedure mpp_get_compute_domain2D
  end interface

  !> Retrieve the entire array of compute domain extents associated with a decomposition.
  !!
  !> @param domain 2D domain
  !> @param[out] xbegin,ybegin x and y domain starting indices
  !> @param[out] xsize,ysize x and y domain sizes
  !! <br>Example usage:
  !!
  !!            call mpp_get_compute_domains( domain, xbegin, xend, xsize, &
  !!                                                ybegin, yend, ysize )
  !> @ingroup mpp_domains_mod
  interface mpp_get_compute_domains
     module procedure mpp_get_compute_domains1D
     module procedure mpp_get_compute_domains2D
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_get_global_domains
     module procedure mpp_get_global_domains1D
     module procedure mpp_get_global_domains2D
  end interface

  !> These routines retrieve the axis specifications associated with the data domains.
  !! The domain is a derived type with private elements. These routines
  !! retrieve the axis specifications associated with the data domains.
  !! The 2D version of these is a simple extension of 1D.
  !! <br>Example usage:
  !!
  !!                    call mpp_get_data_domain(domain_1d, isd, ied)
  !!                    call mpp_get_data_domain(domain_2d, isd, ied, jsd, jed)
  !> @ingroup mpp_domains_mod
  interface mpp_get_data_domain
     module procedure mpp_get_data_domain1D
     module procedure mpp_get_data_domain2D
  end interface

  !> These routines retrieve the axis specifications associated with the global domains.
  !! The domain is a derived type with private elements. These routines
  !! retrieve the axis specifications associated with the global domains.
  !! The 2D version of these is a simple extension of 1D.
  !! <br>Example usage:
  !!
  !!                    call mpp_get_global_domain(domain_1d, isg, ieg)
  !!                    call mpp_get_global_domain(domain_2d, isg, ieg, jsg, jeg)
  !> @ingroup mpp_domains_mod
  interface mpp_get_global_domain
     module procedure mpp_get_global_domain1D
     module procedure mpp_get_global_domain2D
  end interface

  !> These routines retrieve the axis specifications associated with the memory domains.
  !! The domain is a derived type with private elements. These routines
  !! retrieve the axis specifications associated with the memory domains.
  !! The 2D version of these is a simple extension of 1D.
  !! <br>Example usage:
  !!
  !!                    call mpp_get_memory_domain(domain_1d, ism, iem)
  !!                    call mpp_get_memory_domain(domain_2d, ism, iem, jsm, jem)
  !> @ingroup mpp_domains_mod
  interface mpp_get_memory_domain
     module procedure mpp_get_memory_domain1D
     module procedure mpp_get_memory_domain2D
  end interface

  !> @ingroup mpp_domains_mod
  interface mpp_get_domain_extents
     module procedure mpp_get_domain_extents1D
     module procedure mpp_get_domain_extents2D
  end interface

  !> These routines set the axis specifications associated with the compute domains.
  !! The domain is a derived type with private elements. These routines
  !! set the axis specifications associated with the compute domains
  !! The 2D version of these is a simple extension of 1D.
  !! <br>Example usage:
  !!
  !!                    call mpp_get_data_domain(domain_1d, isd, ied)
  !!                    call mpp_get_data_domain(domain_2d, isd, ied, jsd, jed)
  !> @ingroup mpp_domains_mod
  interface mpp_set_compute_domain
     module procedure mpp_set_compute_domain1D
     module procedure mpp_set_compute_domain2D
  end interface

  !> These routines set the axis specifications associated with the data domains.
  !! The domain is a derived type with private elements. These routines
  !! set the axis specifications associated with the data domains.
  !! The 2D version of these is a simple extension of 1D.
  !! <br>Example usage:
  !!
  !!                    call mpp_set_data_domain(domain_1d, isd, ied)
  !!                    call mpp_set_data_domain(domain_2d, isd, ied, jsd, jed)
  !> @ingroup mpp_domains_mod
  interface mpp_set_data_domain
     module procedure mpp_set_data_domain1D
     module procedure mpp_set_data_domain2D
  end interface

  !> These routines set the axis specifications associated with the global domains.
  !! The domain is a derived type with private elements. These routines
  !! set the axis specifications associated with the global domains.
  !! The 2D version of these is a simple extension of 1D.
  !! <br>Example usage:
  !!
  !!                    call mpp_set_global_domain(domain_1d, isg, ieg)
  !!                    call mpp_set_global_domain(domain_2d, isg, ieg, jsg, jeg)
  !> @ingroup mpp_domains_mod
  interface mpp_set_global_domain
     module procedure mpp_set_global_domain1D
     module procedure mpp_set_global_domain2D
  end interface

  !> Retrieve list of PEs associated with a domain decomposition.
  !! The 1D version of this call returns an array of the PEs assigned to
  !! this 1D domain decomposition. In addition the optional argument pos may be
  !! used to retrieve the 0-based position of the domain local to the
  !! calling PE, i.e., <TT> domain\%list(pos)\%pe</TT> is the local PE,
  !! as returned by @ref mpp_pe().
  !! The 2D version of this call is identical to 1D version.
  !> @ingroup mpp_domains_mod
  interface mpp_get_pelist
     module procedure mpp_get_pelist1D
     module procedure mpp_get_pelist2D
  end interface

  !> Retrieve layout associated with a domain decomposition
  !! The 1D version of this call returns the number of divisions that was assigned to this
  !! decomposition axis. The 2D version of this call returns an array of dimension 2 holding the
  !! results on two axes.
  !! <br>Example usage:
  !!
  !!                    call mpp_get_layout( domain, layout )
  !> @ingroup mpp_domains_mod
  interface mpp_get_layout
     module procedure mpp_get_layout1D
     module procedure mpp_get_layout2D
  end interface
  !> Private interface for internal usage, compares two sizes
  !> @ingroup mpp_domains_mod
  interface check_data_size
     module procedure check_data_size_1d
     module procedure check_data_size_2d
  end interface

  !> Nullify domain list. This interface is needed in mpp_domains_test.
  !! 1-D case can be added in if needed.
  !! <br>Example usage:
  !!
  !!                    call mpp_nullify_domain_list(domain)
  !> @ingroup mpp_domains_mod
  interface mpp_nullify_domain_list
     module procedure nullify_domain2d_list
  end interface

  ! Include variable "version" to be written to log file.
#include<file_version.h>
  public version


contains

#include <mpp_define_nest_domains.inc>
#include <mpp_domains_util.inc>
#include <mpp_domains_comm.inc>
#include <mpp_domains_define.inc>
#include <mpp_domains_misc.inc>
#include <mpp_domains_reduce.inc>
#include <mpp_unstruct_domain.inc>

end module mpp_domains_mod
