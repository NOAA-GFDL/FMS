!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @defgroup fms_diag_axis_object_mod fms_diag_axis_object_mod
!> @ingroup diag_manager
!! @brief Modern object-oriented implementation of diagnostic axis management for the FMS diagnostic manager.
!!
!! This module provides type-based classes for managing diagnostic axes in FMS.
!! It serves as the modern, object-oriented replacement for @ref diag_axis_mod, utilizing Fortran 2003+ class structures
!! and polymorphism.

!> @file
!> @brief File for @ref diag_axis_object_mod

!> @addtogroup fms_diag_axis_object_mod
!> @{
module fms_diag_axis_object_mod
#ifdef use_yaml
  use mpp_domains_mod, only:  domain1d, domain2d, domainUG, mpp_get_compute_domain, CENTER, &
                            & mpp_get_global_domain, NORTH, EAST, mpp_get_tile_id, &
                            & mpp_get_ntile_count, mpp_get_io_domain, mpp_get_layout
  use platform_mod,    only:  r8_kind, r4_kind, i4_kind, i8_kind
  use diag_data_mod,   only:  diag_atttype, max_axes, NO_DOMAIN, TWO_D_DOMAIN, UG_DOMAIN, &
                              direction_down, direction_up, fmsDiagAttribute_type, max_axis_attributes, &
                              DIAG_NULL, index_gridtype, latlon_gridtype, pack_size_str, &
                              get_base_year, get_base_month, get_base_day, get_base_hour, get_base_minute,&
                              get_base_second, is_x_axis, is_y_axis
  use mpp_mod,         only:  FATAL, mpp_error, uppercase, mpp_pe, mpp_root_pe, stdout
  use fms2_io_mod,     only:  FmsNetcdfFile_t, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t, &
                            & register_axis, register_field, register_variable_attribute, write_data
  use fms_diag_yaml_mod, only: subRegion_type, diag_yaml, MAX_SUBAXES, diagYamlFilesVar_type
  use diag_grid_mod,       only:  get_local_indices_cubesphere => get_local_indexes
  use axis_utils2_mod,   only: nearest_index
  implicit none

  PRIVATE

  public :: fmsDiagAxis_type, fms_diag_axis_object_init, fms_diag_axis_object_end, &
          & get_domain_and_domain_type, diagDomain_t, &
          & DIAGDOMAIN2D_T, fmsDiagSubAxis_type, fmsDiagAxisContainer_type, fmsDiagFullAxis_type, DIAGDOMAINUG_T
  public :: define_new_axis, parse_compress_att, get_axis_id_from_name, define_diurnal_axis, &
          & fmsDiagDiurnalAxis_type, create_new_z_subaxis, is_parent_axis, define_new_subaxis_latlon, &
          & define_new_subaxis_index, find_z_sub_axis_name

  !> @}

  !> @brief Base type for domain information associated with an axis.
  !!
  !! This type was created to avoid requiring separate "Domain", "Domain2", and "DomainUG" arguments
  !! in subroutines. Instead, a single polymorphic class(diagDomain_t) argument can be used, which is
  !! polymorphically extended to handle different domain types.
  !!
  !! This base provides a unified interface for domain operations regardless of whether
  !! the axis uses 1D domain decomposition, 2D domain decomposition, or an unstructured grid domain.
  !!
  !! @ingroup diag_axis_object_mod
  type diagDomain_t
    contains
      procedure :: set => set_axis_domain
      procedure :: length => get_length
      procedure :: get_ntiles
  end type diagDomain_t

  !> @brief Type to hold 1D domain decomposition information for an axis.
  !!
  !! This type extends the diagDomain_t base type and is used when an axis is
  !! associated with a 1D domain (typically a vertical or time axis).
  !! The 1D domain provides information about how the axis is partitioned across
  !! MPI processes along a single dimension.
  type, extends(diagDomain_t) :: diagDomain1d_t
     type(domain1d) :: Domain !< 1D domain object describing axis decomposition
  end type

  !> @brief Type to hold 2D domain decomposition information for an axis.
  !!
  !! This type extends the diagDomain_t base type and is used when an axis is
  !! part of a 2D domain (typically for horizontal "X" or "Y" axes in atmospheric models).
  !! The 2D domain provides information about how the axis is partitioned across
  !! MPI processes in both the X and Y dimensions.
  type, extends(diagDomain_t) :: diagDomain2d_t
    type(domain2d) :: Domain2 !< 2D domain object describing X-Y decomposition of an axis
  end type

  !> @brief Type to hold unstructured grid domain information for an axis.
  !!
  !! This type extends the diagDomain_t base type and is used when an axis is
  !! associated with an unstructured (irregular) grid domain. Unstructured grids are commonly
  !! used in models with non-uniform spatial decomposition, such as icosahedral grids.
  !! The unstructured domain provides information about cell connectivity and partitioning.
  type, extends(diagDomain_t) :: diagDomainUg_t
    type(domainUG) :: DomainUG !< Unstructured domain object for irregular mesh decomposition
  end type

  !> @brief Base type for diagnostic axis objects.
  !!
  !! This is the base type for all diagnostic axis implementations. It provides
  !! a unified interface for axis operations and is polymorphically extended by more specific
  !! axis types:
  !! - @ref fmsDiagFullAxis_type - A complete axis with coordinates
  !! - @ref fmsDiagSubAxis_type - A subset of a full axis for regional output
  !! - @ref fmsDiagDiurnalAxis_type - A specialized time-sampling axis for diurnal averaging
  !!
  !! The base type defines a common interface for querying axis properties and writing
  !! axis data to output files via the FMS2_IO library.
  !!
  !! @ingroup diag_axis_object_mod
  TYPE :: fmsDiagAxis_type
     INTEGER                        , private :: axis_id         !< Unique identifier for this axis

     contains
       procedure :: get_parent_axis_id
       procedure :: get_subaxes_id
       procedure :: get_axis_name
       procedure :: is_z_axis
       procedure :: write_axis_metadata
       procedure :: write_axis_data
       procedure :: add_structured_axis_ids
       procedure :: get_structured_axis
       procedure :: is_unstructured_grid
       procedure :: get_edges_id
  END TYPE fmsDiagAxis_type

  !> @brief Container type to hold a polymorphic diagnostic axis object.
  !!
  !! This container type allows storage of any type derived from fmsDiagAxis_type
  !! (e.g., fmsDiagFullAxis_type, fmsDiagSubAxis_type, or fmsDiagDiurnalAxis_type)
  !! in an array. This polymorphic storage approach enables dynamic type checking
  !! and method dispatch at runtime.
  !!
  !! @ingroup diag_axis_object_mod
  type :: fmsDiagAxisContainer_type
    class(fmsDiagAxis_type), allocatable :: axis !< Polymorphic axis object (Full, Sub, or Diurnal)
  end type

  !> @brief Type representing a subregion or subset of a parent diagnostic axis.
  !!
  !! A subaxis is created when a user requests output for a limited region of a full axis.
  !! This can occur for regional output (e.g., a subregion of the full domain) or for
  !! dimension compression (e.g., selected depth levels in a vertical axis).
  !!
  !! Each subaxis maintains references to its parent axis and stores the index ranges
  !! that define the subregion on the current PE, as well as globally.
  !! Subaxes may also store Z-axis bounds for identifying equivalent subaxes across files.
  !!
  !! @ingroup diag_axis_object_mod
  TYPE, extends(fmsDiagAxis_type) :: fmsDiagSubAxis_type
    CHARACTER(len=:),  ALLOCATABLE , private  :: subaxis_name   !< Name of the subaxis (typically parent_name_subNN)
    INTEGER                        , private  :: starting_index !< First index of subregion relative to parent axis
    INTEGER                        , private  :: ending_index   !< Last index of subregion relative to parent axis
    INTEGER                        , private  :: parent_axis_id !< Axis ID of the parent full axis
    INTEGER                        , private  :: compute_idx(2) !< [start, end] indices of compute domain on this PE
    INTEGER,            allocatable, private  :: global_idx(:)  !< [start, end] indices in global domain
    real(kind=r4_kind), allocatable, private  :: zbounds(:)     !< Bounds [min, max] if this is a Z-axis subregion
    contains
      procedure :: fill_subaxis
      procedure :: axis_length
      procedure :: get_starting_index
      procedure :: get_ending_index
      procedure :: get_compute_indices
      procedure :: is_same_zbounds
  END TYPE fmsDiagSubAxis_type

  !> @brief Type for diurnal (daily cycle) sampling axes.
  !!
  !! This specialized axis type represents time-of-day sampling for diurnal averaging.
  !! Diurnal axes divide a 24-hour period into regular intervals for accumulating
  !! time-averaged or instantaneous samples at each time-of-day bin. This is commonly
  !! used in climate modeling to analyze the diurnal cycle of variables.
  !!
  !! Each diurnal axis has an associated edges axis that defines the bin boundaries,
  !! and stores the actual diurnal coordinate data for output to NetCDF files.
  !!
  !! @ingroup diag_axis_object_mod
  TYPE, extends(fmsDiagAxis_type) :: fmsDiagDiurnalAxis_type
    INTEGER                      , private :: ndiurnal_samples !< Number of time-of-day samples in 24-hour period
    CHARACTER(len=:), ALLOCATABLE, private :: axis_name        !< Name of the diurnal axis (e.g., time_of_day_06)
    CHARACTER(len=:), ALLOCATABLE, private :: long_name        !< Long name for the diurnal axis
    CHARACTER(len=:), ALLOCATABLE, private :: units            !< Units string (hours since reference time)
    INTEGER                      , private :: edges_id         !< Axis ID of the diurnal edges axis
    CHARACTER(len=:), ALLOCATABLE, private :: edges_name       !< Name of the edges axis (e.g., time_of_day_edges_06)
    CLASS(*),         ALLOCATABLE, private :: diurnal_data(:)  !< Coordinate values: times within 24-hour day

    contains
      procedure :: get_diurnal_axis_samples
      procedure :: write_diurnal_metadata
  END TYPE fmsDiagDiurnalAxis_type

  !> @brief Type representing a complete diagnostic axis with coordinates and metadata.
  !!
  !! This is the primary type for storing axis information. It contains the coordinate
  !! values, metadata (name, units, long name), domain information, and optional attributes.
  !! A full axis can have subaxes defined from it for regional output or selective output.
  !!
  !! The axis stores data in either single or double precision floating point format,
  !! as determined by the input coordinate array. Domain information is stored polymorphically
  !! to support 1D, 2D, and unstructured grid domains.
  !!
  !! @ingroup diag_axis_object_mod
  TYPE, extends(fmsDiagAxis_type) :: fmsDiagFullAxis_type
     CHARACTER(len=:),   ALLOCATABLE, private :: axis_name !< Name identifier for the axis (e.g., "height", "time")
     CHARACTER(len=:),   ALLOCATABLE, private :: units !< Units string for axis values (e.g., "m", "K")
     CHARACTER(len=:),   ALLOCATABLE, private :: long_name !< Descriptive long name for the axis
                                                           !! (e.g., "Height above sea level")
                                                           !! written as "long_name" attribute in output file
     CHARACTER(len=1)               , private :: cart_name !< Cartesian classification: "X", "Y", "Z", "T", "U", "N"
     CLASS(*),           ALLOCATABLE, private :: axis_data(:)    !< Coordinate values as single or double precision
     CHARACTER(len=:),   ALLOCATABLE, private :: type_of_data    !< Data type: "float" (r4) or "double" (r8)
     !< TODO: Consider linked-list implementation to remove MAX_AXES limit
     integer,            ALLOCATABLE, private :: subaxis(:)      !< Array of axis IDs for subaxes derived from this axis
     integer                        , private :: nsubaxis        !< Count of subaxes currently defined
     class(diagDomain_t),ALLOCATABLE, private :: axis_domain     !< Domain decomposition info (1D, 2D, or UG)
     INTEGER                        , private :: type_of_domain  !< Domain type: NO_DOMAIN, TWO_D_DOMAIN, or UG_DOMAIN
     INTEGER                        , private :: length          !< Total number of coordinate points globally
     INTEGER                        , private :: direction       !< Axis direction: -1 (down), 0 (none), +1 (up)
     INTEGER,            ALLOCATABLE, private :: edges_id        !< Axis ID of edges (cell boundaries) if defined
                                                                 !! Written as coordinate in output file
     CHARACTER(len=:),   ALLOCATABLE, private :: edges_name      !< Name of the edges axis
                                                                 !! Written as "edges" attribute in file
     CHARACTER(len=:), ALLOCATABLE,   private :: aux             !< Auxiliary classification: "geolon_t" or "geolat_t"
     CHARACTER(len=128)             , private :: req             !< Required field names (comma-separated list)
     INTEGER                        , private :: tile_count      !< Number of tiles on this axis
     TYPE(fmsDiagAttribute_type),allocatable , private :: attributes(:) !< User-defined custom attributes
     INTEGER                        , private :: num_attributes  !< Number of custom attributes defined
     INTEGER                        , private :: domain_position !< Stagger position: CENTER, NORTH, or EAST
     integer, allocatable           , private :: structured_ids(:) !< Structured axis IDs for unstructured grids
                                                                   !! (maps unstructured to lat/lon axes)
     CHARACTER(len=:), ALLOCATABLE,   private :: set_name        !< Axis set name to distinguish axes with same name
                                                                 !! Used for multiple configurations of same grid

     contains

     PROCEDURE :: add_axis_attribute
     PROCEDURE :: register => register_diag_axis_obj
     PROCEDURE :: axis_length => get_axis_length
     PROCEDURE :: set_edges
     PROCEDURE :: set_axis_id
     PROCEDURE :: get_compute_domain
     PROCEDURE :: get_indices
     PROCEDURE :: get_global_io_domain
     PROCEDURE :: get_aux
     PROCEDURE :: has_aux
     PROCEDURE :: get_set_name
     PROCEDURE :: has_set_name
     PROCEDURE :: is_x_or_y_axis
     PROCEDURE :: get_dim_size_layout
  END TYPE fmsDiagFullAxis_type

  !> @addtogroup fms_diag_yaml_mod
  !> @{
  contains

  !> @brief Initialize and register a diagnostic axis with its coordinates and metadata.
  !!
  !! This subroutine configures a full diagnostic axis object with coordinate data,
  !! units, and optional domain decomposition. The axis must be initialized before
  !! being used to register diagnostic fields or create subaxes.
  !!
  !! The coordinate data is stored as-is (single or double precision), and the type
  !! is inferred from the input array to optimize storage and I/O.
  subroutine register_diag_axis_obj(this, axis_name, axis_data, units, cart_name, long_name, direction,&
  & set_name, Domain, Domain2, DomainU, aux, req, tile_count, domain_position, axis_length )
    class(fmsDiagFullAxis_type),INTENT(inout):: this            !< Diag_axis obj
    CHARACTER(len=*),   INTENT(in)           :: axis_name       !< Name of the axis
    class(*),           INTENT(in)           :: axis_data(:)    !< Array of coordinate values
    CHARACTER(len=*),   INTENT(in)           :: units           !< Units for the axis
    CHARACTER(len=1),   INTENT(in)           :: cart_name       !< Cartesian axis ("X", "Y", "Z", "T", "U", "N")
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: long_name       !< Long name for the axis.
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: set_name        !< Name of the parent axis, if it is a subaxis
    INTEGER,            INTENT(in), OPTIONAL :: direction       !< Indicates the direction of the axis
    TYPE(domain1d),     INTENT(in), OPTIONAL :: Domain          !< 1D domain
    TYPE(domain2d),     INTENT(in), OPTIONAL :: Domain2         !< 2D domain
    TYPE(domainUG),     INTENT(in), OPTIONAL :: DomainU         !< Unstructured domain
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: aux             !< Auxiliary name, can only be <TT>geolon_t</TT>
                                                                !! or <TT>geolat_t</TT>
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: req             !< Required field names.
    INTEGER,            INTENT(in), OPTIONAL :: tile_count      !< Number of tiles
    INTEGER,            INTENT(in), OPTIONAL :: domain_position !< Domain position, "NORTH" or "EAST"
    integer,            intent(in), optional :: axis_length     !< The length of the axis size(axis_data(:))

    this%axis_name = trim(axis_name)
    this%units = trim(units)
    this%cart_name = uppercase(cart_name)
    call check_if_valid_cart_name(this%cart_name)

    if (present(long_name)) this%long_name = trim(long_name)

    select type (axis_data)
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%axis_data(axis_length))
      this%axis_data = axis_data
      this%length = axis_length
      this%type_of_data = "double" !< This is what fms2_io expects in the register_field call
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%axis_data(axis_length))
      this%axis_data = axis_data
      this%length = axis_length
      this%type_of_data = "float" !< This is what fms2_io expects in the register_field call
    class default
      call mpp_error(FATAL, "The axis_data in your diag_axis_init call is not a supported type. &
                          &  Currently only r4 and r8 data is supported.")
    end select

    this%type_of_domain = NO_DOMAIN
    if (present(Domain)) then
      if (present(Domain2) .or. present(DomainU)) call mpp_error(FATAL, &
        "The presence of Domain with any other domain type is prohibited. "//&
        "Check you diag_axis_init call for axis_name:"//trim(axis_name))
      allocate(diagDomain1d_t :: this%axis_domain)
      call this%axis_domain%set(Domain=Domain)
    else if (present(Domain2)) then
        if (present(DomainU)) call mpp_error(FATAL, &
        "The presence of Domain2 with any other domain type is prohibited. "//&
        "Check you diag_axis_init call for axis_name:"//trim(axis_name))
      allocate(diagDomain2d_t :: this%axis_domain)
      call this%axis_domain%set(Domain2=Domain2)
      this%type_of_domain = TWO_D_DOMAIN
    else if (present(DomainU)) then
      allocate(diagDomainUg_t :: this%axis_domain)
      call this%axis_domain%set(DomainU=DomainU)
      this%type_of_domain = UG_DOMAIN
    endif

    this%tile_count = 1
    if (present(tile_count)) this%tile_count = tile_count

    this%domain_position = CENTER
    if (present(domain_position)) this%domain_position = domain_position
    call check_if_valid_domain_position(this%domain_position)

    this%direction = 0
    if (present(direction)) this%direction = direction
    call check_if_valid_direction(this%direction)

    if (present(aux)) this%aux = trim(aux)
    if (present(req)) this%req = trim(req)
    this%set_name = ""
    if (present(set_name)) this%set_name = trim(set_name)

    if (MAX_SUBAXES .gt. 0) then
      allocate(this%subaxis(MAX_SUBAXES))
      this%subaxis = diag_null
    endif

    this%nsubaxis = 0
    this%num_attributes = 0
  end subroutine register_diag_axis_obj

  !> @brief Add a user-defined attribute to this axis.
  !!
  !! User-defined attributes are additional metadata that can be attached to axes.
  !! These attributes will be written to the NetCDF output file as metadata.
  subroutine add_axis_attribute(this, att_name, att_value)
    class(fmsDiagFullAxis_type),INTENT(INOUT) :: this   !< diag_axis obj
    character(len=*), intent(in)    :: att_name     !< Name of the attribute
    class(*),         intent(in)    :: att_value(:) !< The attribute value to add

    integer :: j    !< obj%num_attributes (for less typing)

    if (.not. allocated(this%attributes)) &
      allocate(this%attributes(max_axis_attributes))

    this%num_attributes = this%num_attributes + 1

    j = this%num_attributes
    call this%attributes(j)%add(att_name, att_value)
  end subroutine add_axis_attribute

  !> @brief Write axis metadata (dimensions, coordinates, and attributes) to a NetCDF file.
  !!
  !! This subroutine registers the axis dimension and variable in an open FMS2_IO
  !! file object (FmsNetcdfFile_t, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t), along with all metadata
  !! attributes (units, long_name, etc.).
  !! It handles different axis types (full, sub, diurnal) and domain decomposition types.
  !!
  !! For subaxes, the parent axis information is used to determine coordinate and
  !! attribute values. The domain decomposition attribute is added when applicable.
  subroutine write_axis_metadata(this, fms2io_fileobj, edges_in_file, parent_axis)
    class(fmsDiagAxis_type),          target,  INTENT(IN)    :: this          !< diag_axis obj
    class(FmsNetcdfFile_t),                    INTENT(INOUT) :: fms2io_fileobj!< Fms2_io fileobj to write the data to
    logical,                                   INTENT(IN)    :: edges_in_file !< .True. if the edges to this axis are
                                                                              !! already in the file
    class(fmsDiagAxis_type), OPTIONAL, target, INTENT(IN)    :: parent_axis   !< If the axis is a subaxis, axis object
                                                                              !! for the parent axis (this will be used
                                                                              !! to get some of the metadata info)

    character(len=:), ALLOCATABLE         :: axis_edges_name !< Name of the edges, if it exist
    character(len=:), pointer             :: axis_name       !< Name of the axis
    integer                               :: axis_length     !< Size of the axis
    integer                               :: i               !< For do loops
    type(fmsDiagFullAxis_type), pointer   :: diag_axis       !< Local pointer to the diag_axis

    integer :: type_of_domain !< The type of domain the current axis is in
    logical :: is_subaxis     !< .true. if the axis is a subaxis
    logical :: needs_domain_decomposition !< .True. if the axis needs the domain decomposition attribute
                                          !! (i.e for "X" and "Y" subaxis)
    integer :: domain_decomposition(4) !< indices of the global (1:2) and compute (3:4) domain for a "X" and "Y" subaxis

    is_subaxis = .false.
    needs_domain_decomposition = .false.

    select type(this)
    type is (fmsDiagFullAxis_type)
      axis_name => this%axis_name
      axis_length = this%length
      diag_axis => this
      type_of_domain = this%type_of_domain
    type is (fmsDiagSubAxis_type)
      is_subaxis = .true.
      axis_name => this%subaxis_name
      axis_length = this%ending_index - this%starting_index + 1
      if (allocated(this%global_idx)) then
        needs_domain_decomposition = .true.
        domain_decomposition(1:2) = this%global_idx
        domain_decomposition(3) = this%starting_index
        domain_decomposition(4) = this%ending_index
      endif
      !< Get all the other information from the parent axis (i.e the cart_name, units, etc)
      if (present(parent_axis)) then
        select type(parent_axis)
        type is (fmsDiagFullAxis_type)
          diag_axis => parent_axis
        end select
      endif
      type_of_domain = NO_DOMAIN !< All subaxes are treated as non-domain decomposed (each rank writes it own file)
    type is (fmsDiagDiurnalAxis_type)
      call this%write_diurnal_metadata(fms2io_fileobj)
      return
    end select

    !< Add the axis as a dimension in the netcdf file based on the type of axis_domain and the fileobj type
    select type (fms2io_fileobj)
      !< The register_field calls need to be inside the select type block so that it can go inside the correct
      !! register_field interface
      type is (FmsNetcdfFile_t)
        !< Here the axis is not domain decomposed (i.e z_axis)
        call register_axis(fms2io_fileobj, axis_name, axis_length)
        call register_field(fms2io_fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
        if (needs_domain_decomposition) then
          call register_variable_attribute(fms2io_fileobj, axis_name, "domain_decomposition", &
            domain_decomposition)
        endif
      type is (FmsNetcdfDomainFile_t)
        select case (type_of_domain)
        case (NO_DOMAIN)
          !< Here the fms2io_fileobj is domain decomposed, but the axis is not
          !! Domain decomposed fileobjs can have axis that are not domain decomposed (i.e "Z" axis)
          call register_axis(fms2io_fileobj, axis_name, axis_length)
          call register_field(fms2io_fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
        case (TWO_D_DOMAIN)
          !< Here the axis is domain decomposed
          call register_axis(fms2io_fileobj, axis_name, diag_axis%cart_name, domain_position=diag_axis%domain_position)
          call register_field(fms2io_fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
        end select
      type is (FmsNetcdfUnstructuredDomainFile_t)
        select case (type_of_domain)
        case (UG_DOMAIN)
          !< Here the axis is in a unstructured domain
          call register_axis(fms2io_fileobj, axis_name)
          call register_field(fms2io_fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
        case default
          !< Here the fms2io_fileobj is in the unstructured domain, but the axis is not
          !< Unstructured domain fileobjs can have axis that are not domain decomposed (i.e "Z" axis)
          call register_axis(fms2io_fileobj, axis_name, axis_length)
          call register_field(fms2io_fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
        end select
    end select

    !< Write its metadata
    if(allocated(diag_axis%long_name)) &
      call register_variable_attribute(fms2io_fileobj, axis_name, "long_name", diag_axis%long_name, &
        str_len=len_trim(diag_axis%long_name))

    if (diag_axis%cart_name .NE. "N") &
      call register_variable_attribute(fms2io_fileobj, axis_name, "axis", diag_axis%cart_name, str_len=1)

    if (trim(diag_axis%units) .NE. "none") &
      call register_variable_attribute(fms2io_fileobj, axis_name, "units", diag_axis%units, &
                                       str_len=len_trim(diag_axis%units))

    select case (diag_axis%direction)
    case (direction_up)
      call register_variable_attribute(fms2io_fileobj, axis_name, "positive", "up", str_len=2)
    case (direction_down)
      call register_variable_attribute(fms2io_fileobj, axis_name, "positive", "down", str_len=4)
    end select

    !< Ignore the edges attribute, if the edges are already in the file or if it is subaxis
    if (.not. edges_in_file .and. allocated(diag_axis%edges_name) .and. .not. is_subaxis) then
      call register_variable_attribute(fms2io_fileobj, axis_name, "edges", diag_axis%edges_name, &
        str_len=len_trim(diag_axis%edges_name))
    endif

    if(allocated(diag_axis%attributes)) then
      do i = 1, diag_axis%num_attributes
        select type (att_value => diag_axis%attributes(i)%att_value)
        type is (character(len=*))
          call register_variable_attribute(fms2io_fileobj, axis_name, diag_axis%attributes(i)%att_name, &
                                           trim(att_value(1)), str_len=len_trim(att_value(1)))
        class default
          call register_variable_attribute(fms2io_fileobj, axis_name, diag_axis%attributes(i)%att_name, att_value)
        end select
      enddo
    endif

  end subroutine write_axis_metadata

  !> @brief Write axis coordinate data to a NetCDF file.
  !!
  !! This subroutine writes the actual coordinate values for an axis to the output file.
  !! It handles domain decomposition automatically, writing only the portion of coordinates
  !! relevant to the current MPI process's I/O domain.
  !!
  !! For subaxes, the parent axis data is sliced to the appropriate indices and written.
  !! For diurnal axes, the time-of-day coordinate values are written.
  subroutine write_axis_data(this, fms2io_fileobj, parent_axis)
    class(fmsDiagAxis_type),           target, INTENT(IN)    :: this        !< diag_axis obj
    class(FmsNetcdfFile_t),                    INTENT(INOUT) :: fms2io_fileobj!< Fms2_io fileobj to write the data to
    class(fmsDiagAxis_type), OPTIONAL, target, INTENT(IN)    :: parent_axis !< The parent axis if this is a subaxis

    integer                       :: i                 !< Starting index of a sub_axis
    integer                       :: j                 !< Ending index of a sub_axis
    integer                       :: global_io_index(2)!< Global io domain starting and ending index
    select type(this)
    type is (fmsDiagFullAxis_type)
      call this%get_global_io_domain(global_io_index, fms2io_fileobj%is_file_using_netcdf_mpi())
      call write_data(fms2io_fileobj, this%axis_name, this%axis_data(global_io_index(1):global_io_index(2)))
    type is (fmsDiagSubAxis_type)
      i = this%starting_index
      j = this%ending_index

      if (present(parent_axis)) then
        select type(parent_axis)
        type is (fmsDiagFullAxis_type)
          call write_data(fms2io_fileobj, this%subaxis_name, parent_axis%axis_data(i:j))
        end select
      endif
    type is (fmsDiagDiurnalAxis_type)
      call write_data(fms2io_fileobj, this%axis_name, this%diurnal_data)
    end select
  end subroutine write_axis_data

  !> @brief Create and registers extra axes used when performing diurnal averaging, to capture the midpoints
  !! and bounds of each diurnal sample. This subroutine will be called twice for each diurnal reduction: once to
  !! create the edges axis (time_of_day_edges_<N>) and once to create the center axis (time_of_day_<N>), N being the
  !! number of diurnal samples. The number of diurnal samples is specified by the reduction method name in your
  !! diag_table.yaml, ie. "diurnal3" for 3 diurnal samples, "diurnal24" for 24 diurnal samples, etc.
  !!
  !! The time_of_day_<N> axis will have the midpoint of the sampled segment and the time_of_day_edges_<N> will have
  !! the bounds. This will be written out in hours of a day, so edges will always start at 0 and end at 24, and
  !! minutes will be represented as decimals.
  !! 
  !! For example, if n_diurnal_samples = 3, the time_of_day_03 axis will have the values [4, 12, 20], representing the
  !! time at the midpoint (4:00 am, 12:00 pm, 8:00 pm) for each of the 3 samples and the time_of_day_edges_03 axis
  !! will have the values [0, 8, 16, 24], representing the start/end times of each sample (12:00 am, 8:00 am, 4:00 pm,
  !! 12:00 pm).
  !! For hourly sampling (n_diurnal_samples = 24), the time_of_day_24 axis will have the values
  !! [0.5, 1.5, 2.5, ..., 23.5] and the time_of_day_edges_24 axis will have the values
  !! [0, 1, 2, ..., 23, 24].
  subroutine define_diurnal_axis(diag_axis, naxis, n_diurnal_samples, is_edges)
    class(fmsDiagAxisContainer_type), target, intent(inout) :: diag_axis(:)      !< Array of axis containers
    integer,                                  intent(inout) :: naxis             !< Number of axis that have
                                                                                 !! been defined
    integer,                                  intent(in)    :: n_diurnal_samples !< The number of diurnal samples
                                                                                 !! for the curent axis
    logical,                                  intent(in)    :: is_edges          !< Flag indicating if this is
                                                                                 !! an edge axis

    CHARACTER(32)                   :: axis_name       !< name of the axis
    CHARACTER(32)                   :: long_name       !< long name of the axis
    CHARACTER(32)                   :: edges_name      !< name of the axis edge
    CHARACTER(128)                  :: units           !< units of the axis
    real(kind=r8_kind), allocatable :: diurnal_data(:) !< Data for the axis
    integer                         :: edges_id        !< Id of the axis edge
    integer                         :: i               !< For do loops

    naxis = naxis + 1

    axis_name = ''
    edges_name = ''
    if (is_edges) then
      WRITE (axis_name,'(a,i2.2)') 'time_of_day_edges_', n_diurnal_samples
      long_name = "time of day edges"
      allocate(diurnal_data(n_diurnal_samples + 1))
      diurnal_data(1) = 0.0
      edges_id = diag_null
      do i = 1, n_diurnal_samples
        diurnal_data(i+1) = 24.0* REAL(i)/n_diurnal_samples
      enddo
    else
      WRITE (axis_name,'(a,i2.2)') 'time_of_day_', n_diurnal_samples
      long_name = "time of day"
      allocate(diurnal_data(n_diurnal_samples))
      edges_id = naxis -1 !< The diurnal edges is the last defined axis
      do i = 1, n_diurnal_samples
        diurnal_data(i) = 24.0*(REAL(i)-0.5)/n_diurnal_samples
      enddo
      WRITE (edges_name,'(a,i2.2)') 'time_of_day_edges_', n_diurnal_samples
    endif

    WRITE (units,11) 'hours', get_base_year(), get_base_month(), &
      get_base_day(), get_base_hour(), get_base_minute(), get_base_second()
11  FORMAT(a,' since ',i4.4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2)

    allocate(fmsDiagDiurnalAxis_type :: diag_axis(naxis)%axis)
    select type (diurnal_axis => diag_axis(naxis)%axis)
    type is (fmsDiagDiurnalAxis_type)
      diurnal_axis%axis_id = naxis
      diurnal_axis%ndiurnal_samples = n_diurnal_samples
      diurnal_axis%axis_name = trim(axis_name)
      diurnal_axis%long_name = trim(long_name)
      diurnal_axis%units = trim(units)
      diurnal_axis%diurnal_data = diurnal_data
      diurnal_axis%edges_id = edges_id
      if (is_edges) &
        WRITE (edges_name,'(a,i2.2)') 'time_of_day_edges_', n_diurnal_samples
        diurnal_axis%edges_name = trim(edges_name)
    end select
  end subroutine define_diurnal_axis

  !< @brief Determine if the axis is in the unstructured grid
  !! @return .True. if the axis is in unstructured grid
  pure logical function is_unstructured_grid(this)
    class(fmsDiagAxis_type),           target, INTENT(in)    :: this        !< diag_axis obj

    is_unstructured_grid = .false.
    select type (this)
    type is (fmsDiagFullAxis_type)
      is_unstructured_grid = trim(this%cart_name) .eq. "U"
    end select
  end function is_unstructured_grid

  !< @brief Adds the structured axis ids to the axis object
  subroutine add_structured_axis_ids(this, axis_ids)
    class(fmsDiagAxis_type),           target, INTENT(inout) :: this        !< diag_axis obj
    integer,                                   intent(in)    :: axis_ids(2) !< axis ids to add to the axis object

    select type (this)
    type is (fmsDiagFullAxis_type)
      allocate(this%structured_ids(2))
      this%structured_ids = axis_ids
    end select
  end subroutine add_structured_axis_ids

  !< @brief Get the structured axis ids from the axis object
  !! @return the structured axis ids
  pure function get_structured_axis(this) &
  result(rslt)
    class(fmsDiagAxis_type),           target, INTENT(in) :: this !< diag_axis obj
    integer :: rslt(2)

    rslt = diag_null
    select type (this)
    type is (fmsDiagFullAxis_type)
      rslt = this%structured_ids
    end select
  end function get_structured_axis


  !< @brief Get the edges_id of an axis_object
  !! @return The edges_id of an axis object
  pure integer function get_edges_id(this)
    class(fmsDiagAxis_type), INTENT(in)    :: this        !< diag_axis obj

    get_edges_id = diag_null
    select type (this)
    type is (fmsDiagFullAxis_type)
      if (allocated(this%edges_id)) get_edges_id = this%edges_id
    end select
  end function

  !> @brief Get the starting and ending indices of the global io domain of the axis
  subroutine get_global_io_domain(this, global_io_index, use_collective_writes)
    class(fmsDiagFullAxis_type), target, intent(in)  :: this               !< diag_axis obj
    integer,                     intent(out) :: global_io_index(2) !< Global io domain starting and ending index
    logical,                     intent(in)  :: use_collective_writes !< .True. if using collective writes

    type(domain2d), pointer :: io_domain !< pointer to the io domain

    global_io_index(1) = 1
    global_io_index(2) = this%length

    if (allocated(this%axis_domain)) then
      select type(domain => this%axis_domain)
      type is (diagDomain2d_t)
        if (use_collective_writes) then
          io_domain => domain%domain2
        else
          io_domain => mpp_get_io_domain(domain%domain2)
        endif

        if (this%cart_name .eq. "X") then
          call mpp_get_global_domain(io_domain, xbegin=global_io_index(1), xend=global_io_index(2), &
            position=this%domain_position)
        elseif (this%cart_name .eq. "Y") then
          call mpp_get_global_domain(io_domain, ybegin=global_io_index(1), yend=global_io_index(2), &
            position=this%domain_position)
        endif
      end select
    endif
  end subroutine get_global_io_domain

  !> @brief Get the length of the axis
  !> @return axis length
  function get_axis_length(this) &
  result (axis_length)
    class(fmsDiagFullAxis_type), intent(in) :: this !< diag_axis obj
    integer                                 :: axis_length

    !< If the axis is domain decomposed axis_length will be set to the length for the current PE:
    if (allocated(this%axis_domain)) then
      axis_length = this%axis_domain%length(this%cart_name, this%domain_position, this%length)
    else
      axis_length = this%length
    endif

  end function


  !> @brief Determine if an axis object has an auxiliary name
  !! @return .true. if an axis object has an auxiliary name
  pure function has_aux(this) &
  result(rslt)
    class(fmsDiagFullAxis_type), intent(in)  :: this               !< diag_axis obj
    logical :: rslt

    rslt = .false.
    if (allocated(this%aux)) rslt = trim(this%aux) .ne. ""
  end function has_aux

  !> @brief Determine if an axis object has a set_name
  !! @return .true. if an axis object has a set_name
  pure function has_set_name(this) &
  result(rslt)
    class(fmsDiagFullAxis_type), intent(in)  :: this               !< diag_axis obj
    logical :: rslt

    rslt = .false.
    if (allocated(this%set_name)) rslt = trim(this%set_name) .ne. ""
  end function has_set_name

  !> @brief Determine if an axis object is an x or y axis
  !! @return .true. if an axis object is an x or y axis, optionally return a flag indicating which it is
  function is_x_or_y_axis(this, x_or_y) &
  result(rslt)
    class(fmsDiagFullAxis_type), intent(in)    :: this   !< diag_axis obj
    integer, optional,           intent(inout) :: x_or_y !< returns is_x_axis if it is a x axis
                                                         !! is_y_axis if it is a y axis
    logical :: rslt

    select case (trim(this%cart_name))
    case ("X")
      if (present(x_or_y)) x_or_y = is_x_axis
      rslt = .true.
    case ("Y")
      if (present(x_or_y)) x_or_y = is_y_axis
      rslt = .true.
    case default
      rslt = .false.
      if (present(x_or_y)) x_or_y = diag_null
    end select
  end function is_x_or_y_axis

  !< @brief Get the global size of the axis, and the layout
  !! It is assumed that this function is only called on "X" and "Y" axes
  !! using the `is_x_or_y_axis` function from above
  subroutine get_dim_size_layout(this, dim_size, layout)
    class(fmsDiagFullAxis_type), intent(in)    :: this     !< diag_axis obj
    integer,                     intent(out)   :: dim_size !< Size of the dimension
    integer,                     intent(out)   :: layout   !< Layout of the dimension

    integer :: nx, ny
    integer :: layout_xy(2)

    select type (domain => this%axis_domain)
    type is (diagDomain2d_t)
      call mpp_get_global_domain(domain%Domain2, xsize=nx, ysize=ny)
      call mpp_get_layout(domain%Domain2, layout_xy)

      if (this%cart_name .eq. "X") then
        dim_size = nx
        layout = layout_xy(1)
      else if (this%cart_name .eq. "Y") then
        dim_size = ny
        layout = layout_xy(2)
      endif
    end select
  end subroutine get_dim_size_layout

  !> @brief Get the set name of an axis object
  !! @return the set name of an axis object
  pure function get_set_name(this) &
  result(rslt)
    class(fmsDiagFullAxis_type), intent(in)  :: this               !< diag_axis obj
    character(len=:), allocatable :: rslt

    rslt = this%set_name
  end function get_set_name

  !> @brief Get the auxiliary name of an axis object
  !! @return the auxiliary name of an axis object
  pure function get_aux(this) &
  result(rslt)
    class(fmsDiagFullAxis_type), intent(in)  :: this               !< diag_axis obj
    character(len=:), allocatable :: rslt

    rslt = this%aux
  end function get_aux

  !> @brief Set the axis_id
  subroutine set_axis_id(this, axis_id)
    class(fmsDiagFullAxis_type), intent(inout) :: this    !< diag_axis obj
    integer,                     intent(in)    :: axis_id !< Axis_id

    this%axis_id = axis_id

  end subroutine set_axis_id

    !> @brief Set the name and ids of the edges
  subroutine set_edges(this, edges_name, edges_id)
    class(fmsDiagFullAxis_type), intent(inout) :: this       !< diag_axis obj
    CHARACTER(len=*),            intent(in)    :: edges_name !< Name of the edges
    integer,                     intent(in)    :: edges_id   !< Axis id of the edges

    !< Saving the name and the id of the edges axis because it will make it easier to use
    !! downstream (i.e you need the edges name to write the attribute to the current axis,
    !! and you need the edges id to add to the diag file object so that you can write the edges
    !! to the file)
    this%edges_name = edges_name
    this%edges_id = edges_id
  end subroutine set_edges

  !> @brief Determine if the subRegion is in the current PE.
  !! If it is, determine the starting and ending indices of the current PE that belong to the subRegion
  subroutine get_indices(this, compute_idx, corners_indices, starting_index, ending_index, need_to_define_axis)
    class(fmsDiagFullAxis_type), intent(in)    :: this                !< diag_axis obj
    integer,                     intent(in)    :: compute_idx(:)      !< Current PE's compute domain
    class(*),                    intent(in)    :: corners_indices(:)  !< The indices of the corners of the subRegion
    integer,                     intent(out)   :: starting_index      !< Starting index of the subRegion
                                                                      !! for the current PE
    integer,                     intent(out)   :: ending_index        !< Ending index of the subRegion
                                                                      !! for the current PE
    logical,                     intent(out)   :: need_to_define_axis !< .true. if it is needed to define
                                                                      !! an axis

    integer :: subregion_start !< Starting index of the subRegion
    integer :: subregion_end   !< Ending index of the subRegion

    !< Get the rectangular coordinates of the subRegion
    !! If the subRegion is not rectangular, the points outside of the subRegion will be masked
    !! out later
    select type (corners_indices)
    type is (integer(kind=i4_kind))
      subregion_start = minval(corners_indices)
      subregion_end = maxval(corners_indices)
    end select

    !< Initiliaze the output
    need_to_define_axis = .false.
    starting_index = diag_null
    ending_index = diag_null

    !< If the compute domain of the current PE is outisde of the range of sub_axis, return
    if (compute_idx(1) < subregion_start .and. compute_idx(2) < subregion_start) return
    if (compute_idx(1) > subregion_end   .and. compute_idx(2) > subregion_end) return

    need_to_define_axis = .true.
    if (compute_idx(1)  >= subregion_start .and. compute_idx(2) >= subregion_end) then
      !< In this case all the point of the current PE are inside the range of the sub_axis
      starting_index = compute_idx(1)
      ending_index   = subregion_end
    else if (compute_idx(1)  >= subregion_start .and. compute_idx(2) <= subregion_end) then
      !< In this case all the points of the current PE are valid up to the end point
      starting_index = compute_idx(1)
      ending_index   = compute_idx(2)
    else if (compute_idx(1)  <= subregion_start .and. compute_idx(2) <= subregion_end) then
      !< In this case all the points of the current PE are valid starting with t subregion_start
      starting_index = subregion_start
      ending_index   = compute_idx(2)
    else if (compute_idx(1) <= subregion_start .and. compute_idx(2) >= subregion_end) then
      !< In this case only the points in the current PE ar valid
      starting_index = subregion_start
      ending_index   = subregion_end
    endif

    if (this%domain_position .ne. CENTER) then
      if (subregion_end - subregion_start + 1 .eq. 1) then
          !< If your subregion consitsts of just 1 one, only include 1 PE
          if (ending_index .eq. compute_idx(2)) need_to_define_axis = .false.
      else
        if (ending_index - starting_index + 1 .eq. 1) then
          !< If the PEs section is only 1, only include 1 PE
          if (starting_index .eq. compute_idx(2) .or. ending_index .eq. compute_idx(1)) &
            need_to_define_axis = .false.
        endif
      endif
    endif

  end subroutine get_indices

  !< Get the compute domain of the axis
  subroutine get_compute_domain(this, compute_idx, need_to_define_axis, tile_number)
    class(fmsDiagFullAxis_type), intent(in)    :: this                !< diag_axis obj
    integer,                     intent(inout) :: compute_idx(:)      !< Compute domain of the axis
    logical,                     intent(out)   :: need_to_define_axis !< .true. if it needed to define the axis
    integer, optional,           intent(in)    :: tile_number         !< The tile number of the axis

    !< Initialize the output
    need_to_define_axis = .false.
    compute_idx = diag_null

    if (.not. allocated(this%axis_domain)) then
       !< If the axis is not domain decomposed, use the whole axis as the compute domain
       if (this%cart_name .eq. "X" .or. this%cart_name .eq. "Y") then
         compute_idx(1) = 1
         compute_idx(2) = size(this%axis_data)
         need_to_define_axis = .true.
       endif
      return
    endif

    select type(domain => this%axis_domain)
    type is (diagDomain2d_t)
      if (present(tile_number)) then
        !< If the tile number is present and the current PE is not on the tile, then there is no need
        !! to define the axis
        if (any(mpp_get_tile_id(domain%Domain2) .ne. tile_number)) then
          need_to_define_axis = .false.
          return
        endif
      endif

      !< Get the compute domain for the current PE if it is an "X" or "Y" axis
      select case (this%cart_name)
      case ("X")
        call mpp_get_compute_domain(domain%Domain2, xbegin=compute_idx(1), xend=compute_idx(2), &
                 & position=this%domain_position)
        need_to_define_axis = .true.
      case ("Y")
        call mpp_get_compute_domain(domain%Domain2, ybegin=compute_idx(1), yend=compute_idx(2), &
                 & position=this%domain_position)
        need_to_define_axis = .true.
      end select
    end select

  end subroutine get_compute_domain

  !!!!!!!!!!!!!!!!! SUBAXIS PROCEDURES !!!!!!!!!!!!!!!!!

  !> @brief Initialize a subaxis object with region boundaries and metadata.
  !!
  !! This subroutine populates a subaxis object with all necessary information about
  !! the subregion it represents. The subaxis name is automatically generated from the
  !! parent axis name and a sequential number (e.g., "temperature_sub01").
  !!
  !! Optional parameters allow storage of global domain indices (needed for the
  !! domain_decomposition NetCDF attribute) and Z-axis bounds (for identifying
  !! equivalent subaxes across output files).
  subroutine fill_subaxis(this, starting_index, ending_index, axis_id, parent_id, parent_axis_name, compute_idx, &
                          global_idx, zbounds, nz_subaxis)
    class(fmsDiagSubAxis_type)  , INTENT(INOUT) :: this             !< diag_sub_axis obj
    integer                     , intent(in)    :: starting_index   !< Starting index of the subRegion for the PE
    integer                     , intent(in)    :: ending_index     !< Ending index of the subRegion for the PE
    integer                     , intent(in)    :: axis_id          !< Axis id to assign to the subaxis
    integer                     , intent(in)    :: parent_id        !< The id of the parent axis the subaxis belongs to
    character(len=*)            , intent(in)    :: parent_axis_name !< Name of the parent_axis
    integer                     , intent(in)    :: compute_idx(2)   !< Starting and ending index of
                                                                    !! the axis's compute domain
    integer,            optional, intent(in)    :: global_idx(2)   !< Starting and ending index of
                                                                    !! the axis's compute domain
    real(kind=r4_kind), optional, intent(in)    :: zbounds(2)       !< Bounds of the z-axis
    integer,            optional, intent(in)    :: nz_subaxis       !< The number of z subaxis that have been defined
                                                                    !! in the file

    integer :: nsubaxis !< The subaxis number in the axis name subXX
    character(len=2) :: nsubaxis_char !< nsubaxis converted to a string

    nsubaxis = 1
    if (present(nz_subaxis)) nsubaxis = nz_subaxis

    this%axis_id = axis_id
    this%starting_index = starting_index
    this%ending_index = ending_index
    this%parent_axis_id = parent_id
    write(nsubaxis_char, '(i2.2)')  nsubaxis
    this%subaxis_name = trim(parent_axis_name)//"_sub"//nsubaxis_char
    this%compute_idx = compute_idx

    if (present(zbounds)) then
      ! This is needed to avoid duplicating z sub axis!
      allocate(this%zbounds(2))
      this%zbounds = zbounds
    endif

    if (present(global_idx)) then
      ! This is needed for the "domain_decomposition" attribute which is needed for the combiner
      allocate(this%global_idx(2))
      this%global_idx = global_idx
    endif
  end subroutine fill_subaxis

  !> @brief Get the number of points in a subaxis.
  !!
  !! Calculates the number of coordinate points in this subaxis based on its starting
  !! and ending indices.
  function axis_length(this) &
    result(res)
      class(fmsDiagSubAxis_type)  , INTENT(IN) :: this             !< diag_sub_axis obj
      integer :: res

      res = this%ending_index - this%starting_index + 1
    end function

   !> @brief Get the starting index of this subaxis.
  !!
  !! Returns the first index of the subregion on the current process.
  function get_starting_index(this) result(indx)
    class(fmsDiagSubAxis_type), intent(in) :: this !< diag_sub_axis object
    integer :: indx !< Result to return
    indx = this%starting_index
  end function get_starting_index

  !> @brief Accesses its member ending_index
  !! @return a copy of the ending_index
  function get_ending_index(this) result(indx)
    class(fmsDiagSubAxis_type), intent(in) :: this !< diag_sub_axis object
    integer :: indx !< Result to return
    indx = this%ending_index
  end function get_ending_index

  !> @brief Accesses its member compute_indices
  !! @return a copy of the ending_index
  function get_compute_indices(this) result(indx)
    class(fmsDiagSubAxis_type), intent(in) :: this !< diag_sub_axis object
    integer :: indx(2) !< Result to return
    indx = this%compute_idx
  end function get_compute_indices

  !> @brief Determines if the zbounds passed in are the same as those in the file
  !! @return .True. if the zbounds are the same
  function is_same_zbounds(this, zbounds) result(is_same)
    class(fmsDiagSubAxis_type), intent(in) :: this       !< diag_sub_axis object
    real(kind=r4_kind),         intent(in) :: zbounds(2) !< zbounds to compare with
    logical :: is_same

    is_same = zbounds(1) .eq. this%zbounds(1) .and. zbounds(2) .eq. this%zbounds(2)
  end function

  !> @brief Get the ntiles in a domain
  !> @return the number of tiles in a domain
  function get_ntiles(this) &
  result (ntiles)
    class(diagDomain_t), INTENT(IN)    :: this       !< diag_axis obj

    integer :: ntiles

    select type (this)
    type is (diagDomain2d_t)
      ntiles = mpp_get_ntile_count(this%domain2)
    end select
  end function get_ntiles

  !> @brief Get the length of a 2D domain
  !> @return Length of the 2D domain
  function get_length(this, cart_axis, domain_position, global_length) &
  result (length)
    class(diagDomain_t), INTENT(IN)    :: this       !< diag_axis obj
    character(len=*),    INTENT(IN)    :: cart_axis !< cart_axis of the axis, must be "X" or "Y"
    integer,             INTENT(IN)    :: domain_position !< Domain position (CENTER, NORTH, EAST)
    integer,             INTENT(IN)    :: global_length !< global_length of the axis

    integer :: length

    select type (this)
    type is(diagDomain2d_t)
      if (trim(cart_axis) == "X") call mpp_get_compute_domain(this%Domain2, xsize=length, position=domain_position)
      if (trim(cart_axis) == "Y") call mpp_get_compute_domain(this%Domain2, ysize=length, position=domain_position)
    class default
      !< If domain is 1D or UG, just set it to the global length
      length = global_length
    end select
  end function get_length

  !!!!!!!!!!!!!!!!! FMS_DOMAIN PROCEDURES !!!!!!!!!!!!!!!!!

  !> @brief Assign a specific domain object to this domain wrapper.
  !!
  !! This subroutine assigns the appropriate domain object (1D, 2D, or unstructured)
  !! to the polymorphic domain wrapper based on the type-specific argument provided.
  subroutine set_axis_domain(this, Domain, Domain2, DomainU)
    class(diagDomain_t) :: this !< fms_domain obj
    TYPE(domain1d),     INTENT(in),  OPTIONAL :: Domain  !< 1d domain
    TYPE(domain2d),     INTENT(in),  OPTIONAL :: Domain2 !< 2d domain
    TYPE(domainUG),     INTENT(in),  OPTIONAL :: DomainU !< Unstructured domain

    select type(this)
    type is (diagDomain1d_t)
      this%Domain = Domain
    type is (diagDomain2d_t)
      this%Domain2 = Domain2
    type is (diagDomainUg_t)
      this%DomainUG = DomainU
    end select
  end subroutine set_axis_domain

  !> @brief Allocate the module-level array of diagnostic axis containers.
  !!
  !! This initialization routine must be called before any axes are created.
  !! It allocates a global array to hold all axis objects up to the limit
  !! defined by max_axes. Returns .true. on successful allocation.
  !!
  !! @return .true. if allocation succeeded
  logical function fms_diag_axis_object_init(axis_array)
    class(fmsDiagAxisContainer_type)   , allocatable, intent(inout) :: axis_array(:) !< Array of diag_axis

    if (allocated(axis_array)) call mpp_error(FATAL, "The diag_axis containers is already allocated")
    allocate(axis_array(max_axes))
    !axis_array%axis_id = DIAG_NULL

    fms_diag_axis_object_init = .true.
  end function fms_diag_axis_object_init

  !> @brief Deallocate the module-level array of diagnostic axis containers.
  !!
  !! This cleanup routine should be called when axis management is no longer needed.
  !! It deallocates the global axis array, freeing all associated memory.
  !!
  !! @return .false. after successful deallocation
  logical function fms_diag_axis_object_end(axis_array)
    class(fmsDiagAxisContainer_type)   , allocatable, intent(inout) :: axis_array(:) !< Array of diag_axis

    if (allocated(axis_array)) deallocate(axis_array)
    fms_diag_axis_object_end = .false.

  end function fms_diag_axis_object_end

  !> @brief Determine the name of the axis (or subaxis) represented by the given axis object.
  !!
  !! Returns the axis name, optionally modified for regional subsets.
  !! When is_regional is .true. and this is an X or Y axis, "_sub01" is appended
  !! to indicate a regional subset.
  !!
  !! @note This function may be called from field objects to get dimension names.
  !! Field objects know only parent axis IDs; the is_regional flag distinguishes
  !! whether regional subaxis naming should be applied.
  !!
  !! @return The axis name (possibly modified for regional output)
  pure function get_axis_name(this, is_regional) &
  result(axis_name)
    class(fmsDiagAxis_type), intent(in)           :: this        !< Axis object
    logical,                 intent(in), optional :: is_regional !< Flag indicating if the axis is regional

    character(len=:),       allocatable :: axis_name

    select type (this)
    type is (fmsDiagFullAxis_type)
      axis_name = this%axis_name
      if (present(is_regional)) then
        if (is_regional) then
          if (this%cart_name .eq. "X" .or. this%cart_name .eq. "Y") axis_name = axis_name//"_sub01"
        endif
      endif
    type is (fmsDiagSubAxis_type)
      axis_name = this%subaxis_name
    end select
  end function get_axis_name

  !< @brief Determine if the axis is a Z axis by looking at the cartesian name
  !! @return .True. if the axis is a Z axis
  pure logical function is_z_axis(this)
    class(fmsDiagAxis_type), intent(in)           :: this        !< Axis object
    is_z_axis = .false.
    select type (this)
    type is (fmsDiagFullAxis_type)
      if (this%cart_name .eq. "Z") is_z_axis = .true.
    end select
  end function

  !> @brief Validate a Cartesian axis classification name.
  !!
  !! Checks that the provided Cartesian axis name is one of the recognized values.
  !! Valid values are: X (longitude), Y (latitude), Z (vertical), T (time),
  !! U (unstructured), or N (no axis). Calls mpp_error(FATAL) if invalid.
  subroutine check_if_valid_cart_name(cart_name)
    character(len=*), intent(in) :: cart_name !< axis cartesian classification name to validate

    select case (cart_name)
    case ("X", "Y", "Z", "T", "U", "N")
    case default
      call mpp_error(FATAL, "diag_axit_init: Invalid cart_name: "//cart_name//&
                             "The acceptable values are X, Y, Z, T, U, N.")
    end select
  end subroutine check_if_valid_cart_name

  !> @brief Validate a domain position classification.
  !!
  !! Checks that the provided domain position is one of the recognized values.
  !! Valid values are: CENTER, NORTH, or EAST (from mpp_domains_mod).
  !! Calls mpp_error(FATAL) if invalid.
  subroutine check_if_valid_domain_position(domain_position)
    integer, INTENT(IN) :: domain_position !< domain position to validate (CENTER, NORTH, EAST)

    select case (domain_position)
    case (CENTER, NORTH, EAST)
    case default
        call mpp_error(FATAL, "diag_axit_init: Invalid domain_positon. &
                              &The acceptable values are NORTH, EAST, CENTER")
    end select
  end subroutine check_if_valid_domain_position

  !> @brief Validate an axis direction indicator.
  !!
  !! Checks that the provided direction value is one of the recognized values.
  !! Valid values are: -1 (down/decreasing), 0 (no direction), or 1 (up/increasing).
  !! Calls mpp_error(FATAL) if invalid.
  subroutine check_if_valid_direction(direction)
    integer, INTENT(IN) :: direction !< Direction indicator to validate (-1, 0, 1)

    select case(direction)
    case(-1, 0, 1)
    case default
      call mpp_error(FATAL, "diag_axit_init: Invalid direction. &
                            &The acceptable values are-1 0 1")
    end select
  end subroutine check_if_valid_direction

  !> @brief Determine the domain type and domain object from a list of axis IDs.
  !!
  !! This utility subroutine examines all axes used by a variable and determines
  !! which domain type (and domain object) should be used for output. It handles
  !! the case where a variable has both domain-decomposed and non-decomposed axes
  !! (e.g., X,Y and Z axes on different domain types).
  subroutine get_domain_and_domain_type(diag_axis, axis_id, domain_type, domain, var_name)
    class(fmsDiagAxisContainer_type), target, intent(in)  :: diag_axis(:)  !< Array of diag_axis
    integer,                      INTENT(IN)  :: axis_id(:)    !< Array of axis ids
    integer,                      INTENT(OUT) :: domain_type   !< fileobj_type to use
    CLASS(diagDomain_t), POINTER, INTENT(OUT) :: domain        !< Domain
    character(len=*),             INTENT(IN)  :: var_name      !< Name of the variable (for error messages)

    integer :: i !< For do loops
    integer :: j !< axis_id(i) (for less typing)

    domain_type = NO_DOMAIN
    domain => null()

    do i = 1, size(axis_id)
      j = axis_id(i)
      select type (axis => diag_axis(j)%axis)
      type is (fmsdiagfullaxis_type)
        !< Check that all the axis are in the same domain
        if (domain_type .ne. axis%type_of_domain) then
          !< If they are different domains, one of them can be NO_DOMAIN
          !! i.e a variable can have axis that are domain decomposed (x,y) and an axis that isn't (z)
          if (domain_type .eq. NO_DOMAIN .or. axis%type_of_domain .eq. NO_DOMAIN ) then
            !< Update the domain_type and domain, if needed
            if ((axis%type_of_domain .eq. TWO_D_DOMAIN  .and. size(axis_id) > 1) &
               & .or. axis%type_of_domain .eq. UG_DOMAIN) then
                domain_type = axis%type_of_domain
                domain => axis%axis_domain
            endif
          else
            call mpp_error(FATAL, "The variable:"//trim(var_name)//" has axis that are not in the same domain")
          endif
        endif
      end select
    enddo
  end subroutine get_domain_and_domain_type

  !> @brief Create a subaxis for a subregion defined by explicit index bounds.
  !!
  !! This subroutine creates a new subaxis for output of a limited region specified
  !! by starting and ending indices (corners). It determines which processes
  !! contribute data to the subregion and allocates the subaxis accordingly.
  !!
  !! @note Only processes with compute domains overlapping the subregion will have write_on_this_pe = .true.
  subroutine define_new_subaxis_index(parent_axis, subRegion, diag_axis, naxis, is_x_or_y, write_on_this_pe)
    class(fmsDiagAxisContainer_type), target, intent(inout) :: diag_axis(:)     !< Diag_axis object
    type(fmsDiagFullAxis_type),               intent(inout) :: parent_axis      !< axis object of the parent
    integer,                                  intent(inout) :: naxis            !< Number of axis registered
    type(subRegion_type),                     intent(in)    :: subRegion        !< SubRegion definition from the yaml
    integer,                                  intent(in)    :: is_x_or_y        !< Flag indicating if it is
                                                                                !! a x or y axis
    logical,                                  intent(out)   :: write_on_this_pe !< .true. if the subregion
                                                                                !! is on this PE
    integer :: compute_idx(2) !< Indices of the compute domain
    integer :: global_idx(2)  !< Indices of the "global" domain
    integer :: starting_index !< starting index of the subregion
    integer :: ending_index   !< ending index of the subregion

    call parent_axis%get_compute_domain(compute_idx, write_on_this_pe, tile_number=subRegion%tile)
    if (.not. write_on_this_pe) return

    !< Determine if the PE's compute domain is inside the subRegion
    !! If it is get the starting and ending indices for that PE
    call parent_axis%get_indices(compute_idx, subRegion%corners(:,is_x_or_y), starting_index, ending_index, &
      write_on_this_pe)

    if (.not. write_on_this_pe) return

    select type(corners=> subRegion%corners)
    type is (integer(kind=i4_kind))
      global_idx(1) = minval(corners(:,is_x_or_y))
      global_idx(2) = maxval(corners(:,is_x_or_y))
    end select

    !< If it made it to this point, the current PE is in the subRegion!
    call define_new_axis(diag_axis, parent_axis, naxis, parent_axis%axis_id, &
          starting_index, ending_index, compute_idx, global_idx)

  end subroutine define_new_subaxis_index

  !> @brief Create subaxes for a subregion defined by latitude/longitude bounds.
  !!
  !! This subroutine creates new X and Y subaxes for output of a geographic region
  !! specified by lat/lon bounds. It handles both cubesphere and regular lat/lon grids,
  !! determining the corresponding index ranges for each grid type.
  !!
  !! @note Uses cubesphere index functions if is_cube_sphere is .true., otherwise nearest_index
  subroutine define_new_subaxis_latlon(diag_axis, axis_ids, naxis, subRegion, is_cube_sphere, write_on_this_pe)
    class(fmsDiagAxisContainer_type), target, intent(inout) :: diag_axis(:)     !< Diag_axis object
    integer,                                  INTENT(in)    :: axis_ids(:)      !< Array of axes_ids
    integer,                                  intent(inout) :: naxis            !< Number of axis registered
    type(subRegion_type),                     intent(in)    :: subRegion        !< SubRegion definition from the yaml
    logical,                                  intent(in)    :: is_cube_sphere   !< .true. if this is a cubesphere
    logical,                                  intent(out)   :: write_on_this_pe !< .true. if the subregion
                                                                                !! is on this PE

    real    :: lat(2)                 !< Starting and ending lattiude of the subRegion
    real    :: lon(2)                 !< Starting and ending longitude or the subRegion
    integer :: lat_indices(2)         !< Starting and ending latitude indices of the subRegion
    integer :: lon_indices(2)         !< Starting and ending longitude indices of the subRegion
    integer :: compute_idx(2)         !< Compute domain of the current axis
    integer :: starting_index(2)      !< Starting index of the subRegion for the current PE for the "x" and "y"
                                      !! direction
    integer :: ending_index(2)        !< Ending index of the subRegion for the current PE for the "x" and "y" direction
    logical :: need_to_define_axis(2) !< .true. if it is needed to define the subaxis for the "x" and "y" direction
    integer :: i                      !< For do loops
    integer :: parent_axis_ids(2)     !< The axis id of the parent axis for the "x" and "y" direction
    logical :: is_x_y_axis            !< .true. if the axis is x or y
    integer :: compute_idx_2(2, 2)    !< Starting and ending indices of the compute domain for the "x" and "y" direction
    integer :: global_idx (2, 2)      !< Starting and ending indices of the global domain for the "x" and "y" direction

    write_on_this_pe = .false.
    need_to_define_axis = .true.
    parent_axis_ids = diag_null

    !< Get the rectangular coordinates of the subRegion
    !! If the subRegion is not rectangular, the points outside of the subRegion will be masked
    !! out later
    select type (corners => subRegion%corners)
    type is (real(kind=r4_kind))
      lon(1) = minval(corners(:,1))
      lon(2) = maxval(corners(:,1))
      lat(1) = minval(corners(:,2))
      lat(2) = maxval(corners(:,2))
    end select

    if_is_cube_sphere: if (is_cube_sphere) then
      !< Get the starting and ending indices of the subregion in the cubesphere relative to the global domain
      call get_local_indices_cubesphere(lat(1), lat(2), lon(1), lon(2),&
          & lon_indices(1), lon_indices(2), lat_indices(1), lat_indices(2))
      loop_over_axis_ids: do i = 1, size(axis_ids)
        select_axis_type: select type (parent_axis => diag_axis(axis_ids(i))%axis)
        type is (fmsDiagFullAxis_type)
          !< Get the PEs compute domain
          call parent_axis%get_compute_domain(compute_idx, is_x_y_axis)

          !< If this is not a "X" or "Y" axis go to the next axis
          if (.not. is_x_y_axis) cycle

          !< Determine if the PE's compute domain is inside the subRegion
          !! If it is get the starting and ending indices for that PE
          if (parent_axis%cart_name .eq. "X") then
            call parent_axis%get_indices(compute_idx, lon_indices, starting_index(1), ending_index(1), &
              need_to_define_axis(1))
            parent_axis_ids(1) = axis_ids(i)
            compute_idx_2(1,:) = compute_idx
            global_idx(1,:) = lon_indices
          else if (parent_axis%cart_name .eq. "Y") then
            call parent_axis%get_indices(compute_idx, lat_indices, starting_index(2), ending_index(2), &
              need_to_define_axis(2))
            parent_axis_ids(2) = axis_ids(i)
            compute_idx_2(2,:) = compute_idx
            global_idx(2,:) = lat_indices
          endif
        end select select_axis_type
      enddo loop_over_axis_ids
    else if_is_cube_sphere
      loop_over_axis_ids2: do i = 1, size(axis_ids)
        select type (parent_axis => diag_axis(axis_ids(i))%axis)
        type is (fmsDiagFullAxis_type)
          !< Get the PEs compute domain
          call parent_axis%get_compute_domain(compute_idx, is_x_y_axis)

          !< If this is not a "X" or "Y" axis go to the next axis
          if (.not. is_x_y_axis) cycle

          !< Get the starting and ending indices of the subregion relative to the global grid
          if (parent_axis%cart_name .eq. "X") then
            select type(adata=>parent_axis%axis_data)
            type is (real(kind=r8_kind))
              lon_indices(1) = nearest_index(real(lon(1), kind=r8_kind), adata)
              lon_indices(2) = nearest_index(real(lon(2), kind=r8_kind), adata)
            type is (real(kind=r4_kind))
              lon_indices(1) = nearest_index(real(lon(1), kind=r4_kind), adata)
              lon_indices(2) = nearest_index(real(lon(2), kind=r4_kind), adata)
            end select
            call parent_axis%get_indices(compute_idx, lon_indices, starting_index(1), ending_index(1), &
              need_to_define_axis(1))
            parent_axis_ids(1) = axis_ids(i)
            compute_idx_2(1,:) = compute_idx
            global_idx(1,:) = lon_indices
          else if (parent_axis%cart_name .eq. "Y") then
            select type(adata=>parent_axis%axis_data)
            type is (real(kind=r8_kind))
              lat_indices(1) = nearest_index(real(lat(1), kind=r8_kind), adata)
              lat_indices(2) = nearest_index(real(lat(2), kind=r8_kind), adata)
            type is (real(kind=r4_kind))
              lat_indices(1) = nearest_index(real(lat(1), kind=r4_kind), adata)
              lat_indices(2) = nearest_index(real(lat(2), kind=r4_kind), adata)
            end select
            call parent_axis%get_indices(compute_idx, lat_indices, starting_index(2), ending_index(2), &
              need_to_define_axis(2))
            parent_axis_ids(2) = axis_ids(i)
            compute_idx_2(2,:) = compute_idx
            global_idx(2,:) = lat_indices
          endif
        end select
      enddo loop_over_axis_ids2
    endif if_is_cube_sphere

    !< If the PE's compute is not inside the subRegion move to the next axis
    if (any(.not. need_to_define_axis )) return

    !< If it made it to this point, the current PE is in the subRegion!
    write_on_this_pe = .true.

    do i = 1, size(parent_axis_ids)
      if (parent_axis_ids(i) .eq. diag_null) cycle
      select type (parent_axis => diag_axis(parent_axis_ids(i))%axis)
      type is (fmsDiagFullAxis_type)
        call define_new_axis(diag_axis, parent_axis, naxis, parent_axis_ids(i), &
          starting_index(i), ending_index(i), compute_idx_2(i,:), global_idx(i,:))
     end select
    enddo

  end subroutine define_new_subaxis_latlon

  !> @brief Create a new subaxis and initialize it with index ranges.
  !!
  !! This is the core subroutine for creating subaxes. It allocates a new fmsDiagSubAxis_type
  !! object, fills it with the provided index information, and registers it with the parent axis.
  !! The parent axis's subaxis list is updated to include the new subaxis ID.
  subroutine define_new_axis(diag_axis, parent_axis, naxis, parent_id, &
                             starting_index, ending_index, compute_idx, global_idx, new_axis_id, zbounds, &
                             nz_subaxis)

    class(fmsDiagAxisContainer_type), target, intent(inout) :: diag_axis(:)     !< Diag_axis object
    class(fmsDiagFullAxis_type),              intent(inout) :: parent_axis      !< The parent axis
    integer,                                  intent(inout) :: naxis            !< The number of axis that
                                                                                !! have been defined
    integer,                                  intent(in)    :: parent_id        !< Id of the parent axis
    integer,                                  intent(in)    :: starting_index   !< PE's Starting index
    integer,                                  intent(in)    :: ending_index     !< PE's Ending index
    integer,                                  intent(in)    :: compute_idx(2)   !< Starting and ending index of
                                                                                !! the axis's compute domain
    integer,                        optional, intent(in)    :: global_idx(2)    !< Starting and ending index of
                                                                                !! the axis's global domain
    integer,                        optional, intent(out)   :: new_axis_id      !< Axis id of the axis this is creating
    real(kind=r4_kind),             optional, intent(in)    :: zbounds(2)       !< Bounds of the Z axis
    integer,                        optional, intent(in)    :: nz_subaxis       !< The number of z subaxis that have
                                                                                !! been defined in the file

    naxis = naxis + 1 !< This is the axis id of the new axis!

    !< Add the axis_id of the new subaxis to the parent axis
    parent_axis%nsubaxis = parent_axis%nsubaxis + 1
    parent_axis%subaxis(parent_axis%nsubaxis) = naxis

    !< Allocate the new axis as a subaxis and fill it
    allocate(fmsDiagSubAxis_type :: diag_axis(naxis)%axis)
    diag_axis(naxis)%axis%axis_id = naxis
    if (present(new_axis_id)) new_axis_id = naxis

    select type (sub_axis => diag_axis(naxis)%axis)
    type is (fmsDiagSubAxis_type)
      call sub_axis%fill_subaxis(starting_index, ending_index, naxis, parent_id, &
             parent_axis%axis_name, compute_idx, global_idx=global_idx, zbounds=zbounds, nz_subaxis=nz_subaxis)
    end select
  end subroutine define_new_axis

  !> @brief Get the ID of the parent axis if this is a subaxis.
  !!
  !! For subaxis objects, returns the axis ID of the parent full axis.
  !! For non-subaxis objects (full axes), returns diag_null.
  !!
  !! @return Parent axis ID if this is a subaxis; diag_null otherwise
  pure function get_parent_axis_id(this) &
  result(parent_axis_id)

    class(fmsDiagAxis_type), intent(in) :: this            !< Axis Object
    integer                             :: parent_axis_id

    select type (this)
    type is (fmsDiagFullAxis_type)
      parent_axis_id = diag_null
    type is (fmsDiagSubAxis_type)
      parent_axis_id = this%parent_axis_id
    type is (fmsDiagDiurnalAxis_type)
      parent_axis_id = diag_null
    end select

  end function

  !> @brief Get the ID of the most recently defined subaxis of this axis.
  !!
  !! For non-Z axes with subaxes defined, returns the ID of the latest subaxis.
  !! For Z axes or axes without subaxes, returns the axis's own ID.
  !!
  !! @return ID of the most recent subaxis, or this%axis_id if no subaxes
  pure function get_subaxes_id(this) &
  result(sub_axis_id)

    class(fmsDiagAxis_type), intent(in) :: this !< Axis Object
    integer :: sub_axis_id

    sub_axis_id = this%axis_id
    select type (this)
    type is (fmsDiagFullAxis_type)
      if (this%cart_name .ne. "Z") sub_axis_id = this%subaxis(this%nsubaxis)
    end select

  end function

  !> @brief Parse the "compress" dimension specification into structured axis names.
  !!
  !! The "compress" attribute defines a compressed (packed) dimension by listing
  !! the names of two or more axes that are combined into a single dimension.
  !! This function extracts those axis names from the attribute value.
  !!
  !! @return Array of axis names
  pure function parse_compress_att(compress_att) &
  result(axis_names)
    class(*), intent(in) :: compress_att(:) !< The compress attribute to parse
    character(len=120)   :: axis_names(2)

    integer            :: ios           !< Errorcode after parsing the compress attribute

    select type (compress_att)
      type is (character(len=*))
        read(compress_att(1),*, iostat=ios) axis_names
        if (ios .ne. 0) axis_names = ""
      class default
        axis_names = ""
    end select
  end function parse_compress_att

  !> @brief Look up an axis ID by name and optional axis set name.
  !!
  !! Searches the registered axes for one matching the given name and set name.
  !! The set_name parameter allows disambiguation when multiple axis configurations
  !! with the same name are defined (e.g., different resolutions on the same grid).
  !!
  !! @return Axis ID if found; diag_null if not found
  pure function get_axis_id_from_name(axis_name, diag_axis, naxis, set_name) &
  result(axis_id)
    class(fmsDiagAxisContainer_type), intent(in) :: diag_axis(:) !< Array of axis object
    character(len=*),                 intent(in) :: axis_name    !< Name of the axis
    integer,                          intent(in) :: naxis        !< Number of axis that have been registered
    character(len=*),                 intent(in) :: set_name     !< Name of the axis set
    integer                                      :: axis_id

    integer :: i !< For do loops

    axis_id = diag_null
    do i = 1, naxis
      select type(axis => diag_axis(i)%axis)
      type is (fmsDiagFullAxis_type)
        if (trim(axis%axis_name) .eq. trim(axis_name)) then
          if (trim(axis%set_name) .eq. trim(set_name)) then
            axis_id = i
            return
          endif
        endif
      end select
    enddo

  end function get_axis_id_from_name

  !> @brief Get the number of time-of-day samples for a diurnal axis.
  pure function get_diurnal_axis_samples(this) &
  result(n_diurnal_samples)

    class(fmsDiagDiurnalAxis_type), intent(in) :: this !< Axis Object
    integer :: n_diurnal_samples

    n_diurnal_samples = this%ndiurnal_samples
  end function get_diurnal_axis_samples

  !> @brief Write diurnal axis metadata to a NetCDF file.
  !!
  !! Registers the diurnal axis dimension and variable in the file object,
  !! along with associated attributes (units, long_name, edges reference).
  subroutine write_diurnal_metadata(this, fms2io_fileobj)
    class(fmsDiagDiurnalAxis_type), intent(in)    :: this     !< Diurnal axis Object
    class(FmsNetcdfFile_t),         intent(inout) :: fms2io_fileobj  !< Fms2_io fileobj to write the data to

    call register_axis(fms2io_fileobj, this%axis_name, size(this%diurnal_data))
    call register_field(fms2io_fileobj, this%axis_name, pack_size_str, (/trim(this%axis_name)/))
    call register_variable_attribute(fms2io_fileobj, this%axis_name, "units", &
                                    &trim(this%units), str_len=len_trim(this%units))
    call register_variable_attribute(fms2io_fileobj, this%axis_name, "long_name", &
                                    &trim(this%long_name), str_len=len_trim(this%long_name))
    if (this%edges_id .ne. diag_null) &
      call register_variable_attribute(fms2io_fileobj, this%axis_name, "edges", &
                                      &trim(this%edges_name), str_len=len_trim(this%edges_name))
  end subroutine write_diurnal_metadata

  !> @brief Create or reuse a Z-axis subaxis for the specified Z bounds.
  !!
  !! This subroutine manages Z-axis compression by creating subaxes for specific
  !! depth ranges. If an identical Z subaxis already exists in the file, that
  !! existing subaxis is reused. Otherwise, a new subaxis is created.
  !!
  !! This deduplication avoids creating multiple subaxes with identical ranges.
  subroutine create_new_z_subaxis(zbounds, var_axis_ids, diag_axis, naxis, file_axis_id, nfile_axis, nz_subaxis, &
                                  error_mseg)
    real(kind=r4_kind),                       intent(in)    :: zbounds(2)      !< Bounds of the Z axis
    integer,                                  intent(inout) :: var_axis_ids(:) !< The variable's axis_ids
    class(fmsDiagAxisContainer_type), target, intent(inout) :: diag_axis(:)    !< Array of diag_axis objects
    integer,                                  intent(inout) :: naxis           !< Number of axis that have been
                                                                               !! registered
    integer,                                  intent(inout) :: file_axis_id(:) !< The file's axis_ids
    integer,                                  intent(inout) :: nfile_axis      !< Number of axis that have been
                                                                               !! defined in file
    integer,                                  intent(inout) :: nz_subaxis      !< The number of z subaxis currently
                                                                               !! defined in the file
    character(len=*),                         intent(inout) :: error_mseg      !! Message to include in error message
                                                                               !! if there is an error

    class(*), pointer :: zaxis_data(:)      !< The data of the full zaxis
    integer           :: subaxis_indices(2) !< The starting and ending indices of the subaxis relative to the full
                                            !! axis
    integer           :: i                  !< For do loops
    integer           :: subaxis_id         !< The id of the new z subaxis
    integer           :: parent_axis_id     !< Id of parent axis id
    integer           :: zaxis_index        !< Index of the z axis (i.e 3 if the variable is x,y,z)
    type(fmsDiagFullAxis_type), pointer :: parent_axis !< Pointer to the parent axis

    parent_axis_id = diag_null
    zaxis_index = diag_null

    !< Determine which axis is the z axis:
    do i = 1, size(var_axis_ids)
      select type (parent_axis => diag_axis(var_axis_ids(i))%axis)
      type is (fmsDiagFullAxis_type)
        if (parent_axis%cart_name .eq. "Z") then
          parent_axis_id = var_axis_ids(i)
          zaxis_index = i
        endif
      end select
    enddo

    if (parent_axis_id .eq. DIAG_NULL) then
        call mpp_error(FATAL, "create_new_z_subaxis:: unable to find the zaxis for "//trim(error_mseg))
    endif

    !< Determine if the axis was already created
    do i = 1, nfile_axis
      select type (axis => diag_axis(file_axis_id(i))%axis)
      type is (fmsDiagSubAxis_type)
        if (axis%parent_axis_id .ne. parent_axis_id) cycle
        if (axis%zbounds(1) .eq. zbounds(1) .and. axis%zbounds(2) .eq. zbounds(2)) then
          var_axis_ids(zaxis_index) = file_axis_id(i)
          return
        endif
      end select
    enddo

    select type (axis => diag_axis(parent_axis_id)%axis)
    type is (fmsDiagFullAxis_type)
      zaxis_data => axis%axis_data
      parent_axis => axis
    end select

    select type(zaxis_data)
      type is (real(kind=r4_kind))
        !TODO need to include the conversion to "real" because nearest_index doesn't take r4s and r8s
        subaxis_indices(1) = nearest_index(real(zbounds(1)), real(zaxis_data))
        subaxis_indices(2) = nearest_index(real(zbounds(2)), real(zaxis_data))
      type is (real(kind=r8_kind))
        subaxis_indices(1) = nearest_index(real(zbounds(1)), real(zaxis_data))
        subaxis_indices(2) = nearest_index(real(zbounds(2)), real(zaxis_data))
    end select

    nz_subaxis = nz_subaxis + 1
    call define_new_axis(diag_axis, parent_axis, naxis, parent_axis%axis_id, &
                        &subaxis_indices(1), subaxis_indices(2), (/lbound(zaxis_data,1), ubound(zaxis_data,1)/), &
                        &new_axis_id=subaxis_id, zbounds=zbounds, nz_subaxis=nz_subaxis)
    var_axis_ids(zaxis_index) = subaxis_id

  end subroutine

  !> @brief Check if one axis is the parent of another.
  !!
  !! Determines parent-child relationships between axes by checking if axis_id
  !! is a subaxis of parent_axis_id.
  function is_parent_axis(axis_id, parent_axis_id, diag_axis) &
    result(rslt)
    integer, intent(in) :: axis_id        !< Axis id to check
    integer, intent(in) :: parent_axis_id !< Axis id of the parent to check
    class(fmsDiagAxisContainer_type), target, intent(in) :: diag_axis(:)    !< Array of diag_axis objects

    logical :: rslt

    rslt = .false.
    select type(axis => diag_axis(axis_id)%axis)
    type is (fmsDiagSubAxis_type)
      if (axis%parent_axis_id .eq. parent_axis_id) rslt = .true.
    end select
  end function is_parent_axis

  !> @brief Find the name of a Z-axis subaxis by matching parent and bounds.
  !!
  !! Searches the file's registered axes to find a Z subaxis with the specified
  !! parent axis and Z bounds. The axis name is returned in dim_name.
  subroutine find_z_sub_axis_name(dim_name, parent_axis_id, file_axis_id, field_yaml, diag_axis)
    character(len=*),                intent(inout) :: dim_name         !< Name of z subaxis
    integer,                         intent(in)    :: parent_axis_id   !< Axis id of the parent
    integer,                         intent(in)    :: file_axis_id(:)  !< Axis ids of the file
    type(diagYamlFilesVar_type),     intent(in)    :: field_yaml       !< Field info from diag_table yaml
    class(fmsDiagAxisContainer_type),intent(in)    :: diag_axis(:)     !< Array of axis objections

    integer :: id
    integer :: i

    do i = 1, size(file_axis_id)
      id = file_axis_id(i)
      select type (axis_ptr => diag_axis(id)%axis)
      type is (fmsDiagSubAxis_type)
        if (axis_ptr%parent_axis_id .eq. parent_axis_id) then
          if (axis_ptr%is_same_zbounds(field_yaml%get_var_zbounds())) then
            dim_name = axis_ptr%subaxis_name
            return
          endif
        endif
      end select
    enddo
    call mpp_error(FATAL, "Unable to determine the z subaxis name for field "//&
                          trim(field_yaml%get_var_varname())//" in file: "//&
                          trim(field_yaml%get_var_fname()))
  end subroutine
#endif
end module fms_diag_axis_object_mod
!> @}
! close documentation grouping
