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

!> @defgroup fms_diag_axis_object_mod fms_diag_axis_object_mod
!> @ingroup diag_manager
!! @brief fms_diag_axis_object_mod stores the diag axis object, a diag domain
!! object, and a subaxis object.

!> @file
!> @brief File for @ref diag_axis_object_mod

!> @addtogroup fms_diag_axis_object_mod
!> @{
module fms_diag_axis_object_mod
#ifdef use_yaml
  use mpp_domains_mod, only:  domain1d, domain2d, domainUG, mpp_get_compute_domain, CENTER, &
                            & mpp_get_global_domain, NORTH, EAST, mpp_get_tile_id, &
                            & mpp_get_ntile_count, mpp_get_io_domain
  use platform_mod,    only:  r8_kind, r4_kind, i4_kind, i8_kind
  use diag_data_mod,   only:  diag_atttype, max_axes, NO_DOMAIN, TWO_D_DOMAIN, UG_DOMAIN, &
                              direction_down, direction_up, fmsDiagAttribute_type, max_axis_attributes, &
                              DIAG_NULL, index_gridtype, latlon_gridtype, pack_size_str, &
                              get_base_year, get_base_month, get_base_day, get_base_hour, get_base_minute,&
                              get_base_second, is_x_axis, is_y_axis
  use mpp_mod,         only:  FATAL, mpp_error, uppercase, mpp_pe, mpp_root_pe, stdout
  use fms2_io_mod,     only:  FmsNetcdfFile_t, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t, &
                            & register_axis, register_field, register_variable_attribute, write_data
  use fms_diag_yaml_mod, only: subRegion_type, diag_yaml, MAX_SUBAXES
  use diag_grid_mod,       only:  get_local_indices_cubesphere => get_local_indexes
  use axis_utils2_mod,   only: nearest_index
  implicit none

  PRIVATE

  public :: fmsDiagAxis_type, fms_diag_axis_object_init, fms_diag_axis_object_end, &
          & get_domain_and_domain_type, diagDomain_t, &
          & DIAGDOMAIN2D_T, fmsDiagSubAxis_type, fmsDiagAxisContainer_type, fmsDiagFullAxis_type, DIAGDOMAINUG_T
  public :: define_new_axis, parse_compress_att, get_axis_id_from_name, define_diurnal_axis, &
          & fmsDiagDiurnalAxis_type, create_new_z_subaxis, is_parent_axis, define_new_subaxis_latlon, &
          & define_new_subaxis_index

  !> @}

  !> @brief Type to hold the domain info for an axis
  !! This type was created to avoid having to send in "Domain", "Domain2", "DomainUG" as arguments into subroutines
  !! and instead only 1 class(diagDomain_t) argument can be send
  !> @ingroup diag_axis_object_mod
  type diagDomain_t
    contains
      procedure :: set => set_axis_domain
      procedure :: length => get_length
      procedure :: get_ntiles
  end type diagDomain_t

  !> @brief Type to hold the 1d domain
  type, extends(diagDomain_t) :: diagDomain1d_t
     type(domain1d) :: Domain !< 1d Domain of the axis
  end type

  !> @brief Type to hold the 2d domain
  type, extends(diagDomain_t) :: diagDomain2d_t
    type(domain2d) :: Domain2 !< 2d Domain of an "X" or "Y" axis
  end type

  !> @brief Type to hold the unstructured domain
  type, extends(diagDomain_t) :: diagDomainUg_t
    type(domainUG) :: DomainUG !< Domain of "U" axis
  end type

  !> @brief Type to hold the diagnostic axis description.
  !> @ingroup diag_axis_object_mod
  TYPE :: fmsDiagAxis_type
     INTEGER                        , private :: axis_id         !< ID of the axis

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

  !> @brief Type to hold the diag_axis (either subaxis or a full axis)
  !> @ingroup diag_axis_object_mod
  type :: fmsDiagAxisContainer_type
    class(fmsDiagAxis_type), allocatable :: axis
  end type

  !> @brief Type to hold the subaxis
  !> @ingroup diag_axis_object_mod
  TYPE, extends(fmsDiagAxis_type) :: fmsDiagSubAxis_type
    CHARACTER(len=:),  ALLOCATABLE , private  :: subaxis_name   !< Name of the subaxis
    INTEGER                        , private  :: starting_index !< Starting index of the subaxis relative to the
                                                                !! parent axis
    INTEGER                        , private  :: ending_index   !< Ending index of the subaxis relative to the
                                                                !! parent axis
    INTEGER                        , private  :: parent_axis_id !< Id of the parent_axis
    INTEGER                        , private  :: compute_idx(2) !< Starting and ending index of the compute domain
    INTEGER,            allocatable, private  :: global_idx(:)  !< Starting and ending index of the global domain
    real(kind=r4_kind), allocatable, private  :: zbounds(:)     !< Bounds of the Z axis
    contains
      procedure :: fill_subaxis
      procedure :: axis_length
      procedure :: get_starting_index
      procedure :: get_ending_index
      procedure :: get_compute_indices
  END TYPE fmsDiagSubAxis_type

  !> @brief Type to hold the diurnal axis
  !> @ingroup diag_axis_object_mod
  TYPE, extends(fmsDiagAxis_type) :: fmsDiagDiurnalAxis_type
    INTEGER                      , private :: ndiurnal_samples !< The number of diurnal samples
    CHARACTER(len=:), ALLOCATABLE, private :: axis_name        !< The diurnal axis name
    CHARACTER(len=:), ALLOCATABLE, private :: long_name        !< The longname of the diurnal axis
    CHARACTER(len=:), ALLOCATABLE, private :: units            !< The units
    INTEGER                      , private :: edges_id         !< The id of the diurnal edges
    CHARACTER(len=:), ALLOCATABLE, private :: edges_name       !< The name of the edges axis
    CLASS(*),         ALLOCATABLE, private :: diurnal_data(:)  !< The diurnal data

    contains
      procedure :: get_diurnal_axis_samples
      procedure :: write_diurnal_metadata
  END TYPE fmsDiagDiurnalAxis_type

  !> @brief Type to hold the diagnostic axis description.
  !> @ingroup diag_axis_object_mod
  TYPE, extends(fmsDiagAxis_type) :: fmsDiagFullAxis_type
     CHARACTER(len=:),   ALLOCATABLE, private :: axis_name       !< Name of the axis
     CHARACTER(len=:),   ALLOCATABLE, private :: units           !< Units of the axis
     CHARACTER(len=:),   ALLOCATABLE, private :: long_name       !< Long_name attribute of the axis
     CHARACTER(len=1)               , private :: cart_name       !< Cartesian name "X", "Y", "Z", "T", "U", "N"
     CLASS(*),           ALLOCATABLE, private :: axis_data(:)    !< Data of the axis
     CHARACTER(len=:),   ALLOCATABLE, private :: type_of_data    !< The type of the axis_data ("float" or "double")
     !< TO DO this can be a dlinked to avoid having limits
     integer,            ALLOCATABLE, private :: subaxis(:)      !< Array of subaxis
     integer                        , private :: nsubaxis        !< Number of subaxis
     class(diagDomain_t),ALLOCATABLE, private :: axis_domain     !< Domain
     INTEGER                        , private :: type_of_domain  !< The type of domain ("NO_DOMAIN", "TWO_D_DOMAIN",
                                                                 !! or "UG_DOMAIN")
     INTEGER                        , private :: length          !< Global axis length
     INTEGER                        , private :: direction       !< Direction of the axis 0, 1, -1
     INTEGER,            ALLOCATABLE, private :: edges_id        !< Axis ID for the edges axis
                                                                 !! This axis will be written to the file
     CHARACTER(len=:),   ALLOCATABLE, private :: edges_name      !< Name for the previously defined "edges axis"
                                                                 !! This will be written as an attribute
     CHARACTER(len=:), ALLOCATABLE,   private :: aux             !< Auxiliary name, can only be <TT>geolon_t</TT>
                                                                 !! or <TT>geolat_t</TT>
     CHARACTER(len=128)             , private :: req             !< Required field names.
     INTEGER                        , private :: tile_count      !< The number of tiles
     TYPE(fmsDiagAttribute_type),allocatable , private :: attributes(:) !< Array to hold user definable attributes
     INTEGER                        , private :: num_attributes  !< Number of defined attibutes
     INTEGER                        , private :: domain_position !< The position in the doman (NORTH, EAST or CENTER)
     integer, allocatable           , private :: structured_ids(:) !< If the axis is in the unstructured grid,
                                                                   !! this is the axis ids of the structured axis
     CHARACTER(len=:), ALLOCATABLE,   private :: set_name        !< Name of the axis set. This is to distinguish
                                                                 !! two axis with the same name

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
     ! TO DO:
     ! Get/has/is subroutines as needed
  END TYPE fmsDiagFullAxis_type

  !> @addtogroup fms_diag_yaml_mod
  !> @{
  contains

  !!!!!!!!!!!!!!!!! DIAG AXIS PROCEDURES !!!!!!!!!!!!!!!!!
  !> @brief Initialize the axis
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

  !> @brief Add an attribute to an axis
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

  !> @brief Write the axis meta data to an open fileobj
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

  !> @brief Write the axis data to an open fms2io_fileobj
  subroutine write_axis_data(this, fms2io_fileobj, parent_axis)
    class(fmsDiagAxis_type),           target, INTENT(IN)    :: this        !< diag_axis obj
    class(FmsNetcdfFile_t),                    INTENT(INOUT) :: fms2io_fileobj!< Fms2_io fileobj to write the data to
    class(fmsDiagAxis_type), OPTIONAL, target, INTENT(IN)    :: parent_axis !< The parent axis if this is a subaxis

    integer                       :: i                 !< Starting index of a sub_axis
    integer                       :: j                 !< Ending index of a sub_axis
    integer                       :: global_io_index(2)!< Global io domain starting and ending index
    select type(this)
    type is (fmsDiagFullAxis_type)
      call this%get_global_io_domain(global_io_index)
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


  !> @brief Defined a new diurnal axis
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
  subroutine get_global_io_domain(this, global_io_index)
    class(fmsDiagFullAxis_type), intent(in)  :: this               !< diag_axis obj
    integer,                     intent(out) :: global_io_index(2) !< Global io domain starting and ending index

    type(domain2d), pointer :: io_domain !< pointer to the io domain

    global_io_index(1) = 1
    global_io_index(2) = this%length

    if (allocated(this%axis_domain)) then
      select type(domain => this%axis_domain)
      type is (diagDomain2d_t)
        io_domain => mpp_get_io_domain(domain%domain2)
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
        !< If the the tile number is present and the current PE is not on the tile, then there is no need
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

  !!!!!!!!!!!!!!!!!! SUB AXIS PROCEDURES !!!!!!!!!!!!!!!!!
  !> @brief Fills in the information needed to define a subaxis
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

  !> @brief Get the axis length of a subaxis
  !> @return the axis length
  function axis_length(this) &
    result(res)
      class(fmsDiagSubAxis_type)  , INTENT(IN) :: this             !< diag_sub_axis obj
      integer :: res

      res = this%ending_index - this%starting_index + 1
    end function

   !> @brief Accesses its member starting_index
  !! @return a copy of the starting_index
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
    character(len=*),    INTENT(IN)    :: cart_axis !< cart_axis of the axis
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

  !> @brief Set the axis domain
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

  !< @brief Allocates the array of axis/subaxis objects
  !! @return true if there the aray of axis/subaxis objects is allocated
  logical function fms_diag_axis_object_init(axis_array)
    class(fmsDiagAxisContainer_type)   , allocatable, intent(inout) :: axis_array(:) !< Array of diag_axis

    if (allocated(axis_array)) call mpp_error(FATAL, "The diag_axis containers is already allocated")
    allocate(axis_array(max_axes))
    !axis_array%axis_id = DIAG_NULL

    fms_diag_axis_object_init = .true.
  end function fms_diag_axis_object_init

  !< @brief Deallocates the array of axis/subaxis objects
  !! @return false if the aray of axis/subaxis objects was allocated
  logical function fms_diag_axis_object_end(axis_array)
    class(fmsDiagAxisContainer_type)   , allocatable, intent(inout) :: axis_array(:) !< Array of diag_axis

    if (allocated(axis_array)) deallocate(axis_array)
    fms_diag_axis_object_end = .false.

  end function fms_diag_axis_object_end

  !< @brief Determine the axis name of an axis_object
  !! @return The name of the axis
  !! @note This function may be called from the field object (i.e. to determine the dimension names for io),
  !! The field object only contains the parent axis ids, because the subregion is defined in a per file basis,
  !! so the is_regional flag is needed so that the correct axis name can be used
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

  !> @brief Check if a cart_name is valid and crashes if it isn't
  subroutine check_if_valid_cart_name(cart_name)
    character(len=*), intent(in) :: cart_name

    select case (cart_name)
    case ("X", "Y", "Z", "T", "U", "N")
    case default
      call mpp_error(FATAL, "diag_axit_init: Invalid cart_name: "//cart_name//&
                             "The acceptable values are X, Y, Z, T, U, N.")
    end select
  end subroutine check_if_valid_cart_name

  !> @brief Check if a domain_position is valid and crashes if it isn't
  subroutine check_if_valid_domain_position(domain_position)
    integer, INTENT(IN) :: domain_position

    select case (domain_position)
    case (CENTER, NORTH, EAST)
    case default
        call mpp_error(FATAL, "diag_axit_init: Invalid domain_positon. "&
                             "The acceptable values are NORTH, EAST, CENTER")
    end select
  end subroutine check_if_valid_domain_position

  !> @brief Check if a direction is valid and crashes if it isn't
  subroutine check_if_valid_direction(direction)
    integer, INTENT(IN) :: direction

    select case(direction)
    case(-1, 0, 1)
    case default
      call mpp_error(FATAL, "diag_axit_init: Invalid direction. "&
                             "The acceptable values are-1 0 1")
    end select
  end subroutine check_if_valid_direction

  !> @brief Loop through a variable's axis_id to determine and return the domain type and domain to use
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
      type is (fmsDiagFullAxis_type)
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

  !> @brief Fill in the subaxis object for a subRegion defined by index
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

  !> @brief Fill in the subaxis object for a subRegion defined by lat lon
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

  !> @brief Creates a new subaxis and fills it will all the information it needs
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

  !< @brief Determine the parent_axis_id of a subaxis
  !! @return parent_axis_id if it is a subaxis and diag_null if is not a subaxis
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

  !< @brief Determine the most recent subaxis id in a diag_axis object
  !! @return the most recent subaxis id in a diag_axis object
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

  !< @brief Parses the "compress" attribute to get the names of the two axis
  !! @return the names of the structured axis
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

  !< @brief Determine the axis id of a axis
  !! @return Axis id
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

  !< @brief Get the number of diurnal samples for a diurnal axis
  !! @return The number of diurnal samples
  pure function get_diurnal_axis_samples(this) &
  result(n_diurnal_samples)

    class(fmsDiagDiurnalAxis_type), intent(in) :: this !< Axis Object
    integer :: n_diurnal_samples

    n_diurnal_samples = this%ndiurnal_samples
  end function get_diurnal_axis_samples

  !< @brief Writes out the metadata for a diurnal axis
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

  !> @brief Creates a new z subaxis to use
  subroutine create_new_z_subaxis(zbounds, var_axis_ids, diag_axis, naxis, file_axis_id, nfile_axis, nz_subaxis)
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

    class(*), pointer :: zaxis_data(:)      !< The data of the full zaxis
    integer           :: subaxis_indices(2) !< The starting and ending indices of the subaxis relative to the full
                                            !! axis
    integer           :: i                  !< For do loops
    integer           :: subaxis_id         !< The id of the new z subaxis
    logical           :: axis_found         !< Flag that indicated if the zsubaxis already exists

    !< Determine if the axis was already created
    axis_found = .false.
    do i = 1, nfile_axis
      select type (axis => diag_axis(file_axis_id(i))%axis)
      type is (fmsDiagSubAxis_type)
        if (axis%zbounds(1) .eq. zbounds(1) .and. axis%zbounds(2) .eq. zbounds(2)) then
          axis_found = .true.
          subaxis_id = file_axis_id(i)
          exit
        endif
      end select
    enddo

    !< Determine which of the variable's axis is the zaxis!
    do i = 1, size(var_axis_ids)
      select type (parent_axis => diag_axis(var_axis_ids(i))%axis)
      type is (fmsDiagFullAxis_type)
        if (parent_axis%cart_name .eq. "Z") then
          !< If the axis was previously defined set the var_axis_ids and leave
          if (axis_found) then
            var_axis_ids(i) = subaxis_id
            return
          endif
          zaxis_data => parent_axis%axis_data

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
          var_axis_ids(i) = subaxis_id
          return
        endif
      end select
    enddo

  end subroutine

  !> @brief Determine if the diag_axis(parent_axis_id) is the parent of diag_axis(axis_id)
  !! @return .True. if diag_axis(parent_axis_id) is the parent of diag_axis(axis_id)
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

#endif
end module fms_diag_axis_object_mod
!> @}
! close documentation grouping
