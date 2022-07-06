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
  use mpp_domains_mod, only:  domain1d, domain2d, domainUG, mpp_get_compute_domain, CENTER, &
                            & mpp_get_compute_domain, NORTH, EAST
  use platform_mod,    only:  r8_kind, r4_kind, i4_kind, i8_kind
  use diag_data_mod,   only:  diag_atttype, max_axes, NO_DOMAIN, TWO_D_DOMAIN, UG_DOMAIN, &
                              direction_down, direction_up, fmsDiagAttribute_type, max_axis_attributes
  use mpp_mod,         only:  FATAL, mpp_error, uppercase
  use fms2_io_mod,     only:  FmsNetcdfFile_t, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t, &
                            & register_axis, register_field, register_variable_attribute, write_data
  implicit none

  PRIVATE

  public :: diagAxis_t, set_subaxis, fms_diag_axis_init, fms_diag_axis_object_init, fms_diag_axis_object_end, &
          & get_domain_and_domain_type, axis_obj, diagDomain_t, sub_axis_objs, fms_diag_axis_add_attribute, &
          & DIAGDOMAIN2D_T, fms_get_axis_length
  !> @}

  !> @brief Type to hold the domain info for an axis
  !! This type was created to avoid having to send in "Domain", "Domain2", "DomainUG" as arguments into subroutines
  !! and instead only 1 class(diagDomain_t) argument can be send
  !> @ingroup diag_axis_object_mod
  type diagDomain_t
    contains
      procedure :: set => set_axis_domain
      procedure :: length => get_length
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

  !> @brief Type to hold the subaxis
  !> @ingroup diag_axis_object_mod
  TYPE subaxis_t
    CHARACTER(len=:), ALLOCATABLE :: subaxis_name   !< Name of the subaxis
    INTEGER                       :: starting_index !< Starting index of the subaxis relative to the parent axis
    INTEGER                       :: ending_index   !< Ending index of the subaxis relative to the parent axis
    class(*)        , ALLOCATABLE :: bounds         !< Bounds of the subaxis (lat/lon or indices)
    INTEGER                       :: parent_axis_id !< Id of the parent_axis
    contains
      procedure :: exists => check_if_subaxis_exists
  END TYPE subaxis_t

  !> @brief Type to hold the diagnostic axis description.
  !> @ingroup diag_axis_object_mod
  TYPE diagAxis_t
     CHARACTER(len=:),   ALLOCATABLE, private :: axis_name       !< Name of the axis
     CHARACTER(len=:),   ALLOCATABLE, private :: units           !< Units of the axis
     CHARACTER(len=:),   ALLOCATABLE, private :: long_name       !< Long_name attribute of the axis
     CHARACTER(len=1)               , private :: cart_name       !< Cartesian name "X", "Y", "Z", "T", "U", "N"
     CLASS(*),           ALLOCATABLE, private :: axis_data(:)    !< Data of the axis
     CHARACTER(len=:),   ALLOCATABLE, private :: type_of_data    !< The type of the axis_data ("float" or "double")
     !< TO DO this can be a dlinked to avoid having limits
     type(subaxis_t)                , private :: subaxis(3)      !< Array of subaxis
     integer                        , private :: nsubaxis        !< Number of subaxis
     class(diagDomain_t),ALLOCATABLE, private :: axis_domain     !< Domain
     INTEGER                        , private :: type_of_domain  !< The type of domain ("NO_DOMAIN", "TWO_D_DOMAIN",
                                                                 !! or "UG_DOMAIN")
     INTEGER                        , private :: length          !< Global axis length
     INTEGER                        , private :: direction       !< Direction of the axis 0, 1, -1
     INTEGER                        , private :: edges           !< Axis ID for the previously defined "edges axis"
     CHARACTER(len=128)             , private :: aux             !< Auxiliary name, can only be <TT>geolon_t</TT>
                                                                 !! or <TT>geolat_t</TT>
     CHARACTER(len=128)             , private :: req             !< Required field names.
     INTEGER                        , private :: tile_count      !< The number of tiles
     TYPE(fmsDiagAttribute_type),allocatable , private :: attributes(:) !< Array to hold user definable attributes
     INTEGER                        , private :: num_attributes  !< Number of defined attibutes
     INTEGER                        , private :: domain_position !< The position in the doman (NORTH, EAST or CENTER)

     contains

     PROCEDURE :: add_axis_attribute
     PROCEDURE :: register => register_diag_axis_obj
     PROCEDURE :: axis_length => get_axis_length
     PROCEDURE :: set_subaxis
     PROCEDURE :: write_axis_metadata
     PROCEDURE :: write_axis_data

     ! TO DO:
     ! Get/has/is subroutines as needed
  END TYPE diagAxis_t

  integer                                :: number_of_axis !< Number of axis that has been registered
  type(diagAxis_t), ALLOCATABLE, TARGET  :: axis_obj(:)    !< Diag_axis objects
  logical                                :: module_is_initialized !< Flag indicating if the module is initialized
  integer                                :: nsubaxis_objs  !< Number of sub_axis that has been registered
  type(subaxis_t), ALLOCATABLE, Target   :: sub_axis_objs(:) !< Registered sub_axis objects

  !> @addtogroup fms_diag_yaml_mod
  !> @{
  contains

  !!!!!!!!!!!!!!!!! DIAG AXIS PROCEDURES !!!!!!!!!!!!!!!!!
  !> @brief Initialize the axis
  subroutine register_diag_axis_obj(obj, axis_name, axis_data, units, cart_name, long_name, direction,&
  & set_name, edges, Domain, Domain2, DomainU, aux, req, tile_count, domain_position )
    class(diagAxis_t),  INTENT(out)          :: obj             !< Diag_axis obj
    CHARACTER(len=*),   INTENT(in)           :: axis_name       !< Name of the axis
    class(*),           INTENT(in)           :: axis_data(:)    !< Array of coordinate values
    CHARACTER(len=*),   INTENT(in)           :: units           !< Units for the axis
    CHARACTER(len=1),   INTENT(in)           :: cart_name       !< Cartesian axis ("X", "Y", "Z", "T", "U", "N")
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: long_name       !< Long name for the axis.
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: set_name        !< Name of the parent axis, if it is a subaxis
    INTEGER,            INTENT(in), OPTIONAL :: direction       !< Indicates the direction of the axis
    INTEGER,            INTENT(in), OPTIONAL :: edges           !< Axis ID for the previously defined "edges axis"
    TYPE(domain1d),     INTENT(in), OPTIONAL :: Domain          !< 1D domain
    TYPE(domain2d),     INTENT(in), OPTIONAL :: Domain2         !< 2D domain
    TYPE(domainUG),     INTENT(in), OPTIONAL :: DomainU         !< Unstructured domain
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: aux             !< Auxiliary name, can only be <TT>geolon_t</TT>
                                                                !! or <TT>geolat_t</TT>
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: req             !< Required field names.
    INTEGER,            INTENT(in), OPTIONAL :: tile_count      !< Number of tiles
    INTEGER,            INTENT(in), OPTIONAL :: domain_position !< Domain position, "NORTH" or "EAST"

    obj%axis_name = trim(axis_name)
    obj%units = trim(units)
    obj%cart_name = uppercase(cart_name)
    call check_if_valid_cart_name(obj%cart_name)

    if (present(long_name)) obj%long_name = trim(long_name)

    select type (axis_data)
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: obj%axis_data(size(axis_data)))
      obj%axis_data = axis_data
      obj%type_of_data = "double" !< This is what fms2_io expects in the register_field call
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: obj%axis_data(size(axis_data)))
      obj%axis_data = axis_data
      obj%type_of_data = "float" !< This is what fms2_io expects in the register_field call
    class default
      call mpp_error(FATAL, "The axis_data in your diag_axis_init call is not a supported type. &
                          &  Currently only r4 and r8 data is supported.")
    end select

    obj%type_of_domain = NO_DOMAIN
    if (present(Domain)) then
      if (present(Domain2) .or. present(DomainU)) call mpp_error(FATAL, &
        "The presence of Domain with any other domain type is prohibited. "//&
        "Check you diag_axis_init call for axis_name:"//trim(axis_name))
      allocate(diagDomain1d_t :: obj%axis_domain)
      call obj%axis_domain%set(Domain=Domain)
    else if (present(Domain2)) then
        if (present(DomainU)) call mpp_error(FATAL, &
        "The presence of Domain2 with any other domain type is prohibited. "//&
        "Check you diag_axis_init call for axis_name:"//trim(axis_name))
      allocate(diagDomain2d_t :: obj%axis_domain)
      call obj%axis_domain%set(Domain2=Domain2)
      obj%type_of_domain = TWO_D_DOMAIN
    else if (present(DomainU)) then
      allocate(diagDomainUg_t :: obj%axis_domain)
      call obj%axis_domain%set(DomainU=DomainU)
      obj%type_of_domain = UG_DOMAIN
    endif

    obj%tile_count = 1
    if (present(tile_count)) obj%tile_count = tile_count

    obj%domain_position = CENTER
    if (present(domain_position)) obj%domain_position = domain_position
    call check_if_valid_domain_position(obj%domain_position)

    obj%length = size(axis_data)

    obj%direction = 0
    if (present(direction)) obj%direction = direction
    call check_if_valid_direction(obj%direction)

    obj%edges = 0
    if (present(edges)) obj%edges = edges
    call check_if_valid_edges(obj%edges)

    if (present(aux)) obj%aux = trim(aux)
    if (present(req)) obj%req = trim(req)

    obj%nsubaxis = 0
    obj%num_attributes = 0
  end subroutine register_diag_axis_obj

  !> @brief Add an attribute to an axis
  subroutine add_axis_attribute(obj, att_name, att_value)
    class(diagAxis_t),INTENT(INOUT) :: obj          !< diag_axis obj
    character(len=*), intent(in)    :: att_name     !< Name of the attribute
    class(*),         intent(in)    :: att_value(:) !< The attribute value to add

    integer :: j    !< obj%num_attributes (for less typing)

    if (.not. allocated(obj%attributes)) &
      allocate(obj%attributes(max_axis_attributes))

    obj%num_attributes = obj%num_attributes + 1

    j = obj%num_attributes
    call obj%attributes(j)%add(att_name, att_value)
  end subroutine add_axis_attribute

  !> @brief Write the axis meta data to an open fileobj
  subroutine write_axis_metadata(obj, fileobj, sub_axis_id)
    class(diagAxis_t), target, INTENT(IN)    :: obj         !< diag_axis obj
    class(FmsNetcdfFile_t),    INTENT(INOUT) :: fileobj     !< Fms2_io fileobj to write the data to
    integer, OPTIONAL,         INTENT(IN)    :: sub_axis_id !< ID of the sub_axis, if it exists

    character(len=:), ALLOCATABLE :: axis_edges_name !< Name of the edges, if it exist
    character(len=:), pointer     :: axis_name       !< Name of the axis
    integer                       :: axis_length     !< Size of the axis
    integer                       :: i               !< For do loops

    if (present(sub_axis_id)) then
      axis_name  => obj%subaxis(sub_axis_id)%subaxis_name
      axis_length = obj%subaxis(sub_axis_id)%ending_index - obj%subaxis(sub_axis_id)%starting_index + 1
    else
      axis_name => obj%axis_name
      axis_length = obj%length
    endif

    !< Add the axis as a dimension in the netcdf file based on the type of axis_domain and the fileobj type
    select type (fileobj)
      type is (FmsNetcdfFile_t)
        !< Here the axis is not domain decomposed (i.e z_axis)
        call register_axis(fileobj, axis_name, axis_length)
      type is (FmsNetcdfDomainFile_t)
        select case (obj%type_of_domain)
        case (NO_DOMAIN)
          !< Here the fileobj is domain decomposed, but the axis is not
          !! Domain decomposed fileobjs can have axis that are not domain decomposed (i.e "Z" axis)
          call register_axis(fileobj, axis_name, axis_length)
        case (TWO_D_DOMAIN)
          !< Here the axis is domain decomposed
          call register_axis(fileobj, axis_name, obj%cart_name, domain_position=obj%domain_position)
        end select
      type is (FmsNetcdfUnstructuredDomainFile_t)
        select case (obj%type_of_domain)
        case (NO_DOMAIN)
          !< Here the fileobj is in the unstructured domain, but the axis is not
          !< Unstructured domain fileobjs can have axis that are not domain decomposed (i.e "Z" axis)
          call register_axis(fileobj, axis_name, axis_length)
        case (UG_DOMAIN)
          !< Here the axis is in a unstructured domain
          call register_axis(fileobj, axis_name)
        end select
    end select

    !< Add the axis as a variable and write its metada
    call register_field(fileobj, axis_name, obj%type_of_data, (/axis_name/))
    call register_variable_attribute(fileobj, axis_name, "longname", obj%long_name, &
      str_len=len_trim(obj%long_name))

    if (obj%cart_name .NE. "N") &
      call register_variable_attribute(fileobj, axis_name, "axis", obj%cart_name, str_len=1)

    if (trim(obj%units) .NE. "none") &
      call register_variable_attribute(fileobj, axis_name, "units", obj%units, str_len=len_trim(obj%units))

    select case (obj%direction)
    case (direction_up)
      call register_variable_attribute(fileobj, axis_name, "positive", "up", str_len=2)
    case (direction_down)
      call register_variable_attribute(fileobj, axis_name, "positive", "down", str_len=4)
    end select

    if (obj%edges > 0) then
      axis_edges_name = axis_obj(obj%edges)%axis_name
      call register_variable_attribute(fileobj, axis_name, "edges", axis_edges_name, &
        str_len=len_trim(axis_edges_name))
    endif

    if(allocated(obj%attributes)) then
      do i = 1, size(obj%attributes)
        call register_variable_attribute(fileobj, axis_name, obj%attributes(i)%att_name, &
          & obj%attributes(i)%att_value)
      enddo
    endif

  end subroutine write_axis_metadata

  !> @brief Write the axis data to an open fileobj
  subroutine write_axis_data(obj, fileobj, sub_axis_id)
    class(diagAxis_t),      INTENT(IN)    :: obj       !< diag_axis obj
    class(FmsNetcdfFile_t), INTENT(INOUT) :: fileobj   !< Fms2_io fileobj to write the data to
    integer, OPTIONAL,      INTENT(IN)    :: sub_axis_id !< ID of the sub_axis, if it exists

    integer                       :: i         !< Starting index of a sub_axis
    integer                       :: j         !< Ending index of a sub_axis

    if (present(sub_axis_id)) then
      i = obj%subaxis(sub_axis_id)%starting_index
      j = obj%subaxis(sub_axis_id)%ending_index

      call write_data(fileobj, obj%subaxis(sub_axis_id)%subaxis_name, obj%axis_data(i:j))
    else
      call write_data(fileobj, obj%axis_name, obj%axis_data)
    endif
  end subroutine write_axis_data

  !> @brief Get the length of the axis
  !> @return axis length
  function get_axis_length(obj) &
  result (axis_length)
    class(diagAxis_t), intent(inout) :: obj !< diag_axis obj
    integer                           :: axis_length

    !< If the axis is domain decomposed axis_length will be set to the length for the current PE:
    if (allocated(obj%axis_domain)) then
      axis_length = obj%axis_domain%length(obj%cart_name, obj%domain_position, obj%length)
    else
      axis_length = obj%length
    endif

  end function

  !> @brief Set the subaxis of the axis obj
  !> @return A sub_axis id corresponding to the indices of the sub_axes in the sub_axes_objs array
  function set_subaxis(obj, bounds) &
  result(sub_axes_id)
    class(diagAxis_t),  INTENT(INOUT) :: obj       !< diag_axis obj
    class(*),           INTENT(INOUT) :: bounds(:) !< bound of the subaxis

    integer :: sub_axes_id

    integer :: i !< For do loops

    !< Check if the subaxis for this bouds already exists
    do i = 1, obj%nsubaxis
      if (obj%subaxis(i)%exists(bounds)) return
    enddo

    !< TO DO: everything
    obj%nsubaxis = obj%nsubaxis + 1

    nsubaxis_objs = nsubaxis_objs + 1
    sub_axes_id = nsubaxis_objs
    !< TO DO: set the parent_axis_id
  end function

  !!!!!!!!!!!!!!!!!! SUB AXIS PROCEDURES !!!!!!!!!!!!!!!!!
  !> @brief Check if a subaxis was already defined
  !> @return Flag indicating if a subaxis is already defined
  function check_if_subaxis_exists(obj,bounds) &
  result(exists)
    class(subaxis_t), INTENT(INOUT) :: obj       !< diag_axis obj
    class(*),         INTENT(IN)    :: bounds(:) !< bounds of the subaxis
    logical                         :: exists

    !< TO DO: compare bounds
    exists = .false.
  end function

  !> @brief Get the length of a 2D domain
  !> @return Length of the 2D domain
  function get_length(obj, cart_axis, domain_position, global_length) &
  result (length)
    class(diagDomain_t), INTENT(INOUT) :: obj       !< diag_axis obj
    character(len=*),    INTENT(IN)    :: cart_axis !< cart_axis of the axis
    integer,             INTENT(IN)    :: domain_position !< Domain position (CENTER, NORTH, EAST)
    integer,             INTENT(IN)    :: global_length !< global_length of the axis

    integer :: length

    select type (obj)
    type is(diagDomain2d_t)
      if (trim(cart_axis) == "X") call mpp_get_compute_domain(obj%Domain2, xsize=length, position=domain_position)
      if (trim(cart_axis) == "Y") call mpp_get_compute_domain(obj%Domain2, ysize=length, position=domain_position)
    class default
      !< If domain is 1D or UG, just set it to the global length
      length = global_length
    end select
  end function

  !!!!!!!!!!!!!!!!! FMS_DOMAIN PROCEDURES !!!!!!!!!!!!!!!!!

  !> @brief Set the axis domain
  subroutine set_axis_domain(obj, Domain, Domain2, DomainU)
    class(diagDomain_t) :: obj !< fms_domain obj
    TYPE(domain1d),     INTENT(in),  OPTIONAL :: Domain  !< 1d domain
    TYPE(domain2d),     INTENT(in),  OPTIONAL :: Domain2 !< 2d domain
    TYPE(domainUG),     INTENT(in),  OPTIONAL :: DomainU !< Unstructured domain

    select type(obj)
    type is (diagDomain1d_t)
      obj%Domain = Domain
    type is (diagDomain2d_t)
      obj%Domain2 = Domain2
    type is (diagDomainUg_t)
      obj%DomainUG = DomainU
    end select
  end subroutine set_axis_domain

  subroutine fms_diag_axis_object_init()

    if (module_is_initialized) return

    number_of_axis = 0
    allocate(axis_obj(max_axes))

    module_is_initialized = .true.
  end subroutine fms_diag_axis_object_init

  subroutine fms_diag_axis_object_end()
    deallocate(axis_obj)

    module_is_initialized = .false.
  end subroutine fms_diag_axis_object_end

  !> @brief Wrapper for the register_diag_axis subroutine. This is needed to keep the diag_axis_init
  !! interface the same
  !> @return Axis id
  FUNCTION fms_diag_axis_init(axis_name, axis_data, units, cart_name, long_name, direction,&
    & set_name, edges, Domain, Domain2, DomainU, aux, req, tile_count, domain_position ) &
    & result(id)

    CHARACTER(len=*),   INTENT(in)           :: axis_name       !< Name of the axis
    REAL,               INTENT(in)           :: axis_data(:)    !< Array of coordinate values
    CHARACTER(len=*),   INTENT(in)           :: units           !< Units for the axis
    CHARACTER(len=1),   INTENT(in)           :: cart_name       !< Cartesian axis ("X", "Y", "Z", "T", "U", "N")
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: long_name       !< Long name for the axis.
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: set_name        !< Name of the parent axis, if it is a subaxis
    INTEGER,            INTENT(in), OPTIONAL :: direction       !< Indicates the direction of the axis
    INTEGER,            INTENT(in), OPTIONAL :: edges           !< Axis ID for the previously defined "edges axis"
    TYPE(domain1d),     INTENT(in), OPTIONAL :: Domain          !< 1D domain
    TYPE(domain2d),     INTENT(in), OPTIONAL :: Domain2         !< 2D domain
    TYPE(domainUG),     INTENT(in), OPTIONAL :: DomainU         !< Unstructured domain
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: aux             !< Auxiliary name, can only be <TT>geolon_t</TT>
                                                                !! or <TT>geolat_t</TT>
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: req             !< Required field names.
    INTEGER,            INTENT(in), OPTIONAL :: tile_count      !< Number of tiles
    INTEGER,            INTENT(in), OPTIONAL :: domain_position !< Domain position, "NORTH" or "EAST"
    integer :: id

    number_of_axis = number_of_axis + 1

    if (number_of_axis > max_axes) call mpp_error(FATAL, &
      &"diag_axis_init: max_axes exceeded, increase via diag_manager_nml")

    call axis_obj(number_of_axis)%register(axis_name, axis_data, units, cart_name, long_name=long_name, &
    & direction=direction, set_name=set_name, edges=edges, Domain=Domain, Domain2=Domain2, DomainU=DomainU, aux=aux, &
    & req=req, tile_count=tile_count, domain_position=domain_position)

    id = number_of_axis
  end function

  !> @brief Add an attribute to an axis
  subroutine fms_diag_axis_add_attribute(axis_id, att_name, att_value)
    integer,          intent(in) :: axis_id      !< Id of the axis to add the attribute to
    character(len=*), intent(in) :: att_name     !< Name of the attribute
    class(*),         intent(in) :: att_value(:) !< The attribute value to add

    if (axis_id < 0 .and. axis_id > number_of_axis) &
      call mpp_error(FATAL, "diag_axis_add_attribute: The axis_id is not valid")

    call axis_obj(axis_id)%add_axis_attribute(att_name, att_value)
  end subroutine fms_diag_axis_add_attribute

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

  !> @brief Check if the edges id is valid and crashes if it isn't
  subroutine check_if_valid_edges(edges)
    integer, INTENT(IN) :: edges

    if (edges < 0 .or. edges > number_of_axis) &
       call mpp_error(FATAL, "diag_axit_init: The edge axis has not been defined. "&
                             "Call diag_axis_init for the edge axis first")
  end subroutine check_if_valid_edges

  !> @brief Loop through a variable's axis_id to determine and return the domain type and domain to use
  subroutine get_domain_and_domain_type(axis_id, domain_type, domain, var_name)
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
      !< Check that all the axis are in the same domain
      if (domain_type .ne. axis_obj(j)%type_of_domain) then
        !< If they are different domains, one of them can be NO_DOMAIN
        !! i.e a variable can have axis that are domain decomposed (x,y) and an axis that isn't (z)
        if (domain_type .eq. NO_DOMAIN .or. axis_obj(j)%type_of_domain .eq. NO_DOMAIN ) then
          !< Update the domain_type and domain, if needed
          if ((axis_obj(j)%type_of_domain .eq. TWO_D_DOMAIN  .and. size(axis_id) > 2) &
             & .or. axis_obj(j)%type_of_domain .eq. UG_DOMAIN) then
               domain_type = axis_obj(j)%type_of_domain
               domain => axis_obj(j)%axis_domain
          endif
        else
          call mpp_error(FATAL, "The variable:"//trim(var_name)//" has axis that are not in the same domain")
        endif
      endif
    enddo
  end subroutine get_domain_and_domain_type

  function fms_get_axis_length(axis_id)&
  result(axis_length)
    INTEGER, INTENT(in) :: axis_id !< Axis ID of the axis to the length of
    integer :: axis_length

    if (axis_id < 0 .and. axis_id > number_of_axis) &
      call mpp_error(FATAL, "fms_get_axis_length: The axis_id is not valid")

    axis_length = axis_obj(axis_id)%axis_length()
  end function fms_get_axis_length
end module fms_diag_axis_object_mod
!> @}
! close documentation grouping
