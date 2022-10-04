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
                              direction_down, direction_up, fmsDiagAttribute_type, max_axis_attributes, &
                              MAX_SUBAXES, DIAG_NULL
  use mpp_mod,         only:  FATAL, mpp_error, uppercase
  use fms2_io_mod,     only:  FmsNetcdfFile_t, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t, &
                            & register_axis, register_field, register_variable_attribute, write_data
  implicit none

  PRIVATE

  public :: fmsDiagAxis_type, fms_diag_axis_object_init, fms_diag_axis_object_end, &
          & get_domain_and_domain_type, diagDomain_t, &
          & DIAGDOMAIN2D_T, fmsDiagSubAxis_type, fmsDiagAxisContainer_type, fmsDiagFullAxis_type, DIAGDOMAINUG_T
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

  !> @brief Type to hold the diag_axis (either subaxis or a full axis)
  !> @ingroup diag_axis_object_mod
  type :: fmsDiagAxisContainer_type
    class(fmsDiagAxis_type), allocatable :: axis
  end type

  !> @brief Type to hold the diagnostic axis description.
  !> @ingroup diag_axis_object_mod
  TYPE fmsDiagAxis_type
     INTEGER                        , private :: axis_id         !< ID of the axis
  END TYPE fmsDiagAxis_type

  !> @brief Type to hold the subaxis
  !> @ingroup diag_axis_object_mod
  TYPE, extends(fmsDiagAxis_type) :: fmsDiagSubAxis_type
    INTEGER                      , private  :: subaxis_id     !< ID of the subaxis
    CHARACTER(len=:), ALLOCATABLE, private  :: subaxis_name   !< Name of the subaxis
    INTEGER                      , private  :: starting_index !< Starting index of the subaxis relative to the
                                                              !! parent axis
    INTEGER                      , private  :: ending_index   !< Ending index of the subaxis relative to the
                                                              !! parent axis
    class(*)        , ALLOCATABLE, private  :: bounds         !< Bounds of the subaxis (lat/lon or indices)
    INTEGER                      , private  :: parent_axis_id !< Id of the parent_axis
    contains
      procedure :: exists => check_if_subaxis_exists
  END TYPE fmsDiagSubAxis_type

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
     type(fmsDiagSubAxis_type)      , private :: subaxis(3)      !< Array of subaxis
     integer                        , private :: nsubaxis        !< Number of subaxis
     class(diagDomain_t),ALLOCATABLE, private :: axis_domain     !< Domain
     INTEGER                        , private :: type_of_domain  !< The type of domain ("NO_DOMAIN", "TWO_D_DOMAIN",
                                                                 !! or "UG_DOMAIN")
     INTEGER                        , private :: length          !< Global axis length
     INTEGER                        , private :: direction       !< Direction of the axis 0, 1, -1
     CHARACTER(len=:),   ALLOCATABLE, private :: edges_name      !< Name for the previously defined "edges axis"
                                                                 !! This will be written as an attribute
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
     PROCEDURE :: get_axis_name
     PROCEDURE :: set_edges_name
     PROCEDURE :: set_subaxis
     PROCEDURE :: write_axis_metadata
     PROCEDURE :: write_axis_data

     ! TO DO:
     ! Get/has/is subroutines as needed
  END TYPE fmsDiagFullAxis_type

  !> @addtogroup fms_diag_yaml_mod
  !> @{
  contains

  !!!!!!!!!!!!!!!!! DIAG AXIS PROCEDURES !!!!!!!!!!!!!!!!!
  !> @brief Initialize the axis
  subroutine register_diag_axis_obj(this, axis_name, axis_data, units, cart_name, long_name, direction,&
  & set_name, Domain, Domain2, DomainU, aux, req, tile_count, domain_position )
    class(fmsDiagFullAxis_type),INTENT(out)  :: this            !< Diag_axis obj
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

    this%axis_name = trim(axis_name)
    this%units = trim(units)
    this%cart_name = uppercase(cart_name)
    call check_if_valid_cart_name(this%cart_name)

    if (present(long_name)) this%long_name = trim(long_name)

    select type (axis_data)
    type is (real(kind=r8_kind))
      allocate(real(kind=r8_kind) :: this%axis_data(size(axis_data)))
      this%axis_data = axis_data
      this%type_of_data = "double" !< This is what fms2_io expects in the register_field call
    type is (real(kind=r4_kind))
      allocate(real(kind=r4_kind) :: this%axis_data(size(axis_data)))
      this%axis_data = axis_data
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

    this%length = size(axis_data)

    this%direction = 0
    if (present(direction)) this%direction = direction
    call check_if_valid_direction(this%direction)

    if (present(aux)) this%aux = trim(aux)
    if (present(req)) this%req = trim(req)

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
  subroutine write_axis_metadata(this, fileobj, sub_axis_id)
    class(fmsDiagFullAxis_type), target, INTENT(IN) :: this     !< diag_axis obj
    class(FmsNetcdfFile_t),    INTENT(INOUT) :: fileobj     !< Fms2_io fileobj to write the data to
    integer, OPTIONAL,         INTENT(IN)    :: sub_axis_id !< ID of the sub_axis, if it exists

    character(len=:), ALLOCATABLE :: axis_edges_name !< Name of the edges, if it exist
    character(len=:), pointer     :: axis_name       !< Name of the axis
    integer                       :: axis_length     !< Size of the axis
    integer                       :: i               !< For do loops

    if (present(sub_axis_id)) then
      axis_name  => this%subaxis(sub_axis_id)%subaxis_name
      axis_length = this%subaxis(sub_axis_id)%ending_index - this%subaxis(sub_axis_id)%starting_index + 1
    else
      axis_name => this%axis_name
      axis_length = this%length
    endif

    !< Add the axis as a dimension in the netcdf file based on the type of axis_domain and the fileobj type
    select type (fileobj)
      type is (FmsNetcdfFile_t)
        !< Here the axis is not domain decomposed (i.e z_axis)
        call register_axis(fileobj, axis_name, axis_length)
      type is (FmsNetcdfDomainFile_t)
        select case (this%type_of_domain)
        case (NO_DOMAIN)
          !< Here the fileobj is domain decomposed, but the axis is not
          !! Domain decomposed fileobjs can have axis that are not domain decomposed (i.e "Z" axis)
          call register_axis(fileobj, axis_name, axis_length)
        case (TWO_D_DOMAIN)
          !< Here the axis is domain decomposed
          call register_axis(fileobj, axis_name, this%cart_name, domain_position=this%domain_position)
        end select
      type is (FmsNetcdfUnstructuredDomainFile_t)
        select case (this%type_of_domain)
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
    call register_field(fileobj, axis_name, this%type_of_data, (/axis_name/))
    call register_variable_attribute(fileobj, axis_name, "longname", this%long_name, &
      str_len=len_trim(this%long_name))

    if (this%cart_name .NE. "N") &
      call register_variable_attribute(fileobj, axis_name, "axis", this%cart_name, str_len=1)

    if (trim(this%units) .NE. "none") &
      call register_variable_attribute(fileobj, axis_name, "units", this%units, str_len=len_trim(this%units))

    select case (this%direction)
    case (direction_up)
      call register_variable_attribute(fileobj, axis_name, "positive", "up", str_len=2)
    case (direction_down)
      call register_variable_attribute(fileobj, axis_name, "positive", "down", str_len=4)
    end select

    if (allocated(this%edges_name)) then
      call register_variable_attribute(fileobj, axis_name, "edges", this%edges_name, &
        str_len=len_trim(this%edges_name))
    endif

    if(allocated(this%attributes)) then
      do i = 1, size(this%attributes)
        call register_variable_attribute(fileobj, axis_name, this%attributes(i)%att_name, &
          & this%attributes(i)%att_value)
      enddo
    endif

  end subroutine write_axis_metadata

  !> @brief Write the axis data to an open fileobj
  subroutine write_axis_data(this, fileobj, sub_axis_id)
    class(fmsDiagFullAxis_type),INTENT(IN):: this       !< diag_axis obj
    class(FmsNetcdfFile_t), INTENT(INOUT) :: fileobj   !< Fms2_io fileobj to write the data to
    integer, OPTIONAL,      INTENT(IN)    :: sub_axis_id !< ID of the sub_axis, if it exists

    integer                       :: i         !< Starting index of a sub_axis
    integer                       :: j         !< Ending index of a sub_axis

    if (present(sub_axis_id)) then
      i = this%subaxis(sub_axis_id)%starting_index
      j = this%subaxis(sub_axis_id)%ending_index

      call write_data(fileobj, this%subaxis(sub_axis_id)%subaxis_name, this%axis_data(i:j))
    else
      call write_data(fileobj, this%axis_name, this%axis_data)
    endif
  end subroutine write_axis_data

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

  !> @brief Get the name of the axis
  !> @return axis name
  pure function get_axis_name(this) &
  result (axis_name)
    class(fmsDiagFullAxis_type), intent(in)    :: this !< diag_axis obj
    CHARACTER(len=:),   ALLOCATABLE            :: axis_name

    axis_name = this%axis_name
  end function

  !> @brief Set the name of the edges
  subroutine set_edges_name(this, edges_name)
    class(fmsDiagFullAxis_type), intent(inout) :: this !< diag_axis obj
    CHARACTER(len=*),        intent(in)        :: edges_name !< Name of the edges

    this%edges_name = edges_name
  end subroutine

  !> @brief Set the subaxis of the axis obj
  !> @return A sub_axis id corresponding to the indices of the sub_axes in the sub_axes_objs array
  function set_subaxis(this, bounds) &
  result(sub_axes_id)
    class(fmsDiagFullAxis_type),  INTENT(INOUT) :: this      !< diag_axis obj
    class(*),                 INTENT(INOUT) :: bounds(:) !< bound of the subaxis

    integer :: sub_axes_id

    integer :: i !< For do loops

    !< Check if the subaxis for this bouds already exists
    do i = 1, this%nsubaxis
      if (this%subaxis(i)%exists(bounds)) return
    enddo

    !< TO DO: everything
    this%nsubaxis = this%nsubaxis + 1
    sub_axes_id = -999
  end function

  !!!!!!!!!!!!!!!!!! SUB AXIS PROCEDURES !!!!!!!!!!!!!!!!!
  !> @brief Check if a subaxis was already defined
  !> @return Flag indicating if a subaxis is already defined
  pure function check_if_subaxis_exists(this, bounds) &
  result(exists)
    class(fmsDiagSubAxis_type), INTENT(IN) :: this      !< diag_axis obj
    class(*),                   INTENT(IN)    :: bounds(:) !< bounds of the subaxis
    logical                                   :: exists

    !< TO DO: compare bounds
    exists = .false.
  end function check_if_subaxis_exists

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

end module fms_diag_axis_object_mod
!> @}
! close documentation grouping
