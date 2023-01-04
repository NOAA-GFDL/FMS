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
                              MAX_SUBAXES, DIAG_NULL, index_gridtype, latlon_gridtype
  use mpp_mod,         only:  FATAL, mpp_error, uppercase, mpp_pe, mpp_root_pe, stdout
  use fms2_io_mod,     only:  FmsNetcdfFile_t, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t, &
                            & register_axis, register_field, register_variable_attribute, write_data
  use fms_diag_yaml_mod, only: subRegion_type
  use diag_grid_mod,       only:  get_local_indices_cubesphere => get_local_indexes
  use axis_utils2_mod,   only: nearest_index
  implicit none

  PRIVATE

  public :: fmsDiagAxis_type, fms_diag_axis_object_init, fms_diag_axis_object_end, &
          & get_domain_and_domain_type, diagDomain_t, &
          & DIAGDOMAIN2D_T, fmsDiagSubAxis_type, fmsDiagAxisContainer_type, fmsDiagFullAxis_type, DIAGDOMAINUG_T
  public :: define_new_axis, define_subaxis

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

  !> @brief Type to hold the diag_axis (either subaxis or a full axis)
  !> @ingroup diag_axis_object_mod
  type :: fmsDiagAxisContainer_type
    class(fmsDiagAxis_type), allocatable :: axis
  end type

  !> @brief Type to hold the diagnostic axis description.
  !> @ingroup diag_axis_object_mod
  TYPE fmsDiagAxis_type
     INTEGER                        , private :: axis_id         !< ID of the axis

     contains
       procedure :: get_parent_axis_id
       procedure :: get_subaxes_id
       procedure :: get_axis_name
       procedure :: write_axis_metadata
       procedure :: write_axis_data
  END TYPE fmsDiagAxis_type

  !> @brief Type to hold the subaxis
  !> @ingroup diag_axis_object_mod
  TYPE, extends(fmsDiagAxis_type) :: fmsDiagSubAxis_type
    CHARACTER(len=:), ALLOCATABLE, private  :: subaxis_name   !< Name of the subaxis
    INTEGER                      , private  :: starting_index !< Starting index of the subaxis relative to the
                                                              !! parent axis
    INTEGER                      , private  :: ending_index   !< Ending index of the subaxis relative to the
                                                              !! parent axis
    type(subRegion_type)         , private  :: subRegion      !< Bounds of the subaxis (lat/lon or indices)
    INTEGER                      , private  :: parent_axis_id !< Id of the parent_axis
    contains
      procedure :: fill_subaxis
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
     integer                        , private :: subaxis(MAX_SUBAXES) !< Array of subaxis
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
     PROCEDURE :: set_edges_name
     PROCEDURE :: set_axis_id
     PROCEDURE :: get_compute_domain
     PROCEDURE :: get_indices
     PROCEDURE :: get_global_io_domain
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
  subroutine write_axis_metadata(this, fileobj, parent_axis)
    class(fmsDiagAxis_type),          target,  INTENT(IN)    :: this        !< diag_axis obj
    class(FmsNetcdfFile_t),                    INTENT(INOUT) :: fileobj     !< Fms2_io fileobj to write the data to
    class(fmsDiagAxis_type), OPTIONAL, target, INTENT(IN)    :: parent_axis !< If the axis is a subaxis, axis object
                                                                            !! for the parent axis (this will be used
                                                                            !! to get some of the metadata info)

    character(len=:), ALLOCATABLE         :: axis_edges_name !< Name of the edges, if it exist
    character(len=:), pointer             :: axis_name       !< Name of the axis
    integer                               :: axis_length     !< Size of the axis
    integer                               :: i               !< For do loops
    type(fmsDiagFullAxis_type), pointer   :: diag_axis       !< Local pointer to the diag_axis

    select type(this)
    type is (fmsDiagFullAxis_type)
      axis_name => this%axis_name
      axis_length = this%length
      diag_axis => this
    type is (fmsDiagSubAxis_type)
      axis_name => this%subaxis_name
      axis_length = this%ending_index - this%starting_index + 1
      !< Get all the other information from the parent axis (i.e the cart_name, units, etc)
      if (present(parent_axis)) then
        select type(parent_axis)
        type is (fmsDiagFullAxis_type)
          diag_axis => parent_axis
        end select
      endif
    end select

    !< Add the axis as a dimension in the netcdf file based on the type of axis_domain and the fileobj type
    select type (fileobj)
      type is (FmsNetcdfFile_t)
        !< Here the axis is not domain decomposed (i.e z_axis)
        call register_axis(fileobj, axis_name, axis_length)
        call register_field(fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
      type is (FmsNetcdfDomainFile_t)
        select case (diag_axis%type_of_domain)
        case (NO_DOMAIN)
          !< Here the fileobj is domain decomposed, but the axis is not
          !! Domain decomposed fileobjs can have axis that are not domain decomposed (i.e "Z" axis)
          call register_axis(fileobj, axis_name, axis_length)
          call register_field(fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
        case (TWO_D_DOMAIN)
          !< Here the axis is domain decomposed
          call register_axis(fileobj, axis_name, diag_axis%cart_name, domain_position=diag_axis%domain_position)
          call register_field(fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
        end select
      type is (FmsNetcdfUnstructuredDomainFile_t)
        select case (diag_axis%type_of_domain)
        case (NO_DOMAIN)
          !< Here the fileobj is in the unstructured domain, but the axis is not
          !< Unstructured domain fileobjs can have axis that are not domain decomposed (i.e "Z" axis)
          call register_axis(fileobj, axis_name, axis_length)
          call register_field(fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
        case (UG_DOMAIN)
          !< Here the axis is in a unstructured domain
          call register_axis(fileobj, axis_name)
          call register_field(fileobj, axis_name, diag_axis%type_of_data, (/axis_name/))
        end select
    end select

    !< Write its metadata
    call register_variable_attribute(fileobj, axis_name, "long_name", diag_axis%long_name, &
      str_len=len_trim(diag_axis%long_name))

    if (diag_axis%cart_name .NE. "N") &
      call register_variable_attribute(fileobj, axis_name, "axis", diag_axis%cart_name, str_len=1)

    if (trim(diag_axis%units) .NE. "none") &
      call register_variable_attribute(fileobj, axis_name, "units", diag_axis%units, str_len=len_trim(diag_axis%units))

    select case (diag_axis%direction)
    case (direction_up)
      call register_variable_attribute(fileobj, axis_name, "positive", "up", str_len=2)
    case (direction_down)
      call register_variable_attribute(fileobj, axis_name, "positive", "down", str_len=4)
    end select

    if (allocated(diag_axis%edges_name)) then
      call register_variable_attribute(fileobj, axis_name, "edges", diag_axis%edges_name, &
        str_len=len_trim(diag_axis%edges_name))
    endif

    if(allocated(diag_axis%attributes)) then
      do i = 1, diag_axis%num_attributes
        select type (att_value => diag_axis%attributes(i)%att_value)
        type is (character(len=*))
          call register_variable_attribute(fileobj, axis_name, diag_axis%attributes(i)%att_name, trim(att_value(1)), &
                                           str_len=len_trim(att_value(1)))
        class default
          call register_variable_attribute(fileobj, axis_name, diag_axis%attributes(i)%att_name, att_value)
        end select
      enddo
    endif

  end subroutine write_axis_metadata

  !> @brief Write the axis data to an open fileobj
  subroutine write_axis_data(this, fileobj, parent_axis)
    class(fmsDiagAxis_type),           target, INTENT(IN)    :: this        !< diag_axis obj
    class(FmsNetcdfFile_t),                    INTENT(INOUT) :: fileobj     !< Fms2_io fileobj to write the data to
    class(fmsDiagAxis_type), OPTIONAL, target, INTENT(IN)    :: parent_axis !< The parent axis if this is a subaxis

    integer                       :: i                 !< Starting index of a sub_axis
    integer                       :: j                 !< Ending index of a sub_axis
    integer                       :: global_io_index(2)!< Global io domain starting and ending index
    select type(this)
    type is (fmsDiagFullAxis_type)
      call this%get_global_io_domain(global_io_index)
      call write_data(fileobj, this%axis_name, this%axis_data(global_io_index(1):global_io_index(2)))
    type is (fmsDiagSubAxis_type)
      i = this%starting_index
      j = this%ending_index

      if (present(parent_axis)) then
        select type(parent_axis)
        type is (fmsDiagFullAxis_type)
          call write_data(fileobj, this%subaxis_name, parent_axis%axis_data(i:j))
        end select
      endif
    end select
  end subroutine write_axis_data

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

  !> @brief Set the axis_id
  subroutine set_axis_id(this, axis_id)
    class(fmsDiagFullAxis_type), intent(inout) :: this    !< diag_axis obj
    integer,                     intent(in)    :: axis_id !< Axis_id

    this%axis_id = axis_id

  end subroutine set_axis_id

  !> @brief Set the name of the edges
  subroutine set_edges_name(this, edges_name)
    class(fmsDiagFullAxis_type), intent(inout) :: this !< diag_axis obj
    CHARACTER(len=*),        intent(in)        :: edges_name !< Name of the edges

    this%edges_name = edges_name
  end subroutine

  !> @brief Determine if the subRegion is in the current PE.
  !! If it is, determine the starting and ending indices of the current PE that belong to the subRegion
  subroutine get_indices(this, compute_idx, corners_indices, starting_index, ending_index, need_to_define_axis)
    class(fmsDiagFullAxis_type), intent(inout) :: this                !< diag_axis obj
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
    if (compute_idx(1) > subregion_start .and. compute_idx(2) > subregion_start) return
    if (compute_idx(1) > subregion_end   .and. compute_idx(2) > subregion_end) return

    need_to_define_axis = .true.
    if (compute_idx(1)  >= subregion_start .and. compute_idx(2) >= subregion_end) then
      !< In this case all the point of the current PE are inside the range of the sub_axis
      starting_index = compute_idx(1)
      ending_index   = compute_idx(2)
    else if (compute_idx(1)  >= subregion_start .and. compute_idx(2) <= subregion_end) then
      !< In this case all the points of the current PE are valid up to the end point
      starting_index = compute_idx(1)
      ending_index   = subregion_end
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
  subroutine fill_subaxis(this, starting_index, ending_index, axis_id, parent_id, parent_axis_name, subRegion)
    class(fmsDiagSubAxis_type), INTENT(INOUT) :: this             !< diag_sub_axis obj
    integer                   , intent(in)    :: starting_index   !< Starting index of the subRegion for the PE
    integer                   , intent(in)    :: ending_index     !< Ending index of the subRegion for the PE
    integer                   , intent(in)    :: axis_id          !< Axis id to assign to the subaxis
    integer                   , intent(in)    :: parent_id        !< The id of the parent axis, the subaxis belongs to
    type(subRegion_type)      , intent(in)    :: subRegion        !< SubRegion definition as it is defined in the yaml
    character(len=*)          , intent(in)    :: parent_axis_name !< Name of the parent_axis

    this%axis_id = axis_id
    this%starting_index = starting_index
    this%ending_index = ending_index
    this%parent_axis_id = parent_id
    this%subRegion = subRegion
    this%subaxis_name = trim(parent_axis_name)//"_sub01"
  end subroutine fill_subaxis

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

  !< @brief Determinet the axis name of an axis_object
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

  !> @brief Define a subaxis based on the subRegion defined by the yaml
  subroutine define_subaxis (diag_axis, axis_ids, naxis, subRegion, is_cube_sphere, write_on_this_pe)
    class(fmsDiagAxisContainer_type), target, intent(inout) :: diag_axis(:)     !< Diag_axis object
    integer,                                  INTENT(in)    :: axis_ids(:)      !< Array of axes_ids
    integer,                                  intent(inout) :: naxis            !< Number of axis registered
    type(subRegion_type),                     intent(in)    :: subRegion        !< The subRegion definition from
                                                                                !! the yaml
    logical,                                  intent(in)    :: is_cube_sphere   !< .true. if this is a cubesphere
    logical,                                  intent(out)   :: write_on_this_pe !< .true. if the subregion
                                                                                !! is on this PE

    select case(subRegion%grid_type)
    case (latlon_gridtype)
      call define_subaxis_latlon(diag_axis, axis_ids, naxis, subRegion, is_cube_sphere, write_on_this_pe)
    case (index_gridtype)
      call define_subaxis_index(diag_axis, axis_ids, naxis, subRegion, write_on_this_pe)
    end select
  end subroutine define_subaxis

  !> @brief Fill in the subaxis object for a subRegion defined by index
  subroutine define_subaxis_index(diag_axis, axis_ids, naxis, subRegion, write_on_this_pe)
    class(fmsDiagAxisContainer_type), target, intent(inout) :: diag_axis(:)     !< Diag_axis object
    integer,                                  INTENT(in)    :: axis_ids(:)      !< Array of axes_ids
    integer,                                  intent(inout) :: naxis            !< Number of axis registered
    type(subRegion_type),                     intent(in)    :: subRegion        !< SubRegion definition from the yaml
    logical,                                  intent(out)   :: write_on_this_pe !< .true. if the subregion
                                                                                !! is on this PE
    integer :: i !< For do loops
    integer :: compute_idx(2)
    integer :: starting_index, ending_index
    logical :: need_to_define_axis
    integer :: lat_indices(2), lon_indices(2)


    do i = 1, size(axis_ids)
      select type (parent_axis => diag_axis(axis_ids(i))%axis)
      type is (fmsDiagFullAxis_type)
        !< Get the PEs compute domain
        call parent_axis%get_compute_domain(compute_idx, need_to_define_axis, tile_number=subRegion%tile)

        !< If this is not a "X" or "Y" axis, go to the next axis
        if (.not. need_to_define_axis) then
          cycle
        endif

        !< Determine if the PE's compute domain is inside the subRegion
        !! If it is get the starting and ending indices for that PE
        call parent_axis%get_indices(compute_idx, subRegion%corners(:,i), starting_index, ending_index, &
          need_to_define_axis)

        !< If the PE's compute is not inside the subRegion, define a null subaxis and go to the next axis
        if (.not. need_to_define_axis) then
          call define_new_axis(diag_axis, parent_axis, naxis, axis_ids(i), &
            subRegion, diag_null, diag_null)
          cycle
        endif

        !< If it made it to this point, the current PE is in the subRegion!
        write_on_this_pe = .true.

        call define_new_axis(diag_axis, parent_axis, naxis, axis_ids(i), &
            subRegion, starting_index, ending_index)
        end select
    enddo

  end subroutine define_subaxis_index

  !> @brief Fill in the subaxis object for a subRegion defined by lat lon
  subroutine define_subaxis_latlon(diag_axis, axis_ids, naxis, subRegion, is_cube_sphere, write_on_this_pe)
    class(fmsDiagAxisContainer_type), target, intent(inout) :: diag_axis(:)     !< Diag_axis object
    integer,                                  INTENT(in)    :: axis_ids(:)      !< Array of axes_ids
    integer,                                  intent(inout) :: naxis            !< Number of axis registered
    type(subRegion_type),                     intent(in)    :: subRegion        !< SubRegion definition from the yaml
    logical,                                  intent(in)    :: is_cube_sphere   !< .true. if this is a cubesphere
    logical,                                  intent(out)   :: write_on_this_pe !< .true. if the subregion
                                                                                !! is on this PE

    real    :: lat(2)              !< Starting and ending lattiude of the subRegion
    real    :: lon(2)              !< Starting and ending longitude or the subRegion
    integer :: lat_indices(2)      !< Starting and ending latitude indices of the subRegion
    integer :: lon_indices(2)      !< Starting and ending longitude indices of the subRegion
    integer :: compute_idx(2)      !< Compute domain of the current axis
    integer :: starting_index      !< Starting index of the subRegion for the current PE
    integer :: ending_index        !< Ending index of the subRegion for the current PE
    logical :: need_to_define_axis !< .true. if it is needed to define the subaxis
    integer :: i                   !< For do loops

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
          call parent_axis%get_compute_domain(compute_idx, need_to_define_axis)

          !< If this is not a "X" or "Y" axis go to the next axis
          if (.not. need_to_define_axis) cycle

          !< Determine if the PE's compute domain is inside the subRegion
          !! If it is get the starting and ending indices for that PE
          if (parent_axis%cart_name .eq. "X") then
            call parent_axis%get_indices(compute_idx, lon_indices, starting_index, ending_index, &
              need_to_define_axis)
          else if (parent_axis%cart_name .eq. "Y") then
            call parent_axis%get_indices(compute_idx, lat_indices, starting_index, ending_index, &
              need_to_define_axis)
          endif

          !< If the PE's compute is not inside the subRegion move to the next axis
          if (.not. need_to_define_axis) cycle

          !< If it made it to this point, the current PE is in the subRegion!
          write_on_this_pe = .true.

          call define_new_axis(diag_axis, parent_axis, naxis, axis_ids(i), &
            subRegion, starting_index, ending_index)
        end select select_axis_type
      enddo loop_over_axis_ids
    else if_is_cube_sphere
      loop_over_axis_ids2: do i = 1, size(axis_ids)
        select type (parent_axis => diag_axis(axis_ids(i))%axis)
        type is (fmsDiagFullAxis_type)
          !< Get the PEs compute domain
          call parent_axis%get_compute_domain(compute_idx, need_to_define_axis)

          !< If this is not a "X" or "Y" axis go to the next axis
          if (.not. need_to_define_axis) cycle

          !< Get the starting and ending indices of the subregion relative to the global grid
          if (parent_axis%cart_name .eq. "X") then
            select type(adata=>parent_axis%axis_data)
            type is (real)
              lon_indices(1) = nearest_index(lon(1), adata)
              lon_indices(2) = nearest_index(lon(2), adata) + 1
            end select
            call parent_axis%get_indices(compute_idx, lon_indices, starting_index, ending_index, &
              need_to_define_axis)
          else if (parent_axis%cart_name .eq. "Y") then
            select type(adata=>parent_axis%axis_data)
            type is (real)
              lat_indices(1) = nearest_index(lat(1), adata)
              lat_indices(2) = nearest_index(lat(2), adata) + 1
            end select
            call parent_axis%get_indices(compute_idx, lat_indices, starting_index, ending_index, &
              need_to_define_axis)
          endif

          !< If the PE's compute is not inside the subRegion move to the next axis
          if (.not. need_to_define_axis) cycle

          !< If it made it to this point, the current PE is in the subRegion!
          write_on_this_pe = .true.

          call define_new_axis(diag_axis, parent_axis, naxis, axis_ids(i), &
            subRegion, starting_index, ending_index)
        end select
      enddo loop_over_axis_ids2
    endif if_is_cube_sphere
  end subroutine define_subaxis_latlon

  !< Creates a new subaxis and fills it will all the information it needs
  subroutine define_new_axis(diag_axis, parent_axis, naxis, parent_id, subRegion, &
                             starting_index, ending_index)

    class(fmsDiagAxisContainer_type), target, intent(inout) :: diag_axis(:)     !< Diag_axis object
    class(fmsDiagFullAxis_type),              intent(inout) :: parent_axis      !< The parent axis
    integer,                                  intent(inout) :: naxis            !< The number of axis that
                                                                                !! have been defined
    integer,                                  intent(in)    :: parent_id        !< Id of the parent axis
    type(subRegion_type),                     intent(in)    :: subRegion        !< SubRegion definition from the yaml
    integer,                                  intent(in)    :: starting_index   !< PE's Starting index
    integer,                                  intent(in)    :: ending_index     !< PE's Ending index

    naxis = naxis + 1 !< This is the axis id of the new axis!

    !< Add the axis_id of the new subaxis to the parent axis
    parent_axis%nsubaxis = parent_axis%nsubaxis + 1
    parent_axis%subaxis(parent_axis%nsubaxis) = naxis

    !< Allocate the new axis as a subaxis and fill it
    allocate(fmsDiagSubAxis_type :: diag_axis(naxis)%axis)
    diag_axis(naxis)%axis%axis_id = naxis

    select type (sub_axis => diag_axis(naxis)%axis)
    type is (fmsDiagSubAxis_type)
      call sub_axis%fill_subaxis(starting_index, ending_index, naxis, parent_id, &
             parent_axis%axis_name, subRegion)
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

#endif
end module fms_diag_axis_object_mod
!> @}
! close documentation grouping
