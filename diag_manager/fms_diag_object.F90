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
module fms_diag_object_mod
use mpp_mod, only: fatal, note, warning, mpp_error, mpp_pe, mpp_root_pe, stdout
use diag_data_mod,  only: diag_null, diag_not_found, diag_not_registered, diag_registered_id, &
                         &DIAG_FIELD_NOT_FOUND, diag_not_registered, max_axes, TWO_D_DOMAIN
  USE time_manager_mod, ONLY: set_time, set_date, get_time, time_type, OPERATOR(>=), OPERATOR(>),&
       & OPERATOR(<), OPERATOR(==), OPERATOR(/=), OPERATOR(/), OPERATOR(+), ASSIGNMENT(=), get_date, &
       & get_ticks_per_second
#ifdef use_yaml
use fms_diag_file_object_mod, only: fmsDiagFileContainer_type, fmsDiagFile_type, fms_diag_files_object_init
use fms_diag_field_object_mod, only: fmsDiagField_type, fms_diag_fields_object_init
use fms_diag_yaml_mod, only: diag_yaml_object_init, diag_yaml_object_end, find_diag_field, &
                            & get_diag_files_id, diag_yaml
use fms_diag_axis_object_mod, only: fms_diag_axis_object_init, fmsDiagAxis_type, fmsDiagSubAxis_type, &
                                   &diagDomain_t, get_domain_and_domain_type, diagDomain2d_t, &
                                   &fmsDiagAxisContainer_type, fms_diag_axis_object_end, fmsDiagFullAxis_type
use fms_diag_buffer_mod
#endif
use mpp_domains_mod, only: domain1d, domain2d, domainUG, null_domain2d
implicit none
private

type fmsDiagObject_type
!TODO add container arrays
#ifdef use_yaml
private
!TODO: Remove FMS prefix from variables in this type
  class(fmsDiagFileContainer_type), allocatable :: FMS_diag_files (:) !< array of diag files
  class(fmsDiagField_type), allocatable :: FMS_diag_fields(:) !< Array of diag fields
  type(fmsDiagBufferContainer_type), allocatable :: FMS_diag_buffers(:) !< array of buffer objects
  integer, private :: registered_buffers = 0 !< number of registered buffers, per dimension
  class(fmsDiagAxisContainer_type), allocatable :: diag_axis(:) !< Array of diag_axis
  integer, private :: registered_variables !< Number of registered variables
  integer, private :: registered_axis !< Number of registered axis
  logical, private :: initialized=.false. !< True if the fmsDiagObject is initialized
  logical, private :: files_initialized=.false. !< True if the fmsDiagObject is initialized
  logical, private :: fields_initialized=.false. !< True if the fmsDiagObject is initialized
  logical, private :: buffers_initialized=.false. !< True if the fmsDiagObject is initialized
  logical, private :: axes_initialized=.false. !< True if the fmsDiagObject is initialized
#endif
  contains
    procedure :: init => fms_diag_object_init
    procedure :: diag_end => fms_diag_object_end
    procedure :: fms_register_diag_field_scalar
    procedure :: fms_register_diag_field_array
    procedure :: fms_register_static_field
    procedure :: fms_diag_axis_init
    procedure :: register => fms_register_diag_field_obj !! Merely initialize fields.
    procedure :: fms_diag_field_add_attribute
    procedure :: fms_diag_axis_add_attribute
    procedure :: fms_get_domain2d
    procedure :: fms_get_axis_length
    procedure :: fms_get_diag_field_id_from_name
    procedure :: fms_get_axis_name_from_id
    procedure :: fms_diag_send_complete
    procedure :: fms_diag_do_io
#ifdef use_yaml
    procedure :: get_diag_buffer
#endif
end type fmsDiagObject_type

type (fmsDiagObject_type), target :: fms_diag_object

public :: fms_register_diag_field_obj
public :: fms_register_diag_field_scalar
public :: fms_register_diag_field_array
public :: fms_register_static_field
public :: fms_diag_field_add_attribute
public :: fms_get_diag_field_id_from_name
public :: fms_diag_object
public :: fmsDiagObject_type
integer, private :: registered_variables !< Number of registered variables
public :: dump_diag_obj

contains

!> @brief Initiliazes the  fms_diag_object.
!! Reads the diag_table.yaml and fills in the yaml object
!! Allocates the diag manager object arrays for files, fields, and buffers
!! Initializes variables
subroutine fms_diag_object_init (this,diag_subset_output)
 class(fmsDiagObject_type) :: this !< Diag mediator/controller object
 integer :: diag_subset_output !< Subset of the diag output?
#ifdef use_yaml
 if (this%initialized) return

! allocate(diag_objs(get_num_unique_fields()))
  CALL diag_yaml_object_init(diag_subset_output)
  this%axes_initialized = fms_diag_axis_object_init(this%diag_axis)
  this%files_initialized = fms_diag_files_object_init(this%FMS_diag_files)
  this%fields_initialized = fms_diag_fields_object_init(this%FMS_diag_fields)
  this%buffers_initialized = fms_diag_buffer_init(this%FMS_diag_buffers, SIZE(diag_yaml%get_diag_fields()))
  this%registered_variables = 0
  this%registered_axis = 0
  this%initialized = .true.
#else
  call mpp_error("fms_diag_object_init",&
    "You must compile with -Duse_yaml to use the option use_modern_diag", FATAL)
#endif
end subroutine fms_diag_object_init

!> \description Loops through all files and does one final write.
!! Closes all files
!! Deallocates all buffers, fields, and files
!! Uninitializes the fms_diag_object
subroutine fms_diag_object_end (this, time)
  class(fmsDiagObject_type) :: this
  TYPE(time_type), INTENT(in) :: time

  integer                   :: i
#ifdef use_yaml
  !TODO: loop through files and force write
  if (.not. this%initialized) return

  call this%fms_diag_do_io(time, is_end_of_run=.true.)
  !TODO: Deallocate diag object arrays and clean up all memory
  do i=1, size(this%FMS_diag_buffers)
    if(allocated(this%FMS_diag_buffers(i)%diag_buffer_obj)) then
      call this%FMS_diag_buffers(i)%diag_buffer_obj%flush_buffer()
    endif
  enddo
  deallocate(this%FMS_diag_buffers)
  this%axes_initialized = fms_diag_axis_object_end(this%diag_axis)
  this%initialized = .false.
  call diag_yaml_object_end
#else
  call mpp_error(FATAL, "You can not call fms_diag_object%end without yaml")
#endif
end subroutine fms_diag_object_end

!> @brief Registers a field.
!! @description This to avoid having duplicate code in each of the _scalar, _array and _static register calls
!! @return field index for subsequent call to send_data.
integer function fms_register_diag_field_obj &
       (this, modname, varname, axes, init_time, &
       longname, units, missing_value, varRange, mask_variant, standname, &
       do_not_log, err_msg, interp_method, tile_count, area, volume, realm, static)

 class(fmsDiagObject_type),TARGET,INTENT(inout):: this       !< Diaj_obj to fill
 CHARACTER(len=*),               INTENT(in)    :: modname               !< The module name
 CHARACTER(len=*),               INTENT(in)    :: varname               !< The variable name
 TYPE(time_type),  OPTIONAL,     INTENT(in)    :: init_time             !< Initial time
 INTEGER, TARGET,  OPTIONAL,     INTENT(in)    :: axes(:)               !< The axes indicies
 CHARACTER(len=*), OPTIONAL,     INTENT(in)    :: longname              !< THe variables long name
 CHARACTER(len=*), OPTIONAL,     INTENT(in)    :: units                 !< The units of the variables
 CHARACTER(len=*), OPTIONAL,     INTENT(in)    :: standname             !< The variables stanard name
 class(*),         OPTIONAL,     INTENT(in)    :: missing_value         !< Missing value to add as a attribute
 class(*),         OPTIONAL,     INTENT(in)    :: varRANGE(2)           !< Range to add as a attribute
 LOGICAL,          OPTIONAL,     INTENT(in)    :: mask_variant          !< Mask
 LOGICAL,          OPTIONAL,     INTENT(in)    :: do_not_log            !< if TRUE, field info is not logged
 CHARACTER(len=*), OPTIONAL,     INTENT(out)   :: err_msg               !< Error message to be passed back up
 CHARACTER(len=*), OPTIONAL,     INTENT(in)    :: interp_method         !< The interp method to be used when
                                                                        !! regridding the field in post-processing.
                                                                        !! Valid options are "conserve_order1",
                                                                        !! "conserve_order2", and "none".
 INTEGER,          OPTIONAL,     INTENT(in)    :: tile_count            !< the number of tiles
 INTEGER,          OPTIONAL,     INTENT(in)    :: area                  !< diag_field_id of the cell area field
 INTEGER,          OPTIONAL,     INTENT(in)    :: volume                !< diag_field_id of the cell volume field
 CHARACTER(len=*), OPTIONAL,     INTENT(in)    :: realm                 !< String to set as the value to the
                                                                        !! modeling_realm attribute
 LOGICAL,          OPTIONAL,     INTENT(in)    :: static                !< True if the variable is static
#ifdef use_yaml

 class (fmsDiagFile_type), pointer :: fileptr => null() !< Pointer to the diag_file
 class (fmsDiagField_type), pointer :: fieldptr => null() !< Pointer to the diag_field
 integer, allocatable :: file_ids(:) !< The file IDs for this variable
 integer :: i !< For do loops
 integer, allocatable :: diag_field_indices(:) !< indices where the field was found in the yaml
#endif
#ifndef use_yaml
fms_register_diag_field_obj = DIAG_FIELD_NOT_FOUND
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
 diag_field_indices = find_diag_field(varname, modname)
 if (diag_field_indices(1) .eq. diag_null) then
    !< The field was not found in the table, so return diag_null
    fms_register_diag_field_obj = diag_null
    deallocate(diag_field_indices)
    return
  endif

  this%registered_variables = this%registered_variables + 1
  fms_register_diag_field_obj = this%registered_variables

  call this%FMS_diag_fields(this%registered_variables)%&
    &setID(this%registered_variables)

!> Use pointers for convenience
  fieldptr => this%FMS_diag_fields(this%registered_variables)
!> Register the data for the field
  call fieldptr%register(modname, varname, diag_field_indices, fms_diag_object%diag_axis, &
       axes=axes, longname=longname, units=units, missing_value=missing_value, varRange= varRange, &
       mask_variant= mask_variant, standname=standname, do_not_log=do_not_log, err_msg=err_msg, &
       interp_method=interp_method, tile_count=tile_count, area=area, volume=volume, realm=realm, &
       static=static)
!> Get the file IDs from the field indicies from the yaml
  file_ids = get_diag_files_id(diag_field_indices)
  call fieldptr%set_file_ids(file_ids)
!> Add the axis information, initial time, and field IDs to the files
  if (present(axes) .and. present(init_time)) then
    do i = 1, size(file_ids)
     fileptr => this%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_and_yaml_id(fieldptr%get_id(), diag_field_indices(i))
     call fileptr%set_file_domain(fieldptr%get_domain(), fieldptr%get_type_of_domain())
     call fileptr%add_axes(axes, this%diag_axis, this%registered_axis)
     call fileptr%add_start_time(init_time)
     call fileptr%set_file_time_ops (fieldptr%diag_field(i), fieldptr%is_static())
    enddo
  elseif (present(axes)) then !only axes present
    do i = 1, size(file_ids)
     fileptr => this%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_and_yaml_id(fieldptr%get_id(), diag_field_indices(i))
     call fileptr%set_file_domain(fieldptr%get_domain(), fieldptr%get_type_of_domain())
     call fileptr%add_axes(axes, this%diag_axis, this%registered_axis)
     call fileptr%set_file_time_ops (fieldptr%diag_field(i), fieldptr%is_static())
    enddo
  elseif (present(init_time)) then !only inti time present
    do i = 1, size(file_ids)
     fileptr => this%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_and_yaml_id(fieldptr%get_id(), diag_field_indices(i))
     call fileptr%add_start_time(init_time)
     call fileptr%set_file_time_ops (fieldptr%diag_field(i), fieldptr%is_static())
    enddo
  else !no axis or init time present
    do i = 1, size(file_ids)
     fileptr => this%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_and_yaml_id(fieldptr%get_id(), diag_field_indices(i))
     call fileptr%set_file_time_ops (fieldptr%diag_field(i), fieldptr%is_static())
    enddo
  endif
  nullify (fileptr)
  nullify (fieldptr)
  deallocate(diag_field_indices)
#endif
end function fms_register_diag_field_obj

  !> @brief Registers a scalar field
  !! @return field index for subsequent call to send_data.
INTEGER FUNCTION fms_register_diag_field_scalar(this,module_name, field_name, init_time, &
       & long_name, units, missing_value, var_range, standard_name, do_not_log, err_msg,&
       & area, volume, realm)
    class(fmsDiagObject_type),TARGET,INTENT(inout):: this       !< Diaj_obj to fill
    CHARACTER(len=*),           INTENT(in) :: module_name   !< Module where the field comes from
    CHARACTER(len=*),           INTENT(in) :: field_name    !< Name of the field
    TYPE(time_type),  OPTIONAL, INTENT(in) :: init_time     !< Time to start writing data from
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: long_name     !< Long_name to add as a variable attribute
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: units         !< Units to add as a variable_attribute
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: standard_name !< Standard_name to name the variable in the file
    CLASS(*),         OPTIONAL, INTENT(in) :: missing_value !< Missing value to add as a variable attribute
    CLASS(*),         OPTIONAL, INTENT(in) :: var_range(:)  !< Range to add a variable attribute
    LOGICAL,          OPTIONAL, INTENT(in) :: do_not_log    !< If TRUE, field information is not logged
    CHARACTER(len=*), OPTIONAL, INTENT(out):: err_msg       !< Error_msg from call
    INTEGER,          OPTIONAL, INTENT(in) :: area          !< Id of the area field
    INTEGER,          OPTIONAL, INTENT(in) :: volume        !< Id of the volume field
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: realm         !< String to set as the modeling_realm attribute
#ifndef use_yaml
fms_register_diag_field_scalar=diag_null
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
    fms_register_diag_field_scalar = this%register(&
      & module_name, field_name, init_time=init_time, &
      & longname=long_name, units=units, missing_value=missing_value, varrange=var_range, &
      & standname=standard_name, do_not_log=do_not_log, err_msg=err_msg, &
      & area=area, volume=volume, realm=realm)
#endif
end function fms_register_diag_field_scalar

    !> @brief Registers an array field
  !> @return field index for subsequent call to send_data.
INTEGER FUNCTION fms_register_diag_field_array(this, module_name, field_name, axes, init_time, &
       & long_name, units, missing_value, var_range, mask_variant, standard_name, verbose,&
       & do_not_log, err_msg, interp_method, tile_count, area, volume, realm)
    class(fmsDiagObject_type),TARGET,INTENT(inout):: this       !< Diaj_obj to fill
    CHARACTER(len=*),           INTENT(in) :: module_name   !< Module where the field comes from
    CHARACTER(len=*),           INTENT(in) :: field_name    !< Name of the field
    INTEGER,                    INTENT(in) :: axes(:)       !< Ids corresponding to the variable axis
    TYPE(time_type),  OPTIONAL, INTENT(in) :: init_time     !< Time to start writing data from
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: long_name     !< Long_name to add as a variable attribute
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: units         !< Units to add as a variable_attribute
    CLASS(*),         OPTIONAL, INTENT(in) :: missing_value !< Missing value to add as a variable attribute
    CLASS(*),         OPTIONAL, INTENT(in) :: var_range(:)  !< Range to add a variable attribute
    LOGICAL,          OPTIONAL, INTENT(in) :: mask_variant  !< Mask variant
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: standard_name !< Standard_name to name the variable in the file
    LOGICAL,          OPTIONAL, INTENT(in) :: verbose       !< Print more information
    LOGICAL,          OPTIONAL, INTENT(in) :: do_not_log    !< If TRUE, field information is not logged
    CHARACTER(len=*), OPTIONAL, INTENT(out):: err_msg       !< Error_msg from call
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: interp_method !< The interp method to be used when
                                                            !! regridding the field in post-processing.
                                                            !! Valid options are "conserve_order1",
                                                            !! "conserve_order2", and "none".
    INTEGER,          OPTIONAL, INTENT(in) :: tile_count    !< The current tile number
    INTEGER,          OPTIONAL, INTENT(in) :: area          !< Id of the area field
    INTEGER,          OPTIONAL, INTENT(in) :: volume        !< Id of the volume field
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: realm         !< String to set as the modeling_realm attribute

#ifndef use_yaml
fms_register_diag_field_array=diag_null
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
    fms_register_diag_field_array = this%register( &
      & module_name, field_name, init_time=init_time, &
      & axes=axes, longname=long_name, units=units, missing_value=missing_value, varrange=var_range, &
      & mask_variant=mask_variant, standname=standard_name, do_not_log=do_not_log, err_msg=err_msg, &
      & interp_method=interp_method, tile_count=tile_count, area=area, volume=volume, realm=realm)
#endif
end function fms_register_diag_field_array

!> @brief Return field index for subsequent call to send_data.
!! @return field index for subsequent call to send_data.
INTEGER FUNCTION fms_register_static_field(this, module_name, field_name, axes, long_name, units,&
       & missing_value, range, mask_variant, standard_name, DYNAMIC, do_not_log, interp_method,&
       & tile_count, area, volume, realm)
    class(fmsDiagObject_type),TARGET,INTENT(inout):: this       !< Diaj_obj to fill
    CHARACTER(len=*),                         INTENT(in) :: module_name   !< Name of the module, the field is on
    CHARACTER(len=*),                         INTENT(in) :: field_name    !< Name of the field
    INTEGER,          DIMENSION(:),           INTENT(in) :: axes          !< Axes_id of the field
    CHARACTER(len=*),               OPTIONAL, INTENT(in) :: long_name     !< Longname to be added as a attribute
    CHARACTER(len=*),               OPTIONAL, INTENT(in) :: units         !< Units to be added as a attribute
    CHARACTER(len=*),               OPTIONAL, INTENT(in) :: standard_name !< Standard name to be added as a attribute
    CLASS(*),                       OPTIONAL, INTENT(in) :: missing_value !< Missing value to be added as a attribute
    CLASS(*),                       OPTIONAL, INTENT(in) :: range(:)      !< Range to be added as a attribute
    LOGICAL,                        OPTIONAL, INTENT(in) :: mask_variant  !< Flag indicating if the field is has
                                                                          !! a mask variant
    LOGICAL,                        OPTIONAL, INTENT(in) :: DYNAMIC       !< Flag indicating if the field is dynamic
    LOGICAL,                        OPTIONAL, INTENT(in) :: do_not_log    !< if TRUE, field information is not logged
    CHARACTER(len=*),               OPTIONAL, INTENT(in) :: interp_method !< The interp method to be used when
                                                                          !! regridding the field in post-processing
                                                                          !! Valid options are "conserve_order1",
                                                                          !! "conserve_order2", and "none".
    INTEGER,                        OPTIONAL, INTENT(in) :: tile_count    !! Number of tiles
    INTEGER,                        OPTIONAL, INTENT(in) :: area          !< Field ID for the area field associated
                                                                          !! with this field
    INTEGER,                        OPTIONAL, INTENT(in) :: volume        !< Field ID for the volume field associated
                                                                          !! with this field
    CHARACTER(len=*),               OPTIONAL, INTENT(in) :: realm         !< String to set as the value to the
                                                                          !! modeling_realm attribute

#ifndef use_yaml
fms_register_static_field=diag_null
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
! Include static as optional variable to register here
  fms_register_static_field = this%register( &
      & module_name, field_name, axes=axes, &
      & longname=long_name, units=units, missing_value=missing_value, varrange=range, &
      & mask_variant=mask_variant, do_not_log=do_not_log, interp_method=interp_method, tile_count=tile_count, &
      & standname=standard_name, area=area, volume=volume, realm=realm, &
      & static=.true.)
#endif
end function fms_register_static_field

!> @brief Wrapper for the register_diag_axis subroutine. This is needed to keep the diag_axis_init
!! interface the same
!> @return Axis id
FUNCTION fms_diag_axis_init(this, axis_name, axis_data, units, cart_name, long_name, direction,&
  & set_name, edges, Domain, Domain2, DomainU, aux, req, tile_count, domain_position ) &
  & result(id)

  class(fmsDiagObject_type),TARGET,INTENT(inout):: this       !< Diaj_obj to fill
  CHARACTER(len=*),   INTENT(in)           :: axis_name       !< Name of the axis
  CLASS(*),           INTENT(in)           :: axis_data(:)    !< Array of coordinate values
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

#ifndef use_yaml
id = diag_null
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
  CHARACTER(len=:),   ALLOCATABLE :: edges_name !< Name of the edges

  this%registered_axis = this%registered_axis + 1

  if (this%registered_axis > max_axes) call mpp_error(FATAL, &
    &"diag_axis_init: max_axes exceeded, increase via diag_manager_nml")

  allocate(fmsDiagFullAxis_type :: this%diag_axis(this%registered_axis)%axis)

  select type (axis => this%diag_axis(this%registered_axis)%axis )
  type is (fmsDiagFullAxis_type)
    if(present(edges)) then
      if (edges < 0 .or. edges > this%registered_axis) &
        call mpp_error(FATAL, "diag_axit_init: The edge axis has not been defined. "&
                               "Call diag_axis_init for the edge axis first")
      select type (edges_axis => this%diag_axis(edges)%axis)
      type is (fmsDiagFullAxis_type)
        edges_name = edges_axis%get_axis_name()
        call axis%set_edges_name(edges_name)
      end select
    endif
    call axis%register(axis_name, axis_data, units, cart_name, long_name=long_name, &
      & direction=direction, set_name=set_name, Domain=Domain, Domain2=Domain2, DomainU=DomainU, aux=aux, &
      & req=req, tile_count=tile_count, domain_position=domain_position)

    id = this%registered_axis
    call axis%set_axis_id(id)
  end select
#endif
end function fms_diag_axis_init

!> @brief Loops through all the files, open the file, writes out axis and
!! variable metadata and data when necessary.
subroutine fms_diag_send_complete(this, time_step)
  class(fmsDiagObject_type), target, intent (inout) :: this          !< The diag object
  TYPE (time_type),                  INTENT(in)     :: time_step     !< The current model time

#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else

  call this%fms_diag_do_io(time_step)

#endif

end subroutine fms_diag_send_complete

subroutine fms_diag_do_io(this, time_step, is_end_of_run)
  class(fmsDiagObject_type), target, intent (inout) :: this          !< The diag object
  TYPE (time_type),                  INTENT(in)     :: time_step     !< The current model time
  logical,                 optional, intent(in)     :: is_end_of_run !< If .true. this is the end of the run,
                                                                     !! so force write
#ifdef use_yaml
  integer :: i !< For do loops
  class(fmsDiagFileContainer_type), pointer :: diag_file !< Pointer to this%FMS_diag_files(i) (for convenience)
  logical :: file_is_opened_this_time_step !< True if the file was opened in this time_step
                                           !! If true the metadata will need to be written
  logical :: force_write

  force_write = .false.
  if (present (is_end_of_run)) force_write = .true.

  do i = 1, size(this%FMS_diag_files)
    diag_file => this%FMS_diag_files(i)

    !< Go away if the file is a subregional file and the current PE does not have any data for it
    if (.not. diag_file%writing_on_this_pe()) cycle

    call diag_file%open_diag_file(time_step, file_is_opened_this_time_step)
    if (file_is_opened_this_time_step) then
      call diag_file%write_time_metadata()
      call diag_file%write_axis_metadata(this%diag_axis)
      call diag_file%write_field_metadata(this%FMS_diag_fields, this%diag_axis)
      call diag_file%write_axis_data(this%diag_axis)
    endif

    if (diag_file%is_time_to_write(time_step)) then
      call diag_file%increase_unlimited_dimension()
      call diag_file%write_time_data()
      !TODO call diag_file%add_variable_data()
      call diag_file%update_next_write(time_step)
      call diag_file%update_current_new_file_freq_index(time_step)
      if (diag_file%is_time_to_close_file(time_step)) call diag_file%close_diag_file()
    else if (force_write .and. .not. diag_file%is_file_static()) then
      call diag_file%increase_unlimited_dimension()
      call diag_file%write_time_data()
      call diag_file%close_diag_file()
    endif
  enddo
#endif
end subroutine fms_diag_do_io

!> @brief Add a attribute to the diag_obj using the diag_field_id
subroutine fms_diag_field_add_attribute(this, diag_field_id, att_name, att_value)
  class(fmsDiagObject_type), intent (inout) :: this !< The diag object
  integer,          intent(in) :: diag_field_id      !< Id of the axis to add the attribute to
  character(len=*), intent(in) :: att_name     !< Name of the attribute
  class(*),         intent(in) :: att_value(:) !< The attribute value to add
#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
!TODO: Value for diag not found
  if ( diag_field_id .LE. 0 ) THEN
    RETURN
  else
    if (this%FMS_diag_fields(diag_field_id)%is_registered() ) &
      call this%FMS_diag_fields(diag_field_id)%add_attribute(att_name, att_value)
  endif
#endif
end subroutine fms_diag_field_add_attribute

!> @brief Add an attribute to an axis
subroutine fms_diag_axis_add_attribute(this, axis_id, att_name, att_value)
  class(fmsDiagObject_type), intent (inout) :: this !< The diag object
  integer,          intent(in) :: axis_id      !< Id of the axis to add the attribute to
  character(len=*), intent(in) :: att_name     !< Name of the attribute
  class(*),         intent(in) :: att_value(:) !< The attribute value to add

#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
  if (axis_id < 0 .and. axis_id > this%registered_axis) &
    call mpp_error(FATAL, "diag_axis_add_attribute: The axis_id is not valid")

  select type (axis => this%diag_axis(axis_id)%axis)
  type is (fmsDiagFullAxis_type)
    call axis%add_axis_attribute(att_name, att_value)
  end select
#endif
end subroutine fms_diag_axis_add_attribute

#ifdef use_yaml
!> \brief Gets the diag field ID from the module name and field name.
!> \returns a copy of the ID of the diag field or DIAG_FIELD_NOT_FOUND if the field is not registered
PURE FUNCTION fms_get_diag_field_id_from_name(fms_diag_object, module_name, field_name) &
  result(diag_field_id)
  class(fmsDiagObject_type), intent (in) :: fms_diag_object !< The diag object
  CHARACTER(len=*), INTENT(in) :: module_name !< Module name that registered the variable
  CHARACTER(len=*), INTENT(in) :: field_name !< Variable name
  integer :: diag_field_id
  integer :: i !< For looping
!> Initialize to not found
  diag_field_id = DIAG_FIELD_NOT_FOUND
!> Loop through fields to find it.
  if (fms_diag_object%registered_variables < 1) return
  do i=1,fms_diag_object%registered_variables
   diag_field_id = fms_diag_object%FMS_diag_fields(i)%id_from_name(module_name, field_name)
   if(diag_field_id .ne. DIAG_FIELD_NOT_FOUND) return
  enddo
END FUNCTION fms_get_diag_field_id_from_name
#else
!> \brief This replaces the pure function when not compiled with yaml so that an error can be called
!> \returns Error
FUNCTION fms_get_diag_field_id_from_name(fms_diag_object, module_name, field_name) &
  result(diag_field_id)
  class(fmsDiagObject_type), intent (in) :: fms_diag_object !< The diag object
  CHARACTER(len=*), INTENT(in) :: module_name !< Module name that registered the variable
  CHARACTER(len=*), INTENT(in) :: field_name !< Variable name
  integer :: diag_field_id
  integer :: i !< For looping
!> Initialize to not found
  diag_field_id = DIAG_FIELD_NOT_FOUND
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
END FUNCTION fms_get_diag_field_id_from_name
#endif


#ifdef use_yaml
!> returns the buffer object for the given id
!! actual data comes from %get_buffer_data() on the returned object
function get_diag_buffer(this, bufferid) &
result(rslt)
  class(fmsDiagObject_type), intent(in) :: this
  integer, intent(in)                   :: bufferid
  class(fmsDiagBuffer_class),allocatable:: rslt
  if( (bufferid .gt. UBOUND(this%FMS_diag_buffers, 1)) .or. (bufferid .lt. UBOUND(this%FMS_diag_buffers, 1))) &
    call mpp_error(FATAL, 'get_diag_bufer: invalid bufferid given')
  rslt = fms_diag_object%FMS_diag_buffers(bufferid)%diag_buffer_obj
end function
#endif

!> @brief Return the 2D domain for the axis IDs given.
!! @return 2D domain for the axis IDs given
type(domain2d) FUNCTION fms_get_domain2d(this, ids)
  class(fmsDiagObject_type), intent (in) :: this !< The diag object
  INTEGER, DIMENSION(:), INTENT(in) :: ids !< Axis IDs.

#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
fms_get_domain2d = null_domain2d
#else
  INTEGER                      :: type_of_domain !< The type of domain
  CLASS(diagDomain_t), POINTER :: domain         !< Diag Domain pointer

  call get_domain_and_domain_type(fms_diag_object%diag_axis, ids, type_of_domain, domain, "get_domain2d")
  if (type_of_domain .ne. TWO_D_DOMAIN) &
    call mpp_error(FATAL, 'diag_axis_mod::get_domain2d- The axis do not correspond to a 2d Domain')
  select type(domain)
  type is (diagDomain2d_t)
    fms_get_domain2d = domain%domain2
  end select
#endif
END FUNCTION fms_get_domain2d

 !> @brief Gets the length of the axis based on the axis_id
 !> @return Axis_length
 integer function fms_get_axis_length(this, axis_id)
  class(fmsDiagObject_type), intent (in) :: this !< The diag object
  INTEGER, INTENT(in) :: axis_id !< Axis ID of the axis to the length of

#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
fms_get_axis_length = 0
#else
fms_get_axis_length = 0

  if (axis_id < 0 .and. axis_id > this%registered_axis) &
    call mpp_error(FATAL, "fms_get_axis_length: The axis_id is not valid")

  select type (axis => this%diag_axis(axis_id)%axis)
  type is (fmsDiagFullAxis_type)
    fms_get_axis_length = axis%axis_length()
  end select
#endif
end function fms_get_axis_length

!> @brief Gets the name of the axis based on the axis_id
 !> @return The axis_name
function fms_get_axis_name_from_id (this, axis_id) &
result(axis_name)
  class(fmsDiagObject_type), intent (in) :: this !< The diag object
  INTEGER, INTENT(in) :: axis_id !< Axis ID of the axis to the length of

  character (len=:), allocatable :: axis_name

#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
axis_name=" "
#else
    if (axis_id < 0 .and. axis_id > this%registered_axis) &
    call mpp_error(FATAL, "fms_get_axis_length: The axis_id is not valid")

    select type (axis => this%diag_axis(axis_id)%axis)
    type is (fmsDiagFullAxis_type)
      axis_name = axis%get_axis_name()
    end select
#endif
end function fms_get_axis_name_from_id

!> Dumps as much data as it can from the fmsDiagObject_type.
!! Will dump any fields and files as well (see d)
subroutine dump_diag_obj( filename )
  character(len=*), intent(in), optional :: filename !< optional filename to print to,
                                            !! otherwise prints to stdout
#ifdef use_yaml
  !type(fmsDiagObject_type) :: diag_obj
  type(fmsDiagFile_type), pointer :: fileptr !<  pointer for traversing file list
  type(fmsDiagField_type), pointer :: fieldptr !<  pointer for traversing field list
  integer :: i !< do loops
  integer :: unit_num !< unit num of opened log file or stdout

  if( present(filename) ) then
    open(newunit=unit_num, file=trim(filename), action='WRITE')
  else
    unit_num = stdout()
  endif
  if( mpp_pe() .eq. mpp_root_pe()) then
    write(unit_num, *) '********** dumping diag object ***********'
    write(unit_num, *) 'registered_variables:', fms_diag_object%registered_variables
    write(unit_num, *) 'registered_axis:', fms_diag_object%registered_axis
    write(unit_num, *) 'initialized:', fms_diag_object%initialized
    write(unit_num, *) 'files_initialized:', fms_diag_object%files_initialized
    write(unit_num, *) 'fields_initialized:', fms_diag_object%fields_initialized
    write(unit_num, *) 'buffers_initialized:', fms_diag_object%buffers_initialized
    write(unit_num, *) 'axes_initialized:', fms_diag_object%axes_initialized
    write(unit_num, *) 'Files:'
    if( fms_diag_object%files_initialized ) then
      do i=1, SIZE(fms_diag_object%FMS_diag_files)
        write(unit_num, *) 'File num:', i
        fileptr => fms_diag_object%FMS_diag_files(i)%FMS_diag_file
        call fileptr%dump_file_obj(unit_num)
      enddo
    else
      write(unit_num, *) 'files not initialized'
    endif
    if( fms_diag_object%fields_initialized) then
      do i=1, SIZE(fms_diag_object%FMS_diag_fields)
        write(unit_num, *) 'Field num:', i
        fieldptr => fms_diag_object%FMS_diag_fields(i)
        call fieldptr%dump_field_obj(unit_num)
      enddo
    else
      write(unit_num, *) 'fields not initialized'
    endif
    if( present(filename) ) close(unit_num)
  endif
#else
  call mpp_error( FATAL, "You can not use the modern diag manager without compiling with -Duse_yaml")
#endif
end subroutine
end module fms_diag_object_mod
