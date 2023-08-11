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
                         &DIAG_FIELD_NOT_FOUND, diag_not_registered, max_axes, TWO_D_DOMAIN, &
                         &get_base_time, NULL_AXIS_ID, get_var_type, diag_not_registered

  USE time_manager_mod, ONLY: set_time, set_date, get_time, time_type, OPERATOR(>=), OPERATOR(>),&
       & OPERATOR(<), OPERATOR(==), OPERATOR(/=), OPERATOR(/), OPERATOR(+), ASSIGNMENT(=), get_date, &
       & get_ticks_per_second
#ifdef use_yaml
use fms_diag_file_object_mod, only: fmsDiagFileContainer_type, fmsDiagFile_type, fms_diag_files_object_init
use fms_diag_field_object_mod, only: fmsDiagField_type, fms_diag_fields_object_init, get_default_missing_value
use fms_diag_yaml_mod, only: diag_yaml_object_init, diag_yaml_object_end, find_diag_field, &
                           & get_diag_files_id, diag_yaml, get_diag_field_ids, DiagYamlFilesVar_type
use fms_diag_axis_object_mod, only: fms_diag_axis_object_init, fmsDiagAxis_type, fmsDiagSubAxis_type, &
                                   &diagDomain_t, get_domain_and_domain_type, diagDomain2d_t, &
                                   &fmsDiagAxisContainer_type, fms_diag_axis_object_end, fmsDiagFullAxis_type, &
                                   &parse_compress_att, get_axis_id_from_name
use fms_diag_output_buffer_mod
use fms_mod, only: fms_error_handler
use constants_mod, only: SECONDS_PER_DAY
#endif
#if defined(_OPENMP)
use omp_lib
#endif
use mpp_domains_mod, only: domain1d, domain2d, domainUG, null_domain2d
use platform_mod
implicit none
private

type fmsDiagObject_type
!TODO add container arrays
#ifdef use_yaml
private
!TODO: Remove FMS prefix from variables in this type
  class(fmsDiagFileContainer_type), allocatable :: FMS_diag_files (:) !< array of diag files
  class(fmsDiagField_type), allocatable :: FMS_diag_fields(:) !< Array of diag fields
  type(fmsDiagOutputBufferContainer_type), allocatable :: FMS_diag_output_buffers(:) !< array of output buffer objects
                                                                       !! one for each variable in the diag_table.yaml
  integer, private :: registered_buffers = 0 !< number of registered buffers, per dimension
  class(fmsDiagAxisContainer_type), allocatable :: diag_axis(:) !< Array of diag_axis
  type(time_type)  :: current_model_time !< The current model time
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
    procedure :: fms_diag_accept_data
    procedure :: fms_diag_send_complete
    procedure :: fms_diag_do_io
    procedure :: fms_diag_field_add_cell_measures
    procedure :: allocate_diag_field_output_buffers
    procedure :: fms_diag_compare_window
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
  this%buffers_initialized =fms_diag_output_buffer_init(this%FMS_diag_output_buffers,SIZE(diag_yaml%get_diag_fields()))
  this%registered_variables = 0
  this%registered_axis = 0
  this%current_model_time = get_base_time()
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

  call this%fms_diag_do_io(is_end_of_run=.true.)
  !TODO: Deallocate diag object arrays and clean up all memory
  do i=1, size(this%FMS_diag_output_buffers)
    if(allocated(this%FMS_diag_output_buffers(i)%diag_buffer_obj)) then
      call this%FMS_diag_output_buffers(i)%diag_buffer_obj%flush_buffer()
    endif
  enddo
  deallocate(this%FMS_diag_output_buffers)
  this%axes_initialized = fms_diag_axis_object_end(this%diag_axis)
  this%initialized = .false.
  call diag_yaml_object_end
#else
  call mpp_error(FATAL, "You can not call fms_diag_object%end without yaml")
#endif
end subroutine fms_diag_object_end

!> @brief Registers a field.
!! @description This to avoid having duplicate code in each of the _scalar, _array and _static register calls
!! @return field index to be used in subsequent calls to send_data or DIAG_FIELD_NOT_FOUND if the field is not
!! in the diag_table.yaml
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
    !< The field was not found in the table, so return DIAG_FIELD_NOT_FOUND
    fms_register_diag_field_obj = DIAG_FIELD_NOT_FOUND
    deallocate(diag_field_indices)
    return
  endif

  this%registered_variables = this%registered_variables + 1
  fms_register_diag_field_obj = this%registered_variables

  call this%FMS_diag_fields(this%registered_variables)%&
    &setID(this%registered_variables)

!> Use pointers for convenience
  fieldptr => this%FMS_diag_fields(this%registered_variables)

!> Initialize buffer_ids of this field with the diag_field_indices(diag_field_indices)
!! of the sorted variable list
  fieldptr%buffer_ids = get_diag_field_ids(diag_field_indices)
  do i = 1, size(fieldptr%buffer_ids)
    call this%FMS_diag_output_buffers(fieldptr%buffer_ids(i))%set_field_id(this%registered_variables)
    call this%FMS_diag_output_buffers(fieldptr%buffer_ids(i))%set_yaml_id(fieldptr%buffer_ids(i))
  enddo

!> Allocate and initialize member buffer_allocated of this field
  fieldptr%buffer_allocated = .false.

!> Register the data for the field
  call fieldptr%register(modname, varname, diag_field_indices, this%diag_axis, &
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
     call fileptr%add_buffer_id(fieldptr%buffer_ids(i))
     call fileptr%set_file_domain(fieldptr%get_domain(), fieldptr%get_type_of_domain())
     call fileptr%init_diurnal_axis(this%diag_axis, this%registered_axis, diag_field_indices(i))
     call fileptr%add_axes(axes, this%diag_axis, this%registered_axis, diag_field_indices(i), &
       fieldptr%buffer_ids(i), this%FMS_diag_output_buffers)
     call fileptr%add_start_time(init_time, this%current_model_time)
     call fileptr%set_file_time_ops (fieldptr%diag_field(i), fieldptr%is_static())
    enddo
  elseif (present(axes)) then !only axes present
    do i = 1, size(file_ids)
     fileptr => this%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_and_yaml_id(fieldptr%get_id(), diag_field_indices(i))
     call fileptr%add_buffer_id(fieldptr%buffer_ids(i))
     call fileptr%init_diurnal_axis(this%diag_axis, this%registered_axis, diag_field_indices(i))
     call fileptr%set_file_domain(fieldptr%get_domain(), fieldptr%get_type_of_domain())
     call fileptr%add_axes(axes, this%diag_axis, this%registered_axis, diag_field_indices(i), &
       fieldptr%buffer_ids(i), this%FMS_diag_output_buffers)
     call fileptr%set_file_time_ops (fieldptr%diag_field(i), fieldptr%is_static())
    enddo
  elseif (present(init_time)) then !only inti time present
    do i = 1, size(file_ids)
     fileptr => this%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_and_yaml_id(fieldptr%get_id(), diag_field_indices(i))
     call fileptr%add_buffer_id(fieldptr%buffer_ids(i))
     call fileptr%add_start_time(init_time, this%current_model_time)
     call fileptr%set_file_time_ops (fieldptr%diag_field(i), fieldptr%is_static())
    enddo
  else !no axis or init time present
    do i = 1, size(file_ids)
     fileptr => this%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_and_yaml_id(fieldptr%get_id(), diag_field_indices(i))
     call fileptr%add_buffer_id(fieldptr%buffer_ids(i))
     call fileptr%set_file_time_ops (fieldptr%diag_field(i), fieldptr%is_static())
    enddo
  endif
  nullify (fileptr)
  nullify (fieldptr)
  deallocate(diag_field_indices)
#endif
end function fms_register_diag_field_obj

!> @brief Registers a scalar field
!! @return field index to be used in subsequent calls to send_data or DIAG_FIELD_NOT_FOUND if the field is not
!! in the diag_table.yaml
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
fms_register_diag_field_scalar=DIAG_FIELD_NOT_FOUND
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
!! @return field index to be used in subsequent calls to send_data or DIAG_FIELD_NOT_FOUND if the field is not
!! in the diag_table.yaml
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
fms_register_diag_field_array=DIAG_FIELD_NOT_FOUND
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
!! @return field index to be used in subsequent calls to send_data or DIAG_FIELD_NOT_FOUND if the field is not
!! in the diag_table.yaml
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
fms_register_static_field=DIAG_FIELD_NOT_FOUND
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
  !TODO The register_static_field interface does not have the capabiliy to register a variable as a "scalar"
  !     since the axes argument is required, this forced model code to pass in a null_axis_id as an argument
  if (size(axes) .eq. 1 .and. axes(1) .eq. null_axis_id) then
    ! If they are passing in the null_axis_ids, ignore the `axes` argument
    fms_register_static_field = this%register( &
      & module_name, field_name, &
      & longname=long_name, units=units, missing_value=missing_value, varrange=range, &
      & mask_variant=mask_variant, do_not_log=do_not_log, interp_method=interp_method, tile_count=tile_count, &
      & standname=standard_name, area=area, volume=volume, realm=realm, &
      & static=.true.)
  else
    fms_register_static_field = this%register( &
      & module_name, field_name, axes=axes, &
      & longname=long_name, units=units, missing_value=missing_value, varrange=range, &
      & mask_variant=mask_variant, do_not_log=do_not_log, interp_method=interp_method, tile_count=tile_count, &
      & standname=standard_name, area=area, volume=volume, realm=realm, &
      & static=.true.)
  endif
#endif
end function fms_register_static_field

!> @brief Wrapper for the register_diag_axis subroutine. This is needed to keep the diag_axis_init
!! interface the same
!> @return Axis id
FUNCTION fms_diag_axis_init(this, axis_name, axis_data, units, cart_name, axis_length, long_name, direction,&
  & set_name, edges, Domain, Domain2, DomainU, aux, req, tile_count, domain_position ) &
  & result(id)

  class(fmsDiagObject_type),TARGET,INTENT(inout):: this       !< Diaj_obj to fill
  CHARACTER(len=*),   INTENT(in)           :: axis_name       !< Name of the axis
  CLASS(*),           INTENT(in)           :: axis_data(:)    !< Array of coordinate values
  CHARACTER(len=*),   INTENT(in)           :: units           !< Units for the axis
  CHARACTER(len=1),   INTENT(in)           :: cart_name       !< Cartesian axis ("X", "Y", "Z", "T", "U", "N")
  integer,            intent(in)           :: axis_length     !< The length of the axis size(axis_data(:))
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
        call axis%set_edges(edges_name, edges)
      end select
    endif
    call axis%register(axis_name, axis_data, units, cart_name, long_name=long_name, &
      & direction=direction, set_name=set_name, Domain=Domain, Domain2=Domain2, DomainU=DomainU, aux=aux, &
      & req=req, tile_count=tile_count, domain_position=domain_position, axis_length=axis_length)

    id = this%registered_axis
    call axis%set_axis_id(id)
  end select
#endif
end function fms_diag_axis_init

!> Accepts data from the send_data functions.  If this is in an openmp region with more than
!! one thread, the data is buffered in the field object and processed later.  If only a single thread
!! is being used, then the processing can be done and stored in the buffer object.  The hope is that
!! the increase in memory footprint related to buffering can be handled by the shared memory of the
!! multithreaded case.
!! \note If some of the diag manager is offloaded in the future, then it should be treated similarly
!! to the multi-threaded option for processing later
logical function fms_diag_accept_data (this, diag_field_id, field_data, time, is_in, js_in, ks_in, &
                  mask, rmask, ie_in, je_in, ke_in, weight, err_msg)
  class(fmsDiagObject_type),TARGET,INTENT(inout):: this !< Diaj_obj to fill
  INTEGER, INTENT(in) :: diag_field_id !< The ID of the input diagnostic field
  CLASS(*), DIMENSION(:,:,:,:), INTENT(in) :: field_data !< The data for the input diagnostic
  CLASS(*), INTENT(in), OPTIONAL :: weight !< The weight used for averaging
  TYPE (time_type), INTENT(in), OPTIONAL :: time !< The current time
  INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ks_in,ie_in,je_in, ke_in !< Indicies for the variable
  LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask !< The location of the mask
  CLASS(*), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: rmask !< The masking values
  CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg !< An error message returned
  integer :: is, js, ks !< Starting indicies of the field_data
  integer :: ie, je, ke !< Ending indicied of the field_data
  integer :: n1, n2, n3 !< Size of the 3 indicies of the field data
  integer :: omp_num_threads !< Number of openmp threads
  integer :: omp_level !< The openmp active level
  logical :: buffer_the_data !< True if the user selects to buffer the data and run the calculations
                             !! later.  \note This is experimental
  !TODO logical, allocatable, dimension(:,:,:) :: oor_mask !< Out of range mask
  integer :: sample !< Index along the diurnal time axis
  integer :: day    !< Number of days
  integer :: second !< Number of seconds
  integer :: tick   !< Number of ticks representing fractional second
  integer :: buffer_id !< Index of a buffer
  !TODO: logical :: phys_window
  character(len=128) :: error_string !< Store error text
  integer :: i !< For looping
  logical :: data_buffer_is_allocated !< .true. if the data buffer is allocated

#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
  class(diagYamlFilesVar_type), pointer :: ptr_diag_field_yaml !< Pointer to a field from yaml fields

  !TODO: weight is for time averaging where each time level may have a different weight
  ! call real_copy_set()

  !TODO: oor_mask is only used for checking out of range values.
  ! call init_mask_3d()

  !TODO: Check improper combinations of is, ie, js, and je.
  ! if (check_indices_order()) deallocate(oor_mask)

!> Does the user want to push off calculations until send_diag_complete?
  buffer_the_data = .false.
!> initialize the number of threads and level to be 0
  omp_num_threads = 0
  omp_level = 0
#if defined(_OPENMP)
  omp_num_threads = omp_get_num_threads()
  omp_level = omp_get_level()
  buffer_the_data = (omp_num_threads > 1 .AND. omp_level > 0)
#endif
!If this is true, buffer data
  main_if: if (buffer_the_data) then
!> Calculate the i,j,k start and end
    ! If is, js, or ks not present default them to 1
    is = 1
    js = 1
    ks = 1
    IF ( PRESENT(is_in) ) is = is_in
    IF ( PRESENT(js_in) ) js = js_in
    IF ( PRESENT(ks_in) ) ks = ks_in
    n1 = SIZE(field_data, 1)
    n2 = SIZE(field_data, 2)
    n3 = SIZE(field_data, 3)
    ie = is+n1-1
    je = js+n2-1
    ke = ks+n3-1
    IF ( PRESENT(ie_in) ) ie = ie_in
    IF ( PRESENT(je_in) ) je = je_in
    IF ( PRESENT(ke_in) ) ke = ke_in

!> Only 1 thread allocates the output buffer and sets set_math_needs_to_be_done
!$omp critical
    if (.not. this%FMS_diag_fields(diag_field_id)%is_data_buffer_allocated()) then
      data_buffer_is_allocated = &
        this%FMS_diag_fields(diag_field_id)%allocate_data_buffer(field_data, this%diag_axis)
    endif
    call this%FMS_diag_fields(diag_field_id)%set_data_buffer_is_allocated(.TRUE.)
    call this%FMS_diag_fields(diag_field_id)%set_math_needs_to_be_done(.TRUE.)
!$omp end critical
    call this%FMS_diag_fields(diag_field_id)%set_data_buffer(field_data,&
                                                             is, js, ks, ie, je, ke)
    fms_diag_accept_data = .TRUE.
    return
  else
!!TODO: Loop through fields and do averages/math functions

    call this%allocate_diag_field_output_buffers(field_data, diag_field_id)
    do i = 1, size(this%FMS_diag_fields(diag_field_id)%buffer_ids)
      buffer_id = this%FMS_diag_fields(diag_field_id)%buffer_ids(i)

      !!TODO: Check if the field is a physics window
      !! phys_window = fms_diag_compare_window()

      !!TODO: Get local start and end indices on 3 axes for regional output

      !> Compute the diurnal index
      sample = 1
      if (present(time)) then
        call get_time(time, second, day, tick) !< Current time in days and seconds
        ptr_diag_field_yaml => diag_yaml%get_diag_field_from_id(buffer_id)
        sample = floor((second + real(tick) / get_ticks_per_second()) &
          & * ptr_diag_field_yaml%get_n_diurnal() / SECONDS_PER_DAY) + 1
      end if

      !!TODO: Get the vertical layer start and end indices

      !!TODO: Initialize output time for fields output every time step

      !< Check if time should be present for this field
      if (.not.this%FMS_diag_fields(diag_field_id)%is_static() .and. .not.present(time)) then
        write(error_string, '(a,"/",a)') trim(this%FMS_diag_fields(diag_field_id)%get_modname()),&
          & trim(this%FMS_diag_fields(diag_field_id)%diag_field(i)%get_var_outname())
        if (fms_error_handler('fms_diag_object_mod::fms_diag_accept_data', 'module/output_name: '&
          &//trim(error_string)//', time must be present for nonstatic field', err_msg)) then
            !!TODO: deallocate local pointers/allocatables if needed
          return
        end if
      end if

      !!TODO: Is it time to output for this field? CAREFUL ABOUT > vs >= HERE
      !--- The fields send out within openmp parallel region will be written out in
      !--- diag_send_complete.

      !!TODO: Is check to bounds of current field necessary?

      !!TODO: Take care of submitted field data

    enddo
    call this%FMS_diag_fields(diag_field_id)%set_math_needs_to_be_done(.FALSE.)
    fms_diag_accept_data = .TRUE.
    return
  end if main_if
!> Return false if nothing is done
  fms_diag_accept_data = .FALSE.
  return
#endif
end function fms_diag_accept_data
!! TODO: This entire routine
!> @brief Loops through all the files, open the file, writes out axis and
!! variable metadata and data when necessary.
subroutine fms_diag_send_complete(this, time_step)
  class(fmsDiagObject_type), target, intent (inout) :: this      !< The diag object
  TYPE (time_type),                  INTENT(in)     :: time_step !< The time_step

  integer :: i !< For do loops

  integer :: ifile !< For file loops
  integer :: ifield !< For field loops
#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else

  class(fmsDiagFileContainer_type), pointer :: diag_file !< Pointer to this%FMS_diag_files(i) (for convenience
  class(fmsDiagField_type), pointer :: diag_field !< Pointer to this%FMS_diag_files(i)%diag_field(j)
  logical :: math !< True if the math functions need to be called using the data buffer,
  !! False if the math functions were done in accept_data
  integer, dimension(:), allocatable :: file_field_ids !< Array of field IDs for a file

  !< Update the current model time by adding the time_step
  this%current_model_time = this%current_model_time + time_step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! In the future, this may be parallelized for offloading
  file_loop: do ifile = 1, size(this%FMS_diag_files)
    diag_file => this%FMS_diag_files(ifile)
    field_outer_if: if (size(diag_file%FMS_diag_file%get_field_ids()) .ge. 1) then
      allocate (file_field_ids(size(diag_file%FMS_diag_file%get_field_ids() )))
      file_field_ids = diag_file%FMS_diag_file%get_field_ids()
      field_loop: do ifield = 1, size(file_field_ids)
        ! If the field is not registered go away
        if (.not. diag_file%FMS_diag_file%is_field_registered(ifield)) cycle

        diag_field => this%FMS_diag_fields(file_field_ids(ifield))
        !> Check if math needs to be done
        math = diag_field%get_math_needs_to_be_done()
        calling_math: if (math) then
          call this%allocate_diag_field_output_buffers(diag_field%get_data_buffer(), file_field_ids(ifield))
          !!TODO: call math functions !!
        endif calling_math
        !> Clean up, clean up, everybody everywhere
        if (associated(diag_field)) nullify(diag_field)
      enddo field_loop
      !> Clean up, clean up, everybody do your share
      if (allocated(file_field_ids)) deallocate(file_field_ids)
    endif field_outer_if
  enddo file_loop

  call this%fms_diag_do_io()
#endif

end subroutine fms_diag_send_complete

!> @brief Loops through all the files, open the file, writes out axis and
!! variable metadata and data when necessary.
subroutine fms_diag_do_io(this, is_end_of_run)
  class(fmsDiagObject_type), target, intent(inout)  :: this          !< The diag object
  logical,                 optional, intent(in)     :: is_end_of_run !< If .true. this is the end of the run,
                                                                     !! so force write
#ifdef use_yaml
  integer :: i !< For do loops
  class(fmsDiagFileContainer_type), pointer :: diag_file !< Pointer to this%FMS_diag_files(i) (for convenience)
  TYPE (time_type),                 pointer :: model_time!< The current model time

  logical :: file_is_opened_this_time_step !< True if the file was opened in this time_step
                                           !! If true the metadata will need to be written
  logical :: force_write

  force_write = .false.
  if (present (is_end_of_run)) force_write = .true.

  model_time => this%current_model_time

  do i = 1, size(this%FMS_diag_files)
    diag_file => this%FMS_diag_files(i)

    !< Go away if the file is a subregional file and the current PE does not have any data for it
    if (.not. diag_file%writing_on_this_pe()) cycle

    call diag_file%open_diag_file(model_time, file_is_opened_this_time_step)
    if (file_is_opened_this_time_step) then
      call diag_file%write_global_metadata()
      call diag_file%write_axis_metadata(this%diag_axis)
      call diag_file%write_time_metadata()
      call diag_file%write_field_metadata(this%FMS_diag_fields, this%diag_axis)
      call diag_file%write_axis_data(this%diag_axis)
    endif

    if (diag_file%is_time_to_write(model_time)) then
      call diag_file%increase_unlim_dimension_level()
      call diag_file%write_time_data()
      call diag_file%write_field_data(this%FMS_diag_fields, this%FMS_diag_output_buffers)
      call diag_file%update_next_write(model_time)
      call diag_file%update_current_new_file_freq_index(model_time)
      if (diag_file%is_time_to_close_file(model_time)) call diag_file%close_diag_file()
    else if (force_write) then
      if (diag_file%get_unlim_dimension_level() .eq. 0) then
        call diag_file%increase_unlim_dimension_level()
        call diag_file%write_time_data()
      endif
      call diag_file%close_diag_file()
    endif
  enddo
#endif
end subroutine fms_diag_do_io

!> @brief Adds the diag ids of the Area and or Volume of the diag_field_object
subroutine fms_diag_field_add_cell_measures(this, diag_field_id, area, volume)
  class(fmsDiagObject_type), intent (inout) :: this          !< The diag object
  integer,                   intent(in)     :: diag_field_id !< diag_field to add the are and volume to
  INTEGER, optional,         INTENT(in)     :: area          !< diag ids of area
  INTEGER, optional,         INTENT(in)     :: volume        !< diag ids of volume

#ifndef use_yaml
  CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
  call this%FMS_diag_fields(diag_field_id)%add_area_volume(area, volume)
#endif
end subroutine fms_diag_field_add_cell_measures

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

  character(len=20) :: axis_names(2) !< Names of the uncompress axis
  integer           :: uncmx_ids(2)  !< Ids of the uncompress axis
  integer           :: j             !< For do loops
#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
  if (axis_id < 0 .and. axis_id > this%registered_axis) &
    call mpp_error(FATAL, "diag_axis_add_attribute: The axis_id is not valid")

  select type (axis => this%diag_axis(axis_id)%axis)
  type is (fmsDiagFullAxis_type)
    call axis%add_axis_attribute(att_name, att_value)

    !! Axis that are in the "unstructured" domain require a "compress" attribute for the
    !! combiner and PP. This attribute is passed in via a diag_axis_add_attribute call in the model code
    !! The compress attribute indicates the names of the axis that were compressed
    !! For example grid_index:compress = "grid_yt grid_xt"
    !! The metadata and the data for these axis also needs to be written to the file
    if (trim(att_name) .eq. "compress") then
      !< If the attribute is the "compress" attribute, get the axis names,
      !! and the ids of the axis and add it to the axis object so it can be written to netcdf files
      !! that use this axis
      axis_names = parse_compress_att(att_value)
      do j = 1, size(axis_names)
        uncmx_ids(j) = get_axis_id_from_name(axis_names(j), this%diag_axis, this%registered_axis)
        if (uncmx_ids(j) .eq. diag_null) call mpp_error(FATAL, &
          &"Error parsing the compress attribute for axis: "//trim(axis%get_axis_name())//&
          &". Be sure that the axes in the compress attribute are registered")
      enddo
      call axis%add_structured_axis_ids(uncmx_ids)
    endif
  end select
#endif
end subroutine fms_diag_axis_add_attribute

!> \brief Gets the diag field ID from the module name and field name.
!> \returns a copy of the ID of the diag field or DIAG_FIELD_NOT_FOUND if the field is not registered
FUNCTION fms_get_diag_field_id_from_name(this, module_name, field_name) &
  result(diag_field_id)
  class(fmsDiagObject_type), intent (in) :: this !< The diag object, the caller
  CHARACTER(len=*), INTENT(in) :: module_name !< Module name that registered the variable
  CHARACTER(len=*), INTENT(in) :: field_name !< Variable name
  integer :: diag_field_id

#ifdef use_yaml
  integer              :: i                     !< For looping
  integer, allocatable :: diag_field_indices(:) !< indices where the field was found in the yaml

  diag_field_id = DIAG_FIELD_NOT_FOUND

  !> Loop through fields to find it.
  do i=1, this%registered_variables
    !< Check if the field was registered, if it was return the diag_field_id
    diag_field_id = this%FMS_diag_fields(i)%id_from_name(module_name, field_name)
    if(diag_field_id .ne. DIAG_FIELD_NOT_FOUND) return
  enddo

  !< Check if the field is in the diag_table.yaml. If it is, return DIAG_FIELD_NOT_REGISTERED
  !! Otherwsie it will return DIAG_FIELD_NOT_FOUND
  diag_field_indices = find_diag_field(field_name, module_name)
  if (diag_field_indices(1) .ne. diag_null) then
    diag_field_id = DIAG_NOT_REGISTERED
  endif
  deallocate(diag_field_indices)
#else
  diag_field_id = DIAG_FIELD_NOT_FOUND
  CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#endif
END FUNCTION fms_get_diag_field_id_from_name

#ifdef use_yaml
!> returns the buffer object for the given id
!! actual data comes from %get_buffer_data() on the returned object
function get_diag_buffer(this, bufferid) &
result(rslt)
  class(fmsDiagObject_type), intent(in) :: this
  integer, intent(in)                   :: bufferid
  class(fmsDiagOutputBuffer_class),allocatable:: rslt
  if( (bufferid .gt. UBOUND(this%FMS_diag_output_buffers, 1)) .or. &
      (bufferid .lt. LBOUND(this%FMS_diag_output_buffers, 1))) &
    call mpp_error(FATAL, 'get_diag_bufer: invalid bufferid given')
  rslt = this%FMS_diag_output_buffers(bufferid)%diag_buffer_obj
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
  type is (fmsDiagSubAxis_type)
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

    !! if its a scalar (null axis id) just returns n/a since no axis is defined
    if (axis_id .eq. NULL_AXIS_ID) then
      allocate(character(len=3) :: axis_name)
      axis_name = "n/a"
      return
    endif


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

!> @brief Allocates the output buffers of the fields corresponding to the registered variable
!! Input arguments are the field and its ID passed to routine fms_diag_accept_data()
subroutine allocate_diag_field_output_buffers(this, field_data, field_id)
  class(fmsDiagObject_type), target, intent(inout) :: this !< diag object
  class(*), dimension(:,:,:,:), intent(in) :: field_data !< field data
  integer, intent(in) :: field_id !< Id of the field data
#ifdef use_yaml
  integer :: ndims !< Number of dimensions in the input field data
  integer :: buffer_id !< Buffer index of FMS_diag_buffers
  integer :: num_diurnal_samples !< Number of diurnal samples from diag_yaml
  integer, allocatable :: axes_length(:) !< Length of each axis
  integer :: i, j !< For looping
  class(fmsDiagOutputBuffer_class), pointer :: ptr_diag_buffer_obj !< Pointer to the buffer class
  class(DiagYamlFilesVar_type), pointer :: ptr_diag_field_yaml !< Pointer to a field from yaml fields
  integer, allocatable :: axis_ids(:) !< Pointer to indices of axes of the field variable
  integer :: var_type !< Stores type of the field data (r4, r8, i4, i8, and string) represented as an integer.
  class(*), allocatable :: missing_value !< Missing value to initialize the data to
  character(len=128), allocatable :: var_name !< Field name to initialize output buffers
  logical :: is_scalar !< Flag indicating that the variable is a scalar
  integer :: yaml_id

  if (this%FMS_diag_fields(field_id)%buffer_allocated) return

  ! Determine the type of the field data
  var_type = get_var_type(field_data(1, 1, 1, 1))

  ! Get variable/field name
  var_name = this%Fms_diag_fields(field_id)%get_varname()

  ! Get missing value for the field
  !TODO class (*) is weird missing_value = this%FMS_diag_fields(field_id)%get_missing_value(var_type)
  !!should work ...
  if (this%FMS_diag_fields(field_id)%has_missing_value()) then
    select type (my_type => this%FMS_diag_fields(field_id)%get_missing_value(var_type))
      type is (real(kind=r4_kind))
        missing_value = real(my_type, kind=r4_kind)
      type is (real(kind=r8_kind))
        missing_value = real(my_type, kind=r8_kind)
      class default
        call mpp_error( FATAL, 'fms_diag_object_mod:allocate_diag_field_output_buffers Invalid type')
    end select
  else
    select type (my_type => get_default_missing_value(var_type))
      type is (real(kind=r4_kind))
        missing_value = real(my_type, kind=r4_kind)
      type is (real(kind=r8_kind))
        missing_value = real(my_type, kind=r8_kind)
      class default
        call mpp_error( FATAL, 'fms_diag_object_mod:allocate_diag_field_output_buffers Invalid type')
    end select
  endif

  ! Determine dimensions of the field
  is_scalar = this%FMS_diag_fields(field_id)%is_scalar()

  ! Loop over a number of fields/buffers where this variable occurs
  do i = 1, size(this%FMS_diag_fields(field_id)%buffer_ids)
    buffer_id = this%FMS_diag_fields(field_id)%buffer_ids(i)

    ndims = 0
    if (.not. is_scalar) then
      axis_ids = this%FMS_diag_output_buffers(buffer_id)%get_axis_ids()
      ndims = size(axis_ids)
    endif

    yaml_id = this%FMS_diag_output_buffers(buffer_id)%get_yaml_id()

    ptr_diag_field_yaml => diag_yaml%diag_fields(yaml_id)
    num_diurnal_samples = ptr_diag_field_yaml%get_n_diurnal() !< Get number of diurnal samples

    ! If diurnal axis exists, fill lengths of axes.
    if (num_diurnal_samples .ne. 0) then
      allocate(axes_length(ndims + 1)) !< Include extra length for the diurnal axis
    else
      allocate(axes_length(ndims))
    endif

    do j = 1, ndims
      axes_length(j) = this%fms_get_axis_length(axis_ids(j))
    enddo

    if (num_diurnal_samples .ne. 0) then
      axes_length(ndims + 1) = num_diurnal_samples
      ndims = ndims + 1 !< Add one more dimension for the diurnal axis
    endif

    ! Allocates diag_buffer_obj to the correct outputBuffer type based on the dimension:
    ! outputBuffer0d_type, outputBuffer1d_type, outputBuffer2d_type, outputBuffer3d_type,
    ! outputBuffer4d_type or outputBuffer5d_type.
    if (.not. allocated(this%FMS_diag_output_buffers(buffer_id)%diag_buffer_obj)) then
      call fms_diag_output_buffer_create_container(ndims, this%FMS_diag_output_buffers(buffer_id))
    end if

    ptr_diag_buffer_obj => this%FMS_diag_output_buffers(buffer_id)%diag_buffer_obj

    select type (ptr_diag_buffer_obj)
      type is (outputBuffer0d_type) !< Scalar buffer
        if (allocated(ptr_diag_buffer_obj%buffer)) cycle !< If allocated, loop back
        call ptr_diag_buffer_obj%allocate_buffer(field_data(1, 1, 1, 1), & !< If scalar field variable
          this%FMS_diag_fields(field_id)%get_varname())
        call ptr_diag_buffer_obj%initialize_buffer(missing_value, var_name)
      type is (outputBuffer1d_type) !< 1D buffer
        if (allocated(ptr_diag_buffer_obj%buffer)) cycle !< If allocated, loop back
        call ptr_diag_buffer_obj%allocate_buffer(field_data(1, 1, 1, 1), axes_length(1), &
          this%FMS_diag_fields(field_id)%get_varname(), num_diurnal_samples)
        call ptr_diag_buffer_obj%initialize_buffer(missing_value, var_name)
      type is (outputBuffer2d_type) !< 2D buffer
        if (allocated(ptr_diag_buffer_obj%buffer)) cycle !< If allocated, loop back
        call ptr_diag_buffer_obj%allocate_buffer(field_data(1, 1, 1, 1), axes_length(1:2), &
          this%FMS_diag_fields(field_id)%get_varname(), num_diurnal_samples)
        call ptr_diag_buffer_obj%initialize_buffer(missing_value, var_name)
      type is (outputBuffer3d_type) !< 3D buffer
        if (allocated(ptr_diag_buffer_obj%buffer)) cycle !< If allocated, loop back
        call ptr_diag_buffer_obj%allocate_buffer(field_data(1, 1, 1, 1), axes_length(1:3), &
          this%FMS_diag_fields(field_id)%get_varname(), num_diurnal_samples)
          call ptr_diag_buffer_obj%initialize_buffer(missing_value, var_name)
      type is (outputBuffer4d_type) !< 4D buffer
        if (allocated(ptr_diag_buffer_obj%buffer)) cycle !< If allocated, loop back
        call ptr_diag_buffer_obj%allocate_buffer(field_data(1, 1, 1, 1), axes_length(1:4), &
          this%FMS_diag_fields(field_id)%get_varname(), num_diurnal_samples)
        call ptr_diag_buffer_obj%initialize_buffer(missing_value, var_name)
      type is (outputBuffer5d_type) !< 5D buffer
        if (allocated(ptr_diag_buffer_obj%buffer)) cycle !< If allocated, loop back
        call ptr_diag_buffer_obj%allocate_buffer(field_data(1, 1, 1, 1), axes_length(1:5), &
          this%FMS_diag_fields(field_id)%get_varname(), num_diurnal_samples)
        call ptr_diag_buffer_obj%initialize_buffer(missing_value, var_name)
      class default
        call mpp_error( FATAL, 'allocate_diag_field_output_buffers: invalid buffer type')
    end select

    if (allocated(axis_ids)) deallocate(axis_ids)
    deallocate(axes_length)
  enddo

  this%FMS_diag_fields(field_id)%buffer_allocated = .true.
#else
  call mpp_error( FATAL, "allocate_diag_field_output_buffers: "//&
    "you can not use the modern diag manager without compiling with -Duse_yaml")
#endif
end subroutine allocate_diag_field_output_buffers

!> @brief Determines if the window defined by the input bounds is a physics window.
!> @return TRUE if the window size is less then the actual field size else FALSE.
function fms_diag_compare_window(this, field, field_id, &
  is_in, ie_in, js_in, je_in, ks_in, ke_in) result(is_phys_win)
  class(fmsDiagObject_type), intent(in) :: this !< Diag Object
  class(*), intent(in) :: field(:,:,:,:) !< Field data
  integer, intent(in) :: field_id !< ID of the input field
  integer, intent(in) :: is_in, js_in !< Starting field indices for the first 2 dimensions;
                                      !< pass reconditioned indices fis and fjs
                                      !< which are computed elsewhere.
  integer, intent(in) :: ie_in, je_in !< Ending field indices for the first 2 dimensions;
                                      !< pass reconditioned indices fie and fje
                                      !< which are computed elsewhere.
  integer, intent(in) :: ks_in, ke_in !< Starting and ending indices of the field in 3rd dimension
  logical :: is_phys_win !< Return flag
#ifdef use_yaml
  integer, pointer :: axis_ids(:)
  integer :: total_elements
  integer :: i !< For do loop
  integer :: field_size
  integer, allocatable :: field_shape(:) !< Shape of the field data
  integer :: window_size

  !> Determine shape of the field defined by the input bounds
  field_shape = shape(field(is_in:ie_in, js_in:je_in, ks_in:ke_in, :))

  window_size = field_shape(1) * field_shape(2) * field_shape(3)

  total_elements = 1
  axis_ids => this%FMS_diag_fields(field_id)%get_axis_id()
  do i=1, size(axis_ids)
    total_elements = total_elements * this%fms_get_axis_length(axis_ids(i))
  enddo

  if (total_elements > window_size) then
    is_phys_win = .true.
  else
    is_phys_win = .false.
  end if
#else
  is_phys_win = .false.
  call mpp_error( FATAL, "fms_diag_compare_window: "//&
    "you can not use the modern diag manager without compiling with -Duse_yaml")
#endif
end function fms_diag_compare_window
end module fms_diag_object_mod
