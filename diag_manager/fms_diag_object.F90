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
                         &get_base_time, NULL_AXIS_ID, get_var_type, diag_not_registered, &
                         &time_none, time_max, time_min, time_sum, time_average, time_diurnal, &
                         &time_power, time_rms, r8

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
use fms_diag_reduction_methods_mod, only: check_indices_order, init_mask, set_weight
use constants_mod, only: SECONDS_PER_DAY
#endif
USE fms_diag_bbox_mod, ONLY: fmsDiagIbounds_type, determine_if_block_is_in_region
#if defined(_OPENMP)
use omp_lib
#endif
use mpp_domains_mod, only: domain1d, domain2d, domainUG, null_domain2d
use fms_string_utils_mod, only: string
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
  type(fmsDiagOutputBuffer_type), allocatable :: FMS_diag_output_buffers(:) !< array of output buffer objects
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
    procedure :: fms_get_field_name_from_id
    procedure :: fms_get_axis_name_from_id
    procedure :: fms_diag_accept_data
    procedure :: fms_diag_send_complete
    procedure :: fms_diag_do_io
    procedure :: fms_diag_do_reduction
    procedure :: fms_diag_field_add_cell_measures
    procedure :: allocate_diag_field_output_buffers
    procedure :: fms_diag_compare_window
    procedure :: update_current_model_time
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

  call this%fms_diag_do_io(end_time=time)
  !TODO: Deallocate diag object arrays and clean up all memory
  do i=1, size(this%FMS_diag_output_buffers)
    call this%FMS_diag_output_buffers(i)%flush_buffer()
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
 LOGICAL,          OPTIONAL,     INTENT(in)    :: mask_variant          !< .True. if mask changes over time
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

 class (fmsDiagFile_type), pointer :: fileptr !< Pointer to the diag_file
 class (fmsDiagField_type), pointer :: fieldptr !< Pointer to the diag_field
 class (fmsDiagOutputBuffer_type), pointer :: bufferptr !< Pointer to the output buffer
 class (diagYamlFilesVar_type), pointer :: yamlfptr !< Pointer to yaml object to get the reduction method
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
!> Get the file IDs from the field indicies from the yaml
  file_ids = get_diag_files_id(diag_field_indices)
  call fieldptr%set_file_ids(file_ids)

!> Allocate and initialize member buffer_allocated of this field
  fieldptr%buffer_allocated = .false.
  fieldptr%buffer_ids = get_diag_field_ids(diag_field_indices)

!> Register the data for the field
  call fieldptr%register(modname, varname, diag_field_indices, this%diag_axis, &
       axes=axes, longname=longname, units=units, missing_value=missing_value, varRange= varRange, &
       mask_variant= mask_variant, standname=standname, do_not_log=do_not_log, err_msg=err_msg, &
       interp_method=interp_method, tile_count=tile_count, area=area, volume=volume, realm=realm, &
       static=static)

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

  !> Initialize buffer_ids of this field with the diag_field_indices(diag_field_indices)
!! of the sorted variable list
  do i = 1, size(fieldptr%buffer_ids)
    bufferptr => this%FMS_diag_output_buffers(fieldptr%buffer_ids(i))
    call bufferptr%set_field_id(this%registered_variables)
    call bufferptr%set_yaml_id(fieldptr%buffer_ids(i))
    ! check if diurnal reduction for this buffer and if so set the diurnal sample size
    yamlfptr => diag_yaml%diag_fields(fieldptr%buffer_ids(i))
    if( yamlfptr%get_var_reduction() .eq. time_diurnal) then
      call bufferptr%set_diurnal_sample_size(yamlfptr%get_n_diurnal())
    endif
    call bufferptr%init_buffer_time(init_time)
    call bufferptr%set_next_output(this%FMS_diag_files(file_ids(i))%get_next_output(), fieldptr%is_static())
  enddo

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
    LOGICAL,          OPTIONAL, INTENT(in) :: mask_variant  !< .True. if mask changes over time
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
    LOGICAL,                        OPTIONAL, INTENT(in) :: mask_variant  !< .True. if mask changes over time
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
logical function fms_diag_accept_data (this, diag_field_id, field_data, mask, rmask, &
                                       time, is_in, js_in, ks_in, &
                                       ie_in, je_in, ke_in, weight, err_msg)
  class(fmsDiagObject_type),TARGET,      INTENT(inout)          :: this          !< Diaj_obj to fill
  INTEGER,                               INTENT(in)             :: diag_field_id !< The ID of the diag field
  CLASS(*), DIMENSION(:,:,:,:),          INTENT(in)             :: field_data    !< The data for the diag_field
  LOGICAL,  allocatable,                 INTENT(in)             :: mask(:,:,:,:) !< Logical mask indicating the grid
                                                                                 !! points to mask (null if no mask)
  CLASS(*), allocatable,                 INTENT(in)             :: rmask(:,:,:,:)!< real mask indicating the grid
                                                                                 !! points to mask (null if no mask)
  CLASS(*),                              INTENT(in),   OPTIONAL :: weight        !< The weight used for averaging
  TYPE (time_type),                      INTENT(in),   OPTIONAL :: time          !< The current time
  INTEGER,                               INTENT(in),   OPTIONAL :: is_in, js_in, ks_in !< Starting indices
  INTEGER,                               INTENT(in),   OPTIONAL :: ie_in, je_in, ke_in !< Ending indices
  CHARACTER(len=*),                      INTENT(out),  OPTIONAL :: err_msg       !< An error message returned

  integer                                  :: is, js, ks      !< Starting indicies of the field_data
  integer                                  :: ie, je, ke      !< Ending indicies of the field_data
  integer                                  :: omp_num_threads !< Number of openmp threads
  integer                                  :: omp_level       !< The openmp active level
  logical                                  :: buffer_the_data !< True if the user selects to buffer the data and run
                                                              !! the calculationslater.  \note This is experimental
  character(len=128)                       :: error_string    !< Store error text
  logical                                  :: data_buffer_is_allocated !< .true. if the data buffer is allocated
  character(len=256)                       :: field_info      !< String holding info about the field to append to the
                                                              !! error message
  logical, allocatable, dimension(:,:,:,:) :: oor_mask        !< Out of range mask
  real(kind=r8_kind)                       :: field_weight    !< Weight to use when averaging (it will be converted
                                                              !! based on the type of field_data when doing the math)
  type(fmsDiagIbounds_type)                :: bounds          !< Bounds (starting ending indices) for the field
  logical                                  :: has_halos       !< .True. if field_data contains halos
  logical                                  :: using_blocking  !< .True. if field_data is passed in blocks
#ifndef use_yaml
CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
  field_info = " Check send data call for field:"//trim(this%FMS_diag_fields(diag_field_id)%get_varname())//&
    " and module:"//trim(this%FMS_diag_fields(diag_field_id)%get_modname())

  !< Check if time should be present for this field
  if (.not.this%FMS_diag_fields(diag_field_id)%is_static() .and. .not.present(time)) &
    call mpp_error(FATAL, "Time must be present if the field is not static. "//trim(field_info))

  !< Set the field_weight. If "weight" is not present it will be set to 1.0_r8_kind
  field_weight = set_weight(weight)

  !< Check that the indices are present in the correct combination
  error_string = check_indices_order(is_in, ie_in, js_in, je_in)
  if (trim(error_string) .ne. "") call mpp_error(FATAL, trim(error_string)//". "//trim(field_info))

  using_blocking = .false.
  if ((present(is_in) .and. .not. present(ie_in)) .or. (present(js_in) .and. .not. present(je_in))) &
    using_blocking = .true.

  has_halos = .false.
  if ((present(is_in) .and. present(ie_in)) .or. (present(js_in) .and. present(je_in))) &
    has_halos = .true.

  !< If the field has `mask_variant=.true.`, check that mask OR rmask are present
  if (this%FMS_diag_fields(diag_field_id)%is_mask_variant()) then
    if (.not. allocated(mask) .and. .not. allocated(rmask)) call mpp_error(FATAL, &
      "The field was registered with mask_variant, but mask or rmask are not present in the send_data call. "//&
      trim(field_info))
  endif

  !< Check that mask and rmask are not both present
  if (allocated(mask) .and. allocated(rmask)) call mpp_error(FATAL, &
    "mask and rmask are both present in the send_data call. "//&
    trim(field_info))

  !< Create the oor_mask based on the "mask" and "rmask" arguments
  oor_mask = init_mask(rmask, mask, field_data)

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

  !> Calculate the i,j,k start and end
  ! If is, js, or ks not present default them to 1
  is = 1
  js = 1
  ks = 1
  IF ( PRESENT(is_in) ) is = is_in
  IF ( PRESENT(js_in) ) js = js_in
  IF ( PRESENT(ks_in) ) ks = ks_in
  ie = is+SIZE(field_data, 1)-1
  je = js+SIZE(field_data, 2)-1
  ke = ks+SIZE(field_data, 3)-1
  IF ( PRESENT(ie_in) ) ie = ie_in
  IF ( PRESENT(je_in) ) je = je_in
  IF ( PRESENT(ke_in) ) ke = ke_in

  !If this is true, buffer data
  main_if: if (buffer_the_data) then
!> Only 1 thread allocates the output buffer and sets set_math_needs_to_be_done
!$omp critical

    if (present(time)) call this%update_current_model_time(time)

    !< These set_* calls need to be done inside an omp_critical to avoid any race conditions
    !! and allocation issues
    if(has_halos) call this%FMS_diag_fields(diag_field_id)%set_halo_present()

    !< Set the variable type based off passed in field data
    if(.not. this%FMS_diag_fields(diag_field_id)%has_vartype()) &
      call this%FMS_diag_fields(diag_field_id)%set_type(field_data(1,1,1,1))

    if (allocated(mask) .or. allocated(rmask)) then
      call this%FMS_diag_fields(diag_field_id)%set_var_is_masked(.True.)
    else
      call this%FMS_diag_fields(diag_field_id)%set_var_is_masked(.False.)
    endif

    if (.not. this%FMS_diag_fields(diag_field_id)%is_data_buffer_allocated()) then
      data_buffer_is_allocated = &
        this%FMS_diag_fields(diag_field_id)%allocate_data_buffer(field_data, this%diag_axis)
      if(.not. this%FMS_diag_fields(diag_field_id)%has_mask_allocated()) &
        call this%FMS_diag_fields(diag_field_id)%allocate_mask(oor_mask, this%diag_axis)
    endif
    call this%FMS_diag_fields(diag_field_id)%set_data_buffer_is_allocated(.TRUE.)
    call this%FMS_diag_fields(diag_field_id)%set_math_needs_to_be_done(.TRUE.)
!$omp end critical
    call this%FMS_diag_fields(diag_field_id)%set_data_buffer(field_data, field_weight, &
                                                             is, js, ks, ie, je, ke)
    call this%FMS_diag_fields(diag_field_id)%set_mask(oor_mask, field_info, is, js, ks, ie, je, ke)
    fms_diag_accept_data = .TRUE.
    return
  else
    if (present(time)) call this%update_current_model_time(time)

    !< At this point if we are no longer in an openmp region or running with 1 thread
    !! so it is safe to have these set_* calls
    if(has_halos) call this%FMS_diag_fields(diag_field_id)%set_halo_present()

    !< Set the variable type based off passed in field data
    if(.not. this%FMS_diag_fields(diag_field_id)%has_vartype()) &
      call this%FMS_diag_fields(diag_field_id)%set_type(field_data(1,1,1,1))

    if (allocated(mask) .or. allocated(rmask)) then
      call this%FMS_diag_fields(diag_field_id)%set_var_is_masked(.True.)
    else
      call this%FMS_diag_fields(diag_field_id)%set_var_is_masked(.False.)
    endif

    error_string = bounds%set_bounds(field_data, is, ie, js, je, ks, ke, has_halos)
    if (trim(error_string) .ne. "") call mpp_error(FATAL, trim(error_string)//". "//trim(field_info))

    call this%allocate_diag_field_output_buffers(field_data, diag_field_id)
    error_string = this%fms_diag_do_reduction(field_data, diag_field_id, oor_mask, field_weight, &
      bounds, using_blocking, Time=Time)
    if (trim(error_string) .ne. "") call mpp_error(FATAL, trim(error_string)//". "//trim(field_info))
    call this%FMS_diag_fields(diag_field_id)%set_math_needs_to_be_done(.FALSE.)
    if(.not. this%FMS_diag_fields(diag_field_id)%has_mask_allocated()) &
      call this%FMS_diag_fields(diag_field_id)%allocate_mask(oor_mask)
    call this%FMS_diag_fields(diag_field_id)%set_mask(oor_mask, field_info)
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
  class(*), pointer :: input_data_buffer(:,:,:,:)
  character(len=128) :: error_string
  type(fmsDiagIbounds_type) :: bounds
  integer, dimension(:), allocatable :: file_ids !< Array of file IDs for a field
  logical, parameter :: DEBUG_SC = .false. !< turn on output for debugging

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! In the future, this may be parallelized for offloading
  ! loop through each field
  field_loop: do ifield = 1, size(this%FMS_diag_fields)
    diag_field => this%FMS_diag_fields(ifield)
    if(.not. diag_field%is_registered()) cycle
    if(DEBUG_SC) call mpp_error(NOTE, "fms_diag_send_complete:: var: "//diag_field%get_varname())
    ! get files the field is in
    allocate (file_ids(size(diag_field%get_file_ids() )))
    file_ids = diag_field%get_file_ids()
    math = diag_field%get_math_needs_to_be_done()
    ! if doing math loop through each file for given field
    doing_math: if (size(file_ids) .ge. 1 .and. math) then
      ! Check if buffer alloc'd
      has_input_buff: if (diag_field%has_input_data_buffer()) then
        input_data_buffer => diag_field%get_data_buffer()
        ! reset bounds, allocate output buffer, and update it with reduction
        call bounds%reset_bounds_from_array_4D(input_data_buffer)
        call this%allocate_diag_field_output_buffers(input_data_buffer, ifield)
        error_string = this%fms_diag_do_reduction(input_data_buffer, ifield, &
                              diag_field%get_mask(), diag_field%get_weight(), &
                              bounds, .False., Time=this%current_model_time)
        if (trim(error_string) .ne. "") call mpp_error(FATAL, "Field:"//trim(diag_field%get_varname()//&
                                                       " -"//trim(error_string)))
      else
        call mpp_error(FATAL, "diag_send_complete:: no input buffer allocated for field"//diag_field%get_longname())
      endif has_input_buff
    endif doing_math
    call diag_field%set_math_needs_to_be_done(.False.)
    !> Clean up, clean up, everybody do your share
    if (allocated(file_ids)) deallocate(file_ids)
    if (associated(diag_field)) nullify(diag_field)
  enddo field_loop

call this%fms_diag_do_io()
#endif

end subroutine fms_diag_send_complete

!> @brief Loops through all the files, open the file, writes out axis and
!! variable metadata and data when necessary.
!! TODO: passing in the saved mask from the field obj to diag_reduction_done_wrapper
!! for performance
subroutine fms_diag_do_io(this, end_time)
  class(fmsDiagObject_type), target, intent(inout)  :: this          !< The diag object
  type(time_type), optional, target, intent(in)     :: end_time      !< the model end_time
#ifdef use_yaml
  integer :: i !< For do loops
  class(fmsDiagFileContainer_type), pointer :: diag_file !< Pointer to this%FMS_diag_files(i) (for convenience)
  class(fmsDiagOutputBuffer_type), pointer  :: diag_buff !< pointer to output buffers iterated in buff_loop
  class(fmsDiagField_type), pointer         :: diag_field !< pointer to output buffers iterated in buff_loop
  class(DiagYamlFilesVar_type), pointer     :: field_yaml !< Pointer to a field from yaml fields
  TYPE (time_type),                 pointer :: model_time!< The current model time
  integer, allocatable                      :: buff_ids(:) !< ids for output buffers to loop through
  integer                                   :: ibuff !< buffer index
  logical :: file_is_opened_this_time_step !< True if the file was opened in this time_step
                                           !! If true the metadata will need to be written
  logical :: force_write !< force the last write if at end of run
  logical :: finish_writing !< true if finished writing for all the fields
  logical :: has_mask !< whether we have a mask
  logical, parameter :: DEBUG_REDUCT = .false. !< enables debugging output
  class(*), allocatable :: missing_val !< netcdf missing value for a given field
  real(r8_kind) :: mval !< r8 copy of missing value
  character(len=128) :: error_string !< outputted error string from reducti

  force_write = .false.
  if (present (end_time)) then
    force_write = .true.
    model_time => end_time
  else
    model_time => this%current_model_time
  endif

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
      call diag_file%increase_unlim_dimension_level()
    endif

    finish_writing = diag_file%is_time_to_write(model_time)

    ! finish reduction method if its time to write
    buff_ids = diag_file%FMS_diag_file%get_buffer_ids()
    ! loop through the buffers and finish reduction if needed
    buff_loop: do ibuff=1, SIZE(buff_ids)
      diag_buff => this%FMS_diag_output_buffers(buff_ids(ibuff))
      field_yaml => diag_yaml%diag_fields(diag_buff%get_yaml_id())
      diag_field => this%FMS_diag_fields(diag_buff%get_field_id())

      ! Go away if there is no data to write
      if (.not. diag_buff%is_there_data_to_write()) cycle

      if ( diag_buff%is_time_to_finish_reduction(end_time)) then
        ! sets missing value
        mval = diag_field%find_missing_value(missing_val)
        ! time_average and greater values all involve averaging so need to be "finished" before written
        if( field_yaml%has_var_reduction()) then
          if( field_yaml%get_var_reduction() .ge. time_average) then
            if(DEBUG_REDUCT)call mpp_error(NOTE, "fms_diag_do_io:: finishing reduction for "//diag_field%get_longname())
            error_string = diag_buff%diag_reduction_done_wrapper( &
                                    field_yaml%get_var_reduction(), &
                                   mval, diag_field%get_var_is_masked(), diag_field%get_mask_variant())
          endif
        endif
        call diag_file%write_field_data(diag_field, diag_buff)
        call diag_buff%set_next_output(diag_file%get_next_next_output())
      endif
      nullify(diag_buff)
      nullify(field_yaml)
    enddo buff_loop
    deallocate(buff_ids)

    if (finish_writing) then
      call diag_file%write_time_data()
      call diag_file%update_next_write(model_time)
      call diag_file%update_current_new_file_freq_index(model_time)
      call diag_file%increase_unlim_dimension_level()
      if (diag_file%is_time_to_close_file(model_time)) call diag_file%close_diag_file()
    else if (force_write) then
      call diag_file%write_time_data(is_the_end = .true.)
      call diag_file%close_diag_file()
    endif
  enddo
#endif
end subroutine fms_diag_do_io

!> @brief Computes average, min, max, rms error, etc.
!! based on the specified reduction method for the field.
!> @return Empty string if successful, error message if it fails
function fms_diag_do_reduction(this, field_data, diag_field_id, oor_mask, weight, &
  bounds, using_blocking, time) &
  result(error_msg)
  class(fmsDiagObject_type), intent(inout), target:: this                !< Diag Object
  class(*),                  intent(in)           :: field_data(:,:,:,:) !< Field data
  integer,                   intent(in)           :: diag_field_id       !< ID of the input field
  logical,                   intent(in), target   :: oor_mask(:,:,:,:)   !< mask
  real(kind=r8_kind),        intent(in)           :: weight              !< Must be a updated weight
  type(fmsDiagIbounds_type), intent(in)           :: bounds              !< Bounds for the field
  logical,                   intent(in)           :: using_blocking      !< .True. if field data is passed
                                                                         !! in blocks
  type(time_type),           intent(in), optional :: time                !< Current time

  character(len=50)         :: error_msg          !< Error message to check
  !TODO Mostly everything
#ifdef use_yaml
  type(fmsDiagField_type),          pointer :: field_ptr      !< Pointer to the field's object
  type(fmsDiagOutputBuffer_type),   pointer :: buffer_ptr     !< Pointer to the field's buffer
  class(fmsDiagFileContainer_type), pointer :: file_ptr       !< Pointer to the field's file
  type(diagYamlFilesVar_type),      pointer :: field_yaml_ptr !< Pointer to the field's yaml

  integer                   :: reduction_method   !< Integer representing a reduction method
  integer                   :: ids                !< For looping through buffer ids
  integer                   :: buffer_id          !< Id of the buffer
  integer                   :: file_id            !< File id
  integer, pointer          :: axis_ids(:)        !< Axis ids for the buffer
  logical                   :: is_subregional     !< .True. if the buffer is subregional
  logical                   :: reduced_k_range    !< .True. is the field is only outputing a section
                                                  !! of the z dimension
  type(fmsDiagIbounds_type) :: bounds_in          !< Starting and ending indices of the input field_data
  type(fmsDiagIbounds_type) :: bounds_out         !< Starting and ending indices of the output buffer
  integer                   :: i                  !< For looping through axid ids
  integer                   :: sindex             !< Starting index of a subregion
  integer                   :: eindex             !< Ending index of a subregion
  integer                   :: compute_idx(2)     !< Starting and Ending of the compute domain
  character(len=1)          :: cart_axis          !< Cartesian axis of the axis
  logical                   :: block_in_subregion !< .True. if the current block is part of the subregion
  integer                   :: starting           !< Starting index of the subregion relative to the compute domain
  integer                   :: ending             !< Ending index of the subregion relative to the compute domain
  real(kind=r8_kind)        :: missing_value      !< Missing_value for data points that are masked
                                                  !! This will obtained as r8 and converted to the right type as
                                                  !! needed. This is to avoid yet another select type ...
  logical                   :: new_time           !< .True. if this is a new time (i.e data has not be been
                                                  !! sent for this time)

  !TODO mostly everything
  field_ptr => this%FMS_diag_fields(diag_field_id)
  if (field_ptr%has_missing_value()) then
    select type (missing_val => field_ptr%get_missing_value(r8))
    type is (real(kind=r8_kind))
      missing_value = missing_val
    class default
      call mpp_error(FATAl, "The missing value for the field:"//trim(field_ptr%get_varname())//&
        &" was not allocated to the correct type. This shouldn't have happened")
    end select
  else
    select type (missing_val => get_default_missing_value(r8))
    type is (real(kind=r8_kind))
      missing_value = missing_val
    class default
      call mpp_error(FATAl, "The missing value for the field:"//trim(field_ptr%get_varname())//&
        &" was not allocated to the correct type. This shouldn't have happened")
    end select
  endif

  buffer_loop: do ids = 1, size(field_ptr%buffer_ids)
    error_msg = ""
    buffer_id = this%FMS_diag_fields(diag_field_id)%buffer_ids(ids)
    file_id = this%FMS_diag_fields(diag_field_id)%file_ids(ids)

    !< Gather all the objects needed for the buffer
    field_yaml_ptr => field_ptr%diag_field(ids)
    buffer_ptr     => this%FMS_diag_output_buffers(buffer_id)
    file_ptr       => this%FMS_diag_files(file_id)

    !< Go away if the file is a subregional file and the current PE does not have any data for it
    if (.not. file_ptr%writing_on_this_pe()) cycle

    !< Go away if finished doing math for this buffer
    if (buffer_ptr%is_done_with_math()) cycle

    bounds_out = bounds
    if (.not. using_blocking) then
      !< Set output bounds to start at 1:size(buffer_ptr%buffer)
      call bounds_out%reset_bounds_from_array_4D(buffer_ptr%buffer(:,:,:,:,1))
    endif

    bounds_in = bounds
    if (.not. bounds%has_halos) then
      !< If field_data does not contain halos, set bounds_in to start at 1:size(field_data)
      call bounds_in%reset_bounds_from_array_4D(field_data)
    endif

    is_subregional = file_ptr%is_regional()
    reduced_k_range = field_yaml_ptr%has_var_zbounds()

    !< Reset the bounds based on the reduced k range and subregional
    is_subregional_reduced_k_range: if (is_subregional .or. reduced_k_range) then
      call buffer_ptr%get_axis_ids(axis_ids)
      block_in_subregion = .true.
      axis_loops: do i = 1, size(axis_ids)
        !< Move on if the block does not have any data for the subregion
        if (.not. block_in_subregion) cycle

        select type (diag_axis => this%diag_axis(axis_ids(i))%axis)
        type is (fmsDiagSubAxis_type)
          sindex = diag_axis%get_starting_index()
          eindex = diag_axis%get_ending_index()
          compute_idx = diag_axis%get_compute_indices()
          starting=sindex-compute_idx(1)+1
          ending=eindex-compute_idx(1)+1
          if (using_blocking) then
            block_in_subregion = determine_if_block_is_in_region(starting, ending, bounds, i)
            if (.not. block_in_subregion) cycle

            !< Set bounds_in so that you can the correct section of the data for the block (starting at 1)
            call bounds_in%rebase_input(bounds, starting, ending, i)

            !< Set bounds_out to be the correct section relative to the block starting and ending indices
            call bounds_out%rebase_output(starting, ending, i)
          else
            !< Set bounds_in so that only the subregion section of the data will be used (starting at 1)
            call bounds_in%update_index(starting, ending, i, .false.)

            !< Set bounds_out to 1:size(subregion) for the PE
            call bounds_out%update_index(1, ending-starting+1, i, .true.)
          endif
        end select
      enddo axis_loops
      nullify(axis_ids)
      !< Move on to the next buffer if the block does not have any data for the subregion
      if (.not. block_in_subregion) cycle
    endif is_subregional_reduced_k_range

    !< Determine the reduction method for the buffer
    reduction_method = field_yaml_ptr%get_var_reduction()
    if (present(time)) new_time = buffer_ptr%update_buffer_time(time)
    call buffer_ptr%set_send_data_called()
    select case(reduction_method)
    case (time_none)
      error_msg = buffer_ptr%do_time_none_wrapper(field_data, oor_mask, field_ptr%get_var_is_masked(), &
        bounds_in, bounds_out, missing_value)
      if (trim(error_msg) .ne. "") then
        return
      endif
    case (time_min)
      error_msg = buffer_ptr%do_time_min_wrapper(field_data, oor_mask, field_ptr%get_var_is_masked(), &
        bounds_in, bounds_out, missing_value)
      if (trim(error_msg) .ne. "") then
        return
      endif
    case (time_max)
      error_msg = buffer_ptr%do_time_max_wrapper(field_data, oor_mask, field_ptr%get_var_is_masked(), &
        bounds_in, bounds_out, missing_value)
      if (trim(error_msg) .ne. "") then
        return
      endif
    case (time_sum)
      error_msg = buffer_ptr%do_time_sum_wrapper(field_data, oor_mask, field_ptr%get_var_is_masked(), &
        field_ptr%get_mask_variant(), bounds_in, bounds_out, missing_value, new_time)
      if (trim(error_msg) .ne. "") then
        return
      endif
    case (time_average)
      error_msg = buffer_ptr%do_time_sum_wrapper(field_data, oor_mask, field_ptr%get_var_is_masked(), &
        field_ptr%get_mask_variant(), bounds_in, bounds_out, missing_value, new_time)
      if (trim(error_msg) .ne. "") then
        return
      endif
    case (time_power)
      error_msg = buffer_ptr%do_time_sum_wrapper(field_data, oor_mask, field_ptr%get_var_is_masked(), &
        field_ptr%get_mask_variant(), bounds_in, bounds_out, missing_value, new_time, &
        pow_value=field_yaml_ptr%get_pow_value())
      if (trim(error_msg) .ne. "") then
        return
      endif
    case (time_rms)
      error_msg = buffer_ptr%do_time_sum_wrapper(field_data, oor_mask, field_ptr%get_var_is_masked(), &
        field_ptr%get_mask_variant(), bounds_in, bounds_out, missing_value, new_time, pow_value = 2)
      if (trim(error_msg) .ne. "") then
        return
      endif
    case (time_diurnal)
      if(.not. present(time)) call mpp_error(FATAL, &
                            "fms_diag_do_reduction:: time must be present when using diurnal reductions")
      ! sets the diurnal index for reduction within the buffer object
      call buffer_ptr%set_diurnal_section_index(time)
      error_msg = buffer_ptr%do_time_sum_wrapper(field_data, oor_mask, field_ptr%get_var_is_masked(), &
        field_ptr%get_mask_variant(), bounds_in, bounds_out, missing_value, new_time)
      if (trim(error_msg) .ne. "") then
        return
      endif
    case default
      error_msg = "The reduction method is not supported. "//&
        "Only none, min, max, sum, average, power, rms, and diurnal are supported."
    end select

    if (field_ptr%is_static() .or. file_ptr%FMS_diag_file%is_done_writing_data()) then
      call buffer_ptr%set_done_with_math()
    endif
  enddo buffer_loop
#else
  error_msg = ""
  CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#endif
end function fms_diag_do_reduction

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
  character(len=20) :: set_name      !< Name of the axis set
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
      set_name = ""
      if (axis%has_set_name()) set_name = axis%get_set_name()
      do j = 1, size(axis_names)
        uncmx_ids(j) = get_axis_id_from_name(axis_names(j), this%diag_axis, this%registered_axis, set_name)
        if (uncmx_ids(j) .eq. diag_null) call mpp_error(FATAL, &
          &"Error parsing the compress attribute for axis: "//trim(axis%get_axis_name())//&
          &". Be sure that the axes in the compress attribute are registered")
      enddo
      call axis%add_structured_axis_ids(uncmx_ids)
    endif
  end select
#endif
end subroutine fms_diag_axis_add_attribute

!> \brief Gets the field_name from the diag_field
!> \returns a copy of the field_name
function fms_get_field_name_from_id (this, field_id) &
  result(field_name)

  class(fmsDiagObject_type), intent (in) :: this     !< The diag object, the caller
  integer,                   intent (in) :: field_id !< Field id to get the name for
  character(len=:), allocatable :: field_name
#ifndef use_yaml
  CALL MPP_ERROR(FATAL,"You can not use the modern diag manager without compiling with -Duse_yaml")
#else
  field_name = this%FMS_diag_fields(field_id)%get_varname()
#endif
end function fms_get_field_name_from_id

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
  type(fmsDiagOutputBuffer_type),allocatable:: rslt
  if( (bufferid .gt. UBOUND(this%FMS_diag_output_buffers, 1)) .or. &
      (bufferid .lt. LBOUND(this%FMS_diag_output_buffers, 1))) &
    call mpp_error(FATAL, 'get_diag_bufer: invalid bufferid given')
  rslt = this%FMS_diag_output_buffers(bufferid)
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
  integer :: axes_length(4) !< Length of each axis
  integer :: i, j !< For looping
  class(fmsDiagOutputBuffer_type), pointer :: ptr_diag_buffer_obj !< Pointer to the buffer class
  class(DiagYamlFilesVar_type), pointer :: ptr_diag_field_yaml !< Pointer to a field from yaml fields
  integer, pointer :: axis_ids(:) !< Pointer to indices of axes of the field variable
  integer :: var_type !< Stores type of the field data (r4, r8, i4, i8, and string) represented as an integer.
  character(len=:), allocatable :: var_name !< Field name to initialize output buffers
  logical :: is_scalar !< Flag indicating that the variable is a scalar
  integer :: yaml_id !< Yaml id for the buffer
  integer :: file_id !< File id for the buffer

  if (this%FMS_diag_fields(field_id)%buffer_allocated) return

  ! Determine the type of the field data
  var_type = get_var_type(field_data(1, 1, 1, 1))

  ! Get variable/field name
  var_name = this%FMS_diag_fields(field_id)%get_varname()

  ! Determine dimensions of the field
  is_scalar = this%FMS_diag_fields(field_id)%is_scalar()

  ! Loop over a number of fields/buffers where this variable occurs
  do i = 1, size(this%FMS_diag_fields(field_id)%buffer_ids)
    buffer_id = this%FMS_diag_fields(field_id)%buffer_ids(i)
    file_id = this%FMS_diag_fields(field_id)%file_ids(i)

    !< Go away if the file is a subregional file and the current PE does not have any data for it
    if (.not. this%FMS_diag_files(file_id)%writing_on_this_pe()) cycle

    ndims = 0
    if (.not. is_scalar) then
      call this%FMS_diag_output_buffers(buffer_id)%get_axis_ids(axis_ids)
      ndims = size(axis_ids)
    endif

    yaml_id = this%FMS_diag_output_buffers(buffer_id)%get_yaml_id()

    ptr_diag_field_yaml => diag_yaml%diag_fields(yaml_id)
    num_diurnal_samples = ptr_diag_field_yaml%get_n_diurnal() !< Get number of diurnal samples

    axes_length = 1
    do j = 1, ndims
      axes_length(j) = this%fms_get_axis_length(axis_ids(j))
    enddo

    if (num_diurnal_samples .ne. 0) then
      ndims = ndims + 1 !< Add one more dimension for the diurnal axis
    endif

    ptr_diag_buffer_obj => this%FMS_diag_output_buffers(buffer_id)
    call ptr_diag_buffer_obj%allocate_buffer(field_data(1, 1, 1, 1), ndims, axes_length(1:4), &
      this%FMS_diag_fields(field_id)%get_mask_variant(), var_name, num_diurnal_samples)
    call ptr_diag_buffer_obj%initialize_buffer(ptr_diag_field_yaml%get_var_reduction(), var_name)

  enddo
  nullify(axis_ids)

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

!> @brief Update the current model time in the diag object
subroutine update_current_model_time(this, time)
  class(fmsDiagObject_type), intent(inout) :: this !< Diag Object
  type(time_type),           intent(in)    :: time !< Current diag manager time
#ifdef use_yaml
  if(time > this%current_model_time) this%current_model_time = time
#endif
end subroutine update_current_model_time

end module fms_diag_object_mod
