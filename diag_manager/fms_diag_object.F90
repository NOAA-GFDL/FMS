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
use mpp_mod, only: fatal, note, warning, mpp_error
use diag_data_mod,  only: diag_null, diag_not_found, diag_not_registered, diag_registered_id, &
                         &DIAG_FIELD_NOT_FOUND, diag_not_registered
  USE time_manager_mod, ONLY: set_time, set_date, get_time, time_type, OPERATOR(>=), OPERATOR(>),&
       & OPERATOR(<), OPERATOR(==), OPERATOR(/=), OPERATOR(/), OPERATOR(+), ASSIGNMENT(=), get_date, &
       & get_ticks_per_second
#ifdef use_yaml
use fms_diag_file_object_mod, only: fmsDiagFileContainer_type, fmsDiagFile_type, fms_diag_files_object_init
use fms_diag_field_object_mod, only: fmsDiagField_type, fms_diag_fields_object_init
use fms_diag_yaml_mod, only: diag_yaml_object_init, find_diag_field, get_diag_files_id
use fms_diag_axis_object_mod, only: fms_diag_axis_object_init
#endif
implicit none
private

type fmsDiagObject_type
!TODO add container arrays
#ifdef use_yaml
private
!TODO: Remove FMS prefix from variables in this type 
  class(fmsDiagFileContainer_type), allocatable :: FMS_diag_files (:) !< array of diag files
  class(fmsDiagField_type), allocatable :: FMS_diag_fields(:) !< Array of diag fields
  integer, private :: registered_variables !< Number of registered variables
  logical, private :: initialized=.false. !< True if the fmsDiagObject is initialized
  logical, private :: files_initialized=.false. !< True if the fmsDiagObject is initialized
  logical, private :: fields_initialized=.false. !< True if the fmsDiagObject is initialized
  logical, private :: buffers_initialized=.false. !< True if the fmsDiagObject is initialized
  logical, private :: axes_initialized=.false. !< True if the fmsDiagObject is initialized
#endif
  contains
    procedure :: init => fms_diag_object_init
    procedure :: fms_register_diag_field_scalar
    procedure :: fms_register_diag_field_array
    procedure :: fms_register_static_field
    procedure :: register => fms_register_diag_field_obj !! Merely initialize fields.
    procedure :: fms_diag_field_add_attribute
    procedure :: fms_get_diag_field_id_from_name
    procedure :: diag_end => fms_diag_object_end
end type fmsDiagObject_type

type (fmsDiagObject_type), target :: fms_diag_object
integer, private :: registered_variables !< Number of registered variables
public :: fms_register_diag_field_obj
public :: fms_register_diag_field_scalar
public :: fms_register_diag_field_array 
public :: fms_register_static_field
public :: fms_diag_field_add_attribute
public :: fms_get_diag_field_id_from_name
public :: fms_diag_object 
public :: fmsDiagObject_type

contains

!> @brief Initiliazes the  fms_diag_object.
!! Reads the diag_table.yaml and fills in the yaml object
!! Allocates the diag manager object arrays for files, fields, and buffers
!! Initializes variables
subroutine fms_diag_object_init (obj,diag_subset_output)
 class(fmsDiagObject_type) :: obj !< Diag mediator/controller object
 integer :: diag_subset_output !< Subset of the diag output?
#ifdef use_yaml
 if (obj%initialized) return

!TODO: allocate the file, field, and buffer containers
! allocate(diag_objs(get_num_unique_fields()))
  CALL diag_yaml_object_init(diag_subset_output)
  CALL fms_diag_axis_object_init()
  obj%files_initialized = fms_diag_files_object_init(obj%FMS_diag_files)
  obj%fields_initialized = fms_diag_fields_object_init (obj%FMS_diag_fields)
 registered_variables = 0
 obj%initialized = .true.
#else
  call mpp_error("fms_diag_object_init",&
    "You must compile with -Duse_yaml to use the option use_modern_diag", FATAL)
#endif
end subroutine fms_diag_object_init
!> \description Loops through all files and does one final write.
!! Closes all files
!! Deallocates all buffers, fields, and files
!! Uninitializes the fms_diag_object
subroutine fms_diag_object_end (obj)
  class(fmsDiagObject_type) :: obj
#ifdef use_yaml
  !TODO: loop through files and force write
  !TODO: Close all files
  !TODO: Deallocate diag object arrays and clean up all memory
  obj%initialized = .false.
#endif
end subroutine fms_diag_object_end
!> \Description Fills in and allocates (when necessary) the values in the diagnostic object
subroutine fms_register_diag_field_obj &
                !(field_obj, modname, varname, axes, time, longname, units, missing_value, metadata)
       (fms_diag_object, modname, varname, diag_field_indices, axes, init_time, &
       longname, units, missing_value, varRange, mask_variant, standname, &
       do_not_log, err_msg, interp_method, tile_count, area, volume, realm)

 class(fmsDiagObject_type),TARGET,INTENT(inout):: fms_diag_object       !< Diaj_obj to fill
 CHARACTER(len=*),               INTENT(in)    :: modname               !< The module name
 CHARACTER(len=*),               INTENT(in)    :: varname               !< The variable name
 integer,                        INTENT(in)    :: diag_field_indices(:) !< Array of indices to the field
                                                                        !! in the yaml object
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
#ifdef use_yaml

 class (fmsDiagFile_type), pointer :: fileptr => null()
 class (fmsDiagField_type), pointer :: fieldptr => null()
 integer, allocatable :: file_ids(:) !< The file IDs for this variable
 integer :: i !< For do loops
 integer :: j !< fms_diag_object%FMS_diag_fields%file_ids(i) (for less typing :)
  
!> Use pointers for convenience
  fieldptr => fms_diag_object%FMS_diag_fields(registered_variables)
!> Register the data for the field
  call fieldptr%register(modname, varname, diag_field_indices, &
       axes, init_time, longname, units, missing_value, varRange, mask_variant, standname, &
       do_not_log, err_msg, interp_method, tile_count, area, volume, realm)
!> Get the file IDs from the field indicies from the yaml
  file_ids = get_diag_files_id(diag_field_indices)
!> Add the axis information, initial time, and field IDs to the files
  if (present(axes) .and. present(init_time)) then
    do i = 1, size(file_ids)
     fileptr => fms_diag_object%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_id(fieldptr%get_id())
     call fileptr%set_domain_from_axis(axes)
     call fileptr%add_axes(axes)
     call fileptr%add_start_time(init_time)
    enddo
  elseif (present(axes)) then !only axes present
    do i = 1, size(file_ids)
     fileptr => fms_diag_object%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_id(fieldptr%get_id())
     call fileptr%set_domain_from_axis(axes)
     call fileptr%add_axes(axes)
    enddo
  elseif (present(init_time)) then !only inti time present
    do i = 1, size(file_ids)
     fileptr => fms_diag_object%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_id(fieldptr%get_id())
     call fileptr%add_start_time(init_time)
    enddo
  else !no axis or init time present
    do i = 1, size(file_ids)
     fileptr => fms_diag_object%FMS_diag_files(file_ids(i))%FMS_diag_file
     call fileptr%add_field_id(fieldptr%get_id())
    enddo
  endif
  nullify (fileptr)
  nullify (fieldptr)
#endif
end subroutine fms_register_diag_field_obj

  !> @brief Registers a scalar field
  !! @return field index for subsequent call to send_data.
INTEGER FUNCTION fms_register_diag_field_scalar(fms_diag_object,module_name, field_name, init_time, &
       & long_name, units, missing_value, var_range, standard_name, do_not_log, err_msg,&
       & area, volume, realm)
    class(fmsDiagObject_type),TARGET,INTENT(inout):: fms_diag_object       !< Diaj_obj to fill
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

#ifdef use_yaml
    integer, allocatable :: diag_field_indices(:) !< indices where the field was found

    diag_field_indices = find_diag_field(field_name, module_name)
    if (diag_field_indices(1) .eq. diag_null) then
      !< The field was not found in the table, so return diag_null
      fms_register_diag_field_scalar = diag_null
      deallocate(diag_field_indices)
      return
    endif

    registered_variables = registered_variables + 1
    fms_register_diag_field_scalar = registered_variables

    call fms_diag_object%FMS_diag_fields(registered_variables)%setID(registered_variables)
    call fms_diag_object%FMS_diag_fields(registered_variables)%register(&
      & module_name, field_name, diag_field_indices, init_time=init_time, &
      & longname=long_name, units=units, missing_value=missing_value, varrange=var_range, &
      & standname=standard_name, do_not_log=do_not_log, err_msg=err_msg, &
      & area=area, volume=volume, realm=realm)
    deallocate(diag_field_indices)
#else 
fms_register_diag_field_scalar = diag_not_registered
#endif
end function fms_register_diag_field_scalar

    !> @brief Registers an array field
  !> @return field index for subsequent call to send_data.
INTEGER FUNCTION fms_register_diag_field_array(fms_diag_object, module_name, field_name, axes, init_time, &
       & long_name, units, missing_value, var_range, mask_variant, standard_name, verbose,&
       & do_not_log, err_msg, interp_method, tile_count, area, volume, realm)
    class(fmsDiagObject_type),TARGET,INTENT(inout):: fms_diag_object       !< Diaj_obj to fill
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

#ifdef use_yaml
    integer, allocatable :: diag_field_indices(:) !< indices of diag_field yaml where the field was found

    diag_field_indices = find_diag_field(field_name, module_name)
    if (diag_field_indices(1) .eq. diag_null) then
      !< The field was not found in the table, so return diag_null
      fms_register_diag_field_array = diag_null
      deallocate(diag_field_indices)
      return
    endif

    registered_variables = registered_variables + 1
    fms_register_diag_field_array = registered_variables

    call fms_diag_object%FMS_diag_fields(registered_variables)%setID (registered_variables)
    call fms_diag_object%FMS_diag_fields(registered_variables)%register( &
      & module_name, field_name, diag_field_indices, init_time=init_time, &
      & axes=axes, longname=long_name, units=units, missing_value=missing_value, varrange=var_range, &
      & mask_variant=mask_variant, standname=standard_name, do_not_log=do_not_log, err_msg=err_msg, &
      & interp_method=interp_method, tile_count=tile_count, area=area, volume=volume, realm=realm)
    deallocate(diag_field_indices)
#else
fms_register_diag_field_array = diag_not_registered
#endif
end function fms_register_diag_field_array

!> @brief Return field index for subsequent call to send_data.
!! @return field index for subsequent call to send_data.
INTEGER FUNCTION fms_register_static_field(fms_diag_object, module_name, field_name, axes, long_name, units,&
       & missing_value, range, mask_variant, standard_name, DYNAMIC, do_not_log, interp_method,&
       & tile_count, area, volume, realm)
    class(fmsDiagObject_type),TARGET,INTENT(inout):: fms_diag_object       !< Diaj_obj to fill
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

#ifdef use_yaml
    integer, allocatable :: diag_field_indices(:) !< indices where the field was foun

    diag_field_indices = find_diag_field(field_name, module_name)
    if (diag_field_indices(1) .eq. diag_null) then
      !< The field was not found in the table, so return diag_null
      fms_register_static_field = diag_null
      deallocate(diag_field_indices)
      return
    endif

    registered_variables = registered_variables + 1
    fms_register_static_field = registered_variables

    call fms_diag_object%FMS_diag_fields(registered_variables)%setID(registered_variables)
! Include static as optional variable to register here
    call fms_diag_object%FMS_diag_fields(registered_variables)%register( &
      & module_name, field_name, diag_field_indices, axes=axes, &
      & longname=long_name, units=units, missing_value=missing_value, varrange=range, &
      & standname=standard_name, do_not_log=do_not_log, area=area, volume=volume, realm=realm, &
      & static=.true.)
    deallocate(diag_field_indices)
#else
fms_register_static_field = diag_not_registered
#endif
end function fms_register_static_field

!> @brief Add a attribute to the diag_obj using the diag_field_id
subroutine fms_diag_field_add_attribute(fms_diag_object, diag_field_id, att_name, att_value)
  class(fmsDiagObject_type), intent (inout) :: fms_diag_object !< The diag object
  integer,          intent(in) :: diag_field_id      !< Id of the axis to add the attribute to
  character(len=*), intent(in) :: att_name     !< Name of the attribute
  class(*),         intent(in) :: att_value(:) !< The attribute value to add
#ifdef use_yaml
!TODO: Value for diag not found
  if ( diag_field_id .LE. 0 ) THEN 
    RETURN 
  else
    if (fms_diag_object%FMS_diag_fields(diag_field_id)%is_registered() ) &
      call fms_diag_object%FMS_diag_fields(diag_field_id)%add_attribute(att_name, att_value)
  endif
#endif
end subroutine fms_diag_field_add_attribute
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
#ifdef use_yaml
!> Loop through fields to find it.
  if (registered_variables < 1) return
  do i=1,registered_variables
   diag_field_id = fms_diag_object%FMS_diag_fields(i)%id_from_name(module_name, field_name)
   if(diag_field_id .ne. DIAG_FIELD_NOT_FOUND) return
  enddo
#endif
END FUNCTION fms_get_diag_field_id_from_name
end module fms_diag_object_mod
