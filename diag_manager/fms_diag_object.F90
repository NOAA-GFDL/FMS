module fms_diag_object_mod
use fms_diag_file_object_mod :: FMS_diag_files
use fms_diag_field_object_mod :: fmsDiagField_type

type fmsDiagObject_type
!TODO add container arrays
  TYPE(fmsDiagField_type), private, ALLOCATABLE, target :: diag_fields(:) !< Array of diag objects
                                                                       !! one for each registered variable
  integer, private :: registered_variables !< Number of registered variables
  contains
    procedure :: register => fms_register_diag_field_obj !! Merely initialize fields.
end type fmsDiagObject_type

type (fmsDiagObject_type), target :: fms_diag_object

public :: fms_register_diag_field_obj
public :: fms_register_diag_field_scalar
public :: fms_register_diag_field_array 
public :: fms_register_static_field
public :: fms_diag_object 
contains
!> \Description Fills in and allocates (when necessary) the values in the diagnostic object
subroutine fms_register_diag_field_obj &
                !(field_obj, modname, varname, axes, time, longname, units, missing_value, metadata)
       (fms_diag_object, modname, varname, diag_field_indices, axes, init_time, &
       longname, units, missing_value, varRange, mask_variant, standname, &
       do_not_log, err_msg, interp_method, tile_count, area, volume, realm)

 class(fmsDiagObject_type),      INTENT(inout) :: fms_diag_object       !< Diaj_obj to fill
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

 integer :: i !< For do loops
 integer :: j !< fms_diag_object%field_obj%file_ids(i) (for less typing :)

#ifdef use_yaml
!> Fill in information from the register call
  fms_diag_object%field_obj%varname = trim(varname)
  fms_diag_object%field_obj%modname = trim(modname)

!> Fill in diag_field and find the ids of the files that this variable is in
  fms_diag_object%field_obj%diag_field = get_diag_fields_entries(diag_field_indices)
  fms_diag_object%field_obj%file_ids   = get_diag_files_id(diag_field_indices)

  if (present(axes)) then
    fms_diag_object%field_obj%axis_ids => axes
    call get_domain_and_domain_type(fms_diag_object%field_obj%axis_ids, fms_diag_object%field_obj%type_of_domain, fms_diag_object%field_obj%domain, fms_diag_object%field_obj%varname)
    do i = 1, size(fms_diag_object%field_obj%file_ids)
       j = fms_diag_object%field_obj%file_ids(i)
       call FMS_diag_files(j)%set_file_domain(fms_diag_object%field_obj%domain, fms_diag_object%field_obj%type_of_domain)
       call FMS_diag_files(j)%add_axes(axes)
       if (present(init_time)) call FMS_diag_files(j)%add_start_time(init_time)
    enddo
     !> TO DO:
     !!     Mark the field as registered in the diag_files
  else
     !> The variable is a scalar
    fms_diag_object%field_obj%type_of_domain = NO_DOMAIN
    fms_diag_object%field_obj%domain => null()
  endif

!> get the optional arguments if included and the diagnostic is in the diag table
  if (present(longname))      fms_diag_object%field_obj%longname      = trim(longname)
  if (present(standname))     fms_diag_object%field_obj%standname     = trim(standname)
  if (present(units))         fms_diag_object%field_obj%units         = trim(units)
  if (present(realm))         fms_diag_object%field_obj%realm         = trim(realm)
  if (present(interp_method)) fms_diag_object%field_obj%interp_method = trim(interp_method)
  if (present(tile_count)) then
    allocate(fms_diag_object%field_obj%tile_count)
    fms_diag_object%field_obj%tile_count = tile_count
  endif

  if (present(missing_value)) then
    select type (missing_value)
     type is (integer(kind=i4_kind))
             allocate(integer(kind=i4_kind) :: fms_diag_object%field_obj%missing_value)
             fms_diag_object%field_obj%missing_value = missing_value
     type is (integer(kind=i8_kind))
             allocate(integer(kind=i8_kind) :: fms_diag_object%field_obj%missing_value)
             fms_diag_object%field_obj%missing_value = missing_value
     type is (real(kind=r4_kind))
             allocate(integer(kind=r4_kind) :: fms_diag_object%field_obj%missing_value)
             fms_diag_object%field_obj%missing_value = missing_value
     type is (real(kind=r8_kind))
             allocate(integer(kind=r8_kind) :: fms_diag_object%field_obj%missing_value)
             fms_diag_object%field_obj%missing_value = missing_value
     class default
             call mpp_error("fms_register_diag_field_obj", &
                     "The missing value passed to register a diagnostic is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
  else
      allocate(real :: fms_diag_object%field_obj%missing_value)
      select type (miss => fms_diag_object%field_obj%missing_value)
       type is (real)
        miss = real(CMOR_MISSING_VALUE)
      end select
  endif

  if (present(varRANGE)) then
    select type (varRANGE)
     type is (integer(kind=i4_kind))
             allocate(integer(kind=i4_kind) :: fms_diag_object%field_obj%data_RANGE(2))
             fms_diag_object%field_obj%data_RANGE = varRANGE
     type is (integer(kind=i8_kind))
             allocate(integer(kind=i8_kind) :: fms_diag_object%field_obj%data_RANGE(2))
             fms_diag_object%field_obj%data_RANGE = varRANGE
     type is (real(kind=r4_kind))
             allocate(integer(kind=r4_kind) :: fms_diag_object%field_obj%data_RANGE(2))
             fms_diag_object%field_obj%data_RANGE = varRANGE
     type is (real(kind=r8_kind))
             allocate(integer(kind=r8_kind) :: fms_diag_object%field_obj%data_RANGE(2))
             fms_diag_object%field_obj%data_RANGE = varRANGE
     class default
             call mpp_error("fms_register_diag_field_obj", &
                     "The varRange passed to register a diagnostic is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
  else
      allocate(real :: fms_diag_object%field_obj%data_RANGE(2))
      select type (varRANGE => fms_diag_object%field_obj%data_RANGE)
       type is (real)
        varRANGE = real(CMOR_MISSING_VALUE)
      end select
  endif

  if (present(area)) then
    if (area < 0) call mpp_error("fms_register_diag_field_obj", &
                     "The area id passed with field_name"//trim(varname)//" has not been registered."&
                     "Check that there is a register_diag_field call for the AREA measure and that is in the"&
                     "diag_table.yaml", FATAL)
    allocate(fms_diag_object%field_obj%area)
    fms_diag_object%field_obj%area = area
  endif

  if (present(volume)) then
    if (volume < 0) call mpp_error("fms_register_diag_field_obj", &
                     "The volume id passed with field_name"//trim(varname)//" has not been registered."&
                     "Check that there is a register_diag_field call for the VOLUME measure and that is in the"&
                     "diag_table.yaml", FATAL)
    allocate(fms_diag_object%field_obj%volume)
    fms_diag_object%field_obj%volume = volume
  endif

  if (present(mask_variant)) then
    allocate(fms_diag_object%field_obj%mask_variant)
    fms_diag_object%field_obj%mask_variant = mask_variant
  endif

  if (present(do_not_log)) then
    allocate(fms_diag_object%field_obj%do_not_log)
    fms_diag_object%field_obj%do_not_log = do_not_log
  endif

 !< Allocate space for any additional variable attributes
 !< These will be fill out when calling `diag_field_add_attribute`
 allocate(fms_diag_object%field_obj%attributes(max_field_attributes))
 fms_diag_object%field_obj%num_attributes = 0
 fms_diag_object%field_obj%registered = .true.
#endif
end subroutine fms_register_diag_field_obj

  !> @brief Registers a scalar field
  !! @return field index for subsequent call to send_data.
  INTEGER FUNCTION fms_register_diag_field_scalar(module_name, field_name, init_time, &
       & long_name, units, missing_value, var_range, standard_name, do_not_log, err_msg,&
       & area, volume, realm)
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

    call diag_objs(registered_variables)%setID(registered_variables)
    call diag_objs(registered_variables)%register(module_name, field_name, diag_field_indices, init_time=init_time, &
      & longname=long_name, units=units, missing_value=missing_value, varrange=var_range, &
      & standname=standard_name, do_not_log=do_not_log, err_msg=err_msg, &
      & area=area, volume=volume, realm=realm)
    deallocate(diag_field_indices)
#endif

  end function fms_register_diag_field_scalar

    !> @brief Registers an array field
  !> @return field index for subsequent call to send_data.
  INTEGER FUNCTION fms_register_diag_field_array(module_name, field_name, axes, init_time, &
       & long_name, units, missing_value, var_range, mask_variant, standard_name, verbose,&
       & do_not_log, err_msg, interp_method, tile_count, area, volume, realm)
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

    call diag_objs(registered_variables)%setID(registered_variables)
    call diag_objs(registered_variables)%register(module_name, field_name, diag_field_indices, init_time=init_time, &
      & axes=axes, longname=long_name, units=units, missing_value=missing_value, varrange=var_range, &
      & mask_variant=mask_variant, standname=standard_name, do_not_log=do_not_log, err_msg=err_msg, &
      & interp_method=interp_method, tile_count=tile_count, area=area, volume=volume, realm=realm)
    deallocate(diag_field_indices)
#endif

end function fms_register_diag_field_array

!> @brief Return field index for subsequent call to send_data.
!! @return field index for subsequent call to send_data.
INTEGER FUNCTION fms_register_static_field(module_name, field_name, axes, long_name, units,&
       & missing_value, range, mask_variant, standard_name, DYNAMIC, do_not_log, interp_method,&
       & tile_count, area, volume, realm)
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

    call diag_objs(registered_variables)%setID(registered_variables)
    allocate(diag_objs(registered_variables)%static)
    diag_objs(registered_variables)%static = .true.
    call diag_objs(registered_variables)%register(module_name, field_name, diag_field_indices, axes=axes, &
      & longname=long_name, units=units, missing_value=missing_value, varrange=range, &
      & standname=standard_name, do_not_log=do_not_log, area=area, volume=volume, realm=realm)
    deallocate(diag_field_indices)
#endif
end function fms_register_static_field


end fms_diag_object_mod
