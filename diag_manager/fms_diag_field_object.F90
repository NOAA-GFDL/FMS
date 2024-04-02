module fms_diag_field_object_mod
!> \author Tom Robinson
!> \email thomas.robinson@noaa.gov
!! \brief Contains routines for the diag_objects
!!
!! \description The diag_manager passes an object back and forth between the diag routines and the users.
!! The procedures of this object and the types are all in this module.  The fms_dag_object is a type
!! that contains all of the information of the variable.  It is extended by a type that holds the
!! appropriate buffer for the data for manipulation.
#ifdef use_yaml
use diag_data_mod,  only: prepend_date, diag_null, CMOR_MISSING_VALUE, diag_null_string, MAX_STR_LEN
use diag_data_mod,  only: r8, r4, i8, i4, string, null_type_int, NO_DOMAIN
use diag_data_mod,  only: max_field_attributes, fmsDiagAttribute_type
use diag_data_mod,  only: diag_null, diag_not_found, diag_not_registered, diag_registered_id, &
                         &DIAG_FIELD_NOT_FOUND, avg_name, time_average, time_min, time_max, &
                         &time_none, time_diurnal, time_power, time_rms, time_sum
use fms_string_utils_mod, only: int2str=>string
use mpp_mod, only: fatal, note, warning, mpp_error, mpp_pe, mpp_root_pe
use fms_diag_yaml_mod, only:  diagYamlFilesVar_type, get_diag_fields_entries, get_diag_files_id, &
  & find_diag_field, get_num_unique_fields, diag_yaml
use fms_diag_axis_object_mod, only: diagDomain_t, get_domain_and_domain_type, fmsDiagAxis_type, &
  & fmsDiagAxisContainer_type, fmsDiagFullAxis_Type
use time_manager_mod, ONLY: time_type, get_date
use fms2_io_mod, only: FmsNetcdfFile_t, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t, register_field, &
                       register_variable_attribute
use fms_diag_input_buffer_mod, only: fmsDiagInputBuffer_t
!!!set_time, set_date, get_time, time_type, OPERATOR(>=), OPERATOR(>),&
!!!       & OPERATOR(<), OPERATOR(==), OPERATOR(/=), OPERATOR(/), OPERATOR(+), ASSIGNMENT(=), get_date, &
!!!       & get_ticks_per_second

use platform_mod
use iso_c_binding

implicit none

private

!> \brief Object that holds all variable information
type fmsDiagField_type
     type (diagYamlFilesVar_type), allocatable, dimension(:) :: diag_field !< info from diag_table for this variable
     integer,                      allocatable, dimension(:) :: file_ids   !< Ids of the FMS_diag_files the variable
                                                                           !! belongs to
     integer, allocatable, private                    :: diag_id           !< unique id for varable
     integer, allocatable, dimension(:)               :: buffer_ids        !< index/id for this field's buffers
     type(fmsDiagAttribute_type), allocatable         :: attributes(:)     !< attributes for the variable
     integer,              private                    :: num_attributes    !< Number of attributes currently added
     logical, allocatable, private                    :: static            !< true if this is a static var
     logical, allocatable, private                    :: scalar            !< .True. if the variable is a scalar
     logical, allocatable, private                    :: registered        !< true when registered
     logical, allocatable, private                    :: mask_variant      !< true if the mask changes over time
     logical, allocatable, private                    :: var_is_masked     !< true if the field is masked
     logical, allocatable, private                    :: do_not_log        !< .true. if no need to log the diag_field
     logical, allocatable, private                    :: local             !< If the output is local
     integer,          allocatable, private           :: vartype           !< the type of varaible
     character(len=:), allocatable, private           :: varname           !< the name of the variable
     character(len=:), allocatable, private           :: longname          !< longname of the variable
     character(len=:), allocatable, private           :: standname         !< standard name of the variable
     character(len=:), allocatable, private           :: units             !< the units
     character(len=:), allocatable, private           :: modname           !< the module
     character(len=:), allocatable, private           :: realm             !< String to set as the value
                                                                           !! to the modeling_realm attribute
     character(len=:), allocatable, private           :: interp_method     !< The interp method to be used
                                                            !! when regridding the field in post-processing.
                                                            !! Valid options are "conserve_order1",
                                                            !! "conserve_order2", and "none".
     integer, allocatable, dimension(:), private      :: frequency         !< specifies the frequency
     integer, allocatable, private                    :: tile_count        !< The number of tiles
     integer, allocatable, dimension(:), private      :: axis_ids          !< variable axis IDs
     class(diagDomain_t), pointer,   private          :: domain            !< Domain
     INTEGER                         , private        :: type_of_domain    !< The type of domain ("NO_DOMAIN",
                                                                           !! "TWO_D_DOMAIN", or "UG_DOMAIN")
     integer, allocatable, private                    :: area, volume      !< The Area and Volume
     class(*), allocatable, private                   :: missing_value     !< The missing fill value
     class(*), allocatable, private                   :: data_RANGE(:)     !< The range of the variable data
     type(fmsDiagInputBuffer_t), allocatable          :: input_data_buffer !< Input buffer object for when buffering
                                                                           !! data
     logical, allocatable, private                    :: multiple_send_data!< .True. if send_data is called multiple
                                                                           !! times for the same model time
     logical, allocatable, private                    :: data_buffer_is_allocated !< True if the buffer has
                                                                           !! been allocated
     logical, allocatable, private                    :: math_needs_to_be_done !< If true, do math
                                                                           !! functions. False when done.
     logical, allocatable                             :: buffer_allocated  !< True if a buffer pointed by
                                                                           !! the corresponding index in
                                                                           !! buffer_ids(:) is allocated.
     logical, allocatable                             :: mask(:,:,:,:)     !< Mask passed in send_data
     logical                                          :: halo_present = .false. !< set if any halos are used
  contains
!     procedure :: send_data => fms_send_data  !!TODO
! Get ID functions
     procedure :: get_id => fms_diag_get_id
     procedure :: id_from_name => diag_field_id_from_name
     procedure :: copy => copy_diag_obj
     procedure :: register => fms_register_diag_field_obj !! Merely initialize fields.
     procedure :: setID => set_diag_id
     procedure :: set_type => set_vartype
     procedure :: set_data_buffer => set_data_buffer
     procedure :: prepare_data_buffer
     procedure :: init_data_buffer
     procedure :: set_data_buffer_is_allocated
     procedure :: set_send_data_time
     procedure :: get_send_data_time
     procedure :: is_data_buffer_allocated
     procedure :: allocate_data_buffer
     procedure :: set_math_needs_to_be_done => set_math_needs_to_be_done
     procedure :: add_attribute => diag_field_add_attribute
     procedure :: vartype_inq => what_is_vartype
     procedure :: set_var_is_masked
     procedure :: get_var_is_masked
! Check functions
     procedure :: is_static => diag_obj_is_static
     procedure :: is_scalar
     procedure :: is_registered => get_registered
     procedure :: is_registeredB => diag_obj_is_registered
     procedure :: is_mask_variant => get_mask_variant
     procedure :: is_local => get_local
! Is variable allocated check functions
!TODO     procedure :: has_diag_field
     procedure :: has_diag_id
     procedure :: has_attributes
     procedure :: has_static
     procedure :: has_registered
     procedure :: has_mask_variant
     procedure :: has_local
     procedure :: has_vartype
     procedure :: has_varname
     procedure :: has_longname
     procedure :: has_standname
     procedure :: has_units
     procedure :: has_modname
     procedure :: has_realm
     procedure :: has_interp_method
     procedure :: has_frequency
     procedure :: has_tile_count
     procedure :: has_axis_ids
     procedure :: has_area
     procedure :: has_volume
     procedure :: has_missing_value
     procedure :: has_data_RANGE
     procedure :: has_input_data_buffer
! Get functions
     procedure :: get_attributes
     procedure :: get_static
     procedure :: get_registered
     procedure :: get_mask_variant
     procedure :: get_local
     procedure :: get_vartype
     procedure :: get_varname
     procedure :: get_longname
     procedure :: get_standname
     procedure :: get_units
     procedure :: get_modname
     procedure :: get_realm
     procedure :: get_interp_method
     procedure :: get_frequency
     procedure :: get_tile_count
     procedure :: get_area
     procedure :: get_volume
     procedure :: get_missing_value
     procedure :: get_data_RANGE
     procedure :: get_axis_id
     procedure :: get_data_buffer
     procedure :: get_mask
     procedure :: get_weight
     procedure :: dump_field_obj
     procedure :: get_domain
     procedure :: get_type_of_domain
     procedure :: set_file_ids
     procedure :: get_dimnames
     procedure :: get_var_skind
     procedure :: get_longname_to_write
     procedure :: get_multiple_send_data
     procedure :: write_field_metadata
     procedure :: write_coordinate_attribute
     procedure :: get_math_needs_to_be_done
     procedure :: add_area_volume
     procedure :: append_time_cell_methods
     procedure :: get_file_ids
     procedure :: set_mask
     procedure :: allocate_mask
     procedure :: set_halo_present
     procedure :: is_halo_present
     procedure :: find_missing_value
     procedure :: has_mask_allocated
     procedure :: is_variable_in_file
     procedure :: get_field_file_name
     procedure :: generate_associated_files_att
end type fmsDiagField_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type(fmsDiagField_type) :: null_ob

logical,private :: module_is_initialized = .false. !< Flag indicating if the module is initialized

!type(fmsDiagField_type) :: diag_object_placeholder (10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
public :: fmsDiagField_type
public :: fms_diag_fields_object_init
public :: null_ob
public :: fms_diag_field_object_end
public :: get_default_missing_value
public :: check_for_slices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Deallocates the array of diag_objs
subroutine fms_diag_field_object_end (ob)
  class (fmsDiagField_type), allocatable, intent(inout) :: ob(:) !< diag field object
  if (allocated(ob)) deallocate(ob)
  module_is_initialized = .false.
end subroutine fms_diag_field_object_end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \Description Allocates the diad field object array.
!! Sets the diag_id to the not registered value.
!! Initializes the number of registered variables to be 0
logical function fms_diag_fields_object_init(ob)
  class (fmsDiagField_type), allocatable, intent(inout) :: ob(:) !< diag field object
  integer :: i !< For looping
  allocate(ob(get_num_unique_fields()))
  do i = 1,size(ob)
      ob(i)%diag_id = diag_not_registered !null_ob%diag_id
      ob(i)%registered = .false.
  enddo
  module_is_initialized = .true.
  fms_diag_fields_object_init = .true.
end function fms_diag_fields_object_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \Description Fills in and allocates (when necessary) the values in the diagnostic object
subroutine fms_register_diag_field_obj &
       (this, modname, varname, diag_field_indices, diag_axis, axes, &
       longname, units, missing_value, varRange, mask_variant, standname, &
       do_not_log, err_msg, interp_method, tile_count, area, volume, realm, static, &
       multiple_send_data)

 class(fmsDiagField_type),       INTENT(inout) :: this                  !< Diaj_obj to fill
 CHARACTER(len=*),               INTENT(in)    :: modname               !< The module name
 CHARACTER(len=*),               INTENT(in)    :: varname               !< The variable name
 integer,                        INTENT(in)    :: diag_field_indices(:) !< Array of indices to the field
                                                                        !! in the yaml object
 class(fmsDiagAxisContainer_type),intent(in)   :: diag_axis(:)          !< Array of diag_axis
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
 LOGICAL,          OPTIONAL,     INTENT(in)    :: static                !< Set to true if it is a static field
 LOGICAL,          OPTIONAL,     INTENT(in)    :: multiple_send_data    !< .True. if send data is called, multiple
                                                                        !! times for the same time

!> Fill in information from the register call
  this%varname = trim(varname)
  this%modname = trim(modname)

!> Add the yaml info to the diag_object
  this%diag_field = get_diag_fields_entries(diag_field_indices)

!> Add axis and domain information
  if (present(axes)) then
    this%scalar = .false.
    this%axis_ids = axes
    call get_domain_and_domain_type(diag_axis, this%axis_ids, this%type_of_domain, this%domain, this%varname)
  else
    !> The variable is a scalar
    this%scalar = .true.
    this%type_of_domain = NO_DOMAIN
    this%domain => null()
  endif

!> get the optional arguments if included and the diagnostic is in the diag table
  if (present(longname))      this%longname      = trim(longname)
  if (present(standname))     this%standname     = trim(standname)

  !> Ignore the units if they are set to "none". This is to reproduce previous diag_manager behavior
  if (present(units)) then
    if (trim(units) .ne. "none") this%units = trim(units)
  endif
  if (present(realm))         this%realm         = trim(realm)
  if (present(interp_method)) this%interp_method = trim(interp_method)

  if (present(tile_count)) then
    allocate(this%tile_count)
    this%tile_count = tile_count
  endif
  if (present(static)) then
    this%static = static
  else
    this%static = .false.
  endif

  if (present(missing_value)) then
    select type (missing_value)
     type is (integer(kind=i4_kind))
             allocate(integer(kind=i4_kind) :: this%missing_value)
             this%missing_value = missing_value
     type is (integer(kind=i8_kind))
             allocate(integer(kind=i8_kind) :: this%missing_value)
             this%missing_value = missing_value
     type is (real(kind=r4_kind))
             allocate(real(kind=r4_kind) :: this%missing_value)
             this%missing_value = missing_value
     type is (real(kind=r8_kind))
             allocate(real(kind=r8_kind) :: this%missing_value)
             this%missing_value = missing_value
     class default
             call mpp_error("fms_register_diag_field_obj", &
                     "The missing value passed to register a diagnostic is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
  endif

  if (present(varRANGE)) then
    select type (varRANGE)
     type is (integer(kind=i4_kind))
             allocate(integer(kind=i4_kind) :: this%data_RANGE(2))
             this%data_RANGE = varRANGE
     type is (integer(kind=i8_kind))
             allocate(integer(kind=i8_kind) :: this%data_RANGE(2))
             this%data_RANGE = varRANGE
     type is (real(kind=r4_kind))
             allocate(integer(kind=r4_kind) :: this%data_RANGE(2))
             this%data_RANGE = varRANGE
     type is (real(kind=r8_kind))
             allocate(integer(kind=r8_kind) :: this%data_RANGE(2))
             this%data_RANGE = varRANGE
     class default
             call mpp_error("fms_register_diag_field_obj", &
                     "The varRange passed to register a diagnostic is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
  endif

  if (present(area)) then
    if (area < 0) call mpp_error("fms_register_diag_field_obj", &
                     "The area id passed with field_name"//trim(varname)//" has not been registered."&
                     "Check that there is a register_diag_field call for the AREA measure and that is in the"&
                     "diag_table.yaml", FATAL)
    allocate(this%area)
    this%area = area
  endif

  if (present(volume)) then
    if (volume < 0) call mpp_error("fms_register_diag_field_obj", &
                     "The volume id passed with field_name"//trim(varname)//" has not been registered."&
                     "Check that there is a register_diag_field call for the VOLUME measure and that is in the"&
                     "diag_table.yaml", FATAL)
    allocate(this%volume)
    this%volume = volume
  endif

  this%mask_variant = .false.
  if (present(mask_variant)) then
    this%mask_variant = mask_variant
  endif

  if (present(do_not_log)) then
    allocate(this%do_not_log)
    this%do_not_log = do_not_log
  endif

  if (present(multiple_send_data)) then
    this%multiple_send_data = multiple_send_data
  else
    this%multiple_send_data = .false.
  endif

 !< Allocate space for any additional variable attributes
 !< These will be fill out when calling `diag_field_add_attribute`
 allocate(this%attributes(max_field_attributes))
 this%num_attributes = 0
 this%registered = .true.
end subroutine fms_register_diag_field_obj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief Sets the diag_id.  This can only be done if a variable is unregistered
subroutine set_diag_id(this , id)
 class (fmsDiagField_type) , intent(inout):: this
 integer                                :: id
 if (allocated(this%registered)) then
     if (this%registered) then
          call mpp_error("set_diag_id", "The variable"//this%varname//" is already registered", FATAL)
     else
       this%diag_id = id
     endif
 else
     this%diag_id = id
 endif
end subroutine set_diag_id

!> \brief Find the type of the variable and store it in the object
subroutine set_vartype(objin , var)
 class (fmsDiagField_type) , intent(inout):: objin
 class(*)                               :: var
 select type (var)
     type is (real(kind=8))
          objin%vartype = r8
     type is (real(kind=4))
          objin%vartype = r4
     type is (integer(kind=8))
          objin%vartype = i8
     type is (integer(kind=4))
          objin%vartype = i4
     type is (character(*))
          objin%vartype = string
     class default
          objin%vartype = null_type_int
          call mpp_error("set_vartype", "The variable"//objin%varname//" is not a supported type "// &
          " r8, r4, i8, i4, or string.", warning)
 end select
end subroutine set_vartype

!> @brief Sets the time send data was called last
subroutine set_send_data_time (this, time)
  class (fmsDiagField_type) , intent(inout):: this                !< The field object
  type(time_type),            intent(in)   :: time                !< Current model time

  call this%input_data_buffer%set_send_data_time(time)
end subroutine set_send_data_time

!> @brief Get the time send data was called last
!! @result the time send data was called last
function get_send_data_time(this) &
  result(rslt)
  class (fmsDiagField_type) , intent(in):: this                  !< The field object
  type(time_type) :: rslt

  rslt = this%input_data_buffer%get_send_data_time()
end function get_send_data_time

!> @brief Prepare the input_data_buffer to do the reduction method
subroutine prepare_data_buffer(this)
  class (fmsDiagField_type) , intent(inout):: this                !< The field object

  if (.not. this%multiple_send_data) return
  if (this%mask_variant) return
  call this%input_data_buffer%prepare_input_buffer_object(this%modname//":"//this%varname)
end subroutine prepare_data_buffer

!> @brief Initialize the input_data_buffer
subroutine init_data_buffer(this)
  class (fmsDiagField_type) , intent(inout):: this                !< The field object

  if (.not. this%multiple_send_data) return
  if (this%mask_variant) return
  call this%input_data_buffer%init_input_buffer_object()
end subroutine init_data_buffer

!> @brief Adds the input data to the buffered data.
subroutine set_data_buffer (this, input_data, mask, weight, is, js, ks, ie, je, ke)
  class (fmsDiagField_type) , intent(inout):: this                !< The field object
  class(*),                   intent(in)   :: input_data(:,:,:,:) !< The input array
  logical,                    intent(in)   :: mask(:,:,:,:)       !< Mask that is passed into
                                                                  !! send_data
  real(kind=r8_kind),         intent(in)   :: weight              !< The field weight
  integer,                    intent(in)   :: is, js, ks          !< Starting indicies of the field_data relative
                                                                  !! to the compute domain (1 based)
  integer,                    intent(in)   :: ie, je, ke          !< Ending indicies of the field_data relative
                                                                  !! to the compute domain (1 based)

  character(len=128) :: err_msg !< Error msg
  if (.not.this%data_buffer_is_allocated) &
    call mpp_error ("set_data_buffer", "The data buffer for the field "//trim(this%varname)//" was unable to be "//&
      "allocated.", FATAL)
  if (this%multiple_send_data) then
    err_msg = this%input_data_buffer%update_input_buffer_object(input_data, is, js, ks, ie, je, ke, &
                                                                mask, this%mask, this%mask_variant, this%var_is_masked)
  else
    this%mask(is:ie, js:je, ks:ke, :) = mask
    err_msg = this%input_data_buffer%set_input_buffer_object(input_data, weight, is, js, ks, ie, je, ke)
  endif
  if (trim(err_msg) .ne. "") call mpp_error(FATAL, "Field:"//trim(this%varname)//" -"//trim(err_msg))

end subroutine set_data_buffer
!> Allocates the global data buffer for a given field using a single thread. Returns true when the
!! buffer is allocated
logical function allocate_data_buffer(this, input_data, diag_axis)
  class (fmsDiagField_type), target, intent(inout):: this !< The field object
  class(*), dimension(:,:,:,:), intent(in) :: input_data !< The input array
  class(fmsDiagAxisContainer_type),intent(in)   :: diag_axis(:) !< Array of diag_axis

  character(len=128) :: err_msg !< Error msg
  err_msg = ""

  allocate(this%input_data_buffer)
  err_msg = this%input_data_buffer%allocate_input_buffer_object(input_data, this%axis_ids, diag_axis)
  if (trim(err_msg) .ne. "") then
    call mpp_error(FATAL, "Field:"//trim(this%varname)//" -"//trim(err_msg))
    return
  endif

  allocate_data_buffer = .true.
end function allocate_data_buffer
!> Sets the flag saying that the math functions need to be done
subroutine set_math_needs_to_be_done (this, math_needs_to_be_done)
  class (fmsDiagField_type) , intent(inout):: this
  logical, intent (in) :: math_needs_to_be_done !< Flag saying that the math functions need to be done
  this%math_needs_to_be_done = math_needs_to_be_done
end subroutine set_math_needs_to_be_done

!> @brief Set the mask_variant to .true.
subroutine set_var_is_masked(this, is_masked)
  class (fmsDiagField_type) , intent(inout):: this      !< The diag field object
  logical,                    intent (in)  :: is_masked !< .True. if the field is masked

  this%var_is_masked = is_masked
end subroutine set_var_is_masked

!> @brief Queries a field for the var_is_masked variable
!! @return var_is_masked
function get_var_is_masked(this) &
  result(rslt)
  class (fmsDiagField_type) , intent(inout):: this      !< The diag field object
  logical :: rslt !< .True. if the field is masked

  rslt = this%var_is_masked
end function get_var_is_masked

!> @brief Sets the flag saying that the data buffer is allocated
subroutine set_data_buffer_is_allocated (this, data_buffer_is_allocated)
  class (fmsDiagField_type) , intent(inout) :: this                     !< The field object
  logical,                    intent (in)   :: data_buffer_is_allocated !< .true. if the
                                                                        !! data buffer is allocated
  this%data_buffer_is_allocated = data_buffer_is_allocated
end subroutine set_data_buffer_is_allocated

!> @brief Determine if the data_buffer is allocated
!! @return logical indicating if the data_buffer is allocated
pure logical function is_data_buffer_allocated (this)
  class (fmsDiagField_type) , intent(in) :: this                     !< The field object

  is_data_buffer_allocated = .false.
  if (allocated(this%data_buffer_is_allocated)) is_data_buffer_allocated = this%data_buffer_is_allocated

end function
!> \brief Prints to the screen what type the diag variable is
subroutine what_is_vartype(this)
 class (fmsDiagField_type) , intent(inout):: this
 if (.not. allocated(this%vartype)) then
     call mpp_error("what_is_vartype", "The variable type has not been set prior to this call", warning)
     return
 endif
 select case (this%vartype)
     case (r8)
          call mpp_error("what_is_vartype", "The variable type of "//trim(this%varname)//&
          " is REAL(kind=8)", NOTE)
     case (r4)
          call mpp_error("what_is_vartype", "The variable type of "//trim(this%varname)//&
          " is REAL(kind=4)", NOTE)
     case (i8)
          call mpp_error("what_is_vartype", "The variable type of "//trim(this%varname)//&
          " is INTEGER(kind=8)", NOTE)
     case (i4)
          call mpp_error("what_is_vartype", "The variable type of "//trim(this%varname)//&
          " is INTEGER(kind=4)", NOTE)
     case (string)
          call mpp_error("what_is_vartype", "The variable type of "//trim(this%varname)//&
          " is CHARACTER(*)", NOTE)
     case (null_type_int)
          call mpp_error("what_is_vartype", "The variable type of "//trim(this%varname)//&
          " was not set", WARNING)
     case default
          call mpp_error("what_is_vartype", "The variable type of "//trim(this%varname)//&
          " is not supported by diag_manager", FATAL)
 end select
end subroutine what_is_vartype
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief Copies the calling object into the object that is the argument of the subroutine
subroutine copy_diag_obj(this , objout)
 class (fmsDiagField_type)      , intent(in)                :: this
 class (fmsDiagField_type)      , intent(inout) , allocatable :: objout !< The destination of the copy
select type (objout)
 class is (fmsDiagField_type)

  if (allocated(this%registered)) then
     objout%registered = this%registered
  else
     call mpp_error("copy_diag_obj", "You can only copy objects that have been registered",warning)
  endif
     objout%diag_id = this%diag_id

     if (allocated(this%attributes)) objout%attributes = this%attributes
     objout%static = this%static
     if (allocated(this%frequency)) objout%frequency = this%frequency
     if (allocated(this%varname)) objout%varname = this%varname
end select
end subroutine copy_diag_obj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief Returns the ID integer for a variable
!! \return the diag ID
pure integer function fms_diag_get_id (this) result(diag_id)
 class(fmsDiagField_type)     , intent(in)            :: this
!> Check if the diag_object registration has been done
 if (allocated(this%registered)) then
         !> Return the diag_id if the variable has been registered
         diag_id = this%diag_id
 else
!> If the variable is not regitered, then return the unregistered value
        diag_id = DIAG_NOT_REGISTERED
 endif
end function fms_diag_get_id

!> Function to return a character (string) representation of the most basic
!> object identity info. Intended for debugging and warning. The format produced is:
!> [this: o.varname(string|?), vartype (string|?), o.registered (T|F|?), diag_id (id|?)].
!> A questionmark "?" is set in place of the variable that is not yet allocated
!>TODO: Add diag_id ?
function fms_diag_obj_as_string_basic(this) result(rslt)
    class(fmsDiagField_type), allocatable, intent(in) :: this
    character(:), allocatable :: rslt
    character (len=:), allocatable :: registered, vartype, varname, diag_id
    if ( .not. allocated (this)) then
        varname = "?"
        vartype = "?"
        registered = "?"
        diag_id = "?"
        rslt = "[Obj:" // varname // "," // vartype // "," // registered // "," // diag_id // "]"
        return
     end if

!   if(allocated (this%registered)) then
!       registered = logical_to_cs (this%registered)
!   else
!       registered = "?"
!   end if

!   if(allocated (this%diag_id)) then
!     diag_id = int_to_cs (this%diag_id)
!   else
!       diag_id = "?"
!   end if

!   if(allocated (this%vartype)) then
!       vartype = int_to_cs (this%vartype)
!   else
!       registered = "?"
!   end if

    if(allocated (this%varname)) then
        varname = this%varname
    else
        registered = "?"
    end if

    rslt = "[Obj:" // varname // "," // vartype // "," // registered // "," // diag_id // "]"

end function fms_diag_obj_as_string_basic


function diag_obj_is_registered (this) result (rslt)
    class(fmsDiagField_type), intent(in) :: this
    logical :: rslt
    rslt = this%registered
end function diag_obj_is_registered

function diag_obj_is_static (this) result (rslt)
    class(fmsDiagField_type), intent(in) :: this
    logical :: rslt
    rslt = .false.
    if (allocated(this%static)) rslt = this%static
end function diag_obj_is_static

!> @brief Determine if the field is a scalar
!! @return .True. if the field is a scalar
function is_scalar (this) result (rslt)
  class(fmsDiagField_type), intent(in) :: this !< diag_field object
  logical                              :: rslt
  rslt = this%scalar
end function is_scalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get functions

!> @brief Gets attributes
!! @return A pointer to the attributes of the diag_obj, null pointer if there are no attributes
function get_attributes (this) &
result(rslt)
  class (fmsDiagField_type), target, intent(in) :: this !< diag object
  type(fmsDiagAttribute_type), pointer :: rslt(:)

  rslt => null()
  if (this%num_attributes > 0 ) rslt => this%attributes
end function get_attributes

!> @brief Gets static
!! @return copy of variable static
pure function get_static (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     logical :: rslt
     rslt = this%static
end function get_static

!> @brief Gets regisetered
!! @return copy of registered
pure function get_registered (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     logical :: rslt
     rslt = this%registered
end function get_registered

!> @brief Gets mask variant
!! @return copy of mask variant
pure function get_mask_variant (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     logical :: rslt
     rslt = .false.
     if (allocated(this%mask_variant)) rslt = this%mask_variant
end function get_mask_variant

!> @brief Gets local
!! @return copy of local
pure function get_local (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     logical :: rslt
     rslt = this%local
end function get_local

!> @brief Gets vartype
!! @return copy of The integer related to the variable type
pure function get_vartype (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     integer :: rslt
     rslt = this%vartype
end function get_vartype

!> @brief Gets varname
!! @return copy of the variable name
pure function get_varname (this, to_write) &
result(rslt)
  class (fmsDiagField_type), intent(in) :: this     !< diag object
  logical, optional,         intent(in) :: to_write !< .true. if getting the varname that will be writen to the file
  character(len=:), allocatable :: rslt
  rslt = this%varname

  !< If writing the varname can be the outname which is defined in the yaml
  if (present(to_write)) then
    if (to_write) then
    !TODO this is wrong
    rslt = this%diag_field(1)%get_var_outname()
    endif
  endif

end function get_varname

!> @brief Gets longname
!! @return copy of the variable long name or a single string if there is no long name
pure function get_longname (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     character(len=:), allocatable :: rslt
     if (allocated(this%longname)) then
       rslt = this%longname
     else
       rslt = diag_null_string
     endif
end function get_longname

!> @brief Gets standname
!! @return copy of the standard name or an empty string if standname is not allocated
pure function get_standname (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     character(len=:), allocatable :: rslt
     if (allocated(this%standname)) then
       rslt = this%standname
     else
       rslt = diag_null_string
     endif
end function get_standname

!> @brief Gets units
!! @return copy of the units or an empty string if not allocated
pure function get_units (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     character(len=:), allocatable :: rslt
     if (allocated(this%units)) then
       rslt = this%units
     else
       rslt = diag_null_string
     endif
end function get_units

!> @brief Gets modname
!! @return copy of the module name that the variable is in or an empty string if not allocated
pure function get_modname (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     character(len=:), allocatable :: rslt
     if (allocated(this%modname)) then
       rslt = this%modname
     else
       rslt = diag_null_string
     endif
end function get_modname

!> @brief Gets realm
!! @return copy of the variables modeling realm or an empty string if not allocated
pure function get_realm (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     character(len=:), allocatable :: rslt
     if (allocated(this%realm)) then
       rslt = this%realm
     else
       rslt = diag_null_string
     endif
end function get_realm

!> @brief Gets interp_method
!! @return copy of The interpolation method or an empty string if not allocated
pure function get_interp_method (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     character(len=:), allocatable :: rslt
     if (allocated(this%interp_method)) then
       rslt = this%interp_method
     else
       rslt = diag_null_string
     endif
end function get_interp_method

!> @brief Gets frequency
!! @return copy of the  frequency or DIAG_NULL if obj%frequency is not allocated
pure function get_frequency (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     integer, allocatable, dimension (:) :: rslt
     if (allocated(this%frequency)) then
       allocate (rslt(size(this%frequency)))
       rslt = this%frequency
     else
       allocate (rslt(1))
       rslt = DIAG_NULL
     endif
end function get_frequency

!> @brief Gets tile_count
!! @return copy of the number of tiles or diag_null if tile_count is not allocated
pure function get_tile_count (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     integer :: rslt
     if (allocated(this%tile_count)) then
       rslt = this%tile_count
     else
       rslt = DIAG_NULL
     endif
end function get_tile_count

!> @brief Gets area
!! @return copy of the area or diag_null if not allocated
pure function get_area (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     integer :: rslt
     if (allocated(this%area)) then
       rslt = this%area
     else
       rslt = diag_null
     endif
end function get_area

!> @brief Gets volume
!! @return copy of the volume or diag_null if volume is not allocated
pure function get_volume (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     integer :: rslt
     if (allocated(this%volume)) then
       rslt = this%volume
     else
       rslt = diag_null
     endif
end function get_volume

!> @brief Gets missing_value
!! @return copy of The missing value
!! @note Netcdf requires the type of the variable and the type of the missing_value and _Fillvalue to be the same
!! var_type is the type of the variable which may not be in the same type as the missing_value in the register call
!! For example, if compiling with r8 but the in diag_table.yaml the kind is r4
function get_missing_value (this, var_type) &
result(rslt)
  class (fmsDiagField_type), intent(in) :: this     !< diag object
  integer,                   intent(in) :: var_type !< The type of the variable as it will writen to the netcdf file
                                                    !! and the missing value is return as

  class(*),allocatable :: rslt

  if (.not. allocated(this%missing_value)) then
    call mpp_error ("get_missing_value", &
                 "The missing value is not allocated", FATAL)
  endif

  !< The select types are needed so that the missing_value can be correctly converted and copied as the needed variable
  !! type
  select case (var_type)
  case (r4)
    allocate (real(kind=r4_kind) :: rslt)
    select type (miss => this%missing_value)
    type is (real(kind=r4_kind))
      select type (rslt)
      type is (real(kind=r4_kind))
        rslt = real(miss, kind=r4_kind)
      end select
    type is (real(kind=r8_kind))
      select type (rslt)
      type is (real(kind=r4_kind))
        rslt = real(miss, kind=r4_kind)
      end select
    end select
  case (r8)
    allocate (real(kind=r8_kind) :: rslt)
    select type (miss => this%missing_value)
    type is (real(kind=r4_kind))
      select type (rslt)
      type is (real(kind=r8_kind))
        rslt = real(miss, kind=r8_kind)
      end select
    type is (real(kind=r8_kind))
      select type (rslt)
      type is (real(kind=r8_kind))
        rslt = real(miss, kind=r8_kind)
      end select
    end select
  end select

end function get_missing_value

!> @brief Gets data_range
!! @return copy of the data range
!! @note Netcdf requires the type of the variable and the type of the range to be the same
!! var_type is the type of the variable which may not be in the same type as the range in the register call
!! For example, if compiling with r8 but the in diag_table.yaml the kind is r4
function get_data_RANGE (this, var_type) &
result(rslt)
  class (fmsDiagField_type), intent(in) :: this      !< diag object
  integer,                   intent(in) :: var_type  !< The type of the variable as it will writen to the netcdf file
                                                     !! and the data_range is returned as
  class(*),allocatable :: rslt(:)

  if ( .not. allocated(this%data_RANGE)) call mpp_error ("get_data_RANGE", &
    "The data_RANGE value is not allocated", FATAL)

  !< The select types are needed so that the range can be correctly converted and copied as the needed variable
  !! type
  select case (var_type)
  case (r4)
    allocate (real(kind=r4_kind) :: rslt(2))
    select type (r => this%data_RANGE)
    type is (real(kind=r4_kind))
      select type (rslt)
      type is (real(kind=r4_kind))
        rslt = real(r, kind=r4_kind)
      end select
    type is (real(kind=r8_kind))
      select type (rslt)
      type is (real(kind=r4_kind))
        rslt = real(r, kind=r4_kind)
      end select
    end select
  case (r8)
    allocate (real(kind=r8_kind) :: rslt(2))
    select type (r => this%data_RANGE)
    type is (real(kind=r4_kind))
      select type (rslt)
      type is (real(kind=r8_kind))
        rslt = real(r, kind=r8_kind)
      end select
    type is (real(kind=r8_kind))
      select type (rslt)
      type is (real(kind=r8_kind))
        rslt = real(r, kind=r8_kind)
      end select
    end select
  end select
end function get_data_RANGE

!> @brief Gets axis_ids
!! @return pointer to the axis ids
function get_axis_id (this) &
result(rslt)
  class (fmsDiagField_type), target, intent(in) :: this   !< diag object
  integer, pointer, dimension(:)    :: rslt  !< field's axis_ids

  if(allocated(this%axis_ids)) then
    rslt => this%axis_ids
  else
    rslt => null()
  endif
end function get_axis_id

!> @brief Gets field's domain
!! @return pointer to the domain
function get_domain (this) &
result(rslt)
  class (fmsDiagField_type), target, intent(in) :: this  !< diag field
  class(diagDomain_t),       pointer            :: rslt  !< field's domain

  if (associated(this%domain)) then
    rslt => this%domain
  else
    rslt => null()
  endif

end function get_domain

!> @brief Gets field's type of domain
!! @return integer defining the type of domain (NO_DOMAIN, TWO_D_DOMAIN, UG_DOMAIN)
pure function get_type_of_domain (this) &
result(rslt)
  class (fmsDiagField_type), target, intent(in) :: this  !< diag field
  integer                                       :: rslt  !< field's domain

  rslt = this%type_of_domain
end function get_type_of_domain

!> @brief Set the file ids of the files that the field belongs to
subroutine set_file_ids(this, file_ids)
  class (fmsDiagField_type), intent(inout) :: this        !< diag field
  integer,                   intent(in)    :: file_ids(:) !< File_ids to add

  allocate(this%file_ids(size(file_ids)))
  this%file_ids = file_ids
end subroutine set_file_ids

!> @brief Get the kind of the variable based on the yaml
!! @return A string indicating the kind of the variable (as it is used in fms2_io)
pure function get_var_skind(this, field_yaml) &
result(rslt)
  class (fmsDiagField_type),   intent(in) :: this       !< diag field
  type(diagYamlFilesVar_type), intent(in) :: field_yaml !< The corresponding yaml of the field

  character(len=:), allocatable :: rslt

  integer :: var_kind !< The integer corresponding to the kind of the variable (i4, i8, r4, r8)

  var_kind = field_yaml%get_var_kind()
  select case (var_kind)
  case (r4)
    rslt = "float"
  case (r8)
    rslt = "double"
  case (i4)
    rslt = "int"
  case (i8)
    rslt = "int64"
  end select

end function get_var_skind

!> @brief Get the multiple_send_data member of the field object
!! @return multiple_send_data of the field
pure function get_multiple_send_data(this) &
result(rslt)
  class (fmsDiagField_type),   intent(in) :: this       !< diag field
  logical :: rslt
  rslt = this%multiple_send_data
end function get_multiple_send_data

!> @brief Determine the long name to write for the field
!! @return Long name to write
pure function get_longname_to_write(this, field_yaml) &
result(rslt)
  class (fmsDiagField_type),   intent(in) :: this       !< diag field
  type(diagYamlFilesVar_type), intent(in) :: field_yaml !< The corresponding yaml of the field

  character(len=:), allocatable :: rslt

  rslt = field_yaml%get_var_longname() !! This is the long name defined in the yaml
  if (rslt .eq. "") then !! If the long name is not defined in the yaml, use the long name in the
                         !! register_diag_field
    rslt = this%get_longname()
  else
    return
  endif
  if (rslt .eq. "") then !! If the long name is not defined in the yaml and in the register_diag_field
                         !! use the variable name
    rslt = field_yaml%get_var_varname()
  endif
end function get_longname_to_write

!> @brief Determine the dimension names to use when registering the field to fms2_io
subroutine get_dimnames(this, diag_axis, field_yaml, unlim_dimname, dimnames, is_regional)
  class (fmsDiagField_type),        target, intent(inout) :: this          !< diag field
  class(fmsDiagAxisContainer_type), target, intent(in)    :: diag_axis(:)  !< Diag_axis object
  type(diagYamlFilesVar_type),              intent(in)    :: field_yaml    !< Field info from diag_table yaml
  character(len=*),                         intent(in)    :: unlim_dimname !< The name of unlimited dimension
  character(len=120), allocatable,          intent(out)   :: dimnames(:)   !< Array of the dimension names
                                                                           !! for the field
  logical,                                  intent(in)    :: is_regional   !< Flag indicating if the field is regional

  integer                                   :: i     !< For do loops
  integer                                   :: naxis !< Number of axis for the field
  class(fmsDiagAxisContainer_type), pointer :: axis_ptr !diag_axis(this%axis_ids(i), for convenience
  character(len=23)                         :: diurnal_axis_name !< name of the diurnal axis

  if (this%is_static()) then
    naxis = size(this%axis_ids)
  else
    naxis = size(this%axis_ids) + 1 !< Adding 1 more dimension for the unlimited dimension
  endif

  if (field_yaml%has_n_diurnal()) then
    naxis = naxis + 1 !< Adding 1 more dimension for the diurnal axis
  endif

  allocate(dimnames(naxis))

  !< Duplicated do loops for performance
  if (field_yaml%has_var_zbounds()) then
    do i = 1, size(this%axis_ids)
      axis_ptr => diag_axis(this%axis_ids(i))
      if (axis_ptr%axis%is_z_axis()) then
        dimnames(i) = axis_ptr%axis%get_axis_name(is_regional)//"_sub01"
      else
        dimnames(i) = axis_ptr%axis%get_axis_name(is_regional)
      endif
    enddo
  else
    do i = 1, size(this%axis_ids)
      axis_ptr => diag_axis(this%axis_ids(i))
      dimnames(i) = axis_ptr%axis%get_axis_name(is_regional)
    enddo
  endif

  !< The second to last dimension is always the diurnal axis
  if (field_yaml%has_n_diurnal()) then
    WRITE (diurnal_axis_name,'(a,i2.2)') 'time_of_day_', field_yaml%get_n_diurnal()
    dimnames(naxis - 1) = trim(diurnal_axis_name)
  endif

  !< The last dimension is always the unlimited dimensions
  if (.not. this%is_static()) dimnames(naxis) = unlim_dimname

end subroutine get_dimnames

!> @brief Wrapper for the register_field call. The select types are needed so that the code can go
!! in the correct interface
subroutine register_field_wrap(fms2io_fileobj, varname, vartype, dimensions)
  class(FmsNetcdfFile_t),            INTENT(INOUT) :: fms2io_fileobj!< Fms2_io fileobj to write to
  character(len=*),                  INTENT(IN)    :: varname       !< Name of the variable
  character(len=*),                  INTENT(IN)    :: vartype       !< The type of the variable
  character(len=*), optional,        INTENT(IN)    :: dimensions(:) !< The dimension names of the field

  select type(fms2io_fileobj)
  type is (FmsNetcdfFile_t)
    call register_field(fms2io_fileobj, varname, vartype, dimensions)
  type is (FmsNetcdfDomainFile_t)
    call register_field(fms2io_fileobj, varname, vartype, dimensions)
  type is (FmsNetcdfUnstructuredDomainFile_t)
    call register_field(fms2io_fileobj, varname, vartype, dimensions)
  end select
end subroutine register_field_wrap

!> @brief Write the field's metadata to the file
subroutine write_field_metadata(this, fms2io_fileobj, file_id, yaml_id, diag_axis, unlim_dimname, is_regional, &
                                cell_measures)
  class (fmsDiagField_type), target, intent(inout) :: this          !< diag field
  class(FmsNetcdfFile_t),            INTENT(INOUT) :: fms2io_fileobj!< Fms2_io fileobj to write to
  integer,                           intent(in)    :: file_id       !< File id of the file to write to
  integer,                           intent(in)    :: yaml_id       !< Yaml id of the yaml entry of this field
  class(fmsDiagAxisContainer_type),  intent(in)    :: diag_axis(:)  !< Diag_axis object
  character(len=*),                  intent(in)    :: unlim_dimname !< The name of the unlimited dimension
  logical,                           intent(in)    :: is_regional   !< Flag indicating if the field is regional
  character(len=*),                  intent(in)    :: cell_measures !< The cell measures attribute to write

  type(diagYamlFilesVar_type), pointer     :: field_yaml  !< pointer to the yaml entry
  character(len=:),            allocatable :: var_name    !< Variable name
  character(len=:),            allocatable :: long_name   !< Longname to write
  character(len=:),            allocatable :: units       !< Units of the field to write
  character(len=120),          allocatable :: dimnames(:) !< Dimension names of the field
  character(len=120)                       :: cell_methods!< Cell methods attribute to write
  integer                                  :: i           !< For do loops
  character (len=MAX_STR_LEN), allocatable :: yaml_field_attributes(:,:) !< Variable attributes defined in the yaml

  field_yaml => diag_yaml%get_diag_field_from_id(yaml_id)
  var_name = field_yaml%get_var_outname()

  if (allocated(this%axis_ids)) then
    call this%get_dimnames(diag_axis, field_yaml, unlim_dimname, dimnames, is_regional)
    call register_field_wrap(fms2io_fileobj, var_name, this%get_var_skind(field_yaml), dimnames)
  else
    if (this%is_static()) then
      call register_field_wrap(fms2io_fileobj, var_name, this%get_var_skind(field_yaml))
    else
      !< In this case, the scalar variable is a function of time, so we need to pass in the
      !! unlimited dimension as a dimension
      call register_field_wrap(fms2io_fileobj, var_name, this%get_var_skind(field_yaml), (/unlim_dimname/))
    endif
  endif

  long_name = this%get_longname_to_write(field_yaml)
  call register_variable_attribute(fms2io_fileobj, var_name, "long_name", long_name, str_len=len_trim(long_name))

  units = this%get_units()
  if (units .ne. diag_null_string) &
    call register_variable_attribute(fms2io_fileobj, var_name, "units", units, str_len=len_trim(units))

  if (this%has_missing_value()) then
    call register_variable_attribute(fms2io_fileobj, var_name, "missing_value", &
      this%get_missing_value(field_yaml%get_var_kind()))
    call register_variable_attribute(fms2io_fileobj, var_name, "_FillValue", &
      this%get_missing_value(field_yaml%get_var_kind()))
  else
    call register_variable_attribute(fms2io_fileobj, var_name, "missing_value", &
      get_default_missing_value(field_yaml%get_var_kind()))
      call register_variable_attribute(fms2io_fileobj, var_name, "_FillValue", &
      get_default_missing_value(field_yaml%get_var_kind()))
  endif

  if (this%has_data_RANGE()) then
    call register_variable_attribute(fms2io_fileobj, var_name, "valid_range", &
      this%get_data_range(field_yaml%get_var_kind()))
  endif

  if (this%has_interp_method()) then
    call register_variable_attribute(fms2io_fileobj, var_name, "interp_method", this%get_interp_method(), &
      str_len=len_trim(this%get_interp_method()))
  endif

  cell_methods = ""
  !< Check if any of the attributes defined via a "diag_field_add_attribute" call
  !! are the cell_methods, if so add to the "cell_methods" variable:
  do i = 1, this%num_attributes
    call this%attributes(i)%write_metadata(fms2io_fileobj, var_name, &
      cell_methods=cell_methods)
  enddo

  !< Append the time cell methods based on the variable's reduction
  call this%append_time_cell_methods(cell_methods, field_yaml)
  if (trim(cell_methods) .ne. "") &
    call register_variable_attribute(fms2io_fileobj, var_name, "cell_methods", &
      trim(adjustl(cell_methods)), str_len=len_trim(adjustl(cell_methods)))

  !< Write out the cell_measures attribute (i.e Area, Volume)
  !! The diag field ids for the Area and Volume are sent in the register call
  !! This was defined in file object and passed in here
  if (trim(cell_measures) .ne. "") &
    call register_variable_attribute(fms2io_fileobj, var_name, "cell_measures", &
      trim(adjustl(cell_measures)), str_len=len_trim(adjustl(cell_measures)))

  !< Write out the standard_name (this was defined in the register call)
  if (this%has_standname()) &
  call register_variable_attribute(fms2io_fileobj, var_name, "standard_name", &
    trim(this%get_standname()), str_len=len_trim(this%get_standname()))

  call this%write_coordinate_attribute(fms2io_fileobj, var_name, diag_axis)

  if (field_yaml%has_var_attributes()) then
    yaml_field_attributes = field_yaml%get_var_attributes()
    do i = 1, size(yaml_field_attributes,1)
      call register_variable_attribute(fms2io_fileobj, var_name, trim(yaml_field_attributes(i,1)), &
      trim(yaml_field_attributes(i,2)), str_len=len_trim(yaml_field_attributes(i,2)))
    enddo
    deallocate(yaml_field_attributes)
  endif
end subroutine write_field_metadata

!> @brief Writes the coordinate attribute of a field if any of the field's axis has an
!! auxiliary axis
subroutine write_coordinate_attribute (this, fms2io_fileobj, var_name, diag_axis)
  CLASS(fmsDiagField_type),          intent(in)    :: this         !< The field object
  class(FmsNetcdfFile_t),            INTENT(INOUT) :: fms2io_fileobj!< Fms2_io fileobj to write to
  character(len=*),                  intent(in)    :: var_name     !< Variable name
  class(fmsDiagAxisContainer_type),  intent(in)    :: diag_axis(:) !< Diag_axis object

  integer              :: i         !< For do loops
  character(len = 252) :: aux_coord !< Auxuliary axis name

  !> If the variable is a scalar, go away
  if (.not. allocated(this%axis_ids)) return

  !> Determine if any of the field's axis has an auxiliary axis and the
  !! axis_names as a variable attribute
  aux_coord = ""
  do i = 1, size(this%axis_ids)
    select type (obj => diag_axis(this%axis_ids(i))%axis)
    type is (fmsDiagFullAxis_type)
      if (obj%has_aux()) then
        aux_coord = trim(aux_coord)//" "//obj%get_aux()
      endif
    end select
  enddo

  if (trim(aux_coord) .eq. "") return

  call register_variable_attribute(fms2io_fileobj, var_name, "coordinates", &
    trim(adjustl(aux_coord)), str_len=len_trim(adjustl(aux_coord)))

end subroutine write_coordinate_attribute

!> @brief Gets a fields data buffer
!! @return a pointer to the data buffer
function get_data_buffer (this) &
  result(rslt)
  class (fmsDiagField_type), target, intent(in) :: this  !< diag field
  class(*),dimension(:,:,:,:), pointer      :: rslt !< The field's data buffer

  if (.not. this%data_buffer_is_allocated) &
  call mpp_error(FATAL, "The input data buffer for the field:"&
    //trim(this%varname)//" was never allocated.")

  rslt => this%input_data_buffer%get_buffer()
end function get_data_buffer


!> @brief Gets a fields weight buffer
!! @return a pointer to the weight buffer
function get_weight (this) &
  result(rslt)
  class (fmsDiagField_type), target, intent(in) :: this  !< diag field
  type(real(kind=r8_kind)), pointer :: rslt

  if (.not. this%data_buffer_is_allocated) &
  call mpp_error(FATAL, "The input data buffer for the field:"&
    //trim(this%varname)//" was never allocated.")

  rslt => this%input_data_buffer%get_weight()
end function get_weight

!> Gets the flag telling if the math functions need to be done
!! \return Copy of math_needs_to_be_done flag
pure logical function get_math_needs_to_be_done(this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  get_math_needs_to_be_done = .false.
  if (allocated(this%math_needs_to_be_done)) get_math_needs_to_be_done = this%math_needs_to_be_done
end function get_math_needs_to_be_done
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Allocation checks

!!> @brief Checks if obj%diag_field is allocated
!!! @return true if obj%diag_field is allocated
!logical function has_diag_field (obj)
!  class (fmsDiagField_type), intent(in) :: obj !< diag object
!  has_diag_field = allocated(obj%diag_field)
!end function has_diag_field
!> @brief Checks if obj%diag_id is allocated
!! @return true if obj%diag_id is allocated
pure logical function has_diag_id (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_diag_id = allocated(this%diag_id)
end function has_diag_id

!> @brief Checks if obj%metadata is allocated
!! @return true if obj%metadata is allocated
pure logical function has_attributes (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_attributes = this%num_attributes > 0
end function has_attributes

!> @brief Checks if obj%static is allocated
!! @return true if obj%static is allocated
pure logical function has_static (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_static = allocated(this%static)
end function has_static

!> @brief Checks if obj%registered is allocated
!! @return true if obj%registered is allocated
pure logical function has_registered (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_registered = allocated(this%registered)
end function has_registered

!> @brief Checks if obj%mask_variant is allocated
!! @return true if obj%mask_variant is allocated
pure logical function has_mask_variant (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_mask_variant = allocated(this%mask_variant)
end function has_mask_variant

!> @brief Checks if obj%local is allocated
!! @return true if obj%local is allocated
pure logical function has_local (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_local = allocated(this%local)
end function has_local

!> @brief Checks if obj%vartype is allocated
!! @return true if obj%vartype is allocated
pure logical function has_vartype (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_vartype = allocated(this%vartype)
end function has_vartype

!> @brief Checks if obj%varname is allocated
!! @return true if obj%varname is allocated
pure logical function has_varname (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_varname = allocated(this%varname)
end function has_varname

!> @brief Checks if obj%longname is allocated
!! @return true if obj%longname is allocated
pure logical function has_longname (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_longname = allocated(this%longname)
end function has_longname

!> @brief Checks if obj%standname is allocated
!! @return true if obj%standname is allocated
pure logical function has_standname (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_standname = allocated(this%standname)
end function has_standname

!> @brief Checks if obj%units is allocated
!! @return true if obj%units is allocated
pure logical function has_units (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_units = allocated(this%units)
end function has_units

!> @brief Checks if obj%modname is allocated
!! @return true if obj%modname is allocated
pure logical function has_modname (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_modname = allocated(this%modname)
end function has_modname

!> @brief Checks if obj%realm is allocated
!! @return true if obj%realm is allocated
pure logical function has_realm (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_realm = allocated(this%realm)
end function has_realm

!> @brief Checks if obj%interp_method is allocated
!! @return true if obj%interp_method is allocated
pure logical function has_interp_method (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_interp_method = allocated(this%interp_method)
end function has_interp_method

!> @brief Checks if obj%frequency is allocated
!! @return true if obj%frequency is allocated
pure logical function has_frequency (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_frequency = allocated(this%frequency)
end function has_frequency

!> @brief Checks if obj%tile_count is allocated
!! @return true if obj%tile_count is allocated
pure logical function has_tile_count (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_tile_count = allocated(this%tile_count)
end function has_tile_count

!> @brief Checks if axis_ids of the object is allocated
!! @return true if it is allocated
pure logical function has_axis_ids (this)
  class (fmsDiagField_type), intent(in) :: this !< diag field object
  has_axis_ids = allocated(this%axis_ids)
end function has_axis_ids

!> @brief Checks if obj%area is allocated
!! @return true if obj%area is allocated
pure logical function has_area (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_area = allocated(this%area)
end function has_area

!> @brief Checks if obj%volume is allocated
!! @return true if obj%volume is allocated
pure logical function has_volume (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_volume = allocated(this%volume)
end function has_volume

!> @brief Checks if obj%missing_value is allocated
!! @return true if obj%missing_value is allocated
pure logical function has_missing_value (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_missing_value = allocated(this%missing_value)
end function has_missing_value

!> @brief Checks if obj%data_RANGE is allocated
!! @return true if obj%data_RANGE is allocated
pure logical function has_data_RANGE (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_data_RANGE = allocated(this%data_RANGE)
end function has_data_RANGE

!> @brief Checks if obj%input_data_buffer is allocated
!! @return true if obj%input_data_buffer is allocated
pure logical function has_input_data_buffer (this)
  class (fmsDiagField_type), intent(in) :: this !< diag object
  has_input_data_buffer = allocated(this%input_data_buffer)
end function has_input_data_buffer

!> @brief Add a attribute to the diag_obj using the diag_field_id
subroutine diag_field_add_attribute(this, att_name, att_value)
  class (fmsDiagField_type), intent (inout) :: this !< The field object
  character(len=*), intent(in) :: att_name     !< Name of the attribute
  class(*),         intent(in) :: att_value(:) !< The attribute value to add

  this%num_attributes = this%num_attributes + 1
  if (this%num_attributes > max_field_attributes) &
    call mpp_error(FATAL, "diag_field_add_attribute: Number of attributes exceeds max_field_attributes for field:"&
                           //trim(this%varname)//".  Increase diag_manager_nml:max_field_attributes.")

  call this%attributes(this%num_attributes)%add(att_name, att_value)
end subroutine diag_field_add_attribute

!> @brief Determine the default missing value to use based on the requested variable type
!! @return The missing value
function get_default_missing_value(var_type) &
  result(rslt)

  integer, intent(in) :: var_type !< The type of the variable to return the missing value as
  class(*),allocatable :: rslt

  select case(var_type)
  case (r4)
    allocate(real(kind=r4_kind) :: rslt)
    rslt = real(CMOR_MISSING_VALUE, kind=r4_kind)
  case (r8)
    allocate(real(kind=r8_kind) :: rslt)
    rslt = real(CMOR_MISSING_VALUE, kind=r8_kind)
  case default
  end select
end function

!> @brief Determines the diag_obj id corresponding to a module name and field_name
!> @return diag_obj id
PURE FUNCTION diag_field_id_from_name(this, module_name, field_name) &
  result(diag_field_id)
  CLASS(fmsDiagField_type), INTENT(in) :: this !< The field object
  CHARACTER(len=*), INTENT(in) :: module_name !< Module name that registered the variable
  CHARACTER(len=*), INTENT(in) :: field_name !< Variable name

  integer :: diag_field_id

  diag_field_id = DIAG_FIELD_NOT_FOUND
  if (this%get_varname() .eq. trim(field_name) .and. &
      this%get_modname() .eq. trim(module_name)) then
        diag_field_id = this%get_id()
  endif
end function diag_field_id_from_name

!> @brief Adds the area and volume id to a field object
subroutine add_area_volume(this, area, volume)
  CLASS(fmsDiagField_type),          intent(inout) :: this   !< The field object
  INTEGER,                 optional, INTENT(in)    :: area   !< diag ids of area
  INTEGER,                 optional, INTENT(in)    :: volume !< diag ids of volume

  if (present(area)) then
    if (area > 0) then
      this%area = area
    else
      call mpp_error(FATAL, "diag_field_add_cell_measures: the area id is not valid. "&
                           &"Verify that the area_id passed in to the field:"//this%varname//&
                           &" is valid and that the field is registered and in the diag_table.yaml")
    endif
  endif

  if (present(volume)) then
    if (volume > 0) then
      this%volume = volume
    else
      call mpp_error(FATAL, "diag_field_add_cell_measures: the volume id is not valid. "&
                           &"Verify that the volume_id passed in to the field:"//this%varname//&
                           &" is valid and that the field is registered and in the diag_table.yaml")
    endif
  endif

end subroutine add_area_volume

!> @brief Append the time cell meathods based on the variable's reduction
subroutine append_time_cell_methods(this, cell_methods, field_yaml)
  class (fmsDiagField_type),   target, intent(inout) :: this          !< diag field
  character(len=*),                    intent(inout) :: cell_methods  !< The cell methods var to append to
  type(diagYamlFilesVar_type),         intent(in)    :: field_yaml    !< The field's yaml

  if (this%static) then
    cell_methods = trim(cell_methods)//" time: point "
    return
  endif

  select case (field_yaml%get_var_reduction())
  case (time_none)
    cell_methods = trim(cell_methods)//" time: point "
  case (time_diurnal)
    cell_methods = trim(cell_methods)//" time: mean"
  case (time_power)
    cell_methods = trim(cell_methods)//" time: mean_pow"//int2str(field_yaml%get_pow_value())
  case (time_rms)
    cell_methods = trim(cell_methods)//" time: root_mean_square"
  case (time_max)
    cell_methods = trim(cell_methods)//" time: max"
  case (time_min)
    cell_methods = trim(cell_methods)//" time: min"
  case (time_average)
    cell_methods = trim(cell_methods)//" time: mean"
  case (time_sum)
    cell_methods = trim(cell_methods)//" time: sum"
  end select
end subroutine append_time_cell_methods

!> Dumps any data from a given fmsDiagField_type object
subroutine dump_field_obj (this, unit_num)
  class(fmsDiagField_type), intent(in) :: this
  integer, intent(in) :: unit_num !< passed in from dump_diag_obj if log file is being written to
  integer :: i

  if( mpp_pe() .eq. mpp_root_pe()) then
    if( allocated(this%file_ids)) write(unit_num, *) 'file_ids:' ,this%file_ids
    if( allocated(this%diag_id)) write(unit_num, *) 'diag_id:' ,this%diag_id
    if( allocated(this%static)) write(unit_num, *) 'static:' ,this%static
    if( allocated(this%registered)) write(unit_num, *) 'registered:' ,this%registered
    if( allocated(this%mask_variant)) write(unit_num, *) 'mask_variant:' ,this%mask_variant
    if( allocated(this%do_not_log)) write(unit_num, *) 'do_not_log:' ,this%do_not_log
    if( allocated(this%local)) write(unit_num, *) 'local:' ,this%local
    if( allocated(this%vartype)) write(unit_num, *) 'vartype:' ,this%vartype
    if( allocated(this%varname)) write(unit_num, *) 'varname:' ,this%varname
    if( allocated(this%longname)) write(unit_num, *) 'longname:' ,this%longname
    if( allocated(this%standname)) write(unit_num, *) 'standname:' ,this%standname
    if( allocated(this%units)) write(unit_num, *) 'units:' ,this%units
    if( allocated(this%modname)) write(unit_num, *) 'modname:' ,this%modname
    if( allocated(this%realm)) write(unit_num, *) 'realm:' ,this%realm
    if( allocated(this%interp_method)) write(unit_num, *) 'interp_method:' ,this%interp_method
    if( allocated(this%tile_count)) write(unit_num, *) 'tile_count:' ,this%tile_count
    if( allocated(this%axis_ids)) write(unit_num, *) 'axis_ids:' ,this%axis_ids
    write(unit_num, *) 'type_of_domain:' ,this%type_of_domain
    if( allocated(this%area)) write(unit_num, *) 'area:' ,this%area
    if( allocated(this%missing_value)) then
      select type(missing_val => this%missing_value)
        type is (real(r4_kind))
          write(unit_num, *) 'missing_value:', missing_val
        type is (real(r8_kind))
          write(unit_num, *) 'missing_value:' ,missing_val
       type is(integer(i4_kind))
          write(unit_num, *) 'missing_value:' ,missing_val
       type is(integer(i8_kind))
          write(unit_num, *) 'missing_value:' ,missing_val
      end select
    endif
    if( allocated( this%data_RANGE)) then
      select type(drange => this%data_RANGE)
        type is (real(r4_kind))
          write(unit_num, *) 'data_RANGE:' ,drange
        type is (real(r8_kind))
          write(unit_num, *) 'data_RANGE:' ,drange
       type is(integer(i4_kind))
          write(unit_num, *) 'data_RANGE:' ,drange
       type is(integer(i8_kind))
          write(unit_num, *) 'data_RANGE:' ,drange
      end select
    endif
    write(unit_num, *) 'num_attributes:' ,this%num_attributes
    if( allocated(this%attributes)) then
      do i=1, this%num_attributes
        if( allocated(this%attributes(i)%att_value)) then
          select type( val => this%attributes(i)%att_value)
            type is (real(r8_kind))
              write(unit_num, *) 'attribute name', this%attributes(i)%att_name, 'val:',  val
            type is (real(r4_kind))
              write(unit_num, *) 'attribute name', this%attributes(i)%att_name, 'val:',  val
            type is (integer(i4_kind))
              write(unit_num, *) 'attribute name', this%attributes(i)%att_name, 'val:',  val
            type is (integer(i8_kind))
              write(unit_num, *) 'attribute name', this%attributes(i)%att_name, 'val:', val
          end select
        endif
      enddo
    endif

  endif

end subroutine

!< @brief Get the starting compute domain indices for a set of axis
!! @return compute domain starting indices
function get_starting_compute_domain(axis_ids, diag_axis) &
result(compute_domain)
  integer,                         intent(in) :: axis_ids(:)  !< Array of axis ids
  class(fmsDiagAxisContainer_type),intent(in) :: diag_axis(:) !< Array of axis object

  integer :: compute_domain(4)
  integer :: a              !< For looping through axes
  integer :: compute_idx(2) !< Compute domain indices (starting, ending)
  logical :: dummy          !< Dummy variable for the `get_compute_domain` subroutine

  compute_domain = 1
  axis_loop: do a = 1,size(axis_ids)
    select type (axis => diag_axis(axis_ids(a))%axis)
      type is (fmsDiagFullAxis_type)
        call axis%get_compute_domain(compute_idx, dummy)
        if ( compute_idx(1) .ne. diag_null) compute_domain(a) = compute_idx(1)
    end select
  enddo axis_loop
end function get_starting_compute_domain

!> Get list of field ids
pure function get_file_ids(this)
  class(fmsDiagField_type), intent(in) :: this
  integer, allocatable :: get_file_ids(:) !< Ids of the FMS_diag_files the variable
  get_file_ids = this%file_ids
end function

!> @brief Get the mask from the input buffer object
!! @return a pointer to the mask
function get_mask(this)
  class(fmsDiagField_type), target, intent(in) :: this !< input buffer object
  logical, pointer :: get_mask(:,:,:,:)
  get_mask => this%mask
end function get_mask

!> @brief If in openmp region, omp_axis should be provided in order to allocate to the given axis lengths.
!! Otherwise mask will be allocated to the size of mask_in
subroutine allocate_mask(this, mask_in, omp_axis)
  class(fmsDiagField_type), target, intent(inout) :: this !< input buffer object
  logical, intent(in) :: mask_in(:,:,:,:)
  class(fmsDiagAxisContainer_type), intent(in), optional :: omp_axis(:) !< true if calling from omp region
  integer :: axis_num, length(4)
  integer, pointer :: id_num
  ! if not omp just allocate to whatever is given
  if(.not. present(omp_axis)) then
    allocate(this%mask(size(mask_in,1), size(mask_in,2), size(mask_in,3), &
                     size(mask_in,4)))
  ! otherwise loop through axis and get sizes
  else
    length = 1
    do axis_num=1, size(this%axis_ids)
      id_num => this%axis_ids(axis_num)
      select type(axis => omp_axis(id_num)%axis)
        type is (fmsDiagFullAxis_type)
          length(axis_num) = axis%axis_length()
      end select
    enddo
    allocate(this%mask(length(1), length(2), length(3), length(4)))
  endif
end subroutine allocate_mask

!> Sets previously allocated mask to mask_in at given index ranges
subroutine set_mask(this, mask_in, field_info, is, js, ks, ie, je, ke)
  class(fmsDiagField_type), intent(inout) :: this
  logical, intent(in)                     :: mask_in(:,:,:,:)
  character(len=*), intent(in)            :: field_info !< Field info to add to error message
  integer, optional, intent(in)           :: is, js, ks, ie, je, ke
  if(present(is)) then
    if(is .lt. lbound(this%mask,1) .or. ie .gt. ubound(this%mask,1) .or. &
      js .lt. lbound(this%mask,2) .or. je .gt. ubound(this%mask,2) .or. &
      ks .lt. lbound(this%mask,3) .or. ke .gt. ubound(this%mask,3)) then
        print *, "PE:", int2str(mpp_pe()), "The size of the mask is", &
          SHAPE(this%mask), &
          "But the indices passed in are is=", int2str(is), " ie=", int2str(ie),&
          " js=", int2str(js), " je=", int2str(je), &
          " ks=", int2str(ks), " ke=", int2str(ke), &
          " ", trim(field_info)
        call mpp_error(FATAL,"set_mask:: given indices out of bounds for allocated mask")
    endif
    this%mask(is:ie, js:je, ks:ke, :) = mask_in
  else
    this%mask = mask_in
  endif
end subroutine set_mask

!> sets halo_present to true
subroutine set_halo_present(this)
  class(fmsDiagField_type), intent(inout) :: this !< field object to modify
  this%halo_present = .true.
end subroutine set_halo_present

!> Getter for halo_present
pure function is_halo_present(this)
  class(fmsDiagField_type), intent(in) :: this !< field object to get from
  logical :: is_halo_present
  is_halo_present = this%halo_present
end function is_halo_present

!> Helper routine to find and set the netcdf missing value for a field
!! Always returns r8 due to reduction routine args
!! casts up to r8 from given missing val or default if needed
function find_missing_value(this, missing_val) &
  result(res)
  class(fmsDiagField_type), intent(in) :: this !< field object to get missing value for
  class(*), allocatable, intent(out) :: missing_val !< outputted netcdf missing value (oriignal type)
  real(r8_kind) :: res !< returned r8 copy of missing_val

  if(this%has_missing_value()) then
    missing_val = this%get_missing_value(this%get_vartype())
  else
    missing_val = get_default_missing_value(this%get_vartype())
  endif

  select type(missing_val)
    type is (real(r8_kind))
      res = missing_val
    type is (real(r4_kind))
      res = real(missing_val, r8_kind)
  end select
end function find_missing_value

!> @returns allocation status of logical mask array
!! this just indicates whether the mask array itself has been alloc'd
!! this is different from @ref has_mask_variant, which is set earlier for whether a mask is being used at all
pure logical function has_mask_allocated(this)
  class(fmsDiagField_type),intent(in) :: this !< field object to check mask allocation for
  has_mask_allocated = allocated(this%mask)
end function has_mask_allocated

!> @brief Determine if the variable is in the file
!! @return .True. if the varibale is in the file
pure function is_variable_in_file(this, file_id) &
result(res)
  class(fmsDiagField_type), intent(in) :: this    !< field object to check
  integer,                  intent(in) :: file_id !< File id to check
  logical :: res

  integer :: i

  res = .false.
  if (any(this%file_ids .eq. file_id)) res = .true.
end function is_variable_in_file

!> @brief Determine the name of the first file the variable is in
!! @return filename
function get_field_file_name(this) &
  result(res)
  class(fmsDiagField_type), intent(in) :: this    !< Field object to query
  character(len=:), allocatable :: res

  res = this%diag_field(1)%get_var_fname()
end function get_field_file_name

!> @brief Generate the associated files attribute
subroutine generate_associated_files_att(this, att, start_time)
  class(fmsDiagField_type)        ,  intent(in)            :: this       !< diag_field_object for the area/volume field
  character(len=*),                  intent(inout)         :: att        !< associated_files_att
  type(time_type),                   intent(in)            :: start_time !< The start_time for the field's file

  character(len=:), allocatable :: field_name !< Name of the area/volume field
  character(len=MAX_STR_LEN) :: file_name !< Name of the file the area/volume field is in!
  character(len=128) :: start_date !< Start date to append to the begining of the filename

  integer :: year, month, day, hour, minute, second
  field_name = this%get_varname(to_write = .true.)

  ! Check if the field is already in the associated files attribute (i.e the area can be associated with multiple
  ! fields in the file, but it only needs to be added once)
  if (index(att, field_name) .ne. 0) return

  file_name = this%get_field_file_name()

  if (prepend_date) then
    call get_date(start_time, year, month, day, hour, minute, second)
    write (start_date, '(1I20.4, 2I2.2)') year, month, day
    file_name = TRIM(adjustl(start_date))//'.'//TRIM(file_name)
  endif

  att = trim(att)//" "//trim(field_name)//": "//trim(file_name)//".nc"
end subroutine generate_associated_files_att

!> @brief Determines if the compute domain has been divide further into slices (i.e openmp blocks)
!! @return .True. if the compute domain has been divided furter into slices
function check_for_slices(field, diag_axis, var_size) &
  result(rslt)
  type(fmsDiagField_type),                 intent(in) :: field        !< Field object
  type(fmsDiagAxisContainer_type), target, intent(in) :: diag_axis(:) !< Array of diag axis
  integer,                                 intent(in) :: var_size(:)  !< The size of the buffer pass into send_data

  logical :: rslt
  integer :: i !< For do loops

  rslt = .false.

  if (.not. field%has_axis_ids()) then
    rslt = .false.
    return
  endif
  do i = 1, size(field%axis_ids)
    select type (axis_obj => diag_axis(field%axis_ids(i))%axis)
    type is (fmsDiagFullAxis_type)
      if (axis_obj%axis_length() .ne. var_size(i)) then
        rslt = .true.
        return
      endif
    end select
  enddo
end function
#endif
end module fms_diag_field_object_mod
