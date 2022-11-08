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
use diag_data_mod,  only: diag_null, CMOR_MISSING_VALUE, diag_null_string
use diag_data_mod,  only: r8, r4, i8, i4, string, null_type_int, NO_DOMAIN
use diag_data_mod,  only: max_field_attributes, fmsDiagAttribute_type
use diag_data_mod,  only: diag_null, diag_not_found, diag_not_registered, diag_registered_id, &
                         &DIAG_FIELD_NOT_FOUND
use mpp_mod, only: fatal, note, warning, mpp_error, mpp_pe, mpp_root_pe
use fms_diag_yaml_mod, only:  diagYamlFilesVar_type, get_diag_fields_entries, get_diag_files_id, &
  & find_diag_field, get_num_unique_fields
use fms_diag_axis_object_mod, only: diagDomain_t, get_domain_and_domain_type, fmsDiagAxis_type, &
  & fmsDiagAxisContainer_type
use time_manager_mod, ONLY: time_type
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
     logical, allocatable, private                    :: registered        !< true when registered
     logical, allocatable, private                    :: mask_variant      !< If there is a mask variant
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
    contains
!     procedure :: send_data => fms_send_data  !!TODO
! Get ID functions
     procedure :: get_id => fms_diag_get_id
     procedure :: id_from_name => diag_field_id_from_name
     procedure :: copy => copy_diag_obj
     procedure :: register => fms_register_diag_field_obj !! Merely initialize fields.
     procedure :: setID => set_diag_id
     procedure :: set_type => set_vartype
     procedure :: add_attribute => diag_field_add_attribute
     procedure :: vartype_inq => what_is_vartype
! Check functions
     procedure :: is_static => diag_obj_is_static
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
     procedure :: has_area
     procedure :: has_volume
     procedure :: has_missing_value
     procedure :: has_data_RANGE
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
     procedure :: dump_field_obj
     procedure :: get_domain
     procedure :: get_type_of_domain
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
       do_not_log, err_msg, interp_method, tile_count, area, volume, realm, static)

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

!> Fill in information from the register call
  this%varname = trim(varname)
  this%modname = trim(modname)

!> Add the yaml info to the diag_object
  this%diag_field = get_diag_fields_entries(diag_field_indices)

!> Add axis and domain information
  if (present(axes)) then
    this%axis_ids = axes
    call get_domain_and_domain_type(diag_axis, this%axis_ids, this%type_of_domain, this%domain, this%varname)
  else
     !> The variable is a scalar
    this%type_of_domain = NO_DOMAIN
    this%domain => null()
  endif

!> get the optional arguments if included and the diagnostic is in the diag table
  if (present(longname))      this%longname      = trim(longname)
  if (present(standname))     this%standname     = trim(standname)
  if (present(units))         this%units         = trim(units)
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
             allocate(integer(kind=r4_kind) :: this%missing_value)
             this%missing_value = missing_value
     type is (real(kind=r8_kind))
             allocate(integer(kind=r8_kind) :: this%missing_value)
             this%missing_value = missing_value
     class default
             call mpp_error("fms_register_diag_field_obj", &
                     "The missing value passed to register a diagnostic is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
  else
      allocate(real :: this%missing_value)
      select type (miss => this%missing_value)
       type is (real)
        miss = real(CMOR_MISSING_VALUE)
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
  else
      allocate(real :: this%data_RANGE(2))
      select type (varRANGE => this%data_RANGE)
       type is (real)
        varRANGE = real(CMOR_MISSING_VALUE)
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

  if (present(mask_variant)) then
    allocate(this%mask_variant)
    this%mask_variant = mask_variant
  endif

  if (present(do_not_log)) then
    allocate(this%do_not_log)
    this%do_not_log = do_not_log
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
    rslt = this%static
end function diag_obj_is_static

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
     rslt = this%mask_variant
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
pure function get_varname (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     character(len=:), allocatable :: rslt
     rslt = this%varname
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
function get_missing_value (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     class(*),allocatable :: rslt
     if (allocated(this%missing_value)) then
       select type (miss => this%missing_value)
         type is (integer(kind=i4_kind))
             allocate (integer(kind=i4_kind) :: rslt)
             rslt = miss
         type is (integer(kind=i8_kind))
             allocate (integer(kind=i8_kind) :: rslt)
             rslt = miss
         type is (real(kind=r4_kind))
             allocate (integer(kind=i4_kind) :: rslt)
             rslt = miss
         type is (real(kind=r8_kind))
             allocate (integer(kind=i4_kind) :: rslt)
             rslt = miss
         class default
             call mpp_error ("get_missing_value", &
                     "The missing value is not a r8, r4, i8, or i4",&
                     FATAL)
         end select
       else
         call mpp_error ("get_missing_value", &
                 "The missing value is not allocated", FATAL)
       endif
end function get_missing_value

!> @brief Gets data_range
!! @return copy of the data range
function get_data_RANGE (this) &
result(rslt)
     class (fmsDiagField_type), intent(in) :: this !< diag object
     class(*),allocatable :: rslt(:)
     if (allocated(this%data_RANGE)) then
       select type (r => this%data_RANGE)
         type is (integer(kind=i4_kind))
             allocate (integer(kind=i4_kind) :: rslt(2))
             rslt = r
         type is (integer(kind=i8_kind))
             allocate (integer(kind=i8_kind) :: rslt(2))
             rslt = r
         type is (real(kind=r4_kind))
             allocate (integer(kind=i4_kind) :: rslt(2))
             rslt = r
         type is (real(kind=r8_kind))
             allocate (integer(kind=i4_kind) :: rslt(2))
             rslt = r
         class default
             call mpp_error ("get_data_RANGE", &
                     "The data_RANGE value is not a r8, r4, i8, or i4",&
                     FATAL)
         end select
       else
         call mpp_error ("get_data_RANGE", &
                 "The data_RANGE value is not allocated", FATAL)
       endif
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

#endif
end module fms_diag_field_object_mod
