module fms_diag_object_mod
!> \author Tom Robinson
!> \email thomas.robinson@noaa.gov
!! \brief Contains routines for the diag_objects
!!
!! \description The diag_manager passes an object back and forth between the diag routines and the users.
!! The procedures of this object and the types are all in this module.  The fms_dag_object is a type
!! that contains all of the information of the variable.  It is extended by a type that holds the
!! appropriate buffer for the data for manipulation.
use diag_data_mod,  only: diag_null, CMOR_MISSING_VALUE, diag_null_string
use diag_data_mod,  only: r8, r4, i8, i4, string, null_type_int, NO_DOMAIN
use diag_data_mod,  only: diag_null, diag_not_found, diag_not_registered, diag_registered_id

use diag_axis_mod,  only: diag_axis_type
use mpp_mod, only: fatal, note, warning, mpp_error
#ifdef use_yaml
use fms_diag_yaml_mod, only:  diagYamlFilesVar_type, get_diag_fields_entries, get_diag_files_id, &
  & find_diag_field, get_num_unique_fields
use fms_diag_file_object_mod, only: fmsDiagFile_type, FMS_diag_files
#endif
use fms_diag_axis_object_mod, only: diagDomain_t, get_domain_and_domain_type
use time_manager_mod, ONLY: time_type
!!!set_time, set_date, get_time, time_type, OPERATOR(>=), OPERATOR(>),&
!!!       & OPERATOR(<), OPERATOR(==), OPERATOR(/=), OPERATOR(/), OPERATOR(+), ASSIGNMENT(=), get_date, &
!!!       & get_ticks_per_second

use platform_mod
use iso_c_binding

implicit none

!> \brief Object that holds all variable information
type fmsDiagObject_type
#ifdef use_yaml
     type (diagYamlFilesVar_type), allocatable, dimension(:) :: diag_field !< info from diag_table for this variable
     integer,                      allocatable, dimension(:) :: file_ids   !< Ids of the FMS_diag_files the variable
                                                                           !! belongs to
#endif
     integer, allocatable, private                    :: diag_id           !< unique id for varable
     character(len=:), allocatable, dimension(:)      :: metadata          !< metadata for the variable
     logical, allocatable, private                    :: static            !< true if this is a static var
     logical, allocatable, private                    :: registered        !< true when registered
     logical, allocatable, private                    :: mask_variant      !< If there is a mask variant
     logical, allocatable, private                    :: do_not_log        !< .true. if no need to log the diag_field
     logical, allocatable, private                    :: local             !< If the output is local
     TYPE(time_type), private                         :: init_time         !< The initial time
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
     integer, allocatable, dimension(:), private      :: output_units
     integer, allocatable, private                    :: t
     integer, allocatable, private                    :: tile_count        !< The number of tiles
     integer, pointer, dimension(:), private          :: axis_ids          !< variable axis IDs
     class(diagDomain_t), pointer,   private          :: domain            !< Domain
     INTEGER                         , private        :: type_of_domain    !< The type of domain ("NO_DOMAIN",
                                                                           !! "TWO_D_DOMAIN", or "UG_DOMAIN")
     integer, allocatable, private                    :: area, volume      !< The Area and Volume
     class(*), allocatable, private                   :: missing_value     !< The missing fill value
     class(*), allocatable, private                   :: data_RANGE(:)     !< The range of the variable data
     class(*), allocatable :: vardata0                                     !< Scalar data buffer 
     class(*), allocatable, dimension(:) :: vardata1                       !< 1D data buffer
     class(*), allocatable, dimension(:,:) :: vardata2                     !< 2D data buffer
     class(*), allocatable, dimension(:,:,:) :: vardata3                   !< 3D data buffer
     class(*), allocatable, dimension(:,:,:,:) :: vardata4                 !< 4D data buffer
     class(*), allocatable, dimension(:,:,:,:,:) :: vardata5               !< 5D data buffer
    contains
!     procedure :: send_data => fms_send_data  !!TODO
     procedure :: init_ob => diag_obj_init
     procedure :: get_id => fms_diag_get_id
     procedure :: id => fms_diag_get_id
     procedure :: copy => copy_diag_obj
     procedure :: register => fms_register_diag_field_obj !! Merely initialize fields.
     procedure :: setID => set_diag_id
     procedure :: set_type => set_vartype
     procedure :: vartype_inq => what_is_vartype
! Check functions
     procedure :: is_static => diag_obj_is_static
     procedure :: is_registered => diag_ob_registered
     procedure :: is_registeredB => diag_obj_is_registered
     procedure :: is_mask_variant => get_mask_variant
     procedure :: is_local => get_local
! Is variable allocated check functions
!TODO     procedure :: has_diag_field
     procedure :: has_diag_id
     procedure :: has_metadata
     procedure :: has_static
     procedure :: has_registered
     procedure :: has_mask_variant
     procedure :: has_local
!TODO     procedure :: has_init_time
     procedure :: has_vartype
     procedure :: has_varname
     procedure :: has_longname
     procedure :: has_standname
     procedure :: has_units
     procedure :: has_modname
     procedure :: has_realm
     procedure :: has_interp_method
     procedure :: has_frequency
     procedure :: has_output_units
     procedure :: has_t
     procedure :: has_tile_count
     procedure :: has_area
     procedure :: has_volume
     procedure :: has_missing_value
     procedure :: has_data_RANGE
! Get functions
     procedure :: get_diag_id => fms_diag_get_id
     procedure :: get_metadata
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
     procedure :: get_output_units
     procedure :: get_t
     procedure :: get_tile_count
     procedure :: get_area
     procedure :: get_volume
     procedure :: get_missing_value
     procedure :: get_data_RANGE
!TODO     procedure :: get_init_time
!TODO     procedure :: get_axis
end type fmsDiagObject_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type(fmsDiagObject_type) :: null_ob

integer,private :: MAX_LEN_VARNAME
integer,private :: MAX_LEN_META
logical,private :: module_is_initialized = .false. !< Flag indicating if the module is initialized

TYPE(fmsDiagObject_type), private, ALLOCATABLE, target :: diag_objs(:) !< Array of diag objects
                                                                       !! one for each registered variable
integer, private :: registered_variables !< Number of registered variables

!type(fmsDiagObject_type) :: diag_object_placeholder (10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
public :: fmsDiagObject_type
public :: null_ob
public :: copy_diag_obj, fms_diag_get_id
public :: fms_diag_object_init
public :: fms_diag_object_end
public :: fms_register_diag_field_array
public :: fms_register_diag_field_scalar
public :: fms_register_static_field
public :: get_diag_obj_from_id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> @brief Initiliazes the array of diag_objs based on the number of unique diag_fields in the diag_table
subroutine fms_diag_object_init (mlv,mlm)
 integer, intent(in) :: mlv !< The maximum length of the varname
 integer, intent(in) :: mlm !< The maximum length of the metadata

 if (module_is_initialized) return

!> Get info from the namelist
 MAX_LEN_VARNAME = mlv
 MAX_LEN_META = mlm
!> Initialize the null_d variables
 null_ob%diag_id = DIAG_NULL
#ifdef use_yaml
 allocate(diag_objs(get_num_unique_fields()))
 registered_variables = 0
#endif
 module_is_initialized = .true.
end subroutine fms_diag_object_init

!> @brief Deallocates the array of diag_objs
subroutine fms_diag_object_end ()
  if (.not. module_is_initialized) return

  if (allocated(diag_objs)) deallocate(diag_objs)

  module_is_initialized = .false.
end subroutine fms_diag_object_end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \Description Sets the diag_id to the not registered value.
subroutine diag_obj_init(ob)
 class (fmsDiagObject_type)      , intent(inout)       :: ob
 select type (ob)
  class is (fmsDiagObject_type)
     ob%diag_id = diag_not_registered !null_ob%diag_id
     ob%registered = .false.
 end select
end subroutine diag_obj_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \Description Fills in and allocates (when necessary) the values in the diagnostic object
subroutine fms_register_diag_field_obj &
                !(dobj, modname, varname, axes, time, longname, units, missing_value, metadata)
       (dobj, modname, varname, diag_field_indices, axes, init_time, &
       longname, units, missing_value, varRange, mask_variant, standname, &
       do_not_log, err_msg, interp_method, tile_count, area, volume, realm, metadata)

 class(fmsDiagObject_type),      INTENT(inout) :: dobj                  !< Diaj_obj to fill
 CHARACTER(len=*),               INTENT(in)    :: modname               !< The module name
 CHARACTER(len=*),               INTENT(in)    :: varname               !< The variable name
 integer,                        INTENT(in)    :: diag_field_indices(:) !< Array of indices to the field
                                                                        !! in the yaml object
 TYPE(time_type),  OPTIONAL,     INTENT(in)    :: init_time             !< Initial time !< TO DO
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
 character(len=*), optional,     INTENT(in)    :: metadata(:)           !< metedata for the variable

 integer :: i !< For do loops
 integer :: j !< dobj%file_ids(i) (for less typing :)

#ifdef use_yaml
!> Fill in information from the register call
  dobj%varname = trim(varname)
  dobj%modname = trim(modname)

!> Fill in diag_field and find the ids of the files that this variable is in
  dobj%diag_field = get_diag_fields_entries(diag_field_indices)
  dobj%file_ids   = get_diag_files_id(diag_field_indices)

  if (present(axes)) then
    dobj%axis_ids => axes
    call get_domain_and_domain_type(dobj%axis_ids, dobj%type_of_domain, dobj%domain, dobj%varname)
    do i = 1, size(dobj%file_ids)
       j = dobj%file_ids(i)
       call FMS_diag_files(j)%set_file_domain(dobj%domain, dobj%type_of_domain)
       call FMS_diag_files(j)%add_axes(axes)
    enddo
     !> TO DO:
     !!     Mark the field as registered in the diag_files
  else
     !> The variable is a scalar
    dobj%type_of_domain = NO_DOMAIN
    dobj%domain => null()
  endif

!> get the optional arguments if included and the diagnostic is in the diag table
  if (present(longname))      dobj%longname      = trim(longname)
  if (present(standname))     dobj%standname     = trim(standname)
  if (present(units))         dobj%units         = trim(units)
  if (present(realm))         dobj%realm         = trim(realm)
  if (present(interp_method)) dobj%interp_method = trim(interp_method)
  if (present(tile_count)) then
    allocate(dobj%tile_count)
    dobj%tile_count = tile_count
  endif

  if (present(metadata)) then
     allocate(character(len=MAX_LEN_META) :: dobj%metadata(size(metadata)))
     dobj%metadata = metadata
  endif
  if (present(missing_value)) then
    select type (missing_value)
     type is (integer(kind=i4_kind))
             allocate(integer(kind=i4_kind) :: dobj%missing_value)
             dobj%missing_value = missing_value
     type is (integer(kind=i8_kind))
             allocate(integer(kind=i8_kind) :: dobj%missing_value)
             dobj%missing_value = missing_value
     type is (real(kind=r4_kind))
             allocate(integer(kind=r4_kind) :: dobj%missing_value)
             dobj%missing_value = missing_value
     type is (real(kind=r8_kind))
             allocate(integer(kind=r8_kind) :: dobj%missing_value)
             dobj%missing_value = missing_value
     class default
             call mpp_error("fms_register_diag_field_obj", &
                     "The missing value passed to register a diagnostic is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
  else
      allocate(real :: dobj%missing_value)
      select type (miss => dobj%missing_value)
       type is (real)
        miss = real(CMOR_MISSING_VALUE)
      end select
  endif

  if (present(varRANGE)) then
    select type (varRANGE)
     type is (integer(kind=i4_kind))
             allocate(integer(kind=i4_kind) :: dobj%data_RANGE(2))
             dobj%data_RANGE = varRANGE
     type is (integer(kind=i8_kind))
             allocate(integer(kind=i8_kind) :: dobj%data_RANGE(2))
             dobj%data_RANGE = varRANGE
     type is (real(kind=r4_kind))
             allocate(integer(kind=r4_kind) :: dobj%data_RANGE(2))
             dobj%data_RANGE = varRANGE
     type is (real(kind=r8_kind))
             allocate(integer(kind=r8_kind) :: dobj%data_RANGE(2))
             dobj%data_RANGE = varRANGE
     class default
             call mpp_error("fms_register_diag_field_obj", &
                     "The varRange passed to register a diagnostic is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
  else
      allocate(real :: dobj%data_RANGE(2))
      select type (varRANGE => dobj%data_RANGE)
       type is (real)
        varRANGE = real(CMOR_MISSING_VALUE)
      end select
  endif

  if (present(area)) then
    if (area < 0) call mpp_error("fms_register_diag_field_obj", &
                     "The area id passed with field_name"//trim(varname)//" has not been registered."&
                     "Check that there is a register_diag_field call for the AREA measure and that is in the"&
                     "diag_table.yaml", FATAL)
    allocate(dobj%area)
    dobj%area = area
  endif

  if (present(volume)) then
    if (volume < 0) call mpp_error("fms_register_diag_field_obj", &
                     "The volume id passed with field_name"//trim(varname)//" has not been registered."&
                     "Check that there is a register_diag_field call for the VOLUME measure and that is in the"&
                     "diag_table.yaml", FATAL)
    allocate(dobj%volume)
    dobj%volume = volume
  endif

  if (present(mask_variant)) then
    allocate(dobj%mask_variant)
    dobj%mask_variant = mask_variant
  endif

  if (present(do_not_log)) then
    allocate(dobj%do_not_log)
    dobj%do_not_log = do_not_log
  endif

 dobj%registered = .true.
#endif
end subroutine fms_register_diag_field_obj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \brief Sets the diag_id.  This can only be done if a variable is unregistered
subroutine set_diag_id(objin , id)
 class (fmsDiagObject_type) , intent(inout):: objin
 integer                                :: id
 if (allocated(objin%registered)) then
     if (objin%registered) then
          call mpp_error("set_diag_id", "The variable"//objin%varname//" is already registered", FATAL)
     endif
 else
     objin%diag_id = id
 endif
end subroutine set_diag_id
!> \brief Find the type of the variable and store it in the object
subroutine set_vartype(objin , var)
 class (fmsDiagObject_type) , intent(inout):: objin
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
subroutine what_is_vartype(objin)
 class (fmsDiagObject_type) , intent(inout):: objin
 if (.not. allocated(objin%vartype)) then
     call mpp_error("what_is_vartype", "The variable type has not been set prior to this call", warning)
     return
 endif
 select case (objin%vartype)
     case (r8)
          call mpp_error("what_is_vartype", "The variable type of "//trim(objin%varname)//&
          " is REAL(kind=8)", NOTE)
     case (r4)
          call mpp_error("what_is_vartype", "The variable type of "//trim(objin%varname)//&
          " is REAL(kind=4)", NOTE)
     case (i8)
          call mpp_error("what_is_vartype", "The variable type of "//trim(objin%varname)//&
          " is INTEGER(kind=8)", NOTE)
     case (i4)
          call mpp_error("what_is_vartype", "The variable type of "//trim(objin%varname)//&
          " is INTEGER(kind=4)", NOTE)
     case (string)
          call mpp_error("what_is_vartype", "The variable type of "//trim(objin%varname)//&
          " is CHARACTER(*)", NOTE)
     case (null_type_int)
          call mpp_error("what_is_vartype", "The variable type of "//trim(objin%varname)//&
          " was not set", WARNING)
     case default
          call mpp_error("what_is_vartype", "The variable type of "//trim(objin%varname)//&
          " is not supported by diag_manager", FATAL)
 end select
end subroutine what_is_vartype
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!MZ Is this a TODO. Many problems:
!> \brief Registers the object
subroutine diag_ob_registered(objin , reg)
 class (fmsDiagObject_type)      , intent(inout):: objin
 logical                      , intent(in)   :: reg !< If registering, this is true
 objin%registered = reg
end subroutine diag_ob_registered
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \brief Copies the calling object into the object that is the argument of the subroutine
subroutine copy_diag_obj(objin , objout)
 class (fmsDiagObject_type)      , intent(in)                :: objin
 class (fmsDiagObject_type)      , intent(inout) , allocatable :: objout !< The destination of the copy
select type (objout)
 class is (fmsDiagObject_type)

  if (allocated(objin%registered)) then
     objout%registered = objin%registered
  else
     call mpp_error("copy_diag_obj", "You can only copy objects that have been registered",warning)
  endif
     objout%diag_id = objin%diag_id

     if (allocated(objin%metadata)) objout%metadata = objin%metadata
     objout%static = objin%static
     if (allocated(objin%frequency)) objout%frequency = objin%frequency
     if (allocated(objin%varname)) objout%varname = objin%varname
end select
end subroutine copy_diag_obj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \brief Returns the ID integer for a variable
!! \return the diag ID
integer function fms_diag_get_id (dobj) result(diag_id)
 class(fmsDiagObject_type)     , intent(inout)            :: dobj
! character(*)               , intent(in)               :: varname
!> Check if the diag_object registration has been done
 if (allocated(dobj%registered)) then
         !> Return the diag_id if the variable has been registered
         diag_id = dobj%diag_id
 else
!> If the variable is not regitered, then return the unregistered value
        diag_id = DIAG_NOT_REGISTERED
 endif
end function fms_diag_get_id

!> Function to return a character (string) representation of the most basic
!> object identity info. Intended for debugging and warning. The format produced is:
!> [dobj: o.varname(string|?), vartype (string|?), o.registered (T|F|?), diag_id (id|?)].
!> A questionmark "?" is set in place of the variable that is not yet allocated
!>TODO: Add diag_id ?
function fms_diag_obj_as_string_basic(dobj) result(rslt)
    class(fmsDiagObject_type), allocatable, intent(in) :: dobj
    character(:), allocatable :: rslt
    character (len=:), allocatable :: registered, vartype, varname, diag_id
    if ( .not. allocated (dobj)) then
        varname = "?"
        vartype = "?"
        registered = "?"
        diag_id = "?"
        rslt = "[Obj:" // varname // "," // vartype // "," // registered // "," // diag_id // "]"
        return
     end if

!   if(allocated (dobj%registered)) then
!       registered = logical_to_cs (dobj%registered)
!   else
!       registered = "?"
!   end if

!   if(allocated (dobj%diag_id)) then
!     diag_id = int_to_cs (dobj%diag_id)
!   else
!       diag_id = "?"
!   end if

!   if(allocated (dobj%vartype)) then
!       vartype = int_to_cs (dobj%vartype)
!   else
!       registered = "?"
!   end if

    if(allocated (dobj%varname)) then
        varname = dobj%varname
    else
        registered = "?"
    end if

    rslt = "[Obj:" // varname // "," // vartype // "," // registered // "," // diag_id // "]"

end function fms_diag_obj_as_string_basic


function diag_obj_is_registered (obj) result (rslt)
    class(fmsDiagObject_type), intent(in) :: obj
    logical :: rslt
    rslt = obj%registered
end function diag_obj_is_registered

function diag_obj_is_static (obj) result (rslt)
    class(fmsDiagObject_type), intent(in) :: obj
    logical :: rslt
    rslt = obj%static
end function diag_obj_is_static

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get functions

!> @brief Gets metedata
!! @return copy of metadata string array, or a single space if metadata is not allocated
pure function get_metadata (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable, dimension(:) :: rslt
     if (allocated(obj%metadata)) then
       allocate(character(len=(len(obj%metadata))) :: rslt (size(obj%metadata)) )
       rslt = obj%metadata
     else
       allocate(character(len=1) :: rslt(1:1))
       rslt = diag_null_string
     endif
end function get_metadata
!> @brief Gets static
!! @return copy of variable static
pure function get_static (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     logical :: rslt 
     rslt = obj%static
end function get_static
!> @brief Gets regisetered
!! @return copy of registered
pure function get_registered (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     logical :: rslt 
     rslt = obj%registered
end function get_registered
!> @brief Gets mask variant
!! @return copy of mask variant
pure function get_mask_variant (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     logical :: rslt 
     rslt = obj%mask_variant
end function get_mask_variant
!> @brief Gets local
!! @return copy of local
pure function get_local (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     logical :: rslt 
     rslt = obj%local
end function get_local
!> @brief Gets initial time 
!! @return copy of the initial time
!! TODO
!function get_init_time (obj) &
!result(rslt)
!     class (fmsDiagObject_type), intent(in) :: obj !< diag object
!     TYPE(time_type) :: rslt 
!
!end function get_init_time
!> @brief Gets vartype 
!! @return copy of The integer related to the variable type
pure function get_vartype (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer :: rslt 
     rslt = obj%vartype
end function get_vartype
!> @brief Gets varname
!! @return copy of the variable name
pure function get_varname (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     rslt = obj%varname
end function get_varname
!> @brief Gets longname
!! @return copy of the variable long name or a single string if there is no long name
pure function get_longname (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     if (allocated(obj%longname)) then
       rslt = obj%longname
     else
       rslt = diag_null_string
     endif
end function get_longname
!> @brief Gets standname
!! @return copy of the standard name or an empty string if standname is not allocated
pure function get_standname (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     if (allocated(obj%standname)) then
       rslt = obj%standname
     else
       rslt = diag_null_string
     endif
end function get_standname
!> @brief Gets units
!! @return copy of the units or an empty string if not allocated
pure function get_units (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     if (allocated(obj%units)) then
       rslt = obj%units
     else
       rslt = diag_null_string
     endif
end function get_units
!> @brief Gets modname
!! @return copy of the module name that the variable is in or an empty string if not allocated
pure function get_modname (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     if (allocated(obj%modname)) then
       rslt = obj%modname
     else
       rslt = diag_null_string
     endif
end function get_modname
!> @brief Gets realm
!! @return copy of the variables modeling realm or an empty string if not allocated
pure function get_realm (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     if (allocated(obj%realm)) then
       rslt = obj%realm
     else
       rslt = diag_null_string
     endif
end function get_realm
!> @brief Gets interp_method
!! @return copy of The interpolation method or an empty string if not allocated
pure function get_interp_method (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     if (allocated(obj%interp_method)) then
       rslt = obj%interp_method
     else
       rslt = diag_null_string
     endif
end function get_interp_method
!> @brief Gets frequency
!! @return copy of the  frequency or DIAG_NULL if obj%frequency is not allocated
pure function get_frequency (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer, allocatable, dimension (:) :: rslt 
     if (allocated(obj%frequency)) then
       allocate (rslt(size(obj%frequency)))
       rslt = obj%frequency
     else
       allocate (rslt(1))
       rslt = DIAG_NULL
     endif
end function get_frequency
!> @brief Gets output_units
!! @return copy of The units of the output or DIAG_NULL is output_units is not allocated
pure function get_output_units (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer,allocatable, dimension (:) :: rslt 
     if (allocated(obj%output_units)) then
       allocate (rslt(size(obj%output_units)))
       rslt = obj%output_units
     else
       allocate (rslt(1))
       rslt = DIAG_NULL
     endif
end function get_output_units
!> @brief Gets t
!! @return copy of t 
pure function get_t (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer :: rslt 
     if (allocated(obj%t)) then
       rslt = obj%t
     else
       rslt = -999
     endif
end function get_t
!> @brief Gets tile_count
!! @return copy of the number of tiles or diag_null if tile_count is not allocated
pure function get_tile_count (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer :: rslt 
     if (allocated(obj%tile_count)) then
       rslt = obj%tile_count
     else
       rslt = DIAG_NULL
     endif
end function get_tile_count
!> @brief Gets area
!! @return copy of the area or diag_null if not allocated
pure function get_area (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer :: rslt 
     if (allocated(obj%area)) then
       rslt = obj%area
     else
       rslt = diag_null
     endif
end function get_area
!> @brief Gets volume
!! @return copy of the volume or diag_null if volume is not allocated
pure function get_volume (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer :: rslt
     if (allocated(obj%volume)) then
       rslt = obj%volume
     else
       rslt = diag_null
     endif
end function get_volume
!> @brief Gets missing_value
!! @return copy of The missing value
function get_missing_value (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     class(*),allocatable :: rslt
     if (allocated(obj%missing_value)) then
       select type (miss => obj%missing_value)
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
function get_data_RANGE (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     class(*),allocatable :: rslt(:) 
     if (allocated(obj%data_RANGE)) then
       select type (r => obj%data_RANGE)
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
!> @brief Gets axis
!! @return copy of axis information
!! TODO
!function get_axis (obj) &
!result(rslt)
!     class (fmsDiagObject_type), intent(in) :: obj !< diag object
!     type (diag_axis_type), allocatable, dimension(:) :: rslt 
!
!end function get_axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Allocation checks
!!> @brief Checks if obj%diag_field is allocated
!!! @return true if obj%diag_field is allocated
!logical function has_diag_field (obj)
!  class (fmsDiagObject_type), intent(in) :: obj !< diag object
!  has_diag_field = allocated(obj%diag_field)
!end function has_diag_field
!> @brief Checks if obj%diag_id is allocated
!! @return true if obj%diag_id is allocated
pure logical function has_diag_id (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_diag_id = allocated(obj%diag_id)
end function has_diag_id
!> @brief Checks if obj%metadata is allocated
!! @return true if obj%metadata is allocated
pure logical function has_metadata (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_metadata = allocated(obj%metadata)
end function has_metadata
!> @brief Checks if obj%static is allocated
!! @return true if obj%static is allocated
pure logical function has_static (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_static = allocated(obj%static)
end function has_static
!> @brief Checks if obj%registered is allocated
!! @return true if obj%registered is allocated
pure logical function has_registered (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_registered = allocated(obj%registered)
end function has_registered
!> @brief Checks if obj%mask_variant is allocated
!! @return true if obj%mask_variant is allocated
pure logical function has_mask_variant (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_mask_variant = allocated(obj%mask_variant)
end function has_mask_variant
!> @brief Checks if obj%local is allocated
!! @return true if obj%local is allocated
pure logical function has_local (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_local = allocated(obj%local)
end function has_local
!!> @brief Checks if obj%init_time is allocated
!!! @return true if obj%init_time is allocated
!logical function has_init_time (obj)
!  class (fmsDiagObject_type), intent(in) :: obj !< diag object
!  has_init_time = allocated(obj%init_time)
!end function has_init_time
!> @brief Checks if obj%vartype is allocated
!! @return true if obj%vartype is allocated
pure logical function has_vartype (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_vartype = allocated(obj%vartype)
end function has_vartype
!> @brief Checks if obj%varname is allocated
!! @return true if obj%varname is allocated
pure logical function has_varname (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_varname = allocated(obj%varname)
end function has_varname
!> @brief Checks if obj%longname is allocated
!! @return true if obj%longname is allocated
pure logical function has_longname (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_longname = allocated(obj%longname)
end function has_longname
!> @brief Checks if obj%standname is allocated
!! @return true if obj%standname is allocated
pure logical function has_standname (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_standname = allocated(obj%standname)
end function has_standname
!> @brief Checks if obj%units is allocated
!! @return true if obj%units is allocated
pure logical function has_units (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_units = allocated(obj%units)
end function has_units
!> @brief Checks if obj%modname is allocated
!! @return true if obj%modname is allocated
pure logical function has_modname (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_modname = allocated(obj%modname)
end function has_modname
!> @brief Checks if obj%realm is allocated
!! @return true if obj%realm is allocated
pure logical function has_realm (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_realm = allocated(obj%realm)
end function has_realm
!> @brief Checks if obj%interp_method is allocated
!! @return true if obj%interp_method is allocated
pure logical function has_interp_method (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_interp_method = allocated(obj%interp_method)
end function has_interp_method
!> @brief Checks if obj%frequency is allocated
!! @return true if obj%frequency is allocated
pure logical function has_frequency (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_frequency = allocated(obj%frequency)
end function has_frequency
!> @brief Checks if obj%output_units is allocated
!! @return true if obj%output_units is allocated
pure logical function has_output_units (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_output_units = allocated(obj%output_units)
end function has_output_units
!> @brief Checks if obj%t is allocated
!! @return true if obj%t is allocated
pure logical function has_t (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_t = allocated(obj%t)
end function has_t
!> @brief Checks if obj%tile_count is allocated
!! @return true if obj%tile_count is allocated
pure logical function has_tile_count (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_tile_count = allocated(obj%tile_count)
end function has_tile_count
!> @brief Checks if obj%area is allocated
!! @return true if obj%area is allocated
pure logical function has_area (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_area = allocated(obj%area)
end function has_area
!> @brief Checks if obj%volume is allocated
!! @return true if obj%volume is allocated
pure logical function has_volume (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_volume = allocated(obj%volume)
end function has_volume
!> @brief Checks if obj%missing_value is allocated
!! @return true if obj%missing_value is allocated
pure logical function has_missing_value (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_missing_value = allocated(obj%missing_value)
end function has_missing_value
!> @brief Checks if obj%data_RANGE is allocated
!! @return true if obj%data_RANGE is allocated
pure logical function has_data_RANGE (obj)
  class (fmsDiagObject_type), intent(in) :: obj !< diag object
  has_data_RANGE = allocated(obj%data_RANGE)
end function has_data_RANGE

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
    class(*),                       OPTIONAL, INTENT(in) :: missing_value !< Missing value to be added as a attribute
    class(*),         DIMENSION(:), OPTIONAL, INTENT(in) :: range         !< Range to be added as a attribute
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

!> @brief Get a pointer to the diag_object from the id.
!> @return A pointer to the diag_object or a null pointer if the id is not valid
FUNCTION get_diag_obj_from_id ( id ) result (obj_ptr)
  integer :: id !< Id of the diag_obj to get
  class(fmsDiagObject_type), pointer :: obj_ptr

  obj_ptr => null()
  IF (id >= 1 .and. id <= registered_variables) THEN
    obj_ptr => diag_objs(id)
  END IF
END FUNCTION get_diag_obj_from_id
end module fms_diag_object_mod
