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
use diag_data_mod,  only: r8, r4, i8, i4, string, null_type_int
use diag_data_mod,  only: diag_null, diag_not_found, diag_not_registered, diag_registered_id

use diag_axis_mod,  only: diag_axis_type
use mpp_mod, only: fatal, note, warning, mpp_error
#ifdef use_yaml
use fms_diag_yaml_mod, only: diagYamlFiles_type, diagYamlFilesVar_type
#endif
use time_manager_mod, ONLY: time_type
!!!set_time, set_date, get_time, time_type, OPERATOR(>=), OPERATOR(>),&
!!!       & OPERATOR(<), OPERATOR(==), OPERATOR(/=), OPERATOR(/), OPERATOR(+), ASSIGNMENT(=), get_date, &
!!!       & get_ticks_per_second

!use diag_util_mod,  only: int_to_cs, logical_to_cs
!USE diag_data_mod, ONLY: fileobjU, fileobj, fnum_for_domain, fileobjND

use fms2_io_mod
use platform_mod
use iso_c_binding

implicit none

!> \brief Object that holds all variable information
type fmsDiagObject_type
#ifdef use_yaml
     type (diagYamlFilesVar_type), allocatable, dimension(:) :: diag_field !< info from diag_table
     type (diagYamlFiles_type),     allocatable, dimension(:) :: diag_file  !< info from diag_table
#endif
     integer, allocatable, private                    :: diag_id           !< unique id for varable
     class(FmsNetcdfFile_t), dimension (:), pointer   :: fileob => NULL()  !< A pointer to all of the
                                                                           !! file objects for this variable
     character(len=:), allocatable, dimension(:)      :: metadata          !< metadata for the variable
     logical, allocatable, private                    :: static            !< true if this is a static var
     logical, allocatable, private                    :: registered        !< true when registered
     logical, allocatable, private                    :: mask_variant      !< If there is a mask variant
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
     character(len=:), allocatable, private           :: err_msg           !< An error message
     character(len=:), allocatable, private           :: interp_method     !< The interp method to be used
                                                            !! when regridding the field in post-processing.
                                                            !! Valid options are "conserve_order1",
                                                            !! "conserve_order2", and "none".
     integer, allocatable, dimension(:), private      :: frequency         !< specifies the frequency
     integer, allocatable, dimension(:), private      :: output_units
     integer, allocatable, private                    :: t
     integer, allocatable, private                    :: tile_count        !< The number of tiles
     integer, allocatable, dimension(:), private      :: axis_ids          !< variable axis IDs
     integer, allocatable, private                    :: area, volume      !< The Area and Volume
     class(*), allocatable, private                   :: missing_value     !< The missing fill value
     class(*), allocatable, private                   :: data_RANGE        !< The range of the variable data
     type (diag_axis_type), allocatable, dimension(:) :: axis              !< The axis object
!> \brief Extends the variable object to work with multiple types of data
     class(*), allocatable :: vardata0
     class(*), allocatable, dimension(:) :: vardata1
     class(*), allocatable, dimension(:,:) :: vardata2
     class(*), allocatable, dimension(:,:,:) :: vardata3
     class(*), allocatable, dimension(:,:,:,:) :: vardata4
     class(*), allocatable, dimension(:,:,:,:,:) :: vardata5



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
     procedure :: get_err_msg
     procedure :: get_interp_method
     procedure :: get_frequency
     procedure :: get_output_units
     procedure :: get_t
     procedure :: get_tile_count
     procedure :: get_axis_ids
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

!type(fmsDiagObject_type) :: diag_object_placeholder (10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
public :: fmsDiagObject_type
public :: null_ob
public :: copy_diag_obj, fms_diag_get_id
public :: fms_diag_object_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fms_diag_object_init (mlv,mlm)
 integer, intent(in) :: mlv !< The maximum length of the varname
 integer, intent(in) :: mlm !< The maximum length of the metadata
!> Get info from the namelist
 MAX_LEN_VARNAME = mlv
 MAX_LEN_META = mlm
!> Initialize the null_d variables
 null_ob%diag_id = DIAG_NULL
end subroutine fms_diag_object_init
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
       (dobj, modname, varname, axes, init_time, &
       longname, units, missing_value, varRange, mask_variant, standname, &
       do_not_log, err_msg, interp_method, tile_count, area, volume, realm, metadata)
 class(fmsDiagObject_type)     , intent(inout)            :: dobj
 CHARACTER(len=*), INTENT(in) :: modname !< The module name
 CHARACTER(len=*), INTENT(in) :: varname !< The variable name
 INTEGER, INTENT(in) :: axes(:) !< The axes indicies
 TYPE(time_type), INTENT(in) :: init_time !< Initial time
 CHARACTER(len=*), OPTIONAL, INTENT(in) :: longname !< THe variables long name
 CHARACTER(len=*), OPTIONAL, INTENT(in) :: units !< The units of the variables
 CHARACTER(len=*), OPTIONAL, INTENT(in) :: standname !< The variables stanard name
 class(*), OPTIONAL, INTENT(in) :: missing_value
 class(*), OPTIONAL, INTENT(in) :: varRANGE(2)
 LOGICAL, OPTIONAL, INTENT(in) :: mask_variant
 LOGICAL, OPTIONAL, INTENT(in) :: do_not_log !< if TRUE, field info is not logged
 CHARACTER(len=*), OPTIONAL, INTENT(out):: err_msg !< Error message to be passed back up
 CHARACTER(len=*), OPTIONAL, INTENT(in) :: interp_method !< The interp method to be used when
                                                         !! regridding the field in post-processing.
                                                         !! Valid options are "conserve_order1",
                                                         !! "conserve_order2", and "none".
 INTEGER, OPTIONAL, INTENT(in) :: tile_count !< the number of tiles
 INTEGER, OPTIONAL, INTENT(in) :: area !< diag_field_id containing the cell area field
 INTEGER, OPTIONAL, INTENT(in) :: volume !< diag_field_id containing the cell volume field
 CHARACTER(len=*), OPTIONAL, INTENT(in):: realm !< String to set as the value to the modeling_realm attribute
 character(len=*), optional, intent(in), dimension(:)     :: metadata !< metedata for the variable

!> Fill in information from the register call
  allocate(character(len=MAX_LEN_VARNAME) :: dobj%varname)
  dobj%varname = trim(varname)
  allocate(character(len=len(modname)) :: dobj%modname)
  dobj%modname = trim(modname)
!> Grab the information from the diag_table
!  TO DO:
!  dobj%diag_field = get_diag_table_field(trim(varname))
!  dobj%diag_field = diag_yaml%get_diag_field(
  !! TODO : Discuss design. Is this a premature return that somehow should
  !! indicate a warning or failure to the calling function and/or the log files?
!  if (is_field_type_null(dobj%diag_field)) then
!     dobj%diag_id = diag_not_found
!     dobj%vartype = diag_null
!     return
!  endif

!> get the optional arguments if included and the diagnostic is in the diag table
  if (present(longname)) then
     allocate(character(len=len(longname)) :: dobj%longname)
     dobj%longname = trim(longname)
  endif
  if (present(standname)) then
     allocate(character(len=len(standname)) :: dobj%standname)
     dobj%standname = trim(standname)
  endif
  if (present(units)) then
     allocate(character(len=len(units)) :: dobj%units)
     dobj%units = trim(units)
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

!     write(6,*)"IKIND for diag_fields(1) is",dobj%diag_fields(1)%ikind
!     write(6,*)"IKIND for "//trim(varname)//" is ",dobj%diag_field%ikind
!> Set the registered flag to true
 dobj%registered = .true.
 ! save it in the diag object container.

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
!     type (diag_fields_type)                           :: diag_field         !< info from diag_table
!     type (diag_files_type),allocatable, dimension(:)  :: diag_file          !< info from diag_table

     objout%diag_id = objin%diag_id

!     class (fms_io_obj), allocatable, dimension(:)    :: fms_fileobj        !< fileobjs
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
function get_metadata (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable, dimension(:) :: rslt
     if (allocated(obj%metadata)) then
       allocate(character(len=(len(obj%metadata(1)))) :: rslt (size(obj%metadata)) )
       rslt = obj%metadata
     else
       allocate(character(len=1) :: rslt(1:1))
       rslt = diag_null_string
     endif
end function get_metadata
!> @brief Gets static
!! @return copy of variable static
function get_static (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     logical :: rslt 
     rslt = obj%static
end function get_static
!> @brief Gets regisetered
!! @return copy of registered
function get_registered (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     logical :: rslt 
     rslt = obj%registered
end function get_registered
!> @brief Gets mask variant
!! @return copy of mask variant
function get_mask_variant (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     logical :: rslt 
     rslt = obj%mask_variant
end function get_mask_variant
!> @brief Gets local
!! @return copy of local
function get_local (obj) &
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
function get_vartype (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer :: rslt 
     rslt = obj%vartype
end function get_vartype
!> @brief Gets varname
!! @return copy of the variable name
function get_varname (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     rslt = obj%varname
end function get_varname
!> @brief Gets longname
!! @return copy of the variable long name or a single string if there is no long name
function get_longname (obj) &
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
function get_standname (obj) &
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
function get_units (obj) &
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
function get_modname (obj) &
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
function get_realm (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     if (allocated(obj%realm)) then
       rslt = obj%realm
     else
       rslt = diag_null_string
     endif
end function get_realm
!> @brief Gets err_msg
!! @return copy of The error message stored in err_msg or an empty string if not allocated
function get_err_msg (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     character(len=:), allocatable :: rslt 
     if (allocated(obj%err_msg)) then
       rslt = obj%err_msg
     else
       rslt = diag_null_string
     endif
end function get_err_msg
!> @brief Gets interp_method
!! @return copy of The interpolation method or an empty string if not allocated
function get_interp_method (obj) &
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
function get_frequency (obj) &
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
function get_output_units (obj) &
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
function get_t (obj) &
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
function get_tile_count (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer :: rslt 
     if (allocated(obj%tile_count)) then
       rslt = obj%tile_count
     else
       rslt = DIAG_NULL
     endif
end function get_tile_count
!> @brief Gets axis_ids
!! @return copy of The axis IDs array or a diag_null if no axis IDs are set
function get_axis_ids (obj) &
result(rslt)
     class (fmsDiagObject_type), intent(in) :: obj !< diag object
     integer, allocatable, dimension(:) :: rslt 
     if (allocated(obj%axis_ids)) then
       allocate(rslt(size(obj%axis_ids)))
       rslt = obj%axis_ids
     else
       allocate(rslt(1))
       rslt = diag_null
     endif
end function get_axis_ids
!> @brief Gets area
!! @return copy of the area or diag_null if not allocated
function get_area (obj) &
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
function get_volume (obj) &
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
     class(*),allocatable :: rslt 
     if (allocated(obj%data_RANGE)) then
       select type (r => obj%data_RANGE)
         type is (integer(kind=i4_kind))
             allocate (integer(kind=i4_kind) :: rslt)
             rslt = r
         type is (integer(kind=i8_kind))
             allocate (integer(kind=i8_kind) :: rslt)
             rslt = r
         type is (real(kind=r4_kind))
             allocate (integer(kind=i4_kind) :: rslt)
             rslt = r
         type is (real(kind=r8_kind))
             allocate (integer(kind=i4_kind) :: rslt)
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


end module fms_diag_object_mod
