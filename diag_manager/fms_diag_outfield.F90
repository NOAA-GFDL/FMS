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

!> @defgroup fms_diag_outfield_mod fms_diag_outfield_mod
!> @ingroup diag_manager
!> @brief fms_diag_outfield_mod defines data types and utility or auxiliary routines
!! useful in updating the output buffer.
!!
!> @author Miguel Zuniga
!!
!! <TT>fms_diag_outfield_mod</TT> The output buffer updating routines are passed configuration
!!  and control data with types defined in this module; and some utility functions called by the
!! updating routines are
!! defined here.
!!
!> @file
!> @brief File for @ref fms_diag_outfield_mod
!> @addtogroup fms_diag_outfield_mod
!> @{
MODULE fms_diag_outfield_mod
  USE platform_mod
  USE mpp_mod, only :FATAL, WARNING
  USE fms_mod, only :lowercase, uppercase,  error_mesg, fms_error_handler


  !! TODO: these might need removal or replacement
  USE diag_data_mod, only:Time_zero
  USE diag_data_mod, only: GLO_REG_VAL, GLO_REG_VAL_ALT, region_out_use_alt_value, VERY_LARGE_AXIS_LENGTH, coord_type
  USE diag_data_mod, only: fmsDiagIbounds_type, input_field_type, output_field_type
  USE fms_diag_time_reduction_mod, only: fmsDiagTimeReduction_type, time_none , time_average, time_rms
  USE fms_diag_time_reduction_mod, only:  time_max, time_min, time_sum, time_power

  implicit none

  !> @brief Class fmsDiagOutfield_type (along with class ms_diag_outfield_index_type )
  !! contain information used in updating the output buffers by the diag_manager
  !! send_data routines. In some sense they can be seen as encapsulating related
  !! information in a convenient way (e.g. to pass to functions and for do loop
  !! controls.)
  !!
  !! Class fmsDiagOutfield_type also contains a significant subset of the fields
  !! and routines of of the legacy class output_field_type
  !! TODO: (MDM) This class will need further development for the modern_diag effort.
  !! For its development, consider the legacy diag_util::init_output_field already
  !! in place. Fields added so are used the the field buffer math/dmUpdate functions.
  !! TODO (MDM) : Should the MDM have pow_value be type REAL?
  !> @ingroup fms_diag_outfield_mod
  TYPE, public :: fmsDiagOutfield_type
     PRIVATE
     CHARACTER(len=:), ALLOCATABLE :: module_name !< Module name.
     CHARACTER(len=:), ALLOCATABLE :: field_name  !< Output field name.
     CHARACTER(len=:), ALLOCATABLE :: output_name !< Output name written to file.
     CHARACTER(len=:), ALLOCATABLE :: output_file !< File where field should be written.

     !!Major outer loop controls in send_data functions.
     INTEGER :: pow_value     !< Power value for rms or pow(x) calculations
     LOGICAL :: phys_window   !< TODO: Rename? OMP subsetted data, See output_fields
     LOGICAL :: need_compute  !< True iff is local_output and current PE take part in send_data.
     LOGICAL :: reduced_k_range !< If true, the local start and end indecies are used in k (i.e. 3rd) dim.
     LOGICAL :: missvalue_present !<
     LOGICAL :: mask_variant
     LOGICAL :: mask_present !< True iff mask argument is present in user-facing send function call.
     !< Note this field exists since the actual mask argument in the send
     !< function call may be downstream replaced by a null pointer which
     !< is considered present.

     TYPE(fmsDiagTimeReduction_type) :: time_reduction !< Instance of the fmsDiagTimeTeduction_type.

     !!TODO (Future effort? ) : a pointer for time_min and time_max comparison function
     !! If possible, this can remove the innermost if/then/else construct in the buffer update loops.
     !!       min_max_f_ptr => (should point to < or > operators)

     !! gcc error: Interface ‘addwf’ at (1) must be explicit
     ! procedure (addwf), pointer, nopass :: f_ptr => null () !!A pointer to the field weighing procedure

   CONTAINS
     procedure :: get_module_name
     procedure :: get_field_name
     procedure :: get_output_name
     procedure :: get_output_file
     procedure :: get_pow_value
     procedure :: get_phys_window
     procedure :: get_need_compute
     procedure :: get_reduced_k_range
     procedure :: get_missvalue_present
     procedure :: get_mask_variant
     procedure :: get_mask_present
     procedure :: get_time_reduction
     procedure, public  :: initialize => initialize_outfield_imp
     procedure :: initialize_for_ut

  END TYPE fmsDiagOutfield_type

  !> @brief Class fms_diag_outfield_index_type which (along with class fmsDiagOutfield_type)
  !! encapsulate related information used in updating the output buffers by the diag_manager
  !! send_data routines. This class in particular focuses on do loop index controls or settings.
  !! Note that the index names in this class should be indentical to the names used in the
  !! diag_manager send_data functions and in the "math" buffer update functions. The purpose
  !! of this class is also to allow for a smaller call function signature for the math/buffer
  !! update functions.
  !> @ingroup fms_diag_outfield_mod
  TYPE, public :: fmsDiagOutfieldIndex_type
     PRIVATE
     INTEGER :: f1,f2 !< Indecies used specify 1st dim bounds of field, mask and rmask.
     INTEGER :: f3,f4 !< Indecies used specify 2st dim bounds of field, mask and rmask.
     INTEGER :: is, js, ks  !< Start indecies in each spatial dim of the field_data; and
                                                    !! may be user provided in send_data
     Integer :: ie, je, ke  !< End indecies in each spatial dim of the field_data; and
                                                     !! may be user provided in send_data
     INTEGER :: hi !< halo size in x direction. Same name as in send_data
     INTEGER :: hj !< halo size in y direction. Same
   CONTAINS
     procedure :: initialize => initialize_outfield_index_type
     procedure :: get_f1
     procedure :: get_f2
     procedure :: get_f3
     procedure :: get_f4
     procedure :: get_is
     procedure :: get_js
     procedure :: get_ks
     procedure :: get_ie
     procedure :: get_je
     procedure :: get_ke
     procedure :: get_hi
     procedure :: get_hj
  END TYPE fmsDiagOutfieldIndex_type

CONTAINS

  !> @brief Gets module_name
  !! @return copy of the module_name character array
  pure function get_module_name (this) &
       result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< diag object
    character(len=:), allocatable :: rslt
    rslt = this%module_name
  end function get_module_name

  !> @brief Gets field_name
  !! @return copy of the field_name character array
  pure function get_field_name (this) &
       result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< diag object
    character(len=:), allocatable :: rslt
    rslt = this%field_name
  end function get_field_name

  !> @brief Gets output_name
  !! @return copy of the output_name character array
  pure function get_output_name (this) &
       result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< diag object
    character(len=:), allocatable :: rslt
    rslt = this%output_name
  end function get_output_name

  !> @brief Gets output_file
  !! @return copy of the output_file character array
  pure function get_output_file (this) &
       result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< diag object
    character(len=:), allocatable :: rslt
    rslt = this%output_file
  end function get_output_file

  !> @brief Gets pow_value
  !! @return copy of integer member pow_value
  pure integer function get_pow_value (this) result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< The fmsDiagOutfield_type object
    rslt = this%pow_value
  end function get_pow_value

  !> @brief Gets phys_window
  !! @return copy of integer member phys_window
  pure logical function get_phys_window (this) result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< The fmsDiagOutfield_type object
    rslt = this%phys_window
  end function get_phys_window

  !> @brief Gets need_compute
  !! @return copy of integer member need_compute
  pure logical function get_need_compute (this) result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< The fmsDiagOutfield_type object
    rslt = this%need_compute
  end function get_need_compute

  !> @brief Gets reduced_k_range
  !! @return copy of integer member reduced_k_range
  pure logical function get_reduced_k_range (this) result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< The fmsDiagOutfield_type object
    rslt = this%reduced_k_range
  end function get_reduced_k_range

  !> @brief Gets missvalue_present
  !! @return copy of integer member missvalue_present
  pure logical function get_missvalue_present (this) result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< The fmsDiagOutfield_type object
    rslt = this%missvalue_present
  end function get_missvalue_present

  !> @brief Gets mask_variant
  !! @return copy of integer member mask_variant
  pure logical function get_mask_variant (this) result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< The fmsDiagOutfield_type object
    rslt = this%mask_variant
  end function get_mask_variant

  !> @brief Gets mask_present
  !! @return copy of integer member mask_present
  pure logical function get_mask_present (this) result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< The fmsDiagOutfield_type object
    rslt = this%mask_present
  end function get_mask_present

  !> @brief Gets the time_reduction object
  !! @return copy of the memeber object time_reduction
  function get_time_reduction (this) result(rslt)
    class (fmsDiagOutfield_type), intent(in) :: this !< diag object
    TYPE(fmsDiagTimeReduction_type),  allocatable :: rslt
    allocate( rslt )
    call rslt%copy(this%time_reduction)
  end function get_time_reduction

  !> @brief Gets f1
  !! @return copy of integer member f1
  pure integer function get_f1 (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%f1
  end function get_f1

  !> @brief Gets f2
  !! @return copy of integer member f2
  pure integer function get_f2 (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%f2
  end function get_f2

  !> @brief Gets f3
  !! @return copy of integer member f3
  pure integer function get_f3 (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%f3
  end function get_f3

  !> @brief Gets f4
  !! @return copy of integer member f4
  pure integer function get_f4 (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%f4
  end function get_f4

  !> @brief Gets is
  !! @return copy of integer member is
  pure integer function get_is (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%is
  end function get_is

  !> @brief Gets js
  !! @return copy of integer member js
  pure integer function get_js (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%js
  end function get_js

  !> @brief Gets ks
  !! @return copy of integer member ks
  pure integer function get_ks (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%ks
  end function get_ks

  !> @brief Gets ie
  !! @return copy of integer member ie
  pure integer function get_ie (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%ie
  end function get_ie

  !> @brief Gets je
  !! @return copy of integer member je
  pure integer function get_je (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%je
  end function get_je

  !> @brief Gets ke
  !! @return copy of integer member ke
  pure integer function get_ke (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%ke
  end function get_ke

  !> @brief Gets hi
  !! @return copy of integer member hi
  pure integer function get_hi (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%hi
  end function get_hi

  !> @brief Gets hj
  !! @return copy of integer member hj
  pure integer function get_hj (this) result(rslt)
    class (fmsDiagOutfieldIndex_type), intent(in) :: this !< The fmsDiagOutfieldIndex_type object
    rslt = this%hj
  end function get_hj


  !> #brief initialize all the members of the class.
  SUBROUTINE initialize_outfield_index_type(this, is, js , ks, ie, je, ke, hi, hj, f1, f2, f3, f4)
    CLASS(fmsDiagOutfieldIndex_type), INTENT(inout)  :: this
    INTEGER, INTENT(in) :: is, js, ks !< Variable used to update class member of same names.
    INTEGER, INTENT(in) :: ie, je, ke !< Variable used to update class member of same names.
    INTEGER, INTENT(in) :: hi, hj !< Variable used to update class member of same names.
    INTEGER, INTENT(in) :: f1, f2, f3, f4 !< Variable used to update class member of same names.

    this%is = is
    this%js = js
    this%ks = ks
    this%ie = ie
    this%je = je
    this%ke = ke

    this%hi = hi
    this%hj = hj

    this%f1 = f1
    this%f2 = f2
    this%f3 = f3
    this%f4 = f4
  END SUBROUTINE initialize_outfield_index_type


  !> @brief Update the fmsDiagOutfield_type instance with those fields used in the legacy diag manager.
  !! Note that this is initializing from the legacy structures.
  !! Note that output_frequency  came from file_type;
  SUBROUTINE initialize_outfield_imp(this, input_field,  output_field, mask_present, freq)
    CLASS(fmsDiagOutfield_type), INTENT(inout) :: this !< An instance of the fmsDiagOutfield_type
    TYPE(input_field_type),     INTENT(in) :: input_field !< An instance of the input_field_type
    TYPE(output_field_type),    INTENT(in) :: output_field !< An instance of the output_field_type
    LOGICAL, INTENT(in) :: mask_present !< Was the mask present in the call to send_data?
    INTEGER, INTENT(in) :: freq !< The output frequency.
    INTEGER ::  time_redux !< The time reduction type integer.

    this%module_name = input_field%module_name
    this%field_name = input_field%field_name
    this%output_name = output_field%output_name

    this%pow_value = output_field%pow_value
    this%phys_window = output_field%phys_window
    this%need_compute = output_field%need_compute
    this%reduced_k_range = output_field%reduced_k_range
    this%mask_variant = input_field%mask_variant
    !!Note: in legacy diag manager, presence of missing value vs presence of mask
    !! is determined in different ways (diag table vs send function call)
    this%missvalue_present = input_field%missing_value_present
    this%mask_present = mask_present

    time_redux = get_output_field_time_reduction (output_field)
    call this%time_reduction%initialize( time_redux , freq)

    !!TODO: the time_min and time_max buffer update code is almost the exact same src code, except
    !!  for the compariosn function. Simplify code and set comparison function:
    !!TODO: If possible add to the power function. See issue with pointers and elemental functions

  END SUBROUTINE initialize_outfield_imp

  !> @brief Initialized an fmsDiagOutfield_type as needed for unit tests.
  subroutine initialize_for_ut(this, module_name, field_name, output_name, &
       &  power_val, phys_window, need_compute, mask_variant,  reduced_k_range, num_elems, &
       & time_reduction_type,output_freq)
    CLASS(fmsDiagOutfield_type), intent(inout)  :: this
    CHARACTER(len=*), INTENT(in) :: module_name !< Var with same name in fmsDiagOutfield_type
    CHARACTER(len=*), INTENT(in) :: field_name !< Var with same name in fmsDiagOutfield_type
    CHARACTER(len=*), INTENT(in) :: output_name !< Var with same name in fmsDiagOutfield_type
    INTEGER, INTENT(in) :: power_val    !< Var with same name in fmsDiagOutfield_type
    LOGICAL, INTENT(in) :: phys_window  !< Var with same name in fmsDiagOutfield_type
    LOGICAL, INTENT(in) :: need_compute  !< Var with same name in fmsDiagOutfield_type
    LOGICAL, INTENT(in) :: mask_variant  !< Var with same name in fmsDiagOutfield_type
    LOGICAL, INTENT(in) :: reduced_k_range !< Var with same name in fmsDiagOutfield_type
    INTEGER, INTENT(in) :: num_elems !< Var with same name in fmsDiagOutfield_type
    INTEGER, INTENT(in) :: time_reduction_type !< Var with same name in fmsDiagOutfield_type
    INTEGER, INTENT(in) :: output_freq !< The output_freq need in initaliztion of time_reduction_type

    this%module_name = module_name
    this%field_name = field_name
    this%output_name = output_name
    this%pow_value = power_val
    this%phys_window = phys_window
    this%need_compute = need_compute
    this%reduced_k_range = reduced_k_range
    this%mask_variant = mask_variant
    call this%time_reduction%initialize(time_reduction_type, output_freq)
  end subroutine initialize_for_ut

  !> @brief Reset the time reduction member field. Intended for use in unit tests only.
  SUBROUTINE reset_time_reduction_ut(this,  source )
    CLASS(fmsDiagOutfield_type), INTENT(inout) :: this !< An instance of the fmsDiagOutfield_type
    TYPE(fmsDiagTimeReduction_type) :: source !< The fmsDiagTimeReduction_type to copy from
    call this%time_reduction%copy(source)
  END SUBROUTINE reset_time_reduction_ut



  !> \brief Get the time reduction from a legacy output field.
  !\note Note we do not place this in the time_reduction class to avoid circular dependencies.
  function get_output_field_time_reduction(ofield) result (rslt)
    TYPE(output_field_type), INTENT(in) :: ofield !< An instance of the output_field_type
    INTEGER ::  rslt !< The result integer which is the time reduction.
    if(ofield%time_max) then
       rslt = time_max
    elseif(ofield%time_min)then
       rslt = time_min
    else if (ofield%time_sum) then
       rslt = time_sum
    else if (ofield%time_rms) then
       rslt = time_rms
    else if (ofield%time_average) then
       rslt = time_average
    else
       rslt = time_none
       !if(.NOT. ofield%static) then
       !!TODO: Set error to FATAL. When legacy diag_manager is removed?
       !   CALL error_mesg('fms_diag_outfield:get_output_field_time_reduction', &
       !      & 'result is time_none but out_field%static is not true', WARNING)
       !end if
    endif
  end function get_output_field_time_reduction

END MODULE fms_diag_outfield_mod
!> @}
! close documentation grouping


