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
   USE diag_data_mod, only: fms_diag_ibounds_type, input_field_type, output_field_type
   USE fms_diag_time_reduction_mod, only: time_reduction_type, time_none , time_average, time_rms
   USE fms_diag_time_reduction_mod, only:  time_max, time_min, time_sum, time_power

   !!TODO: for modern diag: if use_yaml then
   !! USE fms_diag_yaml_mod, only : diagYamlFiles_type, diagYamlFilesVar_type
   !!USE fms_diag_field_object_mod, only: fmsDiagField_type
   !!USE diag_data_mod, only: fms_diag_buff_intervals_t, diag_grid
   !!USE time_manager_mod,ONLY: time_type, OPERATOR(==), OPERATOR(>), NO_CALENDAR, increment_date,&
   !!& increment_time, get_calendar_type, get_date, get_time, leap_year, OPERATOR(-),&
   !!& OPERATOR(<), OPERATOR(>=), OPERATOR(<=), OPERATOR(==)

   implicit none

   !> @brief Class fms_diag_outfield_type (along with class ms_diag_outfield_index_type )
   !! contain information used in updating the output buffers by the diag_manager
   !! send_data routines. In some sense they can be seen as encapsulating related
   !! information in a convenient way (e.g. to pass to functions and for do loop
   !! controls.
   !!
   !! Class fms_diag_outfield_type also contains a significant subset of the fields
   !! and routines of of the legacy class output_field_type
   !! TODO: Developemnt of this class is in a seperate and future PR. For its development,
   !! consider the legacy diag_util::init_output_field already in place. Fields added so
   !! are uesd the the field buffer math/dupdate functions.
   !> @ingroup fms_diag_outfield_mod
   TYPE fms_diag_outfield_type
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
      LOGICAL :: mask_present !< True iff mars arguemnt is present in user-facing send function call.
                              !< Note this field exist since the actual mask argument in the send
                              !< function call may be downstream replaced by a null pointer which
                              !< is considered present.

      TYPE(time_reduction_type) :: time_reduction !< Instance of the time_reduction_type.

      !!TODO : a pointer for time_min and time_max comparison function
      !!       min_max_f_ptr => (should point to < or > operators)

      !! gcc error: Interface ‘addwf’ at (1) must be explicit
      ! procedure (addwf), pointer, nopass :: f_ptr => null () !!A pointer to the field weighing procedure

   CONTAINS
      procedure, public  :: initialize => initialize_outfield_imp
   END TYPE fms_diag_outfield_type


   !> @brief Class fms_diag_outfield_index_type which (along with class fms_diag_outfield_type)
   !! encapsulate related information used in updating the output buffers by the diag_manager
   !! send_data routines. This class in particular focuses on do loop index controls or settings.
   !! Note that the index names in this class should be indentical to the names used in the
   !! diag_manager send_data functions and in the "math" buffer update functions. The purpose
   !! of this class is also to allow for a smaller call function signature for the math/buffer
   !! update functions.
   !> @ingroup fms_diag_outfield_mod
   TYPE, public :: fms_diag_outfield_index_type
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
   END TYPE fms_diag_outfield_index_type

CONTAINS
   !!TODO: In the modern diag, the field_val and weight may also be of integer type,
   !!      and so may need to use the pre-processor.


!> #brief initialize all the memebers of the class.
   SUBROUTINE initialize_outfield_index_type(this, is, js , ks, ie, je, ke, hi, hj, f1, f2, f3, f4)
      CLASS(fms_diag_outfield_index_type), INTENT(inout)  :: this
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


   !> @brief Update with those fields used in the legacy diag manager.
  !! Note that this is initializing from the legacy structures.
   !! Note that output_frequency  came from file_type;
   SUBROUTINE initialize_outfield_imp(this, input_field,  output_field, mask_present, freq)
      CLASS(fms_diag_outfield_type), INTENT(inout) :: this
      TYPE(input_field_type),     INTENT(in) :: input_field
      TYPE(output_field_type),    INTENT(in) :: output_field
      LOGICAL, INTENT(in) :: mask_present
      INTEGER, INTENT(in) :: freq
      INTEGER ::  time_redux

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


   !> \brief Get the time reduction from a legacy output field.
    !! Note we do not place this in the time_reduction class to avoid circular dependencies.
   function get_output_field_time_reduction(ofield) result (rslt)
    TYPE(output_field_type), INTENT(in) :: ofield
    INTEGER ::  rslt
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
      if(.NOT. ofield%static) then
         !!TODO: Set error to FATAL. When legacy diag_manager is removed?
         CALL error_mesg('fms_diag_outfield:get_output_field_time_reduction', &
            & 'result is time_none but out_field%static is not true', WARNING)
      end if
    endif
    end function get_output_field_time_reduction

END MODULE fms_diag_outfield_mod
!> @}
! close documentation grouping


