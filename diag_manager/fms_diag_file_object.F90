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
!> @defgroup fms_diag_output_yaml_mod fms_diag_output_yaml_mod
!> @ingroup diag_manager
!! @brief fms_diag_file_object_mod handles the file objects data, functions, and subroutines.
!! @author Tom Robinson
!! @description The fmsDiagFile_type contains the information for each history file to be written.  It has
!! a pointer to the information from the diag yaml, additional metadata that comes from the model, and a
!! list of the variables and their variable IDs that are in the file.
module fms_diag_file_object_mod
#ifdef use_yaml
use fms2_io_mod, only: FmsNetcdfFile_t, FmsNetcdfUnstructuredDomainFile_t, FmsNetcdfDomainFile_t, &
                       get_instance_filename, open_file, close_file, get_mosaic_tile_file, unlimited, &
                       register_axis, register_field, register_variable_attribute, write_data
use diag_data_mod, only: DIAG_NULL, NO_DOMAIN, max_axes, SUB_REGIONAL, get_base_time, DIAG_NOT_REGISTERED, &
                         TWO_D_DOMAIN, UG_DOMAIN, prepend_date, DIAG_DAYS, VERY_LARGE_FILE_FREQ, &
                         get_base_year, get_base_month, get_base_day, get_base_hour, get_base_minute, &
                         get_base_second, time_unit_list, time_average, time_rms, time_max, time_min, time_sum, &
                         time_diurnal, time_power, time_none
use time_manager_mod, only: time_type, operator(>), operator(/=), operator(==), get_date, get_calendar_type, &
                            VALID_CALENDAR_TYPES, operator(>=), date_to_string
use fms_diag_time_utils_mod, only: diag_time_inc, get_time_string, get_date_dif
use fms_diag_yaml_mod, only: diag_yaml, diagYamlObject_type, diagYamlFiles_type, subRegion_type, diagYamlFilesVar_type
use fms_diag_axis_object_mod, only: diagDomain_t, get_domain_and_domain_type, fmsDiagAxis_type, &
                                    fmsDiagAxisContainer_type, DIAGDOMAIN2D_T, DIAGDOMAINUG_T, &
                                    fmsDiagFullAxis_type, define_subaxis
use mpp_mod, only: mpp_get_current_pelist, mpp_npes, mpp_root_pe, mpp_pe, mpp_error, FATAL, stdout, &
                   uppercase, lowercase

implicit none
private

public :: fmsDiagFileContainer_type
public :: fmsDiagFile_type, fms_diag_files_object_init, fms_diag_files_object_initialized

logical :: fms_diag_files_object_initialized = .false.

integer, parameter :: var_string_len = 25

type :: fmsDiagFile_type
 private
  integer :: id !< The number associated with this file in the larger array of files
  TYPE(time_type) :: start_time       !< The start time for the file
  TYPE(time_type) :: last_output      !< Time of the last time output was writen
  TYPE(time_type) :: next_output      !< Time of the next write
  TYPE(time_type) :: next_next_output !< Time of the next next write

  !< This will be used when using the new_file_freq keys in the diag_table.yaml
  TYPE(time_type) :: next_open        !< The next time to open the file
  class(FmsNetcdfFile_t), allocatable :: fileobj !< fms2_io file object for this history file
  type(diagYamlFiles_type), pointer :: diag_yaml_file => null() !< Pointer to the diag_yaml_file data
  integer                                      :: type_of_domain !< The type of domain to use to open the file
                                                                 !! NO_DOMAIN, TWO_D_DOMAIN, UG_DOMAIN, SUB_REGIONAL
  class(diagDomain_t), pointer                 :: domain         !< The domain to use,
                                                                 !! null if NO_DOMAIN or SUB_REGIONAL
  character(len=:) , dimension(:), allocatable :: file_metadata_from_model !< File metadata that comes from
                                                                           !! the model.
  integer, dimension(:), allocatable :: field_ids !< Variable IDs corresponding to file_varlist
  logical, dimension(:), private, allocatable :: field_registered   !< Array corresponding to `field_ids`, .true.
                                                                 !! if the variable has been registered and
                                                                 !! `field_id` has been set for the variable
  integer, allocatable                         :: num_registered_fields !< The number of fields registered
                                                                        !! to the file
  integer, dimension(:), allocatable :: axis_ids !< Array of axis ids in the file
  integer, dimension(:), allocatable :: buffer_ids !< array of buffer ids associated with the file
  integer :: number_of_axis !< Number of axis in the file
  logical :: time_ops !< .True. if file contains variables that are time_min, time_max, time_average or time_sum
  integer :: unlimited_dimension !< The unlimited dimension currently being written
  logical :: is_static !< .True. if the frequency is -1

 contains
  procedure, public :: add_field_id
  procedure, public :: has_file_metadata_from_model
  procedure, public :: has_fileobj
  procedure, public :: has_diag_yaml_file
  procedure, public :: set_domain_from_axis
  procedure, public :: set_file_domain
  procedure, public :: add_axes
  procedure, public :: add_start_time
  procedure, public :: set_file_time_ops
  procedure, public :: has_field_ids
  procedure, public :: get_id
! TODO  procedure, public :: get_fileobj ! TODO
! TODO  procedure, public :: get_diag_yaml_file ! TODO
  procedure, public :: get_file_metadata_from_model
  procedure, public :: get_field_ids
! The following fuctions come will use the yaml inquiry functions
 procedure, public :: get_file_fname
 procedure, public :: get_file_frequnit
 procedure, public :: get_file_freq
 procedure, public :: get_file_timeunit
 procedure, public :: get_file_unlimdim
 procedure, public :: get_file_sub_region
 procedure, public :: get_file_new_file_freq
 procedure, public :: get_file_new_file_freq_units
 procedure, public :: get_file_start_time
 procedure, public :: get_file_duration
 procedure, public :: get_file_duration_units
 procedure, public :: get_file_varlist
 procedure, public :: get_file_global_meta
 procedure, public :: has_file_fname
 procedure, public :: has_file_frequnit
 procedure, public :: has_file_freq
 procedure, public :: has_file_timeunit
 procedure, public :: has_file_unlimdim
 procedure, public :: has_file_sub_region
 procedure, public :: has_file_new_file_freq
 procedure, public :: has_file_new_file_freq_units
 procedure, public :: has_file_start_time
 procedure, public :: has_file_duration
 procedure, public :: has_file_duration_units
 procedure, public :: has_file_varlist
 procedure, public :: has_file_global_meta
 procedure, public :: dump_file_obj
end type fmsDiagFile_type

type, extends (fmsDiagFile_type) :: subRegionalFile_type
  integer, dimension(:), allocatable :: sub_axis_ids !< Array of axis ids in the file
  logical :: write_on_this_pe !< Flag indicating if the subregion is on the current PE
  logical :: is_subaxis_defined !< Flag indicating if the subaxes have already been defined
end type subRegionalFile_type

!> \brief A container for fmsDiagFile_type.  This is used to create the array of files
type fmsDiagFileContainer_type
  class (fmsDiagFile_type),allocatable :: FMS_diag_file !< The individual file object

  contains
  procedure :: open_diag_file
  procedure :: write_time_metadata
  procedure :: write_axis_metadata
  procedure :: write_axis_data
  procedure :: writing_on_this_pe
  procedure :: is_time_to_write
  procedure :: write_time_data
  procedure :: update_next_write
  procedure :: increase_unlimited_dimension
  procedure :: close_diag_file
end type fmsDiagFileContainer_type

!type(fmsDiagFile_type), dimension (:), allocatable, target :: FMS_diag_file !< The array of diag files
!class(fmsDiagFileContainer_type),dimension (:), allocatable, target :: FMS_diag_file

contains

!< @brief Allocates the number of files and sets an ID based for each file
!! @return true if there are files allocated in the YAML object
logical function fms_diag_files_object_init (files_array)
  class(fmsDiagFileContainer_type), allocatable, target, intent(inout) :: files_array (:) !< array of diag files
  class(fmsDiagFile_type), pointer :: obj => null() !< Pointer for each member of the array
  integer :: nFiles !< Number of files in the diag yaml
  integer :: i !< Looping iterator
  if (diag_yaml%has_diag_files()) then
   nFiles = diag_yaml%size_diag_files()
   allocate (files_array(nFiles))
   set_ids_loop: do i= 1,nFiles
     !> If the file has a sub_regional, define it as one and allocate the sub_axis_ids array.
     !! This will be set in a add_axes
     if (diag_yaml%diag_files(i)%has_file_sub_region()) then
       allocate(subRegionalFile_type :: files_array(i)%FMS_diag_file)
       obj => files_array(i)%FMS_diag_file
       select type (obj)
         type is (subRegionalFile_type)
           allocate(obj%sub_axis_ids(max_axes))
           obj%sub_axis_ids = diag_null
           obj%write_on_this_pe = .false.
           obj%is_subaxis_defined = .false.
           obj%number_of_axis = 0
       end select
     else
       allocate(FmsDiagFile_type::files_array(i)%FMS_diag_file)
       obj => files_array(i)%FMS_diag_file
     endif
     !!
     obj%diag_yaml_file => diag_yaml%diag_files(i)
     obj%id = i
     allocate(obj%field_ids(diag_yaml%diag_files(i)%size_file_varlist()))
     allocate(obj%field_registered(diag_yaml%diag_files(i)%size_file_varlist()))
     !! Initialize the integer arrays
     obj%field_ids = DIAG_NOT_REGISTERED
     obj%field_registered = .FALSE.
     obj%num_registered_fields = 0

     !> These will be set in a set_file_domain
     obj%type_of_domain = NO_DOMAIN
     obj%domain => null()

     !> This will be set in a add_axes
     allocate(obj%axis_ids(max_axes))
     obj%number_of_axis = 0

     !> Set the start_time of the file to the base_time and set up the *_output variables
     obj%start_time = get_base_time()
     obj%last_output = get_base_time()
     obj%next_output = diag_time_inc(obj%start_time, obj%get_file_freq(), obj%get_file_frequnit())
     obj%next_next_output = diag_time_inc(obj%next_output, obj%get_file_freq(), obj%get_file_frequnit())
     obj%next_open = get_base_time()
     obj%time_ops = .false.
     obj%unlimited_dimension = 0
     obj%is_static = obj%get_file_freq() .eq. -1

     nullify(obj)
   enddo set_ids_loop
   fms_diag_files_object_init = .true.
  else
   fms_diag_files_object_init = .false.
!  mpp_error("fms_diag_files_object_init: The diag_table.yaml file has not been correctly parsed.",&
!    FATAL)
  endif
end function fms_diag_files_object_init

!> \brief Adds a field ID to the file
subroutine add_field_id (this, new_field_id)
  class(fmsDiagFile_type), intent(inout) :: this !< The file object
  integer, intent(in) :: new_field_id !< The field ID to be added to field_ids
  this%num_registered_fields = this%num_registered_fields + 1
  if (this%num_registered_fields .le. size(this%field_ids)) then
    this%field_ids( this%num_registered_fields ) = new_field_id
    this%field_registered( this%num_registered_fields ) = .true.
  else
    call mpp_error(FATAL, "The file: "//this%get_file_fname()//" has already been assigned its maximum "//&
                 "number of fields.")
  endif
end subroutine add_field_id

!> \brief Set the time_ops variable in the diag_file object
subroutine set_file_time_ops(this, VarYaml, is_static)
  class(fmsDiagFile_type),      intent(inout) :: this      !< The file object
  type (diagYamlFilesVar_type), intent(in)    :: VarYaml   !< The variable's yaml file
  logical,                      intent(in)    :: is_static !< Flag indicating if variable is static

  if (this%time_ops) then
    if (is_static) return
    if (VarYaml%get_var_reduction() .eq. time_none) then
      call mpp_error(FATAL, "The file: "//this%get_file_fname()//&
                            " has variables that are time averaged and instantaneous")
    endif
  else
    select case (VarYaml%get_var_reduction())
      case (time_average, time_rms, time_max, time_min, time_sum, time_diurnal, time_power)
        this%time_ops = .true.
    end select
  endif

end subroutine set_file_time_ops

!> \brief Logical function to determine if the variable file_metadata_from_model has been allocated or associated
!! \return .True. if file_metadata_from_model exists .False. if file_metadata_from_model has not been set
pure logical function has_file_metadata_from_model (this)
  class(fmsDiagFile_type), intent(in) :: this !< The file object
  has_file_metadata_from_model = allocated(this%file_metadata_from_model)
end function has_file_metadata_from_model

!> \brief Logical function to determine if the variable fileobj has been allocated or associated
!! \return .True. if fileobj exists .False. if fileobj has not been set
pure logical function has_fileobj (this)
  class(fmsDiagFile_type), intent(in) :: this !< The file object
  has_fileobj = allocated(this%fileobj)
end function has_fileobj

!> \brief Logical function to determine if the variable diag_yaml_file has been allocated or associated
!! \return .True. if diag_yaml_file exists .False. if diag_yaml has not been set
pure logical function has_diag_yaml_file (this)
  class(fmsDiagFile_type), intent(in) :: this !< The file object
  has_diag_yaml_file = associated(this%diag_yaml_file)
end function has_diag_yaml_file

!> \brief Logical function to determine if the variable field_ids has been allocated or associated
!! \return .True. if field_ids exists .False. if field_ids has not been set
pure logical function has_field_ids (this)
  class(fmsDiagFile_type), intent(in) :: this !< The file object
  has_field_ids = allocated(this%field_ids)
end function has_field_ids

!> \brief Returns a copy of the value of id
!! \return A copy of id
pure function get_id (this) result (res)
  class(fmsDiagFile_type), intent(in) :: this !< The file object
  integer :: res
  res = this%id
end function get_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TODO
!> \brief Returns a copy of the value of fileobj
!! \return A copy of fileobj
!pure function get_fileobj (obj) result (res)
!  class(fmsDiagFile_type), intent(in) :: obj !< The file object
!  class(FmsNetcdfFile_t) :: res
!  res = obj%fileobj
!end function get_fileobj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TODO
!!> \brief Returns a copy of the value of diag_yaml_file
!!! \return A copy of diag_yaml_file
!pure function get_diag_yaml_file (obj) result (res)
!  class(fmsDiagFile_type), intent(in) :: obj !< The file object
!  type(diagYamlFiles_type) :: res
!  res = obj%diag_yaml_file
!end function get_diag_yaml_file

!> \brief Returns a copy of the value of file_metadata_from_model
!! \return A copy of file_metadata_from_model
pure function get_file_metadata_from_model (this) result (res)
  class(fmsDiagFile_type), intent(in) :: this !< The file object
  character(len=:), dimension(:), allocatable :: res
  res = this%file_metadata_from_model
end function get_file_metadata_from_model

!> \brief Returns a copy of the value of field_ids
!! \return A copy of field_ids
pure function get_field_ids (this) result (res)
  class(fmsDiagFile_type), intent(in) :: this !< The file object
  integer, dimension(:), allocatable :: res
  allocate(res(size(this%field_ids)))
  res = this%field_ids
end function get_field_ids

!!!!!!!!! Functions from diag_yaml_file
!> \brief Returns a copy of file_fname from the yaml object
!! \return Copy of file_fname
pure function get_file_fname (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 character (len=:), allocatable :: res
  res = this%diag_yaml_file%get_file_fname()
end function get_file_fname

!> \brief Returns a copy of file_frequnit from the yaml object
!! \return Copy of file_frequnit
pure function get_file_frequnit (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 integer :: res
  res = this%diag_yaml_file%get_file_frequnit()
end function get_file_frequnit

!> \brief Returns a copy of file_freq from the yaml object
!! \return Copy of file_freq
pure function get_file_freq (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 integer :: res
  res = this%diag_yaml_file%get_file_freq()
end function get_file_freq

!> \brief Returns a copy of file_timeunit from the yaml object
!! \return Copy of file_timeunit
pure function get_file_timeunit (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 integer :: res
  res = this%diag_yaml_file%get_file_timeunit()
end function get_file_timeunit

!> \brief Returns a copy of file_unlimdim from the yaml object
!! \return Copy of file_unlimdim
pure function get_file_unlimdim (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 character (len=:), allocatable :: res
  res = this%diag_yaml_file%get_file_unlimdim()
end function get_file_unlimdim

!> \brief Returns a copy of file_sub_region from the yaml object
!! \return Copy of file_sub_region
function get_file_sub_region (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 type(subRegion_type) :: res
  res = obj%diag_yaml_file%get_file_sub_region()
end function get_file_sub_region

!> \brief Returns a copy of file_new_file_freq from the yaml object
!! \return Copy of file_new_file_freq
pure function get_file_new_file_freq (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 integer :: res
  res = this%diag_yaml_file%get_file_new_file_freq()
end function get_file_new_file_freq

!> \brief Returns a copy of file_new_file_freq_units from the yaml object
!! \return Copy of file_new_file_freq_units
pure function get_file_new_file_freq_units (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 integer :: res
  res = this%diag_yaml_file%get_file_new_file_freq_units()
end function get_file_new_file_freq_units

!> \brief Returns a copy of file_start_time from the yaml object
!! \return Copy of file_start_time
pure function get_file_start_time (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 character (len=:), allocatable :: res
  res = this%diag_yaml_file%get_file_start_time()
end function get_file_start_time

!> \brief Returns a copy of file_duration from the yaml object
!! \return Copy of file_duration
pure function get_file_duration (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 integer :: res
  res = this%diag_yaml_file%get_file_duration()
end function get_file_duration

!> \brief Returns a copy of file_duration_units from the yaml object
!! \return Copy of file_duration_units
pure function get_file_duration_units (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 integer :: res
  res = this%diag_yaml_file%get_file_duration_units()
end function get_file_duration_units

!> \brief Returns a copy of file_varlist from the yaml object
!! \return Copy of file_varlist
pure function get_file_varlist (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 character (len=:), allocatable, dimension(:) :: res
  res = this%diag_yaml_file%get_file_varlist()
end function get_file_varlist

!> \brief Returns a copy of file_global_meta from the yaml object
!! \return Copy of file_global_meta
pure function get_file_global_meta (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 character (len=:), allocatable, dimension(:,:) :: res
  res = this%diag_yaml_file%get_file_global_meta()
end function get_file_global_meta

!> \brief Checks if file_fname is allocated in the yaml object
!! \return true if file_fname is allocated
pure function has_file_fname (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_fname()
end function has_file_fname

!> \brief Checks if file_frequnit is allocated in the yaml object
!! \return true if file_frequnit is allocated
pure function has_file_frequnit (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_frequnit()
end function has_file_frequnit

!> \brief Checks if file_freq is allocated in the yaml object
!! \return true if file_freq is allocated
pure function has_file_freq (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_freq()
end function has_file_freq

!> \brief Checks if file_timeunit is allocated in the yaml object
!! \return true if file_timeunit is allocated
pure function has_file_timeunit (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_timeunit()
end function has_file_timeunit

!> \brief Checks if file_unlimdim is allocated in the yaml object
!! \return true if file_unlimdim is allocated
pure function has_file_unlimdim (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_unlimdim()
end function has_file_unlimdim

!> \brief Checks if file_sub_region is allocated in the yaml object
!! \return true if file_sub_region is allocated
pure function has_file_sub_region (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_sub_region()
end function has_file_sub_region

!> \brief Checks if file_new_file_freq is allocated in the yaml object
!! \return true if file_new_file_freq is allocated
pure function has_file_new_file_freq (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_new_file_freq()
end function has_file_new_file_freq

!> \brief Checks if file_new_file_freq_units is allocated in the yaml object
!! \return true if file_new_file_freq_units is allocated
pure function has_file_new_file_freq_units (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_new_file_freq_units()
end function has_file_new_file_freq_units

!> \brief Checks if file_start_time is allocated in the yaml object
!! \return true if file_start_time is allocated
pure function has_file_start_time (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_start_time()
end function has_file_start_time

!> \brief Checks if file_duration is allocated in the yaml object
!! \return true if file_duration is allocated
pure function has_file_duration (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_duration()
end function has_file_duration

!> \brief Checks if file_duration_units is allocated in the yaml object
!! \return true if file_duration_units is allocated
pure function has_file_duration_units (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_duration_units()
end function has_file_duration_units

!> \brief Checks if file_varlist is allocated in the yaml object
!! \return true if file_varlist is allocated
pure function has_file_varlist (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_varlist()
end function has_file_varlist

!> \brief Checks if file_global_meta is allocated in the yaml object
!! \return true if file_global_meta is allocated
pure function has_file_global_meta (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%diag_yaml_file%has_file_global_meta()
end function has_file_global_meta

!> @brief Sets the domain and type of domain from the axis IDs
subroutine set_domain_from_axis(this, diag_axis, axes)
  class(fmsDiagFile_type), intent(inout)       :: this          !< The file object
  class(fmsDiagAxisContainer_type), intent(in) :: diag_axis(:) !< Array of diag_axis
  integer, intent(in) :: axes (:)

  call get_domain_and_domain_type(diag_axis, axes, this%type_of_domain, this%domain, this%get_file_fname())
end subroutine set_domain_from_axis

!> @brief Set the domain and the type_of_domain for a file
!> @details This subroutine is going to be called once by every variable in the file
!! in register_diag_field. It will update the domain and the type_of_domain if needed and verify that
!! all the variables are in the same domain
subroutine set_file_domain(this, domain, type_of_domain)
  class(fmsDiagFile_type), intent(inout)       :: this            !< The file object
  integer,                 INTENT(in)          :: type_of_domain !< fileobj_type to use
  CLASS(diagDomain_t),     INTENT(in), target  :: domain         !< Domain

  if (type_of_domain .ne. this%type_of_domain) then
  !! If the current type_of_domain in the file obj is not the same as the variable calling this subroutine

    if (type_of_domain .eq. NO_DOMAIN .or. this%type_of_domain .eq. NO_DOMAIN) then
    !! If they are not the same then one of them can be NO_DOMAIN
    !! (i.e a file can have variables that are not domain decomposed and variables that are)

      if (type_of_domain .ne. NO_DOMAIN) then
      !! Update the file's type_of_domain and domain if needed
        this%type_of_domain = type_of_domain
        this%domain => domain
      endif

    else
    !! If they are not the same and of them is not NO_DOMAIN, then crash because the variables don't have the
    !! same domain (i.e a file has a variable is that in a 2D domain and one that is in a UG domain)

      call mpp_error(FATAL, "The file: "//this%get_file_fname()//" has variables that are not in the same domain")
    endif
  endif

end subroutine set_file_domain

!> @brief Loops through a variable's axis_ids and adds them to the FMSDiagFile object if they don't exist
subroutine add_axes(this, axis_ids, diag_axis, naxis)
  class(fmsDiagFile_type),          intent(inout)       :: this          !< The file object
  integer,                          INTENT(in)          :: axis_ids(:)   !< Array of axes_ids
  class(fmsDiagAxisContainer_type), intent(inout)       :: diag_axis(:)  !< Diag_axis object
  integer,                          intent(inout)       :: naxis         !< Number of axis that have been registered

  integer :: i, j !< For do loops
  logical :: is_cube_sphere !< Flag indicating if the file's domain is a cubesphere

  is_cube_sphere = .false.

  select type(this)
  type is (subRegionalFile_type)
    if (.not. this%is_subaxis_defined) then
      if (associated(this%domain)) then
        if (this%domain%get_ntiles() .eq. 6) is_cube_sphere = .true.
      endif

      call define_subaxis(diag_axis, axis_ids, naxis, this%get_file_sub_region(), &
        is_cube_sphere, this%write_on_this_pe)
      this%is_subaxis_defined = .true.

      !> add the axis to the list of axis in the file
      if (this%write_on_this_pe) then
        do i = 1, size(axis_ids)
          this%number_of_axis = this%number_of_axis + 1 !< This is the current number of axis in the file
          this%axis_ids(this%number_of_axis) = diag_axis(axis_ids(i))%axis%get_subaxes_id()
        enddo
      else
        this%axis_ids = diag_null
      endif
    endif
    return
  type is (fmsDiagFile_type)
    do i = 1, size(axis_ids)
      do j = 1, this%number_of_axis
        !> Check if the axis already exists, return
        if (axis_ids(i) .eq. this%axis_ids(j)) return
      enddo

      !> If the axis does not exist add it to the list
      this%number_of_axis = this%number_of_axis + 1
      this%axis_ids(this%number_of_axis) = axis_ids(i)
    enddo
  end select
end subroutine add_axes

!> @brief adds the start time to the fileobj
!! @note This should be called from the register field calls. It can be called multiple times (one for each variable)
!! So it needs to make sure that the start_time is the same for each variable. The initial value is the base_time
subroutine add_start_time(this, start_time)
  class(fmsDiagFile_type), intent(inout)       :: this            !< The file object
  TYPE(time_type),         intent(in)          :: start_time     !< Start time to add to the fileobj

  !< If the start_time sent in is equal to the base_time return because
  !! this%start_time was already set to the base_time
  if (start_time .eq. get_base_time()) return

  if (this%start_time .ne. get_base_time()) then
    !> If the this%start_time is not equal to the base_time from the diag_table
    !! this%start_time was already updated so make sure it is the same or error out
    if (this%start_time .ne. start_time)&
      call mpp_error(FATAL, "The variables associated with the file:"//this%get_file_fname()//" have"&
      &" different start_time")
  else
    !> If the this%start_time is equal to the base_time,
    !! simply update it with the start_time and set up the *_output variables
    this%start_time = start_time
    this%last_output = start_time
    this%next_output = diag_time_inc(start_time, this%get_file_freq(), this%get_file_frequnit())
    this%next_next_output = diag_time_inc(this%next_output, this%get_file_freq(), this%get_file_frequnit())
  endif

end subroutine

!> writes out internal values for fmsDiagFile_type object
subroutine dump_file_obj(this, unit_num)
  class(fmsDiagFile_type), intent(in) :: this !< the file object
  integer, intent(in) :: unit_num !< passed in from dump_diag_obj
                                  !! will either be for new log file or stdout
  write( unit_num, *) 'file id:', this%id
  write( unit_num, *) 'start time:', date_to_string(this%start_time)
  write( unit_num, *) 'last_output', date_to_string(this%last_output)
  write( unit_num, *) 'next_output', date_to_string(this%next_output)
  write( unit_num, *)'next_next_output', date_to_string(this%next_next_output)
  write( unit_num, *)'next_open', date_to_string(this%next_open)

  if( allocated(this%fileobj)) write( unit_num, *)'fileobj path', this%fileobj%path

  write( unit_num, *)'type_of_domain', this%type_of_domain
  if( allocated(this%file_metadata_from_model)) write( unit_num, *) 'file_metadata_from_model', &
                                                                    this%file_metadata_from_model
  if( allocated(this%field_ids)) write( unit_num, *)'field_ids', this%field_ids
  if( allocated(this%field_registered)) write( unit_num, *)'field_registered', this%field_registered
  if( allocated(this%num_registered_fields)) write( unit_num, *)'num_registered_fields', this%num_registered_fields
  if( allocated(this%axis_ids)) write( unit_num, *)'axis_ids', this%axis_ids(1:this%number_of_axis)

end subroutine

!< @brief Opens the diag_file if it is time to do so
subroutine open_diag_file(this, time_step, file_is_opened)
  class(fmsDiagFileContainer_type), intent(inout), target :: this            !< The file object
  TYPE(time_type),                  intent(in)            :: time_step       !< Current model step time
  logical,                          intent(out)           :: file_is_opened  !< .true. if the file was opened in this
                                                                             !! time

  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open
  class(diagDomain_t),     pointer     :: domain         !< The domain used in the file
  character(len=:),        allocatable :: diag_file_name !< The file name as defined in the yaml
  character(len=128)                   :: base_name      !< The file name as defined in the yaml
                                                         !! without the wildcard definition
  character(len=128)                   :: file_name      !< The file name as it will be written to disk
  character(len=128)                   :: temp_name      !< Temp variable to store the file_name
  character(len=128)                   :: start_date     !< The start_time as a string that will be added to
                                                         !! the begining of the filename (start_date.filename)
  character(len=128)                   :: suffix         !< The current time as a string that will be added to
                                                         !! the end of filename
  integer                              :: pos            !< Index of the filename with the first "%" in the file name
  INTEGER                              :: year           !< The year of the start_date
  INTEGER                              :: month          !< The month of the start_date
  INTEGER                              :: day            !< The day of the start_date
  INTEGER                              :: hour           !< The hour of the start_date
  INTEGER                              :: minute         !< The minute of the start_date
  INTEGER                              :: second         !< The second of the start_date
  character(len=4)                     :: mype_string    !< The pe as a string
  logical                              :: is_regional    !< Flag indicating if the file is regional
  integer, allocatable                 :: pes(:)         !< Array of the pes in the current pelist

  diag_file => this%FMS_diag_file
  domain => diag_file%domain

  file_is_opened = .false.
  !< Go away if it is not time to open the file
  if (diag_file%next_open > time_step) return

  is_regional = .false.
  !< Figure out what fileobj to use!
  if (.not. allocated(diag_file%fileobj)) then
    select type (diag_file)
    type is (subRegionalFile_type)
      !< In this case each PE is going to write its own file
      allocate(FmsNetcdfFile_t :: diag_file%fileobj)
      is_regional = .true.
    type is (fmsDiagFile_type)
      !< Use the type_of_domain to get the correct fileobj
      select case (diag_file%type_of_domain)
      case (NO_DOMAIN)
        allocate(FmsNetcdfFile_t :: diag_file%fileobj)
      case (TWO_D_DOMAIN)
        allocate(FmsNetcdfDomainFile_t :: diag_file%fileobj)
      case (UG_DOMAIN)
        allocate(FmsNetcdfUnstructuredDomainFile_t :: diag_file%fileobj)
      end select
    end select
  else
    !< In this case, we are opening a new file so close the current the file
    call this%close_diag_file()
  endif

  !< Figure out what to name of the file
  diag_file_name = diag_file%get_file_fname()

  !< If using the new_file_freq figure out what the name is based on the current time
  if (diag_file%has_file_new_file_freq()) then
    !< If using a wildcard file name (i.e ocn%4yr%2mo%2dy%2hr), get the basename (i.e ocn)
    pos = INDEX(diag_file_name, '%')
    if (pos > 0) base_name = diag_file_name(1:pos-1)
    suffix = get_time_string(diag_file_name, time_step) !TODO fname_time?
    base_name = trim(base_name)//trim(suffix)
  else
    base_name = trim(diag_file_name)
  endif

  !< Add the ens number to the file name (if it exists)
  file_name = trim(base_name)
  call get_instance_filename(base_name, file_name)

  !< Prepend the file start_time to the file name if prepend_date == .TRUE. in
  !! the namelist
  IF ( prepend_date ) THEN
    call get_date(diag_file%start_time, year, month, day, hour, minute, second)
    write (start_date, '(1I20.4, 2I2.2)') year, month, day

    file_name = TRIM(adjustl(start_date))//'.'//TRIM(file_name)
  END IF

  file_name = trim(file_name)//".nc"

  !< If this is a regional file add the PE and the tile_number to the filename
  if (is_regional) then
    !< Get the pe number that will be appended to the end of the file
    write(mype_string,'(I0.4)') mpp_pe()

    !< Add the tile number if appropriate
    select type (domain)
    type is (DIAGDOMAIN2D_T)
      temp_name = file_name
      call get_mosaic_tile_file(temp_name, file_name, .true., domain%domain2)
    end select

    file_name = trim(file_name)//"."//trim(mype_string)
  endif

  !< Open the file!
  select type (fileobj => diag_file%fileobj)
  type is (FmsNetcdfFile_t)
    if (is_regional) then
      if (.not. open_file(fileobj, file_name, "overwrite", pelist=(/mpp_pe()/))) &
      &call mpp_error(FATAL, "Error opening the file:"//file_name)
   else
      allocate(pes(mpp_npes()))
      call mpp_get_current_pelist(pes)

      if (.not. open_file(fileobj, file_name, "overwrite", pelist=pes)) &
      &call mpp_error(FATAL, "Error opening the file:"//file_name)
   endif
  type is (FmsNetcdfDomainFile_t)
    select type (domain)
    type is (diagDomain2d_t)
      if (.not. open_file(fileobj, file_name, "overwrite", domain%Domain2)) &
        &call mpp_error(FATAL, "Error opening the file:"//file_name)
    end select
  type is (FmsNetcdfUnstructuredDomainFile_t)
    select type (domain)
    type is (diagDomainUg_t)
      if (.not. open_file(fileobj, file_name, "overwrite", domain%DomainUG)) &
        &call mpp_error(FATAL, "Error opening the file:"//file_name)
    end select
  end select

  if (diag_file%has_file_new_file_freq()) then
    diag_file%next_open = diag_time_inc(diag_file%next_open, diag_file%get_file_new_file_freq(), &
                                        diag_file%get_file_new_file_freq_units())
  else
    diag_file%next_open = diag_time_inc(diag_file%next_open, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
  endif

  file_is_opened = .true.
  domain => null()
  diag_file => null()
end subroutine open_diag_file

!> \brief Write the time metadata to the diag file
subroutine write_time_metadata(this)
  class(fmsDiagFileContainer_type), intent(inout), target :: this !< The file object

  class(fmsDiagFile_type), pointer  :: diag_file      !< Diag_file object to open
  class(FmsNetcdfFile_t),  pointer  :: fileobj        !< The fileobj to write to
  character(len=50)                 :: time_units_str !< Time units written as a string
  character(len=50)                 :: calendar       !< The calendar name

  character(len=:), allocatable :: time_var_name !< The name of the time variable as it is defined in the yaml

  diag_file => this%FMS_diag_file
  fileobj => diag_file%fileobj

  time_var_name = diag_file%get_file_unlimdim()
  call register_axis(fileobj, time_var_name, unlimited)

  WRITE(time_units_str, 11)  &
    TRIM(time_unit_list(diag_file%get_file_timeunit())), get_base_year(),&
         & get_base_month(), get_base_day(), get_base_hour(), get_base_minute(), get_base_second()
11  FORMAT(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)

  !TODO harcodded "double"
  call register_field(fileobj, time_var_name, "double", (/time_var_name/))
  call register_variable_attribute(fileobj, time_var_name, "units", trim(time_units_str), &
    str_len=len_trim(time_units_str))

  call register_variable_attribute(fileobj, time_var_name, "axis", "T", str_len=1 )
  call register_variable_attribute(fileobj, time_var_name, "long_name", trim(time_var_name), &
    str_len=len_trim(time_var_name) )

  !TODO no need to have both attributes, probably?
  calendar = valid_calendar_types(get_calendar_type())
  call register_variable_attribute(fileobj, time_var_name, "calendar_type", &
    uppercase(trim(calendar)), str_len=len_trim(calendar))
  call register_variable_attribute(fileobj, time_var_name, "calendar", &
    lowercase(trim(calendar)), str_len=len_trim(calendar))

  if (diag_file%time_ops) call register_variable_attribute(fileobj, time_var_name, "bounds", &
    trim(time_var_name)//"_bounds", str_len=len_trim(time_var_name//"_bounds"))

end subroutine write_time_metadata

!> \brief Determine if it is time to "write" to the file
logical function is_time_to_write(this, time_step)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object
  TYPE(time_type),                  intent(in)           :: time_step       !< Current model step time

  if (time_step >= this%FMS_diag_file%next_output) then
    is_time_to_write = .true.
    if (this%FMS_diag_file%is_static) return
    if (time_step >= this%FMS_diag_file%next_next_output) &
      call mpp_error(FATAL, this%FMS_diag_file%get_file_fname()//": You skip a time!")
  else
    is_time_to_write = .false.
  endif
end function is_time_to_write

!> \brief Determine if the current PE has data to write
logical function writing_on_this_pe(this)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object

  writing_on_this_pe = .true.

  select type(diag_file => this%FMS_diag_file)
  type is (subRegionalFile_type)
    writing_on_this_pe = diag_file%write_on_this_pe
  end select

end function

!> \brief Write out the time data to the file
subroutine write_time_data(this, time_step)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object
  TYPE(time_type),                  intent(in)           :: time_step       !< Current model step time

  real                                 :: dif            !< The time as a real number
  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open
  class(FmsNetcdfFile_t),  pointer     :: fileobj        !< The fileobj to write to

  diag_file => this%FMS_diag_file
  fileobj => diag_file%fileobj

  !> dif is the time as a real that is evaluated
  dif = get_date_dif(time_step, get_base_time(), diag_file%get_file_timeunit())
  select type (fileobj)
  type is (FmsNetcdfDomainFile_t)
    call write_data(fileobj, diag_file%get_file_unlimdim(), dif, &
      unlim_dim_level=diag_file%unlimited_dimension)
  type is (FmsNetcdfUnstructuredDomainFile_t)
    call write_data(fileobj, diag_file%get_file_unlimdim(), dif, &
      unlim_dim_level=diag_file%unlimited_dimension)
  type is (FmsNetcdfFile_t)
    call write_data(fileobj, diag_file%get_file_unlimdim(), dif, &
      unlim_dim_level=diag_file%unlimited_dimension)
  end select

end subroutine write_time_data

!> \brief Set up the next_output and next_next_output variable in a file obj
subroutine update_next_write(this, time_step)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object
  TYPE(time_type),                  intent(in)           :: time_step       !< Current model step time

  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open

  diag_file => this%FMS_diag_file
  if (diag_file%is_static) then
    diag_file%next_output = diag_time_inc(diag_file%next_output, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
    diag_file%next_next_output = diag_time_inc(diag_file%next_output, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
  else
    diag_file%next_output = diag_time_inc(diag_file%next_output, diag_file%get_file_freq(), &
      diag_file%get_file_frequnit())
    diag_file%next_next_output = diag_time_inc(diag_file%next_output, diag_file%get_file_freq(), &
      diag_file%get_file_frequnit())
  endif

end subroutine update_next_write

!> \brief Increase the unlimited dimension variable that the file is currently being written to
subroutine increase_unlimited_dimension(this)
  class(fmsDiagFileContainer_type), intent(inout), target   :: this            !< The file object

  this%FMS_diag_file%unlimited_dimension = this%FMS_diag_file%unlimited_dimension + 1
end subroutine increase_unlimited_dimension

!< @brief Writes the axis metadata for the file
subroutine write_axis_metadata(this, diag_axis)
  class(fmsDiagFileContainer_type), intent(inout), target :: this            !< The file object
  class(fmsDiagAxisContainer_type), intent(in)            :: diag_axis(:)    !< Diag_axis object

  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open
  class(FmsNetcdfFile_t),  pointer     :: fileobj        !< The fileobj to write to
  integer                              :: i              !< For do loops
  integer                              :: j              !< diag_file%axis_ids(i) (for less typing)
  integer                              :: parent_axis_id !< Id of the parent_axis

  diag_file => this%FMS_diag_file
  fileobj => diag_file%fileobj

  do i = 1, diag_file%number_of_axis
    j = diag_file%axis_ids(i)
    parent_axis_id = diag_axis(j)%axis%get_parent_axis_id()
    if (parent_axis_id .eq. DIAG_NULL) then
      call diag_axis(j)%axis%write_axis_metadata(fileobj)
    else
      call diag_axis(j)%axis%write_axis_metadata(fileobj, diag_axis(parent_axis_id)%axis)
    endif
  enddo

end subroutine write_axis_metadata

!< @brief Writes the axis data for the file
subroutine write_axis_data(this, diag_axis)
  class(fmsDiagFileContainer_type), intent(inout), target :: this            !< The file object
  class(fmsDiagAxisContainer_type), intent(in)            :: diag_axis(:)    !< Diag_axis object

  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open
  class(FmsNetcdfFile_t),  pointer     :: fileobj        !< The fileobj to write to
  integer                              :: i              !< For do loops
  integer                              :: j              !< diag_file%axis_ids(i) (for less typing)
  integer                              :: parent_axis_id !< Id of the parent_axis

  diag_file => this%FMS_diag_file
  fileobj => diag_file%fileobj

  do i = 1, diag_file%number_of_axis
    j = diag_file%axis_ids(i)
    parent_axis_id = diag_axis(j)%axis%get_parent_axis_id()
    if (parent_axis_id .eq. DIAG_NULL) then
      call diag_axis(j)%axis%write_axis_data(fileobj)
    else
      call diag_axis(j)%axis%write_axis_data(fileobj, diag_axis(parent_axis_id)%axis)
    endif
  enddo

end subroutine write_axis_data

!< @brief Closes the diag_file
subroutine close_diag_file(this)
  class(fmsDiagFileContainer_type), intent(inout), target :: this            !< The file object

  !< The select types are needed here because otherwise the code will go to the
  !! wrong close_file routine and things will not close propertly
  select type( fileobj => this%FMS_diag_file%fileobj)
  type is (FmsNetcdfDomainFile_t)
    call close_file(fileobj)
  type is (FmsNetcdfFile_t)
    call close_file(fileobj)
  type is (FmsNetcdfUnstructuredDomainFile_t)
    call close_file(fileobj)
  end select

end subroutine close_diag_file

#endif
end module fms_diag_file_object_mod
