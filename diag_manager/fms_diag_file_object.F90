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
                       register_axis, register_field, register_variable_attribute, write_data, &
                       dimension_exists, register_global_attribute
use diag_data_mod, only: DIAG_NULL, NO_DOMAIN, max_axes, SUB_REGIONAL, get_base_time, DIAG_NOT_REGISTERED, &
                         TWO_D_DOMAIN, UG_DOMAIN, prepend_date, DIAG_DAYS, VERY_LARGE_FILE_FREQ, &
                         get_base_year, get_base_month, get_base_day, get_base_hour, get_base_minute, &
                         get_base_second, time_unit_list, time_average, time_rms, time_max, time_min, time_sum, &
                         time_diurnal, time_power, time_none, avg_name, no_units, pack_size_str, &
                         middle_time, begin_time, end_time, MAX_STR_LEN, index_gridtype, latlon_gridtype
use time_manager_mod, only: time_type, operator(>), operator(/=), operator(==), get_date, get_calendar_type, &
                            VALID_CALENDAR_TYPES, operator(>=), date_to_string, &
                            OPERATOR(/), OPERATOR(+), operator(<)
use fms_diag_time_utils_mod, only: diag_time_inc, get_time_string, get_date_dif
use fms_diag_yaml_mod, only: diag_yaml, diagYamlObject_type, diagYamlFiles_type, subRegion_type, diagYamlFilesVar_type
use fms_diag_axis_object_mod, only: diagDomain_t, get_domain_and_domain_type, fmsDiagAxis_type, &
                                    fmsDiagAxisContainer_type, DIAGDOMAIN2D_T, DIAGDOMAINUG_T, &
                                    fmsDiagFullAxis_type, define_diurnal_axis, &
                                    fmsDiagDiurnalAxis_type, create_new_z_subaxis, is_parent_axis, &
                                    define_new_subaxis_latlon, define_new_subaxis_index, fmsDiagSubAxis_type
use fms_diag_field_object_mod, only: fmsDiagField_type
use fms_diag_output_buffer_mod, only: fmsDiagOutputBuffer_type
use mpp_mod, only: mpp_get_current_pelist, mpp_npes, mpp_root_pe, mpp_pe, mpp_error, FATAL, stdout, &
                   uppercase, lowercase, NOTE

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
  TYPE(time_type) :: no_more_data     !< Time to stop receiving data for this file
  logical         :: done_writing_data!< .True. if finished writing data

  !< This will be used when using the new_file_freq keys in the diag_table.yaml
  TYPE(time_type) :: next_close       !< Time to close the file
  logical         :: is_file_open     !< .True. if the file is opened

  class(FmsNetcdfFile_t), allocatable :: fms2io_fileobj !< fms2_io file object for this history file
  type(diagYamlFiles_type), pointer :: diag_yaml_file => null() !< Pointer to the diag_yaml_file data
  integer                                      :: type_of_domain !< The type of domain to use to open the file
                                                                 !! NO_DOMAIN, TWO_D_DOMAIN, UG_DOMAIN, SUB_REGIONAL
  class(diagDomain_t), pointer                 :: domain         !< The domain to use,
                                                                 !! null if NO_DOMAIN or SUB_REGIONAL
  character(len=:) , dimension(:), allocatable :: file_metadata_from_model !< File metadata that comes from
                                                                           !! the model.
  integer, dimension(:), allocatable :: field_ids !< Variable IDs corresponding to file_varlist
  integer, dimension(:), allocatable :: yaml_ids !< IDs corresponding to the yaml field section
  logical, dimension(:), private, allocatable :: field_registered   !< Array corresponding to `field_ids`, .true.
                                                                 !! if the variable has been registered and
                                                                 !! `field_id` has been set for the variable
  integer, allocatable                         :: num_registered_fields !< The number of fields registered
                                                                        !! to the file
  integer, dimension(:), allocatable :: axis_ids !< Array of axis ids in the file
  integer :: number_of_axis !< Number of axis in the file
  integer, dimension(:), allocatable :: buffer_ids !< array of buffer ids associated with the file
  integer :: number_of_buffers !< Number of buffers that have been added to the file
  logical :: time_ops !< .True. if file contains variables that are time_min, time_max, time_average or time_sum
  integer :: unlim_dimension_level !< The unlimited dimension level currently being written
  logical :: data_has_been_written !< .True. if data has been written for the current unlimited dimension level
  logical :: is_static !< .True. if the frequency is -1
  integer :: nz_subaxis !< The number of Z axis currently added to the file

 contains
  procedure, public :: add_field_and_yaml_id
  procedure, public :: add_buffer_id
  procedure, public :: is_field_registered
  procedure, public :: init_diurnal_axis
  procedure, public :: has_file_metadata_from_model
  procedure, public :: has_fileobj
  procedure, public :: has_diag_yaml_file
  procedure, public :: set_domain_from_axis
  procedure, public :: set_file_domain
  procedure, public :: add_axes
  procedure, public :: add_new_axis
  procedure, public :: update_write_on_this_pe
  procedure, public :: get_write_on_this_pe
  procedure, public :: does_axis_exist
  procedure, public :: define_new_subaxis
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
 procedure, public :: get_file_sub_region_grid_type
 procedure, public :: get_file_new_file_freq
 procedure, public :: get_filename_time
 procedure, public :: get_file_new_file_freq_units
 procedure, public :: get_file_start_time
 procedure, public :: get_file_duration
 procedure, public :: get_file_duration_units
 procedure, public :: get_file_varlist
 procedure, public :: get_file_global_meta
 procedure, public :: is_done_writing_data
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
 procedure, public :: get_buffer_ids
 procedure, public :: get_number_of_buffers
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
  procedure :: is_regional
  procedure :: is_file_static
  procedure :: open_diag_file
  procedure :: write_global_metadata
  procedure :: write_time_metadata
  procedure :: write_field_data
  procedure :: write_axis_metadata
  procedure :: write_field_metadata
  procedure :: write_axis_data
  procedure :: writing_on_this_pe
  procedure :: is_time_to_write
  procedure :: is_time_to_close_file
  procedure :: write_time_data
  procedure :: update_next_write
  procedure :: update_current_new_file_freq_index
  procedure :: increase_unlim_dimension_level
  procedure :: get_unlim_dimension_level
  procedure :: get_next_output
  procedure :: get_next_next_output
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
           obj%write_on_this_pe = .true.
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
     allocate(obj%buffer_ids(diag_yaml%diag_files(i)%size_file_varlist()))
     allocate(obj%yaml_ids(diag_yaml%diag_files(i)%size_file_varlist()))
     allocate(obj%field_registered(diag_yaml%diag_files(i)%size_file_varlist()))
     !! Initialize the integer arrays
     obj%field_ids = DIAG_NOT_REGISTERED
     obj%yaml_ids = DIAG_NOT_REGISTERED
     obj%buffer_ids = DIAG_NOT_REGISTERED
     obj%field_registered = .FALSE.
     obj%num_registered_fields = 0
     obj%number_of_buffers = 0

     !> These will be set in a set_file_domain
     obj%type_of_domain = NO_DOMAIN
     obj%domain => null()

     !> This will be set in a add_axes
     allocate(obj%axis_ids(max_axes))
     obj%number_of_axis = 0

     !> Set the start_time of the file to the base_time and set up the *_output variables
     obj%done_writing_data = .false.
     obj%start_time = get_base_time()
     obj%last_output = get_base_time()
     obj%next_output = diag_time_inc(obj%start_time, obj%get_file_freq(), obj%get_file_frequnit())
     obj%next_next_output = diag_time_inc(obj%next_output, obj%get_file_freq(), obj%get_file_frequnit())

     if (obj%has_file_new_file_freq()) then
       obj%next_close = diag_time_inc(obj%start_time, obj%get_file_new_file_freq(), &
                                        obj%get_file_new_file_freq_units())
     else
       obj%next_close = diag_time_inc(obj%start_time, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
     endif
     obj%is_file_open = .false.

     if(obj%has_file_duration()) then
       obj%no_more_data = diag_time_inc(obj%start_time, obj%get_file_duration(), &
                                          obj%get_file_duration_units())
     else
       obj%no_more_data = diag_time_inc(obj%start_time, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
     endif

     obj%time_ops = .false.
     obj%unlim_dimension_level = 0
     obj%is_static = obj%get_file_freq() .eq. -1
     obj%nz_subaxis = 0

     nullify(obj)
   enddo set_ids_loop
   fms_diag_files_object_init = .true.
  else
   fms_diag_files_object_init = .false.
!  mpp_error("fms_diag_files_object_init: The diag_table.yaml file has not been correctly parsed.",&
!    FATAL)
  endif
end function fms_diag_files_object_init

!< @brief Determine if the field corresponding to the field_id was registered to the file
!! @return .True. if the field was registed to the file
pure logical function is_field_registered(this, field_id)
  class(fmsDiagFile_type), intent(in)    :: this         !< The file object
  integer,                 intent(in)    :: field_id     !< Id of the field to check

  is_field_registered = this%field_registered(field_id)
end function is_field_registered

!> \brief Adds a field and yaml ID to the file
subroutine add_field_and_yaml_id (this, new_field_id, yaml_id)
  class(fmsDiagFile_type), intent(inout) :: this         !< The file object
  integer,                 intent(in)    :: new_field_id !< The field ID to be added to field_ids
  integer,                 intent(in)    :: yaml_id      !< The yaml_id

  this%num_registered_fields = this%num_registered_fields + 1
  if (this%num_registered_fields .le. size(this%field_ids)) then
    this%field_ids( this%num_registered_fields ) = new_field_id
    this%yaml_ids( this%num_registered_fields ) = yaml_id
    this%field_registered( this%num_registered_fields ) = .true.
  else
    call mpp_error(FATAL, "The file: "//this%get_file_fname()//" has already been assigned its maximum "//&
                 "number of fields.")
  endif
end subroutine add_field_and_yaml_id

!> \brief Adds a buffer_id to the file object
subroutine add_buffer_id (this, buffer_id)
  class(fmsDiagFile_type), intent(inout) :: this         !< The file object
  integer,                 intent(in)    :: buffer_id    !< Buffer id to add to the file

  this%number_of_buffers = this%number_of_buffers + 1
  this%buffer_ids(this%number_of_buffers) = buffer_id

end subroutine add_buffer_id

!> \brief Initializes a diurnal axis for a fileobj
!! \note This is going to be called for every variable in the file, if the variable is not a diurnal variable
!! it will do nothing. It only defined a diurnal axis once.
subroutine init_diurnal_axis(this, diag_axis, naxis, yaml_id)
  class(fmsDiagFile_type),          intent(inout) :: this         !< The file object
  class(fmsDiagAxisContainer_type), intent(inout) :: diag_axis(:) !< Array of diag_axis object
  integer,                          intent(inout) :: naxis        !< Number of diag_axis that heve been defined
  integer,                          intent(in)    :: yaml_id      !< The ID to the variable's yaml

  integer                              :: i           !< For do loops
  type(diagYamlFilesVar_type), pointer :: field_yaml  !< pointer to the yaml entry

  field_yaml => diag_yaml%get_diag_field_from_id(yaml_id)

  !< Go away if the file does not need a diurnal axis
  if (.not. field_yaml%has_n_diurnal()) return

  !< Check if the diurnal axis is already defined for this number of diurnal samples
  do i = 1, this%number_of_axis
    select type(axis=>diag_axis(this%axis_ids(i))%axis)
    type is (fmsDiagDiurnalAxis_type)
      if(field_yaml%get_n_diurnal() .eq. axis%get_diurnal_axis_samples()) return
    end select
  end do

  !< If it is not already defined, define it
  call define_diurnal_axis(diag_axis, naxis, field_yaml%get_n_diurnal(), .true.)
  call define_diurnal_axis(diag_axis, naxis, field_yaml%get_n_diurnal(), .False.)

  !< Add it to the list of axis for the file
  this%number_of_axis = this%number_of_axis + 1
  this%axis_ids(this%number_of_axis) = naxis !< This is the diurnal axis edges

  this%number_of_axis = this%number_of_axis + 1
  this%axis_ids(this%number_of_axis) = naxis - 1 !< This the diurnal axis

end subroutine init_diurnal_axis

!> \brief Set the time_ops variable in the diag_file object
subroutine set_file_time_ops(this, VarYaml, is_static)
  class(fmsDiagFile_type),      intent(inout) :: this      !< The file object
  type (diagYamlFilesVar_type), intent(in)    :: VarYaml   !< The variable's yaml file
  logical,                      intent(in)    :: is_static !< Flag indicating if variable is static

  !< Go away if the file is static
  if (this%is_static) return

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
  has_fileobj = allocated(this%fms2io_fileobj)
end function has_fileobj

!> \brief Logical function to determine if the variable diag_yaml_file has been allocated or associated
!! \return .True. if diag_yaml_file exists .False. if diag_yaml has not been set
pure logical function has_diag_yaml_file (this)
  class(fmsDiagFile_type), intent(in) :: this !< The file object
  has_diag_yaml_file = associated(this%diag_yaml_file)
end function has_diag_yaml_file

!> \brief Get the time to use to determine the filename, if using a wildcard file name (i.e ocn%4yr%2mo%2dy%2hr)
!! \return The time to use when determining the filename
function get_filename_time(this) &
  result(res)
    class(fmsDiagFile_type), intent(in) :: this !< The file object
    type(time_type) :: res

    select case (this%diag_yaml_file%get_filename_time())
    case (begin_time)
      res = this%last_output
    case (middle_time)
      res = (this%last_output + this%next_close)/2
    case (end_time)
      res = this%next_close
    end select
end function get_filename_time

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

!< @brief Query for the subregion grid type (latlon or index)
!! @return subregion grid type
function get_file_sub_region_grid_type(this) &
  result(res)
  class(fmsDiagFile_type), intent(in) :: this !< Diag file object
  integer :: res

  type(subRegion_type) :: subregion !< Subregion type

  subregion = this%diag_yaml_file%get_file_sub_region()
  res = subregion%grid_type
end function get_file_sub_region_grid_type

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
 character (len=MAX_STR_LEN), allocatable, dimension(:,:) :: res
  res = this%diag_yaml_file%get_file_global_meta()
end function get_file_global_meta

!> \brief Determines if done writing data
!! \return .True. if done writing data
pure function is_done_writing_data (this) result(res)
 class(fmsDiagFile_type), intent(in) :: this !< The file object
 logical :: res
  res = this%done_writing_data
end function is_done_writing_data

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
subroutine add_axes(this, axis_ids, diag_axis, naxis, yaml_id, buffer_id, output_buffers)
  class(fmsDiagFile_type),                 intent(inout) :: this              !< The file object
  integer,                                 INTENT(in)    :: axis_ids(:)       !< Array of axes_ids
  class(fmsDiagAxisContainer_type),        intent(inout) :: diag_axis(:)      !< Diag_axis object
  integer,                                 intent(inout) :: naxis             !< Number of axis that have been
                                                                              !! registered
  integer,                                 intent(in)    :: yaml_id           !< Yaml id of the field section for
                                                                              !! this var
  integer,                                 intent(in)    :: buffer_id        !< ID of the buffer
  type(fmsDiagOutputBuffer_type),          intent(inout) :: output_buffers(:) !< Array of output buffers

  type(diagYamlFilesVar_type), pointer     :: field_yaml  !< pointer to the yaml entry

  integer              :: i, j             !< For do loops
  logical              :: is_cube_sphere   !< Flag indicating if the file's domain is a cubesphere
  logical              :: axis_found       !< Flag indicating that the axis was already to the file obj
  integer, allocatable :: var_axis_ids(:)  !< Array of the variable's axis ids
  integer              :: x_y_axis_id(2)   !< Ids of the x and y axis
  integer              :: x_or_y           !< integer indicating if the axis is x or y
  logical              :: is_x_or_y        !< flag indicating if the axis is x or y
  integer              :: subregion_gridtype !< The type of the subregion (latlon or index)
  logical              :: write_on_this_pe !< Flag indicating if the current pe is in the subregion

  is_cube_sphere = .false.
  subregion_gridtype = this%get_file_sub_region_grid_type()

  field_yaml => diag_yaml%get_diag_field_from_id(yaml_id)

  !< Created a copy here, because if the variable has a z subaxis var_axis_ids will be modified in
  !! `create_new_z_subaxis` to contain the id of the new z subaxis instead of the parent axis,
  !! which will be added to the the list of axis in the file object (axis_ids is intent(in),
  !! which is why the copy was needed)
  var_axis_ids = axis_ids

  if (field_yaml%has_var_zbounds()) then
    call create_new_z_subaxis(field_yaml%get_var_zbounds(), var_axis_ids, diag_axis, naxis, &
                              this%axis_ids, this%number_of_axis, this%nz_subaxis)
  endif

  select type(this)
  type is (subRegionalFile_type)
    if (associated(this%domain)) then
      if (this%domain%get_ntiles() .eq. 6) is_cube_sphere = .true.
    endif
    if (.not. this%get_write_on_this_pe()) return
    subaxis_defined: if (this%is_subaxis_defined) then
      do i = 1, size(var_axis_ids)
        select type (parent_axis => diag_axis(var_axis_ids(i))%axis)
        type is (fmsDiagFullAxis_type)
          axis_found = .false.
          is_x_or_y = parent_axis%is_x_or_y_axis()
          do j = 1, this%number_of_axis
            if (is_x_or_y) then
              if(is_parent_axis(this%axis_ids(j), var_axis_ids(i), diag_axis)) then
                axis_found = .true.
                var_axis_ids(i) = this%axis_ids(j) !Set the var_axis_id to the sub axis_id
                cycle
              endif
            elseif (var_axis_ids(i) .eq. this%axis_ids(j)) then
              axis_found = .true.
            endif
          enddo

          if (.not. axis_found) then
            if (is_x_or_y) then
              if (subregion_gridtype .eq. latlon_gridtype .and. is_cube_sphere) &
                call mpp_error(FATAL, "If using the cube sphere and defining the subregion with latlon "//&
                "the variable need to have the same x and y axis. Please check the variables in the file "//&
                trim(this%get_file_fname())//" or use indices to define the subregion.")

              select case (subregion_gridtype)
              case (index_gridtype)
                call define_new_subaxis_index(parent_axis, this%get_file_sub_region(), diag_axis, naxis, &
                  i, write_on_this_pe)
              case (latlon_gridtype)
                call define_new_subaxis_latlon(diag_axis, var_axis_ids(i:i), naxis, this%get_file_sub_region(), &
                  .false., write_on_this_pe)
              end select
              call this%update_write_on_this_pe(write_on_this_pe)
              if (.not. this%get_write_on_this_pe()) cycle
              call this%add_new_axis(naxis)
              var_axis_ids(i) = naxis
            else
              call this%add_new_axis(var_axis_ids(i))
            endif
          endif
        type is (fmsDiagSubAxis_type)
          axis_found = this%does_axis_exist(var_axis_ids(i))
          if (.not. axis_found) call this%add_new_axis(var_axis_ids(i))
        end select
      enddo
    else
      x_y_axis_id = diag_null
      do i = 1, size(var_axis_ids)
        select type (parent_axis => diag_axis(var_axis_ids(i))%axis)
        type is (fmsDiagFullAxis_type)
          if (.not. parent_axis%is_x_or_y_axis(x_or_y)) then
            axis_found = this%does_axis_exist(var_axis_ids(i))
            if (.not. axis_found) call this%add_new_axis(var_axis_ids(i))
          else
            x_y_axis_id(x_or_y) = var_axis_ids(i)
          endif
        type is (fmsDiagSubAxis_type)
          axis_found = this%does_axis_exist(var_axis_ids(i))
          if (.not. axis_found) call this%add_new_axis(var_axis_ids(i))
        end select
      enddo

      call this%define_new_subaxis(var_axis_ids, x_y_axis_id, is_cube_sphere, diag_axis, naxis)
      this%is_subaxis_defined = .true.
    endif subaxis_defined
  type is (fmsDiagFile_type)
    do i = 1, size(var_axis_ids)
      axis_found = this%does_axis_exist(var_axis_ids(i))
      if (.not. axis_found) call this%add_new_axis(var_axis_ids(i))
    enddo
  end select
  !> Add the axis to the buffer object
  call output_buffers(buffer_id)%add_axis_ids(var_axis_ids)
end subroutine add_axes

!> @brief Adds a new axis the list of axis in the diag file object
subroutine add_new_axis(this, var_axis_id)
  class(fmsDiagFile_type),                 intent(inout) :: this        !< The file object
  integer,                                 intent(in)    :: var_axis_id !< Axis id of the variable

  this%number_of_axis = this%number_of_axis + 1
  this%axis_ids(this%number_of_axis) = var_axis_id
end subroutine add_new_axis

!> @brief This updates write on this pe
subroutine update_write_on_this_pe(this, write_on_this_pe)
  class(fmsDiagFile_type),                 intent(inout) :: this             !< The file object
  logical,                                 intent(in)    :: write_on_this_pe !< .True. if the current PE is in
                                                                             !! subregion

  select type (this)
  type is (subRegionalFile_type)
    if (this%write_on_this_pe) this%write_on_this_pe = write_on_this_pe
  end select
end subroutine update_write_on_this_pe

!> @brief Query for the write_on_this_pe member of the diag file object
!! @return the write_on_this_pe member of the diag file object
function get_write_on_this_pe(this) &
  result(rslt)
  class(fmsDiagFile_type),                 intent(inout) :: this             !< The file object
  logical :: rslt
  rslt = .true.
  select type (this)
  type is (subRegionalFile_type)
    rslt= this%write_on_this_pe
  end select
end function get_write_on_this_pe

!< @brief Determine if an axis is already in the list of axis for a diag file
!! @return .True. if the axis is already in the list of axis for a diag file
function does_axis_exist(this, var_axis_id) &
  result(rslt)
  class(fmsDiagFile_type), intent(inout) :: this         !< The file object
  integer,                 intent(in)    :: var_axis_id  !< Variable axis id to check

  logical :: rslt
  integer :: j !< For do loops

  rslt = .false.
  do j = 1, this%number_of_axis
    !> Check if the axis already exists, move on
    if (var_axis_id .eq. this%axis_ids(j)) then
      rslt = .true.
      return
    endif
  enddo
end function

!> @brief Define a new sub axis
subroutine define_new_subaxis(this, var_axis_ids, x_y_axis_id, is_cube_sphere, diag_axis, naxis)
  class(fmsDiagFile_type),                 intent(inout) :: this              !< The file object
  integer,                                 INTENT(inout) :: var_axis_ids(:)   !< Original variable axis ids
  integer,                                 INTENT(in)    :: x_y_axis_id(:)    !< The ids of the x and y axis
  logical,                                 intent(in)    :: is_cube_sphere    !< .True. if the axis is in the cubesphere
  integer,                                 intent(inout) :: naxis             !< Number of axis current registered
  class(fmsDiagAxisContainer_type),        intent(inout) :: diag_axis(:)      !< Diag_axis object

  logical :: write_on_this_pe !< .True. if the current PE is in the subregion
  integer :: i, j             !< For do loop

  select case (this%get_file_sub_region_grid_type())
  case(latlon_gridtype)
    call define_new_subaxis_latlon(diag_axis, x_y_axis_id, naxis, this%get_file_sub_region(), is_cube_sphere, &
      write_on_this_pe)
    call this%update_write_on_this_pe(write_on_this_pe)
    if (.not. this%get_write_on_this_pe()) return
    call this%add_new_axis(naxis)
    call this%add_new_axis(naxis-1)
    do j = 1, size(var_axis_ids)
      if (x_y_axis_id(1) .eq. var_axis_ids(j)) var_axis_ids(j) = naxis - 1
      if (x_y_axis_id(2) .eq. var_axis_ids(j)) var_axis_ids(j) = naxis
    enddo
  case (index_gridtype)
    do i = 1, size(x_y_axis_id)
      select type (parent_axis => diag_axis(x_y_axis_id(i))%axis)
      type is (fmsDiagFullAxis_type)
        call define_new_subaxis_index(parent_axis, this%get_file_sub_region(), diag_axis, naxis, i, &
          write_on_this_pe)
        call this%update_write_on_this_pe(write_on_this_pe)
        if (.not. this%get_write_on_this_pe()) return
        call this%add_new_axis(naxis)
        do j = 1, size(var_axis_ids)
          if (x_y_axis_id(i) .eq. var_axis_ids(j)) var_axis_ids(j) = naxis
        enddo
      end select
    enddo
  end select
end subroutine define_new_subaxis

!> @brief adds the start time to the fileobj
!! @note This should be called from the register field calls. It can be called multiple times (one for each variable)
!! So it needs to make sure that the start_time is the same for each variable. The initial value is the base_time
subroutine add_start_time(this, start_time, model_time)
  class(fmsDiagFile_type), intent(inout)       :: this           !< The file object
  TYPE(time_type),         intent(in)          :: start_time     !< Start time to add to the fileobj
  TYPE(time_type),         intent(out)         :: model_time     !< The current model time
                                                                 !! this will be set to the start_time
                                                                 !! at the begining of the run

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
    model_time = start_time
    this%start_time = start_time
    this%last_output = start_time
    this%next_output = diag_time_inc(start_time, this%get_file_freq(), this%get_file_frequnit())
    this%next_next_output = diag_time_inc(this%next_output, this%get_file_freq(), this%get_file_frequnit())
    if (this%has_file_new_file_freq()) then
       this%next_close = diag_time_inc(this%start_time, this%get_file_new_file_freq(), &
                                        this%get_file_new_file_freq_units())
    else
      if (this%is_static) then
        ! If the file is static, set the close time to be equal to the start_time, so that it can be closed
        ! after the first write!
        this%next_close = this%start_time
        this%next_next_output = diag_time_inc(this%start_time, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
      else
        this%next_close = diag_time_inc(this%start_time, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
      endif
     endif

    if(this%has_file_duration()) then
       this%no_more_data = diag_time_inc(this%start_time, this%get_file_duration(), &
                                          this%get_file_duration_units())
    else
       this%no_more_data = diag_time_inc(this%start_time, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
    endif

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
  write( unit_num, *)'next_close', date_to_string(this%next_close)

  if( allocated(this%fms2io_fileobj)) write( unit_num, *)'fileobj path', this%fms2io_fileobj%path

  write( unit_num, *)'type_of_domain', this%type_of_domain
  if( allocated(this%file_metadata_from_model)) write( unit_num, *) 'file_metadata_from_model', &
                                                                    this%file_metadata_from_model
  if( allocated(this%field_ids)) write( unit_num, *)'field_ids', this%field_ids
  if( allocated(this%field_registered)) write( unit_num, *)'field_registered', this%field_registered
  if( allocated(this%num_registered_fields)) write( unit_num, *)'num_registered_fields', this%num_registered_fields
  if( allocated(this%axis_ids)) write( unit_num, *)'axis_ids', this%axis_ids(1:this%number_of_axis)

end subroutine

!> @brief Determine if a file is regional
!! @return Flag indicating if the file is regional or not
logical pure function is_regional(this)
  class(fmsDiagFileContainer_type), intent(in) :: this            !< The file object

  select type (wut=>this%FMS_diag_file)
  type is (subRegionalFile_type)
    is_regional = .true.
  type is (fmsDiagFile_type)
    is_regional = .false.
  end select

end function is_regional

!> @brief Determine if a file is static
!! @return Flag indicating if the file is static or not
logical pure function is_file_static(this)
class(fmsDiagFileContainer_type), intent(in) :: this            !< The file object

is_file_static = .false.

select type (fileptr=>this%FMS_diag_file)
type is (fmsDiagFile_type)
  is_file_static = fileptr%is_static
end select

end function is_file_static

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
  !< Go away if it the file is already open
  if (diag_file%is_file_open) return

  is_regional = .false.
  !< Figure out what fms2io_fileobj to use!
  if (.not. allocated(diag_file%fms2io_fileobj)) then
    select type (diag_file)
    type is (subRegionalFile_type)
      !< In this case each PE is going to write its own file
      allocate(FmsNetcdfFile_t :: diag_file%fms2io_fileobj)
      is_regional = .true.
    type is (fmsDiagFile_type)
      !< Use the type_of_domain to get the correct fms2io_fileobj
      select case (diag_file%type_of_domain)
      case (NO_DOMAIN)
        allocate(FmsNetcdfFile_t :: diag_file%fms2io_fileobj)
      case (TWO_D_DOMAIN)
        allocate(FmsNetcdfDomainFile_t :: diag_file%fms2io_fileobj)
      case (UG_DOMAIN)
        allocate(FmsNetcdfUnstructuredDomainFile_t :: diag_file%fms2io_fileobj)
      end select
    end select
  endif

  !< Figure out what to name of the file
  diag_file_name = diag_file%get_file_fname()

  !< If using the new_file_freq figure out what the name is based on the current time
  if (diag_file%has_file_new_file_freq()) then
    !< If using a wildcard file name (i.e ocn%4yr%2mo%2dy%2hr), get the basename (i.e ocn)
    pos = INDEX(diag_file_name, '%')
    if (pos > 0) base_name = diag_file_name(1:pos-1)
    suffix = get_time_string(diag_file_name, diag_file%get_filename_time())
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
  select type (fms2io_fileobj => diag_file%fms2io_fileobj)
  type is (FmsNetcdfFile_t)
    if (is_regional) then
      if (.not. open_file(fms2io_fileobj, file_name, "overwrite", pelist=(/mpp_pe()/))) &
      &call mpp_error(FATAL, "Error opening the file:"//file_name)
      call register_global_attribute(fms2io_fileobj, "is_subregional", "True", str_len=4)
   else
      allocate(pes(mpp_npes()))
      call mpp_get_current_pelist(pes)

      if (.not. open_file(fms2io_fileobj, file_name, "overwrite", pelist=pes)) &
      &call mpp_error(FATAL, "Error opening the file:"//file_name)
   endif
  type is (FmsNetcdfDomainFile_t)
    select type (domain)
    type is (diagDomain2d_t)
      if (.not. open_file(fms2io_fileobj, file_name, "overwrite", domain%Domain2)) &
        &call mpp_error(FATAL, "Error opening the file:"//file_name)
    end select
  type is (FmsNetcdfUnstructuredDomainFile_t)
    select type (domain)
    type is (diagDomainUg_t)
      if (.not. open_file(fms2io_fileobj, file_name, "overwrite", domain%DomainUG)) &
        &call mpp_error(FATAL, "Error opening the file:"//file_name)
    end select
  end select

  file_is_opened = .true.
  diag_file%is_file_open = file_is_opened
  domain => null()
  diag_file => null()
end subroutine open_diag_file

!< @brief Write global attributes in the diag_file
subroutine write_global_metadata(this)
  class(fmsDiagFileContainer_type), intent(inout), target :: this !< The file object

  class(FmsNetcdfFile_t),  pointer  :: fms2io_fileobj !< The fileobj to write to
  integer                           :: i              !< For do loops
  character (len=MAX_STR_LEN), allocatable :: yaml_file_attributes(:,:) !< Global attributes defined in the yaml

  type(diagYamlFiles_type), pointer :: diag_file_yaml !< The diag_file yaml

  diag_file_yaml => this%FMS_diag_file%diag_yaml_file
  fms2io_fileobj => this%FMS_diag_file%fms2io_fileobj

  if (diag_file_yaml%has_file_global_meta()) then
    yaml_file_attributes = diag_file_yaml%get_file_global_meta()
    do i = 1, size(yaml_file_attributes,1)
      call register_global_attribute(fms2io_fileobj, trim(yaml_file_attributes(i,1)), &
      trim(yaml_file_attributes(i,2)), str_len=len_trim(yaml_file_attributes(i,2)))
    enddo
    deallocate(yaml_file_attributes)
  endif
end subroutine write_global_metadata

!< @brief Writes a variable's metadata in the netcdf file
subroutine write_var_metadata(fms2io_fileobj, variable_name, dimensions, long_name, units)
  class(FmsNetcdfFile_t), intent(inout) :: fms2io_fileobj !< The file object to write into
  character(len=*)      , intent(in)    :: variable_name  !< The name of the time variables
  character(len=*)      , intent(in)    :: dimensions(:)  !< The dimensions of the variable
  character(len=*)      , intent(in)    :: long_name      !< The long_name of the variable
  character(len=*)      , intent(in)    :: units          !< The units of the variable

  call register_field(fms2io_fileobj, variable_name, pack_size_str, dimensions)
  call register_variable_attribute(fms2io_fileobj, variable_name, "long_name", &
                                  trim(long_name), str_len=len_trim(long_name))
  if (trim(units) .ne. no_units) &
    call register_variable_attribute(fms2io_fileobj, variable_name, "units", &
                                    trim(units), str_len=len_trim(units))
end subroutine write_var_metadata

!> \brief Write the time metadata to the diag file
subroutine write_time_metadata(this)
  class(fmsDiagFileContainer_type), intent(inout), target :: this !< The file object

  class(fmsDiagFile_type), pointer  :: diag_file      !< Diag_file object to open
  class(FmsNetcdfFile_t),  pointer  :: fms2io_fileobj !< The fileobj to write to
  character(len=50)                 :: time_units_str !< Time units written as a string
  character(len=50)                 :: calendar       !< The calendar name

  character(len=:), allocatable :: time_var_name !< The name of the time variable as it is defined in the yaml
  character(len=50)             :: dimensions(2) !< Array of dimensions names for the variable

  diag_file => this%FMS_diag_file
  fms2io_fileobj => diag_file%fms2io_fileobj

  time_var_name = diag_file%get_file_unlimdim()
  call register_axis(fms2io_fileobj, time_var_name, unlimited)

  WRITE(time_units_str, 11)  &
    TRIM(time_unit_list(diag_file%get_file_timeunit())), get_base_year(),&
         & get_base_month(), get_base_day(), get_base_hour(), get_base_minute(), get_base_second()
11  FORMAT(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)

  dimensions(1) = "nv"
  dimensions(2) = trim(time_var_name)

  call write_var_metadata(fms2io_fileobj, time_var_name, dimensions(2:2), &
    time_var_name, time_units_str)

  !< Add additional variables to the time variable
  call register_variable_attribute(fms2io_fileobj, time_var_name, "axis", "T", str_len=1 )

  !TODO no need to have both attributes, probably?
  calendar = valid_calendar_types(get_calendar_type())
  call register_variable_attribute(fms2io_fileobj, time_var_name, "calendar_type", &
    uppercase(trim(calendar)), str_len=len_trim(calendar))
  call register_variable_attribute(fms2io_fileobj, time_var_name, "calendar", &
    lowercase(trim(calendar)), str_len=len_trim(calendar))

  if (diag_file%time_ops) then
    call register_variable_attribute(fms2io_fileobj, time_var_name, "bounds", &
      trim(time_var_name)//"_bnds", str_len=len_trim(time_var_name//"_bnds"))

    !< It is possible that the "nv" "axis" was registered via "diag_axis_init" call
    !! so only adding it if it doesn't exist already
    if ( .not. dimension_exists(fms2io_fileobj, "nv")) then
      call register_axis(fms2io_fileobj, "nv", 2) !< Time bounds need a vertex number
      call write_var_metadata(fms2io_fileobj, "nv", dimensions(1:1), &
        "vertex number", no_units)
    endif
    call write_var_metadata(fms2io_fileobj, time_var_name//"_bnds", dimensions, &
      trim(time_var_name)//" axis boundaries", time_units_str)
  endif

end subroutine write_time_metadata

!> \brief Write out the field data to the file
subroutine write_field_data(this, field_obj, buffer_obj)
  class(fmsDiagFileContainer_type),        intent(in),    target :: this           !< The diag file object to write to
  type(fmsDiagField_type),                 intent(in),    target :: field_obj      !< The field object to write from
  type(fmsDiagOutputBuffer_type),          intent(inout), target :: buffer_obj     !< The buffer object with the data

  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open
  class(FmsNetcdfFile_t),  pointer     :: fms2io_fileobj !< Fileobj to write to
  logical                              :: has_diurnal    !< indicates if theres a diurnal axis to adjust for

  diag_file => this%FMS_diag_file
  fms2io_fileobj => diag_file%fms2io_fileobj

  !TODO This may be offloaded in the future
  if (diag_file%is_static) then
    !< Here the file is static so there is no need for the unlimited dimension
    !! as a variables are static
    call buffer_obj%write_buffer(fms2io_fileobj)
    diag_file%data_has_been_written = .true.
  else
    if (field_obj%is_static()) then
      !< If the variable is static, only write it the first time
      if (diag_file%unlim_dimension_level .eq. 1) then
        call buffer_obj%write_buffer(fms2io_fileobj)
        diag_file%data_has_been_written = .true.
      endif
    else
     diag_file%data_has_been_written = .true.
     has_diurnal = buffer_obj%get_diurnal_sample_size() .gt. 1
      if (.not. buffer_obj%is_there_data_to_write()) then
        ! Only print the error message once
        if (diag_file%unlim_dimension_level .eq. 1) &
          call mpp_error(NOTE, "Send data was never called. Writing fill values for variable "//&
            field_obj%get_varname()//" in mod "//field_obj%get_modname())
      endif
      call buffer_obj%write_buffer(fms2io_fileobj, &
                        unlim_dim_level=diag_file%unlim_dimension_level, is_diurnal=has_diurnal)
    endif
  endif

end subroutine write_field_data

!> \brief Determine if it is time to close the file
!! \return .True. if it is time to close the file
logical function is_time_to_close_file (this, time_step)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object
  TYPE(time_type),                  intent(in)           :: time_step       !< Current model step time

  if (time_step >= this%FMS_diag_file%next_close) then
    is_time_to_close_file = .true.
  else
    is_time_to_close_file = .false.
  endif
end function

!> \brief Determine if it is time to "write" to the file
logical function is_time_to_write(this, time_step)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object
  TYPE(time_type),                  intent(in)           :: time_step       !< Current model step time

  if (time_step > this%FMS_diag_file%next_output) then
    is_time_to_write = .true.
    if (this%FMS_diag_file%is_static) return
    if (time_step > this%FMS_diag_file%next_next_output) &
      call mpp_error(FATAL, this%FMS_diag_file%get_file_fname()//&
        &": Diag_manager_mod:: You skipped a time_step. Be sure that diag_send_complete is called at every time step "&
        &" needed by the file.")
  else
    is_time_to_write = .false.
    if (this%FMS_diag_file%is_static) then
      ! This is to ensure that static files get finished in the begining of the run
      if (this%FMS_diag_file%unlim_dimension_level .eq. 1) is_time_to_write = .true.
    endif
  endif
end function is_time_to_write

!> \brief Determine if the current PE has data to write
logical function writing_on_this_pe(this)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object

  select type(diag_file => this%FMS_diag_file)
  type is (subRegionalFile_type)
    writing_on_this_pe = diag_file%write_on_this_pe
  class default
    writing_on_this_pe = .true.
  end select

end function

!> \brief Write out the time data to the file
subroutine write_time_data(this)
  class(fmsDiagFileContainer_type), intent(in), target   :: this !< The file object

  real                                 :: dif            !< The time as a real number
  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open
  class(FmsNetcdfFile_t),  pointer     :: fms2io_fileobj !< The fileobj to write to
  TYPE(time_type)                      :: middle_time    !< The middle time of the averaging period

  real :: T1 !< The beginning time of the averaging period
  real :: T2 !< The ending time of the averaging period

  diag_file => this%FMS_diag_file
  fms2io_fileobj => diag_file%fms2io_fileobj

  !< If data has not been written for the current unlimited dimension
  !! ignore this
  if (.not. diag_file%data_has_been_written) return

  if (diag_file%time_ops) then
    middle_time = (diag_file%last_output+diag_file%next_output)/2
    dif = get_date_dif(middle_time, get_base_time(), diag_file%get_file_timeunit())
  else
    dif = get_date_dif(diag_file%next_output, get_base_time(), diag_file%get_file_timeunit())
  endif

  call write_data(fms2io_fileobj, diag_file%get_file_unlimdim(), dif, &
    unlim_dim_level=diag_file%unlim_dimension_level)

  if (diag_file%time_ops) then
    T1 = get_date_dif(diag_file%last_output, get_base_time(), diag_file%get_file_timeunit())
    T2 = get_date_dif(diag_file%next_output, get_base_time(), diag_file%get_file_timeunit())

    call write_data(fms2io_fileobj, trim(diag_file%get_file_unlimdim())//"_bnds", &
                    (/T1, T2/), unlim_dim_level=diag_file%unlim_dimension_level)

    if (diag_file%unlim_dimension_level .eq. 1) then
      call write_data(fms2io_fileobj, "nv", (/1, 2/))
    endif
  endif

end subroutine write_time_data

!> \brief Updates the current_new_file_freq_index if using a new_file_freq
subroutine update_current_new_file_freq_index(this, time_step)
  class(fmsDiagFileContainer_type), intent(inout), target   :: this            !< The file object
  TYPE(time_type),                  intent(in)           :: time_step       !< Current model step time

  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open

  diag_file => this%FMS_diag_file

  if (time_step >= diag_file%no_more_data) then
    call diag_file%diag_yaml_file%increase_new_file_freq_index()

    if (diag_file%has_file_duration()) then
      diag_file%no_more_data = diag_time_inc(diag_file%no_more_data, diag_file%get_file_duration(), &
                                          diag_file%get_file_duration_units())
    else
      !< At this point you are done writing data
       diag_file%done_writing_data = .true.
       diag_file%no_more_data = diag_time_inc(diag_file%no_more_data, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
       diag_file%next_output = diag_file%no_more_data
       diag_file%next_next_output = diag_file%no_more_data
       diag_file%last_output = diag_file%no_more_data
       diag_file%next_close = diag_file%no_more_data
    endif
  endif
end subroutine update_current_new_file_freq_index

!> \brief Set up the next_output and next_next_output variable in a file obj
subroutine update_next_write(this, time_step)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object
  TYPE(time_type),                  intent(in)           :: time_step       !< Current model step time

  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open

  diag_file => this%FMS_diag_file
  if (diag_file%is_static) then
    diag_file%last_output = diag_file%next_output
    diag_file%next_output = diag_time_inc(diag_file%next_output, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
    diag_file%next_next_output = diag_time_inc(diag_file%next_output, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
  else
    diag_file%last_output = diag_file%next_output
    diag_file%next_output = diag_time_inc(diag_file%next_output, diag_file%get_file_freq(), &
      diag_file%get_file_frequnit())
    diag_file%next_next_output = diag_time_inc(diag_file%next_output, diag_file%get_file_freq(), &
      diag_file%get_file_frequnit())
  endif

end subroutine update_next_write

!> \brief Increase the unlimited dimension level that the file is currently being written to
subroutine increase_unlim_dimension_level(this)
  class(fmsDiagFileContainer_type), intent(inout), target   :: this            !< The file object

  this%FMS_diag_file%unlim_dimension_level = this%FMS_diag_file%unlim_dimension_level + 1
  this%FMS_diag_file%data_has_been_written = .false.
end subroutine increase_unlim_dimension_level

!> \brief Get the unlimited dimension level that is in the file
!! \return The unlimited dimension
pure function get_unlim_dimension_level(this) &
result(res)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object
  integer :: res

  res = this%FMS_diag_file%unlim_dimension_level
end function

!> \brief Get the next_output for the file object
!! \return The next_output
pure function get_next_output(this) &
result(res)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object
  type(time_type) :: res

  res = this%FMS_diag_file%next_output
end function get_next_output

!> \brief Get the next_output for the file object
!! \return The next_output
pure function get_next_next_output(this) &
result(res)
  class(fmsDiagFileContainer_type), intent(in), target   :: this            !< The file object
  type(time_type) :: res

  res = this%FMS_diag_file%next_next_output
  if (this%FMS_diag_file%is_static) then
    res = this%FMS_diag_file%no_more_data
  endif
end function get_next_next_output

!< @brief Writes the axis metadata for the file
subroutine write_axis_metadata(this, diag_axis)
  class(fmsDiagFileContainer_type), intent(inout), target :: this            !< The file object
  class(fmsDiagAxisContainer_type), intent(in),    target :: diag_axis(:)    !< Diag_axis object

  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open
  class(FmsNetcdfFile_t),  pointer     :: fms2io_fileobj !< The fileobj to write to
  integer                              :: i,k            !< For do loops
  integer                              :: parent_axis_id !< Id of the parent_axis
  integer                              :: structured_ids(2) !< Ids of the uncompress axis
  integer                              :: edges_id       !< Id of the axis edge

  class(fmsDiagAxisContainer_type), pointer :: axis_ptr      !< pointer to the axis object currently writing
  logical                                   :: edges_in_file !< .true. if the edges are already in the file

  diag_file => this%FMS_diag_file
  fms2io_fileobj => diag_file%fms2io_fileobj

  do i = 1, diag_file%number_of_axis
    edges_in_file = .false.
    axis_ptr => diag_axis(diag_file%axis_ids(i))
    parent_axis_id = axis_ptr%axis%get_parent_axis_id()

    edges_id = axis_ptr%axis%get_edges_id()
    if (edges_id .ne. diag_null) then
      !< write the edges if is not in the list of axis in the file, otherwrise ignore
      if (any(diag_file%axis_ids(1:diag_file%number_of_axis) .eq. edges_id)) then
        edges_in_file = .true.
     else
        call diag_axis(edges_id)%axis%write_axis_metadata(fms2io_fileobj, .true.)
        call diag_file%add_new_axis(edges_id)
      endif
    endif

    if (parent_axis_id .eq. DIAG_NULL) then
      call axis_ptr%axis%write_axis_metadata(fms2io_fileobj, edges_in_file)
    else
      call axis_ptr%axis%write_axis_metadata(fms2io_fileobj, edges_in_file, diag_axis(parent_axis_id)%axis)
    endif

    if (axis_ptr%axis%is_unstructured_grid()) then
      structured_ids = axis_ptr%axis%get_structured_axis()
      do k = 1, size(structured_ids)
        call diag_axis(structured_ids(k))%axis%write_axis_metadata(fms2io_fileobj, .false.)
      enddo
    endif

  enddo

end subroutine write_axis_metadata

!< @brief Writes the field metadata for the file
subroutine write_field_metadata(this, diag_field, diag_axis)
  class(fmsDiagFileContainer_type), intent(inout), target :: this            !< The file object
  class(fmsDiagField_type)        , intent(inout), target :: diag_field(:)   !<
  class(fmsDiagAxisContainer_type), intent(in)            :: diag_axis(:)    !< Diag_axis object

  class(FmsNetcdfFile_t),           pointer :: fms2io_fileobj !< The fileobj to write to
  class(fmsDiagFile_type),          pointer :: diag_file  !< Diag_file object to open
  class(fmsDiagField_type),         pointer :: field_ptr  !< diag_field(diag_file%field_ids(i)), for convenience

  integer            :: i             !< For do loops
  logical            :: is_regional   !< Flag indicating if the field is in a regional file
  character(len=255) :: cell_measures !< cell_measures attributes for the field
  logical            :: need_associated_files !< .True. if the 'associated_files' global attribute is needed
  character(len=255) :: associated_files !< Associated files attribute to add

  is_regional = this%is_regional()

  diag_file => this%FMS_diag_file
  fms2io_fileobj => diag_file%fms2io_fileobj

  associated_files = ""
  need_associated_files = .false.
  do i = 1, size(diag_file%field_ids)
    if (.not. diag_file%field_registered(i)) cycle !TODO do something else here
    field_ptr => diag_field(diag_file%field_ids(i))

    cell_measures = ""
    if (field_ptr%has_area()) then
      cell_measures = "area: "//diag_field(field_ptr%get_area())%get_varname(to_write=.true.)

      !! Determine if the area field is already in the file. If it is not create the "associated_files" attribute
      !! which contains the file name of the file the area field is in. This is needed for PP/fregrid.
      if (.not. diag_field(field_ptr%get_area())%is_variable_in_file(diag_file%id)) then
        need_associated_files = .true.
        call diag_field(field_ptr%get_area())%generate_associated_files_att(associated_files, diag_file%start_time)
      endif
    endif

    if (field_ptr%has_volume()) then
      cell_measures = trim(cell_measures)//" volume: "//diag_field(field_ptr%get_volume())%get_varname(to_write=.true.)

      !! Determine if the volume field is already in the file. If it is not create the "associated_files" attribute
      !! which contains the file name of the file the volume field is in. This is needed for PP/fregrid.
      if (.not. diag_field(field_ptr%get_volume())%is_variable_in_file(diag_file%id)) then
        need_associated_files = .true.
        call diag_field(field_ptr%get_volume())%generate_associated_files_att(associated_files, diag_file%start_time)
      endif
    endif

    call field_ptr%write_field_metadata(fms2io_fileobj, diag_file%id, diag_file%yaml_ids(i), diag_axis, &
      this%FMS_diag_file%get_file_unlimdim(), is_regional, cell_measures)
  enddo

  if (need_associated_files) &
    call register_global_attribute(fms2io_fileobj, "associated_files", trim(ADJUSTL(associated_files)), &
      str_len=len_trim(ADJUSTL(associated_files)))

end subroutine write_field_metadata

!< @brief Writes the axis data for the file
subroutine write_axis_data(this, diag_axis)
  class(fmsDiagFileContainer_type), intent(inout), target :: this            !< The file object
  class(fmsDiagAxisContainer_type), intent(in)            :: diag_axis(:)    !< Diag_axis object

  class(fmsDiagFile_type), pointer     :: diag_file      !< Diag_file object to open
  class(FmsNetcdfFile_t),  pointer     :: fms2io_fileobj !< The fileobj to write to
  integer                              :: i, k           !< For do loops
  integer                              :: j              !< diag_file%axis_ids(i) (for less typing)
  integer                              :: parent_axis_id !< Id of the parent_axis
  integer                              :: structured_ids(2) !< Ids of the uncompress axis

  diag_file => this%FMS_diag_file
  fms2io_fileobj => diag_file%fms2io_fileobj

  do i = 1, diag_file%number_of_axis
    j = diag_file%axis_ids(i)
    parent_axis_id = diag_axis(j)%axis%get_parent_axis_id()
    if (parent_axis_id .eq. DIAG_NULL) then
      call diag_axis(j)%axis%write_axis_data(fms2io_fileobj)
    else
      call diag_axis(j)%axis%write_axis_data(fms2io_fileobj, diag_axis(parent_axis_id)%axis)
    endif

    if (diag_axis(j)%axis%is_unstructured_grid()) then
      structured_ids = diag_axis(j)%axis%get_structured_axis()
      do k = 1, size(structured_ids)
        call diag_axis(structured_ids(k))%axis%write_axis_data(fms2io_fileobj)
      enddo
    endif
  enddo

end subroutine write_axis_data

!< @brief Closes the diag_file
subroutine close_diag_file(this)
  class(fmsDiagFileContainer_type), intent(inout), target :: this            !< The file object

  if (.not. this%FMS_diag_file%is_file_open) return

  !< The select types are needed here because otherwise the code will go to the
  !! wrong close_file routine and things will not close propertly
  select type( fms2io_fileobj => this%FMS_diag_file%fms2io_fileobj)
  type is (FmsNetcdfDomainFile_t)
    call close_file(fms2io_fileobj)
  type is (FmsNetcdfFile_t)
    call close_file(fms2io_fileobj)
  type is (FmsNetcdfUnstructuredDomainFile_t)
    call close_file(fms2io_fileobj)
  end select

  !< Reset the unlimited dimension level back to 0, in case the fms2io_fileobj is re-used
  this%FMS_diag_file%unlim_dimension_level = 0
  this%FMS_diag_file%is_file_open = .false.

  if (this%FMS_diag_file%has_file_new_file_freq()) then
    this%FMS_diag_file%next_close = diag_time_inc(this%FMS_diag_file%next_close, &
                                        this%FMS_diag_file%get_file_new_file_freq(), &
                                        this%FMS_diag_file%get_file_new_file_freq_units())
  else
    this%FMS_diag_file%next_close = diag_time_inc(this%FMS_diag_file%next_close, VERY_LARGE_FILE_FREQ, DIAG_DAYS)
  endif
end subroutine close_diag_file

!> \brief Gets the buffer_id list from the file object
pure function get_buffer_ids (this)
  class(fmsDiagFile_type), intent(in) :: this !< The file object
  integer, allocatable :: get_buffer_ids(:) !< returned buffer ids for this file

  allocate(get_buffer_ids(this%number_of_buffers))
  get_buffer_ids = this%buffer_ids(1:this%number_of_buffers)
end function get_buffer_ids

!> Gets the stored number of buffers from the file object
pure function get_number_of_buffers(this)
  class(fmsDiagFile_type), intent(in) :: this !< file object
  integer :: get_number_of_buffers !< returned number of buffers
  get_number_of_buffers = this%number_of_buffers
end function get_number_of_buffers

#endif
end module fms_diag_file_object_mod
