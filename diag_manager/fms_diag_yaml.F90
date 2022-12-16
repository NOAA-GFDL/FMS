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

!> @defgroup fms_diag_yaml_mod fms_diag_yaml_mod
!> @ingroup diag_manager
!! @brief fms_diag_yaml_mod is an integral part of
!!   diag_manager_mod. Its function is to read the diag_table.yaml to fill in
!!   the diag_yaml_object

!> @file
!> @brief File for @ref diag_yaml_mod

!> @addtogroup fms_diag_yaml_mod
!> @{
module fms_diag_yaml_mod

#define DEBUG .false.
#ifdef use_yaml
use diag_data_mod,   only: DIAG_NULL, DIAG_OCEAN, DIAG_ALL, DIAG_OTHER, set_base_time, latlon_gridtype, &
                           index_gridtype, null_gridtype, DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, &
                           DIAG_MONTHS, DIAG_YEARS, time_average, time_rms, time_max, time_min, time_sum, &
                           time_diurnal, time_power, time_none, r8, i8, r4, i4
use yaml_parser_mod, only: open_and_parse_file, get_value_from_key, get_num_blocks, get_nkeys, &
                           get_block_ids, get_key_value, get_key_ids, get_key_name
use mpp_mod,         only: mpp_error, FATAL, mpp_pe, mpp_root_pe, stdout
use, intrinsic :: iso_c_binding, only : c_ptr, c_null_char
use fms_string_utils_mod, only: fms_array_to_pointer, fms_find_my_string, fms_sort_this, fms_find_unique, &
                                fms_f2c_string
use fms_yaml_output_mod
use platform_mod

implicit none

private

public :: diag_yaml
public :: diag_yaml_object_init, diag_yaml_object_end
public :: diagYamlObject_type, get_diag_yaml_obj, subRegion_type
public :: diagYamlFiles_type, diagYamlFilesVar_type
public :: get_num_unique_fields, find_diag_field, get_diag_fields_entries, get_diag_files_id
public :: dump_diag_yaml_obj
public :: fms_diag_yaml_out
!> @}

integer, parameter :: basedate_size = 6
integer, parameter :: NUM_SUB_REGION_ARRAY = 8
integer, parameter :: MAX_STR_LEN = 255
integer, parameter :: MAX_FREQ = 12


!> @brief type to hold an array of sorted diag_fiels
type varList_type
  character(len=255), allocatable :: var_name(:) !< Array of diag_field
  type(c_ptr), allocatable :: var_pointer(:) !< Array of pointers
  integer, allocatable :: diag_field_indices(:) !< Index of the field in the diag_field array
end type

!> @brief type to hold an array of sorted diag_files
type fileList_type
  character(len=255), allocatable :: file_name(:) !< Array of diag_field
  type(c_ptr), allocatable :: file_pointer(:) !< Array of pointers
  integer, allocatable :: diag_file_indices(:)  !< Index of the file in the diag_file array
end type

!> @brief type to hold the sub region information about a file
type subRegion_type
  INTEGER                        :: grid_type   !< Flag indicating the type of region,
                                                !! acceptable values are latlon_gridtype, index_gridtype,
                                                !! null_gridtype
  class(*),          allocatable :: corners(:,:)!< (x, y) coordinates of the four corner of the region
  integer                        :: tile        !< Tile number of the sub region
                                                !! required if using the "index" grid type

end type subRegion_type

!> @brief type to hold the diag_file information
type diagYamlFiles_type
  private
  character (len=:),    allocatable :: file_fname                        !< file name
  integer                           :: file_frequnit(MAX_FREQ)           !< the frequency unit (DIAG_SECONDS,
                                                                         !! DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS,
                                                                         !! DIAG_YEARS)
  integer                           :: file_freq(MAX_FREQ)               !< the frequency of data
  integer                           :: file_timeunit                     !< The unit of time (DIAG_SECONDS,
                                                                         !! DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS,
                                                                         !! DIAG_YEARS)
  character (len=:),   allocatable :: file_unlimdim                      !< The name of the unlimited dimension
  type(subRegion_type)             :: file_sub_region                    !< type containing info about the subregion
  integer                          :: file_new_file_freq(MAX_FREQ)       !< Frequency for closing the existing file
  integer                          :: file_new_file_freq_units(MAX_FREQ) !< Time units for creating a new file.
                                                                         !! Required if “new_file_freq” used
                                                                         !! (DIAG_SECONDS, DIAG_MINUTES, &
                                                                         !! DIAG_HOURS, DIAG_DAYS, DIAG_YEARS)
  character (len=:),   allocatable :: file_start_time                    !< Time to start the file for the
                                                                         !! first time. Requires “new_file_freq”
  integer                          :: file_duration(MAX_FREQ)            !< How long the file should receive data
                                                                         !! after start time in file_duration_units.
                                                                         !! This optional field can only be used if
                                                                         !! the start_time field is present.  If this
                                                                         !! field is absent, then the file duration
                                                                         !! will be equal to the frequency for
                                                                         !! creating new files. NOTE: The
                                                                         !! file_duration_units field must also
                                                                         !! be present if this field is present.
  integer                          :: file_duration_units(MAX_FREQ)      !< The file duration units
                                                                         !! (DIAG_SECONDS, DIAG_MINUTES, &
                                                                         !! DIAG_HOURS, DIAG_DAYS, DIAG_YEARS)
  integer                          :: current_new_file_freq_index        !< The index of the new_file_freq array
  !< Need to use `MAX_STR_LEN` because not all filenames/global attributes are the same length
  character (len=MAX_STR_LEN), allocatable :: file_varlist(:)            !< An array of variable names
                                                                         !! within a file
  character (len=MAX_STR_LEN), allocatable :: file_global_meta(:,:)      !< Array of key(dim=1)
                                                                         !! and values(dim=2) to be
                                                                         !! added as global meta data to
                                                                         !! the file
 contains

 !> All getter functions (functions named get_x(), for member field named x)
 !! return copies of the member variables unless explicitly noted.
 procedure, public :: size_file_varlist
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
 procedure, public :: is_global_meta
 !> Has functions to determine if allocatable variables are true.  If a variable is not an allocatable
 !! then is will always return .true.
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
 procedure, public :: increase_new_file_freq_index
end type diagYamlFiles_type

!> @brief type to hold the info a diag_field
type diagYamlFilesVar_type
  character (len=:), private, allocatable :: var_fname !< The field/diagnostic name
  character (len=:), private, allocatable :: var_varname !< The name of the variable
  integer          , private, allocatable :: var_reduction !< Reduction to be done on var
                                                           !! time_average, time_rms, time_max,
                                                           !! time_min, time_sum, time_diurnal, time_power
  character (len=:), private, allocatable :: var_module !< The module that th variable is in
  integer          , private, allocatable :: var_kind !< The type/kind of the variable
  character (len=:), private, allocatable :: var_outname !< Name of the variable as written to the file
  character (len=:), private, allocatable :: var_longname !< Overwrites the long name of the variable
  character (len=:), private, allocatable :: var_units !< Overwrites the units
  real(kind=r4_kind), private             :: var_zbounds(2)  !< The z axis limits [vert_min, vert_max]
  integer          , private              :: n_diurnal !< Number of diurnal samples
                                                       !! 0 if var_reduction is not "diurnalXX"
  integer          , private              :: pow_value !< The power value
                                                       !! 0 if pow_value is not "powXX"

  !< Need to use `MAX_STR_LEN` because not all filenames/global attributes are the same length
  character (len=MAX_STR_LEN), dimension (:, :), private, allocatable :: var_attributes !< Attributes to overwrite or
                                                                                        !! add from diag_yaml
 contains
 !> All getter functions (functions named get_x(), for member field named x)
 !! return copies of the member variables unless explicitly noted.
  procedure :: get_var_fname
  procedure :: get_var_varname
  procedure :: get_var_reduction
  procedure :: get_var_module
  procedure :: get_var_kind
  procedure :: get_var_outname
  procedure :: get_var_longname
  procedure :: get_var_units
  procedure :: get_var_zbounds
  procedure :: get_var_attributes
  procedure :: get_n_diurnal
  procedure :: get_pow_value
  procedure :: is_var_attributes

  procedure :: has_var_fname
  procedure :: has_var_varname
  procedure :: has_var_reduction
  procedure :: has_var_module
  procedure :: has_var_kind
  procedure :: has_var_outname
  procedure :: has_var_longname
  procedure :: has_var_units
  procedure :: has_var_zbounds
  procedure :: has_var_attributes
  procedure :: has_n_diurnal
  procedure :: has_pow_value

end type diagYamlFilesVar_type

!> @brief Object that holds the information of the diag_yaml
!> @ingroup fms_diag_yaml_mod
type diagYamlObject_type
  character(len=:), allocatable, private :: diag_title                   !< Experiment name
  integer, private, dimension (basedate_size) :: diag_basedate           !< basedate array
  type(diagYamlFiles_type), allocatable, public, dimension (:) :: diag_files!< History file info
  type(diagYamlFilesVar_type), allocatable, private, dimension (:) :: diag_fields !< Diag fields info
  contains
  procedure :: size_diag_files

  procedure :: get_title        !< Returns the title
  procedure :: get_basedate     !< Returns the basedate array
  procedure :: get_diag_files   !< Returns the diag_files array
  procedure :: get_diag_fields  !< Returns the diag_field array

  procedure :: has_diag_title
  procedure :: has_diag_basedate
  procedure :: has_diag_files
  procedure :: has_diag_fields

end type diagYamlObject_type

type (diagYamlObject_type), target :: diag_yaml  !< Obj containing the contents of the diag_table.yaml
type (varList_type), save :: variable_list !< List of all the variables in the diag_table.yaml
type (fileList_type), save :: file_list !< List of all files in the diag_table.yaml

logical, private :: diag_yaml_module_initialized = .false.


!> @addtogroup fms_diag_yaml_mod
!> @{
contains

!> @brief gets the diag_yaml module variable
!! @return a copy of the diag_yaml module variable
function get_diag_yaml_obj() &
result(res)
  type (diagYamlObject_type) :: res

  res= diag_yaml
end function get_diag_yaml_obj

!> @brief get the basedate of a diag_yaml type
!! @return the basedate as an integer array
pure function get_basedate (diag_yaml) &
result (diag_basedate)
  class (diagYamlObject_type), intent(in) :: diag_yaml               !< The diag_yaml
  integer, dimension (basedate_size) :: diag_basedate !< Basedate array result to return

  diag_basedate = diag_yaml%diag_basedate
end function get_basedate

!> @brief Find the number of files listed in the diag yaml
!! @return the number of files in the diag yaml
pure integer function size_diag_files(diag_yaml)
  class (diagYamlObject_type), intent(in) :: diag_yaml      !< The diag_yaml
  if (diag_yaml%has_diag_files()) then
    size_diag_files = size(diag_yaml%diag_files)
  else
    size_diag_files = 0
  endif
end function size_diag_files

!> @brief get the title of a diag_yaml type
!! @return the title of the diag table as an allocated string
pure function get_title (diag_yaml) &
  result (diag_title)
  class (diagYamlObject_type), intent(in) :: diag_yaml      !< The diag_yaml
  character(len=:),allocatable :: diag_title !< Basedate array result to return

  diag_title = diag_yaml%diag_title
end function get_title

!> @brief get the diag_files of a diag_yaml type
!! @return the diag_files
function get_diag_files(diag_yaml) &
result(diag_files)
  class (diagYamlObject_type), intent(in) :: diag_yaml               !< The diag_yaml
  type(diagYamlFiles_type), allocatable, dimension (:) :: diag_files!< History file info

  diag_files = diag_yaml%diag_files
end function get_diag_files

!> @brief get the diag_fields of a diag_yaml type
!! @return the diag_fields
pure function get_diag_fields(diag_yaml) &
result(diag_fields)
  class (diagYamlObject_type), intent(in) :: diag_yaml               !< The diag_yaml
  type(diagYamlFilesVar_type), allocatable, dimension (:) :: diag_fields !< Diag fields info

  diag_fields = diag_yaml%diag_fields
end function get_diag_fields

!> @brief Uses the yaml_parser_mod to read in the diag_table and fill in the
!! diag_yaml object
subroutine diag_yaml_object_init(diag_subset_output)
  integer, intent(in)  :: diag_subset_output !< DIAG_ALL   - Current PE is in the one and only pelist
                                             !! DIAG_OTHER - Current PE is not in the ocean pelist
                                             !! and there are multiple pelists
                                             !! DIAG_OCEAN - Current PE is in the ocean pelist
                                             !! and there are multiple pelists
  integer              :: diag_yaml_id     !< Id for the diag_table yaml
  integer              :: nfiles           !< Number of files in the diag_table yaml
  integer, allocatable :: diag_file_ids(:) !< Ids of the files in the diag_table yaml
  integer              :: i, j             !< For do loops
  integer              :: total_nvars      !< The total number of variables in the diag_table yaml
  integer              :: var_count        !< The current number of variables added to the diag_yaml obj
  integer              :: file_var_count   !< The current number of variables added in the diag_file
  integer              :: nvars            !< The number of variables in the current file
  integer, allocatable :: var_ids(:)       !< Ids of the variables in diag_table yaml
  logical              :: is_ocean         !< Flag indicating if it is an ocean file
  logical, allocatable :: ignore(:)        !< Flag indicating if the diag_file is going to be ignored
  integer              :: actual_num_files !< The actual number of files that were saved
  integer              :: file_count       !! The current number of files added to the diag_yaml obj
  logical              :: write_file       !< Flag indicating if the user wants the file to be written
  logical              :: write_var        !< Flag indicating if the user wants the variable to be written

  if (diag_yaml_module_initialized) return

  diag_yaml_id = open_and_parse_file("diag_table.yaml")

  call diag_get_value_from_key(diag_yaml_id, 0, "title", diag_yaml%diag_title)
  call get_value_from_key(diag_yaml_id, 0, "base_date", diag_yaml%diag_basedate)
  call set_base_time(diag_yaml%diag_basedate)

  nfiles = get_num_blocks(diag_yaml_id, "diag_files")
  allocate(diag_file_ids(nfiles))
  allocate(ignore(nfiles))

  call get_block_ids(diag_yaml_id, "diag_files", diag_file_ids)

  ignore = .false.
  total_nvars = 0
  !< If you are on two seperate pelists
  if(diag_subset_output .ne. DIAG_ALL) then
    do i = 1, nfiles
      is_ocean = .false.
      call get_value_from_key(diag_yaml_id, diag_file_ids(i), "is_ocean", is_ocean, is_optional=.true.)
      !< If you are on the ocean pelist and the file is not an ocean file, skip the file
      if (diag_subset_output .eq. DIAG_OCEAN .and. .not. is_ocean) ignore(i) = .true.

      !< If you are not on the ocean pelist and the file is ocean, skip the file
      if(diag_subset_output .eq. DIAG_OTHER .and. is_ocean) ignore(i) = .true.
    enddo
  endif

  !< Determine how many files are in the diag_yaml, ignoring those with write_file = False
  actual_num_files = 0
  do i = 1, nfiles
    write_file = .true.
    call get_value_from_key(diag_yaml_id, diag_file_ids(i), "write_file", write_file, is_optional=.true.)
    if(.not. write_file) ignore(i) = .true.

    if (.not. ignore(i)) then
        actual_num_files = actual_num_files + 1
        !< If ignoring the file, ignore the fields in that file too!
        total_nvars = total_nvars + get_total_num_vars(diag_yaml_id, diag_file_ids(i))
    endif
  enddo

  allocate(diag_yaml%diag_files(actual_num_files))
  allocate(diag_yaml%diag_fields(total_nvars))
  allocate(variable_list%var_name(total_nvars))
  allocate(variable_list%diag_field_indices(total_nvars))
  allocate(file_list%file_name(actual_num_files))
  allocate(file_list%diag_file_indices(actual_num_files))

  var_count = 0
  file_count = 0
  !> Loop through the number of nfiles and fill in the diag_yaml obj
  nfiles_loop: do i = 1, nfiles
    if(ignore(i)) cycle
    file_count = file_count + 1
    call diag_yaml_files_obj_init(diag_yaml%diag_files(file_count))
    call fill_in_diag_files(diag_yaml_id, diag_file_ids(i), diag_yaml%diag_files(file_count))

    !> Save the file name in the file_list
    file_list%file_name(file_count) = trim(diag_yaml%diag_files(file_count)%file_fname)//c_null_char
    file_list%diag_file_indices(file_count) = file_count

    nvars = 0
    nvars = get_num_blocks(diag_yaml_id, "varlist", parent_block_id=diag_file_ids(i))
    allocate(var_ids(nvars))
    call get_block_ids(diag_yaml_id, "varlist", var_ids, parent_block_id=diag_file_ids(i))
    file_var_count = 0
    allocate(diag_yaml%diag_files(file_count)%file_varlist(get_total_num_vars(diag_yaml_id, diag_file_ids(i))))
    nvars_loop: do j = 1, nvars
      write_var = .true.
      call get_value_from_key(diag_yaml_id, var_ids(j), "write_var", write_var, is_optional=.true.)
      if (.not. write_var) cycle

      var_count = var_count + 1
      file_var_count = file_var_count + 1

      !> Save the filename in the diag_field type
      diag_yaml%diag_fields(var_count)%var_fname = diag_yaml%diag_files(file_count)%file_fname

      call fill_in_diag_fields(diag_yaml_id, var_ids(j), diag_yaml%diag_fields(var_count))

      !> Save the variable name in the diag_file type
      diag_yaml%diag_files(file_count)%file_varlist(file_var_count) = diag_yaml%diag_fields(var_count)%var_varname

      !> Save the variable name and the module name in the variable_list
      variable_list%var_name(var_count) = trim(diag_yaml%diag_fields(var_count)%var_varname)//&
                                        ":"//trim(diag_yaml%diag_fields(var_count)%var_module)//c_null_char
      variable_list%diag_field_indices(var_count) = var_count
    enddo nvars_loop
    deallocate(var_ids)
  enddo nfiles_loop

  !> Sort the file list in alphabetical order
  file_list%file_pointer = fms_array_to_pointer(file_list%file_name)
  call fms_sort_this(file_list%file_pointer, actual_num_files, file_list%diag_file_indices)

  variable_list%var_pointer = fms_array_to_pointer(variable_list%var_name)
  call fms_sort_this(variable_list%var_pointer, total_nvars, variable_list%diag_field_indices)

  deallocate(diag_file_ids)
  diag_yaml_module_initialized = .true.
end subroutine

!> @brief Destroys the diag_yaml object
subroutine diag_yaml_object_end()
  integer :: i !< For do loops

  do i = 1, size(diag_yaml%diag_files, 1)
    if(allocated(diag_yaml%diag_files(i)%file_varlist)) deallocate(diag_yaml%diag_files(i)%file_varlist)
    if(allocated(diag_yaml%diag_files(i)%file_global_meta)) deallocate(diag_yaml%diag_files(i)%file_global_meta)
    if(allocated(diag_yaml%diag_files(i)%file_sub_region%corners)) &
      deallocate(diag_yaml%diag_files(i)%file_sub_region%corners)
  enddo
  if(allocated(diag_yaml%diag_files)) deallocate(diag_yaml%diag_files)

  do i = 1, size(diag_yaml%diag_fields, 1)
    if(allocated(diag_yaml%diag_fields(i)%var_attributes)) deallocate(diag_yaml%diag_fields(i)%var_attributes)
  enddo
  if(allocated(diag_yaml%diag_fields)) deallocate(diag_yaml%diag_fields)

  if(allocated(file_list%file_pointer)) deallocate(file_list%file_pointer)
  if(allocated(file_list%file_name)) deallocate(file_list%file_name)
  if(allocated(file_list%diag_file_indices)) deallocate(file_list%diag_file_indices)

  if(allocated(variable_list%var_pointer)) deallocate(variable_list%var_pointer)
  if(allocated(variable_list%var_name)) deallocate(variable_list%var_name)
  if(allocated(variable_list%diag_field_indices)) deallocate(variable_list%diag_field_indices)

end subroutine diag_yaml_object_end

!> @brief Fills in a diagYamlFiles_type with the contents of a file block in diag_table.yaml
subroutine fill_in_diag_files(diag_yaml_id, diag_file_id, fileobj)
  integer,                    intent(in)    :: diag_yaml_id !< Id of the diag_table.yaml
  integer,                    intent(in)    :: diag_file_id !< Id of the file block to read
  type(diagYamlFiles_type), intent(inout) :: fileobj      !< diagYamlFiles_type obj to read the contents into

  integer :: nsubregion       !< Flag indicating of there any regions (0 or 1)
  integer :: sub_region_id(1) !< Id of the sub_region block
  integer :: natt             !< Number of global attributes in the current file
  integer :: global_att_id(1) !< Id of the global attributes block
  integer :: nkeys            !< Number of key/value global attributes pair
  integer :: j                !< For do loops

  integer, allocatable :: key_ids(:) !< Id of the gloabl atttributes key/value pairs
  character(len=:), ALLOCATABLE :: grid_type !< grid_type as it is read in from the yaml
  character(len=:), ALLOCATABLE :: freq_buffer !< buffer to store any freq as it is read from the yaml
  character(len=:), ALLOCATABLE :: buffer      !< buffer to store any *_units as it is read from the yaml

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "file_name", fileobj%file_fname)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "freq_units", buffer)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "freq", freq_buffer)
  call set_file_freq(fileobj, freq_buffer, buffer)

  deallocate(freq_buffer, buffer)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "unlimdim", fileobj%file_unlimdim)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "time_units", buffer)
  call set_file_time_units(fileobj, buffer)

  deallocate(buffer)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "new_file_freq", freq_buffer, is_optional=.true.)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "new_file_freq_units", buffer, &
         is_optional=.true.)
  call set_new_file_freq(fileobj, freq_buffer, buffer)

  deallocate(freq_buffer, buffer)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "start_time", fileobj%file_start_time, is_optional=.true.)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "file_duration", freq_buffer, is_optional=.true.)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "file_duration_units", buffer, &
    is_optional=.true.)
  call set_file_duration(fileobj, freq_buffer, buffer)

  nsubregion = 0
  nsubregion = get_num_blocks(diag_yaml_id, "sub_region", parent_block_id=diag_file_id)
  if (nsubregion .eq. 1) then
    call get_block_ids(diag_yaml_id, "sub_region", sub_region_id, parent_block_id=diag_file_id)
    call diag_get_value_from_key(diag_yaml_id, sub_region_id(1), "grid_type", grid_type)
    call get_sub_region(diag_yaml_id, sub_region_id(1), fileobj%file_sub_region, grid_type, fileobj%file_fname)
  elseif (nsubregion .eq. 0) then
    fileobj%file_sub_region%grid_type = null_gridtype
  else
    call mpp_error(FATAL, "diag_yaml_object_init: file "//trim(fileobj%file_fname)//" has multiple region blocks")
  endif

  natt = 0
  natt = get_num_blocks(diag_yaml_id, "global_meta", parent_block_id=diag_file_id)
  if (natt .eq. 1) then
    call get_block_ids(diag_yaml_id, "global_meta", global_att_id, parent_block_id=diag_file_id)
    nkeys = get_nkeys(diag_yaml_id, global_att_id(1))
    allocate(key_ids(nkeys))
    call get_key_ids(diag_yaml_id, global_att_id(1), key_ids)

    allocate(fileobj%file_global_meta(nkeys, 2))
    do j = 1, nkeys
      call get_key_name(diag_yaml_id,  key_ids(j), fileobj%file_global_meta(j, 1))
      call get_key_value(diag_yaml_id, key_ids(j), fileobj%file_global_meta(j, 2))
    enddo
    deallocate(key_ids)
  elseif (natt .ne. 0) then
    call mpp_error(FATAL, "diag_yaml_object_init: file "//trim(fileobj%file_fname)//&
                         &" has multiple global_meta blocks")
  endif

end subroutine

!> @brief Fills in a diagYamlFilesVar_type with the contents of a variable block in
!! diag_table.yaml
subroutine fill_in_diag_fields(diag_file_id, var_id, field)
  integer,                        intent(in)  :: diag_file_id !< Id of the file block in the yaml file
  integer,                        intent(in)  :: var_id       !< Id of the variable block in the yaml file
  type(diagYamlFilesVar_type), intent(inout)  :: field        !< diagYamlFilesVar_type obj to read the contents into

  integer :: natt          !< Number of attributes in variable
  integer :: var_att_id(1) !< Id of the variable attribute block
  integer :: nkeys         !< Number of key/value pairs of attributes
  integer :: j             !< For do loops

  integer, allocatable :: key_ids(:) !< Id of each attribute key/value pair
  character(len=:), ALLOCATABLE :: buffer    !< buffer to store the reduction method as it is read from the yaml

  call diag_get_value_from_key(diag_file_id, var_id, "var_name", field%var_varname)
  call diag_get_value_from_key(diag_file_id, var_id, "reduction", buffer)
  call set_field_reduction(field, buffer)

  call diag_get_value_from_key(diag_file_id, var_id, "module", field%var_module)
  deallocate(buffer)
  call diag_get_value_from_key(diag_file_id, var_id, "kind", buffer)
  call set_field_kind(field, buffer)

  call diag_get_value_from_key(diag_file_id, var_id, "output_name", field%var_outname, is_optional=.true.)
  call diag_get_value_from_key(diag_file_id, var_id, "long_name", field%var_longname, is_optional=.true.)
  !! VAR_UNITS !!

  natt = 0
  natt = get_num_blocks(diag_file_id, "attributes", parent_block_id=var_id)
  if (natt .eq. 1) then
    call get_block_ids(diag_file_id, "attributes", var_att_id, parent_block_id=var_id)
    nkeys = get_nkeys(diag_file_id, var_att_id(1))
    allocate(key_ids(nkeys))
    call get_key_ids(diag_file_id, var_att_id(1), key_ids)

    allocate(field%var_attributes(nkeys, 2))
    do j = 1, nkeys
      call get_key_name(diag_file_id,  key_ids(j), field%var_attributes(j, 1))
      call get_key_value(diag_file_id, key_ids(j), field%var_attributes(j, 2))
    enddo
    deallocate(key_ids)
  elseif (natt .ne. 0) then
      call mpp_error(FATAL, "diag_yaml_object_init: variable "//trim(field%var_varname)//&
                            " has multiple attribute blocks")
  endif

  !> Set the zbounds if they exist
  field%var_zbounds = DIAG_NULL
  call get_value_from_key(diag_file_id, var_id, "zbounds", field%var_zbounds, is_optional=.true.)
end subroutine

!> @brief diag_manager wrapper to get_value_from_key to use for allocatable
!! string variables
subroutine diag_get_value_from_key(diag_file_id, par_id, key_name, value_name, is_optional)
  integer,           intent(in) :: diag_file_id!< Id of the file block in the yaml file
  integer,           intent(in) :: par_id      !< Id of the parent block in the yaml file
  character(len=*),  intent(in) :: key_name    !< Key to look for in the parent block
  character(len=:), allocatable :: value_name  !< Value of the key
  logical, intent(in), optional :: is_optional !< Flag indicating if the key is optional

  character(len=255) :: buffer !< String buffer to read in to

  buffer = "" !< Needs to be initialized for optional keys that are not present
  call get_value_from_key(diag_file_id, par_id, trim(key_name), buffer, is_optional= is_optional)
  allocate(character(len=len_trim(buffer)) :: value_name)
  value_name = trim(buffer)

end subroutine diag_get_value_from_key

!> @brief gets the lat/lon of the sub region to use in a diag_table yaml
subroutine get_sub_region(diag_yaml_id, sub_region_id, sub_region, grid_type, fname)
  integer,             intent(in)    :: diag_yaml_id       !< Id of the diag_table yaml file
  integer,             intent(in)    :: sub_region_id      !< Id of the region block to read from
  type(subRegion_type),intent(inout) :: sub_region         !< Type that stores the sub_region
  character(len=*),    intent(in)    :: grid_type          !< The grid_type as it is read from the file
  character(len=*),    intent(in)    :: fname              !< filename of the subregion (for error messages)

  select case (trim(grid_type))
  case ("latlon")
    sub_region%grid_type = latlon_gridtype
    allocate(real(kind=r4_kind) :: sub_region%corners(4,2))
  case ("index")
    sub_region%grid_type = index_gridtype
    allocate(integer(kind=i4_kind) :: sub_region%corners(4,2))

    call get_value_from_key(diag_yaml_id, sub_region_id, "tile", sub_region%tile, is_optional=.true.)
    if (sub_region%tile .eq. DIAG_NULL) call mpp_error(FATAL, &
      "The tile number is required when defining a "//&
      "subregion. Check your subregion entry for "//trim(fname))
  case default
    call mpp_error(FATAL, trim(grid_type)//" is not a valid region type. &
    &The acceptable values are latlon and index. &
    &Check your entry for file:"//trim(fname))
  end select

  call get_value_from_key(diag_yaml_id, sub_region_id, "corner1", sub_region%corners(1,:))
  call get_value_from_key(diag_yaml_id, sub_region_id, "corner2", sub_region%corners(2,:))
  call get_value_from_key(diag_yaml_id, sub_region_id, "corner3", sub_region%corners(3,:))
  call get_value_from_key(diag_yaml_id, sub_region_id, "corner4", sub_region%corners(4,:))

end subroutine get_sub_region

!> @brief gets the total number of variables in the diag_table yaml file
!! @return total number of variables
function get_total_num_vars(diag_yaml_id, diag_file_id) &
result(total_nvars)

  integer, intent(in) :: diag_yaml_id     !< Id for the diag_table yaml
  integer, intent(in) :: diag_file_id     !< Id of the file in the diag_table yaml
  integer :: total_nvars

  integer :: i !< For do loop
  integer :: nvars !< Number of variables in a file
  integer, allocatable :: var_ids(:) !< Id of the variables in the file block of the yaml file
  logical :: var_write !< Flag indicating if the user wants the variable to be written

  nvars = get_num_blocks(diag_yaml_id, "varlist", parent_block_id=diag_file_id)
  allocate(var_ids(nvars))
  call get_block_ids(diag_yaml_id, "varlist", var_ids, parent_block_id=diag_file_id)

  !< Loop through all the variables in the diag_file block and only count those that don't have write_var=false
  total_nvars = 0
  do i = 1, nvars
    var_write = .true.
    call get_value_from_key(diag_yaml_id, var_ids(i), "write_var", var_write, is_optional=.true.)
    if (var_write) total_nvars = total_nvars + 1
  end do
end function

!> @brief This checks if the file frequency and file frequency units in a diag file are valid and
!! sets the integer equivalent
subroutine set_file_freq(fileobj, file_freq, file_frequnit)
  type(diagYamlFiles_type), intent(inout) :: fileobj         !< diagYamlFiles_type obj to check
  character(len=*),         intent(in)    :: file_freq       !< File_freq as it is read from the diag_table
  character(len=*),         intent(in)    :: file_frequnit   !< File_freq_units as it is read from the diag_table

  integer           :: i                         !< For do loops
  character(len=10) :: file_freq_units(MAX_FREQ) !< Array of file frequencies as a string
  integer           :: err_unit                  !< Dummy error unit

  file_freq_units = ""
  read(file_freq, *, iostat=err_unit) fileobj%file_freq
  read(file_frequnit, *, iostat=err_unit) file_freq_units

  do i = 1, MAX_FREQ
    if (fileobj%file_freq(i) >= -1) then
      if (trim(file_freq_units(i)) .eq. "") &
        call mpp_error(FATAL, "file_freq_units is required. &
        &Check your entry for file:"//trim(fileobj%file_fname))

      fileobj%file_frequnit(i) = set_valid_time_units(file_freq_units(i), &
      "file_freq_units for file:"//trim(fileobj%file_fname))
    else
      return
    endif
  enddo
end subroutine set_file_freq

!> @brief This checks if the time unit in a diag file is valid and sets the integer equivalent
subroutine set_file_time_units (fileobj, file_timeunit)
  type(diagYamlFiles_type), intent(inout) :: fileobj       !< diagYamlFiles_type obj to checK
  character(len=*),         intent(in)    :: file_timeunit !< file_timeunit as it is read from the diag_table

 fileobj%file_timeunit = set_valid_time_units(file_timeunit, "timeunit for file:"//trim(fileobj%file_fname))
end subroutine set_file_time_units

!> @brief This checks if the new file frequency and the new file frequency units in a diag file are valid
!! and sets the integer equivalent
subroutine set_new_file_freq(fileobj, new_file_freq, new_file_freq_units)
  type(diagYamlFiles_type), intent(inout) :: fileobj                  !< diagYamlFiles_type obj to check
  character(len=*),         intent(in)    :: new_file_freq            !< new file freq units as it is read from
                                                                      !! the diag_table
  character(len=*),         intent(in)    :: new_file_freq_units      !< new file freq units as it is read from
                                                                      !! the diag_table
  integer           :: i                                  !< For do loops
  character(len=10) :: file_new_file_freq_units(MAX_FREQ) !< Array of new file frequencies as string
  integer           :: err_unit                           !< Dummy error unit

  file_new_file_freq_units = ""
  read(new_file_freq, *, iostat=err_unit) fileobj%file_new_file_freq
  read(new_file_freq_units, *, iostat=err_unit) file_new_file_freq_units

  do i = 1, MAX_FREQ
    if (fileobj%file_new_file_freq(i) > 0) then
      if (trim(file_new_file_freq_units(i)) .eq. "") &
        call mpp_error(FATAL, "new_file_freq_units is required if using new_file_freq. &
        &Check your entry for file:"//trim(fileobj%file_fname))

      fileobj%file_new_file_freq_units(i) = set_valid_time_units(file_new_file_freq_units(i), &
      "new_file_freq_units for file:"//trim(fileobj%file_fname))
    else
      return
    endif
  enddo
end subroutine set_new_file_freq

!> @brief This checks if the file duration and the file duration units in a diag file are valid
!! and sets the integer equivalent
subroutine set_file_duration(fileobj, file_duration, file_duration_units)
  type(diagYamlFiles_type), intent(inout) :: fileobj             !< diagYamlFiles_type obj to check
  character(len=*),         intent(in)    :: file_duration       !< file_duration as it is read from the yaml
  character(len=*),         intent(in)    :: file_duration_units !< file_duration units as it is read from the yaml

  integer           :: i                                   !< For do loops
  character(len=10) :: file_duration_units_array(MAX_FREQ) !< Array of file_duration_units as string
  integer           :: err_unit                            !< Dummy error unit
  logical           :: mask(MAX_FREQ)                      !< Array of logical
  integer           :: nfile_duration                      !< Number of file durations defined
  integer           :: nfile_freq                          !< Number of file frequencies defined
  integer           :: nnew_file_freq                      !< Number of new file frequencies defined

  file_duration_units_array = ""
  read(file_duration, *, iostat=err_unit) fileobj%file_duration
  read(file_duration_units, *, iostat=err_unit) file_duration_units_array

  nfile_duration = 0
  do i = 1, MAX_FREQ
    if (fileobj%file_duration(i) > 0) then
      if(trim(file_duration_units_array(i)) .eq. "") &
      call mpp_error(FATAL, "file_duration_units is required if using file_duration. &
        &Check your entry for file:"//trim(fileobj%file_fname))

      fileobj%file_duration_units(i) = set_valid_time_units(file_duration_units_array(i), &
      "file_duration_units for file:"//trim(fileobj%file_fname))
      nfile_duration = nfile_duration + 1
    else
      exit
    endif
  enddo

  !< Make sure the user send in the correct number of freq, new_file_freq, and file_duration
  mask = .FALSE.
  mask = fileobj%file_freq .ne. DIAG_NULL
  nfile_freq = count(mask)

  mask = .FALSE.
  mask = fileobj%file_new_file_freq .ne. DIAG_NULL
  nnew_file_freq = count(mask)

  if (nfile_freq .ne. nfile_duration .and. nfile_freq-1 .ne. nfile_duration) &
    call mpp_error(FATAL, "freq and file_duration do not have consistent size. &
      &Check your entry for file:"//trim(fileobj%file_fname))

  if (nfile_freq .ne. nnew_file_freq .and. nfile_freq-1 .ne. nnew_file_freq) &
    call mpp_error(FATAL, "freq and new_file_freq do not have consistent size. &
      &Check your entry for file:"//trim(fileobj%file_fname))

  if (nnew_file_freq .ne. nfile_duration .and. nnew_file_freq-1 .ne. nfile_duration) &
    call mpp_error(FATAL, "new_file_freq and file_duration do not have consistent size. &
      &Check your entry for file:"//trim(fileobj%file_fname))

end subroutine set_file_duration

!> @brief This checks if the kind of a diag field is valid and sets it
subroutine set_field_kind(field, skind)
  type(diagYamlFilesVar_type), intent(inout) :: field        !< diagYamlFilesVar_type obj to read the contents into
  character(len=*),            intent(in)    :: skind        !< The variable kind as read from diag_yaml

  select case (TRIM(skind))
  case ("r4")
    field%var_kind = r4
  case ("r8")
    field%var_kind = r8
  case ("i4")
    field%var_kind = i4
  case ("i8")
    field%var_kind = i8
  case default
    call mpp_error(FATAL, trim(skind)//" is an invalid kind! &
      &The acceptable values are r4, r8, i4, i8. &
      &Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
  end select

end subroutine set_field_kind

!> @brief This checks if the reduction of a diag field is valid and sets it
!! If the reduction method is diurnalXX or powXX, it gets the number of diurnal sample and the power value
subroutine set_field_reduction(field, reduction_method)
  type(diagYamlFilesVar_type), intent(inout) :: field           !< diagYamlFilesVar_type obj to read the contents into
  character(len=*)           , intent(in)    :: reduction_method!< reduction method as read from the yaml

  integer :: n_diurnal !< number of diurnal samples
  integer :: pow_value !< The power value
  integer :: ioerror   !< io error status after reading in the diurnal samples

  n_diurnal = 0
  pow_value = 0
  ioerror = 0
  if (index(reduction_method, "diurnal") .ne. 0) then
    READ (UNIT=reduction_method(8:LEN_TRIM(reduction_method)), FMT=*, IOSTAT=ioerror) n_diurnal
    if (ioerror .ne. 0) &
      call mpp_error(FATAL, "Error getting the number of diurnal samples from "//trim(reduction_method))
    if (n_diurnal .le. 0) &
      call mpp_error(FATAL, "Diurnal samples should be greater than 0. &
        & Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
    field%var_reduction = time_diurnal
  elseif (index(reduction_method, "pow") .ne. 0) then
    READ (UNIT=reduction_method(4:LEN_TRIM(reduction_method)), FMT=*, IOSTAT=ioerror) pow_value
    if (ioerror .ne. 0) &
      call mpp_error(FATAL, "Error getting the power value from "//trim(reduction_method))
      if (pow_value .le. 0) &
      call mpp_error(FATAL, "The power value should be greater than 0. &
        & Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
    field%var_reduction = time_power
  else
    select case (reduction_method)
    case ("none")
      field%var_reduction = time_none
    case ("average")
      field%var_reduction = time_average
    case ("min")
      field%var_reduction = time_min
    case ("max")
      field%var_reduction = time_max
    case ("rms")
      field%var_reduction = time_rms
    case ("sum")
      field%var_reduction = time_sum
    case default
      call mpp_error(FATAL, trim(reduction_method)//" is an invalid reduction method! &
        &The acceptable values are none, average, pow##, diurnal##, min, max, and rms. &
        &Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
    end select
  endif

  field%n_diurnal = n_diurnal
  field%pow_value = pow_value
end subroutine set_field_reduction

!> @brief This checks if a time unit is valid and if it is, it assigns the integer equivalent
!! @return The integer equivalent to the time units
function set_valid_time_units(time_units, error_msg) &
result(time_units_int)

  character(len=*), intent(in)  :: time_units      !< The time_units as a string
  character(len=*), intent(in)  :: error_msg       !< Error message to append

  integer                       :: time_units_int  !< The integer equivalent of the time_units

  select case (TRIM(time_units))
  case ("seconds")
    time_units_int = DIAG_SECONDS
  case ("minutes")
    time_units_int = DIAG_MINUTES
  case ("hours")
    time_units_int = DIAG_HOURS
  case ("days")
    time_units_int = DIAG_DAYS
  case ("months")
    time_units_int = DIAG_MONTHS
  case ("years")
    time_units_int = DIAG_YEARS
  case default
    time_units_int =DIAG_NULL
    call mpp_error(FATAL, trim(error_msg)//" is not valid. Acceptable values are "&
                   "seconds, minutes, hours, days, months, years")
  end select
end function set_valid_time_units

!!!!!!! YAML FILE INQUIRIES !!!!!!!
!> @brief Finds the number of variables in the file_varlist
!! @return the size of the diag_files_obj%file_varlist array
integer pure function size_file_varlist (diag_files_obj)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 size_file_varlist = size(diag_files_obj%file_varlist)
end function size_file_varlist

!> @brief Inquiry for diag_files_obj%file_fname
!! @return file_fname of a diag_yaml_file obj
pure function get_file_fname (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_files_obj%file_fname
end function get_file_fname
!> @brief Inquiry for diag_files_obj%file_frequnit
!! @return file_frequnit of a diag_yaml_file_obj
pure function get_file_frequnit (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_frequnit(diag_files_obj%current_new_file_freq_index)
end function get_file_frequnit
!> @brief Inquiry for diag_files_obj%file_freq
!! @return file_freq of a diag_yaml_file_obj
pure function get_file_freq(diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_freq(diag_files_obj%current_new_file_freq_index)
end function get_file_freq
!> @brief Inquiry for diag_files_obj%file_timeunit
!! @return file_timeunit of a diag_yaml_file_obj
pure function get_file_timeunit (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_timeunit
end function get_file_timeunit
!> @brief Inquiry for diag_files_obj%file_unlimdim
!! @return file_unlimdim of a diag_yaml_file_obj
pure function get_file_unlimdim(diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_files_obj%file_unlimdim
end function get_file_unlimdim
!> @brief Inquiry for diag_files_obj%file_subregion
!! @return file_sub_region of a diag_yaml_file_obj
function get_file_sub_region (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 type(subRegion_type) :: res !< What is returned
  res = diag_files_obj%file_sub_region
end function get_file_sub_region
!> @brief Inquiry for diag_files_obj%file_new_file_freq
!! @return file_new_file_freq of a diag_yaml_file_obj
pure function get_file_new_file_freq(diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_new_file_freq(diag_files_obj%current_new_file_freq_index)
end function get_file_new_file_freq
!> @brief Inquiry for diag_files_obj%file_new_file_freq_units
!! @return file_new_file_freq_units of a diag_yaml_file_obj
pure function get_file_new_file_freq_units (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_new_file_freq_units(diag_files_obj%current_new_file_freq_index)
end function get_file_new_file_freq_units
!> @brief Inquiry for diag_files_obj%file_start_time
!! @return file_start_time of a diag_yaml_file_obj
pure function get_file_start_time (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_files_obj%file_start_time
end function get_file_start_time
!> @brief Inquiry for diag_files_obj%file_duration
!! @return file_duration of a diag_yaml_file_obj
pure function get_file_duration (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_duration(diag_files_obj%current_new_file_freq_index)
end function get_file_duration
!> @brief Inquiry for diag_files_obj%file_duration_units
!! @return file_duration_units of a diag_yaml_file_obj
pure function get_file_duration_units (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
  integer :: res !< What is returned
  res = diag_files_obj%file_duration_units(diag_files_obj%current_new_file_freq_index)
end function get_file_duration_units
!> @brief Inquiry for diag_files_obj%file_varlist
!! @return file_varlist of a diag_yaml_file_obj
pure function get_file_varlist (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (:), allocatable :: res(:) !< What is returned
  res = diag_files_obj%file_varlist
end function get_file_varlist
!> @brief Inquiry for diag_files_obj%file_global_meta
!! @return file_global_meta of a diag_yaml_file_obj
pure function get_file_global_meta (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (:), allocatable :: res(:,:) !< What is returned
  res = diag_files_obj%file_global_meta
end function get_file_global_meta
!> @brief Inquiry for whether file_global_meta is allocated
!! @return Flag indicating if file_global_meta is allocated
function is_global_meta(diag_files_obj) &
  result(res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 logical :: res
 res = .false.
 if (allocated(diag_files_obj%file_global_meta)) &
   res = .true.
end function

!> @brief Increate the current_new_file_freq_index by 1
subroutine increase_new_file_freq_index(this)
  class(diagYamlFiles_type), intent(inout) :: this !< The file object
  this%current_new_file_freq_index = this%current_new_file_freq_index + 1
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! VARIABLES ROUTINES AND FUNCTIONS !!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! YAML VAR INQUIRIES !!!!!!!
!> @brief Inquiry for diag_yaml_files_var_obj%var_fname
!! @return var_fname of a diag_yaml_files_var_obj
pure function get_var_fname (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_fname
end function get_var_fname
!> @brief Inquiry for diag_yaml_files_var_obj%var_varname
!! @return var_varname of a diag_yaml_files_var_obj
pure function get_var_varname (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_varname
end function get_var_varname
!> @brief Inquiry for diag_yaml_files_var_obj%var_reduction
!! @return var_reduction of a diag_yaml_files_var_obj
pure function get_var_reduction (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 integer, allocatable :: res !< What is returned
  res = diag_var_obj%var_reduction
end function get_var_reduction
!> @brief Inquiry for diag_yaml_files_var_obj%var_module
!! @return var_module of a diag_yaml_files_var_obj
pure function get_var_module (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_module
end function get_var_module
!> @brief Inquiry for diag_yaml_files_var_obj%var_kind
!! @return var_kind of a diag_yaml_files_var_obj
pure function get_var_kind (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 integer, allocatable :: res !< What is returned
  res = diag_var_obj%var_kind
end function get_var_kind
!> @brief Inquiry for diag_yaml_files_var_obj%var_outname
!! @return var_outname of a diag_yaml_files_var_obj
pure function get_var_outname (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_outname
end function get_var_outname
!> @brief Inquiry for diag_yaml_files_var_obj%var_longname
!! @return var_longname of a diag_yaml_files_var_obj
pure function get_var_longname (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_longname
end function get_var_longname
!> @brief Inquiry for diag_yaml_files_var_obj%var_units
!! @return var_units of a diag_yaml_files_var_obj
pure function get_var_units (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_units
end function get_var_units
!> @brief Inquiry for diag_yaml_files_var_obj%var_zbounds
!! @return var_zbounds of a diag_yaml_files_var_obj
pure function get_var_zbounds (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 real(kind=r4_kind) :: res(2) !< What is returned
  res = diag_var_obj%var_zbounds
end function get_var_zbounds
!> @brief Inquiry for diag_yaml_files_var_obj%var_attributes
!! @return var_attributes of a diag_yaml_files_var_obj
pure function get_var_attributes(diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=MAX_STR_LEN), allocatable :: res (:,:) !< What is returned
 res = diag_var_obj%var_attributes
end function get_var_attributes
!> @brief Inquiry for diag_yaml_files_var_obj%n_diurnal
!! @return the number of diurnal samples of a diag_yaml_files_var_obj
pure function get_n_diurnal(diag_var_obj) &
result (res)
  class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
  integer :: res !< What is returned
  res = diag_var_obj%n_diurnal
end function get_n_diurnal
!> @brief Inquiry for diag_yaml_files_var_obj%pow_value
!! @return the pow_value of a diag_yaml_files_var_obj
pure function get_pow_value(diag_var_obj) &
result (res)
  class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
  integer :: res !< What is returned
  res = diag_var_obj%pow_value
end function get_pow_value
!> @brief Inquiry for whether var_attributes is allocated
!! @return Flag indicating if var_attributes is allocated
function is_var_attributes(diag_var_obj) &
result(res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 logical :: res
 res = .false.
 if (allocated(diag_var_obj%var_attributes)) &
   res = .true.
end function is_var_attributes

!> @brief Initializes the non string values of a diagYamlFiles_type to its
!! default values
subroutine diag_yaml_files_obj_init(obj)
  type(diagYamlFiles_type), intent(out) :: obj !< diagYamlFiles_type object to initialize

  obj%file_freq           = DIAG_NULL
  obj%file_sub_region%tile = DIAG_NULL
  obj%file_new_file_freq = DIAG_NULL
  obj%file_duration = DIAG_NULL
  obj%file_new_file_freq_units = DIAG_NULL
  obj%file_duration_units = DIAG_NULL
  obj%current_new_file_freq_index = 1
end subroutine diag_yaml_files_obj_init

!> @brief Checks if obj%file_fname is allocated
!! @return true if obj%file_fname is allocated
pure logical function has_file_fname (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_fname = allocated(obj%file_fname)
end function has_file_fname
!> @brief Checks if obj%file_frequnit is allocated
!! @return true if obj%file_frequnit is allocated
pure logical function has_file_frequnit (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_frequnit = obj%file_frequnit(obj%current_new_file_freq_index) .NE. DIAG_NULL
end function has_file_frequnit
!> @brief obj%file_freq is on the stack, so the object always has it
!! @return true if obj%file_freq is allocated
pure logical function has_file_freq (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_freq = .true.
end function has_file_freq
!> @brief Checks if obj%file_timeunit is allocated
!! @return true if obj%file_timeunit is allocated
pure logical function has_file_timeunit (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_timeunit = obj%file_timeunit .ne. diag_null
end function has_file_timeunit
!> @brief Checks if obj%file_unlimdim is allocated
!! @return true if obj%file_unlimdim is allocated
pure logical function has_file_unlimdim (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_unlimdim = allocated(obj%file_unlimdim)
end function has_file_unlimdim
!> @brief Checks if obj%file_write is on the stack, so this will always be true
!! @return true
pure logical function has_file_write (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_write = .true.
end function has_file_write
!> @brief Checks if obj%file_sub_region is being used and has the sub region variables allocated
!! @return true if obj%file_sub_region sub region variables are allocated
pure logical function has_file_sub_region (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  if ( obj%file_sub_region%grid_type .eq. latlon_gridtype .or. obj%file_sub_region%grid_type .eq. index_gridtype) then
       has_file_sub_region = .true.
  else
       has_file_sub_region = .false.
  endif
end function has_file_sub_region
!> @brief obj%file_new_file_freq is defined on the stack, so this will return true
!! @return true
pure logical function has_file_new_file_freq (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_new_file_freq = obj%file_new_file_freq(obj%current_new_file_freq_index) .ne. DIAG_NULL
end function has_file_new_file_freq
!> @brief Checks if obj%file_new_file_freq_units is allocated
!! @return true if obj%file_new_file_freq_units is allocated
pure logical function has_file_new_file_freq_units (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_new_file_freq_units = obj%file_new_file_freq_units(obj%current_new_file_freq_index) .ne. diag_null
end function has_file_new_file_freq_units
!> @brief Checks if obj%file_start_time is allocated
!! @return true if obj%file_start_time is allocated
pure logical function has_file_start_time (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_start_time = allocated(obj%file_start_time)
end function has_file_start_time
!> @brief obj%file_duration is allocated on th stack, so this is always true
!! @return true
pure logical function has_file_duration (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_duration = obj%file_duration(obj%current_new_file_freq_index) .ne. DIAG_NULL
end function has_file_duration
!> @brief obj%file_duration_units is on the stack, so this will retrun true
!! @return true
pure logical function has_file_duration_units (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_duration_units = obj%file_duration_units(obj%current_new_file_freq_index) .ne. diag_null
end function has_file_duration_units
!> @brief Checks if obj%file_varlist is allocated
!! @return true if obj%file_varlist is allocated
pure logical function has_file_varlist (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_varlist = allocated(obj%file_varlist)
end function has_file_varlist
!> @brief Checks if obj%file_global_meta is allocated
!! @return true if obj%file_global_meta is allocated
pure logical function has_file_global_meta (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_global_meta = allocated(obj%file_global_meta)
end function has_file_global_meta

!> @brief Checks if obj%var_fname is allocated
!! @return true if obj%var_fname is allocated
pure logical function has_var_fname (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_fname = allocated(obj%var_fname)
end function has_var_fname
!> @brief Checks if obj%var_varname is allocated
!! @return true if obj%var_varname is allocated
pure logical function has_var_varname (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_varname = allocated(obj%var_varname)
end function has_var_varname
!> @brief Checks if obj%var_reduction is allocated
!! @return true if obj%var_reduction is allocated
pure logical function has_var_reduction (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_reduction = allocated(obj%var_reduction)
end function has_var_reduction
!> @brief Checks if obj%var_module is allocated
!! @return true if obj%var_module is allocated
pure logical function has_var_module (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_module = allocated(obj%var_module)
end function has_var_module
!> @brief Checks if obj%var_kind is allocated
!! @return true if obj%var_kind is allocated
pure logical function has_var_kind (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_kind = allocated(obj%var_kind)
end function has_var_kind
!> @brief obj%var_write is on the stack, so this returns true
!! @return true
pure logical function has_var_write (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_write = .true.
end function has_var_write
!> @brief Checks if obj%var_outname is allocated
!! @return true if obj%var_outname is allocated
pure logical function has_var_outname (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_outname = allocated(obj%var_outname)
end function has_var_outname
!> @brief Checks if obj%var_longname is allocated
!! @return true if obj%var_longname is allocated
pure logical function has_var_longname (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_longname = allocated(obj%var_longname)
end function has_var_longname
!> @brief Checks if obj%var_units is allocated
!! @return true if obj%var_units is allocated
pure logical function has_var_units (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_units = allocated(obj%var_units)
end function has_var_units
!> @brief Checks if obj%var_zbounds is allocated
!! @return true if obj%var_zbounds is allocated
pure logical function has_var_zbounds (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_zbounds = any(obj%var_zbounds .eq. diag_null)
end function has_var_zbounds
!> @brief Checks if obj%var_attributes is allocated
!! @return true if obj%var_attributes is allocated
pure logical function has_var_attributes (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_attributes = allocated(obj%var_attributes)
end function has_var_attributes
!> @brief Checks if obj%n_diurnal is set
!! @return true if obj%n_diurnal is set
pure logical function has_n_diurnal(obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to inquire
  has_n_diurnal = (obj%n_diurnal .ne. 0)
end function has_n_diurnal
!> @brief Checks if obj%pow_value is set
!! @return true if obj%pow_value is set
pure logical function has_pow_value(obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to inquire
  has_pow_value = (obj%pow_value .ne. 0)
end function has_pow_value

!> @brief Checks if obj%diag_title is allocated
!! @return true if obj%diag_title is allocated
pure logical function has_diag_title (obj)
  class(diagYamlObject_type), intent(in) :: obj !< diagYamlObject_type object to initialize
  has_diag_title = allocated(obj%diag_title)
end function has_diag_title
!> @brief obj%diag_basedate is on the stack, so this is always true
!! @return true
pure logical function has_diag_basedate (obj)
  class(diagYamlObject_type), intent(in) :: obj !< diagYamlObject_type object to initialize
  has_diag_basedate = .true.
end function has_diag_basedate
!> @brief Checks if obj%diag_files is allocated
!! @return true if obj%diag_files is allocated
pure logical function has_diag_files (obj)
  class(diagYamlObject_type), intent(in) :: obj !< diagYamlObject_type object to initialize
  has_diag_files = allocated(obj%diag_files)
end function has_diag_files
!> @brief Checks if obj%diag_fields is allocated
!! @return true if obj%diag_fields is allocated
pure logical function has_diag_fields (obj)
  class(diagYamlObject_type), intent(in) :: obj !< diagYamlObject_type object to initialize
  has_diag_fields = allocated(obj%diag_fields)
end function has_diag_fields

!> @brief Determine the number of unique diag_fields in the diag_yaml_object
!! @return The number of unique diag_fields
function get_num_unique_fields() &
  result(nfields)
  integer :: nfields
  nfields = fms_find_unique(variable_list%var_pointer, size(variable_list%var_pointer))

end function get_num_unique_fields

!> @brief Determines if a diag_field is in the diag_yaml_object
!! @return Indices of the locations where the field was found
function find_diag_field(diag_field_name, module_name) &
result(indices)

  character(len=*), intent(in) :: diag_field_name !< diag_field name to search for
  character(len=*), intent(in) :: module_name     !< Name of the module, the variable is in

  integer, allocatable :: indices(:)

  indices = fms_find_my_string(variable_list%var_pointer, size(variable_list%var_pointer), &
                               & trim(diag_field_name)//":"//trim(module_name)//c_null_char)
end function find_diag_field

!> @brief Gets the diag_field entries corresponding to the indices of the sorted variable_list
!! @return Array of diag_fields
function get_diag_fields_entries(indices) &
  result(diag_field)

  integer, intent(in) :: indices(:) !< Indices of the field in the sorted variable_list array
  type(diagYamlFilesVar_type), dimension (:), allocatable :: diag_field

  integer :: i !< For do loops
  integer :: field_id !< Indices of the field in the diag_yaml array

  allocate(diag_field(size(indices)))

  do i = 1, size(indices)
    field_id = variable_list%diag_field_indices(indices(i))
    diag_field(i) = diag_yaml%diag_fields(field_id)
  end do

end function get_diag_fields_entries

!> @brief Finds the indices of the diag_yaml%diag_files(:) corresponding to fields in variable_list(indices)
!! @return indices of the diag_yaml%diag_files(:)
function get_diag_files_id(indices) &
  result(file_id)

  integer, intent(in) :: indices(:) !< Indices of the field in the sorted variable_list
  integer, allocatable :: file_id(:)

  integer :: field_id !< Indices of the field in the diag_yaml field array
  integer :: i !< For do loops
  character(len=120) :: filename !< Filename of the field
  integer, allocatable :: file_indices(:) !< Indices of the file in the sorted variable_list

  allocate(file_id(size(indices)))

  do i = 1, size(indices)
    field_id = variable_list%diag_field_indices(indices(i))
    !< Get the filename of the field
    filename = diag_yaml%diag_fields(field_id)%var_fname

    !< File indice of that file in the array of list of sorted files
    file_indices = fms_find_my_string(file_list%file_pointer, size(file_list%file_pointer), &
      & trim(filename)//c_null_char)

    if (size(file_indices) .ne. 1) &
      & call mpp_error(FATAL, "get_diag_files_id: Error getting the correct number of file indices!")

    if (file_indices(1) .eq. diag_null) &
      & call mpp_error(FATAL, "get_diag_files_id: Error finding the filename in the diag_files yaml")

    !< Get the index of the file in the diag_yaml file
    file_id(i) = file_list%diag_file_indices(file_indices(1))
  end do

end function get_diag_files_id

!> Prints out values from diag_yaml object for debugging.
!! Only writes on root.
subroutine dump_diag_yaml_obj( filename )
  character(len=*), optional, intent(in)        :: filename !< optional name of logfile to write to, otherwise
                                                            !! prints to stdout
  type(diagyamlfilesvar_type), allocatable      :: fields(:)
  type(diagyamlfiles_type), allocatable         :: files(:)
  integer                                       :: i, unit_num
  if( present(filename)) then
    open(newunit=unit_num, file=trim(filename), action='WRITE')
  else
    unit_num = stdout()
  endif
  !! TODO write to log
  if( mpp_pe() .eq. mpp_root_pe()) then
    write(unit_num, *) '**********Dumping diag_yaml object**********'
    if( diag_yaml%has_diag_title())    write(unit_num, *) 'Title:', diag_yaml%diag_title
    if( diag_yaml%has_diag_basedate()) write(unit_num, *) 'basedate array:', diag_yaml%diag_basedate
    write(unit_num, *) 'FILES'
    allocate(fields(SIZE(diag_yaml%get_diag_fields())))
    allocate(files(SIZE(diag_yaml%get_diag_files())))
    files = diag_yaml%get_diag_files()
    fields = diag_yaml%get_diag_fields()
    do i=1, SIZE(files)
      write(unit_num, *) 'File: ', files(i)%get_file_fname()
      if(files(i)%has_file_frequnit()) write(unit_num, *) 'file_frequnit:', files(i)%get_file_frequnit()
      if(files(i)%has_file_freq()) write(unit_num, *) 'freq:', files(i)%get_file_freq()
      if(files(i)%has_file_timeunit()) write(unit_num, *) 'timeunit:', files(i)%get_file_timeunit()
      if(files(i)%has_file_unlimdim()) write(unit_num, *) 'unlimdim:', files(i)%get_file_unlimdim()
      !if(files(i)%has_file_sub_region()) write(unit_num, *) 'sub_region:', files(i)%get_file_sub_region()
      if(files(i)%has_file_new_file_freq()) write(unit_num, *) 'new_file_freq:', files(i)%get_file_new_file_freq()
      if(files(i)%has_file_new_file_freq_units()) write(unit_num, *) 'new_file_freq_units:', &
                                                         & files(i)%get_file_new_file_freq_units()
      if(files(i)%has_file_start_time()) write(unit_num, *) 'start_time:', files(i)%get_file_start_time()
      if(files(i)%has_file_duration()) write(unit_num, *) 'duration:', files(i)%get_file_duration()
      if(files(i)%has_file_duration_units()) write(unit_num, *) 'duration_units:', files(i)%get_file_duration_units()
      if(files(i)%has_file_varlist()) write(unit_num, *) 'varlist:', files(i)%get_file_varlist()
      if(files(i)%has_file_global_meta()) write(unit_num, *) 'global_meta:', files(i)%get_file_global_meta()
      if(files(i)%is_global_meta()) write(unit_num, *) 'global_meta:', files(i)%is_global_meta()
      write(unit_num, *) ''
    enddo
    write(unit_num, *) 'FIELDS'
    do i=1, SIZE(fields)
      write(unit_num, *) 'Field: ', fields(i)%get_var_fname()
      if(fields(i)%has_var_fname()) write(unit_num, *) 'fname:', fields(i)%get_var_fname()
      if(fields(i)%has_var_varname()) write(unit_num, *) 'varname:', fields(i)%get_var_varname()
      if(fields(i)%has_var_reduction()) write(unit_num, *) 'reduction:', fields(i)%get_var_reduction()
      if(fields(i)%has_var_module()) write(unit_num, *) 'module:', fields(i)%get_var_module()
      if(fields(i)%has_var_kind()) write(unit_num, *) 'kind:', fields(i)%get_var_kind()
      if(fields(i)%has_var_outname()) write(unit_num, *) 'outname:', fields(i)%get_var_outname()
      if(fields(i)%has_var_longname()) write(unit_num, *) 'longname:', fields(i)%get_var_longname()
      if(fields(i)%has_var_units()) write(unit_num, *) 'units:', fields(i)%get_var_units()
      if(fields(i)%has_var_zbounds()) write(unit_num, *) 'zbounds:', fields(i)%get_var_zbounds()
      if(fields(i)%has_var_attributes()) write(unit_num, *) 'attributes:', fields(i)%get_var_attributes()
      if(fields(i)%has_n_diurnal()) write(unit_num, *) 'n_diurnal:', fields(i)%get_n_diurnal()
      if(fields(i)%has_pow_value()) write(unit_num, *) 'pow_value:', fields(i)%get_pow_value()
      if(fields(i)%has_var_attributes()) write(unit_num, *) 'is_var_attributes:', fields(i)%is_var_attributes()
    enddo
    deallocate(files, fields)
    if( present(filename)) then
      close(unit_num)
    endif
  endif
end subroutine

!> Writes an output yaml with all available information on the written files.
!! Will only write with root pe.
!! Global attributes are limited to 16 per file.
subroutine fms_diag_yaml_out()
  type(diagYamlFiles_type), pointer :: fileptr !< pointer for individual variables
  type(diagYamlFilesVar_type), pointer :: varptr !< pointer for individual variables
  type (fmsyamloutkeys_type), allocatable :: keys(:), keys2(:), keys3(:)
  type (fmsyamloutvalues_type), allocatable :: vals(:), vals2(:), vals3(:)
  integer :: i, j
  character(len=128) :: tmpstr1, tmpstr2 !< string to store output fields
  integer, parameter :: tier1size = 3 !< size of first tier, will always be 3 for basedate, title and diag_files
  integer :: tier2size, tier3size !< size of each 'tier'(based one numbers of tabs) in the yaml
  integer, allocatable :: tier3each(:) !< tier 3 list sizes corresponding to where they are in the second tier
  integer, dimension(basedate_size) :: basedate_loc !< local copy of basedate to loop through
  integer :: varnum_i, key3_i, gm
  character(len=32), allocatable :: st_vals(:) !< start times for gcc bug

  if( mpp_pe() .ne. mpp_root_pe()) return

  allocate(tier3each(SIZE(diag_yaml%diag_files) * 3))
  tier3size = 0; tier3each = 0

  !! allocations for key+val structs
  allocate(keys(1))
  allocate(vals(1))
  allocate(keys2(SIZE(diag_yaml%diag_files)))
  allocate(vals2(SIZE(diag_yaml%diag_files)))
  allocate(st_vals(SIZE(diag_yaml%diag_files)))
  do i=1, SIZE(diag_yaml%diag_files)
    call initialize_key_struct(keys2(i))
    call initialize_val_struct(vals2(i))
    if (allocated(diag_yaml%diag_files(i)%file_varlist) ) then
      do j=1, SIZE(diag_yaml%diag_files(i)%file_varlist)
        tier3size = tier3size + 1
      enddo
    endif
    tier3size = tier3size + 2
  enddo
  allocate(keys3(tier3size))
  allocate(vals3(tier3size))

  !! tier 1 - title, basedate, diag_files
  call initialize_key_struct(keys(1))
  call initialize_val_struct(vals(1))
  call fms_f2c_string( keys(1)%key1, 'title')
  call fms_f2c_string( vals(1)%val1, diag_yaml%diag_title)
  call fms_f2c_string( keys(1)%key2, 'base_date')
  basedate_loc = diag_yaml%get_basedate()
  tmpstr1 = ''; tmpstr2 = ''
  write (tmpstr1, '(I0)') basedate_loc(1)
  tmpstr2 = trim(tmpstr1)
  do i=2, basedate_size
    write (tmpstr1, '(I0)') basedate_loc(i)
    tmpstr2 = trim(tmpstr2) // ' ' // trim(tmpstr1)
  enddo
  call fms_f2c_string(vals(1)%val2, trim(tmpstr2))
  call yaml_out_add_level2key('diag_files', keys(1))
  key3_i = 0
  !! tier 2 - diag files
  do i=1, SIZE(diag_yaml%diag_files)
    fileptr => diag_yaml%diag_files(i)

    call fms_f2c_string(keys2(i)%key1, 'file_name')
    call fms_f2c_string(keys2(i)%key2, 'freq')
    call fms_f2c_string(keys2(i)%key3, 'freq_units')
    call fms_f2c_string(keys2(i)%key4, 'time_units')
    call fms_f2c_string(keys2(i)%key5, 'unlimdim')
    call fms_f2c_string(keys2(i)%key6, 'new_file_freq')
    call fms_f2c_string(keys2(i)%key7, 'new_file_freq_units')
    !! gcc bug causes this key to sometimes not show up
    call fms_f2c_string(keys2(i)%key8, 'start_time')
    call fms_f2c_string(keys2(i)%key9, 'file_duration')
    call fms_f2c_string(keys2(i)%key10, 'file_duration_units')
    if (fileptr%has_file_fname())     call fms_f2c_string(vals2(i)%val1, fileptr%file_fname)
    if (fileptr%has_file_unlimdim())  call fms_f2c_string(vals2(i)%val5, fileptr%file_unlimdim)
    if (fileptr%has_file_timeunit()) then
      if ( get_diag_unit_string(fileptr%file_timeunit) .eq. 'N/A' ) call mpp_error(FATAL, &
            "fms_diag_yaml_out: invalid time unit for file:"//fileptr%file_fname)
      call fms_f2c_string(vals2(i)%val4, get_diag_unit_string(fileptr%file_timeunit))
    endif
    if (fileptr%has_file_freq()) then
      write(tmpstr1, '(I0)') fileptr%file_freq
      call fms_f2c_string(vals2(i)%val2, trim(tmpstr1))
    endif
    if (fileptr%has_file_frequnit()) then
      if ( get_diag_unit_string(fileptr%file_frequnit) .eq. 'N/A' ) call mpp_error(FATAL, &
            "fms_diag_yaml_out: invalid frequency unit for file:"//fileptr%file_fname)
      call fms_f2c_string(vals2(i)%val3, get_diag_unit_string(fileptr%file_frequnit))
    endif
    if(fileptr%has_file_new_file_freq()) then
      write(tmpstr1, '(I0)') fileptr%file_new_file_freq
      call fms_f2c_string(vals2(i)%val6, tmpstr1)
    endif
    if(fileptr%has_file_new_file_freq_units()) then
      if ( get_diag_unit_string(fileptr%file_new_file_freq_units) .eq. 'N/A' ) call mpp_error(FATAL, &
            "fms_diag_yaml_out: invalid new file frequency unit for file:"//fileptr%file_fname)
      call fms_f2c_string(vals2(i)%val7, get_diag_unit_string(fileptr%file_new_file_freq_units))
    endif
    if(fileptr%has_file_start_time()) then
      call fms_f2c_string(vals(i)%val8, trim(fileptr%get_file_start_time()))
      st_vals(i) = fileptr%get_file_start_time()
    endif
    if(fileptr%has_file_duration()) then
      write(tmpstr1, '(I0)') fileptr%file_duration
      call fms_f2c_string(vals2(i)%val9, tmpstr1)
    endif
    if(fileptr%has_file_duration_units()) then
      if ( get_diag_unit_string(fileptr%file_duration_units) .eq. 'N/A' ) call mpp_error(FATAL, &
            "fms_diag_yaml_out: invalid new file duration unit for file:"//fileptr%file_fname)
      call fms_f2c_string(vals2(i)%val10, get_diag_unit_string(fileptr%file_duration_units))
    endif
    !! tier 3 - varlists, subregion, global metadata
    call yaml_out_add_level2key('varlist', keys2(i))
    j = 0
    if( SIZE(fileptr%file_varlist) .gt. 0) then
      do j=1, SIZE(fileptr%file_varlist)
        key3_i = key3_i + 1
        call initialize_key_struct(keys3(key3_i))
        call initialize_val_struct(vals3(key3_i))
        !! find the variable object from the list
        varptr => NULL()
        do varnum_i=1, SIZE(diag_yaml%diag_fields)
          if(DEBUG) print *, 'diag_obj', trim(diag_yaml%diag_fields(varnum_i)%var_varname)
          if(DEBUG) print *, 'file obj', trim(diag_yaml%diag_fields(varnum_i)%var_varname)
          if( trim(diag_yaml%diag_fields(varnum_i)%var_varname ) .eq. trim(fileptr%file_varlist(j)) .and. &
              trim(diag_yaml%diag_fields(varnum_i)%var_fname) .eq. trim(fileptr%file_fname)) then
            if(DEBUG) print *, 'match'
            varptr => diag_yaml%diag_fields(varnum_i)
            exit
          endif
        enddo
        if( .not. associated(varptr)) call mpp_error(FATAL, "diag_yaml_output: could not find variable in list."//&
                                                            " var: "// trim(fileptr%file_varlist(j)))
        call fms_f2c_string(keys3(key3_i)%key1, 'module')
        call fms_f2c_string(keys3(key3_i)%key2, 'var_name')
        call fms_f2c_string(keys3(key3_i)%key3, 'reduction')
        call fms_f2c_string(keys3(key3_i)%key4, 'kind')
        call fms_f2c_string(keys3(key3_i)%key5, 'output_name')
        call fms_f2c_string(keys3(key3_i)%key6, 'long_name')
        call fms_f2c_string(keys3(key3_i)%key7, 'units')
        call fms_f2c_string(keys3(key3_i)%key8, 'zbounds')
        call fms_f2c_string(keys3(key3_i)%key9, 'n_diurnal')
        call fms_f2c_string(keys3(key3_i)%key10, 'pow_value')
        if (varptr%has_var_module())   call fms_f2c_string(vals3(key3_i)%val1, varptr%var_module)
        if (varptr%has_var_varname())  call fms_f2c_string(vals3(key3_i)%val2, varptr%var_varname)
        if( get_diag_reduction_string(varptr%var_reduction) .eq. 'N/A') call mpp_error(FATAL,  &
                    "fms_diag_yaml_out: invalid reduction type for variable:" //varptr%var_varname)
        if (varptr%has_var_reduction())call fms_f2c_string(vals3(key3_i)%val3, &
                                                           get_diag_reduction_string(varptr%var_reduction))
        if (varptr%has_var_outname())  call fms_f2c_string(vals3(key3_i)%val5, varptr%var_outname)
        if (varptr%has_var_longname()) call fms_f2c_string(vals3(key3_i)%val6, varptr%var_longname)
        if (varptr%has_var_units()) call fms_f2c_string(vals3(key3_i)%val7, varptr%var_units)
        if (varptr%has_var_kind()) then
          select case(varptr%var_kind)
            case(i4)
              call fms_f2c_string(vals3(key3_i)%val4, 'i4')
            case(i8)
              call fms_f2c_string(vals3(key3_i)%val4, 'i8')
            case(r4)
              call fms_f2c_string(vals3(key3_i)%val4, 'r4')
            case(r8)
              call fms_f2c_string(vals3(key3_i)%val4, 'r8')
          end select
        endif
        if (varptr%has_var_zbounds()) then
          tmpstr1 = ''
          !! trims don't seem to work with just (F) for the format code
          !! these are just abritrary lengths
          write (tmpstr1, '(F10.5)') varptr%var_zbounds(1)
          tmpstr2 = trim(tmpstr1)
          write (tmpstr1, '(F10.5)') varptr%var_zbounds(2)
          tmpstr2 = trim(tmpstr2) // ', ' // trim(tmpstr1)
          call fms_f2c_string(vals3(key3_i)%val8, trim(tmpstr2))
        endif
        if (varptr%has_n_diurnal()) then
          tmpstr1 = ''; write(tmpstr1, '(I0)') varptr%n_diurnal
          call fms_f2c_string(vals3(key3_i)%val9, tmpstr1)
        endif
        if (varptr%has_pow_value()) then
          tmpstr1 = ''; write(tmpstr1, '(I0)') varptr%pow_value
          call fms_f2c_string(vals3(key3_i)%val10, tmpstr1)
        endif
      enddo
    endif

    key3_i = key3_i + 1
    tier3each(i*3-2) = j-1 ! j-1 structs to print for varlist keys
    tier3each(i*3-1) = 1   ! 1 struct per sub_region key
    tier3each(i*3) = 1     ! 1 struct per global metadata key
    call initialize_key_struct(keys3(key3_i))
    call initialize_val_struct(vals3(key3_i))
    !! sub region
    call yaml_out_add_level2key('sub_region', keys2(i))
    call fms_f2c_string(keys3(key3_i)%key1, 'grid_type')
    call fms_f2c_string(keys3(key3_i)%key2, 'tile')
    call fms_f2c_string(keys3(key3_i)%key3, 'corner1')
    call fms_f2c_string(keys3(key3_i)%key4, 'corner2')
    call fms_f2c_string(keys3(key3_i)%key5, 'corner3')
    call fms_f2c_string(keys3(key3_i)%key6, 'corner4')
    select case (fileptr%file_sub_region%grid_type)
      case(latlon_gridtype)
        call fms_f2c_string(vals3(key3_i)%val1, 'latlon')
      case(index_gridtype)
        call fms_f2c_string(vals3(key3_i)%val1, 'index')
      case default
        call fms_f2c_string(vals3(key3_i)%val1, 'null')
    end select
    tmpstr1 = ''; write(tmpstr1, '(I0)') fileptr%file_sub_region%tile
    call fms_f2c_string(vals3(key3_i)%val2, tmpstr1)
    if( allocated(fileptr%file_sub_region%corners)) then
      select type (corners => fileptr%file_sub_region%corners)
      type is (real(r8_kind))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(1,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(1,2)
        call fms_f2c_string(vals3(key3_i)%val3, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(2,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(2,2)
        call fms_f2c_string(vals3(key3_i)%val4, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(3,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(3,2)
        call fms_f2c_string(vals3(key3_i)%val5, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(4,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(4,2)
        call fms_f2c_string(vals3(key3_i)%val6, trim(tmpstr1)//' '//trim(tmpstr2))
      type is (real(r4_kind))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(1,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(1,2)
        call fms_f2c_string(vals3(key3_i)%val3, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(2,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(2,2)
        call fms_f2c_string(vals3(key3_i)%val4, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(3,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(3,2)
        call fms_f2c_string(vals3(key3_i)%val5, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(4,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(4,2)
        call fms_f2c_string(vals3(key3_i)%val6, trim(tmpstr1)//' '//trim(tmpstr2))
      type is (integer(i4_kind))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(1,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(1,2)
        call fms_f2c_string(vals3(key3_i)%val3, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(2,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(2,2)
        call fms_f2c_string(vals3(key3_i)%val4, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(3,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(3,2)
        call fms_f2c_string(vals3(key3_i)%val5, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(4,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(4,2)
        call fms_f2c_string(vals3(key3_i)%val6, trim(tmpstr1)//' '//trim(tmpstr2))
      type is (integer(i8_kind))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(1,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(1,2)
        call fms_f2c_string(vals3(key3_i)%val3, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(2,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(2,2)
        call fms_f2c_string(vals3(key3_i)%val4, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(3,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(3,2)
        call fms_f2c_string(vals3(key3_i)%val5, trim(tmpstr1)//' '//trim(tmpstr2))
        tmpstr1 = ''; write(tmpstr1, '(I0)') corners(4,1)
        tmpstr2 = ''; write(tmpstr2, '(I0)') corners(4,2)
        call fms_f2c_string(vals3(key3_i)%val6, trim(tmpstr1)//' '//trim(tmpstr2))
      end select
    endif
    !! global metadata
    key3_i = key3_i + 1
    call initialize_key_struct(keys3(key3_i))
    call initialize_val_struct(vals3(key3_i))
    call yaml_out_add_level2key('global_meta', keys2(i))
    if ( fileptr%has_file_global_meta()) then
      do gm=1, SIZE(fileptr%file_global_meta, 1)
        select case(gm)
          case (1)
            call fms_f2c_string(keys3(key3_i)%key1, fileptr%file_global_meta(1,1))
            call fms_f2c_string(vals3(key3_i)%val1, fileptr%file_global_meta(1,2))
          case (2)
            call fms_f2c_string(keys3(key3_i)%key2, fileptr%file_global_meta(2,1))
            call fms_f2c_string(vals3(key3_i)%val2, fileptr%file_global_meta(2,2))
          case (3)
            call fms_f2c_string(keys3(key3_i)%key3, fileptr%file_global_meta(3,1))
            call fms_f2c_string(vals3(key3_i)%val3, fileptr%file_global_meta(3,2))
          case (4)
            call fms_f2c_string(keys3(key3_i)%key4, fileptr%file_global_meta(4,1))
            call fms_f2c_string(vals3(key3_i)%val4, fileptr%file_global_meta(4,2))
          case (5)
            call fms_f2c_string(keys3(key3_i)%key5, fileptr%file_global_meta(5,1))
            call fms_f2c_string(vals3(key3_i)%val5, fileptr%file_global_meta(5,2))
          case (6)
            call fms_f2c_string(keys3(key3_i)%key6, fileptr%file_global_meta(6,1))
            call fms_f2c_string(vals3(key3_i)%val6, fileptr%file_global_meta(6,2))
          case (7)
            call fms_f2c_string(keys3(key3_i)%key7, fileptr%file_global_meta(7,1))
            call fms_f2c_string(vals3(key3_i)%val7, fileptr%file_global_meta(7,2))
          case (8)
            call fms_f2c_string(keys3(key3_i)%key8, fileptr%file_global_meta(8,1))
            call fms_f2c_string(vals3(key3_i)%val8, fileptr%file_global_meta(8,2))
          case (9)
            call fms_f2c_string(keys3(key3_i)%key9, fileptr%file_global_meta(9,1))
            call fms_f2c_string(vals3(key3_i)%val9, fileptr%file_global_meta(9,2))
          case (10)
            call fms_f2c_string(keys3(key3_i)%key10, fileptr%file_global_meta(10,1))
            call fms_f2c_string(vals3(key3_i)%val10, fileptr%file_global_meta(10,2))
          case (11)
            call fms_f2c_string(keys3(key3_i)%key11, fileptr%file_global_meta(11,1))
            call fms_f2c_string(vals3(key3_i)%val11, fileptr%file_global_meta(11,2))
          case (12)
            call fms_f2c_string(keys3(key3_i)%key12, fileptr%file_global_meta(12,1))
            call fms_f2c_string(vals3(key3_i)%val12, fileptr%file_global_meta(12,2))
          case (13)
            call fms_f2c_string(keys3(key3_i)%key13, fileptr%file_global_meta(13,1))
            call fms_f2c_string(vals3(key3_i)%val13, fileptr%file_global_meta(13,2))
          case (14)
            call fms_f2c_string(keys3(key3_i)%key14, fileptr%file_global_meta(14,1))
            call fms_f2c_string(vals3(key3_i)%val14, fileptr%file_global_meta(14,2))
          case (15)
            call fms_f2c_string(keys3(key3_i)%key15, fileptr%file_global_meta(15,1))
            call fms_f2c_string(vals3(key3_i)%val15, fileptr%file_global_meta(15,2))
          case (16)
            call fms_f2c_string(keys3(key3_i)%key16, fileptr%file_global_meta(16,1))
            call fms_f2c_string(vals3(key3_i)%val16, fileptr%file_global_meta(16,2))
        end select
      enddo
    endif
  enddo
  tier2size = i
  ! gcc bug seems to occasionally ignore the start_time key when its set in the loop above
  ! so we'll just set it again
  do i=1, SIZE(keys2)
    call fms_f2c_string(keys2(i)%key8, 'start_time')
    call fms_f2c_string(vals2(i)%val8, st_vals(i))
  enddo
  if (DEBUG) print *, 'tier1size', 1, 'tier2size', SIZE(diag_yaml%diag_files), &
                                                      & 'tier3size', tier3size, 'tier3each', tier3each
  call write_yaml_from_struct_3( 'diag_out.yaml'//c_null_char,  1, keys, vals,          &
                                 SIZE(diag_yaml%diag_files), keys2, vals2, &
                                 tier3size, tier3each, keys3, vals3,       &
                                 (/size(diag_yaml%diag_files), 0, 0, 0, 0, 0, 0, 0/))
  deallocate( keys, keys2, keys3, vals, vals2, vals3)

end subroutine

!> private function for getting unit string from diag_data parameter values
character(len=7) pure function get_diag_unit_string( unit_param )
  integer, intent(in) :: unit_param !< A diag unit parameter value from diag_data_mod.
                                    !! <br>eg. DIAG_SECONDS, DIAG_MINUTES,DIAG_HOURS, DIAG_DAYS, DIAG_YEARS
  select case(unit_param)
    case (DIAG_SECONDS)
      get_diag_unit_string = 'seconds'
    case (DIAG_MINUTES)
      get_diag_unit_string = 'minutes'
    case (DIAG_HOURS)
      get_diag_unit_string = 'hours'
    case (DIAG_DAYS)
      get_diag_unit_string = 'days'
    case (DIAG_YEARS)
      get_diag_unit_string = 'years'
    case default
      get_diag_unit_string = 'N/A'
  end select
end function

!> private function for getting reduction type string from parameter values
character(len=7) pure function get_diag_reduction_string( reduction_val ) &
  result(rslt)
  integer, intent(in) :: reduction_val !< reduction type (eg. time_average)
  select case (reduction_val)
    case (time_none)
      rslt = 'none'
    case (time_average)
      rslt = 'average'
    case (time_min)
      rslt = 'min'
    case (time_max)
      rslt = 'max'
    case (time_rms)
      rslt = 'rms'
    case (time_sum)
      rslt = 'sum'
    case default
      rslt = 'N/A'
  end select
end function

#endif
end module fms_diag_yaml_mod
!> @}
! close documentation grouping
