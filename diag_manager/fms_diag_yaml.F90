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
#ifdef use_yaml
use diag_data_mod,   only: DIAG_NULL, DIAG_OCEAN, DIAG_ALL, DIAG_OTHER, set_base_time, latlon_gridtype, &
                           index_gridtype, null_gridtype, DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, &
                           DIAG_MONTHS, DIAG_YEARS, time_average, time_rms, time_max, time_min, time_sum, &
                           time_diurnal, time_power, time_none, r8, i8, r4, i4, DIAG_NOT_REGISTERED, &
                           middle_time, begin_time, end_time, MAX_STR_LEN
use yaml_parser_mod, only: open_and_parse_file, get_value_from_key, get_num_blocks, get_nkeys, &
                           get_block_ids, get_key_value, get_key_ids, get_key_name
use mpp_mod,         only: mpp_error, FATAL, mpp_pe, mpp_root_pe, stdout
use, intrinsic :: iso_c_binding, only : c_ptr, c_null_char
use fms_string_utils_mod, only: fms_array_to_pointer, fms_find_my_string, fms_sort_this, fms_find_unique
use platform_mod, only: r4_kind, i4_kind
use fms_mod, only: lowercase

implicit none

private

public :: diag_yaml
public :: diag_yaml_object_init, diag_yaml_object_end
public :: diagYamlObject_type, get_diag_yaml_obj, subRegion_type
public :: diagYamlFiles_type, diagYamlFilesVar_type
public :: get_num_unique_fields, find_diag_field, get_diag_fields_entries, get_diag_files_id
public :: get_diag_field_ids
public :: dump_diag_yaml_obj
!> @}

integer, parameter :: basedate_size = 6
integer, parameter :: NUM_SUB_REGION_ARRAY = 8
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
  integer                          :: filename_time                      !< The time to use when setting the name of
                                                                         !! new files: begin, middle, or end of the
                                                                         !! time_bounds
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
 procedure, public :: get_filename_time
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
  type(diagYamlFilesVar_type), allocatable, public, dimension (:) :: diag_fields !< Diag fields info
  contains
  procedure :: size_diag_files

  procedure :: get_title        !< Returns the title
  procedure :: get_basedate     !< Returns the basedate array
  procedure :: get_diag_files   !< Returns the diag_files array
  procedure :: get_diag_fields  !< Returns the diag_field array
  procedure :: get_diag_field_from_id

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
pure function get_basedate (this) &
result (diag_basedate)
  class (diagYamlObject_type), intent(in) :: this          !< The diag_yaml
  integer, dimension (basedate_size)      :: diag_basedate !< Basedate array result to return

  diag_basedate = this%diag_basedate
end function get_basedate

!> @brief Find the number of files listed in the diag yaml
!! @return the number of files in the diag yaml
pure integer function size_diag_files(this)
  class (diagYamlObject_type), intent(in) :: this      !< The diag_yaml
  if (this%has_diag_files()) then
    size_diag_files = size(this%diag_files)
  else
    size_diag_files = 0
  endif
end function size_diag_files

!> @brief get the title of a diag_yaml type
!! @return the title of the diag table as an allocated string
pure function get_title (this) &
  result (diag_title)
  class (diagYamlObject_type), intent(in) :: this       !< The diag_yaml
  character(len=:),allocatable            :: diag_title !< Basedate array result to return

  diag_title = this%diag_title
end function get_title

!> @brief get the diag_files of a diag_yaml type
!! @return the diag_files
function get_diag_files(this) &
result(diag_files)
  class (diagYamlObject_type), intent(in)              :: this       !< The diag_yaml
  type(diagYamlFiles_type), allocatable, dimension (:) :: diag_files !< History file info

  diag_files = this%diag_files
end function get_diag_files

!> @brief Get the diag_field yaml corresponding to a yaml_id
!! @return Pointer to the diag_field yaml entry
function get_diag_field_from_id(this, yaml_id) &
  result(diag_field)
    class (diagYamlObject_type), target, intent(in) :: this      !< The diag_yaml
    integer,                             intent(in) :: yaml_id   !< Yaml id

    type(diagYamlFilesVar_type), pointer :: diag_field !< Diag fields info

    if (yaml_id .eq. DIAG_NOT_REGISTERED) call mpp_error(FATAL, &
      "Diag_manager: The yaml id for this field is not is not set")

    diag_field => this%diag_fields(variable_list%diag_field_indices(yaml_id))

end function get_diag_field_from_id

!> @brief get the diag_fields of a diag_yaml type
!! @return the diag_fields
pure function get_diag_fields(this) &
result(diag_fields)
  class (diagYamlObject_type), intent(in)                 :: this        !< The diag_yaml
  type(diagYamlFilesVar_type), allocatable, dimension (:) :: diag_fields !< Diag fields info

  diag_fields = this%diag_fields
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
    !! The diag_table is not case sensitive (so we are saving it as lowercase)
    file_list%file_name(file_count) = lowercase(trim(diag_yaml%diag_files(file_count)%file_fname)//c_null_char)
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
      !! The diag_table is not case sensitive (so we are saving it as lowercase)
      variable_list%var_name(var_count) = lowercase(variable_list%var_name(var_count))
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
  character(len=:), ALLOCATABLE :: buffer      !< buffer to store any *_units as it is read from the yaml

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "file_name", fileobj%file_fname)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "freq", buffer)
  call parse_key(fileobj%file_fname, buffer, fileobj%file_freq, fileobj%file_frequnit, "freq")
  deallocate(buffer)

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "unlimdim", fileobj%file_unlimdim)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "time_units", buffer)
  call set_file_time_units(fileobj, buffer)
  deallocate(buffer)

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "new_file_freq", buffer, is_optional=.true.)
  call parse_key(fileobj%file_fname, buffer, fileobj%file_new_file_freq, fileobj%file_new_file_freq_units, &
    "new_file_freq")
  deallocate(buffer)

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "filename_time", buffer, is_optional=.true.)
  call set_filename_time(fileobj, buffer)
  deallocate(buffer)

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "start_time", fileobj%file_start_time, is_optional=.true.)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "file_duration", buffer, is_optional=.true.)
  call parse_key(fileobj%file_fname, buffer, fileobj%file_duration, fileobj%file_duration_units, &
    "file_duration")

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

!> @brief This parses the freq, new_file_freq, or file_duration keys which are read in as a comma list
subroutine parse_key(filename, buffer, file_freq, file_frequnit, var)
  character(len=*),         intent(in)     :: filename         !< The name of the file (for error messages)
  character(len=*),         intent(inout)  :: buffer           !< Buffer that was read in from the yaml
  integer,                  intent(out)    :: file_freq(:)     !< buffer to store the freq, new_file_freq, or
                                                               !! file_duration after it is parsed
  integer,                  intent(out)    :: file_frequnit(:) !< buffer to store the freq units, new_file_freq units,
                                                               !! or file_duration units after it is parsed
  character(len=*),         intent(in)     :: var              !< Name of the key parsing

  integer            :: j           !< location of the ",' in the buffer
  integer            :: k           !< location of the " " that seperated the units
  logical            :: finished    !< .true. if the parsing is complete
  integer            :: count       !< Number of keys that have been parsed
  character(len=255) :: str         !< Member of the comma seperated list
  character(len=10)  :: units       !< String to hold the units
  integer            :: err_unit    !< Error key

  if (buffer .eq. "") return

  finished = .false.
  j = 0
  count = 0
  do while (.not. finished)
    count = count + 1
    buffer = buffer(j+1:len_trim(buffer))
    j = index(buffer, ",")
    if (j == 0) then
      !< There is only 1 member in the list
      j = len_trim(buffer)+1
      finished = .true.
    endif

    str = adjustl(buffer(1:j-1))

    k = index(str, " ")
    read(str(1:k-1), *, iostat=err_unit) file_freq(count)
    units = str(k+1:len_trim(str))

    if (err_unit .ne. 0) &
      call mpp_error(FATAL, "Error parsing "//trim(var)//". Check your entry for file"//&
        trim(filename))

    if (file_freq(count) .lt. -1) &
        call mpp_error(FATAL, trim(var)//" is not valid. &
        &Check your entry for file:"//trim(filename))

    if (file_freq(count) .eq. -1 .or. file_freq(count) .eq. 0) then
      !! The file is static so no need to read the units
      file_frequnit(count) = DIAG_DAYS
    else
      if (trim(units) .eq. "") &
        call mpp_error(FATAL, trim(var)//" units is required. &
        &Check your entry for file:"//trim(filename))

      file_frequnit(count) = set_valid_time_units(units, &
        trim(var)//" for file:"//trim(filename))
    endif
  enddo
end subroutine parse_key

!> @brief This checks if the time unit in a diag file is valid and sets the integer equivalent
subroutine set_file_time_units (fileobj, file_timeunit)
  type(diagYamlFiles_type), intent(inout) :: fileobj       !< diagYamlFiles_type obj to checK
  character(len=*),         intent(in)    :: file_timeunit !< file_timeunit as it is read from the diag_table

 fileobj%file_timeunit = set_valid_time_units(file_timeunit, "timeunit for file:"//trim(fileobj%file_fname))
end subroutine set_file_time_units

!> @brief This checks if the filename_time in a diag file is correct and sets the integer equivalent
subroutine set_filename_time(fileobj, filename_time)
  type(diagYamlFiles_type), intent(inout) :: fileobj             !< diagYamlFiles_type obj to check
  character(len=*),         intent(in)    :: filename_time       !< filename_time as it is read from the yaml

  select case (trim(filename_time))
  case ("")
    fileobj%filename_time = middle_time !< This is the default
  case ("begin")
    fileobj%filename_time = begin_time
  case ("middle")
    fileobj%filename_time = middle_time
  case ("end")
    fileobj%filename_time = end_time
  case default
    call mpp_error(FATAL, trim(filename_time)//" is an invalid filename_time &
    &The acceptable values are begin, middle, and end. &
    &Check your entry for file "//trim(fileobj%file_fname))
  end select
end subroutine set_filename_time

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
integer pure function size_file_varlist (this)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 size_file_varlist = size(this%file_varlist)
end function size_file_varlist

!> @brief Inquiry for diag_files_obj%file_fname
!! @return file_fname of a diag_yaml_file obj
pure function get_file_fname (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = this%file_fname
end function get_file_fname
!> @brief Inquiry for diag_files_obj%file_frequnit
!! @return file_frequnit of a diag_yaml_file_obj
pure function get_file_frequnit (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 integer :: res !< What is returned
  res = this%file_frequnit(this%current_new_file_freq_index)
end function get_file_frequnit
!> @brief Inquiry for diag_files_obj%file_freq
!! @return file_freq of a diag_yaml_file_obj
pure function get_file_freq(this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 integer :: res !< What is returned
  res = this%file_freq(this%current_new_file_freq_index)
end function get_file_freq
!> @brief Inquiry for diag_files_obj%file_timeunit
!! @return file_timeunit of a diag_yaml_file_obj
pure function get_file_timeunit (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 integer :: res !< What is returned
  res = this%file_timeunit
end function get_file_timeunit
!> @brief Inquiry for diag_files_obj%file_unlimdim
!! @return file_unlimdim of a diag_yaml_file_obj
pure function get_file_unlimdim(this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = this%file_unlimdim
end function get_file_unlimdim
!> @brief Inquiry for diag_files_obj%file_subregion
!! @return file_sub_region of a diag_yaml_file_obj
function get_file_sub_region (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 type(subRegion_type) :: res !< What is returned
  res = this%file_sub_region
end function get_file_sub_region
!> @brief Inquiry for diag_files_obj%file_new_file_freq
!! @return file_new_file_freq of a diag_yaml_file_obj
pure function get_file_new_file_freq(this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 integer :: res !< What is returned
  res = this%file_new_file_freq(this%current_new_file_freq_index)
end function get_file_new_file_freq
!> @brief Inquiry for diag_files_obj%file_new_file_freq_units
!! @return file_new_file_freq_units of a diag_yaml_file_obj
pure function get_file_new_file_freq_units (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 integer :: res !< What is returned
  res = this%file_new_file_freq_units(this%current_new_file_freq_index)
end function get_file_new_file_freq_units
!> @brief Inquiry for diag_files_obj%file_start_time
!! @return file_start_time of a diag_yaml_file_obj
pure function get_file_start_time (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = this%file_start_time
end function get_file_start_time
!> @brief Inquiry for diag_files_obj%file_duration
!! @return file_duration of a diag_yaml_file_obj
pure function get_file_duration (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 integer :: res !< What is returned
  res = this%file_duration(this%current_new_file_freq_index)
end function get_file_duration
!> @brief Inquiry for diag_files_obj%file_duration_units
!! @return file_duration_units of a diag_yaml_file_obj
pure function get_file_duration_units (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
  integer :: res !< What is returned
  res = this%file_duration_units(this%current_new_file_freq_index)
end function get_file_duration_units
!> @brief Inquiry for diag_files_obj%file_varlist
!! @return file_varlist of a diag_yaml_file_obj
pure function get_file_varlist (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 character (:), allocatable :: res(:) !< What is returned
  res = this%file_varlist
end function get_file_varlist
!> @brief Inquiry for diag_files_obj%file_global_meta
!! @return file_global_meta of a diag_yaml_file_obj
pure function get_file_global_meta (this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 character (len=MAX_STR_LEN), allocatable :: res(:,:) !< What is returned
  res = this%file_global_meta
end function get_file_global_meta
!> @brief Get the integer equivalent of the time to use to determine the filename,
!! if using a wildcard file name (i.e ocn%4yr%2mo%2dy%2hr)
!! @return the integer equivalent of the time to use to determine the filename
pure function get_filename_time(this) &
result (res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 integer :: res !< What is returned
  res = this%filename_time
end function
!> @brief Inquiry for whether file_global_meta is allocated
!! @return Flag indicating if file_global_meta is allocated
function is_global_meta(this) &
  result(res)
 class (diagYamlFiles_type), intent(in) :: this !< The object being inquiried
 logical :: res
 res = .false.
 if (allocated(this%file_global_meta)) &
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
pure function get_var_fname (this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = this%var_fname
end function get_var_fname
!> @brief Inquiry for diag_yaml_files_var_obj%var_varname
!! @return var_varname of a diag_yaml_files_var_obj
pure function get_var_varname (this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = this%var_varname
end function get_var_varname
!> @brief Inquiry for diag_yaml_files_var_obj%var_reduction
!! @return var_reduction of a diag_yaml_files_var_obj
pure function get_var_reduction (this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 integer, allocatable :: res !< What is returned
  res = this%var_reduction
end function get_var_reduction
!> @brief Inquiry for diag_yaml_files_var_obj%var_module
!! @return var_module of a diag_yaml_files_var_obj
pure function get_var_module (this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = this%var_module
end function get_var_module
!> @brief Inquiry for diag_yaml_files_var_obj%var_kind
!! @return var_kind of a diag_yaml_files_var_obj
pure function get_var_kind (this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 integer, allocatable :: res !< What is returned
  res = this%var_kind
end function get_var_kind
!> @brief Inquiry for diag_yaml_files_var_obj%var_outname
!! @return var_outname of a diag_yaml_files_var_obj
pure function get_var_outname (this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned

 if (this%has_var_outname()) then
   res = this%var_outname
 else
   res = this%var_varname !< If outname is not set, the variable name will be used
 endif
end function get_var_outname
!> @brief Inquiry for diag_yaml_files_var_obj%var_longname
!! @return var_longname of a diag_yaml_files_var_obj
pure function get_var_longname (this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = this%var_longname
end function get_var_longname
!> @brief Inquiry for diag_yaml_files_var_obj%var_units
!! @return var_units of a diag_yaml_files_var_obj
pure function get_var_units (this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = this%var_units
end function get_var_units
!> @brief Inquiry for diag_yaml_files_var_obj%var_zbounds
!! @return var_zbounds of a diag_yaml_files_var_obj
pure function get_var_zbounds (this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 real(kind=r4_kind) :: res(2) !< What is returned
  res = this%var_zbounds
end function get_var_zbounds
!> @brief Inquiry for diag_yaml_files_var_obj%var_attributes
!! @return var_attributes of a diag_yaml_files_var_obj
pure function get_var_attributes(this) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 character (len=MAX_STR_LEN), allocatable :: res (:,:) !< What is returned
 res = this%var_attributes
end function get_var_attributes
!> @brief Inquiry for diag_yaml_files_var_obj%n_diurnal
!! @return the number of diurnal samples of a diag_yaml_files_var_obj
pure function get_n_diurnal(this) &
result (res)
  class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
  integer :: res !< What is returned
  res = this%n_diurnal
end function get_n_diurnal
!> @brief Inquiry for diag_yaml_files_var_obj%pow_value
!! @return the pow_value of a diag_yaml_files_var_obj
pure function get_pow_value(this) &
result (res)
  class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
  integer :: res !< What is returned
  res = this%pow_value
end function get_pow_value
!> @brief Inquiry for whether var_attributes is allocated
!! @return Flag indicating if var_attributes is allocated
function is_var_attributes(this) &
result(res)
 class (diagYamlFilesVar_type), intent(in) :: this !< The object being inquiried
 logical :: res
 res = .false.
 if (allocated(this%var_attributes)) &
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

!> @brief Checks if diag_file_obj%file_fname is allocated
!! @return true if diag_file_obj%file_fname is allocated
pure logical function has_file_fname (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_fname = allocated(this%file_fname)
end function has_file_fname
!> @brief Checks if diag_file_obj%file_frequnit is allocated
!! @return true if diag_file_obj%file_frequnit is allocated
pure logical function has_file_frequnit (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_frequnit = this%file_frequnit(this%current_new_file_freq_index) .NE. DIAG_NULL
end function has_file_frequnit
!> @brief diag_file_obj%file_freq is on the stack, so the object always has it
!! @return true if diag_file_obj%file_freq is allocated
pure logical function has_file_freq (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_freq = .true.
end function has_file_freq
!> @brief Checks if diag_file_obj%file_timeunit is allocated
!! @return true if diag_file_obj%file_timeunit is allocated
pure logical function has_file_timeunit (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_timeunit = this%file_timeunit .ne. diag_null
end function has_file_timeunit
!> @brief Checks if diag_file_obj%file_unlimdim is allocated
!! @return true if diag_file_obj%file_unlimdim is allocated
pure logical function has_file_unlimdim (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_unlimdim = allocated(this%file_unlimdim)
end function has_file_unlimdim
!> @brief Checks if diag_file_obj%file_write is on the stack, so this will always be true
!! @return true
pure logical function has_file_write (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_write = .true.
end function has_file_write
!> @brief Checks if diag_file_obj%file_sub_region is being used and has the sub region variables allocated
!! @return true if diag_file_obj%file_sub_region sub region variables are allocated
pure logical function has_file_sub_region (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  if ( this%file_sub_region%grid_type .eq. latlon_gridtype .or. this%file_sub_region%grid_type .eq. index_gridtype) then
       has_file_sub_region = .true.
  else
       has_file_sub_region = .false.
  endif
end function has_file_sub_region
!> @brief diag_file_obj%file_new_file_freq is defined on the stack, so this will return true
!! @return true
pure logical function has_file_new_file_freq (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_new_file_freq = this%file_new_file_freq(this%current_new_file_freq_index) .ne. DIAG_NULL
end function has_file_new_file_freq
!> @brief Checks if diag_file_obj%file_new_file_freq_units is allocated
!! @return true if diag_file_obj%file_new_file_freq_units is allocated
pure logical function has_file_new_file_freq_units (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_new_file_freq_units = this%file_new_file_freq_units(this%current_new_file_freq_index) .ne. diag_null
end function has_file_new_file_freq_units
!> @brief Checks if diag_file_obj%file_start_time is allocated
!! @return true if diag_file_obj%file_start_time is allocated
pure logical function has_file_start_time (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_start_time = allocated(this%file_start_time)
end function has_file_start_time
!> @brief diag_file_obj%file_duration is allocated on th stack, so this is always true
!! @return true
pure logical function has_file_duration (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_duration = this%file_duration(this%current_new_file_freq_index) .ne. DIAG_NULL
end function has_file_duration
!> @brief diag_file_obj%file_duration_units is on the stack, so this will retrun true
!! @return true
pure logical function has_file_duration_units (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_duration_units = this%file_duration_units(this%current_new_file_freq_index) .ne. diag_null
end function has_file_duration_units
!> @brief Checks if diag_file_obj%file_varlist is allocated
!! @return true if diag_file_obj%file_varlist is allocated
pure logical function has_file_varlist (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_varlist = allocated(this%file_varlist)
end function has_file_varlist
!> @brief Checks if diag_file_obj%file_global_meta is allocated
!! @return true if diag_file_obj%file_global_meta is allocated
pure logical function has_file_global_meta (this)
  class(diagYamlFiles_type), intent(in) :: this !< diagYamlFiles_type object to initialize
  has_file_global_meta = allocated(this%file_global_meta)
end function has_file_global_meta

!> @brief Checks if diag_file_obj%var_fname is allocated
!! @return true if diag_file_obj%var_fname is allocated
pure logical function has_var_fname (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_fname = allocated(this%var_fname)
end function has_var_fname
!> @brief Checks if diag_file_obj%var_varname is allocated
!! @return true if diag_file_obj%var_varname is allocated
pure logical function has_var_varname (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_varname = allocated(this%var_varname)
end function has_var_varname
!> @brief Checks if diag_file_obj%var_reduction is allocated
!! @return true if diag_file_obj%var_reduction is allocated
pure logical function has_var_reduction (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_reduction = allocated(this%var_reduction)
end function has_var_reduction
!> @brief Checks if diag_file_obj%var_module is allocated
!! @return true if diag_file_obj%var_module is allocated
pure logical function has_var_module (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_module = allocated(this%var_module)
end function has_var_module
!> @brief Checks if diag_file_obj%var_kind is allocated
!! @return true if diag_file_obj%var_kind is allocated
pure logical function has_var_kind (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_kind = allocated(this%var_kind)
end function has_var_kind
!> @brief diag_file_obj%var_write is on the stack, so this returns true
!! @return true
pure logical function has_var_write (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_write = .true.
end function has_var_write
!> @brief Checks if diag_file_obj%var_outname is allocated
!! @return true if diag_file_obj%var_outname is allocated
pure logical function has_var_outname (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  if (allocated(this%var_outname)) then
    if (trim(this%var_outname) .ne. "") then
      has_var_outname = .true.
    else
      has_var_outname = .false.
    endif
  else
    has_var_outname = .true.
  endif
end function has_var_outname
!> @brief Checks if diag_file_obj%var_longname is allocated
!! @return true if diag_file_obj%var_longname is allocated
pure logical function has_var_longname (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_longname = allocated(this%var_longname)
end function has_var_longname
!> @brief Checks if diag_file_obj%var_units is allocated
!! @return true if diag_file_obj%var_units is allocated
pure logical function has_var_units (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_units = allocated(this%var_units)
end function has_var_units
!> @brief Checks if diag_file_obj%var_zbounds is allocated
!! @return true if diag_file_obj%var_zbounds is allocated
pure logical function has_var_zbounds (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_zbounds = any(this%var_zbounds .ne. diag_null)
end function has_var_zbounds
!> @brief Checks if diag_file_obj%var_attributes is allocated
!! @return true if diag_file_obj%var_attributes is allocated
pure logical function has_var_attributes (this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to initialize
  has_var_attributes = allocated(this%var_attributes)
end function has_var_attributes
!> @brief Checks if diag_file_obj%n_diurnal is set
!! @return true if diag_file_obj%n_diurnal is set
pure logical function has_n_diurnal(this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to inquire
  has_n_diurnal = (this%n_diurnal .ne. 0)
end function has_n_diurnal
!> @brief Checks if diag_file_obj%pow_value is set
!! @return true if diag_file_obj%pow_value is set
pure logical function has_pow_value(this)
  class(diagYamlFilesVar_type), intent(in) :: this !< diagYamlvar_type object to inquire
  has_pow_value = (this%pow_value .ne. 0)
end function has_pow_value

!> @brief Checks if diag_file_obj%diag_title is allocated
!! @return true if diag_file_obj%diag_title is allocated
pure logical function has_diag_title (this)
  class(diagYamlObject_type), intent(in) :: this !< diagYamlObject_type object to inquire
  has_diag_title = allocated(this%diag_title)
end function has_diag_title
!> @brief diag_file_obj%diag_basedate is on the stack, so this is always true
!! @return true
pure logical function has_diag_basedate (this)
  class(diagYamlObject_type), intent(in) :: this !< diagYamlObject_type object to initialize
  has_diag_basedate = .true.
end function has_diag_basedate
!> @brief Checks if diag_file_obj%diag_files is allocated
!! @return true if diag_file_obj%diag_files is allocated
pure logical function has_diag_files (this)
  class(diagYamlObject_type), intent(in) :: this !< diagYamlObject_type object to initialize
  has_diag_files = allocated(this%diag_files)
end function has_diag_files
!> @brief Checks if diag_file_obj%diag_fields is allocated
!! @return true if diag_file_obj%diag_fields is allocated
pure logical function has_diag_fields (this)
  class(diagYamlObject_type), intent(in) :: this !< diagYamlObject_type object to initialize
  has_diag_fields = allocated(this%diag_fields)
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
                               & lowercase(trim(diag_field_name))//":"//lowercase(trim(module_name)//c_null_char))
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

!> @brief Gets field indices corresponding to the indices (input argument) in the sorted variable_list
!! @return Copy of array of field indices
function get_diag_field_ids(indices) result(field_ids)

  integer, intent(in) :: indices(:) !< Indices of the fields in the sorted variable_list array
  integer, allocatable :: field_ids(:)
  integer :: i !< For do loop

  allocate(field_ids(size(indices)))

  do i = 1, size(indices)
    field_ids(i) = variable_list%diag_field_indices(indices(i))
  end do

end function get_diag_field_ids

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

#endif
end module fms_diag_yaml_mod
!> @}
! close documentation grouping
