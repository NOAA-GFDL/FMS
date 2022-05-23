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
                           index_gridtype, null_gridtype
use yaml_parser_mod, only: open_and_parse_file, get_value_from_key, get_num_blocks, get_nkeys, &
                           get_block_ids, get_key_value, get_key_ids, get_key_name
use mpp_mod,         only: mpp_error, FATAL
use, intrinsic :: iso_c_binding, only : c_ptr, c_null_char
use fms_string_utils_mod, only: fms_array_to_pointer, fms_find_my_string, fms_sort_this, fms_find_unique
use platform_mod, only: r4_kind, i4_kind

implicit none

private

public :: diag_yaml_object_init, diag_yaml_object_end
public :: diagYamlObject_type, get_diag_yaml_obj, get_title, get_basedate, get_diag_files, get_diag_fields
public :: diagYamlFiles_type, diagYamlFilesVar_type
public :: get_num_unique_fields, find_diag_field, get_diag_fields_entries, get_diag_files_entries

!> @}

integer, parameter :: basedate_size = 6
integer, parameter :: NUM_SUB_REGION_ARRAY = 8
integer, parameter :: MAX_STR_LEN = 255

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
  integer                        :: zbounds(2)  !< indices of the z axis limits (zbegin, zend)
  integer                        :: tile        !< Tile number of the sub region
                                                !! required if using the "index" grid type

end type subRegion_type

!> @brief type to hold the diag_file information
type diagYamlFiles_type
  character (len=:), private, allocatable :: file_fname !< file name
  character (len=:), private, allocatable :: file_frequnit !< the frequency unit
  integer, private    :: file_freq !< the frequency of data
  character (len=:), private, allocatable :: file_timeunit !< The unit of time
  character (len=:), private, allocatable :: file_unlimdim !< The name of the unlimited dimension
  type(subRegion_type), private :: file_sub_region !< type containing info about the subregion, if any
  integer, private :: file_new_file_freq !< Frequency for closing the existing file
  character (len=:), private, allocatable :: file_new_file_freq_units !< Time units for creating a new file.
                                                                      !! Required if “new_file_freq” used
  character (len=:), private, allocatable :: file_start_time !< Time to start the file for the first time. Requires
                                                             !! “new_file_freq”
  integer, private :: file_duration !< How long the file should receive data after start time
                                    !! in “file_duration_units”.  This optional field can only
                                    !! be used if the start_time field is present.  If this field
                                    !! is absent, then the file duration will be equal to the
                                    !! frequency for creating new files.
                                    !! NOTE: The file_duration_units field must also be present if
                                    !! this field is present.
  character (len=:), private, allocatable :: file_duration_units !< The file duration units
  !< Need to use `MAX_STR_LEN` because not all filenames/global attributes are the same length
  character (len=MAX_STR_LEN), dimension(:), private, allocatable :: file_varlist !< An array of variable names
                                                                                  !! within a file
  character (len=MAX_STR_LEN), dimension(:,:), private, allocatable :: file_global_meta !< Array of key(dim=1)
                                                                                        !! and values(dim=2) to be
                                                                                        !! added as global meta data to
                                                                                        !! the file

 contains

 !> All getter functions (functions named get_x(), for member field named x)
 !! return copies of the member variables unless explicitly noted.
 procedure :: get_file_fname
 procedure :: get_file_frequnit
 procedure :: get_file_freq
 procedure :: get_file_timeunit
 procedure :: get_file_unlimdim
 procedure :: get_file_sub_region
 procedure :: get_file_new_file_freq
 procedure :: get_file_new_file_freq_units
 procedure :: get_file_start_time
 procedure :: get_file_duration
 procedure :: get_file_duration_units
 procedure :: get_file_varlist
 procedure :: get_file_global_meta
 procedure :: is_global_meta
 !> Has functions to determine if allocatable variables are true.  If a variable is not an allocatable
 !! then is will always return .true.
 procedure :: has_file_fname 
 procedure :: has_file_frequnit 
 procedure :: has_file_freq 
 procedure :: has_file_timeunit 
 procedure :: has_file_unlimdim 
 procedure :: has_file_sub_region 
 procedure :: has_file_new_file_freq 
 procedure :: has_file_new_file_freq_units 
 procedure :: has_file_start_time 
 procedure :: has_file_duration 
 procedure :: has_file_duration_units 
 procedure :: has_file_varlist 
 procedure :: has_file_global_meta 

end type diagYamlFiles_type

!> @brief type to hold the info a diag_field
type diagYamlFilesVar_type
  character (len=:), private, allocatable :: var_fname !< The field/diagnostic name
  character (len=:), private, allocatable :: var_varname !< The name of the variable
  character (len=:), private, allocatable :: var_reduction !< Reduction to be done on var
  character (len=:), private, allocatable :: var_module !< The module that th variable is in
  character (len=:), private, allocatable :: var_skind !< The type/kind of the variable
  character (len=:), private, allocatable :: var_outname !< Name of the variable as written to the file
  character (len=:), private, allocatable :: var_longname !< Overwrites the long name of the variable
  character (len=:), private, allocatable :: var_units !< Overwrites the units
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
  procedure :: get_var_skind
  procedure :: get_var_outname
  procedure :: get_var_longname
  procedure :: get_var_units
  procedure :: get_var_attributes
  procedure :: get_n_diurnal
  procedure :: get_pow_value
  procedure :: is_var_attributes

  procedure :: has_var_fname 
  procedure :: has_var_varname 
  procedure :: has_var_reduction 
  procedure :: has_var_module 
  procedure :: has_var_skind 
  procedure :: has_var_outname 
  procedure :: has_var_longname 
  procedure :: has_var_units 
  procedure :: has_var_attributes 
  procedure :: has_n_diurnal
  procedure :: has_pow_value

end type diagYamlFilesVar_type

!> @brief Object that holds the information of the diag_yaml
!> @ingroup fms_diag_yaml_mod
type diagYamlObject_type
  character(len=:), allocatable, private :: diag_title                   !< Experiment name
  integer, private, dimension (basedate_size) :: diag_basedate           !< basedate array
  type(diagYamlFiles_type), allocatable, private, dimension (:) :: diag_files!< History file info
  type(diagYamlFilesVar_type), allocatable, private, dimension (:) :: diag_fields !< Diag fields info
  contains
  procedure :: get_title        !< Returns the title
  procedure :: get_basedate     !< Returns the basedate array
  procedure :: get_diag_files   !< Returns the diag_files array
  procedure :: get_diag_fields  !< Returns the diag_field array

  procedure :: has_diag_title                   
  procedure :: has_diag_basedate           
  procedure :: has_diag_files
  procedure :: has_diag_fields 

end type diagYamlObject_type

type (diagYamlObject_type) :: diag_yaml  !< Obj containing the contents of the diag_table.yaml
type (varList_type), save :: variable_list !< List of all the variables in the diag_table.yaml
type (fileList_type), save :: file_list !< List of all files in the diag_table.yaml

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

      !> Save the variable name in the variable_list
      variable_list%var_name(var_count) = trim(diag_yaml%diag_fields(var_count)%var_varname)//c_null_char
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

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "file_name", fileobj%file_fname)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "freq_units", fileobj%file_frequnit)
  call get_value_from_key(diag_yaml_id, diag_file_id, "freq", fileobj%file_freq)
  call check_file_freq(fileobj)

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "unlimdim", fileobj%file_unlimdim)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "time_units", fileobj%file_timeunit)
  call check_file_time_units(fileobj)

  call get_value_from_key(diag_yaml_id, diag_file_id, "new_file_freq", fileobj%file_new_file_freq, is_optional=.true.)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "new_file_freq_units", fileobj%file_new_file_freq_units, &
         is_optional=.true.)
  call check_new_file_freq(fileobj)

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "start_time", fileobj%file_start_time, is_optional=.true.)
  call get_value_from_key(diag_yaml_id, diag_file_id, "file_duration", fileobj%file_duration, is_optional=.true.)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "file_duration_units", fileobj%file_duration_units, &
    is_optional=.true.)
  call check_file_duration(fileobj)

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
    call mpp_error(FATAL, "diag_yaml_object_init: file "//trim(fileobj%file_fname)//" has multiple global_meta blocks")
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

  call diag_get_value_from_key(diag_file_id, var_id, "var_name", field%var_varname)
  call diag_get_value_from_key(diag_file_id, var_id, "reduction", field%var_reduction)
  call check_field_reduction(field)

  call diag_get_value_from_key(diag_file_id, var_id, "module", field%var_module)
  call diag_get_value_from_key(diag_file_id, var_id, "kind", field%var_skind)
  call check_field_kind(field)

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

  sub_region%zbounds = DIAG_NULL
  call get_value_from_key(diag_yaml_id, sub_region_id, "zbounds", sub_region%zbounds, is_optional=.true.)

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

!> @brief This checks if the file frequency in a diag file is valid and crashes if it isn't
subroutine check_file_freq(fileobj)
  type(diagYamlFiles_type), intent(inout) :: fileobj      !< diagYamlFiles_type obj to check

  if (.not. (fileobj%file_freq >= -1) ) &
    call mpp_error(FATAL, "freq must be greater than or equal to -1. &
      &Check you entry for"//trim(fileobj%file_fname))
  if(.not. is_valid_time_units(fileobj%file_frequnit)) &
    call mpp_error(FATAL, trim(fileobj%file_frequnit)//" is not a valid file_frequnit. &
      &The acceptable values are seconds, minuts, hours, days, months, years. &
      &Check your entry for file:"//trim(fileobj%file_fname))
end subroutine check_file_freq

!> @brief This checks if the time unit in a diag file is valid and crashes if it isn't
subroutine check_file_time_units (fileobj)
  type(diagYamlFiles_type), intent(inout) :: fileobj      !< diagYamlFiles_type obj to checK

  if(.not. is_valid_time_units(fileobj%file_timeunit)) &
    call mpp_error(FATAL, trim(fileobj%file_timeunit)//" is not a valid time_unit. &
      &The acceptable values are seconds, minuts, hours, days, months, years. &
      &Check your entry for file:"//trim(fileobj%file_fname))
end subroutine check_file_time_units

!> @brief This checks if the new file frequency in a diag file is valid and crashes if it isn't
subroutine check_new_file_freq(fileobj)
  type(diagYamlFiles_type), intent(inout) :: fileobj      !< diagYamlFiles_type obj to check

  if (fileobj%file_new_file_freq > 0) then
    if (trim(fileobj%file_new_file_freq_units) .eq. "") &
      call mpp_error(FATAL, "new_file_freq_units is required if using new_file_freq. &
        &Check your entry for file:"//trim(fileobj%file_fname))

    if (.not. is_valid_time_units(fileobj%file_new_file_freq_units)) &
      call mpp_error(FATAL, trim(fileobj%file_new_file_freq_units)//" is not a valid new_file_freq_units. &
        &The acceptable values are seconds, minuts, hours, days, months, years. &
        &Check your entry for file:"//trim(fileobj%file_fname))
  endif
end subroutine check_new_file_freq

!> @brief This checks if the file duration in a diag file is valid and crashes if it isn't
subroutine check_file_duration(fileobj)
  type(diagYamlFiles_type), intent(inout) :: fileobj      !< diagYamlFiles_type obj to check

  if (fileobj%file_duration > 0) then
    if(trim(fileobj%file_duration_units) .eq. "") &
      call mpp_error(FATAL, "file_duration_units is required if using file_duration. &
        &Check your entry for file:"//trim(fileobj%file_fname))

    if (.not. is_valid_time_units(fileobj%file_duration_units)) &
      call mpp_error(FATAL, trim(fileobj%file_duration_units)//" is not a valid file_duration_units. &
        &The acceptable values are seconds, minuts, hours, days, months, years. &
        &Check your entry for file:"//trim(fileobj%file_duration_units))
  endif
end subroutine check_file_duration

!> @brief This checks if the kind of a diag field is valid and crashes if it isn't
subroutine check_field_kind(field)
  type(diagYamlFilesVar_type), intent(in) :: field        !< diagYamlFilesVar_type obj to read the contents into

  select case (TRIM(field%var_skind))
  case ("r4", "r8", "i4", "i8")
  case default
    call mpp_error(FATAL, trim(field%var_skind)//" is an invalid kind! &
      &The acceptable values are r4, r8, i4, i8. &
      &Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
  end select

end subroutine check_field_kind

!> @brief This checks if the reduction of a diag field is valid and crashes if it isn't
!! If the reduction method is diurnalXX or powXX, it gets the number of diurnal sample and the power value
subroutine check_field_reduction(field)
  type(diagYamlFilesVar_type), intent(inout) :: field        !< diagYamlFilesVar_type obj to read the contents into

  integer :: n_diurnal !< number of diurnal samples
  integer :: pow_value !< The power value
  integer :: ioerror   !< io error status after reading in the diurnal samples

  n_diurnal = 0
  pow_value = 0
  ioerror = 0
  if (index(field%var_reduction, "diurnal") .ne. 0) then
    READ (UNIT=field%var_reduction(8:LEN_TRIM(field%var_reduction)), FMT=*, IOSTAT=ioerror) n_diurnal
    if (ioerror .ne. 0) &
      call mpp_error(FATAL, "Error getting the number of diurnal samples from "//trim(field%var_reduction))
    if (n_diurnal .le. 0) &
      call mpp_error(FATAL, "Diurnal samples should be greater than 0. &
        & Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
  elseif (index(field%var_reduction, "pow") .ne. 0) then
    READ (UNIT=field%var_reduction(4:LEN_TRIM(field%var_reduction)), FMT=*, IOSTAT=ioerror) pow_value
    if (ioerror .ne. 0) &
      call mpp_error(FATAL, "Error getting the power value from "//trim(field%var_reduction))
      if (pow_value .le. 0) &
      call mpp_error(FATAL, "The power value should be greater than 0. &
        & Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
  else
    select case (TRIM(field%var_reduction))
    case ("none", "average", "min", "max", "rms")
    case default
      call mpp_error(FATAL, trim(field%var_reduction)//" is an invalid reduction method! &
        &The acceptable values are none, average, pow##, diurnal##, min, max, and rms. &
        &Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
    end select
  endif

  field%n_diurnal = n_diurnal
  field%pow_value = pow_value
end subroutine check_field_reduction

!> @brief This checks if a time unit is valid
!! @return Flag indicating if the time units are valid
pure function is_valid_time_units(time_units) &
result(is_valid)
  character(len=*), intent(in) :: time_units
  logical :: is_valid

  select case (TRIM(time_units))
  case ("seconds", "minutes", "hours", "days", "months", "years")
    is_valid = .true.
  case default
    is_valid = .false.
  end select
end function is_valid_time_units

!!!!!!! YAML FILE INQUIRIES !!!!!!!
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
 character (len=:), allocatable :: res !< What is returned
  res = diag_files_obj%file_frequnit
end function get_file_frequnit
!> @brief Inquiry for diag_files_obj%file_freq
!! @return file_freq of a diag_yaml_file_obj
pure function get_file_freq(diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_freq
end function get_file_freq
!> @brief Inquiry for diag_files_obj%file_timeunit
!! @return file_timeunit of a diag_yaml_file_obj
pure function get_file_timeunit (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
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
  res = diag_files_obj%file_new_file_freq
end function get_file_new_file_freq
!> @brief Inquiry for diag_files_obj%file_new_file_freq_units
!! @return file_new_file_freq_units of a diag_yaml_file_obj
pure function get_file_new_file_freq_units (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (:), allocatable :: res !< What is returned
  res = diag_files_obj%file_new_file_freq_units
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
  res = diag_files_obj%file_duration
end function get_file_duration
!> @brief Inquiry for diag_files_obj%file_duration_units
!! @return file_duration_units of a diag_yaml_file_obj
pure function get_file_duration_units (diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
  character (:), allocatable :: res !< What is returned
  res = diag_files_obj%file_duration_units
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
 character (len=:), allocatable :: res !< What is returned
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
!> @brief Inquiry for diag_yaml_files_var_obj%var_skind
!! @return var_skind of a diag_yaml_files_var_obj
pure function get_var_skind (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_skind
end function get_var_skind
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
  obj%file_duration       = DIAG_NULL
  obj%file_new_file_freq  = DIAG_NULL
  obj%file_sub_region%tile = DIAG_NULL
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
  has_file_frequnit = allocated(obj%file_frequnit)
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
  has_file_timeunit = allocated(obj%file_timeunit)
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
  has_file_new_file_freq = .true.
end function has_file_new_file_freq
!> @brief Checks if obj%file_new_file_freq_units is allocated
!! @return true if obj%file_new_file_freq_units is allocated
pure logical function has_file_new_file_freq_units (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_new_file_freq_units = allocated(obj%file_new_file_freq_units)
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
  has_file_duration = .true.
end function has_file_duration
!> @brief obj%file_duration_units is on the stack, so this will retrun true
!! @return true 
pure logical function has_file_duration_units (obj)
  class(diagYamlFiles_type), intent(in) :: obj !< diagYamlFiles_type object to initialize
  has_file_duration_units = .true.
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
!> @brief Checks if obj%var_skind is allocated
!! @return true if obj%var_skind is allocated
pure logical function has_var_skind (obj)
  class(diagYamlFilesVar_type), intent(in) :: obj !< diagYamlvar_type object to initialize
  has_var_skind = allocated(obj%var_skind)
end function has_var_skind
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
function find_diag_field(diag_field_name) &
result(indices)

  character(len=*), intent(in) :: diag_field_name !< diag_field name to search for

  integer, allocatable :: indices(:)

  indices = fms_find_my_string(variable_list%var_pointer, size(variable_list%var_pointer), &
                               & diag_field_name//c_null_char)
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

!> @brief Gets the diag_files entries corresponding to the indices of the sorted variable_list
!! @return Array of diag_files
function get_diag_files_entries(indices) &
  result(diag_file)

  integer, intent(in) :: indices(:) !< Indices of the field in the sorted variable_list
  type(diagYamlFiles_type), dimension (:), allocatable :: diag_file

  integer :: i !< For do loops
  integer :: field_id !< Indices of the field in the diag_yaml array
  integer :: file_id !< Indices of the file in the diag_yaml array
  character(len=120) :: filename !< Filename of the field
  integer, allocatable :: file_indices(:) !< Indices of the file in the sorted variable_list

  allocate(diag_file(size(indices)))

  do i = 1, size(indices)
    field_id = variable_list%diag_field_indices(indices(i))
    filename = diag_yaml%diag_fields(field_id)%var_fname

    file_indices = fms_find_my_string(file_list%file_pointer, size(file_list%file_pointer), &
      & trim(filename)//c_null_char)

    if (size(file_indices) .ne. 1) &
      & call mpp_error(FATAL, "get_diag_files_entries: Error getting the correct number of file indices!")

    if (file_indices(1) .eq. diag_null) &
      & call mpp_error(FATAL, "get_diag_files_entries: Error finding the filename in the diag_files")

    file_id = file_list%diag_file_indices(file_indices(1))
    diag_file(i) = diag_yaml%diag_files(file_id)
  end do

end function get_diag_files_entries
#endif
end module fms_diag_yaml_mod
!> @}
! close documentation grouping
