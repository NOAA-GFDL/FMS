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
use diag_data_mod, only: DIAG_NULL
use yaml_parser_mod
use mpp_mod

implicit none

private

public :: diag_yaml_object_init, diag_yaml_object_end
public :: diagYamlObject_type, get_diag_yaml_obj, get_title, get_basedate, get_diag_files, get_diag_fields
public :: diagYamlFiles_type, diagYamlFilesVar_type
!> @}

integer, parameter :: basedate_size = 6
integer, parameter :: NUM_SUB_REGION_ARRAY = 8
integer, parameter :: MAX_STR_LEN = 255

!> @brief type to hold the sub region information about a file
type subRegion_type
  character (len=:), allocatable :: grid_type !< Flag indicating the type of region,
                                              !! acceptable values are "latlon" and "index"
  real, allocatable :: lat_lon_sub_region (:) !< Array that stores the grid point bounds for the sub region
                                              !! to use if grid_type is set to "latlon"
                                              !! [dim1_begin, dim1_end, dim2_begin, dim2_end,
                                              !!  dim3_begin, dim3_end, dim4_begin, dim4_end]
  integer, allocatable :: index_sub_region (:) !< Array that stores the index bounds for the sub region to
                                               !! to use if grid_type is set to "index"
                                               !! [dim1_begin, dim1_end, dim2_begin, dim2_end,
                                               !!  dim3_begin, dim3_end, dim4_begin, dim4_end]
  integer :: tile !< Tile number of the sub region, required if using the "index" grid type

end type subRegion_type

!> @brief type to hold the diag_file information
type diagYamlFiles_type
  character (len=:), private, allocatable :: file_fname !< file name
  character (len=:), private, allocatable :: file_frequnit !< the frequency unit
  integer, private    :: file_freq !< the frequency of data
  character (len=:), private, allocatable :: file_timeunit !< The unit of time
  character (len=:), private, allocatable :: file_unlimdim !< The name of the unlimited dimension
  logical, private :: file_write !< false if the user doesn't want to the file to be created
  character (len=:), private, allocatable :: string_file_write !< false if the user doesn’t want the file to be
                                                               !! created (default is true).
  character (len=:), private, allocatable :: file_realm !< The modeling realm that the variables come from
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
                                                                                        !! and values(dim=2) to be added as global
                                                                                        !! meta data to the file

 contains
 procedure :: get_file_fname
 procedure :: get_file_frequnit
 procedure :: get_file_freq
 procedure :: get_file_timeunit
 procedure :: get_file_unlimdim
 procedure :: get_file_write
 procedure :: get_file_realm
 procedure :: get_file_sub_region
 procedure :: get_file_new_file_freq
 procedure :: get_file_new_file_freq_units
 procedure :: get_file_start_time
 procedure :: get_file_duration
 procedure :: get_file_duration_units
 procedure :: get_file_varlist
 procedure :: get_file_global_meta
 procedure :: is_global_meta

end type diagYamlFiles_type

!> @brief type to hold the info a diag_field
type diagYamlFilesVar_type
  character (len=:), private, allocatable :: var_fname !< The field/diagnostic name
  character (len=:), private, allocatable :: var_varname !< The name of the variable
  character (len=:), private, allocatable :: var_reduction !< Reduction to be done on var
  character (len=:), private, allocatable :: var_module !< The module that th variable is in
  character (len=:), private, allocatable :: var_skind !< The type/kind of the variable
  character (len=:), private, allocatable :: string_var_write !< false if the user doesn’t want the variable to be
                                                              !! written to the file (default: true).
  logical, private :: var_write !< false if the user doesn’t want the variable to be
                                !! written to the file (default: true).
  character (len=:), private, allocatable :: var_outname !< Name of the variable as written to the file
  character (len=:), private, allocatable :: var_longname !< Overwrites the long name of the variable
  character (len=:), private, allocatable :: var_units !< Overwrites the units

  !< Need to use `MAX_STR_LEN` because not all filenames/global attributes are the same length
  character (len=MAX_STR_LEN), dimension (:, :), private, allocatable :: var_attributes !< Attributes to overwrite or
                                                                                        !! add from diag_yaml
 contains
  procedure :: get_var_fname
  procedure :: get_var_varname
  procedure :: get_var_reduction
  procedure :: get_var_module
  procedure :: get_var_skind
  procedure :: get_var_outname
  procedure :: get_var_longname
  procedure :: get_var_units
  procedure :: get_var_write
  procedure :: get_var_attributes
  procedure :: is_var_attributes
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
end type diagYamlObject_type

type (diagYamlObject_type) :: diag_yaml  !< Obj containing the contents of the diag_table.yaml

!> @addtogroup fms_diag_yaml_mod
!> @{
contains

!> @brief gets the diag_yaml module variable
!! @return a copy of the diag_yaml module variable
function get_diag_yaml_obj() &
result(res)
  type (diagYamlObject_type) :: res

   res = diag_yaml
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
pure function get_diag_files(diag_yaml) &
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
subroutine diag_yaml_object_init
  integer              :: diag_yaml_id     !< Id for the diag_table yaml
  integer              :: nfiles           !< Number of files in the diag_table yaml
  integer, allocatable :: diag_file_ids(:) !< Ids of the files in the diag_table yaml
  integer              :: i, j             !< For do loops
  integer              :: total_nvars      !< The total number of variables in the diag_table yaml
  integer              :: var_count        !< The current number of variables added to the diag_yaml obj
  integer              :: nvars            !< The number of variables in the current file
  integer, allocatable :: var_ids(:)       !< Ids of the variables in diag_table yaml

  diag_yaml_id = open_and_parse_file("diag_table.yaml")

  call diag_get_value_from_key(diag_yaml_id, 0, "title", diag_yaml%diag_title)
  call get_value_from_key(diag_yaml_id, 0, "base_date", diag_yaml%diag_basedate)

  nfiles = get_num_blocks(diag_yaml_id, "diag_files")
  allocate(diag_yaml%diag_files(nfiles))
  allocate(diag_file_ids(nfiles))
  call get_block_ids(diag_yaml_id, "diag_files", diag_file_ids)

  total_nvars = get_total_num_vars(diag_yaml_id, diag_file_ids)
  allocate(diag_yaml%diag_fields(total_nvars))

  var_count = 0
  nfiles_loop: do i = 1, nfiles
    call diag_yaml_files_obj_init(diag_yaml%diag_files(i))
    call fill_in_diag_files(diag_yaml_id, diag_file_ids(i), diag_yaml%diag_files(i))

    nvars = 0
    nvars = get_num_blocks(diag_yaml_id, "varlist", parent_block_id=diag_file_ids(i))
    allocate(var_ids(nvars))
    call get_block_ids(diag_yaml_id, "varlist", var_ids, parent_block_id=diag_file_ids(i))
    allocate(diag_yaml%diag_files(i)%file_varlist(nvars))
    nvars_loop: do j = 1, nvars
      var_count = var_count + 1
      !> Save the filename in the diag_field type
      diag_yaml%diag_fields(var_count)%var_fname = diag_yaml%diag_files(i)%file_fname

      call fill_in_diag_fields(diag_yaml_id, var_ids(j), diag_yaml%diag_fields(var_count))

      !> Save the variable name in the diag_file type
      diag_yaml%diag_files(i)%file_varlist(j) = diag_yaml%diag_fields(var_count)%var_varname
    enddo nvars_loop
    deallocate(var_ids)
  enddo nfiles_loop

  deallocate(diag_file_ids)
end subroutine

!> @brief Destroys the diag_yaml object
subroutine diag_yaml_object_end()
  integer :: i !< For do loops

  do i = 1, size(diag_yaml%diag_files, 1)
    if(allocated(diag_yaml%diag_files(i)%file_varlist)) deallocate(diag_yaml%diag_files(i)%file_varlist)
    if(allocated(diag_yaml%diag_files(i)%file_global_meta)) deallocate(diag_yaml%diag_files(i)%file_global_meta)
    if(allocated(diag_yaml%diag_files(i)%file_sub_region%lat_lon_sub_region)) &
      deallocate(diag_yaml%diag_files(i)%file_sub_region%lat_lon_sub_region)
    if(allocated(diag_yaml%diag_files(i)%file_sub_region%index_sub_region)) &
      deallocate(diag_yaml%diag_files(i)%file_sub_region%index_sub_region)
  enddo
  if(allocated(diag_yaml%diag_files)) deallocate(diag_yaml%diag_files)

  do i = 1, size(diag_yaml%diag_fields, 1)
    if(allocated(diag_yaml%diag_fields(i)%var_attributes)) deallocate(diag_yaml%diag_fields(i)%var_attributes)
  enddo
  if(allocated(diag_yaml%diag_fields)) deallocate(diag_yaml%diag_fields)

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

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "file_name", fileobj%file_fname)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "freq_units", fileobj%file_frequnit)
  call get_value_from_key(diag_yaml_id, diag_file_id, "freq", fileobj%file_freq)
  call check_file_freq(fileobj)

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "unlimdim", fileobj%file_unlimdim)
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "time_units", fileobj%file_timeunit)
  call check_file_time_units(fileobj)

  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "write_file", fileobj%string_file_write, is_optional=.true.)
  if (fileobj%string_file_write .eq. "false") fileobj%file_write = .false.
  call diag_get_value_from_key(diag_yaml_id, diag_file_id, "realm", fileobj%file_realm, is_optional=.true.)
  call check_file_realm(fileobj)

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
    call diag_get_value_from_key(diag_yaml_id, sub_region_id(1), "grid_type", fileobj%file_sub_region%grid_type)
    if (trim(fileobj%file_sub_region%grid_type) .eq. "latlon") then
      allocate(fileobj%file_sub_region%lat_lon_sub_region(8))
      fileobj%file_sub_region%lat_lon_sub_region = DIAG_NULL
      call get_sub_region(diag_yaml_id, sub_region_id(1), fileobj%file_sub_region%lat_lon_sub_region)
    elseif (trim(fileobj%file_sub_region%grid_type) .eq. "index") then
      allocate(fileobj%file_sub_region%index_sub_region(8))
      fileobj%file_sub_region%index_sub_region = DIAG_NULL
      call get_sub_region(diag_yaml_id, sub_region_id(1), fileobj%file_sub_region%index_sub_region)
      call get_value_from_key(diag_yaml_id, sub_region_id(1), "tile", fileobj%file_sub_region%tile, is_optional=.true.)
      if (fileobj%file_sub_region%tile .eq. DIAG_NULL) call mpp_error(FATAL, "The tile number is required when defining a "//&
        "subregion. Check your subregion entry for "//trim(fileobj%file_fname))
    else
      call mpp_error(FATAL, trim(fileobj%file_sub_region%grid_type)//" is not a valid region type. &
      &The acceptable values are latlon and index. &
      &Check your entry for file:"//trim(fileobj%file_fname))
    endif
  elseif (nsubregion .ne. 0) then
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

  field%var_write = .true.
  call diag_get_value_from_key(diag_file_id, var_id, "var_name", field%var_varname)
  call diag_get_value_from_key(diag_file_id, var_id, "reduction", field%var_reduction)
  call check_field_reduction(field)

  call diag_get_value_from_key(diag_file_id, var_id, "module", field%var_module)
  call diag_get_value_from_key(diag_file_id, var_id, "kind", field%var_skind)
  call check_field_kind(field)

  call diag_get_value_from_key(diag_file_id, var_id, "write_var", field%string_var_write, is_optional=.true.)
  if (trim(field%string_var_write) .eq. "false") field%var_write = .false.

  call diag_get_value_from_key(diag_file_id, var_id, "output_name", field%var_outname)
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
      call mpp_error(FATAL, "diag_yaml_object_init: variable "//trim(field%var_varname)//" has multiple attribute blocks")
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
subroutine get_sub_region(diag_yaml_id, sub_region_id, sub_region)
  integer, intent(in)  :: diag_yaml_id       !< Id of the diag_table yaml file
  integer, intent(in)  :: sub_region_id      !< Id of the region block to read from
  class(*),intent(out) :: sub_region (NUM_SUB_REGION_ARRAY) !< Array storing the bounds of the sub region

  call get_value_from_key(diag_yaml_id, sub_region_id, "dim1_begin", sub_region(1), is_optional=.true.)
  call get_value_from_key(diag_yaml_id, sub_region_id, "dim1_end", sub_region(2), is_optional=.true.)
  call get_value_from_key(diag_yaml_id, sub_region_id, "dim2_begin", sub_region(3), is_optional=.true.)
  call get_value_from_key(diag_yaml_id, sub_region_id, "dim2_end", sub_region(4), is_optional=.true.)
  call get_value_from_key(diag_yaml_id, sub_region_id, "dim3_begin", sub_region(5), is_optional=.true.)
  call get_value_from_key(diag_yaml_id, sub_region_id, "dim3_end", sub_region(6), is_optional=.true.)
  call get_value_from_key(diag_yaml_id, sub_region_id, "dim4_begin", sub_region(7), is_optional=.true.)
  call get_value_from_key(diag_yaml_id, sub_region_id, "dim4_end", sub_region(8), is_optional=.true.)

end subroutine get_sub_region

!> @brief gets the total number of variables in the diag_table yaml file
!! @return total number of variables
function get_total_num_vars(diag_yaml_id, diag_file_ids) &
result(total_nvars)

  integer, intent(in) :: diag_yaml_id     !< Id for the diag_table yaml
  integer, intent(in) :: diag_file_ids(:) !< Ids of the files in the diag_table yaml
  integer :: total_nvars

  integer :: i !< For do loop

  total_nvars = 0
  do i = 1, size(diag_file_ids,1)
    total_nvars = total_nvars + get_num_blocks(diag_yaml_id, "varlist", parent_block_id=diag_file_ids(i))
  end do
end function

!> @brief This checks if the file frequency in a diag file is valid and crashes if it isn't
subroutine check_file_freq(fileobj)
  type(diagYamlFiles_type), intent(inout) :: fileobj      !< diagYamlFiles_type obj to check

  if (fileobj%file_freq < 1 ) &
    call mpp_error(FATAL, "freq must be greater than 0. &
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

!> @brief This checks if the realm in a diag file is valid and crashes if it isn't
subroutine check_file_realm(fileobj)
  type(diagYamlFiles_type), intent(inout) :: fileobj      !< diagYamlFiles_type obj to checK

  select case (TRIM(fileobj%file_realm))
  case ("ATM", "OCN", "LND", "ICE", "")
  case default
    call mpp_error(FATAL, trim(fileobj%file_realm)//" is an invalid realm! &
      &The acceptable values are ATM, OCN, LND, ICE. &
      &Check your entry for file:"//trim(fileobj%file_fname))
  end select

end subroutine check_file_realm

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
  case ("double", "float")
  case default
    call mpp_error(FATAL, trim(field%var_skind)//" is an invalid kind! &
      &The acceptable values are double and float. &
      &Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
  end select

end subroutine check_field_kind

!> @brief This checks if the reduction of a diag field is valid and crashes if it isn't
subroutine check_field_reduction(field)
  type(diagYamlFilesVar_type), intent(in) :: field        !< diagYamlFilesVar_type obj to read the contents into

  integer :: n_diurnal !< number of diurnal samples
  integer :: pow_value !< The power value
  integer :: ioerror   !< io error status after reading in the diurnal samples

  n_diurnal = 0
  pow_value = 0
  ioerror = 0
  if (field%var_reduction(1:7) .eq. "diurnal") then
    READ (UNIT=field%var_reduction(8:LEN_TRIM(field%var_reduction)), FMT=*, IOSTAT=ioerror) n_diurnal
    if (ioerror .ne. 0) &
      call mpp_error(FATAL, "Error getting the number of diurnal samples from "//trim(field%var_reduction))
    if (n_diurnal .le. 0) &
      call mpp_error(FATAL, "Diurnal samples should be greater than 0. &
        & Check your entry for file:"//trim(field%var_varname)//" in file "//trim(field%var_fname))
  elseif (field%var_reduction(1:3) .eq. "pow") then
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
!> @brief Inquiry for diag_files_obj%file_write
!! @return file_write of a diag_yaml_file_obj
pure function get_file_write(diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 logical :: res !< What is returned
  res = diag_files_obj%file_write
end function get_file_write
!> @brief Inquiry for diag_files_obj%file_realm
!! @return file_realm of a diag_yaml_file_obj
pure function get_file_realm(diag_files_obj) &
result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (:), allocatable :: res !< What is returned
  res = diag_files_obj%file_realm
end function get_file_realm
!> @brief Inquiry for diag_files_obj%file_subregion
!! @return file_sub_region of a diag_yaml_file_obj
pure function get_file_sub_region (diag_files_obj) &
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
!> @brief Inquiry for diag_yaml_files_var_obj%var_write
!! @return var_write of a diag_yaml_files_var_obj
pure function get_var_write (diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 logical :: res !< What is returned
  res = diag_var_obj%var_write
end function get_var_write
!> @brief Inquiry for diag_yaml_files_var_obj%var_attributes
!! @return var_attributes of a diag_yaml_files_var_obj
pure function get_var_attributes(diag_var_obj) &
result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=MAX_STR_LEN), allocatable :: res (:,:) !< What is returned
 res = diag_var_obj%var_attributes
end function get_var_attributes
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
  obj%file_write          = .true.
  obj%file_duration       = DIAG_NULL
  obj%file_new_file_freq  = DIAG_NULL
  obj%file_sub_region%tile = DIAG_NULL
end subroutine diag_yaml_files_obj_init

#endif
end module fms_diag_yaml_mod
!> @}
! close documentation grouping
