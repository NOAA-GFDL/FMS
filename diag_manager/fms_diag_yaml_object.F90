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
!> @defgroup fms_diag_yaml_object_mod fms_diag_yaml_object_mod
!> @ingroup diag_manager
!! @brief The diag yaml objects are handled here, with variables the correspond to
!! entries in the diag yaml file.  The actual parsing of the yaml is handled in
!! @ref fms_diag_yaml_mod.
!! @author Tom Robinson, Uriel Ramirez

!> @file
!> @brief File for @ref fms_diag_yaml_object_mod

!> @addtogroup fms_diag_yaml_object_mod
!> @{
module fms_diag_yaml_object_mod

use fms_mod , only: fms_c2f_string
use iso_c_binding
      implicit none
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

type diagYamlFiles_type
     character (len=:), allocatable :: file_fname !< file name
     character (len=:), allocatable :: file_frequnit !< the frequency unit
     integer (c_int)    :: file_freq !< the frequency of data
     character (len=:), allocatable :: file_timeunit !< The unit of time
     character (len=:), allocatable :: file_unlimdim !< The name of the unlimited dimension
     logical :: file_write
     character (len=:), allocatable :: string_file_write !< false if the user doesn’t want the file to be
                                               !! created (default is true).
     character (len=:), allocatable :: file_realm !< The modeling realm that the variables come from
     type(subRegion_type) :: file_sub_region !< type containing info about the subregion, if any
     integer :: file_new_file_freq !< Frequency for closing the existing file
     character (len=:), allocatable :: file_new_file_freq_units !< Time units for creating a new file.
                                                        !! Required if “new_file_freq” used
     character (len=:), allocatable :: file_start_time !< Time to start the file for the first time. Requires
                                                       !! “new_file_freq”
     integer :: file_duration !< How long the file should receive data after start time
                                      !! in “file_duration_units”.  This optional field can only
                                      !! be used if the start_time field is present.  If this field
                                      !! is absent, then the file duration will be equal to the
                                      !! frequency for creating new files.
                                      !! NOTE: The file_duration_units field must also be present if
                                      !! this field is present.
    character (len=:), allocatable :: file_duration_units !< The file duration units
    !< Need to use `MAX_STR_LEN` because not all filenames/global attributes are the same length
    character (len=MAX_STR_LEN), dimension(:), allocatable :: file_varlist !< An array of variable names
                                                             !! within a file
    character (len=MAX_STR_LEN), dimension(:,:), allocatable :: file_global_meta !< Array of key(dim=1)
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

type diagYamlFilesVar_type
     character (len=:), allocatable :: var_fname !< The field/diagnostic name
     character (len=:), allocatable :: var_varname !< The name of the variable
     character (len=:), allocatable :: var_reduction !< Reduction to be done on var
     character (len=:), allocatable :: var_module !< The module that th variable is in
     character (len=:), allocatable :: var_skind !< The type/kind of the variable
     character (len=:), allocatable :: string_var_write !< false if the user doesn’t want the variable to be
                                             !! written to the file (default: true).
     logical :: var_write !< false if the user doesn’t want the variable to be
                          !! written to the file (default: true).
     character (len=:), allocatable :: var_outname !< Name of the variable as written to the file
     character (len=:), allocatable :: var_longname !< Overwrites the long name of the variable
     character (len=:), allocatable :: var_units !< Overwrites the units
     !< Need to use `MAX_STR_LEN` because not all filenames/global attributes are the same length
     character (len=MAX_STR_LEN), dimension (:, :), allocatable :: var_attributes !< Attributes to overwrite or
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

contains
!!!!!!! YAML FILE INQUIRIES !!!!!!!
!> @brief Inquiry for diag_files_obj%file_fname
!! @return file_fname of a diag_yaml_file obj
pure function get_file_fname (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_files_obj%file_fname
end function get_file_fname
!> @brief Inquiry for diag_files_obj%file_frequnit
!! @return file_frequnit of a diag_yaml_file_obj
pure function get_file_frequnit (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_files_obj%file_frequnit
end function get_file_frequnit
!> @brief Inquiry for diag_files_obj%file_freq
!! @return file_freq of a diag_yaml_file_obj
pure function get_file_freq(diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_freq
end function get_file_freq
!> @brief Inquiry for diag_files_obj%file_timeunit
!! @return file_timeunit of a diag_yaml_file_obj
pure function get_file_timeunit (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_files_obj%file_timeunit
end function get_file_timeunit
!> @brief Inquiry for diag_files_obj%file_unlimdim
!! @return file_unlimdim of a diag_yaml_file_obj
pure function get_file_unlimdim(diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_files_obj%file_unlimdim
end function get_file_unlimdim
!> @brief Inquiry for diag_files_obj%file_write
!! @return file_write of a diag_yaml_file_obj
pure function get_file_write(diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 logical :: res !< What is returned
  res = diag_files_obj%file_write
end function get_file_write
!> @brief Inquiry for diag_files_obj%file_realm
!! @return file_realm of a diag_yaml_file_obj
pure function get_file_realm(diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (:), allocatable :: res !< What is returned
  res = diag_files_obj%file_realm
end function get_file_realm
!> @brief Inquiry for diag_files_obj%file_subregion
!! @return file_sub_region of a diag_yaml_file_obj
pure function get_file_sub_region (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 type(subRegion_type) :: res !< What is returned
  res = diag_files_obj%file_sub_region
end function get_file_sub_region
!> @brief Inquiry for diag_files_obj%file_new_file_freq
!! @return file_new_file_freq of a diag_yaml_file_obj
pure function get_file_new_file_freq(diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_new_file_freq
end function get_file_new_file_freq
!> @brief Inquiry for diag_files_obj%file_new_file_freq_units
!! @return file_new_file_freq_units of a diag_yaml_file_obj
pure function get_file_new_file_freq_units (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (:), allocatable :: res !< What is returned
  res = diag_files_obj%file_new_file_freq_units
end function get_file_new_file_freq_units
!> @brief Inquiry for diag_files_obj%file_start_time
!! @return file_start_time of a diag_yaml_file_obj
pure function get_file_start_time (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_files_obj%file_start_time
end function get_file_start_time
!> @brief Inquiry for diag_files_obj%file_duration
!! @return file_duration of a diag_yaml_file_obj
pure function get_file_duration (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_duration
end function get_file_duration
!> @brief Inquiry for diag_files_obj%file_duration_units
!! @return file_duration_units of a diag_yaml_file_obj
pure function get_file_duration_units (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
  character (:), allocatable :: res !< What is returned
  res = diag_files_obj%file_duration_units
end function get_file_duration_units
!> @brief Inquiry for diag_files_obj%file_varlist
!! @return file_varlist of a diag_yaml_file_obj
pure function get_file_varlist (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (:), allocatable :: res(:) !< What is returned
  res = diag_files_obj%file_varlist
end function get_file_varlist
!> @brief Inquiry for diag_files_obj%file_global_meta
!! @return file_global_meta of a diag_yaml_file_obj
pure function get_file_global_meta (diag_files_obj) result (res)
 class (diagYamlFiles_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (:), allocatable :: res(:,:) !< What is returned
  res = diag_files_obj%file_global_meta
end function get_file_global_meta
!> @brief Inquiry for whether file_global_meta is allocated
!! @return Flag indicating if file_global_meta is allocated
function is_global_meta(diag_files_obj) result(res)
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
pure function get_var_fname (diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_fname
end function get_var_fname
!> @brief Inquiry for diag_yaml_files_var_obj%var_varname
!! @return var_varname of a diag_yaml_files_var_obj
pure function get_var_varname (diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_varname
end function get_var_varname
!> @brief Inquiry for diag_yaml_files_var_obj%var_reduction
!! @return var_reduction of a diag_yaml_files_var_obj
pure function get_var_reduction (diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_reduction
end function get_var_reduction
!> @brief Inquiry for diag_yaml_files_var_obj%var_module
!! @return var_module of a diag_yaml_files_var_obj
pure function get_var_module (diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_module
end function get_var_module
!> @brief Inquiry for diag_yaml_files_var_obj%var_skind
!! @return var_skind of a diag_yaml_files_var_obj
pure function get_var_skind (diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_skind
end function get_var_skind
!> @brief Inquiry for diag_yaml_files_var_obj%var_outname
!! @return var_outname of a diag_yaml_files_var_obj
pure function get_var_outname (diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_outname
end function get_var_outname
!> @brief Inquiry for diag_yaml_files_var_obj%var_longname
!! @return var_longname of a diag_yaml_files_var_obj
pure function get_var_longname (diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_longname
end function get_var_longname
!> @brief Inquiry for diag_yaml_files_var_obj%var_units
!! @return var_units of a diag_yaml_files_var_obj
pure function get_var_units (diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=:), allocatable :: res !< What is returned
  res = diag_var_obj%var_units
end function get_var_units
!> @brief Inquiry for diag_yaml_files_var_obj%var_write
!! @return var_write of a diag_yaml_files_var_obj
pure function get_var_write (diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 logical :: res !< What is returned
  res = diag_var_obj%var_write
end function get_var_write
!> @brief Inquiry for diag_yaml_files_var_obj%var_attributes
!! @return var_attributes of a diag_yaml_files_var_obj
pure function get_var_attributes(diag_var_obj) result (res)
 class (diagYamlFilesVar_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=MAX_STR_LEN), allocatable :: res (:,:) !< What is returned
 res = diag_var_obj%var_attributes
end function get_var_attributes
!> @brief Inquiry for whether var_attributes is allocated
!! @return Flag indicating if var_attributes is allocated
function is_var_attributes(diag_var_obj) result(res)
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

  obj%file_freq           = 0
  obj%file_write          = .true.
  obj%file_duration       = 0
  obj%file_new_file_freq  = 0
  obj%file_sub_region%tile = 0
end subroutine diag_yaml_files_obj_init

end module fms_diag_yaml_object_mod
!> @}
! close documentation grouping

