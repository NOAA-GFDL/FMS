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
!! @author Tom Robinson

!> @file
!> @brief File for @ref fms_diag_yaml_object_mod

!> @addtogroup fms_diag_yaml_object_mod
!> @{
module fms_diag_yaml_object_mod

use fms_mod , only: fms_c2f_string
use iso_c_binding
      implicit none
integer, parameter :: NUM_REGION_ARRAY = 8
  !> @brief The files type matching a C struct containing diag_yaml information
    !> @ingroup fms_diag_files_mod
type, bind(c) :: diag_yaml_files_struct
     character (kind=c_char) :: file_fname (20) !< file name
     character (kind=c_char) :: file_frequnit (7) !< the frequency unit
     integer (c_int)    :: file_freq !< the frequency of data
     character (kind=c_char) :: file_timeunit(7) !< The unit of time
     character (kind=c_char) :: file_unlimdim(8) !< The name of the unlimited dimension
     character (kind=c_char) :: file_write (5) !< false if the user doesn’t want the file to be 
                                               !! created (default is true).
     character (kind=c_char) :: file_realm (3) !< The modeling realm that the variables come from
     real (c_float) :: file_region (NUM_REGION_ARRAY) !< Bounds of the regional section to capture
     integer (c_int) :: file_new_file_freq !< Frequency for closing the existing file
     character (kind=c_char) :: file_new_file_freq_units (3)!< Time units for creating a new file. 
                                                        !! Required if “new_file_freq” used
     integer (c_int) :: file_start_time !< Time to start the file for the first time. Requires “new_file_freq”
     integer (c_int) :: file_duration !< How long the file should receive data after start time 
                                      !! in “file_duration_units”.  This optional field can only 
                                      !! be used if the start_time field is present.  If this field
                                      !! is absent, then the file duration will be equal to the
                                      !! frequency for creating new files.  
                                      !! NOTE: The file_duration_units field must also be present if
                                      !! this field is present.
    character (kind=c_char) :: file_duration_units (3)!< The file duration units
end type diag_yaml_files_struct

type diag_yaml_files_type
     character (len=:), allocatable :: file_fname !< file name
     character (len=:), allocatable :: file_frequnit !< the frequency unit
     integer (c_int)    :: file_freq !< the frequency of data
     character (len=:), allocatable :: file_timeunit !< The unit of time
     character (len=:), allocatable :: file_unlimdim !< The name of the unlimited dimension
     logical :: file_write
     character (len=:), allocatable :: string_file_write !< false if the user doesn’t want the file to be 
                                               !! created (default is true).
     character (len=:), allocatable :: file_realm !< The modeling realm that the variables come from
     real :: file_region (NUM_REGION_ARRAY) !< Bounds of the regional section to capture
     integer :: file_new_file_freq !< Frequency for closing the existing file
     character (len=:), allocatable :: file_new_file_freq_units !< Time units for creating a new file. 
                                                        !! Required if “new_file_freq” used
     integer :: file_start_time !< Time to start the file for the first time. Requires “new_file_freq”
     integer :: file_duration !< How long the file should receive data after start time 
                                      !! in “file_duration_units”.  This optional field can only 
                                      !! be used if the start_time field is present.  If this field
                                      !! is absent, then the file duration will be equal to the
                                      !! frequency for creating new files.  
                                      !! NOTE: The file_duration_units field must also be present if
                                      !! this field is present.
    character (len=:), allocatable :: file_duration_units !< The file duration units
    character (len=:), dimension(:), allocatable :: file_varlist !< An array of variable names
                                                             !! within a file
    character (len=:), dimension(:,:), allocatable :: file_global_meta !< Array of key(dim=1) 
                                                        !! and values(dim=2) to be added as global 
                                                        !! meta data to the file

 contains
 procedure :: copy_struct => copy_file_struct_to_object
 procedure :: fname => get_file_fname
 procedure :: frequnit => get_file_frequnit
 procedure :: freq => get_file_freq
 procedure :: timeunit => get_file_timeunit
 procedure :: unlimdim => get_file_unlimdim
 procedure :: write_file => get_file_write
 procedure :: realm => get_file_realm
 procedure :: region => get_file_region
 procedure :: new_file_freq => get_file_new_file_freq
 procedure :: new_file_freq_units => get_file_new_file_freq_units
 procedure :: start_time => get_file_start_time
 procedure :: duration => get_file_duration
 procedure :: duration_units => get_file_duration_units
 procedure :: varlist => get_file_varlist
 procedure :: global_meta => get_file_global_meta

end type diag_yaml_files_type

!> @brief The field type matching the C struct for diag_yaml information
  !> @ingroup fms_diag_files_mod
type, bind(c) :: diag_yaml_files_var_struct
     character (kind=c_char) :: var_fname (20) !< The field/diagnostic name
     character (kind=c_char) :: var_varname(20) !< The name of the variable
     character (kind=c_char) :: var_reduction(20) !< Reduction to be done on var
     character (kind=c_char) :: var_module(20) !< The module that th variable is in
     character (kind=c_char) :: var_skind(8) !< The type/kind of the variable
     character (kind=c_char) :: var_write(5) !< false if the user doesn’t want the variable to be
                                             !! written to the file (default: true).
     character (kind=c_char) :: var_outname(20) !< Name of the variable as written to the file
     character (kind=c_char) :: var_longname(100) !< Overwrites the long name of the variable
     character (kind=c_char) :: var_units(10) !< Overwrites the units
end type diag_yaml_files_var_struct

type diag_yaml_files_var_type
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
     character (len=:), dimension (:), allocatable :: var_attributes !< Attributes to overwrite or
                                                                     !! add from diag_yaml
 contains
  procedure :: copy_struct => copy_variable_struct_to_object 
  procedure :: fname => get_var_fname
  procedure :: varname => get_var_varname
  procedure :: reduction => get_var_reduction
  procedure :: module_var => get_var_module
  procedure :: skind => get_var_skind
  procedure :: outname => get_var_outname
  procedure :: longname => get_var_longname
  procedure :: units => get_var_units
  procedure :: write_var => get_var_write
  procedure :: attr => get_var_attributes

end type diag_yaml_files_var_type

contains
!!!!!!!! YAML FILE ROUTINES !!!!!!!!
!< \brief Copies the information of the yaml struct to the fortran object holding the info
subroutine copy_file_struct_to_object(diag_files_obj, diag_files_struct)
 class(diag_yaml_files_type) :: diag_files_obj !< Fortran-side object with diag_yaml info
 type(diag_yaml_files_struct) :: diag_files_struct !< The C struct that has the diag_yaml
                                                   !! info
 integer :: i !< For looping 
!< Convert the C strings to Fortran strings
 diag_files_obj%file_fname = fms_c2f_string (diag_files_struct%file_fname)
 diag_files_obj%file_frequnit = fms_c2f_string (diag_files_struct%file_frequnit)
 diag_files_obj%file_timeunit = fms_c2f_string (diag_files_struct%file_timeunit)
 diag_files_obj%file_unlimdim = fms_c2f_string (diag_files_struct%file_unlimdim)
 diag_files_obj%file_realm = fms_c2f_string (diag_files_struct%file_realm)
 diag_files_obj%file_new_file_freq_units = fms_c2f_string (diag_files_struct%file_new_file_freq_units)
 diag_files_obj%file_duration_units = fms_c2f_string (diag_files_struct%file_duration_units)
!< Set the file_write to be true or false
 diag_files_obj%string_file_write = fms_c2f_string (diag_files_struct%file_write)
 diag_files_obj%file_write = .true.
 if (diag_files_obj%string_file_write(1:1)=="f" .or. &
     diag_files_obj%string_file_write(1:1)=="F") &
        diag_files_obj%file_write = .false.
 deallocate (diag_files_obj%string_file_write)
!< Store the numbers
 diag_files_obj%file_freq = diag_files_struct%file_freq
!$omp simd
 do i = 1, NUM_REGION_ARRAY
        diag_files_obj%file_region(i) = diag_files_struct%file_region(i)
 enddo
 diag_files_obj%file_new_file_freq = diag_files_struct%file_new_file_freq
 diag_files_obj%file_start_time = diag_files_struct%file_start_time
 diag_files_obj%file_duration = diag_files_struct%file_duration

end subroutine copy_file_struct_to_object
!!!!!!! YAML FILE INQUIRIES !!!!!!!
!> \brief Inquiry for diag_files_obj%file_fname
pure function get_file_fname (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_files_obj%file_fname
end function get_file_fname
!> \brief Inquiry for diag_files_obj%file_frequnit
pure function get_file_frequnit (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_files_obj%file_frequnit
end function get_file_frequnit
!> \brief Inquiry for diag_files_obj%file_freq
pure function get_file_freq(diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_freq
end function get_file_freq
!> \brief Inquiry for diag_files_obj%file_timeunit
pure function get_file_timeunit (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_files_obj%file_timeunit
end function get_file_timeunit
!> \brief Inquiry for diag_files_obj%file_unlimdim
pure function get_file_unlimdim(diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_files_obj%file_unlimdim
end function get_file_unlimdim
!> \brief Inquiry for diag_files_obj%file_write
pure function get_file_write(diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 logical :: res !< What is returned
  res = diag_files_obj%file_write
end function get_file_write
!> \brief Inquiry for diag_files_obj%file_realm
pure function get_file_realm(diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (*) :: res !< What is returned
  res = diag_files_obj%file_realm
end function get_file_realm
!> \brief Inquiry for diag_files_obj%file_region
pure function get_file_region (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 real :: res (NUM_REGION_ARRAY) !< What is returned
  res = diag_files_obj%file_region
end function get_file_region
!> \brief Inquiry for diag_files_obj%file_new_file_freq
pure function get_file_new_file_freq(diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_new_file_freq
end function get_file_new_file_freq
!> \brief Inquiry for diag_files_obj%file_new_file_freq_units
pure function get_file_new_file_freq_units (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (*) :: res !< What is returned
  res = diag_files_obj%file_new_file_freq_units
end function get_file_new_file_freq_units
!> \brief Inquiry for diag_files_obj%file_start_time
pure function get_file_start_time (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_start_time
end function get_file_start_time
!> \brief Inquiry for diag_files_obj%file_duration
pure function get_file_duration (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 integer :: res !< What is returned
  res = diag_files_obj%file_duration
end function get_file_duration
!> \brief Inquiry for diag_files_obj%file_duration_units
pure function get_file_duration_units (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
  character (*) :: res !< What is returned
  res = diag_files_obj%file_duration_units
end function get_file_duration_units
!> \brief Inquiry for diag_files_obj%file_varlist
pure function get_file_varlist (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (*) :: res(:) !< What is returned
  res = diag_files_obj%file_varlist
end function get_file_varlist
!> \brief Inquiry for diag_files_obj%file_global_meta
pure function get_file_global_meta (diag_files_obj) result (res)
 class (diag_yaml_files_type), intent(in) :: diag_files_obj !< The object being inquiried
 character (*) :: res(:,:) !< What is returned
  res = diag_files_obj%file_global_meta
end function get_file_global_meta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! VARIABLES ROUTINES AND FUNCTIONS !!!!!!!
!< \brief Copies the information of the yaml struct to the fortran object holding the var info
subroutine copy_variable_struct_to_object(diag_var_obj, diag_var_struct)
 class(diag_yaml_files_var_type) :: diag_var_obj !< Fortran-side object with diag_yaml var info
 type(diag_yaml_files_var_struct) :: diag_var_struct !< The C struct that has the diag_yaml
                                                   !! var info
!< Convert the C strings to Fortran strings
 diag_var_obj%var_fname = fms_c2f_string (diag_var_struct%var_fname)
 diag_var_obj%var_varname = fms_c2f_string (diag_var_struct%var_varname)
 diag_var_obj%var_reduction = fms_c2f_string (diag_var_struct%var_reduction)
 diag_var_obj%var_module = fms_c2f_string (diag_var_struct%var_module)
 diag_var_obj%var_skind = fms_c2f_string (diag_var_struct%var_skind)
 diag_var_obj%var_outname = fms_c2f_string (diag_var_struct%var_outname)
 diag_var_obj%var_longname = fms_c2f_string (diag_var_struct%var_longname)
 diag_var_obj%var_units = fms_c2f_string (diag_var_struct%var_units)
!< Set the file_write to be true or false
 diag_var_obj%string_var_write= fms_c2f_string (diag_var_struct%var_write)
 diag_var_obj%var_write= .true.
 if (diag_var_obj%string_var_write(1:1)=="f" .or. &
     diag_var_obj%string_var_write(1:1)=="F") &
        diag_var_obj%var_write= .false.
 deallocate (diag_var_obj%string_var_write)
end subroutine copy_variable_struct_to_object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! YAML VAR INQUIRIES !!!!!!!
!> \brief Inquiry for diag_yaml_files_var_obj%var_fname
pure function get_var_fname (diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_var_objn%var_fname
end function get_var_fname
!> \brief Inquiry for diag_yaml_files_var_obj%var_varname
pure function get_var_varname (diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_var_obj%var_varname
end function get_var_varname
!> \brief Inquiry for diag_yaml_files_var_obj%var_reduction
pure function get_var_reduction (diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_var_obj%var_reduction
end function get_var_reduction
!> \brief Inquiry for diag_yaml_files_var_obj%var_module
pure function get_var_module (diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_var_obj%var_module
end function get_var_module
!> \brief Inquiry for diag_yaml_files_var_obj%var_skind
pure function get_var_skind (diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_var_obj%var_skind
end function get_var_skind
!> \brief Inquiry for diag_yaml_files_var_obj%var_outname
pure function get_var_outname (diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_var_obj%var_outname
end function get_var_outname
!> \brief Inquiry for diag_yaml_files_var_obj%var_longname
pure function get_var_longname (diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_var_obj%var_longname
end function get_var_longname
!> \brief Inquiry for diag_yaml_files_var_obj%var_units
pure function get_var_units (diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=*) :: res !< What is returned
  res = diag_var_obj%var_units
end function get_var_units
!> \brief Inquiry for diag_yaml_files_var_obj%var_write
pure function get_var_write (diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 logical :: res !< What is returned
  res = diag_var_obj%var_write
end function get_var_write
!> \brief Inquiry for diag_yaml_files_var_obj%var_attributes
pure function get_var_attributes(diag_var_obj) result (res)
 class (diag_yaml_files_var_type), intent(in) :: diag_var_obj !< The object being inquiried
 character (len=*) :: res (:) !< What is returned
 res = diag_var_obj%var_attributes
end function get_var_attributes

end module fms_diag_yaml_object_mod
!> @}
! close documentation grouping

