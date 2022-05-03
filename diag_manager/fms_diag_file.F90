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
!! @brief fms_diag_file_mod handles the file objects data, functions, and subroutines.
!! @author Tom Robinson
!! @description The fmsDiagFile_type contains the information for each history file to be written.  It has
!! a pointer to the information from the diag yaml, additional metadata that comes from the model, and a
!! list of the variables and their variable IDs that are in the file.
module fms_diag_file_mod
!use mpp_mod
use fms2_io_mod, only: FmsNetcdfFile_t, FmsNetcdfUnstructuredDomainFile_t, FmsNetcdfDomainFile_t
#ifdef use_yaml
use fms_diag_yaml_mod, only: diagYamlObject_type
#endif

implicit none
private

public :: fmsDiagFile_type, FMS_diag_files

integer, parameter :: var_string_len = 25

type :: fmsDiagFile_type
 private
  integer :: id !< The number associated with this file in the larger array of files
  class(FmsNetcdfFile_t), allocatable :: fileobj !< fms2_io file object for this history file 
  character(len=1) :: file_domain_type !< (I don't think we will need this)
#ifdef use_yaml
  type(diagYamlObject_type), pointer :: diag_yaml => null() !< Pointer to the diag_yaml data
#endif
  character(len=:) , dimension(:), allocatable :: file_metadata_from_model !< File metadata that comes from
                                                                           !! the model.
  character(len=var_string_len), dimension(:), allocatable :: var_list !< List of the variables by name
  integer, dimension(:), allocatable :: var_ids !< Variable IDs corresponding to var_list
  integer, dimension(:), private, allocatable :: var_index !< An array of the variable indicies in the 
                                                                 !! diag_object.  This should be the same size as
                                                                 !! `file_varlist`
  logical, dimension(:), private, allocatable :: var_reg   !< Array corresponding to `file_varlist`, .true. 
                                                                 !! if the variable has been registered and 
                                                                 !! `file_var_index` has been set for the variable

 contains
  procedure, public :: has_file_metadata_from_model
  procedure, public :: has_fileobj
#ifdef use_yaml
  procedure, public :: has_diag_yaml
#endif
  procedure, public :: has_var_ids
  procedure, public :: has_var_list
  procedure, public :: get_id
! TODO  procedure, public :: get_fileobj ! TODO
  procedure, public :: get_file_domain_type
! TODO  procedure, public :: get_diag_yaml ! TODO
  procedure, public :: get_file_metadata_from_model
  procedure, public :: get_var_list
  procedure, public :: get_var_ids

end type fmsDiagFile_type

type(fmsDiagFile_type), dimension (:), allocatable, target :: FMS_diag_files !< The array of diag files

contains

!> \brief Logical function to determine if the variable file_metadata_from_model has been allocated or associated
!! \return .True. if file_metadata_from_model exists .False. if file_metadata_from_model has not been set
pure logical function has_file_metadata_from_model (obj)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  has_file_metadata_from_model = allocated(obj%file_metadata_from_model)
end function has_file_metadata_from_model
!> \brief Logical function to determine if the variable fileobj has been allocated or associated
!! \return .True. if fileobj exists .False. if fileobj has not been set
pure logical function has_fileobj (obj)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  has_fileobj = allocated(obj%fileobj)
end function has_fileobj
#ifdef use_yaml
!> \brief Logical function to determine if the variable diag_yaml has been allocated or associated
!! \return .True. if diag_yaml exists .False. if diag_yaml has not been set
pure logical function has_diag_yaml (obj)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  has_diag_yaml = associated(obj%diag_yaml)
end function has_diag_yaml
#endif
!> \brief Logical function to determine if the variable var_ids has been allocated or associated
!! \return .True. if var_ids exists .False. if var_ids has not been set
pure logical function has_var_ids (obj)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  has_var_ids = allocated(obj%var_ids)
end function has_var_ids
!> \brief Logical function to determine if the variable var_list has been allocated or associated
!! \return .True. if var_list exists .False. if var_list has not been set
pure logical function has_var_list (obj)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  has_var_list = allocated(obj%var_list)
end function has_var_list
!> \brief Returns a copy of the value of id
!! \return A copy of id
pure function get_id (obj) result (res)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  integer :: res
  res = obj%id
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
!> \brief Returns a copy of the value of file_domain_type
!! \return A copy of file_domain_type
pure function get_file_domain_type (obj) result (res)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  character(1) :: res
  res = obj%file_domain_type
end function get_file_domain_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TODO
!!> \brief Returns a copy of the value of diag_yaml
!!! \return A copy of diag_yaml
!#ifdef use_yaml
!pure function get_diag_yaml (obj) result (res)
!  class(fmsDiagFile_type), intent(in) :: obj !< The file object
!  type(diagYamlFiles_type) :: res
!  res = obj%diag_yaml
!end function get_diag_yaml
!#endif
!> \brief Returns a copy of the value of file_metadata_from_model
!! \return A copy of file_metadata_from_model
pure function get_file_metadata_from_model (obj) result (res)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  character(len=:), dimension(:), allocatable :: res
  res = obj%file_metadata_from_model
end function get_file_metadata_from_model
!> \brief Returns a copy of the value of var_list
!! \return A copy of var_list
pure function get_var_list (obj) result (res)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  character(len=:), dimension(:), allocatable :: res
!  allocate(character(len=len(obj%var_list(1))), dimension (size(obj%var_list(1))) :: res) 
  res = obj%var_list
end function get_var_list
!> \brief Returns a copy of the value of var_ids
!! \return A copy of var_ids
pure function get_var_ids (obj) result (res)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  integer, dimension(:), allocatable :: res
  allocate(res(size(obj%var_ids)))
  res = obj%var_ids
end function get_var_ids


end module fms_diag_file_mod
