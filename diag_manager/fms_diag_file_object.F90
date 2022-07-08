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
use fms2_io_mod, only: FmsNetcdfFile_t, FmsNetcdfUnstructuredDomainFile_t, FmsNetcdfDomainFile_t
use diag_data_mod, only: DIAG_NULL, NO_DOMAIN, max_axes, SUB_REGIONAL, get_base_time
use diag_util_mod, only: diag_time_inc
use time_manager_mod, only: time_type, operator(/=)
#ifdef use_yaml
use fms_diag_yaml_mod, only: diag_yaml, diagYamlObject_type, diagYamlFiles_type
#endif
use fms_diag_axis_object_mod, only: diagDomain_t
use mpp_mod, only: mpp_error, FATAL
implicit none
private

public :: fmsDiagFile_type, FMS_diag_files, fms_diag_files_object_init, fms_diag_files_object_initialized

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
#ifdef use_yaml
  type(diagYamlFiles_type), pointer :: diag_yaml_file => null() !< Pointer to the diag_yaml_file data
#endif
  integer                                      :: type_of_domain !< The type of domain to use to open the file
                                                                 !! NO_DOMAIN, TWO_D_DOMAIN, UG_DOMAIN, SUB_REGIONAL
  class(diagDomain_t), pointer                 :: domain         !< The domain to use,
                                                                 !! null if NO_DOMAIN or SUB_REGIONAL
  character(len=:) , dimension(:), allocatable :: file_metadata_from_model !< File metadata that comes from
                                                                           !! the model.
  integer, dimension(:), allocatable :: var_ids !< Variable IDs corresponding to file_varlist
  integer, dimension(:), private, allocatable :: var_index !< An array of the variable indicies in the
                                                                 !! diag_object.  This should be the same size as
                                                                 !! `file_varlist`
  logical, dimension(:), private, allocatable :: var_reg   !< Array corresponding to `file_varlist`, .true.
                                                                 !! if the variable has been registered and
                                                                 !! `file_var_index` has been set for the variable
  integer, dimension(:), allocatable :: axis_ids !< Array of axis ids in the file
  integer, dimension(:), allocatable :: sub_axis_ids !< Array of axis ids in the file
  integer :: number_of_axis !< Number of axis in the file

 contains
  procedure, public :: has_file_metadata_from_model
  procedure, public :: has_fileobj
#ifdef use_yaml
  procedure, public :: has_diag_yaml_file
  procedure, public :: set_file_domain
  procedure, public :: add_axes
  procedure, public :: add_start_time
#endif
  procedure, public :: has_var_ids
  procedure, public :: get_id
! TODO  procedure, public :: get_fileobj ! TODO
! TODO  procedure, public :: get_diag_yaml_file ! TODO
  procedure, public :: get_file_metadata_from_model
  procedure, public :: get_var_ids
! The following fuctions come will use the yaml inquiry functions
#ifdef use_yaml
 procedure, public :: get_file_fname
 procedure, public :: get_file_frequnit
 procedure, public :: get_file_freq
 procedure, public :: get_file_timeunit
 procedure, public :: get_file_unlimdim
!! TODO get functions for sub region stuff
! procedure, public :: get_file_sub_region
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
#endif

end type fmsDiagFile_type

type(fmsDiagFile_type), dimension (:), allocatable, target :: FMS_diag_files !< The array of diag files

contains

!< @brief Allocates the number of files and sets an ID based for each file
!! @return true if there are files allocated in the YAML object
logical function fms_diag_files_object_init ()
#ifdef use_yaml
  integer :: nFiles !< Number of files in the diag yaml
  integer :: i !< Looping iterator
  type(fmsDiagFile_type), pointer :: obj !< FMS_diag_files(i) (for less typing)
  if (diag_yaml%has_diag_files()) then
   nFiles = diag_yaml%size_diag_files()
   allocate (FMS_diag_files(nFiles))
   set_ids_loop: do i= 1,nFiles
     obj => FMS_diag_files(i)
     obj%diag_yaml_file => diag_yaml%diag_files(i)
     obj%id = i
     allocate(obj%var_ids(diag_yaml%diag_files(i)%size_file_varlist()))
     allocate(obj%var_index(diag_yaml%diag_files(i)%size_file_varlist()))
     allocate(obj%var_reg(diag_yaml%diag_files(i)%size_file_varlist()))
     !! Initialize the integer arrays
     obj%var_ids = DIAG_NULL
     obj%var_reg = .FALSE.
     obj%var_index = DIAG_NULL

     !> These will be set in a set_file_domain
     obj%type_of_domain = NO_DOMAIN
     obj%domain => null()

     !> This will be set in a add_axes
     allocate(obj%axis_ids(max_axes))

     !> If the file has a sub_regional, define it as one and allocate the sub_axis_ids array.
     !! This will be set in a add_axes
     if (obj%has_file_sub_region()) then
      obj%type_of_domain = SUB_REGIONAL
      allocate(obj%sub_axis_ids(max_axes))
      obj%sub_axis_ids = diag_null
     endif

     obj%number_of_axis = 0

     !> Set the start_time of the file to the base_time and set up the *_output variables
     obj%start_time = get_base_time()
     obj%last_output = get_base_time()
     obj%next_output = diag_time_inc(obj%start_time, obj%get_file_freq(), obj%get_file_frequnit())
     obj%next_next_output = diag_time_inc(obj%next_output, obj%get_file_freq(), obj%get_file_frequnit())

     nullify(obj)
   enddo set_ids_loop
   fms_diag_files_object_init = .true.
  else
   fms_diag_files_object_init = .false.
!  mpp_error("fms_diag_files_object_init: The diag_table.yaml file has not been correctly parsed.",&
!    FATAL)
  endif
#else
  fms_diag_files_object_init = .false.
#endif
end function fms_diag_files_object_init
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
!> \brief Logical function to determine if the variable diag_yaml_file has been allocated or associated
!! \return .True. if diag_yaml_file exists .False. if diag_yaml has not been set
pure logical function has_diag_yaml_file (obj)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  has_diag_yaml_file = associated(obj%diag_yaml_file)
end function has_diag_yaml_file
#endif
!> \brief Logical function to determine if the variable var_ids has been allocated or associated
!! \return .True. if var_ids exists .False. if var_ids has not been set
pure logical function has_var_ids (obj)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  has_var_ids = allocated(obj%var_ids)
end function has_var_ids
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TODO
!!> \brief Returns a copy of the value of diag_yaml_file
!!! \return A copy of diag_yaml_file
!#ifdef use_yaml
!pure function get_diag_yaml_file (obj) result (res)
!  class(fmsDiagFile_type), intent(in) :: obj !< The file object
!  type(diagYamlFiles_type) :: res
!  res = obj%diag_yaml_file
!end function get_diag_yaml_file
!#endif
!> \brief Returns a copy of the value of file_metadata_from_model
!! \return A copy of file_metadata_from_model
pure function get_file_metadata_from_model (obj) result (res)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  character(len=:), dimension(:), allocatable :: res
  res = obj%file_metadata_from_model
end function get_file_metadata_from_model
!> \brief Returns a copy of the value of var_ids
!! \return A copy of var_ids
pure function get_var_ids (obj) result (res)
  class(fmsDiagFile_type), intent(in) :: obj !< The file object
  integer, dimension(:), allocatable :: res
  allocate(res(size(obj%var_ids)))
  res = obj%var_ids
end function get_var_ids
!!!!!!!!! Functions from diag_yaml_file
#ifdef use_yaml
!> \brief Returns a copy of file_fname from the yaml object
!! \return Copy of file_fname
pure function get_file_fname (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 character (len=:), allocatable :: res
  res = obj%diag_yaml_file%get_file_fname()
end function get_file_fname
!> \brief Returns a copy of file_frequnit from the yaml object
!! \return Copy of file_frequnit
pure function get_file_frequnit (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 integer :: res
  res = obj%diag_yaml_file%get_file_frequnit()
end function get_file_frequnit
!> \brief Returns a copy of file_freq from the yaml object
!! \return Copy of file_freq
pure function get_file_freq (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 integer :: res
  res = obj%diag_yaml_file%get_file_freq()
end function get_file_freq
!> \brief Returns a copy of file_timeunit from the yaml object
!! \return Copy of file_timeunit
pure function get_file_timeunit (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 integer :: res
  res = obj%diag_yaml_file%get_file_timeunit()
end function get_file_timeunit
!> \brief Returns a copy of file_unlimdim from the yaml object
!! \return Copy of file_unlimdim
pure function get_file_unlimdim (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 character (len=:), allocatable :: res
  res = obj%diag_yaml_file%get_file_unlimdim()
end function get_file_unlimdim
!! TODO - get functions for sub region stuff
!> \brief Returns a copy of file_sub_region from the yaml object
!! \return Copy of file_sub_region
!pure function get_file_sub_region (obj) result(res)
! class(fmsDiagFile_type), intent(in) :: obj !< The file object
! integer :: res
!  res = obj%diag_yaml_file%get_file_sub_region()
!end function get_file_sub_region
!> \brief Returns a copy of file_new_file_freq from the yaml object
!! \return Copy of file_new_file_freq
pure function get_file_new_file_freq (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 integer :: res
  res = obj%diag_yaml_file%get_file_new_file_freq()
end function get_file_new_file_freq
!> \brief Returns a copy of file_new_file_freq_units from the yaml object
!! \return Copy of file_new_file_freq_units
pure function get_file_new_file_freq_units (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 integer :: res
  res = obj%diag_yaml_file%get_file_new_file_freq_units()
end function get_file_new_file_freq_units
!> \brief Returns a copy of file_start_time from the yaml object
!! \return Copy of file_start_time
pure function get_file_start_time (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 character (len=:), allocatable :: res
  res = obj%diag_yaml_file%get_file_start_time()
end function get_file_start_time
!> \brief Returns a copy of file_duration from the yaml object
!! \return Copy of file_duration
pure function get_file_duration (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 integer :: res
  res = obj%diag_yaml_file%get_file_duration()
end function get_file_duration
!> \brief Returns a copy of file_duration_units from the yaml object
!! \return Copy of file_duration_units
pure function get_file_duration_units (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 integer :: res
  res = obj%diag_yaml_file%get_file_duration_units()
end function get_file_duration_units
!> \brief Returns a copy of file_varlist from the yaml object
!! \return Copy of file_varlist
pure function get_file_varlist (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 character (len=:), allocatable, dimension(:) :: res
  res = obj%diag_yaml_file%get_file_varlist()
end function get_file_varlist
!> \brief Returns a copy of file_global_meta from the yaml object
!! \return Copy of file_global_meta
pure function get_file_global_meta (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 character (len=:), allocatable, dimension(:,:) :: res
  res = obj%diag_yaml_file%get_file_global_meta()
end function get_file_global_meta
!> \brief Checks if file_fname is allocated in the yaml object
!! \return true if file_fname is allocated
pure function has_file_fname (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_fname()
end function has_file_fname
!> \brief Checks if file_frequnit is allocated in the yaml object
!! \return true if file_frequnit is allocated
pure function has_file_frequnit (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_frequnit()
end function has_file_frequnit
!> \brief Checks if file_freq is allocated in the yaml object
!! \return true if file_freq is allocated
pure function has_file_freq (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_freq()
end function has_file_freq
!> \brief Checks if file_timeunit is allocated in the yaml object
!! \return true if file_timeunit is allocated
pure function has_file_timeunit (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_timeunit()
end function has_file_timeunit
!> \brief Checks if file_unlimdim is allocated in the yaml object
!! \return true if file_unlimdim is allocated
pure function has_file_unlimdim (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_unlimdim()
end function has_file_unlimdim
!> \brief Checks if file_sub_region is allocated in the yaml object
!! \return true if file_sub_region is allocated
pure function has_file_sub_region (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_sub_region()
end function has_file_sub_region
!> \brief Checks if file_new_file_freq is allocated in the yaml object
!! \return true if file_new_file_freq is allocated
pure function has_file_new_file_freq (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_new_file_freq()
end function has_file_new_file_freq
!> \brief Checks if file_new_file_freq_units is allocated in the yaml object
!! \return true if file_new_file_freq_units is allocated
pure function has_file_new_file_freq_units (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_new_file_freq_units()
end function has_file_new_file_freq_units
!> \brief Checks if file_start_time is allocated in the yaml object
!! \return true if file_start_time is allocated
pure function has_file_start_time (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_start_time()
end function has_file_start_time
!> \brief Checks if file_duration is allocated in the yaml object
!! \return true if file_duration is allocated
pure function has_file_duration (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_duration()
end function has_file_duration
!> \brief Checks if file_duration_units is allocated in the yaml object
!! \return true if file_duration_units is allocated
pure function has_file_duration_units (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_duration_units()
end function has_file_duration_units
!> \brief Checks if file_varlist is allocated in the yaml object
!! \return true if file_varlist is allocated
pure function has_file_varlist (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_varlist()
end function has_file_varlist
!> \brief Checks if file_global_meta is allocated in the yaml object
!! \return true if file_global_meta is allocated
pure function has_file_global_meta (obj) result(res)
 class(fmsDiagFile_type), intent(in) :: obj !< The file object
 logical :: res
  res = obj%diag_yaml_file%has_file_global_meta()
end function has_file_global_meta

!> @brief Set the domain and the type_of_domain for a file
!> @details This subroutine is going to be called once by every variable in the file
!! in register_diag_field. It will update the domain and the type_of_domain if needed and verify that
!! all the variables are in the same domain
subroutine set_file_domain(obj, domain, type_of_domain)
  class(fmsDiagFile_type), intent(inout)       :: obj            !< The file object
  integer,                 INTENT(in)          :: type_of_domain !< fileobj_type to use
  CLASS(diagDomain_t),     INTENT(in), target  :: domain         !< Domain

  !! If this a sub_regional, don't do anything here
  if (obj%type_of_domain .eq. SUB_REGIONAL) return

  if (type_of_domain .ne. obj%type_of_domain) then
  !! If the current type_of_domain in the file obj is not the same as the variable calling this subroutine

    if (type_of_domain .eq. NO_DOMAIN .or. obj%type_of_domain .eq. NO_DOMAIN) then
    !! If they are not the same then one of them can be NO_DOMAIN
    !! (i.e a file can have variables that are not domain decomposed and variables that are)

      if (type_of_domain .ne. NO_DOMAIN) then
      !! Update the file's type_of_domain and domain if needed
        obj%type_of_domain = type_of_domain
        obj%domain => domain
      endif

    else
    !! If they are not the same and of them is not NO_DOMAIN, then crash because the variables don't have the
    !! same domain (i.e a file has a variable is that in a 2D domain and one that is in a UG domain)

      call mpp_error(FATAL, "The file: "//obj%get_file_fname()//" has variables that are not in the same domain")
    endif
  endif

end subroutine set_file_domain

!> @brief Loops through a variable's axis_ids and adds them to the FMSDiagFile object if they don't exist
subroutine add_axes(obj, axis_ids)
  class(fmsDiagFile_type), intent(inout)       :: obj            !< The file object
  integer,                 INTENT(in)          :: axis_ids(:)    !< Array of axes_ids

  integer :: i, j !< For do loops

  do i = 1, size(axis_ids)
    do j = 1, obj%number_of_axis
      !> Check if the axis already exists, if it does leave this do loop
      if (axis_ids(i) .eq. obj%axis_ids(j)) exit
    enddo

    !> If the axis does not exist add it to the list
    obj%number_of_axis = obj%number_of_axis + 1
    obj%axis_ids(obj%number_of_axis) = axis_ids(i)

    !> If this is a sub_regional file, set up the sub_axes
    !> TO DO:
    !!
  enddo

end subroutine add_axes

!> @brief adds the start time to the fileobj
!! @note This should be called from the register field calls. It can be called multiple times (one for each variable)
!! So it needs to make sure that the start_time is the same for each variable. The initial value is the base_time
subroutine add_start_time(obj, start_time)
  class(fmsDiagFile_type), intent(inout)       :: obj            !< The file object
  TYPE(time_type),         intent(in)          :: start_time     !< Start time to add to the fileobj

  if (obj%start_time .ne. get_base_time()) then
    !> If the obj%start_time is not equal to the base_time from the diag_table
    !! obj%start_time was already updated so make sure it is the same or error out
    if (obj%start_time .ne. start_time)&
      call mpp_error(FATAL, "The variables associated with the file:"//obj%get_file_fname()//" have"&
      &" different start_time")
  else
    !> If the obj%start_time is equal to the base_time,
    !! simply update it with the start_time and set up the *_output variables
    obj%start_time = start_time
    obj%last_output = start_time
    obj%next_output = diag_time_inc(start_time, obj%get_file_freq(), obj%get_file_frequnit())
    obj%next_next_output = diag_time_inc(obj%next_output, obj%get_file_freq(), obj%get_file_frequnit())
  endif

end subroutine
#endif
end module fms_diag_file_object_mod
