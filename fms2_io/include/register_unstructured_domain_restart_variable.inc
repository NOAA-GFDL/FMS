!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> @file
!> @brief Adds restart variable for unstructured domains for different dimension to be used in
!! the @ref register_restart_field interface

!> @addtogroup fms2_io_mod
!> @{

!> @brief Add a domain decomposed variable.
subroutine register_unstructured_domain_restart_variable_0d(fileobj, variable_name, vdata, &
                                                            dimensions, is_optional,  &
                                                            chunksizes)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.
  integer, intent(in), optional :: chunksizes(:) !< netcdf chunksize to use for this variable
                                                 !! This feature is only
                                                 !! available for netcdf4 files

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional, &
                                    chunksizes=chunksizes)
end subroutine register_unstructured_domain_restart_variable_0d


!> @brief Add a domain decomposed variable.
subroutine register_unstructured_domain_restart_variable_1d(fileobj, variable_name, vdata, &
                                                            dimensions, is_optional,  &
                                                            chunksizes)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.
  integer, intent(in), optional :: chunksizes(:) !< netcdf chunksize to use for this variable
                                                 !! This feature is only
                                                 !! available for netcdf4 files

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional, &
                                    chunksizes=chunksizes)
end subroutine register_unstructured_domain_restart_variable_1d


!> @brief Add a domain decomposed variable.
subroutine register_unstructured_domain_restart_variable_2d(fileobj, variable_name, vdata, &
                                                            dimensions, is_optional,  &
                                                            chunksizes )

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.
  integer, intent(in), optional :: chunksizes(:) !< netcdf chunksize to use for this variable
                                                 !! This feature is only
                                                 !! available for netcdf4 files

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional, &
                                    chunksizes=chunksizes)
end subroutine register_unstructured_domain_restart_variable_2d


!> @brief Add a domain decomposed variable.
subroutine register_unstructured_domain_restart_variable_3d(fileobj, variable_name, vdata, &
                                                            dimensions, is_optional,  &
                                                            chunksizes)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.
  integer, intent(in), optional :: chunksizes(:) !< netcdf chunksize to use for this variable
                                                 !! This feature is only
                                                 !! available for netcdf4 files

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional, &
                                    chunksizes=chunksizes)
end subroutine register_unstructured_domain_restart_variable_3d


!> @brief Add a domain decomposed variable.
subroutine register_unstructured_domain_restart_variable_4d(fileobj, variable_name, vdata, &
                                                            dimensions, is_optional,  &
                                                            chunksizes)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.
  integer, intent(in), optional :: chunksizes(:) !< netcdf chunksize to use for this variable
                                                 !! This feature is only
                                                 !! available for netcdf4 files

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional, &
                                    chunksizes=chunksizes)
end subroutine register_unstructured_domain_restart_variable_4d


!> @brief Add a domain decomposed variable.
subroutine register_unstructured_domain_restart_variable_5d(fileobj, variable_name, vdata, &
                                                            dimensions, is_optional,  &
                                                            chunksizes)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.
  integer, intent(in), optional :: chunksizes(:) !< netcdf chunksize to use for this variable
                                                 !! This feature is only
                                                 !! available for netcdf4 files

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional, &
                                    chunksizes=chunksizes)
end subroutine register_unstructured_domain_restart_variable_5d
!> @}
