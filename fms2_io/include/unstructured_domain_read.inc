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
!> @brief Routines for reading unstructured domain data, used in @ref read_data interface

!> @addtogroup fms2_io_mod
!> @{

!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_read_0d(fileobj, &
                                                 variable_name, &
                                                 buf, &
                                                 unlim_dim_level, &
                                                 corner, &

                                                 broadcast)
  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(inout) :: buf !< Array that the data
                                                       !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be read from.

  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  call netcdf_read_data(fileobj, variable_name, buf, unlim_dim_level, corner, &

                        broadcast)
end subroutine unstructured_domain_read_0d
!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_read_1d(fileobj, &
                                                 variable_name, &
                                                 buf, &
                                                 unlim_dim_level, &
                                                 corner, &
                                                 edge_lengths, &
                                                 broadcast)
  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(inout) :: buf !< Array that the data
                                                       !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                               !! elements that
                                                               !! will be read
                                                               !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  call netcdf_read_data(fileobj, variable_name, buf, unlim_dim_level, corner, &
                        edge_lengths, &
                        broadcast)
end subroutine unstructured_domain_read_1d
!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_read_2d(fileobj, &
                                                 variable_name, &
                                                 buf, &
                                                 unlim_dim_level, &
                                                 corner, &
                                                 edge_lengths, &
                                                 broadcast)
  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(inout) :: buf !< Array that the data
                                                       !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                               !! elements that
                                                               !! will be read
                                                               !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  call netcdf_read_data(fileobj, variable_name, buf, unlim_dim_level, corner, &
                        edge_lengths, &
                        broadcast)
end subroutine unstructured_domain_read_2d
!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_read_3d(fileobj, &
                                                 variable_name, &
                                                 buf, &
                                                 unlim_dim_level, &
                                                 corner, &
                                                 edge_lengths, &
                                                 broadcast)
  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(inout) :: buf !< Array that the data
                                                       !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                               !! elements that
                                                               !! will be read
                                                               !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  call netcdf_read_data(fileobj, variable_name, buf, unlim_dim_level, corner, &
                        edge_lengths, &
                        broadcast)
end subroutine unstructured_domain_read_3d
!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_read_4d(fileobj, &
                                                 variable_name, &
                                                 buf, &
                                                 unlim_dim_level, &
                                                 corner, &
                                                 edge_lengths, &
                                                 broadcast)
  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(inout) :: buf !< Array that the data
                                                       !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                               !! elements that
                                                               !! will be read
                                                               !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  call netcdf_read_data(fileobj, variable_name, buf, unlim_dim_level, corner, &
                        edge_lengths, &
                        broadcast)
end subroutine unstructured_domain_read_4d
!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_read_5d(fileobj, &
                                                 variable_name, &
                                                 buf, &
                                                 unlim_dim_level, &
                                                 corner, &
                                                 edge_lengths, &
                                                 broadcast)
  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(inout) :: buf !< Array that the data
                                                       !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                               !! elements that
                                                               !! will be read
                                                               !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  call netcdf_read_data(fileobj, variable_name, buf, unlim_dim_level, corner, &
                        edge_lengths, &
                        broadcast)
end subroutine unstructured_domain_read_5d
!> @}
