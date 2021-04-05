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

!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_write_0d(fileobj, variable_name, variable_data, &
                                        unlim_dim_level, corner)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(in) :: variable_data !< Data that will be
                                         !! written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.

  call compressed_write(fileobj, variable_name, variable_data, unlim_dim_level, corner)
end subroutine unstructured_domain_write_0d


!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_write_1d(fileobj, variable_name, variable_data, &
                                        unlim_dim_level, corner, edge_lengths)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in) :: variable_data !< Data that will be
                                                      !! written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(1), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(1), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call compressed_write(fileobj, variable_name, variable_data, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)
end subroutine unstructured_domain_write_1d


!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_write_2d(fileobj, variable_name, variable_data, &
                                        unlim_dim_level, corner, edge_lengths)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in) :: variable_data !< Data that will be
                                                        !! written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(2), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(2), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call compressed_write(fileobj, variable_name, variable_data, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)
end subroutine unstructured_domain_write_2d


!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_write_3d(fileobj, variable_name, variable_data, &
                                        unlim_dim_level, corner, edge_lengths)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in) :: variable_data !< Data that will be
                                                          !! written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(3), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(3), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call compressed_write(fileobj, variable_name, variable_data, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)
end subroutine unstructured_domain_write_3d


!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_write_4d(fileobj, variable_name, variable_data, &
                                        unlim_dim_level, corner, edge_lengths)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in) :: variable_data !< Data that will be
                                                            !! written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(4), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(4), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call compressed_write(fileobj, variable_name, variable_data, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)
end subroutine unstructured_domain_write_4d


!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_domain_write_5d(fileobj, variable_name, variable_data, &
                                        unlim_dim_level, corner, edge_lengths)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in) :: variable_data !< Data that will be
                                                              !! written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(5), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(5), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call compressed_write(fileobj, variable_name, variable_data, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)
end subroutine unstructured_domain_write_5d
