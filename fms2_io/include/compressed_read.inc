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
!> @brief Compressed read routines for @ref read_data interface

!> @addtogroup fms2_io_mod
!> @{

!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_0d(fileobj, variable_name, cdata, unlim_dim_level, &
                              corner)

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), intent(inout) :: cdata !< Buffer where the
                                   !! read in data will
                                   !! be stored.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.


  call netcdf_read_data(fileobj, variable_name, cdata, &
                        unlim_dim_level=unlim_dim_level, corner=corner, &
                        broadcast=.true.)
end subroutine compressed_read_0d


!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_1d(fileobj, variable_name, cdata, unlim_dim_level, &
                              corner, edge_lengths)

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(inout) :: cdata !< Buffer where the
                                                 !! read in data will
                                                 !! be stored.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call netcdf_read_data(fileobj, variable_name, cdata, &
                        unlim_dim_level=unlim_dim_level, corner=corner, &
                        edge_lengths=edge_lengths, broadcast=.true.)
end subroutine compressed_read_1d


!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_2d(fileobj, variable_name, cdata, unlim_dim_level, &
                              corner, edge_lengths)

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(inout) :: cdata !< Buffer where the
                                                   !! read in data will
                                                   !! be stored.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call netcdf_read_data(fileobj, variable_name, cdata, &
                        unlim_dim_level=unlim_dim_level, corner=corner, &
                        edge_lengths=edge_lengths, broadcast=.true.)
end subroutine compressed_read_2d


!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_3d(fileobj, variable_name, cdata, unlim_dim_level, &
                              corner, edge_lengths)

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(inout) :: cdata !< Buffer where the
                                                     !! read in data will
                                                     !! be stored.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call netcdf_read_data(fileobj, variable_name, cdata, &
                        unlim_dim_level=unlim_dim_level, corner=corner, &
                        edge_lengths=edge_lengths, broadcast=.true.)
end subroutine compressed_read_3d


!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_4d(fileobj, variable_name, cdata, unlim_dim_level, &
                              corner, edge_lengths)

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(inout) :: cdata !< Buffer where the
                                                         !! read in data will
                                                         !! be stored.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call netcdf_read_data(fileobj, variable_name, cdata, &
                        unlim_dim_level=unlim_dim_level, corner=corner, &
                        edge_lengths=edge_lengths, broadcast=.true.)
end subroutine compressed_read_4d


!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_5d(fileobj, variable_name, cdata, unlim_dim_level, &
                              corner, edge_lengths)

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(inout) :: cdata !< Buffer where the
                                                         !! read in data will
                                                         !! be stored.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call netcdf_read_data(fileobj, variable_name, cdata, &
                        unlim_dim_level=unlim_dim_level, corner=corner &
                       ,edge_lengths=edge_lengths &
                        ,broadcast=.true.)
end subroutine compressed_read_5d
!> @}
