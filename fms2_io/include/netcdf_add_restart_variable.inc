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

!> @brief Add a restart variable to a netcdf file.
subroutine netcdf_add_restart_variable_0d(fileobj, variable_name, vdata, dimensions, &
                                          is_optional)

  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(1), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  character(len=8) :: buf

  if (fileobj%is_readonly .and. present(is_optional)) then
    if (is_optional .and. .not. variable_exists(fileobj, variable_name)) then
      return
    endif
  endif
  call add_restart_var_to_array(fileobj, variable_name)
  fileobj%restart_vars(fileobj%num_restart_vars)%data0d => vdata
  if (.not. fileobj%is_readonly) then
    call get_data_type_string(vdata, buf)
    if (present(dimensions)) then
      if (.not. is_dimension_unlimited(fileobj, dimensions(1), .true.)) then
        call error("a scalar input variable can only have an unlimited dimension.")
      endif
    endif
    call netcdf_add_variable(fileobj, variable_name, buf, dimensions)
  endif
end subroutine netcdf_add_restart_variable_0d


!> @brief Wrapper to distinguish interfaces.
subroutine netcdf_add_restart_variable_0d_wrap(fileobj, variable_name, vdata, &
                                               dimensions, is_optional)
  type(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(1), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
end subroutine netcdf_add_restart_variable_0d_wrap


!> @brief Add a restart variable to a netcdf file.
subroutine netcdf_add_restart_variable_1d(fileobj, variable_name, vdata, &
                                          dimensions, is_optional)
  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  character(len=8) :: buf
  integer :: ndims
  integer :: vdata_rank

  if (fileobj%is_readonly .and. present(is_optional)) then
    if (is_optional .and. .not. variable_exists(fileobj, variable_name)) then
      return
    endif
  endif
  call add_restart_var_to_array(fileobj, variable_name)
  fileobj%restart_vars(fileobj%num_restart_vars)%data1d => vdata
  if (.not. fileobj%is_readonly) then
    call get_data_type_string(vdata, buf)
    if (.not. present(dimensions)) then
      call error("dimension names required if the file is not read-only.")
    endif
    ndims = size(dimensions)
    vdata_rank = size(shape(vdata))
    if (ndims .eq. vdata_rank+1) then
      if (.not. is_dimension_unlimited(fileobj, dimensions(ndims), .true.)) then
        call error("the slowest dimension must be unlimited.")
      endif
    elseif (ndims .ne. vdata_rank) then
      call error("rank mismatch between vdata and dimensions arrays.")
    endif
    call netcdf_add_variable(fileobj, variable_name, buf, dimensions)
  endif
end subroutine netcdf_add_restart_variable_1d


!> @brief Wrapper to distinguish interfaces.
subroutine netcdf_add_restart_variable_1d_wrap(fileobj, variable_name, vdata, &
                                               dimensions, is_optional)
  type(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
end subroutine netcdf_add_restart_variable_1d_wrap


!> @brief Add a restart variable to a netcdf file.
subroutine netcdf_add_restart_variable_2d(fileobj, variable_name, vdata, &
                                          dimensions, is_optional)
  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  character(len=8) :: buf
  integer :: ndims
  integer :: vdata_rank

  if (fileobj%is_readonly .and. present(is_optional)) then
    if (is_optional .and. .not. variable_exists(fileobj, variable_name)) then
      return
    endif
  endif
  call add_restart_var_to_array(fileobj, variable_name)
  fileobj%restart_vars(fileobj%num_restart_vars)%data2d => vdata
  if (.not. fileobj%is_readonly) then
    call get_data_type_string(vdata, buf)
    if (.not. present(dimensions)) then
      call error("dimension names required if the file is not read-only.")
    endif
    ndims = size(dimensions)
    vdata_rank = size(shape(vdata))
    if (ndims .eq. vdata_rank+1) then
      if (.not. is_dimension_unlimited(fileobj, dimensions(ndims), .true.)) then
        call error("the slowest dimension must be unlimited.")
      endif
    elseif (ndims .ne. vdata_rank) then
      call error("rank mismatch between vdata and dimensions arrays.")
    endif
    call netcdf_add_variable(fileobj, variable_name, buf, dimensions)
  endif
end subroutine netcdf_add_restart_variable_2d


!> @brief Wrapper to distinguish interfaces.
subroutine netcdf_add_restart_variable_2d_wrap(fileobj, variable_name, vdata, &
                                               dimensions, is_optional)
  type(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
end subroutine netcdf_add_restart_variable_2d_wrap


!> @brief Add a restart variable to a netcdf file.
subroutine netcdf_add_restart_variable_3d(fileobj, variable_name, vdata, &
                                          dimensions, is_optional)
  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  character(len=8) :: buf
  integer :: ndims
  integer :: vdata_rank

  if (fileobj%is_readonly .and. present(is_optional)) then
    if (is_optional .and. .not. variable_exists(fileobj, variable_name)) then
      return
    endif
  endif
  call add_restart_var_to_array(fileobj, variable_name)
  fileobj%restart_vars(fileobj%num_restart_vars)%data3d => vdata
  if (.not. fileobj%is_readonly) then
    call get_data_type_string(vdata, buf)
    if (.not. present(dimensions)) then
      call error("dimension names required if the file is not read-only.")
    endif
    ndims = size(dimensions)
    vdata_rank = size(shape(vdata))
    if (ndims .eq. vdata_rank+1) then
      if (.not. is_dimension_unlimited(fileobj, dimensions(ndims), .true.)) then
        call error("the slowest dimension must be unlimited.")
      endif
    elseif (ndims .ne. vdata_rank) then
      call error("rank mismatch between vdata and dimensions arrays.")
    endif
    call netcdf_add_variable(fileobj, variable_name, buf, dimensions)
  endif
end subroutine netcdf_add_restart_variable_3d


!> @brief Wrapper to distinguish interfaces.
subroutine netcdf_add_restart_variable_3d_wrap(fileobj, variable_name, vdata, &
                                               dimensions, is_optional)
  type(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
end subroutine netcdf_add_restart_variable_3d_wrap


!> @brief Add a restart variable to a netcdf file.
subroutine netcdf_add_restart_variable_4d(fileobj, variable_name, vdata, &
                                          dimensions, is_optional)
  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  character(len=8) :: buf
  integer :: ndims
  integer :: vdata_rank

  if (fileobj%is_readonly .and. present(is_optional)) then
    if (is_optional .and. .not. variable_exists(fileobj, variable_name)) then
      return
    endif
  endif
  call add_restart_var_to_array(fileobj, variable_name)
  fileobj%restart_vars(fileobj%num_restart_vars)%data4d => vdata
  if (.not. fileobj%is_readonly) then
    call get_data_type_string(vdata, buf)
    if (.not. present(dimensions)) then
      call error("dimension names required if the file is not read-only.")
    endif
    ndims = size(dimensions)
    vdata_rank = size(shape(vdata))
    if (ndims .eq. vdata_rank+1) then
      if (.not. is_dimension_unlimited(fileobj, dimensions(ndims), .true.)) then
        call error("the slowest dimension must be unlimited.")
      endif
    elseif (ndims .ne. vdata_rank) then
      call error("rank mismatch between vdata and dimensions arrays.")
    endif
    call netcdf_add_variable(fileobj, variable_name, buf, dimensions)
  endif
end subroutine netcdf_add_restart_variable_4d


!> @brief Wrapper to distinguish interfaces.
subroutine netcdf_add_restart_variable_4d_wrap(fileobj, variable_name, vdata, &
                                               dimensions, is_optional)
  type(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
end subroutine netcdf_add_restart_variable_4d_wrap


!> @brief Add a restart variable to a netcdf file.
subroutine netcdf_add_restart_variable_5d(fileobj, variable_name, vdata, &
                                          dimensions, is_optional)
  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  character(len=8) :: buf
  integer :: ndims
  integer :: vdata_rank

  if (fileobj%is_readonly .and. present(is_optional)) then
    if (is_optional .and. .not. variable_exists(fileobj, variable_name)) then
      return
    endif
  endif
  call add_restart_var_to_array(fileobj, variable_name)
  fileobj%restart_vars(fileobj%num_restart_vars)%data5d => vdata
  if (.not. fileobj%is_readonly) then
    call get_data_type_string(vdata, buf)
    if (.not. present(dimensions)) then
      call error("dimension names required if the file is not read-only.")
    endif
    ndims = size(dimensions)
    vdata_rank = size(shape(vdata))
    if (ndims .eq. vdata_rank+1) then
      if (.not. is_dimension_unlimited(fileobj, dimensions(ndims), .true.)) then
        call error("the slowest dimension must be unlimited.")
      endif
    elseif (ndims .ne. vdata_rank) then
      call error("rank mismatch between vdata and dimensions arrays.")
    endif
    call netcdf_add_variable(fileobj, variable_name, buf, dimensions)
  endif
end subroutine netcdf_add_restart_variable_5d


!> @brief Wrapper to distinguish interfaces.
subroutine netcdf_add_restart_variable_5d_wrap(fileobj, variable_name, vdata, &
                                               dimensions, is_optional)
  type(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in), target :: vdata !< Pointer to
                                                              !! variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
end subroutine netcdf_add_restart_variable_5d_wrap
