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
!> @file
!> @ingroup fms2_io

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

!> @brief Registers a regional 2D variable and stores the information needed
subroutine register_restart_region_2d(fileobj, variable_name, vdata, indices, global_size, pelist, is_root_pe, &
                         x_halo, y_halo, jshift, ishift, is_optional)
  type(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in), target :: vdata !< Pointer to variable data
  integer, dimension(4), intent(in)                 :: indices!< Indices for the halo region for the variable
                                                              !! (starting x, ending x, starting y, ending y)
  integer, dimension(:), intent(in)                 :: global_size !> Global_size(1:2) == size of the variable in (x,y)
  integer, dimension(:), intent(in)                 :: pelist !> List of pelist that have the data for
                                                              !! the variable
  logical, intent(in)                               :: is_root_pe !> Flag indicating if this is the root_pe
                                                                  !! from the pelist
  integer, intent(in), optional :: x_halo !< Number of halos in x
  integer, intent(in), optional :: y_halo !< Number of halos in y
  integer, intent(in), optional :: jshift !< Shift in the x axis (from center)
  integer, intent(in), optional :: ishift !< Shift in the y axis (from center)
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  integer :: current_var !> Current number of registered variables
  character(len=7) :: dimnames(4) !> Names of the dimensions for the variable_name

  !> If you reading the file and the variable is optional and not in the file, leave
  if (fileobj%is_readonly .and. present(is_optional)) then
    if (is_optional .and. .not. variable_exists(fileobj, variable_name)) then
      return
    endif
  endif

  !> Register the axis for the variable dimension and get the names of the dimensions
  call set_dimensions(fileobj, fileobj%bc_dimensions, dimnames, global_size)
  call netcdf_add_restart_variable_2d(fileobj, variable_name, vdata, dimnames(1:3), is_optional)

  !> Save the boundary conditions information to the bc_info type for the given variable
  current_var = fileobj%num_restart_vars
  fileobj%restart_vars(current_var)%is_bc_variable = .true.
  fileobj%restart_vars(current_var)%bc_info%is_root_pe = is_root_pe

  allocate(fileobj%restart_vars(current_var)%bc_info%indices(4))
  fileobj%restart_vars(current_var)%bc_info%indices = indices

  allocate(fileobj%restart_vars(current_var)%bc_info%global_size(2))
  fileobj%restart_vars(current_var)%bc_info%global_size = global_size

  if (any(fileobj%io_root .eq. pelist(:))) then
      allocate(fileobj%restart_vars(current_var)%bc_info%pelist(size(pelist)))
      fileobj%restart_vars(current_var)%bc_info%pelist = pelist
      fileobj%restart_vars(current_var)%bc_info%data_on_file_root = .true.
  else
      allocate(fileobj%restart_vars(current_var)%bc_info%pelist(size(pelist)+1))
      fileobj%restart_vars(current_var)%bc_info%pelist = (/ fileobj%io_root, pelist(:) /)
      fileobj%restart_vars(current_var)%bc_info%data_on_file_root = .false.
  endif

  fileobj%restart_vars(current_var)%bc_info%x_halo = 0
  fileobj%restart_vars(current_var)%bc_info%y_halo = 0
  fileobj%restart_vars(current_var)%bc_info%ishift = 0
  fileobj%restart_vars(current_var)%bc_info%jshift = 0

  if (present(x_halo)) fileobj%restart_vars(current_var)%bc_info%x_halo = x_halo
  if (present(y_halo)) fileobj%restart_vars(current_var)%bc_info%y_halo = y_halo
  if (present(ishift)) fileobj%restart_vars(current_var)%bc_info%ishift = ishift
  if (present(jshift)) fileobj%restart_vars(current_var)%bc_info%jshift = jshift

end subroutine register_restart_region_2d

!> @brief Registers a regional 3D variable and stores the information needed
subroutine register_restart_region_3d(fileobj, variable_name, vdata, indices, global_size, pelist, is_root_pe, &
                         x_halo, y_halo, jshift, ishift, is_optional)
  type(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in), target :: vdata !< Pointer to variable data
  integer, dimension(4), intent(in)                 :: indices !< Indices for the halo region for the variable
                                                               !! (starting x, ending x, starting y, ending y)
  integer, dimension(:), intent(in)                 :: global_size !> Global_size(1:3) == size of the variable in (x,y,z)
  integer, dimension(:), intent(in)                 :: pelist !> List of pelist that have the data for
                                                              !! the variable
  logical, intent(in)                               :: is_root_pe !> Flag indicating if this is the root_pe
                                                                  !! from the pelist
  integer, intent(in), optional :: x_halo !< Number of halos in x
  integer, intent(in), optional :: y_halo !< Number of halos in y
  integer, intent(in), optional :: jshift !< Shift in the y axis (from center)
  integer, intent(in), optional :: ishift !< Shift in the x axis (from center)
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  integer :: current_var !> Current number of registered variables
  character(len=7) :: dimnames(4) !> Names of the dimensions for the variable_name

  !> If you reading the file and the variable is optional and not in the file, leave
  if (fileobj%is_readonly .and. present(is_optional)) then
    if (is_optional .and. .not. variable_exists(fileobj, variable_name)) then
      return
    endif
  endif

  !> Register the axis for the variable dimension and get the names of the dimensions
  call set_dimensions(fileobj, fileobj%bc_dimensions, dimnames, global_size)
  call netcdf_add_restart_variable_3d(fileobj, variable_name, vdata, dimnames(1:4), is_optional)

  !> Save the boundary conditions information to the bc_info type for the given variable
  current_var = fileobj%num_restart_vars
  fileobj%restart_vars(current_var)%is_bc_variable = .true.
  fileobj%restart_vars(current_var)%bc_info%is_root_pe = is_root_pe

  allocate(fileobj%restart_vars(current_var)%bc_info%indices(4))
  fileobj%restart_vars(current_var)%bc_info%indices = indices

  allocate(fileobj%restart_vars(current_var)%bc_info%global_size(3))
  fileobj%restart_vars(current_var)%bc_info%global_size = global_size

  if (any(fileobj%io_root .eq. pelist(:))) then
      allocate(fileobj%restart_vars(current_var)%bc_info%pelist(size(pelist)))
      fileobj%restart_vars(current_var)%bc_info%pelist = pelist
      fileobj%restart_vars(current_var)%bc_info%data_on_file_root = .true.
  else
      allocate(fileobj%restart_vars(current_var)%bc_info%pelist(size(pelist)+1))
      fileobj%restart_vars(current_var)%bc_info%pelist = (/ fileobj%io_root, pelist(:) /)
      fileobj%restart_vars(current_var)%bc_info%data_on_file_root = .false.
  endif

  fileobj%restart_vars(current_var)%bc_info%x_halo = 0
  fileobj%restart_vars(current_var)%bc_info%y_halo = 0
  fileobj%restart_vars(current_var)%bc_info%ishift = 0
  fileobj%restart_vars(current_var)%bc_info%jshift = 0

  if (present(x_halo)) fileobj%restart_vars(current_var)%bc_info%x_halo = x_halo
  if (present(y_halo)) fileobj%restart_vars(current_var)%bc_info%y_halo = y_halo
  if (present(ishift)) fileobj%restart_vars(current_var)%bc_info%ishift = ishift
  if (present(jshift)) fileobj%restart_vars(current_var)%bc_info%jshift = jshift

end subroutine register_restart_region_3d

!> @brief Updates the dimension information (bc_dimensions) for the give fileobj and outputs the
!! dimensions name of the variable
subroutine set_dimensions(fileobj, bc_dimensions, dimnames, global_size)
  type(FmsNetcdfFile_t), intent(inout) :: fileobj !> File object.
  type(dimension_information), intent(inout) :: bc_dimensions !> Dimension information for the file
  character(len=7), intent(inout) :: dimnames(4) !> List of the dimension names for the variable
                                                 !! dimnames(1): Name of the x dimension
                                                 !! dimnames(2): Name of the y dimension
                                                 !! dimnames(3): Name of the z dimension (or time if 2d)
                                                 !! dimnames(4): Name of the Time dimension (if 3d variable)
  integer, dimension(:), intent(in) :: global_size !> Size of the variable for each dimension
                                                   !! global_size(1): Size of the x dimension
                                                   !! global_size(2): Size of the y dimension
                                                   !! global_size(3): Size of the z dimension (if it exists)
  character(len=1) :: suffix   !> Axis number converted to a string
  character(len=7) :: axisname !> Name of the axis
  integer          :: i, found_index !> no description

  if (.not. dimension_exists(fileobj, "Time")) call netcdf_add_dimension(fileobj, "Time", unlimited)
  dimnames(size(global_size)+1) = "Time"

  !> Check to see if a x dimension with the correct size is already registered
  found_index = 0
  do i = 1, bc_dimensions%cur_dim_len(1)
     if (bc_dimensions%xlen(i) .eq. global_size(1)) then
        found_index= i
        write(suffix,'(I1)') i
        axisname = 'xaxis_'//suffix
        dimnames(1) = axisname
        exit
      endif
  enddo

  !> If the x dimension doesn't already exist, register it
  if (found_index .eq. 0) then
     bc_dimensions%cur_dim_len(1) = bc_dimensions%cur_dim_len(1) + 1
     bc_dimensions%xlen(bc_dimensions%cur_dim_len(1)) = global_size(1)
     write(suffix,'(I1)') bc_dimensions%cur_dim_len(1)
     axisname = 'xaxis_'//suffix
     call netcdf_add_dimension(fileobj, axisname, global_size(1))
     dimnames(1) = axisname
  endif

  !> Check to see if a y dimension with the correct size is already registered
  found_index = 0
  do i = 1, bc_dimensions%cur_dim_len(2)
     if (bc_dimensions%ylen(i) .eq. global_size(2)) then
        found_index= i
        write(suffix,'(I1)') i
        axisname = 'yaxis_'//suffix
        dimnames(2) = axisname
        exit
      endif
  enddo

  !> If the y dimension doesn't already exist, register it
  if (found_index .eq. 0) then
     bc_dimensions%cur_dim_len(2) = bc_dimensions%cur_dim_len(2) + 1
     bc_dimensions%ylen(bc_dimensions%cur_dim_len(2)) = global_size(2)
     write(suffix,'(I1)') bc_dimensions%cur_dim_len(2)
     axisname = 'yaxis_'//suffix
     call netcdf_add_dimension(fileobj, axisname, global_size(2))
     dimnames(2) = axisname
  endif

  !> If the variable is 3d:
  if (size(global_size) .eq. 3) then
     !> Check to see if a z dimension with the correct size is already registered
     found_index = 0
     do i = 1, bc_dimensions%cur_dim_len(3)
        if (bc_dimensions%zlen(i) .eq. global_size(3)) then
           found_index= i
           write(suffix,'(I1)') i
           axisname = 'zaxis_'//suffix
           dimnames(3) = axisname
           exit
        endif
     enddo

     !> If the z dimension doesn't already exist, register it
     if (found_index .eq. 0) then
        bc_dimensions%cur_dim_len(3) = bc_dimensions%cur_dim_len(3) + 1
        bc_dimensions%zlen(bc_dimensions%cur_dim_len(3)) = global_size(3)
        write(suffix,'(I1)') bc_dimensions%cur_dim_len(3)
        axisname = 'zaxis_'//suffix
        call netcdf_add_dimension(fileobj, axisname, global_size(3))
        dimnames(3) = axisname
      endif
  endif
end subroutine set_dimensions
