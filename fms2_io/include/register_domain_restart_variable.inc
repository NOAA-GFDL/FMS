!> @brief Add a domain decomposed variable.
subroutine register_domain_restart_variable_0d(fileobj, variable_name, vdata, &
                                               dimensions, domain_position, is_optional)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  integer, intent(in), optional :: domain_position !< Domain position.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
  call append_domain_decomposed_variable(fileobj, variable_name, domain_position)
end subroutine register_domain_restart_variable_0d


!> @brief Add a domain decomposed variable.
subroutine register_domain_restart_variable_1d(fileobj, variable_name, vdata, &
                                               dimensions, domain_position, is_optional)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  integer, intent(in), optional :: domain_position !< Domain position.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
  call append_domain_decomposed_variable(fileobj, variable_name, domain_position)
end subroutine register_domain_restart_variable_1d


!> @brief Add a domain decomposed variable.
subroutine register_domain_restart_variable_2d(fileobj, variable_name, vdata, &
                                               dimensions, domain_position, is_optional)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  integer, intent(in), optional :: domain_position !< Domain position.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
  call append_domain_decomposed_variable(fileobj, variable_name, domain_position)
end subroutine register_domain_restart_variable_2d


!> @brief Add a domain decomposed variable.
subroutine register_domain_restart_variable_3d(fileobj, variable_name, vdata, &
                                               dimensions, domain_position, is_optional)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  integer, intent(in), optional :: domain_position !< Domain position.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
  call append_domain_decomposed_variable(fileobj, variable_name, domain_position)
end subroutine register_domain_restart_variable_3d


!> @brief Add a domain decomposed variable.
subroutine register_domain_restart_variable_4d(fileobj, variable_name, vdata, &
                                               dimensions, domain_position, is_optional)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  integer, intent(in), optional :: domain_position !< Domain position.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions, is_optional)
  call append_domain_decomposed_variable(fileobj, variable_name, domain_position)
end subroutine register_domain_restart_variable_4d


!> @brief Add a domain decomposed variable.
subroutine register_domain_restart_variable_5d(fileobj, variable_name, vdata, &
                                               dimensions, domain_position, is_optional)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in), target :: vdata !< Variable data.
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.
  integer, intent(in), optional :: domain_position !< Domain position.
  logical, intent(in), optional :: is_optional !< Prevent errors in read-only files
                                               !! if a variable does not exist.

  call netcdf_add_restart_variable(fileobj, variable_name, vdata, dimensions,is_optional)
  call append_domain_decomposed_variable(fileobj, variable_name, domain_position)
end subroutine register_domain_restart_variable_5d
