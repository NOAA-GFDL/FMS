!> @brief Add a restart variable to a file.
subroutine register_restart_variable_0d(fileobj, &
                                                  variable_name, &
                                                  vdata, &
                                                  dimensions)
    type(FmsNetcdfFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*), intent(in),target :: vdata !< Pointer to
                                                              !! variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_restart_variable_0d
!> @brief Add a restart variable to a file.
subroutine register_restart_variable_1d(fileobj, &
                                                  variable_name, &
                                                  vdata, &
                                                  dimensions)
    type(FmsNetcdfFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:), intent(in),target :: vdata !< Pointer to
                                                              !! variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_restart_variable_1d
!> @brief Add a restart variable to a file.
subroutine register_restart_variable_2d(fileobj, &
                                                  variable_name, &
                                                  vdata, &
                                                  dimensions)
    type(FmsNetcdfFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:), intent(in),target :: vdata !< Pointer to
                                                              !! variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_restart_variable_2d
!> @brief Add a restart variable to a file.
subroutine register_restart_variable_3d(fileobj, &
                                                  variable_name, &
                                                  vdata, &
                                                  dimensions)
    type(FmsNetcdfFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:), intent(in),target :: vdata !< Pointer to
                                                              !! variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_restart_variable_3d
!> @brief Add a restart variable to a file.
subroutine register_restart_variable_4d(fileobj, &
                                                  variable_name, &
                                                  vdata, &
                                                  dimensions)
    type(FmsNetcdfFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:), intent(in),target :: vdata !< Pointer to
                                                              !! variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_restart_variable_4d
!> @brief Add a restart variable to a file.
subroutine register_restart_variable_5d(fileobj, &
                                                  variable_name, &
                                                  vdata, &
                                                  dimensions)
    type(FmsNetcdfFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:,:), intent(in),target :: vdata !< Pointer to
                                                              !! variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_restart_variable_5d
