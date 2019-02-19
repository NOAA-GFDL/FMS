!> @brief Add a compressed restart variable.
subroutine register_compressed_restart_variable_0d(fileobj, &
                                                             variable_name, &
                                                             vdata, &
                                                             dimensions)
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*), intent(in),target :: vdata !< Variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call append_compressed_variable(fileobj, &
                                    variable_name, &
                                    dimensions)
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_compressed_restart_variable_0d
!> @brief Add a compressed restart variable.
subroutine register_compressed_restart_variable_1d(fileobj, &
                                                             variable_name, &
                                                             vdata, &
                                                             dimensions)
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:), intent(in),target :: vdata !< Variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call append_compressed_variable(fileobj, &
                                    variable_name, &
                                    dimensions)
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_compressed_restart_variable_1d
!> @brief Add a compressed restart variable.
subroutine register_compressed_restart_variable_2d(fileobj, &
                                                             variable_name, &
                                                             vdata, &
                                                             dimensions)
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:), intent(in),target :: vdata !< Variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call append_compressed_variable(fileobj, &
                                    variable_name, &
                                    dimensions)
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_compressed_restart_variable_2d
!> @brief Add a compressed restart variable.
subroutine register_compressed_restart_variable_3d(fileobj, &
                                                             variable_name, &
                                                             vdata, &
                                                             dimensions)
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:), intent(in),target :: vdata !< Variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call append_compressed_variable(fileobj, &
                                    variable_name, &
                                    dimensions)
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_compressed_restart_variable_3d
!> @brief Add a compressed restart variable.
subroutine register_compressed_restart_variable_4d(fileobj, &
                                                             variable_name, &
                                                             vdata, &
                                                             dimensions)
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:), intent(in),target :: vdata !< Variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call append_compressed_variable(fileobj, &
                                    variable_name, &
                                    dimensions)
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_compressed_restart_variable_4d
!> @brief Add a compressed restart variable.
subroutine register_compressed_restart_variable_5d(fileobj, &
                                                             variable_name, &
                                                             vdata, &
                                                             dimensions)
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:,:), intent(in),target :: vdata !< Variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call append_compressed_variable(fileobj, &
                                    variable_name, &
                                    dimensions)
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
end subroutine register_compressed_restart_variable_5d
