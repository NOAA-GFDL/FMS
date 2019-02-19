include(`macros.m4')dnl
!> @brief Add a restart variable to a file.
`subroutine register_restart_variable_'NUM_DIMS`d(fileobj, &'
                                                  variable_name, &
                                                  vdata, &
                                                  dimensions)
    type(FmsNetcdfFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dim_declare(NUM_DIMS) intent(in),target :: vdata !< Pointer to
                                                              !! variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
`end subroutine register_restart_variable_'NUM_DIMS`d'
