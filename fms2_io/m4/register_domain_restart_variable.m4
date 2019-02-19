include(`macros.m4')dnl
!> @brief Add a domain decomposed variable.
`subroutine register_domain_restart_variable_'NUM_DIMS`d(fileobj, &'
                                                         variable_name, &
                                                         vdata, &
                                                         dimensions, &
                                                         domain_position)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dim_declare(NUM_DIMS) intent(in),target :: vdata !< Variable data.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    integer,intent(in),optional :: domain_position !< Domain position.

    call netcdf_add_restart_variable(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     dimensions)
    call append_domain_decomposed_variable(fileobj, &
                                           variable_name, &
                                           domain_position)
`end subroutine register_domain_restart_variable_'NUM_DIMS`d'
