include(`macros.m4')dnl
!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
`subroutine compressed_read_'NUM_DIMS`d(fileobj, &'
                                        variable_name, &
                                        cdata, &
                                        unlim_dim_level)

    !Inputs/outputs
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dim_declare(NUM_DIMS) intent(inout) :: cdata !< Buffer where the
                                                      !! read in data will
                                                      !! be stored.
    integer,intent(in),optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          cdata, &
                          unlim_dim_level=unlim_dim_level)
`end subroutine compressed_read_'NUM_DIMS`d'
