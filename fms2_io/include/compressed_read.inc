!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_0d(fileobj, &
                                        variable_name, &
                                        cdata, &
                                        unlim_dim_level)

    !Inputs/outputs
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*), intent(inout) :: cdata !< Buffer where the
                                                      !! read in data will
                                                      !! be stored.
    integer,intent(in),optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          cdata, &
                          unlim_dim_level=unlim_dim_level)
end subroutine compressed_read_0d
!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_1d(fileobj, &
                                        variable_name, &
                                        cdata, &
                                        unlim_dim_level)

    !Inputs/outputs
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:), intent(inout) :: cdata !< Buffer where the
                                                      !! read in data will
                                                      !! be stored.
    integer,intent(in),optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          cdata, &
                          unlim_dim_level=unlim_dim_level)
end subroutine compressed_read_1d
!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_2d(fileobj, &
                                        variable_name, &
                                        cdata, &
                                        unlim_dim_level)

    !Inputs/outputs
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:), intent(inout) :: cdata !< Buffer where the
                                                      !! read in data will
                                                      !! be stored.
    integer,intent(in),optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          cdata, &
                          unlim_dim_level=unlim_dim_level)
end subroutine compressed_read_2d
!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_3d(fileobj, &
                                        variable_name, &
                                        cdata, &
                                        unlim_dim_level)

    !Inputs/outputs
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:), intent(inout) :: cdata !< Buffer where the
                                                      !! read in data will
                                                      !! be stored.
    integer,intent(in),optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          cdata, &
                          unlim_dim_level=unlim_dim_level)
end subroutine compressed_read_3d
!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_4d(fileobj, &
                                        variable_name, &
                                        cdata, &
                                        unlim_dim_level)

    !Inputs/outputs
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:), intent(inout) :: cdata !< Buffer where the
                                                      !! read in data will
                                                      !! be stored.
    integer,intent(in),optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          cdata, &
                          unlim_dim_level=unlim_dim_level)
end subroutine compressed_read_4d
!> @brief I/O domain reads in data from the netcdf file and broadcasts the
!!        data to the rest of the ranks.  This routine may only be used with
!!        variables that are "compressed".
subroutine compressed_read_5d(fileobj, &
                                        variable_name, &
                                        cdata, &
                                        unlim_dim_level)

    !Inputs/outputs
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:,:), intent(inout) :: cdata !< Buffer where the
                                                      !! read in data will
                                                      !! be stored.
    integer,intent(in),optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          cdata, &
                          unlim_dim_level=unlim_dim_level)
end subroutine compressed_read_5d
