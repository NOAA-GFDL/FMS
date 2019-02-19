!> @brief Read in data from a variable in a netcdf file.
subroutine read_data_0d(fileobj, &
                                  variable_name, &
                                  buf, &
                                  unlim_dim_level, &
                                  corner, &

                                  broadcast)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.

    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          buf, &
                          unlim_dim_level, &
                          corner, &

                          broadcast)
end subroutine read_data_0d
!> @brief Read in data from a variable in a netcdf file.
subroutine read_data_1d(fileobj, &
                                  variable_name, &
                                  buf, &
                                  unlim_dim_level, &
                                  corner, &
edge_lengths, &
                                  broadcast)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          buf, &
                          unlim_dim_level, &
                          corner, &
edge_lengths, &
                          broadcast)
end subroutine read_data_1d
!> @brief Read in data from a variable in a netcdf file.
subroutine read_data_2d(fileobj, &
                                  variable_name, &
                                  buf, &
                                  unlim_dim_level, &
                                  corner, &
edge_lengths, &
                                  broadcast)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          buf, &
                          unlim_dim_level, &
                          corner, &
edge_lengths, &
                          broadcast)
end subroutine read_data_2d
!> @brief Read in data from a variable in a netcdf file.
subroutine read_data_3d(fileobj, &
                                  variable_name, &
                                  buf, &
                                  unlim_dim_level, &
                                  corner, &
edge_lengths, &
                                  broadcast)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          buf, &
                          unlim_dim_level, &
                          corner, &
edge_lengths, &
                          broadcast)
end subroutine read_data_3d
!> @brief Read in data from a variable in a netcdf file.
subroutine read_data_4d(fileobj, &
                                  variable_name, &
                                  buf, &
                                  unlim_dim_level, &
                                  corner, &
edge_lengths, &
                                  broadcast)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          buf, &
                          unlim_dim_level, &
                          corner, &
edge_lengths, &
                          broadcast)
end subroutine read_data_4d
!> @brief Read in data from a variable in a netcdf file.
subroutine read_data_5d(fileobj, &
                                  variable_name, &
                                  buf, &
                                  unlim_dim_level, &
                                  corner, &
edge_lengths, &
                                  broadcast)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:,:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    call netcdf_read_data(fileobj, &
                          variable_name, &
                          buf, &
                          unlim_dim_level, &
                          corner, &
edge_lengths, &
                          broadcast)
end subroutine read_data_5d
