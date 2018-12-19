subroutine write_data_0d(fileobj, &
                                   variable_name, &
                                   variable_data, &
                                   unlim_dim_level, &
                                   corner &

                                  )
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*), intent(in) :: variable_data !< Data that will be
                                                               !! written.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(0),intent(in),optional :: corner !< Array of starting
                                                              !! indices describing
                                                              !! where the data
                                                              !! will be written to.

    call netcdf_write_data(fileobj, &
                           variable_name, &
                           variable_data, &
                           unlim_dim_level, &
                           corner &

                          )
end subroutine write_data_0d
subroutine write_data_1d(fileobj, &
                                   variable_name, &
                                   variable_data, &
                                   unlim_dim_level, &
                                   corner &
,edge_lengths &
                                  )
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:), intent(in) :: variable_data !< Data that will be
                                                               !! written.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(1),intent(in),optional :: corner !< Array of starting
                                                              !! indices describing
                                                              !! where the data
                                                              !! will be written to.
integer,dimension(1),intent(in),optional :: edge_lengths !< The number of
                                                                    !! elements that
                                                                    !! will be written
                                                                    !! in each dimension.
    call netcdf_write_data(fileobj, &
                           variable_name, &
                           variable_data, &
                           unlim_dim_level, &
                           corner &
,edge_lengths &
                          )
end subroutine write_data_1d
subroutine write_data_2d(fileobj, &
                                   variable_name, &
                                   variable_data, &
                                   unlim_dim_level, &
                                   corner &
,edge_lengths &
                                  )
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:), intent(in) :: variable_data !< Data that will be
                                                               !! written.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(2),intent(in),optional :: corner !< Array of starting
                                                              !! indices describing
                                                              !! where the data
                                                              !! will be written to.
integer,dimension(2),intent(in),optional :: edge_lengths !< The number of
                                                                    !! elements that
                                                                    !! will be written
                                                                    !! in each dimension.
    call netcdf_write_data(fileobj, &
                           variable_name, &
                           variable_data, &
                           unlim_dim_level, &
                           corner &
,edge_lengths &
                          )
end subroutine write_data_2d
subroutine write_data_3d(fileobj, &
                                   variable_name, &
                                   variable_data, &
                                   unlim_dim_level, &
                                   corner &
,edge_lengths &
                                  )
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:), intent(in) :: variable_data !< Data that will be
                                                               !! written.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(3),intent(in),optional :: corner !< Array of starting
                                                              !! indices describing
                                                              !! where the data
                                                              !! will be written to.
integer,dimension(3),intent(in),optional :: edge_lengths !< The number of
                                                                    !! elements that
                                                                    !! will be written
                                                                    !! in each dimension.
    call netcdf_write_data(fileobj, &
                           variable_name, &
                           variable_data, &
                           unlim_dim_level, &
                           corner &
,edge_lengths &
                          )
end subroutine write_data_3d
subroutine write_data_4d(fileobj, &
                                   variable_name, &
                                   variable_data, &
                                   unlim_dim_level, &
                                   corner &
,edge_lengths &
                                  )
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:), intent(in) :: variable_data !< Data that will be
                                                               !! written.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(4),intent(in),optional :: corner !< Array of starting
                                                              !! indices describing
                                                              !! where the data
                                                              !! will be written to.
integer,dimension(4),intent(in),optional :: edge_lengths !< The number of
                                                                    !! elements that
                                                                    !! will be written
                                                                    !! in each dimension.
    call netcdf_write_data(fileobj, &
                           variable_name, &
                           variable_data, &
                           unlim_dim_level, &
                           corner &
,edge_lengths &
                          )
end subroutine write_data_4d
subroutine write_data_5d(fileobj, &
                                   variable_name, &
                                   variable_data, &
                                   unlim_dim_level, &
                                   corner &
,edge_lengths &
                                  )
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:,:), intent(in) :: variable_data !< Data that will be
                                                               !! written.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(5),intent(in),optional :: corner !< Array of starting
                                                              !! indices describing
                                                              !! where the data
                                                              !! will be written to.
integer,dimension(5),intent(in),optional :: edge_lengths !< The number of
                                                                    !! elements that
                                                                    !! will be written
                                                                    !! in each dimension.
    call netcdf_write_data(fileobj, &
                           variable_name, &
                           variable_data, &
                           unlim_dim_level, &
                           corner &
,edge_lengths &
                          )
end subroutine write_data_5d
