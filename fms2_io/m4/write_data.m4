include(`macros.m4')dnl
`subroutine write_data_'NUM_DIMS`d(fileobj, &'
                                   variable_name, &
                                   variable_data, &
                                   unlim_dim_level, &
                                   corner &
ifelse(NUM_DIMS,0,,
                                 `,edge_lengths &')
                                  )
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dim_declare(NUM_DIMS) intent(in) :: variable_data !< Data that will be
                                                               !! written.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(NUM_DIMS),intent(in),optional :: corner !< Array of starting
                                                              !! indices describing
                                                              !! where the data
                                                              !! will be written to.
ifelse(NUM_DIMS,0,,
   `integer,dimension(NUM_DIMS),intent(in),optional :: edge_lengths !< The number of
                                                                    !! elements that
                                                                    !! will be written
                                                                    !! in each dimension.')
    call netcdf_write_data(fileobj, &
                           variable_name, &
                           variable_data, &
                           unlim_dim_level, &
                           corner &
ifelse(NUM_DIMS,0,,
                         `,edge_lengths &')
                          )
`end subroutine write_data_'NUM_DIMS`d'
