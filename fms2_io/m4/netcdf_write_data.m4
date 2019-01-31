include(`macros.m4')dnl
!> @brief Write data to a variable in a netcdf file.
`subroutine netcdf_write_data_'NUM_DIMS`d(fileobj, &'
                                          variable_name, &
                                          variable_data, &
                                          unlim_dim_level, &
                                          corner &
ifelse(NUM_DIMS,0,,
                                        `,edge_lengths &')
                                         )
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
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
    integer :: err
    integer :: varid
    integer :: unlim_dim_index
    integer,dimension(incr(NUM_DIMS)) :: c
ifelse(NUM_DIMS,0,,
   `integer,dimension(incr(NUM_DIMS)) :: e')
    if (fileobj%is_root) then
        c(:) = 1
        if (present(corner)) then
            c(1:NUM_DIMS) = corner(:)
        endif
ifelse(NUM_DIMS,0,,
       `e(:) = 1
        if (present(edge_lengths)) then
            e(1:NUM_DIMS) = edge_lengths(:)
        else
            e(1:NUM_DIMS) = shape(variable_data)
        endif')
        if (present(unlim_dim_level)) then
            unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, &
                                                                     variable_name, &
                                                                     broadcast=.false.)
            if (unlim_dim_index .ne. incr(NUM_DIMS)) then
                call error("unlimited dimension must be the slowest varying" &
                           //" dimension.")
            endif
            c(unlim_dim_index) = unlim_dim_level
        endif
        call register_variable_attribute(fileobj, &
                                         variable_name, &
                                         "checksum", &
                                         get_checksum(variable_data))
        call set_netcdf_mode(fileobj%ncid, &
                             data_mode)
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        select type(variable_data)
            type is (integer(kind=int32))
                err = nf90_put_var(fileobj%ncid, &
                                   varid, &
                                   variable_data, &
                                   start=c &
ifelse(NUM_DIMS,0,,
                                 `,count=e &')
                                  )
            type is (integer(kind=int64))
                err = nf90_put_var(fileobj%ncid, &
                                   varid, &
                                   variable_data, &
                                   start=c &
ifelse(NUM_DIMS,0,,
                                 `,count=e &')
                                  )
            type is (real(kind=real32))
                err = nf90_put_var(fileobj%ncid, &
                                   varid, &
                                   variable_data, &
                                   start=c &
ifelse(NUM_DIMS,0,,
                                 `,count=e &')
                                  )
            type is (real(kind=real64))
                err = nf90_put_var(fileobj%ncid, &
                                   varid, &
                                   variable_data, &
                                   start=c &
ifelse(NUM_DIMS,0,,
                                 `,count=e &')
                                  )
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
`end subroutine netcdf_write_data_'NUM_DIMS`d'
