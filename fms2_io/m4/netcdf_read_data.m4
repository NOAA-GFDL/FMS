include(`macros.m4')dnl
!> @brief Read in data from a variable in a netcdf file.
`subroutine netcdf_read_data_'NUM_DIMS`d(fileobj, &'
                                         variable_name, &
                                         buf, &
                                         unlim_dim_level, &
                                         corner, &
ifelse(NUM_DIMS,0,,
                                        `edge_lengths, &')
                                         broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dim_declare(NUM_DIMS) intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
ifelse(NUM_DIMS,0,,
   `integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.')
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: bcast
    integer :: err
    integer :: varid
    integer :: unlim_dim_index
    integer,dimension(incr(NUM_DIMS)) :: c
ifelse(NUM_DIMS,0,,
   `integer,dimension(incr(NUM_DIMS)) :: e')
    if (present(broadcast)) then
        bcast = broadcast
    else
        bcast = .true.
    endif
    c(:) = 1
    if (present(corner)) then
        c(1:NUM_DIMS) = corner(:)
    endif
ifelse(NUM_DIMS,0,,
   `e(:) = 1
    if (present(edge_lengths)) then
        e(1:NUM_DIMS) = edge_lengths(:)
    else
        e(1:NUM_DIMS) = shape(buf)
    endif')
    if (present(unlim_dim_level)) then
        unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, &
                                                                 variable_name, &
                                                                 broadcast=bcast)
        if (unlim_dim_index .ne. incr(NUM_DIMS)) then
            call error("unlimited dimension must be the slowest varying" &
                       //" dimension.")
        endif
        c(unlim_dim_index) = unlim_dim_level
    endif
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        select type(buf)
            type is (integer(kind=int32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
ifelse(NUM_DIMS,0,,
                                 `,count=e &')
                                  )
            type is (integer(kind=int64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
ifelse(NUM_DIMS,0,,
                                 `,count=e &')
                                  )
            type is (real(kind=real32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
ifelse(NUM_DIMS,0,,
                                 `,count=e &')
                                  )
            type is (real(kind=real64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
ifelse(NUM_DIMS,0,,
                                 `,count=e &')
                                  )
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
    if (bcast) then
        select type(buf)
            type is (integer(kind=int32))
                call mpp_broadcast(buf, &
ifelse(NUM_DIMS,0,,
                                  `size(buf), &')
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (integer(kind=int64))
                call mpp_broadcast(buf, &
ifelse(NUM_DIMS,0,,
                                  `size(buf), &')
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real32))
                call mpp_broadcast(buf, &
ifelse(NUM_DIMS,0,,
                                  `size(buf), &')
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real64))
                call mpp_broadcast(buf, &
ifelse(NUM_DIMS,0,,
                                  `size(buf), &')
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            class default
                call error("unsupported type.")
        end select
    endif
`end subroutine netcdf_read_data_'NUM_DIMS`d'
