subroutine MPP_TYPE_CREATE_(field, array_of_subsizes, array_of_starts, &
                            dtype_out)
    MPP_TYPE_, intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype_out

    type(mpp_type), pointer :: dtype
    integer :: newtype      ! MPI datatype ID

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE_: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call SYSTEM_CLOCK(start_tick)

    if (verbose) &
        call mpp_error(NOTE, 'MPP_TYPE_CREATE_: &
                             &using MPI_Type_create_subarray...')

    dtype => datatypes%head
    ! TODO: Check mpp_byte

    ! Check if mpp_type already exists
    do while (.not. associated(dtype))
        dtype => dtype%next

        if (dtype%ndims /= rank(field)) cycle
        if (any(dtype%sizes /= shape(field))) cycle
        if (any(dtype%subsizes /= array_of_subsizes)) cycle
        if (any(dtype%starts /= array_of_starts)) cycle
        if (dtype%etype /= MPI_TYPE_) cycle

        ! If all parameters match, then the datatype exists and return dtype
        dtype%counter = dtype%counter + 1
        dtype_out = dtype
        return
    end do

    ! The type does not exist; create a new internal type
    call MPI_Type_create_subarray( &
        rank(field), &
        shape(field), &
        array_of_subsizes, &
        array_of_starts, &
        MPI_ORDER_FORTRAN, &
        MPI_TYPE_, &
        newtype, &
        error &
    )

    ! Register on the MPI runtime
    call MPI_Type_commit(newtype, error)

    ! Create new entry
    allocate(dtype)
    allocate(dtype%sizes(rank(field)))
    allocate(dtype%subsizes(rank(field)))
    allocate(dtype%starts(rank(field)))

    ! Populate values
    dtype%counter = 1
    dtype%ndims = rank(field)
    dtype%sizes = shape(field)
    dtype%subsizes = array_of_subsizes
    dtype%starts = array_of_starts
    dtype%etype = MPI_TYPE_
    dtype%id = newtype

    ! Add dtype to the list
    dtype%prev => datatypes%tail
    dtype%prev%next => dtype
    datatypes%tail => dtype
    datatypes%length = datatypes%length + 1

    ! Copy dtype to output
    dtype_out = dtype

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, MPP_TYPE_BYTELEN_)

end subroutine MPP_TYPE_CREATE_
