subroutine MPP_TYPE_CREATE_(field, array_of_subsizes, array_of_starts, dtype)
    MPP_TYPE_, intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call SYSTEM_CLOCK(start_tick)

    call mpp_error(NOTE, 'MPP_TYPE_CREATE: &
                         &This function is not used in serial mode.')

    ! For consistency with the MPI version, we return a valid mpp_type
    dtype = mpp_byte

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, MPP_TYPE_BYTELEN_)

end subroutine MPP_TYPE_CREATE_
