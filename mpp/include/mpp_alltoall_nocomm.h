subroutine MPP_ALLTOALL_(sbuf, scount, rbuf, rcount, pelist)

    MPP_TYPE_, dimension(:), intent(in) :: sbuf
    MPP_TYPE_, dimension(:), intent(inout) :: rbuf
    integer,   intent(in) :: scount, rcount

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call SYSTEM_CLOCK(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, MPP_TYPE_BYTELEN_)

end subroutine MPP_ALLTOALL_


subroutine MPP_ALLTOALLV_(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)

    MPP_TYPE_, intent(in) :: sbuf(:)
    MPP_TYPE_, intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call SYSTEM_CLOCK(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, MPP_TYPE_BYTELEN_)

end subroutine MPP_ALLTOALLV_
