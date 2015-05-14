subroutine MPP_ALLTOALL_(sbuf, scount, rbuf, rcount, pelist)

    MPP_TYPE_, dimension(:), intent(in) :: sbuf
    MPP_TYPE_, dimension(:), intent(inout) :: rbuf
    integer,   intent(in) :: scount, rcount

    integer, intent(in), optional :: pelist(0:)

    call mpp_error(FATAL, 'MPP_ALLTOALL: No SHMEM implementation.')

end subroutine MPP_ALLTOALL_


subroutine MPP_ALLTOALLV_(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)

    MPP_TYPE_, intent(in) :: sbuf(:)
    MPP_TYPE_, intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)

    call mpp_error(FATAL, 'MPP_ALLTOALLV: No SHMEM implementation.')

end subroutine MPP_ALLTOALLV_
