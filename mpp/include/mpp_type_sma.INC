subroutine MPP_TYPE_CREATE_(field, array_of_subsizes, array_of_starts, dtype)
    MPP_TYPE_, intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    call mpp_error(FATAL, 'MPP_TYPE_CREATE_: Unsupported in SHMEM.')
end subroutine MPP_TYPE_CREATE_
