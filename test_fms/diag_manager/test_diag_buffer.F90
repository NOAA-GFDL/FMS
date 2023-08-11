program test_diag_buffer
#ifdef use_yaml

    use fms_diag_output_buffer_mod
    use platform_mod
    use fms_mod, only: string, fms_init, fms_end
    use mpp_mod, only: mpp_error, FATAL
    use diag_data_mod, only: i4, i8, r4, r8

    implicit none

    type(outputBuffer_type) :: buffobj(6)
    integer :: buff_sizes(5)
    class(*),allocatable       :: p_val(:,:,:,:,:)
    integer :: i, j
    real(r8_kind)  :: r8_data
    real(r4_kind)  :: r4_data
    integer(i8_kind)  :: i8_data
    integer(i4_kind)  :: i4_data
    character(len=4) :: fname = 'test'

    call fms_init
    buff_sizes = 1
    do i=0, 5
        if (i < 5) buff_sizes(i+1) = i+5
        call buffobj(i+1)%allocate_buffer(r8_data, i, buff_sizes, fname)
        call buffobj(i+1)%initialize_buffer( real(i, kind=r8_kind) , fname)
        call buffobj(i+1)%get_buffer(p_val, fname)
        select type(p_val)
        type is (real(kind=r8_kind))
          if (any(p_val .ne. real(i, kind=r8_kind))) &
            call mpp_error(FATAL, "The "//string(i)//"d buffer was not initialized to the correct value")
        do j = 1, 5
          if (size(p_val, j) .ne. buff_sizes(j)) &
            call mpp_error(FATAL, "The "//string(i)//"d buffer was not allocated to the size")
        enddo
        class default
          call mpp_error(FATAL, "The "//string(i)//"d buffer was not allocated to the correct type")
        end select
        deallocate(p_val)
        call buffobj(i+1)%flush_buffer()
    end do

    buff_sizes = 1
    do i=0, 5
        if (i < 5) buff_sizes(i+1) = i+5
        call buffobj(i+1)%allocate_buffer(r4_data, i, buff_sizes, fname)
        call buffobj(i+1)%initialize_buffer( real(i, kind=r4_kind) , fname)
        call buffobj(i+1)%get_buffer(p_val, fname)
        select type(p_val)
        type is (real(kind=r4_kind))
          if (any(p_val .ne. real(i, kind=r4_kind))) &
            call mpp_error(FATAL, "The "//string(i)//"d buffer was not initialized to the correct value")
        do j = 1, 5
          if (size(p_val, j) .ne. buff_sizes(j)) &
            call mpp_error(FATAL, "The "//string(i)//"d buffer was not allocated to the size")
        enddo
        class default
          call mpp_error(FATAL, "The "//string(i)//"d buffer was not allocated to the correct type")
        end select
        deallocate(p_val)
        call buffobj(i+1)%flush_buffer()
    end do

    buff_sizes = 1
    do i=0, 5
        if (i < 5) buff_sizes(i+1) = i+5
        call buffobj(i+1)%allocate_buffer(i8_data, i, buff_sizes, fname)
        call buffobj(i+1)%initialize_buffer( int(i, kind=i8_kind) , fname)
        call buffobj(i+1)%get_buffer(p_val, fname)
        select type(p_val)
        type is (integer(kind=i8_kind))
          if (any(p_val .ne. int(i, kind=i8_kind))) &
            call mpp_error(FATAL, "The "//string(i)//"d buffer was not initialized to the correct value")
        do j = 1, 5
          if (size(p_val, j) .ne. buff_sizes(j)) &
            call mpp_error(FATAL, "The "//string(i)//"d buffer was not allocated to the size")
        enddo
        class default
          call mpp_error(FATAL, "The "//string(i)//"d buffer was not allocated to the correct type")
        end select
        deallocate(p_val)
        call buffobj(i+1)%flush_buffer()
    end do

    buff_sizes = 1
    do i=0, 5
        if (i < 5) buff_sizes(i+1) = i+5
        call buffobj(i+1)%allocate_buffer(i4_data, i, buff_sizes, fname)
        call buffobj(i+1)%initialize_buffer( int(i, kind=i4_kind) , fname)
        call buffobj(i+1)%get_buffer(p_val, fname)
        select type(p_val)
        type is (integer(kind=i4_kind))
          if (any(p_val .ne. int(i, kind=i4_kind))) &
            call mpp_error(FATAL, "The "//string(i)//"d buffer was not initialized to the correct value")
        do j = 1, 5
          if (size(p_val, j) .ne. buff_sizes(j)) &
            call mpp_error(FATAL, "The "//string(i)//"d buffer was not allocated to the size")
        enddo
        class default
          call mpp_error(FATAL, "The "//string(i)//"d buffer was not allocated to the correct type")
        end select
        deallocate(p_val)
        call buffobj(i+1)%flush_buffer()
    end do

    call fms_end()
#endif
end program
