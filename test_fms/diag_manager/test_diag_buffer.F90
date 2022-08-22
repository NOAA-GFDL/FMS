program test_diag_bufer

    use fms_diag_buffer_mod
    use platform_mod

    implicit none

    type(buffer0d) :: buffobj0(10), tmp0
    type(buffer1d) :: buffobj1, tmp1
    type(buffer2d) :: buffobj2, tmp2
    type(buffer3d) :: buffobj3, tmp3
    type(buffer4d) :: buffobj4, tmp4
    type(buffer5d) :: buffobj5, tmp5
    real(r8_kind)  :: r8_data
    real(r4_kind)  :: r4_data
    integer(i8_kind)  :: i8_data
    integer(i4_kind)  :: i4_data
    integer        :: buff_id
    class(*), allocatable       :: buffer_out, buffer_out1(:), buffer_out2(:,:), buffer_out3(:,:,:)
    class(*), allocatable       :: buffer_out4(:,:,:,:), buffer_out5(:,:,:,:,:)
    class(*), pointer :: remap_buffer_out(:,:,:,:,:)
    integer :: i
    real(4) :: arr(9)


    !! 0d
    ! allocate some buffers
    do i=1, 10
        call buffobj0(i)%allocate_buffer(r8_data)
        call buffobj0(i)%initialize_buffer( real(i, kind=r8_kind) )
    end do
    ! add some values
    call buffobj0(5)%add_to_buffer(real(-1, kind=r8_kind))
    ! get the buffer data
    buffer_out = buffobj0(5)%get_buffer_data()
    ! get the 5d remapped buffer data
    ! TODO scalar buffer remapping
    !remap_buffer_out => buffobj0(5)%get_remapped_buffer_pointer()
    ! check output from object and remapped buffer
    select type (buffer_out)
      type is(real(8))
        print *, buffer_out
    end select
    !call print_5d(remap_buffer_out)

    !! 1d
    ! allocate a buffer to the given type and get it's id
    call buffobj1%allocate_buffer(r4_data, 10)
    ! init to given value
    call buffobj1%initialize_buffer( real(0.1, kind=r4_kind) )
    ! add some values to the buffer
    arr = 4.0
    call buffobj1%add_to_buffer(arr)
    ! get the buffer
    buffer_out1= buffobj1%get_buffer_data()
    ! get the remapped buffer
    remap_buffer_out => buffobj1%get_remapped_buffer_pointer()
    ! check output
    select type (buffer_out1)
      type is(real(4))
        print *, buffer_out1
    end select
    call print_5d(remap_buffer_out)


    !! 2d
    ! allocate a buffer to the given type and get it's id
    call buffobj2%allocate_buffer(i4_data, (/ 5, 10 /) )
    !! init to given value
    call buffobj2%initialize_buffer( int(2, kind=i4_kind) )
    !! get the buffer
    buffer_out2= buffobj2%get_buffer_data()
    !! get the remapped buffer
    remap_buffer_out => buffobj2%get_remapped_buffer_pointer()
    !! check output
    select type (buffer_out2)
        type is(integer(4))
            print *, buffer_out2
    end select
    call print_5d(remap_buffer_out)


    !! 3d
    ! allocate a buffer to the given type and get it's id
    call buffobj3%allocate_buffer(i8_data, (/ 5, 10, 5/) )
    !! init to given value
    call buffobj3%initialize_buffer( int(3, kind=i8_kind) )
    !! get the buffer
    buffer_out3= buffobj3%get_buffer_data()
    !! get the remapped buffer
    remap_buffer_out => buffobj3%get_remapped_buffer_pointer()
    !! check output
    select type (buffer_out3)
        type is(integer(8))
            print *, buffer_out3
    end select
    call print_5d(remap_buffer_out)

    !! 4d
    ! allocate a buffer to the given type and get it's id
    call buffobj4%allocate_buffer(i8_data, (/ 5, 5, 5, 5/) )
    !! init to given value
    call buffobj4%initialize_buffer( int(4, kind=i8_kind) )
    !! get the buffer
    buffer_out4= buffobj4%get_buffer_data()
    !! get the remapped buffer
    remap_buffer_out => buffobj4%get_remapped_buffer_pointer()
    !! check output
    select type (buffer_out4)
        type is(integer(8))
            print *, buffer_out4
    end select
    call print_5d(remap_buffer_out)


    !! 5d
    ! allocate a buffer to the given type and get it's id
    call buffobj5%allocate_buffer(i8_data, (/ 5, 5, 5, 5, 5/) )
    !! init to given value
    call buffobj5%initialize_buffer( int(5, kind=i8_kind) )
    !! get the buffer
    buffer_out5= buffobj5%get_buffer_data()
    !! get the remapped buffer
    remap_buffer_out => buffobj5%get_remapped_buffer_pointer()
    !! check output
    select type (buffer_out5)
        type is(integer(8))
            print *, buffer_out5
    end select
    call print_5d(remap_buffer_out)

    contains
    ! just prints polymorphic data types
    subroutine print_5d(val)
        class(*), intent(in) :: val(:,:,:,:,:)

        select type (val)
        type is (real(r4_kind))
            print *, val
        type is (real(r8_kind))
            print *, val
        type is (integer(i4_kind))
            print *, val
        end select
    end subroutine
end program
