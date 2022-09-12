program test_diag_buffer
#ifdef use_yaml

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
    class(*), pointer :: remap_buffer_out(:,:,:,:,:)
    integer :: i
    real(4) :: arr(9)
    real(4), allocatable :: arr1d(:)
    integer(4), allocatable :: arr2d(:,:)
    real(8), allocatable :: r8val
    integer(8), allocatable :: i8arr3d(:,:,:), i8arr4d(:,:,:,:), i8arr5d(:,:,:,:,:)
    logical :: test_5d = .true.

    !! 0d
    ! allocate some buffers
    do i=1, 10
        call buffobj0(i)%allocate_buffer(r8_data)
        call buffobj0(i)%initialize_buffer( real(i, kind=r8_kind) )
    end do
    ! add some values
    call buffobj0(5)%add_to_buffer(real(-1, kind=r8_kind))
    ! get the buffer data
    call buffobj0(5)%get_buffer(r8val)
    ! get the 5d remapped buffer data
    remap_buffer_out => buffobj0(5)%remap_buffer()
    ! check output from object and remapped buffer
    print *, r8val
    call print_5d(remap_buffer_out)
    do i=1, 10
      call buffobj0(i)%flush_buffer()
    enddo

    !! 1d
    ! allocate a buffer to the given type and get it's id
    call buffobj1%allocate_buffer(r4_data, 10)
    !! init to given value
    call buffobj1%initialize_buffer( real(0.1, kind=r4_kind) )
    !! add some values to the buffer
    arr = 4.0
    call buffobj1%add_to_buffer(arr)
    !! get the buffer
    allocate(arr1d(10))
    call buffobj1%get_buffer(arr1d)
    !! get the remapped buffer
    remap_buffer_out => buffobj1%remap_buffer()
    !! check output
    print *, arr1d
    call print_5d(remap_buffer_out)
    call buffobj1%flush_buffer()


    !! 2d
    ! allocate a buffer to the given type and get it's id
    call buffobj2%allocate_buffer(i4_data, (/ 5, 10 /) )
    !!! init to given value
    call buffobj2%initialize_buffer( int(2, kind=i4_kind) )
    !! set some values in the buffer
    allocate(arr2d(5,10))
    arr2d = 1
    call buffobj2%add_to_buffer(arr2d)
    !!! get the buffer
    call buffobj2%get_buffer(arr2d)
    !!! get the remapped buffer
    remap_buffer_out => buffobj2%remap_buffer()
    !!! check output
    print *, arr2d
    call print_5d(remap_buffer_out)
    call buffobj2%flush_buffer()

    !! 3d
    ! allocate a buffer to the given type and get it's id
    call buffobj3%allocate_buffer(i8_data, (/ 2, 2, 2/) )
    !! init to given value
    call buffobj3%initialize_buffer( int(3, kind=i8_kind) )
    !! set some values in the buffer
    allocate(i8arr3d(2,2,2))
    i8arr3d = 6
    call buffobj3%add_to_buffer(i8arr3d)
    !! get the buffer
    call buffobj3%get_buffer(i8arr3d)
    !! get the remapped buffer
    remap_buffer_out => buffobj3%remap_buffer()
    !! check output
    print *, i8arr3d
    call print_5d(remap_buffer_out)
    call buffobj3%flush_buffer()

    !! 4d
    ! allocate a buffer to the given type and get it's id
    call buffobj4%allocate_buffer(i8_data, (/ 2, 2, 2, 2/) )
    !! init to given value
    call buffobj4%initialize_buffer( int(4, kind=i8_kind) )
    !! set some values in the buffer
    allocate(i8arr4d(2,2,2,2))
    i8arr4d = 8
    call buffobj4%add_to_buffer(i8arr4d)
    !! get the buffer
    call buffobj4%get_buffer(i8arr4d)
    !! get the remapped buffer
    remap_buffer_out => buffobj4%remap_buffer()
    !! check output
    print *, i8arr4d
    call print_5d(remap_buffer_out)
    call buffobj4%flush_buffer()

    !! 5d
    call buffobj5%allocate_buffer(i8_data, (/ 2, 2, 2, 2, 2/) )
    !! init to given value
    call buffobj5%initialize_buffer( int(5, kind=i8_kind) )
    !! get the remapped buffer
    remap_buffer_out => buffobj5%remap_buffer()
    !! set some values in the buffer
    allocate(i8arr5d(2,2,2,2,2))
    i8arr5d = 10
    call buffobj5%add_to_buffer(i8arr5d)
    !! get the buffer
    call buffobj5%get_buffer(i8arr5d)
    !! check output
    print *, i8arr5d
    call print_5d(remap_buffer_out)
    call buffobj5%flush_buffer()

  contains

  ! just prints polymorphic data types
  subroutine print_5d(val)
    class(*), intent(in) :: val(:,:,:,:,:)

    select type (val)
      type is (real(r4_kind))
        print *, "5d:", val
      type is (real(r8_kind))
        print *, "5d:", val
      type is (integer(i4_kind))
        print *, "5d:",val
      type is (integer(i8_kind))
        print *, "5d:",val
    end select
  end subroutine



#endif
end program
