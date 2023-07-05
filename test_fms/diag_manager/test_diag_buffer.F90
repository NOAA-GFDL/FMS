program test_diag_buffer
#ifdef use_yaml

    use fms_diag_output_buffer_mod
    use platform_mod
    use diag_data_mod, only: i4, i8, r4, r8

    implicit none

    type(outputBuffer0d_type) :: buffobj0(10)
    type(outputBuffer1d_type) :: buffobj1
    type(outputBuffer2d_type) :: buffobj2
    type(outputBuffer3d_type) :: buffobj3
    type(outputBuffer4d_type) :: buffobj4
    type(outputBuffer5d_type) :: buffobj5
    class(*),allocatable       :: p_val, p_data1(:), p_data2(:,:)
    real(r8_kind)  :: r8_data
    real(r4_kind)  :: r4_data
    integer(i8_kind)  :: i8_data
    integer(i4_kind)  :: i4_data
    integer        :: buff_id
    class(*), pointer :: remap_buffer_out(:,:,:,:,:)
    integer :: i
    real(4) :: arr(9)
    real(4), allocatable :: arr1d(:)
    class(*), allocatable :: arr2d(:,:)
    integer(8), allocatable :: i8arr2d(:,:)
    real(8), allocatable :: r8val
    class(*), allocatable :: arr3d(:,:,:), arr4d(:,:,:,:), arr5d(:,:,:,:,:)
    integer(8), allocatable :: i8arr3d(:,:,:), i8arr4d(:,:,:,:), i8arr5d(:,:,:,:,:)
    logical :: test_5d = .true.
    character(len=4) :: fname = 'test'

    !! 0d
    ! allocate some buffers
    do i=1, 10
        call buffobj0(i)%allocate_buffer(r8_data, fname)
        call buffobj0(i)%initialize_buffer( real(i, kind=r8_kind) , fname)
    end do
    ! add some values
    call buffobj0(5)%add_to_buffer(real(-1, kind=r8_kind), fname)
    ! get the buffer data
    !allocate(real(8) :: p_val)
    !allocate(r8val)
    call buffobj0(5)%get_buffer(p_val, fname)
    select type(p_val)
      type is(real(r8_kind))
        print *, p_val
        r8val = p_val
    end select
    ! get the 5d remapped buffer data
    remap_buffer_out => buffobj0(5)%remap_buffer(fname, .false.)
    ! check output from object and remapped buffer
    print *, r8val
    call print_5d(remap_buffer_out)
    do i=1, 10
      call buffobj0(i)%flush_buffer()
    enddo

    !! 1d
    ! allocate a buffer to the given type and get it's id
    call buffobj1%allocate_buffer(r4_data, 10, fname)
    !! init to given value
    call buffobj1%initialize_buffer( real(0.1, kind=r4_kind), fname )
    !! add some values to the buffer
    arr = 4.0
    call buffobj1%add_to_buffer(arr, fname)
    !! get the buffer
    allocate(real(8) :: p_data1(10))
    allocate(arr1d(10))
    call buffobj1%get_buffer(p_data1, fname)
    select type(p_data1)
      type is(real(4))
        print *, p_data1
        arr1d = p_data1
    end select
    !! get the remapped buffer
    remap_buffer_out => buffobj1%remap_buffer(fname, .false.)
    !! check output
    print *, arr1d
    call print_5d(remap_buffer_out)
    call buffobj1%flush_buffer()
    print *, '********** 2d **********'

    !! 2d
    ! allocate a buffer to the given type and get it's id
    call buffobj2%allocate_buffer(i4_data, (/ 5, 10 /), fname )
    !!! init to given value
    call buffobj2%initialize_buffer( int(2, kind=i4_kind), fname )
    !! set some values in the buffer
    allocate(integer(4) :: arr2d(5,10))
    arr2d = 1
    call buffobj2%add_to_buffer(arr2d, fname)
    !!! get the buffer
    call buffobj2%get_buffer(arr2d, fname)
    !!! get the remapped buffer
    remap_buffer_out => buffobj2%remap_buffer(fname, .false.)
    !!! check output
    select type(arr2d)
      type is(integer(i4_kind))
        print *, arr2d
    end select
    call print_5d(remap_buffer_out)
    call buffobj2%flush_buffer()

    !! 3d
    ! allocate a buffer to the given type and get it's id
    call buffobj3%allocate_buffer(i8_data, (/ 2, 2, 2/), fname )
    !! init to given value
    call buffobj3%initialize_buffer( int(3, kind=i8_kind), fname )
    !! set some values in the buffer
    allocate(i8arr3d(2,2,2))
    i8arr3d = 6
    call buffobj3%add_to_buffer(i8arr3d, fname)
    !! get the buffer
    call buffobj3%get_buffer(arr3d, fname)
    !! get the remapped buffer
    remap_buffer_out => buffobj3%remap_buffer(fname, .false.)
    !! check output
    select type (arr3d)
      type is(integer(i8_kind))
        print *, arr3d
    end select
    call print_5d(remap_buffer_out)
    call buffobj3%flush_buffer()

    !! 4d
    ! allocate a buffer to the given type and get it's id
    call buffobj4%allocate_buffer(i8_data, (/ 2, 2, 2, 2/), fname)
    !! init to given value
    call buffobj4%initialize_buffer( int(4, kind=i8_kind), fname )
    !! set some values in the buffer
    allocate(i8arr4d(2,2,2,2))
    i8arr4d = 8
    call buffobj4%add_to_buffer(i8arr4d, fname)
    !! get the buffer
    call buffobj4%get_buffer(arr4d, fname)
    !! get the remapped buffer
    remap_buffer_out => buffobj4%remap_buffer(fname, .false.)
    !! check output
    select type (arr4d)
      type is(integer(i8_kind))
        print *, arr4d
    end select
    call print_5d(remap_buffer_out)
    call buffobj4%flush_buffer()

    !! 5d
    call buffobj5%allocate_buffer(i8_data, (/ 2, 2, 2, 2, 2/), fname )
    !! init to given value
    call buffobj5%initialize_buffer( int(5, kind=i8_kind), fname )
    !! get the remapped buffer
    remap_buffer_out => buffobj5%remap_buffer(fname, .false.)
    !! set some values in the buffer
    allocate(i8arr5d(2,2,2,2,2))
    i8arr5d = 10
    call buffobj5%add_to_buffer(i8arr5d, fname)
    !! get the buffer
    call buffobj5%get_buffer(arr5d, fname)
    !! check output
    select type (arr4d)
      type is(integer(i8_kind))
        print *, arr4d
    end select
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
