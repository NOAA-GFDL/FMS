program test_diag_bufer

    use fms_diag_buffer_object_mod
    use platform_mod

    implicit none

    type(buffer0d) :: buffobj0(10)
    type(buffer1d) :: buffobj1, tmp1
    type(buffer2d) :: buffobj2, tmp2
    type(buffer3d) :: buffobj3, tmp3
    type(buffer4d) :: buffobj4, tmp4
    type(buffer5d) :: buffobj5, tmp5
    class(fmsDiagBuffer_class), allocatable :: tmp
    real(r8_kind)  :: r8_data
    real(r4_kind)  :: r4_data
    integer(i8_kind)  :: i8_data
    integer(i4_kind)  :: i4_data
    integer        :: buff_id
    class(*), allocatable       :: buffer_out
    class(*), allocatable       :: buffer_out1(:)
    class(*), allocatable       :: buffer_out2(:,:)
    class(*), pointer :: remap_buffer_out(:,:,:,:,:)
    integer :: i

    !procedure :: get_remapped_buffer_pointer
    !procedure :: get_area
    !procedure :: get_volume
    !procedure :: get_missing_value
    !procedure :: get_data_RANGE
    !cedure :: allocate_buffer => allocate_buffer_0d
    !procedure :: get_buffer => get_buffer_0d
    !procedure :: initialize_buffer => initialize_buffer_0d
         
    !! 0d
    do i=1, 10 
        ! allocate a buffer to the given type and get it's id
        print *, i, '1'
        buff_id = buffobj0(i)%allocate_buffer(r8_data)
        call buffobj0(i)%initialize_buffer( real(i, kind=r8_kind) )
        print *, i, '2'
        ! check object data 
        buffer_out = buffobj0(i)%get_buffer()  
        print *, i, '3'
        select type (buffer_out)
            type is(real(8))
                print *, buffer_out
        end select
        ! check data from list 
        select type (obj => buffer_list_0d(i)%buffer_obj)
            type is(buffer0d)
                select type(bdata => obj%get_buffer())
                    type is(real(8))
                        print *, bdata
                end select
        end select
    end do
    ! get the buffer
    buffer_out = buffobj0(5)%get_buffer()  
    ! get the remapped buffer
    remap_buffer_out => buffobj0(5)%get_remapped_buffer_pointer()
    ! check output
    select type (buffer_out)
      type is(real(8))
        print *, buffer_out
    end select
    call print_5d(remap_buffer_out)
    ! get buffer from list with id
    allocate(buffer0d :: tmp)
    tmp = get_buffer_object(1, 0)
    select type (tmp)
        type is (buffer0d)
            buffer_out = tmp%get_buffer()
    end select
    select type (buffer_out)
      type is(real(8))
        print *, buffer_out
    end select

    !! 1d
    ! allocate a buffer to the given type and get it's id
    buff_id = buffobj1%allocate_buffer(r4_data, 10)
    ! init to given value 
    call buffobj1%initialize_buffer( real(0.2, kind=r4_kind) )
    ! get the buffer
    buffer_out1= buffobj1%get_buffer()  
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
    buff_id = buffobj2%allocate_buffer(i4_data, (/ 5, 10 /) )
    !! init to given value 
    call buffobj2%initialize_buffer( int(3, kind=i4_kind) )
    !! get the buffer
    buffer_out2= buffobj2%get_buffer()  
    !! get the remapped buffer
    remap_buffer_out => buffobj2%get_remapped_buffer_pointer()
    !! check output
    select type (buffer_out2)
        type is(integer(4))
            print *, buffer_out2
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
