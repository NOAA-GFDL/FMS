program test_metadata_transfer
  use fms_mod, only: fms_init, fms_end, string
  use mpp_mod
  use metadata_transfer_mod
  use platform_mod
  use mpi
  use, intrinsic :: iso_c_binding

  implicit none

  type(metadata_r8_type) :: file_metadata(3)

  logical :: debug = .true.
  integer :: i

  call fms_init()

  ! all PEs need to initialize the metadata object with a datatype
  do i=1, 3
    call file_metadata(i)%fms_metadata_transfer_init(real8_type)
  enddo

  ! set metadata only on root PE
  if (mpp_pe() .eq. mpp_root_pe()) then
    call file_metadata(1)%set_attribute_name("_FillValue"//c_null_char)
    call file_metadata(1)%set_attribute_type(real8_type)
    call file_metadata(1)%set_attribute_value([666.0_r8_kind])

    call file_metadata(2)%set_attribute_name("missing_value"//c_null_char)
    call file_metadata(2)%set_attribute_type(real8_type)
    call file_metadata(2)%set_attribute_value([-100.0_r8_kind, 100.0_r8_kind])

    call file_metadata(3)%set_attribute_name("a_third_name"//c_null_char)
    call file_metadata(3)%set_attribute_type(real8_type)
    call file_metadata(3)%set_attribute_value([-200.0_r8_kind, -50.0_r8_kind, 0.0_r8_kind, 50.0_r8_kind, 200.0_r8_kind])
  endif
  ! Broadcast the metadata to all PEs
  call fms_metadata_broadcast_all(file_metadata)

  if(debug) call dump_metadata(file_metadata)

  call check_metadata(file_metadata)

  call fms_end()

  contains

  subroutine dump_metadata(this)
    type(metadata_r8_type), intent(inout) :: this(:)
    real(r8_kind), allocatable :: arr(:)

    integer :: i
    do i = 1, size(this)
      arr = this(i)%get_attribute_value()
      print *, "pe: ", mpp_pe(), "i: ", i, " attribute_name is ", trim(adjustl(this(i)%get_attribute_name()))
      print *, "pe: ", mpp_pe(), "i: ", i, " attribute_type is ", string(this(i)%get_attribute_type())
      print *, "pe: ", mpp_pe(), "i: ", i, " attribute_value is ", arr 
    enddo
  end subroutine

  subroutine check_metadata(this)
    type(metadata_r8_type), intent(inout) :: this(:)
    real(r8_kind), allocatable :: arr(:)
    logical :: is_correct
    character(len=32) :: attr_names(3)
    real(r8_kind) :: attr_vals(8)
    integer :: i, j, last_j =1

    attr_names = (/"_FillValue"//c_null_char, "missing_value"//c_null_char, "a_third_name"//c_null_char/)
    attr_vals = (/ 666.0_r8_kind, -100.0_r8_kind, 100.0_r8_kind, -200.0_r8_kind, &
                   -50.0_r8_kind, 0.0_r8_kind, 50.0_r8_kind, 200.0_r8_kind /) 

    do i = 1, size(this)
      arr = this(i)%get_attribute_value()
      if (trim(this(i)%get_attribute_name()) .ne. attr_names(i)) then
        call mpp_error(FATAL, "incorrect metadata name")
      endif 
      if( this(i)%get_attribute_type() .ne. real8_type) then
        call mpp_error(FATAL, "incorrect metadata type")
      endif

      do j=1, size(arr) 
        if (arr(j) .ne. attr_vals(last_j)) then
          print *, "got ", arr(j), " expected ", attr_vals(last_j) 
          call mpp_error(FATAL, "incorrect metadata value")
        endif
        last_j = last_j + 1
      enddo

    enddo
  end subroutine

end program