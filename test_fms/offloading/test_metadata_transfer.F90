program test_metadata_transfer
  use fms_mod, only: fms_init, fms_end, string
  use mpp_mod
  use metadata_transfer_mod
  use platform_mod
  use mpi
  use, intrinsic :: iso_c_binding

  implicit none

  type(metadata_type) :: file_metadata(3)

  logical :: debug = .false.

  call fms_init()

  call fms_metadata_transfer_init()

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
    type(metadata_type), intent(in) :: this(:)
    real(r8_kind), allocatable :: arr(:)

    integer :: i
    do i = 1, size(this)
      arr = this(i)%get_attribute_value()
      print *, "pe: ", mpp_pe(), " attribute_name is ", trim(adjustl(this(i)%get_attribute_name()))
      print *, "pe: ", mpp_pe(), " attribute_type is ", string(this(i)%get_attribute_type())
      print *, "pe: ", mpp_pe(), " attribute_value is ", arr 
    enddo
  end subroutine

  subroutine check_metadata(this)
    type(metadata_type), intent(in) :: this(:)
    real(r8_kind), allocatable :: arr(:)
    logical :: is_correct
    character(len=32) :: attr_names(3)
    real(r8_kind) :: attr_vals(8)
    integer :: i, i_expected = 1
    integer :: val_inds(3) = (/ 0, 1, 4/) 

    attr_names = (/"_FillValue"//c_null_char, "missing_value"//c_null_char, "a_third_name"//c_null_char/)
    attr_vals = (/ 666.0_r8_kind, -100.0_r8_kind, 100.0_r8_kind, -200.0_r8_kind, &
                   -50.0_r8_kind, 0.0_r8_kind, 50.0_r8_kind, 200.0_r8_kind /) 

    do i = 1, size(this)
      arr = this(i)%get_attribute_value()
      is_correct = (trim(this(i)%get_attribute_name()) .eq. attr_names(i)) .and. &
                   (this(i)%get_attribute_type() .eq. real8_type) .and. & 
                   ALL(abs(this(i)%get_attribute_value() - attr_vals(i_expected:i_expected+val_inds(i))) .lt. 1.0e-8_r8_kind)

      i_expected = i_expected + val_inds(i) + 1 
      if(.not. is_correct) call mpp_error(FATAL, "incorrect metadata value")
    enddo
  end subroutine

end program