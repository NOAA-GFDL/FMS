program test_metadata_transfer
  use fms_mod, only: fms_init, fms_end, string
  use mpp_mod
  use metadata_transfer_mod
  use platform_mod
  use mpi
  use, intrinsic :: iso_c_binding

  implicit none

  type(metadata_type) :: file_metadata(2)

  call fms_init()

  call fms_metadata_transfer_init()

  ! set metadata only on root PE
  if (mpp_pe() .eq. mpp_root_pe()) then
    call file_metadata(1)%set_attribute_name("_FillValue"//c_null_char)
    call file_metadata(1)%set_attribute_type(real8_type)
    call file_metadata(1)%set_attribute_value([-666.0_r8_kind, 666.0_r8_kind])

    call file_metadata(2)%set_attribute_name("missing_value"//c_null_char)
    call file_metadata(2)%set_attribute_type(real8_type)
    call file_metadata(2)%set_attribute_value([-100.0_r8_kind, 100.0_r8_kind])
  endif
  ! Broadcast the metadata to all PEs
  call file_metadata(1)%fms_metadata_broadcast()
  call file_metadata(2)%fms_metadata_broadcast()

  call dump_metadata(file_metadata)

  call fms_end()

  contains

  subroutine dump_metadata(this)
    type(metadata_type), intent(in) :: this(:)

    integer :: i
    do i = 1, size(this)
      print *, mpp_pe(), " knows that the attribute_name is ", trim(adjustl(this(i)%get_attribute_name())), &
        " and the attribute_type is ", string(this(i)%get_attribute_type()), &
        " and the value is ", this(i)%get_attribute_value()
    enddo
  end subroutine
end program