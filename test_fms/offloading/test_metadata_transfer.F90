!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
program test_metadata_transfer
  use fms_mod, only: fms_init, fms_end, string
  use mpp_mod
  use metadata_transfer_mod
  use platform_mod
  use, intrinsic :: iso_c_binding

  implicit none

  class(metadata_class), allocatable :: file_metadata(:)

  logical :: debug = .true.

  call fms_init()

  allocate(metadata_r8_type :: file_metadata(3))
  ! all PEs need to initialize the metadata object with a datatype
  call file_metadata(1)%fms_metadata_transfer_init(real8_type)
  call file_metadata(2)%fms_metadata_transfer_init(real8_type)
  call file_metadata(3)%fms_metadata_transfer_init(real8_type)
  ! set metadata only on root PE
  if (mpp_pe() .eq. mpp_root_pe()) then
    call file_metadata(1)%set_attribute_name("_FillValue"//c_null_char)
    select type(obj => file_metadata(1))
    type is(metadata_r8_type)
      call obj%set_attribute_value([666.0_r8_kind])
    end select
    call file_metadata(2)%set_attribute_name("missing_value"//c_null_char)
    select type(obj => file_metadata(2))
    type is(metadata_r8_type)
      call obj%set_attribute_value([-100.0_r8_kind, 100.0_r8_kind])
    end select
    call file_metadata(3)%set_attribute_name("a_third_name"//c_null_char)
    select type(obj => file_metadata(3))
    type is(metadata_r8_type)
      call obj%set_attribute_value([-200.0_r8_kind, -50.0_r8_kind, 0.0_r8_kind, 50.0_r8_kind, 200.0_r8_kind])
    end select
  endif
  ! Broadcast the metadata to all PEs
  call fms_metadata_broadcast_all(file_metadata)
  ! Check data on all PEs
  select type(file_metadata)
  type is(metadata_r8_type)
      if(debug) call dump_metadata_r8(file_metadata)
      call check_metadata_r8(file_metadata)
  end select
  call mpp_sync()
  deallocate(file_metadata) !! fails here on every PE besides root

  !! test with real4_type metadata
  allocate(metadata_r4_type :: file_metadata(1))
  call file_metadata(1)%fms_metadata_transfer_init(real4_type)
  if(mpp_pe() .eq. mpp_root_pe()) then
    call file_metadata(1)%set_attribute_name("Valuez_r4"//c_null_char)
    select type(obj => file_metadata(1))
    type is(metadata_r4_type)
      call obj%set_attribute_value([666.0_r4_kind, -100.0_r4_kind, 100.0_r4_kind, -200.0_r4_kind, &
                                    -50.0_r4_kind, 0.0_r4_kind, 50.0_r4_kind, 200.0_r4_kind])
    end select
  endif
  call fms_metadata_broadcast_all(file_metadata)
  select type(file_metadata)
  type is(metadata_r4_type)
      print *, "PE: ", mpp_pe(), " metadata name is ", trim(adjustl(file_metadata(1)%get_attribute_name()))
      print *, "PE: ", mpp_pe(), " metadata value is ", file_metadata(1)%get_attribute_value()
      call check_metadata_r4(file_metadata)
  end select
  deallocate(file_metadata)

  !! test with int4_type metadata
  allocate(metadata_i4_type :: file_metadata(1))
  call file_metadata(1)%fms_metadata_transfer_init(int4_type)
  if(mpp_pe() .eq. mpp_root_pe()) then
    call file_metadata(1)%set_attribute_name("Valuez_int4"//c_null_char)
    select type(obj => file_metadata(1))
    type is(metadata_i4_type)
      call obj%set_attribute_value([666, -100, 100, -200, &
                                    -50, 0, 50, 200])
    end select
  endif
  call fms_metadata_broadcast_all(file_metadata)
  select type(file_metadata)
  type is(metadata_i4_type)
      print *, "PE: ", mpp_pe(), " metadata name is ", trim(adjustl(file_metadata(1)%get_attribute_name()))
      print *, "PE: ", mpp_pe(), " metadata value is ", file_metadata(1)%get_attribute_value()
  end select
  deallocate(file_metadata)

  !! test with int8_type metadata
  allocate(metadata_i8_type :: file_metadata(1))
  call file_metadata(1)%fms_metadata_transfer_init(int8_type)
  if(mpp_pe() .eq. mpp_root_pe()) then
    call file_metadata(1)%set_attribute_name("Valuez_int8"//c_null_char)
    select type(obj => file_metadata(1))
    type is(metadata_i8_type)
      call obj%set_attribute_value([666_i8_kind, -100_i8_kind, 100_i8_kind, -200_i8_kind, &
                                    -50_i8_kind, 0_i8_kind, 50_i8_kind, 200_i8_kind])
    end select
  endif
  call fms_metadata_broadcast_all(file_metadata)
  select type(file_metadata)
  type is(metadata_i8_type)
      print *, "PE: ", mpp_pe(), " metadata name is ", trim(adjustl(file_metadata(1)%get_attribute_name()))
      print *, "PE: ", mpp_pe(), " metadata value is ", file_metadata(1)%get_attribute_value()
  end select
  deallocate(file_metadata)

  !! test with str_type metadata
  allocate(metadata_str_type :: file_metadata(1))
  call file_metadata(1)%fms_metadata_transfer_init(int4_type)
  if(mpp_pe() .eq. mpp_root_pe()) then
    call file_metadata(1)%set_attribute_name("foo"//c_null_char)
    select type(obj => file_metadata(1))
    type is(metadata_str_type)
      call obj%set_attribute_value("bar")
    end select
  endif
  call fms_metadata_broadcast_all(file_metadata)
  select type(file_metadata)
  type is(metadata_str_type)
      print *, "PE: ", mpp_pe(), " metadata name is ", trim(adjustl(file_metadata(1)%get_attribute_name()))
      print *, "PE: ", mpp_pe(), " metadata value is ", file_metadata(1)%get_attribute_value()
      if(trim(file_metadata(1)%get_attribute_name()) .ne. "foo"//c_null_char .or. &
         trim(file_metadata(1)%get_attribute_value()) .ne. "bar") then
        call mpp_error(FATAL, "incorrect metadata name")
      endif
  end select
  deallocate(file_metadata)

  call fms_end()

  contains

  subroutine dump_metadata_r8(this)
    type(metadata_r8_type), intent(inout) :: this(:)
    real(r8_kind), allocatable :: arr(:)

    integer :: i
    do i = 1, size(this)
      arr = this(i)%get_attribute_value()
      print *, "pe: ", mpp_pe(), "i: ", i, " attribute_name is ", trim(adjustl(this(i)%get_attribute_name()))
      print *, "pe: ", mpp_pe(), "i: ", i, " attribute_value is ", arr
    enddo
  end subroutine

  subroutine check_metadata_r8(this)
    type(metadata_r8_type), intent(inout) :: this(:)
    real(r8_kind), allocatable :: arr(:)
    character(len=32) :: attr_names(3)
    real(r8_kind) :: attr_vals(8)
    integer :: i, j, last_j =1

    attr_names(1) = "_FillValue"//c_null_char
    attr_names(2) = "missing_value"//c_null_char
    attr_names(3) = "a_third_name"//c_null_char
    attr_vals = (/ 666.0_r8_kind, -100.0_r8_kind, 100.0_r8_kind, -200.0_r8_kind, &
                   -50.0_r8_kind, 0.0_r8_kind, 50.0_r8_kind, 200.0_r8_kind /)

    do i = 1, size(this)
      arr = this(i)%get_attribute_value()
      if (trim(this(i)%get_attribute_name()) .ne. attr_names(i)) then
        call mpp_error(FATAL, "incorrect metadata name")
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

  subroutine check_metadata_r4(this)
    type(metadata_r4_type), intent(inout) :: this(:)
    real(r4_kind), allocatable :: arr(:)
    character(len=32) :: attr_name
    real(r4_kind) :: attr_vals(8)
    integer :: j, last_j =1

    attr_name = "Valuez_r4"//c_null_char
    attr_vals = (/ 666.0_r4_kind, -100.0_r4_kind, 100.0_r4_kind, -200.0_r4_kind, &
                   -50.0_r4_kind, 0.0_r4_kind, 50.0_r4_kind, 200.0_r4_kind /)

    arr = this(1)%get_attribute_value()
    if (trim(this(1)%get_attribute_name()) .ne. attr_name) then
      print *, "got ", trim(this(1)%get_attribute_name()), " expected ", trim(attr_name)
      call mpp_error(FATAL, "incorrect metadata name")
    endif

    do j=1, size(arr)
      if (arr(j) .ne. attr_vals(j)) then
          print *, "got ", arr(j), " expected ", attr_vals(last_j)
          call mpp_error(FATAL, "incorrect metadata value")
        endif
    enddo
  end subroutine

end program
