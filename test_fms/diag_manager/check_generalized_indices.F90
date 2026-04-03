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

!> @brief Checker for test_generalized_indicies output.
!!        Verifies swapped-axis variables match identity variables under transpose:
!!          var2_id(x,y)   == var2_swap(y,x)
!!          var3_id(x,y,z) == var3_swap(y,x,z)
program check_generalized_indices
  use fms_mod,           only: fms_init, fms_end, string
  use fms2_io_mod,       only: FmsNetcdfFile_t, read_data, open_file, close_file, get_global_attribute
  use mpp_mod,           only: mpp_error, FATAL, mpp_pe
  use platform_mod,      only: r4_kind

  implicit none

  type(FmsNetcdfFile_t)           :: fileobj
  integer                         :: nx, ny, nz
  integer                         :: i

  real(kind=r4_kind), allocatable :: var2_id(:,:)     ! (x,y)
  real(kind=r4_kind), allocatable :: var2_swap(:,:)   ! (y,x)
  real(kind=r4_kind), allocatable :: var3_id(:,:,:)   ! (x,y,z)
  real(kind=r4_kind), allocatable :: var3_swap(:,:,:) ! (y,x,z)

  call fms_init()

  nx = 96
  ny = 96
  nz = 5

  if (.not. open_file(fileobj, "test_gen.nc", "read")) &
    call mpp_error(FATAL, "unable to open test_gen.nc")

  call check_global_attribute(fileobj, "test_generalized_indices")

  allocate(var2_id(nx,ny), var2_swap(ny,nx))
  allocate(var3_id(nx,ny,nz), var3_swap(ny,nx,nz))

  ! Output every 6 hours over 48 hours => 8 records
  do i = 1, 8
    var2_id   = -999._r4_kind
    var2_swap = -999._r4_kind
    var3_id   = -999._r4_kind
    var3_swap = -999._r4_kind

    print *, "Checking var2_swap vs var2_id - time_level:", string(i)
    call read_data(fileobj, "var2_id",   var2_id,   unlim_dim_level=i)
    call read_data(fileobj, "var2_swap", var2_swap, unlim_dim_level=i)
    call check_var2_relation(var2_id, var2_swap)

    print *, "Checking var3_swap vs var3_id - time_level:", string(i)
    call read_data(fileobj, "var3_id",   var3_id,   unlim_dim_level=i)
    call read_data(fileobj, "var3_swap", var3_swap, unlim_dim_level=i)
    call check_var3_relation(var3_id, var3_swap)
  enddo

  call close_file(fileobj)
  call fms_end()

contains

  subroutine check_global_attribute(fileobj, expected_title)
    type(FmsNetcdfFile_t), intent(in) :: fileobj
    character(len=*),      intent(in) :: expected_title

    character(len=100) :: attribute_value

    call get_global_attribute(fileobj, "title", attribute_value)
    if (trim(attribute_value) .ne. trim(expected_title)) then
      call mpp_error(FATAL, "Global attribute 'title' not expected value.")
    endif
  end subroutine check_global_attribute

  subroutine check_var2_relation(v_id, v_sw)
    real(kind=r4_kind), intent(in) :: v_id(:,:)  ! (x,y)
    real(kind=r4_kind), intent(in) :: v_sw(:,:)  ! (y,x)

    integer :: x, y

    if (size(v_id,1) /= size(v_sw,2) .or. size(v_id,2) /= size(v_sw,1)) then
      call mpp_error(FATAL, "check_var2_relation: dimension mismatch between var2_id and var2_swap")
    endif

    do x = 1, size(v_id,1)
      do y = 1, size(v_id,2)
        if (abs(v_id(x,y) - v_sw(y,x)) > 0) then
          print *, mpp_pe(), "var2 mismatch at (x,y)=", x, y, " id=", v_id(x,y), " swap(y,x)=", v_sw(y,x)
          call mpp_error(FATAL, "check_var2_relation: var2_swap != transpose(var2_id)")
        endif
      enddo
    enddo
  end subroutine check_var2_relation

  subroutine check_var3_relation(v_id, v_sw)
    real(kind=r4_kind), intent(in) :: v_id(:,:,:)  ! (x,y,z)
    real(kind=r4_kind), intent(in) :: v_sw(:,:,:)  ! (y,x,z)

    integer :: x, y, z

    if (size(v_id,1) /= size(v_sw,2) .or. size(v_id,2) /= size(v_sw,1) .or. size(v_id,3) /= size(v_sw,3)) then
      call mpp_error(FATAL, "check_var3_relation: dimension mismatch between var3_id and var3_swap")
    endif

    do x = 1, size(v_id,1)
      do y = 1, size(v_id,2)
        do z = 1, size(v_id,3)
          if (abs(v_id(x,y,z) - v_sw(y,x,z)) > 0) then
            print *, mpp_pe(), "var3 mismatch at (x,y,z)=", x, y, z, &
                               " id=", v_id(x,y,z), " swap(y,x,z)=", v_sw(y,x,z)
            call mpp_error(FATAL, "check_var3_relation: var3_swap != var3_id with x/y swapped")
          endif
        enddo
      enddo
    enddo
  end subroutine check_var3_relation

end program check_generalized_indices

