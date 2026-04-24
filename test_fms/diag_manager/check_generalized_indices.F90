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

!> @brief Checker for test_generalized_indices output.
!!        Verifies permuted-axis variables match identity variables under axis permutations
program check_generalized_indices
  use fms_mod,           only: fms_init, fms_end, string
  use testing_utils,     only: check_perm
  use fms2_io_mod,       only: FmsNetcdfFile_t, read_data, open_file, close_file, get_global_attribute
  use mpp_mod,           only: mpp_error, FATAL, mpp_pe
  use platform_mod,      only: r8_kind

  implicit none

  type(FmsNetcdfFile_t)           :: fileobj
  integer                         :: nx, ny, nz
  integer                         :: i

  real(kind=r8_kind), allocatable :: var2_id(:,:)     ! (x,y)
  real(kind=r8_kind), allocatable :: var2_yx(:,:)     ! (y,x)
  real(kind=r8_kind), allocatable :: var3_id(:,:,:)   ! (x,y,z)
  real(kind=r8_kind), allocatable :: var3_zx(:,:,:)   ! (z,y,x)
  real(kind=r8_kind), allocatable :: var3_yzx(:,:,:)  ! (y,z,x)
  real(kind=r8_kind), allocatable :: var3_zxy(:,:,:)  ! (z,x,y)

  call fms_init()

  nx = 96
  ny = 96
  nz = 5

  if (.not. open_file(fileobj, "test_gen.nc", "read")) &
    call mpp_error(FATAL, "unable to open test_gen.nc")

  call check_global_attribute(fileobj, "test_generalized_indices")

  allocate(var2_id(nx,ny),    var2_yx(ny,nx))
  allocate(var3_id(nx,ny,nz), var3_zx(nz,ny,nx), var3_yzx(ny,nz,nx), var3_zxy(nz,nx,ny))

  ! Output every 6 hours over 48 hours => 8 records
  do i = 1, 8
    var2_id  = -999._r8_kind
    var2_yx  = -999._r8_kind
    var3_id  = -999._r8_kind
    var3_zx  = -999._r8_kind
    var3_yzx = -999._r8_kind
    var3_zxy = -999._r8_kind

    print *, "Checking var2_yx vs var2_id - time_level:", i
    call read_data(fileobj, "var2_id", var2_id, unlim_dim_level=i)
    call read_data(fileobj, "var2_yx", var2_yx, unlim_dim_level=i)
    call check_perm(var2_id, var2_yx, [2,1])

    print *, "Checking var3_zx vs var3_id - time_level:", i
    call read_data(fileobj, "var3_id", var3_id, unlim_dim_level=i)
    call read_data(fileobj, "var3_zx", var3_zx, unlim_dim_level=i)
    call check_perm(var3_id, var3_zx, [3,2,1])

    print *, "Checking var3_yzx vs var3_id - time_level:", i
    call read_data(fileobj, "var3_yzx", var3_yzx, unlim_dim_level=i)
    call check_perm(var3_id, var3_yzx, [2,3,1])

    print *, "Checking var3_zxy vs var3_id - time_level:", i
    call read_data(fileobj, "var3_zxy", var3_zxy, unlim_dim_level=i)
    call check_perm(var3_id, var3_zxy, [3,1,2])
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
end program check_generalized_indices
