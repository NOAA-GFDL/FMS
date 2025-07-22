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

program test_generalized_indices
use fms_mod,              only: fms_init, fms_end
use mpp_mod,              only: mpp_pe
use fms_string_utils_mod, only: string
use fms2_io_mod,          only: FmsNetcdfDomainFile_t, register_field, register_axis, unlimited, &
                                read_data
use fms2_io_mod,          only: open_file, close_file, write_data
use mpp_domains_mod,      only: mpp_define_domains, mpp_define_io_domain, mpp_get_compute_domain, &
                                domain2d, center
use random_numbers_mod,   only: randomNumberStream, initializeRandomNumberStream, getRandomNumbers
use platform_mod,         only: r8_kind

implicit none

integer, parameter :: nx = 16
integer, parameter :: ny = 16
integer, parameter :: nz = 10
integer, parameter :: nt = 5
integer, dimension(2), parameter :: layout = [4, 4]
integer, dimension(2), parameter :: io_layout = [2, 2]
character(*), parameter :: var_name = "rand"
character(1), parameter :: axis_labels(3) = ["x", "y", "z"]

type (domain2d) :: domain

real(r8_kind), dimension(:, :, :, :, :), allocatable :: data
integer :: lb(3), ub(3)

call fms_init

call define_domain

call test_permutation(1, 2, 3)
call test_permutation(1, 3, 2)
call test_permutation(2, 1, 3)
call test_permutation(2, 3, 1)
call test_permutation(3, 1, 2)
call test_permutation(3, 2, 1)

call fms_end

contains

subroutine define_domain
  integer :: isc, iec, jsc, jec !< Compute domain

  call mpp_define_domains( [1, nx, 1, ny], layout, domain)
  call mpp_define_io_domain(domain, io_layout)
  call mpp_get_compute_domain(domain, isc, iec, jsc, jec)

  lb = [isc, jsc, 1]
  ub = [iec, jec, nz]
end subroutine define_domain

subroutine set_data(arr)
  real(r8_kind), intent(out) :: arr(:, :, :, :)
  type(randomNumberStream) :: random_stream
  integer :: i3, i4, n3, n4

  random_stream = initializeRandomNumberStream(mpp_pe())

  n4 = size(arr, 4)
  n3 = size(arr, 3)

  do i4=1, n4
    do i3=1, n3
      call getRandomNumbers(random_stream, arr(:, :, i3, i4))
    enddo
  enddo
end subroutine set_data

subroutine test_permutation(n1, n2, n3)
  integer, intent(in) :: n1, n2, n3
  character(1) :: axes(4)

  axes(1) = axis_labels(n1)
  axes(2) = axis_labels(n2)
  axes(3) = axis_labels(n3)
  axes(4) = "t"

  allocate(data(lb(n1):ub(n1), lb(n2):ub(n2), lb(n3):ub(n3), nt, 2))

  call set_data(data(:, :, :, :, 1))
  call write_permutation(axes, data(:, :, :, :, 1))
  call read_permutation(axes, data(:, :, :, :, 2))

  if (any(data(:, :, :, :, 1).ne.data(:, :, :, :, 2))) then
    print "(A)", "The data written out do not match the data read back in"
    stop 1
  endif

  deallocate(data)
end subroutine test_permutation

subroutine open_netcdf_file(fileobj, axes, mode)
  type(FmsNetcdfDomainFile_t), intent(out) :: fileobj
  character(1), intent(in) :: axes(4)
  character(*), intent(in) :: mode
  character(:), allocatable :: filename

  filename = "permutation_" // axes(1) // axes(2) // axes(3) // ".nc"

  if (.not.open_file(fileobj, filename, mode, domain)) then
    print "(A)", "Error: Could not open "//filename//" in "//mode//" mode"
    stop 1
  endif

  call register_axis(fileobj, "x", "x")
  call register_axis(fileobj, "y", "y")
  call register_axis(fileobj, "z", nz)
  call register_axis(fileobj, "t", unlimited)
end subroutine open_netcdf_file

subroutine write_permutation(axes, arr)
  character(1), intent(in) :: axes(4)
  real(r8_kind), intent(in) :: arr(:, :, :, :)
  type(FmsNetcdfDomainFile_t) :: fileobj

  call open_netcdf_file(fileobj, axes, "write")

  call register_field(fileobj, var_name, "double", axes)
  call write_data(fileobj, var_name, arr)
  call close_file(fileobj)
end subroutine write_permutation

subroutine read_permutation(axes, arr)
  character(1), intent(in) :: axes(4)
  real(r8_kind), intent(out) :: arr(:, :, :, :)
  type(FmsNetcdfDomainFile_t) :: fileobj

  call open_netcdf_file(fileobj, axes, "read")
  call read_data(fileobj, var_name, arr)
  call close_file(fileobj)
end subroutine read_permutation

end program test_generalized_indices
