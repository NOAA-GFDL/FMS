!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!> @brief Checks the output file after running test_subregional
program check_subregional
  use fms_mod,           only: fms_init, fms_end, string
  use fms2_io_mod,       only: FmsNetcdfFile_t, read_data, close_file, open_file, get_dimension_size
  use mpp_mod,           only: mpp_npes, mpp_error, FATAL, mpp_pe
  use platform_mod,      only: r4_kind, r8_kind

  implicit none

  call fms_init()

  call check_zsubaxis_file("test_subZaxis.nc")
  ! The files are in the same subregion, one of them is defined using latlon and another one indices
  call check_subregional_file("test_subregional.nc")
  call check_subregional_file("test_subregional2.nc")

  call fms_end()

  contains

  !> @brief Check dimension data
  subroutine check_dims(err_msg, actual_data, expected_data)
    character(len=*), intent(in) :: err_msg          !< Error message to append
    real,             intent(in) :: actual_data(:)   !< Dimension data from file
    real,             intent(in) :: expected_data(:) !< Expected data

    integer :: i

    do i = 1, size(actual_data)
      if (actual_data(i) .ne. expected_data(i)) &
        call mpp_error(FATAL, "The data is not expected for "//trim(err_msg))
    enddo
  end subroutine check_dims

  !> @brief Check the data for the Z subaxis
  subroutine check_zsubaxis_file(file_name)
    character(len=*), intent(in) :: file_name !< Name of the file to check

    type(FmsNetcdfFile_t)              :: fileobj            !< FMS2 fileobj
    integer                            :: dim_size           !< dim_size as read in from the file
    real, allocatable                  :: dims(:)            !< dimension data as read in from the file
    real, allocatable                  :: dims_exp(:)        !< dimensions data expected

    if (.not. open_file(fileobj, file_name, "read")) &
      call mpp_error(FATAL, "unable to open "//trim(file_name))

    call get_dimension_size(fileobj, "z_sub01", dim_size)
    if (dim_size .ne. 3) call mpp_error(FATAL, "z_sub01 is not the correct size!")
    allocate(dims(dim_size), dims_exp(dim_size))
    call read_data(fileobj, "z_sub01", dims)
    dims_exp = (/3., 4., 5. /)
    call check_dims("z_sub01",dims, dims_exp)
    deallocate(dims, dims_exp)

    call get_dimension_size(fileobj, "z_sub02", dim_size)
    if (dim_size .ne. 2) call mpp_error(FATAL, "z_sub02 is not the correct size!")
    allocate(dims(dim_size), dims_exp(dim_size))
    call read_data(fileobj, "z_sub02", dims)
    dims_exp = (/2., 3./)
    call check_dims("z_sub01",dims, dims_exp)
    deallocate(dims, dims_exp)

    call close_file(fileobj)

  end subroutine check_zsubaxis_file

  !> @brief Check the data for the subregional file
  subroutine check_subregional_file(file_name)
    character(len=*), intent(in) :: file_name !< Name of the file to check

    type(FmsNetcdfFile_t)              :: fileobj            !< FMS2 fileobj
    integer                            :: dim_size           !< dim_size as read in from the file
    real, allocatable                  :: dims(:)            !< dimension data as read in from the file
    real, allocatable                  :: dims_exp(:)        !< dimensions data expected

    if (.not. open_file(fileobj, trim(file_name)//".0003", "read")) &
      call mpp_error(FATAL, "unable to open "//trim(file_name))

    call get_dimension_size(fileobj, "x_sub01", dim_size)
    if (dim_size .ne. 6) call mpp_error(FATAL, "x_sub01 is not the correct size!")
    allocate(dims(dim_size), dims_exp(dim_size))
    call read_data(fileobj, "x_sub01", dims)
    dims_exp = (/60., 61., 62., 63., 64., 65. /)
    call check_dims("x_sub01",dims, dims_exp)
    deallocate(dims, dims_exp)

    call get_dimension_size(fileobj, "y_sub01", dim_size)
    if (dim_size .ne. 5) call mpp_error(FATAL, "y_sub01 is not the correct size!")
    allocate(dims(dim_size), dims_exp(dim_size))
    call read_data(fileobj, "y_sub01", dims)
    dims_exp = (/60., 61., 62., 63., 64./)
    call check_dims("y_sub01",dims, dims_exp)
    deallocate(dims, dims_exp)

    call close_file(fileobj)

    if (.not. open_file(fileobj, trim(file_name)//".0004", "read")) &
      call mpp_error(FATAL, "unable to open "//trim(file_name))

    call get_dimension_size(fileobj, "x_sub01", dim_size)
    if (dim_size .ne. 6) call mpp_error(FATAL, "x_sub01 is not the correct size!")
    allocate(dims(dim_size), dims_exp(dim_size))
    call read_data(fileobj, "x_sub01", dims)
    dims_exp = (/60., 61., 62., 63., 64., 65. /)
    call check_dims("x_sub01",dims, dims_exp)
    deallocate(dims, dims_exp)

    call get_dimension_size(fileobj, "y_sub01", dim_size)
    if (dim_size .ne. 1) call mpp_error(FATAL, "y_sub01 is not the correct size!")
    allocate(dims(dim_size), dims_exp(dim_size))
    call read_data(fileobj, "y_sub01", dims)
    dims_exp = (/65./)
    call check_dims("y_sub01",dims, dims_exp)
    deallocate(dims, dims_exp)

    call close_file(fileobj)

  end subroutine check_subregional_file

end program
