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

!> @brief  Checks the output file after running test_reduction_methods using the "min" reduction method
program check_time_min
  use fms_mod,           only: fms_init, fms_end, string
  use fms2_io_mod,       only: FmsNetcdfFile_t, read_data, close_file, open_file
  use mpp_mod,           only: mpp_npes, mpp_error, FATAL, mpp_pe, input_nml_file
  use platform_mod,      only: r4_kind, r8_kind
  use testing_utils,     only: allocate_buffer, test_normal, test_openmp, test_halos, no_mask, logical_mask, real_mask

  implicit none

  type(FmsNetcdfFile_t)              :: fileobj            !< FMS2 fileobj
  type(FmsNetcdfFile_t)              :: fileobj1           !< FMS2 fileobj for subregional file 1
  type(FmsNetcdfFile_t)              :: fileobj2           !< FMS2 fileobj for subregional file 2
  real(kind=r4_kind), allocatable    :: cdata_out(:,:,:,:) !< Data in the compute domain
  integer                            :: nx                 !< Number of points in the x direction
  integer                            :: ny                 !< Number of points in the y direction
  integer                            :: nz                 !< Number of points in the z direction
  integer                            :: nw                 !< Number of points in the 4th dimension
  integer                            :: i                  !< For looping
  integer                            :: io_status          !< Io status after reading the namelist
  logical                            :: use_mask           !< .true. if using masks

  integer :: test_case = test_normal !< Indicates which test case to run
  integer :: mask_case = no_mask     !< Indicates which masking option to run

  namelist / test_reduction_methods_nml / test_case, mask_case

  call fms_init()

  read (input_nml_file, test_reduction_methods_nml, iostat=io_status)

  select case(mask_case)
  case (no_mask)
    use_mask = .false.
  case (logical_mask, real_mask)
    use_mask = .true.
  end select
  nx = 96
  ny = 96
  nz = 5
  nw = 2

  if (.not. open_file(fileobj, "test_min.nc", "read")) &
    call mpp_error(FATAL, "unable to open file")

  if (.not. open_file(fileobj1, "test_min_regional.nc.0004", "read")) &
    call mpp_error(FATAL, "unable to open file")

  if (.not. open_file(fileobj2, "test_min_regional.nc.0005", "read")) &
    call mpp_error(FATAL, "unable to open file")

  cdata_out = allocate_buffer(1, nx, 1, ny, nz, nw)

  do i = 1, 8
    cdata_out = -999_r4_kind
    print *, "Checking answers for var0_min - time_level:", string(i)
    call read_data(fileobj, "var0_min", cdata_out(1,1,1,1), unlim_dim_level=i)
    call check_data_0d(cdata_out(1,1,1,1), i)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var1_min - time_level:", string(i)
    call read_data(fileobj, "var1_min", cdata_out(:,1,1,1), unlim_dim_level=i)
    call check_data_1d(cdata_out(:,1,1,1), i)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var2_min - time_level:", string(i)
    call read_data(fileobj, "var2_min", cdata_out(:,:,1,1), unlim_dim_level=i)
    call check_data_2d(cdata_out(:,:,1,1), i)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var3_min - time_level:", string(i)
    call read_data(fileobj, "var3_min", cdata_out(:,:,:,1), unlim_dim_level=i)
    call check_data_3d(cdata_out(:,:,:,1), i, .false.)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var3_Z_min - time_level:", string(i)
    call read_data(fileobj, "var3_Z_min", cdata_out(:,:,1:2,1), unlim_dim_level=i)
    call check_data_3d(cdata_out(:,:,1:2,1), i, .true., nz_offset=1)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var3_min in the first regional file- time_level:", string(i)
    call read_data(fileobj1, "var3_min", cdata_out(1:4,1:3,1:2,1), unlim_dim_level=i)
    call check_data_3d(cdata_out(1:4,1:3,1:2,1), i, .true., nx_offset=77, ny_offset=77, nz_offset=1)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var3_min in the second regional file- time_level:", string(i)
    call read_data(fileobj2, "var3_min", cdata_out(1:4,1:1,1:2,1), unlim_dim_level=i)
    call check_data_3d(cdata_out(1:4,1:1,1:2,1), i, .true., nx_offset=77, ny_offset=80, nz_offset=1)
  enddo

  call fms_end()

contains

  !> @brief Check that the 0d data read in is correct
  subroutine check_data_0d(buffer, time_level)
    real(kind=r4_kind), intent(inout) :: buffer        !< Buffer read from the table
    integer,            intent(in)    :: time_level    !< Time level read in

    real(kind=r4_kind)                :: buffer_exp    !< Expected result

    buffer_exp = real(1000_r8_kind+10_r8_kind+1_r8_kind + &
      real(6*(time_level-1)+1, kind=r8_kind)/100_r8_kind, kind=r4_kind)

    if (abs(buffer - buffer_exp) > 0) then
      print *, mpp_pe(), time_level, buffer, buffer_exp
      call mpp_error(FATAL, "Check_time_min::check_data_0d:: Data is not correct")
    endif
  end subroutine check_data_0d

  !> @brief Check that the 1d data read in is correct
  subroutine check_data_1d(buffer, time_level)
    real(kind=r4_kind), intent(in)    :: buffer(:)     !< Buffer read from the table
    integer,            intent(in)    :: time_level    !< Time level read in
    real(kind=r4_kind)                :: buffer_exp    !< Expected result

    integer ii, j, k, l !< For looping

    do ii = 1, size(buffer, 1)
      buffer_exp = real(real(ii, kind=r8_kind)* 1000_r8_kind+10_r8_kind+1_r8_kind + &
                   real(6*(time_level-1)+1, kind=r8_kind)/100_r8_kind, kind=r4_kind)
      if (use_mask .and. ii .eq. 1) buffer_exp = -666_r4_kind
      if (abs(buffer(ii) - buffer_exp) > 0) then
        print *, mpp_pe(), ii, buffer(ii), buffer_exp
        call mpp_error(FATAL, "Check_time_min::check_data_1d:: Data is not correct")
      endif
    enddo
  end subroutine check_data_1d

  !> @brief Check that the 2d data read in is correct
  subroutine check_data_2d(buffer, time_level)
    real(kind=r4_kind), intent(in)    :: buffer(:,:)   !< Buffer read from the table
    integer,            intent(in)    :: time_level    !< Time level read in
    real(kind=r4_kind)                :: buffer_exp    !< Expected result

    integer ii, j, k, l !< For looping

    do ii = 1, size(buffer, 1)
      do j = 1, size(buffer, 2)
        buffer_exp = real(real(ii, kind=r8_kind)* 1000_r8_kind+ &
                     10_r8_kind*real(j, kind=r8_kind)+1_r8_kind + &
                     real(6*(time_level-1)+1, kind=r8_kind)/100_r8_kind, kind=r4_kind)
        if (use_mask .and. ii .eq. 1 .and. j .eq. 1) buffer_exp = -666_r4_kind
        if (abs(buffer(ii, j) - buffer_exp) > 0) then
          print *, mpp_pe(), ii, j, buffer(ii, j), buffer_exp
          call mpp_error(FATAL, "Check_time_min::check_data_2d:: Data is not correct")
        endif
      enddo
    enddo
  end subroutine check_data_2d

  !> @brief Check that the 3d data read in is correct
  subroutine check_data_3d(buffer, time_level, is_regional, nx_offset, ny_offset, nz_offset)
    real(kind=r4_kind), intent(in)    :: buffer(:,:,:) !< Buffer read from the table
    integer,            intent(in)    :: time_level    !< Time level read in
    logical,            intent(in)    :: is_regional   !< .True. if the variable is subregional
    real(kind=r4_kind)                :: buffer_exp    !< Expected result
    integer, optional,  intent(in)    :: nx_offset     !< Offset in the x direction
    integer, optional,  intent(in)    :: ny_offset     !< Offset in the y direction
    integer, optional,  intent(in)    :: nz_offset     !< Offset in the z direction

    integer :: ii, j, k, l !< For looping
    integer :: nx_oset !< Offset in the x direction (local variable)
    integer :: ny_oset !< Offset in the y direction (local variable)
    integer :: nz_oset !< Offset in the z direction (local variable)

    nx_oset = 0
    if (present(nx_offset)) nx_oset = nx_offset

    ny_oset = 0
    if (present(ny_offset)) ny_oset = ny_offset

    nz_oset = 0
    if (present(nz_offset)) nz_oset = nz_offset

    do ii = 1, size(buffer, 1)
      do j = 1, size(buffer, 2)
        do k = 1, size(buffer, 3)
          buffer_exp = real(real(ii+nx_oset, kind=r8_kind)* 1000_r8_kind + &
                       10_r8_kind*real(j+ny_oset, kind=r8_kind) + &
                       1_r8_kind*real(k+nz_oset, kind=r8_kind) + &
                       real(6*(time_level-1)+1, kind=r8_kind)/100_r8_kind, kind=r4_kind)
          if (use_mask .and. ii .eq. 1 .and. j .eq. 1 .and. k .eq. 1 .and. .not. is_regional) buffer_exp = -666_r4_kind
          if (abs(buffer(ii, j, k) - buffer_exp) > 0) then
            print *, mpp_pe(), ii, j, k, buffer(ii, j, k), buffer_exp
            call mpp_error(FATAL, "Check_time_min::check_data_3d:: Data is not correct")
          endif
        enddo
      enddo
    enddo
  end subroutine check_data_3d
end program