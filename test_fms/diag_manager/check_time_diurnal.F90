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
!! TODO more complicated cases and data
!> @brief  Checks the output file after running test_reduction_methods using the "time_diurnal" reduction method
program check_time_diurnal
  use fms_mod,           only: fms_init, fms_end, string
  use fms2_io_mod,       only: FmsNetcdfFile_t, read_data, close_file, open_file
  use mpp_mod,           only: mpp_npes, mpp_error, FATAL, mpp_pe, input_nml_file, NOTE
  use platform_mod,      only: r4_kind, r8_kind
  use testing_utils,     only: allocate_buffer, test_normal, test_openmp, test_halos, no_mask, logical_mask, real_mask

  implicit none

  type(FmsNetcdfFile_t)              :: fileobj            !< FMS2 fileobj
  type(FmsNetcdfFile_t)              :: fileobj1           !< FMS2 fileobj for subregional file 1
  type(FmsNetcdfFile_t)              :: fileobj2           !< FMS2 fileobj for subregional file 2
  real(kind=r4_kind), allocatable    :: cdata_out(:,:,:,:) !< Data in the compute domain
  real(kind=r4_kind), allocatable    :: cdata_out_5d(:,:,:,:,:) !< Data in the compute domain
  integer                            :: nx                 !< Number of points in the x direction
  integer                            :: ny                 !< Number of points in the y direction
  integer                            :: nz                 !< Number of points in the z direction
  integer                            :: nw                 !< Number of points in the w direction
  integer                            :: nd                 !< Number of points in the diurnal axis
  integer                            :: ti                  !< For looping through time levels
  integer                            :: io_status          !< Io status after reading the namelist
  logical                            :: use_mask           !< .true. if using masks
  integer, parameter :: file_freq = 6 !< file frequency as set in diag_table.yaml

  integer :: test_case = test_normal !< Indicates which test case to run
  integer :: mask_case = no_mask     !< Indicates which masking option to run
  integer, parameter :: kindl = KIND(0.0) !< compile-time default kind size
  integer :: nmonths !< number of months the test ran for
  namelist / test_diag_diurnal_nml / test_case, mask_case

  call fms_init()

  read (input_nml_file, test_diag_diurnal_nml, iostat=io_status)

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
  nmonths = 3
  nd = 3 !< diurnal sample size

  if (.not. open_file(fileobj, "test_diurnal.nc", "read")) &
    call mpp_error(FATAL, "unable to open test_diurnal.nc")

  if (.not. open_file(fileobj1, "test_diurnal_regional.nc.0004", "read")) &
    call mpp_error(FATAL, "unable to open test_diurnal_regional.nc.0004")

  if (.not. open_file(fileobj2, "test_diurnal_regional.nc.0005", "read")) &
    call mpp_error(FATAL, "unable to open test_diurnal_regional.nc.0005")

  !cdata_out = allocate_buffer(1, nx, 1, ny, nz, nd)
  allocate(cdata_out(nx, ny, nz, nd))
  allocate(cdata_out_5d(nx, ny, nz, nw, nd))

  do ti = 1, nmonths
    cdata_out = -999_r4_kind
    print *, "Checking answers for var1 - time_level:", string(ti)
    call read_data(fileobj, "var1", cdata_out(:,1:nd,1,1), unlim_dim_level=ti)
    call check_data_1d(cdata_out(:,1:nd,1,1), ti, sample_size=nd)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var2 - time_level:", string(ti)
    call read_data(fileobj, "var2", cdata_out(:,:,1:nd,1), unlim_dim_level=ti)
    call check_data_2d(cdata_out(:,:,1:nd,1), ti, sample_size=nd)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var3 - time_level:", string(ti)
    call read_data(fileobj, "var3", cdata_out(:,:,:,:), unlim_dim_level=ti)
    call check_data_3d(cdata_out(:,:,:,:), ti, .false., sample_size=nd)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var4_diurnal - time_level:", string(ti)
    call read_data(fileobj, "var4", cdata_out_5d, unlim_dim_level=ti)
    call check_data_4d(cdata_out_5d(:,:,:,:,:), ti, .false., sample_size=nd)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var3_diurnal in the first regional file- time_level:", string(ti)
    call read_data(fileobj1, "var3_diurnal", cdata_out(1:4,1:3,1:2,1:1), unlim_dim_level=ti)
    call check_data_3d(cdata_out(1:4,1:3,1:2,1:1), ti, .true., sample_size=nd, nx_offset=77, ny_offset=77, nz_offset=1)

    cdata_out = -999_r4_kind
    print *, "Checking answers for var3_diurnal in the second regional file- time_level:", string(ti)
    call read_data(fileobj2, "var3_diurnal", cdata_out(1:4,1:1,1:2,1:1), unlim_dim_level=ti)
    call check_data_3d(cdata_out(1:4,1:1,1:2,1:1), ti, .true., sample_size=nd, nx_offset=77, ny_offset=80, nz_offset=1)
  enddo

  call fms_end()

contains


  !> @brief Check that the 2d data read in is correct
  subroutine check_data_2d(buffer, time_level, sample_size)
    real(kind=r4_kind), intent(in)    :: buffer(:,:,:)   !< Buffer read from the table (2d + the diurnal axis)
    integer,            intent(in)    :: time_level    !< Time level read in
    integer,            intent(in)    :: sample_size   !< diurnal sample size of variable to check
    real(kind=r4_kind)                :: buffer_exp    !< Expected result

    integer :: ii,i, j, k, l, d!< For looping
    integer :: step_avg !< avg of time step increments to use in generating reference data
    integer :: d_index
    real(r8_kind) :: hrly_sums(sample_size)

    ! sum of hours in diurnal section
    hrly_sums = 0
    do i=1, 23
      d_index = i / (24/sample_size) + 1
      hrly_sums(d_index) = hrly_sums(d_index) + i
    enddo
    hrly_sums = hrly_sums / (24/sample_size)

    ! 2d answer is the
    do ii = 1, size(buffer, 1)
      do j = 1, size(buffer, 2)
        do d = 1, sample_size
          buffer_exp = hrly_sums(d)
          if (use_mask .and. ii .eq. 1 .and. j .eq. 1) buffer_exp = -666_r4_kind
          if (abs(buffer(ii, j, d) - buffer_exp) > 0.0) then
            print *, "indices:", ii, j, d, "expected:", buffer_exp, "read in:",buffer(ii, j, d)
            call mpp_error(FATAL, "Check_time_diurnal::check_data_2d:: Data is not correct")
          endif
        enddo
      enddo
    enddo
  end subroutine check_data_2d

  !> @brief Check that the 3d data read in is correct
  subroutine check_data_3d(buffer, time_level, is_regional, sample_size, nx_offset, ny_offset, nz_offset)
    real(kind=r4_kind), intent(in)    :: buffer(:,:,:,:) !< Buffer read from the table
    integer,            intent(in)    :: time_level    !< Time level read in
    logical,            intent(in)    :: is_regional   !< .True. if the variable is subregional
    integer, intent(in)               :: sample_size   !< diurnal sample size
    real(kind=r4_kind)                :: buffer_exp    !< Expected result
    integer, optional,  intent(in)    :: nx_offset     !< Offset in the x direction
    integer, optional,  intent(in)    :: ny_offset     !< Offset in the y direction
    integer, optional,  intent(in)    :: nz_offset     !< Offset in the z direction

    integer :: ii, i, j, k, l, d!< For looping
    integer :: nx_oset !< Offset in the x direction (local variable)
    integer :: ny_oset !< Offset in the y direction (local variable)
    integer :: nz_oset !< Offset in the z direction (local variable)
    integer :: step_avg!< avg of time step increments to use in generating reference data
    real(r8_kind) :: hrly_sums(24/sample_size) !< can i even do this (yes)
    integer :: d_index !< diurnal index

    ! data is just the hour it was sent at
    ! sum of hours in each diurnal section
    hrly_sums = 0
    do i=1, 23
      d_index = i / (24/sample_size) + 1
      hrly_sums(d_index) = hrly_sums(d_index) + i
    enddo
    hrly_sums = hrly_sums / (24/sample_size)

    nx_oset = 0
    if (present(nx_offset)) nx_oset = nx_offset

    ny_oset = 0
    if (present(ny_offset)) ny_oset = ny_offset

    nz_oset = 0
    if (present(nz_offset)) nz_oset = nz_offset

    ! 3d answer is
    !
    do ii = 1, size(buffer, 1)
      do j = 1, size(buffer, 2)
        do k = 1, size(buffer, 3)
          do d=1, size(buffer, 4)
            buffer_exp = hrly_sums(d)
            if (use_mask .and. ii .eq. 1 .and. j .eq. 1 .and. k .eq. 1 .and. .not. is_regional) &
              buffer_exp = -666_r4_kind
            if (abs(buffer(ii, j, k, d) - buffer_exp) > 0.0) then
              print *, mpp_pe(),'indices:',ii, j, k, d, "read in:", buffer(ii, j, k, d), "expected:",buffer_exp
              call mpp_error(FATAL, "Check_time_diurnal::check_data_3d:: Data is not correct")
            endif
          enddo
        enddo
      enddo
    enddo
  end subroutine check_data_3d

  !> @brief Check that the 1d data read in is correct
  subroutine check_data_1d(buffer, time_level, sample_size)
    real(kind=r4_kind), intent(in)    :: buffer(:,:)   !< Buffer read from the table
    integer,            intent(in)    :: time_level    !< Time level read in
    integer,            intent(in)    :: sample_size   !< diurnal sample size of variable to check
    real(kind=r4_kind)                :: buffer_exp    !< Expected result
    integer :: ii,i, j, k, l, d!< For looping
    integer :: step_avg !< avg of time step increments to use in generating reference data
    integer :: d_index
    real(r8_kind) :: hrly_sums(sample_size)

    ! sum of hours in diurnal section
    hrly_sums = 0
    do i=1, 23
      d_index = i / (24/sample_size) + 1
      hrly_sums(d_index) = hrly_sums(d_index) + i
    enddo
    hrly_sums = hrly_sums / (24/sample_size)

    do ii = 1, size(buffer, 1)
      do d = 1, sample_size
        buffer_exp = hrly_sums(d)
        if (use_mask .and. ii .eq. 1) buffer_exp = -666_r4_kind
        if (abs(buffer(ii,d) - buffer_exp) > 0.0) then
          print *, "indices:", ii, d, "expected:", buffer_exp, "read in:",buffer(ii,d)
          call mpp_error(FATAL, "Check_time_diurnal::check_data_1d:: Data is not correct")
        endif
      enddo
    enddo
  end subroutine check_data_1d

  !> @brief Check that the 4d data read in is correct
  subroutine check_data_4d(buffer, time_level, is_regional, sample_size, nx_offset, ny_offset, nz_offset)
    real(kind=r4_kind), intent(in)    :: buffer(:,:,:,:,:) !< Buffer read from the table
    integer,            intent(in)    :: time_level    !< Time level read in
    logical,            intent(in)    :: is_regional   !< .True. if the variable is subregional
    integer, intent(in)               :: sample_size   !< diurnal sample size
    real(kind=r4_kind)                :: buffer_exp    !< Expected result
    integer, optional,  intent(in)    :: nx_offset     !< Offset in the x direction
    integer, optional,  intent(in)    :: ny_offset     !< Offset in the y direction
    integer, optional,  intent(in)    :: nz_offset     !< Offset in the z direction

    integer :: ii, i, j, k, l, d, w!< For looping
    integer :: nx_oset !< Offset in the x direction (local variable)
    integer :: ny_oset !< Offset in the y direction (local variable)
    integer :: nz_oset !< Offset in the z direction (local variable)
    integer :: step_avg!< avg of time step increments to use in generating reference data
    real(r8_kind) :: hrly_sums(24/sample_size) !< calculated hourly sums for each diurnal section
    integer :: d_index !< diurnal index

    ! data is just the hour it was sent at
    ! sum of hours in each diurnal section
    hrly_sums = 0
    do i=1, 23
      d_index = i / (24/sample_size) + 1
      hrly_sums(d_index) = hrly_sums(d_index) + i
    enddo
    hrly_sums = hrly_sums / (24/sample_size)

    nx_oset = 0
    if (present(nx_offset)) nx_oset = nx_offset

    ny_oset = 0
    if (present(ny_offset)) ny_oset = ny_offset

    nz_oset = 0
    if (present(nz_offset)) nz_oset = nz_offset

    do ii = 1, size(buffer, 1)
      do j = 1, size(buffer, 2)
        do k = 1, size(buffer, 3)
          do w = 1, size(buffer, 4)
            do d = 1, sample_size
              buffer_exp = hrly_sums(d)
              if (use_mask .and. ii .eq. 1 .and. j .eq. 1 .and. k .eq. 1 .and. &
                 .not. is_regional) then
                buffer_exp = -666_r4_kind
              endif
              if (abs(buffer(ii, j, k, w, d) - buffer_exp) > 0.0) then
                print *, mpp_pe(),'indices:',ii, j, k, w, d, "read in:", buffer(ii, j, k, w, d), "expected:",buffer_exp
                call mpp_error(FATAL, "Check_time_diurnal::check_data_4d:: Data is not correct")
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine check_data_4d
end program
