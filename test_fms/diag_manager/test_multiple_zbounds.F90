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

program test_multiple_zbounds
    use fms_mod, only: fms_init, fms_end
    use fms2_io_mod
    use mpp_mod
    use time_manager_mod, only: time_type, set_calendar_type, set_date, JULIAN, set_time, OPERATOR(+)
    use diag_manager_mod

    implicit none

    type(time_type)                    :: Time             !< Time of the simulation
    type(time_type)                    :: Time_step        !< Time_step of the simulation
    integer                            :: id_z1            !< Axis id for the z dimension
    integer                            :: id_z2            !< Axis id for the z dimension
    integer                            :: id_z_reverse     !< Axis id for the z dimension (decreasing)
    integer                            :: id_var1          !< var id for the first variable
    integer                            :: id_var2          !< var_id for the second variable
    integer                            :: id_var3          !< var_id for the third variable
    real,                  allocatable :: z(:)             !< z axis data
    real,                  allocatable :: z_reverse(:)     !< z axis data (decreasing)
    integer                            :: nz               !< Size of the z dimension
    integer                            :: i
    logical                            :: used             !< Dummy argument to send_data

    call fms_init
    call set_calendar_type(JULIAN)
    call diag_manager_init

    nz = 10
    allocate(z(nz))
    allocate(z_reverse(nz))
    do i=1, nz
        z(i) = i
        z_reverse(i) = nz - i + 1
    enddo

    Time = set_date(2,1,1,0,0,0)
    Time_step = set_time (3600,0)

    id_z1 =  diag_axis_init('zaxis1',  z,  'z', 'z', long_name='Z1')
    id_z2 =  diag_axis_init('zaxis2',  z,  'z', 'z', long_name='Z2')
    id_z_reverse =  diag_axis_init('zaxis3',  z_reverse,  'z_reverse', 'z', long_name='Z3')
    id_var1 = register_diag_field  ('atmos', 'ua_1', (/id_z1/), Time)
    id_var2 = register_diag_field  ('atmos', 'ua_2', (/id_z2/), Time)
    id_var3 = register_diag_field  ('atmos', 'ua_3', (/id_z_reverse/), Time)

    call diag_manager_set_time_end(set_date(2,1,2,0,0,0))
    do i = 1, 24
        Time = Time + Time_step
        used = send_data(id_var1, z, Time)
        used = send_data(id_var2, z, Time)
        used = send_data(id_var3, z, Time)
        call diag_send_complete(Time_step)
    end do

    call diag_manager_end(Time)

    call check_output()
    call fms_end

    contains

    subroutine check_output()
        type(FmsNetcdfFile_t) :: fileobj
        integer :: EXPECTED_NTIMES = 24
        integer :: SUB1_SIZE = 3
        integer :: SUB2_SIZE = 1
        character(len=20) :: EXPECTED_DIM_NAMES(2)
        real :: EXPECTED_ZSUBAXIS_1(3)
        real :: EXPECTED_ZSUBAXIS_2(1)
        real :: EXPECTED_ZSUBAXIS_3(3)

        EXPECTED_DIM_NAMES(2) = "time"
        if (.not. open_file(fileobj, "test_multiple_zbounds.nc", "read")) then
            call mpp_error(FATAL, "Unable to open the expected output file: test_var_masks.nc")
        endif

        call check_dimension(fileobj, "time", EXPECTED_NTIMES)

        EXPECTED_ZSUBAXIS_1 = (/3., 4., 5./)
        call check_dimension(fileobj, "zaxis1_sub01", SUB1_SIZE, EXPECTED_ZSUBAXIS_1)

        EXPECTED_ZSUBAXIS_2 = (/1./)
        call check_dimension(fileobj, "zaxis2_sub02", SUB2_SIZE, EXPECTED_ZSUBAXIS_2)

        EXPECTED_ZSUBAXIS_3 = (/5., 4., 3./)
        call check_dimension(fileobj, "zaxis3_sub03", SUB1_SIZE, EXPECTED_ZSUBAXIS_3)

        EXPECTED_DIM_NAMES(1) = "zaxis1_sub01"
        call check_variable(fileobj, "ua_1", EXPECTED_DIM_NAMES)

        EXPECTED_DIM_NAMES(1) = "zaxis2_sub02"
        call check_variable(fileobj, "ua_2", EXPECTED_DIM_NAMES)
        call close_file(fileobj)

    end subroutine check_output

    subroutine check_variable(fileobj, variable_name, expected_dimnames)
        type(FmsNetcdfFile_t), intent(in) :: fileobj
        character(len=*),      intent(in) :: variable_name
        character(len=*),      intent(in) :: expected_dimnames(:)

        character(len=20) :: dim_names(2)
        integer :: i

        call get_variable_dimension_names(fileobj, variable_name, dim_names)
        do i = 1, size(dim_names)
            if (trim(dim_names(i)) .ne. trim(expected_dimnames(i))) then
                print *, trim(dim_names(i)), " vs ", trim(expected_dimnames(i))
                call mpp_error(FATAL, "The dimension names are correct for "//trim(variable_name))
            endif
        enddo

    end subroutine check_variable

      !> @brief Check dimension data
  subroutine check_data(err_msg, actual_data, expected_data)
    character(len=*), intent(in) :: err_msg          !< Error message to append
    real,             intent(in) :: actual_data(:)   !< Dimension data from file
    real,             intent(in) :: expected_data(:) !< Expected data

    integer :: i

    do i = 1, size(actual_data)
      if (actual_data(i) .ne. expected_data(i)) &
        call mpp_error(FATAL, "The data is not expected for "//trim(err_msg))
    enddo
  end subroutine check_data

    subroutine check_dimension(fileobj, dimension_name, expected_size, expected_data)
        type(FmsNetcdfFile_t), intent(in) :: fileobj
        character(len=*),      intent(in) :: dimension_name
        integer,               intent(in) :: expected_size
        real, optional,        intent(in) :: expected_data(:)

        integer :: dim_size
        real, allocatable :: z_data(:)

        call get_dimension_size(fileobj, dimension_name, dim_size)
        if (dim_size .ne. expected_size) then
            call mpp_error(FATAL, trim(dimension_name)//" is not the expected size!")
        endif

        if (present(expected_data)) then
            allocate(z_data(dim_size))
            call read_data(fileobj, dimension_name, z_data)
            call check_data(dimension_name, z_data, expected_data)
        endif


    end subroutine
end program test_multiple_zbounds