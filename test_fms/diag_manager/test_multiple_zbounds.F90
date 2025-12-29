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
    integer                            :: id_var1          !< var id for the first variable
    integer                            :: id_var2          !< var_id for the second variable
    real,                  allocatable :: z(:)             !< z axis data
    integer                            :: nz               !< Size of the z dimension
    integer                            :: i
    logical                            :: used             !< Dummy argument to send_data

    call fms_init
    call set_calendar_type(JULIAN)
    call diag_manager_init

    nz = 10
    allocate(z(nz))
    do i=1, nz
        z(i) = i
    enddo

    Time = set_date(2,1,1,0,0,0)
    Time_step = set_time (3600,0)

    id_z1 =  diag_axis_init('zaxis1',  z,  'z', 'z', long_name='Z1')
    id_z2 =  diag_axis_init('zaxis2',  z,  'z', 'z', long_name='Z2')
    id_var1 = register_diag_field  ('atmos', 'ua_1', (/id_z1/), Time)
    id_var2 = register_diag_field  ('atmos', 'ua_2', (/id_z2/), Time)

    call diag_manager_set_time_end(set_date(2,1,2,0,0,0))
    do i = 1, 24
        Time = Time + Time_step
        z = real(i)
        used = send_data(id_var1, z, Time)
        used = send_data(id_var2, z, Time)

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

        EXPECTED_DIM_NAMES(2) = "time"
        if (.not. open_file(fileobj, "test_multiple_zbounds.nc", "read")) then
            call mpp_error(FATAL, "Unable to open the expected output file: test_var_masks.nc")
        endif

        call check_dimension(fileobj, "time", EXPECTED_NTIMES)
        call check_dimension(fileobj, "zaxis1_sub01", SUB1_SIZE)
        call check_dimension(fileobj, "zaxis2_sub02", SUB2_SIZE)

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

    subroutine check_dimension(fileobj, dimension_name, expected_size)
        type(FmsNetcdfFile_t), intent(in) :: fileobj
        character(len=*),      intent(in) :: dimension_name
        integer,               intent(in) :: expected_size

        integer :: dim_size

        call get_dimension_size(fileobj, dimension_name, dim_size)
        if (dim_size .ne. expected_size) then
            call mpp_error(FATAL, trim(dimension_name)//" is not the expected size!")
        endif

    end subroutine
end program test_multiple_zbounds