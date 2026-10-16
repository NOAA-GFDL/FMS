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

!> @brief  This program tests diag_manager's wildcard (%) file name suffix, and confirms that
!! diag_manager_nml::wildcard_filename_prefix/wildcard_filename_separator control how the
!! substituted time fields are joined. It expects a diag_table.yaml with:
!! - file_name: test_wildcard%4yr%2mo%2dy%2hr
!!   freq: 6 hours
!!   new_file_freq: 12 hours
!!   file_duration: 18 hours
!!   base_date: 2 1 1 0 0 0
!! which, over a 48 hour run starting at 2-1-1 0:0:0, produces two files at hour 06 and hour 15.
!! The expected suffixes for those two files are passed in via the test's own namelist so that the
!! same program can be re-run with different wildcard_filename_prefix/separator settings.
program test_wildcard_filename

  use fms_mod,          only: fms_init, fms_end
  use diag_manager_mod, only: send_data, diag_send_complete, diag_manager_set_time_end, &
                              register_diag_field, diag_manager_init, diag_manager_end
  use time_manager_mod, only: time_type, operator(+), JULIAN, set_time, set_calendar_type, set_date
  use mpp_mod,          only: FATAL, mpp_error, input_nml_file
  use fms2_io_mod,      only: FmsNetcdfFile_t, open_file, close_file

  implicit none

  integer         :: id_var0                       !< diag field id
  logical         :: used                           !< for send_data calls
  integer         :: ntimes = 48                    !< Number of hourly time steps (2 days)
  type(time_type) :: Time                           !< "Model" time
  type(time_type) :: Time_step                      !< Time step for the "simulation"
  integer         :: i                              !< For do loops
  integer         :: io_status                      !< Status when reading the namelist

  !< Expected suffixes of the two files created by the 12 hour new_file_freq/18 hour file_duration
  !! diag_table.yaml described above, given the wildcard_filename_prefix/separator that was set in
  !! diag_manager_nml for this run. Defaults match the historical "_"-joined behavior.
  character(len=32) :: expected_suffix1 = '_0002_01_01_06'
  character(len=32) :: expected_suffix2 = '_0002_01_01_15'

  namelist / test_wildcard_filename_nml / expected_suffix1, expected_suffix2

  call fms_init

  read (input_nml_file, test_wildcard_filename_nml, iostat=io_status)
  if (io_status > 0) call mpp_error(FATAL, '=>test_wildcard_filename: Error reading input.nml')

  call set_calendar_type(JULIAN)
  call diag_manager_init()

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time(3600, 0) !< 1 hour
  call diag_manager_set_time_end(set_date(2,1,3,0,0,0))

  id_var0 = register_diag_field('ocn_mod', 'var0', Time)

  do i = 1, ntimes
    Time = Time + Time_step
    used = send_data(id_var0, real(i), Time)
    call diag_send_complete(Time_step)
  enddo

  call diag_manager_end(Time)

  call check_output()
  call fms_end

  contains

  !< @brief Confirm the two expected wildcard-named files were created
  subroutine check_output()
    type(FmsNetcdfFile_t) :: fileobj !< Fms2io fileobj

    if (.not. open_file(fileobj, "test_wildcard"//trim(expected_suffix1)//".nc", "read")) &
      call mpp_error(FATAL, "Error opening file: test_wildcard"//trim(expected_suffix1)//".nc")
    call close_file(fileobj)

    if (.not. open_file(fileobj, "test_wildcard"//trim(expected_suffix2)//".nc", "read")) &
      call mpp_error(FATAL, "Error opening file: test_wildcard"//trim(expected_suffix2)//".nc")
    call close_file(fileobj)
  end subroutine check_output

end program test_wildcard_filename
