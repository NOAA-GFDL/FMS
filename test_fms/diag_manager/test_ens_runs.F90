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

!> @brief  This programs tests diag manager when the file frequency is set to 0 days
program test_ens_runs

  use fms_mod,          only: fms_init, fms_end, string
  use diag_manager_mod, only: diag_axis_init, send_data, diag_send_complete, diag_manager_set_time_end, &
                              register_diag_field, diag_manager_init, diag_manager_end, register_static_field, &
                              diag_axis_init
  use time_manager_mod, only: time_type, operator(+), JULIAN, set_time, set_calendar_type, set_date
  use mpp_mod,          only: FATAL, mpp_error, mpp_npes, mpp_pe, mpp_get_current_pelist, mpp_set_current_pelist, &
                              input_nml_file
  use fms2_io_mod,      only: FmsNetcdfFile_t, open_file, close_file, read_data, get_dimension_size, &
                              set_filename_appendix, get_instance_filename
  use ensemble_manager_mod, only: get_ensemble_size, ensemble_manager_init

  implicit none

  integer         :: id_var0          !< diag field ids
  integer         :: id_axis1         !< Id for axis
  logical         :: used             !< for send_data calls
  integer         :: ntimes = 48      !< Number of time steps
  real            :: vdata            !< Buffer to store the data
  type(time_type) :: Time             !< "Model" time
  type(time_type) :: Time_step        !< Time step for the "simulation"
  integer         :: i                !< For do loops
  integer         :: npes             !< Number of pes in the current pelist
  integer, allocatable :: pelist(:)   !< Full pelist
  integer         :: ensemble_id      !< The ensemble id
  integer         :: ens_siz(6)       !< The size of the ensemble
  character(len=10) :: text           !< The filename appendix
  integer :: expected_ntimes

  integer                            :: io_status       !< Status after reading the namelist

  !< Configuration parameters
  logical :: using_ens_yaml = .true.      !< Indicates whether or not ensembles yamls are used in the test

  namelist / test_ens_runs_nml / using_ens_yaml

  call fms_init

  read (input_nml_file, test_ens_runs_nml, iostat=io_status)
  if (io_status > 0) call mpp_error(FATAL,'=>test_ens_runs: Error reading input.nml')

  call ensemble_manager_init
  npes = mpp_npes()
  if (npes .ne. 2) &
    call mpp_error(FATAL, "This test requires two pes to run")

  allocate(pelist(npes))
  call mpp_get_current_pelist(pelist)

  ens_siz = get_ensemble_size()
  if (ens_siz(1) .ne. 2) &
    call mpp_error(FATAL, "This test requires 2 ensembles")

  if (mpp_pe() < 1) then
    !< PE 0 is the first ensemble
    ensemble_id = 1
    call mpp_set_current_pelist((/0/))
    expected_ntimes = 48
  else
    ensemble_id = 2
    call mpp_set_current_pelist((/1/))
    expected_ntimes = 24
    if (.not. using_ens_yaml) expected_ntimes = 48
  endif

  write( text,'(a,i2.2)' ) 'ens_', ensemble_id
  call set_filename_appendix(trim(text))

  call set_calendar_type(JULIAN)
  call diag_manager_init

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time (3600,0) !< 1 hour
  call diag_manager_set_time_end(set_date(2,1,3,0,0,0))

  id_var0 = register_diag_field  ('ocn_mod', 'var0', Time)

  do i = 1, ntimes
    Time = Time + Time_step
    vdata = real(i)

    used = send_data(id_var0, vdata, Time)
    call diag_send_complete(Time_step)
  enddo

  call diag_manager_end(Time)

  call check_output()
  call fms_end

  contains

  !< @brief Check the diag manager output
  subroutine check_output()
    type(FmsNetcdfFile_t) :: fileobj     !< Fms2io fileobj
    integer               :: var_size    !< Size of the variable reading
    real, allocatable     :: var_data(:) !< Buffer to read variable data to
    integer               :: j           !< For looping
    character(len=255)    :: filename    !< Name of the diag file
    integer               :: scale_fac   !< Scale factor to use when determining the expected answer

    call get_instance_filename("test_ens.nc", filename)
    if (.not. open_file(fileobj, filename, "read")) &
      call mpp_error(FATAL, "Error opening file:"//trim(filename)//" to read")

    call get_dimension_size(fileobj, "time", var_size)
    if (var_size .ne. expected_ntimes) call mpp_error(FATAL, "The dimension of time in the file:"//&
                                                             "test_ens is not the correct size!")
    allocate(var_data(var_size))
    var_data = -999.99

    scale_fac = 1
    if (using_ens_yaml) scale_fac = ensemble_id

    call read_data(fileobj, "var0", var_data)
    do j = 1, var_size
      if (var_data(j) .ne. real(j * scale_fac))&
        call mpp_error(FATAL, "The variable data for var1 at time level:"//&
                              string(j)//" is not the correct value!")
    enddo

    call close_file(fileobj)
  end subroutine check_output
end program test_ens_runs
