program test_dm_weights
  use fms_mod, only: fms_init, fms_end
  use diag_manager_mod, only: diag_manager_init, diag_axis_init, register_diag_field, send_data, &
                              diag_send_complete, diag_manager_set_time_end, diag_manager_end
  use mpp_mod, only: mpp_error, mpp_pe, mpp_root_pe, FATAL
  use time_manager_mod, only: time_type, set_calendar_type, JULIAN, set_time, set_date, operator(+)
  use fms2_io_mod

  implicit none

  type(time_type)                    :: Time             !< Time of the simulation
  type(time_type)                    :: Time_step        !< Time_step of the simulation
  integer                            :: nx               !< Number of x points
  integer                            :: ny               !< Number of y points
  integer                            :: id_x             !< Axis id for the x dimension
  integer                            :: id_y             !< Axis id for the y dimension
  integer                            :: id_var1          !< Field id for 1st variable
  integer                            :: id_var2          !< Field id for 2nd variable
  logical                            :: used             !< Dummy argument to send_data
  real,                  allocatable :: x(:)             !< X axis data
  real,                  allocatable :: y(:)             !< Y axis_data
  real,                  allocatable :: var1_data(:,:)   !< Data for variable 1
  integer                            :: i                !< For do loops
  integer                            :: ntimes           !< Number of times to run the simulation for

  call fms_init()
  call set_calendar_type(JULIAN)
  call diag_manager_init()

  nx = 10
  ny = 15
  ntimes = 6

  allocate(x(nx), y(ny))
  allocate(var1_data(nx,ny))
  do i=1,nx
    x(i) = i
  enddo
  do i=1,ny
    y(i) = -91 + i
  enddo

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time (3600,0) !< 1 hour

  id_x  = diag_axis_init('x',  x,  'point_E', 'x', long_name='point_E')
  id_y  = diag_axis_init('y',  y,  'point_N', 'y', long_name='point_N')

  id_var1 = register_diag_field  ('atmos', 'ua', (/id_x, id_y/), Time)

  call diag_manager_set_time_end(set_date(2,1,1,ntimes,0,0))
  do i = 1, ntimes
    Time = Time + Time_step
    var1_data = real(i)
    used = send_data(id_var1, var1_data, Time, weight=real(i/10.))
    call diag_send_complete(Time_step)
  enddo

  call diag_manager_end(Time)

  call check_answers()
  call fms_end()

  contains

  subroutine check_answers()
    type(FmsNetcdfFile_t) :: fileobj
    integer :: j
    real :: ans_var
    real :: vardata_out(nx, ny)

    if (.not. open_file(fileobj, "test_weights.nc", "read")) &
        call mpp_error(FATAL, "unable to open test_var_masks.nc for reading")

    ans_var = 0
    do j = 1, ntimes
      if (mod(j,2) .eq. 0) then
        print *, "Checking answers for time = ", j/2
        ans_var = (j*(j/10.)+(j-1.)*(j-1.)/10.)/(j/10. + (j-1)/10.)
        call read_data(fileobj, "ua", vardata_out, unlim_dim_level=j/2)
        if (any(abs(ans_var-vardata_out) > 0.0000001)) &
          call mpp_error(FATAL, "The answer is not the expected result!")
      endif
    enddo
    call close_file(fileobj)
  end subroutine check_answers
end program test_dm_weights