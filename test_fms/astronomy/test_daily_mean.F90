program test_daily_mean

  use fms_mod,           only: fms_init, fms_end
  use mpp_mod,           only : mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use astronomy_mod,     only: astronomy_init, set_orbital_parameters, get_orbital_parameters, &
                               diurnal_solar, daily_mean_solar, annual_mean_solar
  use time_manager_mod,  only: JULIAN, set_calendar_type
  use platform_mod,      only: r4_kind, r8_kind

  implicit none

 
  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init
  call fms_end

  call test_daily_mean_solar_2d

  contains
  
  subroutine test_daily_mean_solar_2d
    
    !print *, "NEW TEST "


  end subroutine test_daily_mean_solar_2d


end program test_daily_mean