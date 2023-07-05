program test_orbital_parameters

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

    call test_set_orbital_parameters

    contains
    
    !Test set_orbital_parameters breaks when arguments are not in correct bounds (XFAIL)
    subroutine test_set_orbital_parameters
        real :: ecc_in = 0.01671d0, ecc_in_check
        real :: obliq_in = 23.439d0, obliq_in_check
        real :: per_in = 102.932d0, per_in_check
        
        ! Eccentricity of Earth's orbit not in between 0.0 and 0.99
        ecc_in_check = -1.0
        call set_orbital_parameters(ecc_in_check, obliq_in, per_in)
        ecc_in_check = 1.0 
        call set_orbital_parameters(ecc_in_check, obliq_in, per_in)

        ! Obliquity not in between -90 and 90
        obliq_in_check = -91
        call set_orbital_parameters(ecc_in, obliq_in_check, per_in)
        obliq_in_check = 91
        call set_orbital_parameters(ecc_in, obliq_in_check, per_in)

        ! Perihelion not in between 0.0 and 360.0
        per_in_check = -1.0
        call set_orbital_parameters(ecc_in, obliq_in, per_in_check)
        per_in_check = 361.0
        call set_orbital_parameters(ecc_in, obliq_in, per_in_check)


    end subroutine test_set_orbital_parameters


end program test_orbital_parameters