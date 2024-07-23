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
!> @file
!! @brief unit tests for atmos_ocean_fluxes_mod
!! @author MiKyung.Lee
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program tests the two main subroutines in atmos_ocean_fluxes.
program test_atmos_ocean_fluxes

  use fms_mod, only: fms_init, fms_end
  use coupler_types_mod, only: coupler_1d_bc_type
  use field_manager_mod, only: fm_exists, fm_get_value
  use fm_util_mod, only: fm_util_get_real_array
  use mpp_mod, only: mpp_error, FATAL
  use platform_mod, only: r4_kind, r8_kind

  use atmos_ocean_fluxes_mod

  implicit none

  type(coupler_1d_bc_type) :: gas_fluxes
  type(coupler_1d_bc_type) :: gas_fields_atm
  type(coupler_1d_bc_type) :: gas_fields_ice

  integer, parameter :: num_bcs=4  !< number of fields
  integer, dimension(num_bcs) :: coupler_index !< coupler_index; index returned after calling aof_set_coupler_flux
  integer, dimension(num_bcs) :: atm_tr_index  !< tracer_index
  integer, dimension(num_bcs) :: nparameter    !< number of parameters for a given flux_type/implementation
  character(8), dimension(num_bcs) :: flux_name !< flux name
  character(24), dimension(num_bcs) :: flux_type !< flux_type for the given flux name
  character(6), dimension(num_bcs) :: impl !< implementation type for a given flux_type
  real(FMS_CP_TEST_KIND_), dimension(num_bcs) :: mol_wt      !< value of mol_wt
  real(FMS_CP_TEST_KIND_), dimension(num_bcs) :: param_array !< parameter array

  !> The flux names are made up.
  flux_name=["vampires", &
             "eric    ", &
             "grapes  ", &
             "coffee  "]
  !> Flux types
  flux_type=["air_sea_gas_flux_generic", &
             "air_sea_gas_flux        ", &
             "air_sea_deposition      ", &
             "land_sea_runoff         "]
  !> Implmentation type for the corresponding flux type
  impl=["ocmip2", &
        "linear", &
        "dry   ", &
        "river "]
  !> made up parameters
  param_array =[1.0, 2.0, 3.0, 4.0]
  !> made up atm_tr_index
  atm_tr_index=[1, 2, 3, 4]
  !> made up mol_wt
  mol_wt      =[1.0, 2.0, 3.0, 4.0]
  !> number of parameters associated with the corresponding implementation
  nparameter  =[2, 3, 1, 1]

  call fms_init
  !> initialize atmos_ocean_type_fluxes to set the /coupler_mod/types field tree
  call atmos_ocean_type_fluxes_init
  !> setting fake flux fields
  call test_aof_set_coupler_flux
  !> initialize gas_fluxes, gas_fields_atm, gas_fields_ice
  call test_atmos_ocean_fluxes_init
  !> checking gas_fluxes, gas_fields_atm, and gas_fields_ice have been initialized correctly
  call test_coupler_1d_bc_type
  call fms_end

contains
  !--------------------------------------
  subroutine test_atmos_ocean_fluxes_init

    !> This subroutine calls atmos_ocean_fluxes to initialize
    !! gas_fluxes, gas_fields_atm, gas_fields_ice

    implicit none

    write(*,*) "*** TEST_ATMOS_OCEAN_FLUXES_INIT ***"
    call atmos_ocean_fluxes_init(gas_fluxes, gas_fields_atm, gas_fields_ice, &
                                 use_r4_kind= FMS_CP_TEST_KIND_ .eq. r4_kind)

  end subroutine test_atmos_ocean_fluxes_init
  !--------------------------------------
  subroutine test_aof_set_coupler_flux

    !> This subroutine sets up the flux fields via aof_set_coupler_flux

    implicit none

    character(100) :: cresults, thelist
    real(FMS_CP_TEST_KIND_) :: rresults, rresults2(num_bcs)
    integer :: i, success, n

    write(*,*) "*** TEST_AOF_SET_COUPLER_FLUX ***"

    !> set flux fields
    do i=1, num_bcs
       coupler_index(i)=aof_set_coupler_flux(name=flux_name(i), &
                                             flux_type=flux_type(i), &
                                             implementation=impl(i), &
                                             atm_tr_index=atm_tr_index(i), &
                                             param=real(param_array(1:nparameter(i)), r8_kind), &
                                             mol_wt=real(mol_wt(i), r8_kind))
    end do

    !> check answers
    do i=1, num_bcs
       !> check to see that the flux field has been created
       thelist="coupler_mod/fluxes/"//trim(flux_name(i))
       if(.not. fm_exists( trim(thelist)) ) call mpp_error(FATAL, 'coupler')

       !> check to see that the flux_type value has been set correctly
       success=fm_get_value( trim(thelist)//"/flux_type", cresults )
       call check_answers(flux_type(i), cresults, 'test_aof_set_coupler_flux with flux_type')

       !> check to see the implementation value has been set correctly
       success=fm_get_value( trim(thelist)//"/implementation", cresults )
       call check_answers(impl(i), cresults, 'test_aof_set_coupler_flux with implementation')

       !> check to see the param array has been set correctly
       rresults2(1:nparameter(i))=fm_util_get_real_array( trim(thelist)//"/param" )
       do n=1, nparameter(i)
          call check_answers(param_array(n), rresults2(n), 'test_aof_set_coupler_fluxes ERROR with param')
       end do

       !> check to see that the mol_wt value has been set correctly
       success=fm_get_value( trim(thelist)//'/mol_wt', rresults )
       call check_answers(mol_wt(i), rresults, 'test_aof_set_coupler_fluxes ERROR with mol_wt')
    end do

  end subroutine test_aof_set_coupler_flux
  !--------------------------------------
  subroutine test_coupler_1d_bc_type

    !> This subroutine checks that the three coupler_1d_bc_types have been
    !! initialized correctly in atmos_ocean_fluxes_init

    implicit none

    integer :: i, n

    write(*,*) "*** TEST_COUPLER_1D_BC_TYPE ***"

    if( gas_fluxes%num_bcs.ne.num_bcs ) call mpp_error(FATAL, 'error')

    do i=1, num_bcs
       !> check fluxes name
       call check_answers(flux_name(i), gas_fluxes%FMS_TEST_BC_TYPE_(i)%name, 'gas_fluxes flux name')
       call check_answers(flux_name(i), gas_fields_atm%FMS_TEST_BC_TYPE_(i)%name, 'gas_fields_atms flux name')
       call check_answers(flux_name(i), gas_fields_ice%FMS_TEST_BC_TYPE_(i)%name, 'gas_fields_ice flux name')

       !> check implementation
       call check_answers(impl(i), gas_fluxes%FMS_TEST_BC_TYPE_(i)%implementation, 'gas_fluxes impl')
       call check_answers(impl(i), gas_fields_atm%FMS_TEST_BC_TYPE_(i)%implementation, 'gas_fields_atm impl')
       call check_answers(impl(i), gas_fields_ice%FMS_TEST_BC_TYPE_(i)%implementation, 'gas_fields_ice impl')

       !> check param
       do n=1, nparameter(i)
          call check_answers(param_array(n), gas_fluxes%FMS_TEST_BC_TYPE_(i)%param(n), 'gas_fluxes param')
       end do

       !> check mol_wt
       call check_answers(mol_wt(i), gas_fluxes%FMS_TEST_BC_TYPE_(i)%mol_wt, 'gas_fluxes mol_wt')
       call check_answers(mol_wt(i), gas_fields_atm%FMS_TEST_BC_TYPE_(i)%mol_wt, 'gas_fields_atm mol_wt')
       call check_answers(mol_wt(i), gas_fields_ice%FMS_TEST_BC_TYPE_(i)%mol_wt, 'gas_fields_ice mol_wt')

    end do

  end subroutine test_coupler_1d_bc_type
  !--------------------------------------
  subroutine check_answers(answers, myresults, whoami)

    implicit none
    class(*), intent(in) :: answers, myresults
    character(*) :: whoami

    select type(answers) ; type is(character(*))
       select type(myresults) ; type is(character(*))
          if(trim(answers).ne.trim(myresults)) then
             write(*,*) 'EXPECTED '//trim(answers)//' but got '//trim(myresults)
             call mpp_error(FATAL, trim(whoami))
          end if
       end select
    type is(real(r4_kind))
       select type(myresults) ; type is(real(r4_kind))
          if(answers.ne.myresults) then
             write(*,*) 'EXPECTED ', answers, ' but got ', myresults
             call mpp_error(FATAL, trim(whoami))
          end if
       end select
    type is(real(r8_kind))
       select type(myresults) ; type is(real(r8_kind))
          if(answers.ne.myresults) then
             write(*,*) 'EXPECTED ', answers, ' but got ', myresults
             call mpp_error(FATAL, trim(whoami))
          end if
       end select
    end select

  end subroutine check_answers
  !--------------------------------------
end program test_atmos_ocean_fluxes
