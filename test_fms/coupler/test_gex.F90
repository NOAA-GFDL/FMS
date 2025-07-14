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

!> This program tests the public API of gex_mod.
program test_coupler_gex

use fms_mod, only: fms_init, fms_end, check_nml_error
use mpp_mod, only: mpp_error, FATAL, input_nml_file
use fms_string_utils_mod, only: string
use tracer_manager_mod,  only: NO_TRACER
use field_manager_mod,   only: MODEL_LAND, MODEL_ATMOS, MODEL_OCEAN, NUM_MODELS
use gex_mod

implicit none

interface assert_eq
  procedure assert_int_eq
  procedure assert_str_eq
end interface assert_eq

! String indicating which test to run
character(50) :: test_name
namelist /test_gex_nml/ test_name

call fms_init
call gex_init
call read_test_nml

select case (trim(test_name))
  case ("get_n_ex_invalid_model_src")
    call get_n_ex_invalid_model_src
  case ("get_n_ex_invalid_model_rec")
    call get_n_ex_invalid_model_rec
  case ("get_property_invalid_tracer")
    call get_property_invalid_tracer
  case ("get_property_invalid_property")
    call get_property_invalid_property
  case default
    call atm_to_land
end select

call fms_end

contains

subroutine atm_to_land
  integer :: i, n
  character(:), allocatable :: name, units

  ! Number of atmosphere-to-land tracers (there should only be one)
  n = gex_get_n_ex(MODEL_ATMOS, MODEL_LAND)
  call assert_eq(n, 1, "gex_get_n_ex returned " // string(n) // &
                           " rather than the expected value of 1")

  ! Attempt to lookup a nonexistent tracer
  i = gex_get_index(MODEL_ATMOS, MODEL_LAND, "does_not_exist")
  call assert_eq(i, NO_TRACER, "gex_get_index returned " // string(i) // &
                           " rather than the expected value of " // string(NO_TRACER))

  ! Look up a tracer defined in the field table
  i = gex_get_index(MODEL_ATMOS, MODEL_LAND, "dryoa")

  name = gex_get_property(MODEL_ATMOS, MODEL_LAND, i, gex_name)
  call assert_eq(name, "dryoa", "gex_get_property returned tracer name '" // name // &
                                "' rather than the expected tracer name 'dryoa'")

  units = gex_get_property(MODEL_ATMOS, MODEL_LAND, i, gex_units)
  call assert_eq(units, "kg/m2/s", "gex_get_property returned units '" // units // &
                                  "' rather than the expected units 'kg/m2/s'")
end subroutine atm_to_land

subroutine get_n_ex_invalid_model_src
  integer :: n

  print "(A)", "Testing gex_get_n_ex with an invalid MODEL_SRC index"
  n = gex_get_n_ex(NUM_MODELS+1, 1)
  print "(A)", "Result: " // string(n)
end subroutine get_n_ex_invalid_model_src

subroutine get_n_ex_invalid_model_rec
  integer :: n

  print "(A)", "Testing gex_get_n_ex with an invalid MODEL_REC index"
  n = gex_get_n_ex(1, -1)
  print "(A)", "Result: " // string(n)
end subroutine get_n_ex_invalid_model_rec

subroutine get_property_invalid_tracer
  integer :: n
  character(:), allocatable :: property

  print "(A)", "Testing gex_get_property with an invalid tracer index"
  n = gex_get_n_ex(MODEL_ATMOS, MODEL_LAND)
  property = gex_get_property(MODEL_ATMOS, MODEL_LAND, n + 1, gex_name)
  print "(A)", "Result: " // property
end subroutine get_property_invalid_tracer

subroutine get_property_invalid_property
  integer :: i
  character(:), allocatable :: property

  print "(A)", "Testing gex_get_property with an invalid property index"
  i = gex_get_index(MODEL_ATMOS, MODEL_LAND, "dryoa")
  property = gex_get_property(MODEL_ATMOS, MODEL_LAND, i, 3)
  print "(A)", "Result: " // property
end subroutine get_property_invalid_property

subroutine read_test_nml
  integer :: ierr

  read (input_nml_file, nml=test_gex_nml, iostat=ierr)
  ierr = check_nml_error(ierr, 'test_gex_nml')
end subroutine

subroutine assert_int_eq(i1, i2, msg)
  integer, intent(in) :: i1, i2
  character(*), intent(in) :: msg

  if (i1.ne.i2) then
    call mpp_error(FATAL, msg)
  endif
end subroutine assert_int_eq

subroutine assert_str_eq(s1, s2, msg)
  character(*), intent(in) :: s1, s2
  character(*), intent(in) :: msg

  if (s1.ne.s2) then
    call mpp_error(FATAL, msg)
  endif
end subroutine assert_str_eq

end program test_coupler_gex
