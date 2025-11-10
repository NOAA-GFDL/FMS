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
!> @file
!! @brief unit test for set_tracer_profile function
!! @author MiKyung Lee
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program only tests set_tracer_profile for cases where
!! profile_type is fixed; and for cases where profile_type is profile
!! TODO:  Unit tests for the remaining subroutines and functions in tracer_manager_mod

program test_tracer_manager

  use fms_mod, only: fms_init, fms_end
  use mpp_mod, only: mpp_error, FATAL
  use field_manager_mod, only: field_manager_init, MODEL_ATMOS, MODEL_OCEAN, MODEL_LAND, &
                               fm_change_list, fm_get_value, fm_get_current_list
  use tracer_manager_mod
  use platform_mod, only: r4_kind, r8_kind

  implicit none

  call fms_init
  call test_set_tracer_profile
  call fms_end

contains

  subroutine test_set_tracer_profile
  integer, parameter :: numlevels=10
  integer, parameter :: npoints=5

  integer :: tracer_index, i, j, k
  logical :: success
  real(TEST_FMS_KIND_) :: top_value, bottom_value, surf_value, multiplier
  real(TEST_FMS_KIND_) :: tracer_out1(1,1,1), tracer_out2(npoints,npoints,numlevels)
  real(TEST_FMS_KIND_) :: answer1(1,1,1), answer2(npoints,npoints,numlevels)

  integer, parameter :: lkind=TEST_FMS_KIND_

  character(128) :: err_message

  call fms_init
  call tracer_manager_init

  !-- profile_type=fixed --!

  !> the tracer 'radon' profile type is 'fixed' (see field_table)
  !> the tracer field value should be zero.
  tracer_index=get_tracer_index(MODEL_ATMOS, 'radon')
  call set_tracer_profile(MODEL_ATMOS, tracer_index, tracer_out1,err_message)
  !> answer
  answer1(1,1,1)=0.0_lkind
  !> check results
  if(tracer_out1(1,1,1).ne.answer1(1,1,1)) call mpp_error(FATAL,'ATM tracer field value should be 0.0')

  !-- ATM profile_type=profile --!
  tracer_index=get_tracer_index(MODEL_ATMOS, 'immadeup')
  call set_tracer_profile(MODEL_ATMOS,tracer_index,tracer_out2,err_message)
  !> answer
  success=fm_get_value("/atmos_mod/tracer/immadeup/profile_type/profile/top_value", top_value)
  success=fm_get_value("/atmos_mod/tracer/immadeup/profile_type/profile/surface_value", surf_value)
  multiplier = exp( log (top_value/surf_value) /real(numlevels-1,lkind) )
  answer2(:,:,1)=surf_value
  do i=2,numlevels
     answer2(:,:,i) = answer2(:,:,i-1)*multiplier
  end do
  !> check results
  do k=1, numlevels
     do j=1, npoints
        do i=1, npoints
           if( tracer_out2(i,j,k) .ne. answer2(i,j,k)) &
                call mpp_error(FATAL, 'ATM tracer field value error for profile_type=profile')
        end do
     end do
  end do

  !-- OCEAN profile_type=profile --!
  tracer_index=get_tracer_index(MODEL_OCEAN, 'immadeup2')
  call set_tracer_profile(MODEL_OCEAN,tracer_index,tracer_out2,err_message)
  !> answer
  success=fm_get_value("/ocean_mod/tracer/immadeup2/profile_type/profile/bottom_value", bottom_value)
  success=fm_get_value("/ocean_mod/tracer/immadeup2/profile_type/profile/surface_value", surf_value)
  multiplier = exp( log (bottom_value/surf_value) /real(numlevels-1,lkind))
  answer2(:,:,numlevels)=surf_value
  do i=numlevels-1, 1, -1
     answer2(:,:,i) = answer2(:,:,i+1)*multiplier
  end do
  !> check results
  do k=1, numlevels
     do j=1, npoints
        do i=1, npoints
           if( tracer_out2(i,j,k) .ne. answer2(i,j,k)) &
                call mpp_error(FATAL, 'OCEAN tracer field value error for profile_type=profile')
        end do
     end do
  end do


end subroutine test_set_tracer_profile

end program test_tracer_manager
