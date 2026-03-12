! -*-f90-*-

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
!> @ingroup mpp_data_mod
!> @brief Holds dummy constants and stack data for @ref mpp_mod and @ref mpp_domains_mod.
!! Accessible through @ref mpp_data_mod

!----------------------------------------------------------------!
! The following data is used in mpp_mod and its components       !
!----------------------------------------------------------------!
real(r8_kind), allocatable :: mpp_stack(:)

!--- some dummy variables with dummy values that will never be used
integer, parameter :: stat=-999
integer, parameter :: ptr_stack = -999
integer, parameter :: status=-999, ptr_status = -999
integer, parameter :: remote_data_loc=-999, ptr_remote = -999
integer, parameter :: sync=-999, ptr_sync = -999
integer, parameter :: mpp_from_pe = -999, ptr_from = -999

!-------------------------------------------------------------------!
! The following data is used in mpp_domains_mod and its components  !
!-------------------------------------------------------------------!
real(r8_kind), allocatable :: mpp_domains_stack(:)
real(r8_kind), allocatable :: mpp_domains_stack_nonblock(:)
!--- some dummy variables with dummy values that will never be used
integer, parameter :: ptr_domains_stack = -999
integer, parameter :: ptr_domains_stack_nonblock = -999
