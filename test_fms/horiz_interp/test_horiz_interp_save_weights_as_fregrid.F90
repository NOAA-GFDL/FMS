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
program main

  use constants_mod, only: DEG_TO_RAD
  use fms_mod, only: fms_init, fms_end
  use horiz_interp_mod, only: horiz_interp_type, horiz_interp_new
  use mpp_mod, only: FATAL, mpp_error
  use platform_mod, only: r8_kind

  implicit none

  integer, parameter :: nlon = 12
  integer, parameter :: nlat = 8

  type(horiz_interp_type) :: interp
  real(r8_kind),allocatable :: lon(:,:), lat(:,:), answer_area(:,:)
  integer :: i, j

  allocate(lon(nlon+1, nlat+1))
  allocate(lat(nlon+1, nlat+1))
  allocate(answer_area(nlon, nlat))

  do j=1, nlat+1
    do i=1, nlon+1
        lon(i,j) = real(i, r8_kind)*DEG_TO_RAD
        lat(i,j) = real(j, r8_kind)*DEG_TO_RAD
    end do
  end do

  call fms_init()

  call get_grid_area(nlon, nlat,lon, lat, answer_area)

  call horiz_interp_new(interp, lon, lat, lon, lat, interp_method="conservative", &
    is_latlon_in=.false., is_latlon_out=.false., save_weights_as_fregrid=.true.)

  if(any(interp%xgrid_area /= pack(answer_area, .true.))) then
    call mpp_error(FATAL, "saved xgrid_area does not match answers")
  end if

  call fms_end()

end program main
