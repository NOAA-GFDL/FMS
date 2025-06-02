!***********************************************************************
!                   GNU Lesser General Public License
!
! This file is part of the GFDL Flexible Modeling System (FMS).
!
! FMS is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! FMS is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!**********************************************************************/
!> This program ensures all necessary functions exist in grid_utils for
!! creating the exchange grid for second order remapping.
!! The first order exchange grid creation is tested in test_horiz_interp.
!! This test is rudimentary and only checks that create_xgrid_order2 returns
!! without failure.

program test_create_xgrid_order2

  use horiz_interp_mod
  use constants_mod, only: DEG_TO_RAD
  implicit none

  integer, parameter :: lkind = HI_TEST_KIND_

  integer, parameter :: nlon_in = 10   !< number of input grid cells in lon direction
  integer, parameter :: nlat_in = 10   !< number of input grid cells in the lat direction
  integer, parameter :: nlon_out = nlon_in * 2 !< number of output grid cells in lon direction
  integer, parameter :: nlat_out = nlat_in * 2 !< number of output grid cells in lat direction
  integer, parameter :: ngridpts_in = (nlon_in+1)*(nlat_in+1)    !< number of input gridpoints
  integer, parameter :: ngridpts_out = (nlon_out+1)*(nlat_out+1) !< number of output gridpoints
  integer, parameter :: nxgrid = nlon_out *nlat_out !< expected number of exchange grid cells

  real(HI_TEST_KIND_) :: lon_in(ngridpts_in)   !< longitudinal values of input grid cell vertices
  real(HI_TEST_KIND_) :: lat_in(ngridpts_in)   !< latitudinal values of input grid cell vertices
  real(HI_TEST_KIND_) :: lon_out(ngridpts_out) !< longitudinal values of output grid cell vertices
  real(HI_TEST_KIND_) :: lat_out(ngridpts_out) !< latitudinal values of output grid cell vertices
  real(HI_TEST_KIND_) :: mask(nlon_in*nlat_in) !< mask to skip input grid cell

  integer :: i_in(nxgrid)  !< input parent cell indices
  integer :: j_in(nxgrid)  !< input parent cell indices
  integer :: i_out(nxgrid) !< output parent cell indices
  integer :: j_out(nxgrid) !< output parent cell indices
  real(HI_TEST_KIND_) :: xgrid_area(nxgrid) !< exchange grid cell areas
  real(HI_TEST_KIND_) :: xgrid_clon(nxgrid) !< longitudinal values of exchange grid cell centroid point
  real(HI_TEST_KIND_) :: xgrid_clat(nxgrid) !< latitudinal values of exchange grid cell centroid point

  mask = 1.0_lkind

  call get_grid(nlon_in, nlat_in, lon_in, lat_in)
  call get_grid(nlon_out, nlat_out, lon_out, lat_out)

  call test_create_xgrid_2dx2d_order2(nlon_in, nlat_in, nlon_out, nlat_out, nxgrid, &
                                      mask, lon_in, lat_in, lon_out, lat_out,       &
                                      i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat)

contains

  !> Returns lon and lat arrays holding grid point values
  subroutine get_grid(nlon, nlat, lon, lat)

    implicit none
    integer, intent(in) :: nlon, nlat                  !< number of cell in lon and lat direction
    real(HI_TEST_KIND_), intent(out) :: lon(:), lat(:) !< lon and lat values at cell vertices

    integer :: ilon, ilat, igridpt
    real :: dlat=0.0_lkind, dlon=0.0_lkind
    real :: lon_start=0.0_lkind, lat_start=-90.0_lkind*DEG_TO_RAD

    dlat = 180._lkind/real(nlat, HI_TEST_KIND_) * DEG_TO_RAD
    dlon = 360._lkind/real(nlon, HI_TEST_KIND_) * DEG_TO_RAD

    igridpt = 1
    do ilat=1, nlat+1
      do ilon=1, nlon+1
        lon(igridpt) = lon_start + real(ilon-1, HI_TEST_KIND_)*dlon
        lat(igridpt) = lat_start + real(ilat-1, HI_TEST_KIND_)*dlat
        igridpt = igridpt + 1
      end do
    end do

  end subroutine get_grid


  !> Calls create_xgrid_2dx2d_order2.  This subroutine also checks the returned value of nxgrid
  subroutine test_create_xgrid_2dx2d_order2(nlon_inl, nlat_inl, nlon_outl, nlat_outl, nxgridl, &
                                            maskl, lon_inl, lat_inl, lon_outl, lat_outl,       &
                                            i_inl, j_inl, i_outl, j_outl, xgrid_areal, xgrid_clonl, xgrid_clatl)

    implicit none
    integer, intent(in)    :: nlon_inl, nlat_inl, nlon_outl, nlat_outl, nxgridl !< number of grid cells
    integer, intent(inout) :: i_inl(:), j_inl(:), i_outl(:), j_outl(:)  !< parent cell indices
    real(HI_TEST_KIND_), intent(in) :: lon_inl(:), lat_inl(:), lon_outl(:), lat_outl(:) !< lon and lat
    real(HI_TEST_KIND_), intent(in) :: maskl(:) !< input grid cell mask
    real(HI_TEST_KIND_), intent(out) :: xgrid_areal(:), xgrid_clonl(:), xgrid_clatl(:) !< returned xgrid info

    integer :: create_xgrid_2dx2d_order2
    integer :: nxgrid_out

    nxgrid_out = create_xgrid_2dx2d_order2(nlon_inl, nlat_inl, nlon_outl, nlat_outl, lon_inl, lat_inl, &
                                           lon_outl, lat_outl, maskl, i_inl, j_inl, i_outl, j_outl, xgrid_areal, &
                                           xgrid_clonl, xgrid_clatl)

    if(nxgrid_out /= nxgridl) then
      write(*,*) 'Expected', nxgrid, 'but got', nxgrid_out
      error stop
    end if

  end subroutine test_create_xgrid_2dx2d_order2

end program test_create_xgrid_order2
