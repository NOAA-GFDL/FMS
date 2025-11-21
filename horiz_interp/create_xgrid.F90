module create_xgrid_mod

  use, intrinsic :: iso_c_binding

  implicit none

  interface 

    function create_xgrid_1dx2d_order1(nlon_in, nlat_in, nlon_out, nlat_out, &
                                       lon_in, lat_in, lon_out, lat_out, mask_in, &
                                       i_in, j_in, i_out, j_out, xgrid_area) BIND(C, NAME="create_xgrid_1dx2d_order1")
      import :: C_INT, C_DOUBLE
      integer(C_INT), value       :: nlon_in
      integer(C_INT), value       :: nlat_in
      integer(C_INT), value       :: nlon_out
      integer(C_INT), value       :: nlat_out
      real(C_DOUBLE), dimension(*) :: lon_in
      real(C_DOUBLE), dimension(*) :: lat_in
      real(C_DOUBLE), dimension(*) :: lon_out
      real(C_DOUBLE), dimension(*) :: lat_out
      real(C_DOUBLE), dimension(*) :: mask_in
      integer(C_INT), dimension(*) :: i_in
      integer(C_INT), dimension(*) :: j_in
      integer(C_INT), dimension(*) :: i_out
      integer(C_INT), dimension(*) :: j_out
      real(C_DOUBLE), dimension(*) :: xgrid_area
    end function

    function create_xgrid_2dx1d_order1(nlon_in, nlat_in, nlon_out, nlat_out, &
                                       lon_in, lat_in, lon_out, lat_out, mask_in, &
                                       i_in, j_in, i_out, j_out, xgrid_area) BIND(C, NAME="create_xgrid_2dx1d_order1")
      import :: C_INT, C_DOUBLE
      integer(C_INT), value       :: nlon_in
      integer(C_INT), value       :: nlat_in
      integer(C_INT), value       :: nlon_out
      integer(C_INT), value       :: nlat_out
      real(C_DOUBLE), dimension(*) :: lon_in
      real(C_DOUBLE), dimension(*) :: lat_in
      real(C_DOUBLE), dimension(*) :: lon_out
      real(C_DOUBLE), dimension(*) :: lat_out
      real(C_DOUBLE), dimension(*) :: mask_in
      integer(C_INT), dimension(*) :: i_in
      integer(C_INT), dimension(*) :: j_in
      integer(C_INT), dimension(*) :: i_out
      integer(C_INT), dimension(*) :: j_out
      real(C_DOUBLE), dimension(*) :: xgrid_area
    end function

    function create_xgrid_2dx2d_order1(nlon_in, nlat_in, nlon_out, nlat_out, &
                                       lon_in, lat_in, lon_out, lat_out, mask_in, &
                                       i_in, j_in, i_out, j_out, xgrid_area) BIND(C, NAME="create_xgrid_2dx2d_order1")
      import :: C_INT, C_DOUBLE
      integer(C_INT), value       :: nlon_in
      integer(C_INT), value       :: nlat_in
      integer(C_INT), value       :: nlon_out
      integer(C_INT), value       :: nlat_out
      real(C_DOUBLE), dimension(*) :: lon_in
      real(C_DOUBLE), dimension(*) :: lat_in
      real(C_DOUBLE), dimension(*) :: lon_out
      real(C_DOUBLE), dimension(*) :: lat_out
      real(C_DOUBLE), dimension(*) :: mask_in
      integer(C_INT), dimension(*) :: i_in
      integer(C_INT), dimension(*) :: j_in
      integer(C_INT), dimension(*) :: i_out
      integer(C_INT), dimension(*) :: j_out
      real(C_DOUBLE), dimension(*) :: xgrid_area
    end function


    function get_grid_area(nlon, nlat, lon, lat, area) BIND(C, NAME="get_grid_area")
      import :: C_INT, C_DOUBLE
      integer(C_INT), value :: nlon
      integer(C_INT), value :: nlat
      real(C_DOUBLE), dimension(*) :: lon
      real(C_DOUBLE), dimension(*) :: lat
      real(C_DOUBLE), dimension(*) :: area
    end function

  end interface

! the routines used in fms:
!../../libFMS/.libs/libFMS.so: undefined reference to `get_grid_area_'
!../../libFMS/.libs/libFMS.so: undefined reference to `create_xgrid_great_circle_'
!../../libFMS/.libs/libFMS.so: undefined reference to `create_xgrid_2dx2d_order1_'
!../../libFMS/.libs/libFMS.so: undefined reference to `create_xgrid_2dx1d_order1_'
!../../libFMS/.libs/libFMS.so: undefined reference to `get_grid_great_circle_area_'
!../../libFMS/.libs/libFMS.so: undefined reference to `get_maxxgrid_'
!../../libFMS/.libs/libFMS.so: undefined reference to `create_xgrid_1dx2d_order1_'


end module
