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
!> @defgroup axis_utils2_mod axis_utils2_mod
!> @ingroup axis_utils
!> @brief A set of utilities for manipulating axes and extracting axis attributes.
!! FMS2_IO equivalent version of @ref axis_utils_mod.
!> @author M.J. Harrison

!> @addtogroup axis_utils2_mod
!> @{
module axis_utils2_mod
  use mpp_mod,      only: mpp_error, FATAL, stdout
  use fms_mod,      only: lowercase, uppercase, string_array_index, fms_error_handler
  use fms2_io_mod,  only: FmsNetcdfDomainFile_t, variable_att_exists, FmsNetcdfFile_t, &
                          get_variable_num_dimensions, get_variable_attribute,  &
                          get_variable_size, read_data, variable_exists
  use platform_mod, only: r4_kind, r8_kind

  implicit none

  public get_axis_cart, get_axis_modulo, lon_in_range, &
         tranlon, frac_index, nearest_index, interp_1d, get_axis_modulo_times, axis_edges

  private

  integer, parameter :: maxatts = 100
  real(r8_kind), parameter    :: epsln = 1.e-10_r8_kind
  real(r8_kind), parameter    :: fp5 = 0.5_r8_kind, f360 = 360.0_r8_kind

!> @}
! Include variable "version" to be written to log file.
#include<file_version.h>

  !> Perform 1D interpolation between grids.
  !!
  !> Data and grids can have 1, 2, or 3 dimensions.
  !! @param grid1 grid for data1
  !! @param grid2 grid for data2
  !! @param data1 Data to interpolate
  !! @param [inout] data2 Interpolated data
  !! @param method Either "linear" or "cubic_spline" interpolation method, default="linear"
  !! @ingroup axis_utils2_mod

  interface axis_edges
    module procedure axis_edges_r4, axis_edges_r8
  end interface axis_edges

  interface lon_in_range
    module procedure lon_in_range_r4, lon_in_range_r8
  end interface lon_in_range

  interface frac_index
    module procedure frac_index_r4, frac_index_r8
  end interface frac_index

  interface nearest_index
      module procedure nearest_index_r4, nearest_index_r8
  end interface nearest_index

  interface tranlon
      module procedure tranlon_r4, tranlon_r8
  end interface tranlon

  interface interp_1d_linear
      module procedure interp_1d_linear_r4, interp_1d_linear_r8
  end interface interp_1d_linear

  interface interp_1d_cubic_spline
        module procedure interp_1d_cubic_spline_r4, interp_1d_cubic_spline_r8
  end interface interp_1d_cubic_spline

  interface interp_1d
     module procedure interp_1d_1d_r4, interp_1d_1d_r8
     module procedure interp_1d_2d_r4, interp_1d_2d_r8
     module procedure interp_1d_3d_r4, interp_1d_3d_r8
  end interface interp_1d

  interface find_index
      module procedure find_index_r4, find_index_r8
  end interface find_index

!> @addtogroup axis_utils2_mod
!> @{

contains

  !> @brief Returns X,Y,Z or T cartesian attribute
  subroutine get_axis_cart(fileobj, axisname, cart)
    type(FmsNetcdfFile_t), intent(in) :: fileobj !< file object to read from
    character(len=*), intent(in) :: axisname !< name of axis to retrieve
    character(len=1), intent(out) :: cart !< Returned attribute axis

    character(len=1) :: axis_cart
    character(len=16), dimension(2) :: lon_names, lat_names
    character(len=16), dimension(3) :: z_names
    character(len=16), dimension(2) :: t_names
    character(len=16), dimension(3) :: lon_units, lat_units
    character(len=8) , dimension(4) :: z_units
    character(len=3) , dimension(6) :: t_units
    character(len=32) :: name
    integer :: i

    lon_names = (/'lon','x  '/)
    lat_names = (/'lat','y  '/)
    z_names = (/'depth ','height','z     '/)
    t_names = (/'time','t   '/)
    lon_units = (/'degrees_e   ', 'degrees_east', 'degreese    '/)
    lat_units = (/'degrees_n    ', 'degrees_north', 'degreesn     '/)
    z_units = (/'cm ','m  ','pa ','hpa'/)
    t_units = (/'sec', 'min','hou','day','mon','yea'/)

    cart = "N"
    if (variable_exists(fileobj, axisname)) then
      if (variable_att_exists(fileobj, axisname, "cartesian_axis")) then
        call get_variable_attribute(fileobj, axisname, "cartesian_axis", cart(1:1))
      elseif (variable_att_exists(fileobj, axisname, "axis")) then
        call get_variable_attribute(fileobj, axisname, "axis", cart(1:1))
      endif
      axis_cart = uppercase(cart)
      if (axis_cart .eq. 'X' .or. axis_cart .eq. 'Y' .or. axis_cart .eq. 'Z' &
          .or. axis_cart .eq. 'T') then
        cart = axis_cart
        return
      endif
    endif

    if (cart /= 'X' .and. cart /= 'Y' .and. cart /= 'Z' .and. cart /= 'T') then
       name = lowercase(axisname)
       do i=1,size(lon_names(:))
          if (trim(name(1:3)) == trim(lon_names(i))) cart = 'X'
       enddo
       do i=1,size(lat_names(:))
          if (trim(name(1:3)) == trim(lat_names(i))) cart = 'Y'
       enddo
       do i=1,size(z_names(:))
          if (trim(name) == trim(z_names(i))) cart = 'Z'
       enddo
       do i=1,size(t_names(:))
          if (trim(name) == t_names(i)) cart = 'T'
       enddo
    end if

    if (cart /= 'X' .and. cart /= 'Y' .and. cart /= 'Z' .and. cart /= 'T') then
       name = lowercase(axisname)
       do i=1,size(lon_units(:))
          if (trim(name) == trim(lon_units(i))) cart = 'X'
       enddo
       do i=1,size(lat_units(:))
          if (trim(name) == trim(lat_units(i))) cart = 'Y'
       enddo
       do i=1,size(z_units(:))
          if (trim(name) == trim(z_units(i))) cart = 'Z'
       enddo
       do i=1,size(t_units(:))
          if (name(1:3) == trim(t_units(i))) cart = 'T'
       enddo
    end if
  end subroutine get_axis_cart

  !> @brief Checks if 'modulo' variable exists for a given axis.
  !!
  !> @return true if modulo variable exists in fileobj for the given axis name.
  function get_axis_modulo(fileobj, axisname)
    type(FmsNetcdfFile_t), intent(in) :: fileobj
    character(len=*), intent(in) :: axisname
    logical :: get_axis_modulo

    get_axis_modulo = variable_att_exists(fileobj, axisname, "modulo")
  end function get_axis_modulo

  !> @return true if modulo_beg and modulo_end exist in fileobj with the given
  !! axis, and returns their values in tbeg and tend.
  function get_axis_modulo_times(fileobj, axisname, tbeg, tend)
    type(FmsNetcdfFile_t), intent(in) :: fileobj
    character(len=*), intent(in) :: axisname
    character(len=*), intent(out) :: tbeg, tend
    logical :: get_axis_modulo_times
    logical :: found_tbeg, found_tend

    found_tbeg = variable_att_exists(fileobj, axisname, "modulo_beg")
    found_tend = variable_att_exists(fileobj, axisname, "modulo_end")

    if (found_tbeg .and. .not. found_tend) then
      call mpp_error(FATAL,'error in get: Found modulo_beg but not modulo_end')
    endif
    if (.not. found_tbeg .and. found_tend) then
      call mpp_error(FATAL,'error in get: Found modulo_end but not modulo_beg')
    endif

    if (found_tbeg) then
      call get_variable_attribute(fileobj, axisname, "modulo_beg", tbeg)
      call get_variable_attribute(fileobj, axisname, "modulo_end", tend)
    else
      tbeg = ""
      tend = ""
    endif
    get_axis_modulo_times = found_tbeg
  end function get_axis_modulo_times

#include "axis_utils2_r4.fh"
#include "axis_utils2_r8.fh"

end module axis_utils2_mod
!> @}
! close documentation grouping
