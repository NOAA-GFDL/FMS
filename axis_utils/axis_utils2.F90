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
!!
!> FMS2_IO equivalent version of @ref axis_utils_mod.
!> @author M.J. Harrison

!> @file
!> @brief File for @ref axis_utils2_mod

!> @addtogroup axis_utils2_mod
!> @{
module axis_utils2_mod
  use, intrinsic :: iso_fortran_env
  use mpp_mod,    only: mpp_error, FATAL, stdout
  use fms_mod,    only: lowercase, uppercase, string_array_index, fms_error_handler
  use fms2_io_mod, only: FmsNetcdfDomainFile_t, variable_att_exists, FmsNetcdfFile_t, &
                         get_variable_num_dimensions, get_variable_attribute,  &
                         get_variable_size, read_data, variable_exists
  use platform_mod

  implicit none

  public get_axis_cart, get_axis_modulo, lon_in_range, &
         tranlon, frac_index, nearest_index, interp_1d, get_axis_modulo_times, axis_edges

  private

  integer, parameter :: maxatts = 100
  real, parameter    :: epsln= 1.e-10
  real, parameter    :: fp5 = 0.5, f360 = 360.0
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
  interface interp_1d
     module procedure interp_1d_1d
     module procedure interp_1d_2d
     module procedure interp_1d_3d
  end interface

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

  !> get axis edge data from a given file
  subroutine axis_edges(fileobj, name, edge_data, reproduce_null_char_bug_flag)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object to read from
  character(len=*), intent(in) :: name !< Name of a given axis
  class(*), dimension(:), intent(out) :: edge_data !< Returned edge data from given axis name
  logical, intent(in), optional :: reproduce_null_char_bug_flag !< Flag indicating to reproduce
                                     !! the mpp_io bug where the null characters were not removed
                                     !! after reading a string attribute

  integer :: ndims
  character(len=128) :: buffer
  integer, dimension(:), allocatable :: dim_sizes
  real(kind=r4_kind), dimension(:), allocatable :: r32
  real(kind=r4_kind), dimension(:,:), allocatable :: r322d
  real(kind=r8_kind), dimension(:), allocatable :: r64
  real(kind=r8_kind), dimension(:,:), allocatable :: r642d
  integer :: i
  integer :: n
  logical :: reproduce_null_char_bug !< Local flag indicating to reproduce the mpp_io bug where
                                     !! the null characters were not removed after reading a string attribute

  ndims = get_variable_num_dimensions(fileobj, name)
  allocate(dim_sizes(ndims))
  call get_variable_size(fileobj, name, dim_sizes)
  n = dim_sizes(1)
  if (size(edge_data) .ne. n+1) then
    call mpp_error(FATAL, "axis_edge: incorrect size of edge_data array.")
  endif
  deallocate(dim_sizes)

  reproduce_null_char_bug = .false.
  if (present(reproduce_null_char_bug_flag)) reproduce_null_char_bug = reproduce_null_char_bug_flag

  buffer = ""
  if (variable_att_exists(fileobj, name, "edges")) then
    !! If the reproduce_null_char_bug flag is turned on fms2io will not remove the null character
    call get_variable_attribute(fileobj, name, "edges", buffer(1:128), &
        reproduce_null_char_bug_flag=reproduce_null_char_bug)

    !! Check for a null character here, if it exists *_bnds will be calculated instead of read in
    if (reproduce_null_char_bug) then
        i = 0
        i = index(buffer, char(0))
        if (i > 0) buffer = ""
    endif
  elseif (variable_att_exists(fileobj, name, "bounds")) then
    !! If the reproduce_null_char_bug flag is turned on fms2io will not remove the null character
    call get_variable_attribute(fileobj, name, "bounds", buffer(1:128), &
        reproduce_null_char_bug_flag=reproduce_null_char_bug)

    !! Check for a null character here, if it exists *_bnds will be calculated instead of read in
    if (reproduce_null_char_bug) then
        i = 0
        i = index(buffer, char(0))
        if (i > 0) buffer = ""
    endif
  endif
  if (trim(buffer) .ne. "") then
    ndims = get_variable_num_dimensions(fileobj, buffer)
    allocate(dim_sizes(ndims))
    call get_variable_size(fileobj, buffer, dim_sizes)
    if (size(dim_sizes) .eq. 1) then
      if (dim_sizes(1) .ne. n+1) then
        call mpp_error(FATAL, "axis_edges: incorrect size of edge data.")
      endif
      call read_data(fileobj, buffer, edge_data)
    elseif (size(dim_sizes) .eq. 2) then
      if (dim_sizes(1) .ne. 2) then
        call mpp_error(FATAL, "axis_edges: first dimension of edge must be of size 2")
      endif
      if (dim_sizes(2) .ne. n) then
        call mpp_error(FATAL, "axis_edges: incorrect size of edge data.")
      endif
      select type (edge_data)
        type is (real(kind=r4_kind))
          allocate(r322d(dim_sizes(1), dim_sizes(2)))
          call read_data(fileobj, buffer, r322d)
          edge_data(1:dim_sizes(2)) = r322d(1,:)
          edge_data(dim_sizes(2)+1) = r322d(2,dim_sizes(2))
          deallocate(r322d)
        type is (real(kind=r8_kind))
          allocate(r642d(dim_sizes(1), dim_sizes(2)))
          call read_data(fileobj, buffer, r642d)
          edge_data(1:dim_sizes(2)) = r642d(1,:)
          edge_data(dim_sizes(2)+1) = r642d(2,dim_sizes(2))
          deallocate(r642d)
        class default
          call mpp_error(FATAL, "axis_edges: unsupported kind.")
      end select
    endif
    deallocate(dim_sizes)
  else
    select type (edge_data)
      type is (real(kind=r4_kind))
        allocate(r32(n))
        call read_data(fileobj, name, r32)
        do i = 2, n
          edge_data(i) = r32(i-1) + 0.5_r4_kind*(r32(i) - r32(i-1))
        enddo
        edge_data(1) = r32(1) - 0.5_r4_kind*(r32(2) - r32(1))
        if (abs(edge_data(1)) .lt. 1.e-10) then
          edge_data(1) = 0._r4_kind
        endif
        edge_data(n+1) = r32(n) + 0.5_r4_kind*(r32(n) - r32(n-1))
        deallocate(r32)
      type is (real(kind=r8_kind))
        allocate(r64(n))
        call read_data(fileobj, name, r64)
        do i = 2, n
          edge_data(i) = r64(i-1) + 0.5_r8_kind*(r64(i) - r64(i-1))
        enddo
        edge_data(1) = r64(1) - 0.5_r8_kind*(r64(2) - r64(1))
        if (abs(edge_data(1)) .lt. 1.d-10) then
          edge_data(1) = 0._r8_kind
        endif
        edge_data(n+1) = r64(n) + 0.5_r8_kind*(r64(n) - r64(n-1))
        deallocate(r64)
      class default
        call mpp_error(FATAL, "axis_edges: unsupported kind.")
    end select
  endif
end subroutine axis_edges

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

  !> @brief Returns lon_strt <= longitude <= lon_strt+360
  !! @return real lon_in_range
  function lon_in_range(lon, l_strt)
    real, intent(in) :: lon, l_strt
    real :: lon_in_range
    real :: l_end

    lon_in_range = lon
    l_end = l_strt+360.

    if (abs(lon_in_range - l_strt) < 1.e-4) then
      lon_in_range = l_strt
      return
    endif

    if (abs(lon_in_range - l_end) < 1.e-4) then
      lon_in_range = l_strt
      return
    endif

    do
      if (lon_in_range < l_strt) then
        lon_in_range = lon_in_range +  f360
      else if (lon_in_range  >  l_end) then
        lon_in_range  = lon_in_range - f360
      else
        exit
      end if
    end do

  end function lon_in_range

  !> @brief Returns monotonic array of longitudes s.t., lon_strt <= lon(:) <= lon_strt+360.
  !!
  !> <br>The first istrt-1 entries are moved to the end of the array:
  !!
  !! e.g.
  !!        lon =      0 1 2 3 4 5  ...  358 359; lon_strt = 3 ==>
  !!        tranlon =  3 4 5 6 7 8  ...  359 360 361 362; istrt = 4
  subroutine tranlon(lon, lon_start, istrt)

    ! returns array of longitudes s.t.  lon_strt <= lon < lon_strt+360.
    ! also, the first istrt-1 entries are moved to the end of the array
    !
    ! e.g.
    !        lon =      0 1 2 3 4 5  ...  358 359; lon_strt = 3 ==>
    !        tranlon =  3 4 5 6 7 8  ...  359 360 361 362; istrt = 4

    real, intent(inout), dimension(:) :: lon
    real, intent(in) :: lon_start
    integer, intent(out) :: istrt


    integer :: len, i
    real :: lon_strt, tmp(size(lon(:))-1)

    len = size(lon(:))

    do i=1,len
       lon(i) = lon_in_range(lon(i),lon_start)
    enddo

    istrt=0
    do i=1,len-1
       if (lon(i+1) < lon(i)) then
          istrt=i+1
          exit
       endif
    enddo

    if (istrt>1) then ! grid is not monotonic
       if (abs(lon(len)-lon(1)) < epsln) then
          tmp = cshift(lon(1:len-1),istrt-1)
          lon(1:len-1) = tmp
          lon(len) = lon(1)
       else
          lon = cshift(lon,istrt-1)
       endif
       lon_strt = lon(1)
       do i=2,len+1
          lon(i) = lon_in_range(lon(i),lon_strt)
          lon_strt = lon(i)
       enddo
    endif

    return
  end subroutine tranlon

  !>     nearest_index = index of nearest data point within "array" corresponding to
  !!            "value".
  !!
  !!     inputs:
  !!
  !!     value  = arbitrary data...same units as elements in "array"
  !!     array  = array of data points  (must be monotonically increasing)
  !!
  !!     output:
  !!
  !!     nearest_index =  index of nearest data point to "value"
  !!             if "value" is outside the domain of "array" then nearest_index = 1
  !!             or "ia" depending on whether array(1) or array(ia) is
  !!             closest to "value"
  !!
  !!             note: if "array" is dimensioned array(0:ia) in the calling
  !!                   program, then the returned index should be reduced
  !!                   by one to account for the zero base.
  !!
  !!     example:
  !!
  !!     let model depths be defined by the following:
  !!     parameter (km=5)
  !!     dimension z(km)
  !!     data z /5.0, 10.0, 50.0, 100.0, 250.0/
  !!
  !!     k1 = nearest_index (12.5, z, km)
  !!     k2 = nearest_index (0.0, z, km)
  !!
  !!     k1 would be set to 2, and k2 would be set to 1 so that
  !!     z(k1) would be the nearest data point to 12.5 and z(k2) would
  !!   be the nearest data point to 0.0
  !!
  !!   @return real frac_index
  function frac_index (value, array)

    integer :: ia, i, ii, unit
    real :: value !< arbitrary data...same units as elements in "array"
    real :: frac_index
    real, dimension(:) :: array !< array of data points  (must be monotonically increasing)
    logical keep_going

    ia = size(array(:))

    do i=2,ia
       if (array(i) < array(i-1)) then
          unit = stdout()
          write (unit,*) '=> Error: "frac_index" array must be monotonically increasing when searching for nearest value to ',&
                              value
          write (unit,*) '          array(i) < array(i-1) for i=',i
          write (unit,*) '          array(i) for i=1..ia follows:'
          do ii=1,ia
             write (unit,*) 'i=',ii, ' array(i)=',array(ii)
          enddo
          call mpp_error(FATAL,' "frac_index" array must be monotonically increasing.')
       endif
    enddo
    if (value < array(1) .or. value > array(ia)) then
!       if (value < array(1))  frac_index = 1.
!       if (value > array(ia)) frac_index = float(ia)
        frac_index = -1.0
    else
       i=1
       keep_going = .true.
       do while (i <= ia .and. keep_going)
          i = i+1
          if (value <= array(i)) then
             frac_index = float(i-1) + (value-array(i-1))/(array(i)-array(i-1))
             keep_going = .false.
          endif
       enddo
    endif
  end function frac_index

  !> @brief Return index of nearest point along axis
  !!
  !>     nearest_index = index of nearest data point within "array" corresponding to
  !!            "value".
  !!
  !!     inputs:
  !!
  !!     value  = arbitrary data...same units as elements in "array"
  !!     array  = array of data points  (must be monotonically increasing)
  !!     ia     = dimension of "array"
  !!
  !!     output:
  !!
  !!     nearest_index =  index of nearest data point to "value"
  !!             if "value" is outside the domain of "array" then nearest_index = 1
  !!             or "ia" depending on whether array(1) or array(ia) is
  !!             closest to "value"
  !!
  !!             note: if "array" is dimensioned array(0:ia) in the calling
  !!                   program, then the returned index should be reduced
  !!                   by one to account for the zero base.
  !!
  !!     example:
  !!
  !!     let model depths be defined by the following:
  !!     parameter (km=5)
  !!     dimension z(km)
  !!     data z /5.0, 10.0, 50.0, 100.0, 250.0/
  !!
  !!     k1 = nearest_index (12.5, z, km)
  !!     k2 = nearest_index (0.0, z, km)
  !!
  !!     k1 would be set to 2, and k2 would be set to 1 so that
  !!     z(k1) would be the nearest data point to 12.5 and z(k2) would
  !!     be the nearest data point to 0.0
  !! @return integer nearest_index
  function nearest_index (value, array)

    integer :: nearest_index
    integer :: ia !< dimension of "array"
    integer :: i, ii, unit
    real :: value !< arbitrary data...same units as elements in "array"
    real, dimension(:) :: array !< array of data points  (must be monotonically increasing)
    logical keep_going

    ia = size(array(:))

    do i=2,ia
       if (array(i) < array(i-1)) then
          unit = stdout()
          write (unit,*) '=> Error: "nearest_index" array must be monotonically increasing &
                         &when searching for nearest value to ',value
          write (unit,*) '          array(i) < array(i-1) for i=',i
          write (unit,*) '          array(i) for i=1..ia follows:'
          do ii=1,ia
             write (unit,*) 'i=',ii, ' array(i)=',array(ii)
          enddo
          call mpp_error(FATAL,' "nearest_index" array must be monotonically increasing.')
       endif
    enddo
    if (value < array(1) .or. value > array(ia)) then
       if (value < array(1))  nearest_index = 1
       if (value > array(ia)) nearest_index = ia
    else
       i=1
       keep_going = .true.
       do while (i <= ia .and. keep_going)
          i = i+1
          if (value <= array(i)) then
             nearest_index = i
             if (array(i)-value > value-array(i-1)) nearest_index = i-1
             keep_going = .false.
          endif
       enddo
    endif
  end function nearest_index

  !#############################################################################

  subroutine interp_1d_linear(grid1,grid2,data1,data2)

    real, dimension(:),    intent(in) :: grid1, data1, grid2
    real, dimension(:), intent(inout) :: data2

    integer :: n1, n2, i, n, ext
    real :: w

    n1 = size(grid1(:))
    n2 = size(grid2(:))


    do i=2,n1
       if (grid1(i) <= grid1(i-1)) call mpp_error(FATAL, 'grid1 not monotonic')
    enddo

    do i=2,n2
       if (grid2(i) <= grid2(i-1)) call mpp_error(FATAL, 'grid2 not monotonic')
    enddo

    if (grid1(1) > grid2(1) ) call mpp_error(FATAL, 'grid2 lies outside grid1')
    if (grid1(n1) < grid2(n2) ) call mpp_error(FATAL, 'grid2 lies outside grid1')

    do i=1,n2
       n = nearest_index(grid2(i),grid1)

       if (grid1(n) < grid2(i)) then
          w = (grid2(i)-grid1(n))/(grid1(n+1)-grid1(n))
          data2(i) = (1.-w)*data1(n) + w*data1(n+1)
       else
          if(n==1) then
             data2(i) = data1(n)
          else
             w = (grid2(i)-grid1(n-1))/(grid1(n)-grid1(n-1))
             data2(i) = (1.-w)*data1(n-1) + w*data1(n)
          endif
       endif
    enddo


    return

  end subroutine interp_1d_linear

  !###################################################################
  subroutine interp_1d_cubic_spline(grid1, grid2, data1, data2, yp1, ypn)

    real, dimension(:),    intent(in) :: grid1, grid2, data1
    real, dimension(:), intent(inout) :: data2
    real,                  intent(in) :: yp1, ypn

    real, dimension(size(grid1))      :: y2, u
    real                              :: sig, p, qn, un, h, a ,b
    integer                           :: n, m, i, k, klo, khi

    n = size(grid1(:))
    m = size(grid2(:))

    do i=2,n
       if (grid1(i) <= grid1(i-1)) call mpp_error(FATAL, 'grid1 not monotonic')
    enddo

    do i=2,m
       if (grid2(i) <= grid2(i-1)) call mpp_error(FATAL, 'grid2 not monotonic')
    enddo

    if (grid1(1) > grid2(1) ) call mpp_error(FATAL, 'grid2 lies outside grid1')
    if (grid1(n) < grid2(m) ) call mpp_error(FATAL, 'grid2 lies outside grid1')

    if (yp1 >.99e30) then
       y2(1)=0.
       u(1)=0.
    else
       y2(1)=-0.5
       u(1)=(3./(grid1(2)-grid1(1)))*((data1(2)-data1(1))/(grid1(2)-grid1(1))-yp1)
    endif

    do i=2,n-1
       sig=(grid1(i)-grid1(i-1))/(grid1(i+1)-grid1(i-1))
       p=sig*y2(i-1)+2.
       y2(i)=(sig-1.)/p
       u(i)=(6.*((data1(i+1)-data1(i))/(grid1(i+1)-grid1(i))-(data1(i)-data1(i-1)) &
             /(grid1(i)-grid1(i-1)))/(grid1(i+1)-grid1(i-1))-sig*u(i-1))/p
    enddo

    if (ypn > .99e30) then
       qn=0.
       un=0.
    else
       qn=0.5
       un=(3./(grid1(n)-grid1(n-1)))*(ypn-(data1(n)-data1(n-1))/(grid1(n)-grid1(n-1)))
    endif

    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

    do  k=n-1,1,-1
       y2(k)=y2(k)*y2(k+1)+u(k)
    enddo

    do k = 1, m
       n = nearest_index(grid2(k),grid1)
       if (grid1(n) < grid2(k)) then
          klo = n
       else
          if(n==1) then
            klo = n
          else
            klo = n -1
          endif
       endif
       khi = klo+1
       h   = grid1(khi)-grid1(klo)
       a   = (grid1(khi) - grid2(k))/h
       b   = (grid2(k) - grid1(klo))/h
       data2(k) = a*data1(klo) + b*data1(khi)+ ((a**3-a)*y2(klo) + (b**3-b)*y2(khi))*(h**2)/6.
    enddo

  end subroutine interp_1d_cubic_spline

  !###################################################################

  subroutine interp_1d_1d(grid1,grid2,data1,data2, method, yp1, yp2)

    real, dimension(:),      intent(in)    :: grid1, data1, grid2
    real, dimension(:),      intent(inout) :: data2
    character(len=*), optional, intent(in) :: method
    real,             optional, intent(in) :: yp1, yp2

    real              :: y1, y2
    character(len=32) :: interp_method
    integer           :: k2, ks, ke

    k2 = size(grid2(:))

    interp_method = "linear"
    if(present(method)) interp_method = method
    y1 = 1.0e30
    if(present(yp1)) y1 = yp1
    y2 = 1.0e30
    if(present(yp2)) y2 = yp2
    call find_index(grid1, grid2(1), grid2(k2), ks, ke)
    select case(trim(interp_method))
    case("linear")
       call interp_1d_linear(grid1(ks:ke),grid2,data1(ks:ke),data2)
    case("cubic_spline")
       call interp_1d_cubic_spline(grid1(ks:ke),grid2,data1(ks:ke),data2, y1, y2)
    case default
       call mpp_error(FATAL,"axis_utils: interp_method should be linear or cubic_spline")
    end select

    return

  end subroutine interp_1d_1d

  !###################################################################


  subroutine interp_1d_2d(grid1,grid2,data1,data2)

    real, dimension(:,:),    intent(in) :: grid1, data1, grid2
    real, dimension(:,:), intent(inout) :: data2

    integer :: n1, n2, i, n, k2, ks, ke
    real :: w

    n1 = size(grid1,1)
    n2 = size(grid2,1)
    k2 = size(grid2,2)

    if (n1 /= n2) call mpp_error(FATAL,'grid size mismatch')

    do n=1,n1
       call find_index(grid1(n,:), grid2(n,1), grid2(n,k2), ks, ke)
       call interp_1d_linear(grid1(n,ks:ke),grid2(n,:),data1(n,ks:ke),data2(n,:))
    enddo

    return

  end subroutine interp_1d_2d

  !###################################################################

  subroutine interp_1d_3d(grid1,grid2,data1,data2, method, yp1, yp2)

    real, dimension(:,:,:),  intent(in)    :: grid1, data1, grid2
    real, dimension(:,:,:),  intent(inout) :: data2
    character(len=*), optional, intent(in) :: method
    real,             optional, intent(in) :: yp1, yp2

    integer           :: n1, n2, m1, m2, k2, i, n, m
    real              :: w, y1, y2
    character(len=32) :: interp_method
    integer           :: ks, ke
    n1 = size(grid1,1)
    n2 = size(grid2,1)
    m1 = size(grid1,2)
    m2 = size(grid2,2)
    k2 = size(grid2,3)

    interp_method = "linear"
    if(present(method)) interp_method = method
    y1 = 1.0e30
    if(present(yp1)) y1 = yp1
    y2 = 1.0e30
    if(present(yp2)) y2 = yp2

    if (n1 /= n2 .or. m1 /= m2) call mpp_error(FATAL,'grid size mismatch')

    select case(trim(interp_method))
    case("linear")
       do m=1,m1
          do n=1,n1
            call find_index(grid1(n,m,:), grid2(n,m,1), grid2(n,m,k2), ks, ke)
             call interp_1d_linear(grid1(n,m,ks:ke),grid2(n,m,:),data1(n,m,ks:ke),data2(n,m,:))
          enddo
       enddo
    case("cubic_spline")
       do m=1,m1
          do n=1,n1
            call find_index(grid1(n,m,:), grid2(n,m,1), grid2(n,m,k2), ks, ke)
            call interp_1d_cubic_spline(grid1(n,m,ks:ke),grid2(n,m,:), data1(n,m,ks:ke),data2(n,m,:), y1, y2)
          enddo
       enddo
    case default
       call mpp_error(FATAL,"axis_utils: interp_method should be linear or cubic_spline")
    end select

    return

  end subroutine interp_1d_3d


  !#####################################################################
  subroutine find_index(grid1, xs, xe, ks, ke)
    real, dimension(:), intent(in) :: grid1
    real,               intent(in) :: xs, xe
    integer,           intent(out) :: ks, ke

    integer :: k, nk

    nk = size(grid1(:))

    ks = 0; ke = 0
    do k = 1, nk-1
       if(grid1(k) <= xs .and. grid1(k+1) > xs ) then
          ks = k
          exit
       endif
    enddo
    do k = nk, 2, -1
       if(grid1(k) >= xe .and. grid1(k-1) < xe ) then
          ke = k
          exit
       endif
    enddo

    if(ks == 0 ) call mpp_error(FATAL,' xs locate outside of grid1')
    if(ke == 0 ) call mpp_error(FATAL,' xe locate outside of grid1')

  end subroutine find_index

end module axis_utils2_mod
!> @}
! close documentation grouping
