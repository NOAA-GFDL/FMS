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

program horiz_interp_test

use mpp_mod,          only : mpp_init, mpp_exit, mpp_error, FATAL, stdout, mpp_npes
use mpp_mod,          only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,          only : mpp_pe, mpp_root_pe, NOTE, MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
use mpp_mod,          only : input_nml_file
use mpp_domains_mod,  only : mpp_define_layout, mpp_define_domains, mpp_get_compute_domain
use mpp_domains_mod,  only : mpp_domains_init, domain2d
use fms_mod,          only : check_nml_error, fms_init
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_new, horiz_interp_del
use horiz_interp_mod, only : horiz_interp, horiz_interp_type
use constants_mod,    only : constants_init, PI

implicit none

  integer :: ni_src = 360, nj_src = 180
  integer :: ni_dst = 144, nj_dst = 72

  namelist /test_horiz_interp_nml/ ni_src, nj_src, ni_dst, nj_dst

  real :: lon_src_beg = 0,    lon_src_end = 360
  real :: lat_src_beg = -90,  lat_src_end = 90
  real :: lon_dst_beg = -280, lon_dst_end = 80
  real :: lat_dst_beg = -90,  lat_dst_end = 90
  real :: D2R = PI/180.
  real, parameter :: SMALL = 1.0e-10

  type(domain2d)                    :: domain
  type(horiz_interp_type)           :: Interp
  integer                           :: id1, id2, id3, id4
  integer                           :: isc, iec, jsc, jec, i, j
  integer                           :: io, ierr, layout(2)
  real                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
  real, allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
  real, allocatable, dimension(:,:) :: lon2D_src, lat2D_src, lon2D_dst, lat2D_dst
  real, allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst, data3_dst, data4_dst

  call constants_init
  call mpp_init
  call mpp_domains_init
  call fms_init
  call horiz_interp_init

  !--- read namelist
  read (input_nml_file, test_horiz_interp_nml, iostat=io)
  ierr = check_nml_error(io, 'test_horiz_interp_nml')

  !--- define domains
  call mpp_define_layout( (/1, ni_dst, 1, nj_dst/), mpp_npes(), layout)
  call mpp_define_domains((/1, ni_dst, 1, nj_dst/), layout, domain)
  call mpp_get_compute_domain(domain,isc,iec,jsc,jec)

  !--- test conservative horiz_interp with a simple test. the source grid is the region
  !    (0:360,-90:90) with grid size ni_src, nj_src ( default 360X180). and the destination
  !    is the region (-280:80, -90:90) with grid size ni_dstXnj_dst( default 144X72).
  !    integer checksum and global sum will be printed out for both the 1D and 2D version.

  allocate(lon2D_src(ni_src+1, nj_src+1), lat2D_src(ni_src+1, nj_src+1) )
  allocate(lon1D_src(ni_src+1), lat1D_src(nj_src+1), data_src(ni_src, nj_src) )

  allocate(lon2D_dst(isc:iec+1, jsc:jec+1), lat2D_dst(isc:iec+1, jsc:jec+1) )
  allocate(lon1D_dst(isc:iec+1), lat1D_dst(jsc:jec+1) )
  allocate(data1_dst(isc:iec, jsc:jec), data2_dst(isc:iec, jsc:jec) )
  allocate(data3_dst(isc:iec, jsc:jec), data4_dst(isc:iec, jsc:jec) )

  ! set up longitude and latitude of source/destination grid.
  dlon_src = (lon_src_end-lon_src_beg)/ni_src
  dlat_src = (lat_src_end-lat_src_beg)/nj_src
  dlon_dst = (lon_dst_end-lon_dst_beg)/ni_dst
  dlat_dst = (lat_dst_end-lat_dst_beg)/nj_dst

  do i = 1, ni_src+1
     lon1D_src(i) = lon_src_beg + (i-1)*dlon_src
  end do

  do j = 1, nj_src+1
     lat1D_src(j) = lat_src_beg + (j-1)*dlat_src
  end do

  do i = isc, iec+1
     lon1D_dst(i) = lon_dst_beg + (i-1)*dlon_dst
  end do

  do j = jsc, jec+1
     lat1D_dst(j) = lat_dst_beg + (j-1)*dlat_dst
  end do

  ! scale grid to radians.
  lon1D_src = lon1D_src * D2R
  lat1D_src = lat1D_src * D2R
  lon1D_dst = lon1D_dst * D2R
  lat1D_dst = lat1D_dst * D2R

  do i = 1, ni_src+1
     lon2D_src(i,:) = lon1D_src(i)
  end do

  do j = 1, nj_src+1
     lat2D_src(:,j) = lat1D_src(j)
  end do

  do i = isc, iec+1
     lon2D_dst(i,:) = lon1D_dst(i)
  end do

  do j = jsc, jec+1
     lat2D_dst(:,j) = lat1D_dst(j)
  end do

  !--- set up the source data
  do j = 1, nj_src
     do i = 1, ni_src
        data_src(i,j) = i + j*0.001
     end do
  end do

  id1 = mpp_clock_id( 'horiz_interp_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  id2 = mpp_clock_id( 'horiz_interp_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  id3 = mpp_clock_id( 'horiz_interp_2dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  id4 = mpp_clock_id( 'horiz_interp_2dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

  ! --- 1dx1d version conservative interpolation
  call mpp_clock_begin(id1)
  call horiz_interp_new(Interp, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, interp_method = "conservative")
  call horiz_interp(Interp, data_src, data1_dst)
  call horiz_interp_del(Interp)
  call mpp_clock_end(id1)

  ! --- 1dx2d version conservative interpolation
  call mpp_clock_begin(id2)
  call horiz_interp_new(Interp, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, interp_method = "conservative")
  call horiz_interp(Interp, data_src, data2_dst)
  call horiz_interp_del(Interp)
  call mpp_clock_end(id2)

  ! --- 2dx1d version conservative interpolation
  call mpp_clock_begin(id3)
  call horiz_interp_new(Interp, lon2D_src, lat2D_src, lon1D_dst, lat1D_dst, interp_method = "conservative")
  call horiz_interp(Interp, data_src, data3_dst)
  call horiz_interp_del(Interp)
  call mpp_clock_end(id3)

  ! --- 2dx2d version conservative interpolation
  call mpp_clock_begin(id4)
  call horiz_interp_new(Interp, lon2D_src, lat2D_src, lon2D_dst, lat2D_dst, interp_method = "conservative")
  call horiz_interp(Interp, data_src, data4_dst)
  call horiz_interp_del(Interp)
  call mpp_clock_end(id4)

  !--- compare the data after interpolation between 1-D and 2-D version interpolation
  do j = jsc, jsc
     do i = isc, iec

        if( abs(data1_dst(i,j)-data2_dst(i,j)) > SMALL ) then
           print*, "After interpolation At point (i,j) = (", i, ",", j, "), data1 = ", data1_dst(i,j), &
           ", data2 = ", data2_dst(i,j), ", data1-data2 = ",  data1_dst(i,j) - data2_dst(i,j)
           call mpp_error(FATAL,"horiz_interp_test: data1_dst does not approxiamate data2_dst")
        end if
     end do
  end do

  if(mpp_pe() == mpp_root_pe()) call mpp_error(NOTE,   &
       "The test that verify 1dx2d version horiz_interp can reproduce 1dx1d version of horiz_interp is succesful")

  do j = jsc, jsc
     do i = isc, iec

        if( abs(data1_dst(i,j)-data3_dst(i,j)) > SMALL ) then
           print*, "After interpolation At point (i,j) = (", i, ",", j, "), data1 = ", data1_dst(i,j), &
           ", data2 = ", data3_dst(i,j), ", data1-data2 = ",  data1_dst(i,j) - data3_dst(i,j)
           call mpp_error(FATAL,"horiz_interp_test: data1_dst does not approxiamate data3_dst")
        end if
     end do
  end do

  if(mpp_pe() == mpp_root_pe()) call mpp_error(NOTE,   &
       "The test that verify 2dx1d version horiz_interp can reproduce 1dx1d version of horiz_interp is succesful")

  do j = jsc, jsc
     do i = isc, iec

        if( abs(data1_dst(i,j)-data4_dst(i,j)) > SMALL ) then
           print*, "After interpolation At point (i,j) = (", i, ",", j, "), data1 = ", data1_dst(i,j), &
           ", data2 = ", data4_dst(i,j), ", data1-data2 = ",  data1_dst(i,j) - data4_dst(i,j)
           call mpp_error(FATAL,"horiz_interp_test: data1_dst does not approxiamate data4_dst")
        end if
     end do
  end do

  if(mpp_pe() == mpp_root_pe()) call mpp_error(NOTE,   &
       "The test that verify 2dx2d version horiz_interp can reproduce 1dx1d version of horiz_interp is succesful")

  call mpp_exit

end program horiz_interp_test
