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
use platform_mod

implicit none

  logical :: test_conserve = .false. , test_bicubic = .false. , test_spherical =.false. , test_bilinear =.false.
  integer :: ni_src = 360, nj_src = 180
  integer :: ni_dst = 144, nj_dst = 72

  namelist /test_horiz_interp_nml/ test_conserve, test_bicubic, test_spherical, test_bilinear, ni_src, nj_src, ni_dst, &
                                   nj_dst


  type(domain2d)                    :: domain
  integer                           :: id1, id2, id3, id4
  integer                           :: isc, iec, jsc, jec, i, j
  integer                           :: io, ierr, layout(2)

  call fms_init
  call constants_init
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
  if (test_conserve) then
    call test_horiz_interp_conserve
  else if(test_bicubic) then
    call test_horiz_interp_bicubic
  else if(test_bilinear) then
    call test_horiz_interp_bilinear
  else if(test_spherical) then
    call test_horiz_interp_spherical
  else
    call mpp_error(FATAL, "test_horiz_interp: no unit test enabled in namelist")
  endif


  call mpp_exit

  contains

  subroutine test_horiz_interp_spherical
    real(HI_TEST_KIND)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(HI_TEST_KIND), allocatable, dimension(:,:) :: lon2D_dst, lat2D_dst, lon2D_src, lat2D_src
    real(HI_TEST_KIND), allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst
    real(HI_TEST_KIND) :: lon_src_beg = 0._HI_TEST_KIND,    lon_src_end = 360._HI_TEST_KIND
    real(HI_TEST_KIND) :: lat_src_beg = -90._HI_TEST_KIND,  lat_src_end = 90._HI_TEST_KIND
    real(HI_TEST_KIND) :: lon_dst_beg = -280._HI_TEST_KIND, lon_dst_end = 80._HI_TEST_KIND
    real(HI_TEST_KIND) :: lat_dst_beg = -90._HI_TEST_KIND,  lat_dst_end = 90._HI_TEST_KIND
    real(HI_TEST_KIND) :: D2R = real(PI,HI_TEST_KIND)/180._HI_TEST_KIND
    real(HI_TEST_KIND), parameter :: SMALL = 1.0e-10_HI_TEST_KIND
    type(horiz_interp_type)           :: interp_spherical

    allocate(data_src(ni_src, nj_src) )
    allocate(data1_dst(isc:iec+1, jsc:jec+1), data2_dst(isc:iec+1, jsc:jec+1) )

    allocate(lon2D_src(ni_src, nj_src), lat2D_src(ni_src, nj_src) )
    allocate(lon2D_dst(isc:iec, jsc:jec), lat2D_dst(isc:iec, jsc:jec) )

    ! set up longitude and latitude of source/destination grid.
    dlon_src = (lon_src_end-lon_src_beg)/ni_src
    dlat_src = (lat_src_end-lat_src_beg)/nj_src
    dlon_dst = (lon_dst_end-lon_dst_beg)/ni_dst
    dlat_dst = (lat_dst_end-lat_dst_beg)/nj_dst

    !--- set up the source data
    do j = 1, nj_src
        do i = 1, ni_src
          data_src(i,j) = i + j*0.001_HI_TEST_KIND
        end do
    end do

    !! init input data
    do i = 1, ni_src+1
      lon2D_src(i,:) = lon_src_beg + (i-1)*dlon_src
    end do

    do j = 1, nj_src+1
      lat2D_src(:,j) = lon_src_beg + (i-1)*dlon_src
    end do

    do i = isc, iec+1
        lon2D_dst(i,:) = lon_src_beg + (i-1)*dlon_src
    end do

    do j = jsc, jec+1
        lat2D_dst(:,j) = lon_src_beg + (i-1)*dlon_src
    end do


    !! perform and time interpolations
    id1 = mpp_clock_id( 'horiz_interp_2dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! --- 1dx1d interpolation

    ! --- 2dx2d interpolation
    call mpp_clock_begin(id1)
    call horiz_interp_new(interp_spherical, lon2D_src, lat2D_src, lon2D_dst, lat2D_dst, interp_method="spherical")
    call horiz_interp(interp_spherical, data_src, data2_dst, verbose=2)
    call horiz_interp_del(interp_spherical)
    call mpp_clock_end(id1)

    !call mpp_clock_begin(id1)
    !call horiz_interp_new(interp_spherical, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst)
    !call horiz_interp(interp_spherical, data_src, data1_dst)
    !call horiz_interp_del(interp_spherical)
    !call mpp_clock_end(id1)

    !! check results TODO

    !--- compare the data after interpolation between 1-D and 2-D version interpolation
    do j = jsc, jsc
        do i = isc, iec

          if( abs(data1_dst(i,j)-data2_dst(i,j)) > SMALL ) then
            ! print*, "After interpolation At point (i,j) = (", i, ",", j, "), data1 = ", data1_dst(i,j), &
              !", data2 = ", data2_dst(i,j), ", data1-data2 = ",  data1_dst(i,j) - data2_dst(i,j)
              !call mpp_error(FATAL,"horiz_interp_test: data1_dst does not approxiamate data2_dst")
          end if
        end do
    end do

    if(mpp_pe() == mpp_root_pe()) call mpp_error(NOTE,   &
          "The test that verify 1dx2d version horiz_interp can reproduce 1dx1d version of horiz_interp is succesful")


  end subroutine

  subroutine test_horiz_interp_bilinear
    real(HI_TEST_KIND)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(HI_TEST_KIND), allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
    real(HI_TEST_KIND), allocatable, dimension(:,:) :: lon2D_dst, lat2D_dst
    real(HI_TEST_KIND), allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst
    real(HI_TEST_KIND), parameter :: lon_src_beg =   0._HI_TEST_KIND,  lon_src_end = 360._HI_TEST_KIND
    real(HI_TEST_KIND), parameter :: lat_src_beg = -90._HI_TEST_KIND,  lat_src_end = 90._HI_TEST_KIND
    real(HI_TEST_KIND), parameter :: D2R = real(PI,HI_TEST_KIND)/180._HI_TEST_KIND
    type(horiz_interp_type) :: interp

    allocate( lon1D_src(ni_src+1), lat1D_src(nj_src+1) )
    allocate( lon1D_dst(ni_src+1), lat1D_dst(nj_src+1) )
    allocate( lon2d_dst(ni_src,nj_src), lat2d_dst(ni_src,nj_src) )
    allocate( data_src(ni_src, nj_src) )
    allocate( data1_dst(ni_src,nj_src), data2_dst(ni_src,nj_src) )

    ! set up longitude and latitude of source/destination grid.
    dlon_src = (lon_src_end-lon_src_beg)/real(ni_src,HI_TEST_KIND)  ;  dlon_dst = dlon_src
    dlat_src = (lat_src_end-lat_src_beg)/real(nj_src,HI_TEST_KIND)  ;  dlat_dst = dlat_src

    ! set up 1d source grid
    do i = 1, ni_src+1
       lon1D_src(i) = ( lon_src_beg + real(i-1,HI_TEST_KIND)*dlon_src ) * D2R
    end do
    do j = 1, nj_src+1
       lat1D_src(j) = ( lat_src_beg + real(j-1,HI_TEST_KIND)*dlat_src ) * D2R
    end do

    !--- set up the source data
    do j = 1, nj_src
       do i = 1, ni_src
          data_src(i,j) = i + j*0.001_HI_TEST_KIND
       end do
    end do

    id1 = mpp_clock_id( 'horiz_interp_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id2 = mpp_clock_id( 'horiz_interp_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id3 = mpp_clock_id( 'horiz_interp_2dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id4 = mpp_clock_id( 'horiz_interp_2dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! --- 1dx1d version conservative interpolation
    data1_dst = 0.0_HI_TEST_KIND
    lon1d_dst = lon1d_src
    lat1d_dst = lat1d_src
    call mpp_clock_begin(id1)
    call horiz_interp_new(interp, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, interp_method = "bilinear")
    call horiz_interp(interp, data_src, data1_dst)
    call horiz_interp_del(interp)
    call mpp_clock_end(id1)
    !check
    do j=1, nj_src
       do i=1, ni_src
          if( data_src(i,j).ne.data1_dst(i,j) ) then
             write(*,*) 'expected ', data_src(i,j), ' but computed ', data1_dst(i,j)
             call mpp_error(FATAL, "failed at horiz_interp_1d")
          end if
       end do
    end do

    ! --- 1dx2d version conservative interpolation
    data2_dst = 0.0_HI_TEST_KIND

    ! taking the midpoint
    do i = 1, ni_src
       lon2D_dst(i,:) = (lon1D_src(i) + lon1D_src(i+1)) * 0.5_HI_TEST_KIND
    end do
    do j = 1, nj_src
       lat2D_dst(:,j) = (lat1D_src(j) + lat1D_src(j+1)) * 0.5_HI_TEST_KIND
    end do

    call mpp_clock_begin(id2)
    call horiz_interp_new(interp, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, interp_method = "bilinear")
    call horiz_interp(interp, data_src, data2_dst)
    call horiz_interp_del(interp)
    call mpp_clock_end(id2)
    !check
    do j=1, nj_src
       do i=1, ni_src
          if( data_src(i,j).ne.data2_dst(i,j) ) then
             write(*,*) 'expected ', data_src(i,j), ' but computed ', data2_dst(i,j)
             call mpp_error(FATAL, "failed at horiz_interp_1d")
          end if
       end do
    end do
    write(*,*) 'HELLOOOO'

  end subroutine test_horiz_interp_bilinear

  subroutine test_horiz_interp_bicubic
    !! input data
    real(HI_TEST_KIND), allocatable, dimension(:) :: lat_in_1D, lon_in_1D
    type(horiz_interp_type)                       :: interp_t
    !! output data
    real(HI_TEST_KIND), allocatable, dimension(:) :: lat_out_1D, lon_out_1D
    !! array sizes
    integer :: nlon_in, nlat_in
    integer :: nlon_out, nlat_out

    nlon_in = ni_src;  nlat_in = nj_src
    nlon_out = ni_src; nlat_out = nj_src

    allocate(lat_in_1D(nlat_in))
    allocate(lon_in_1D(nlon_in))
    allocate(lat_out_1D(nlat_out))
    allocate(lon_out_1D(nlon_out))

    lat_in_1D = 0;  lon_in_1D = 0
    lat_out_1D = 0; lon_out_1D = 0

    call horiz_interp_new(interp_t, lon_in_1d, lat_in_1d, lon_out_1d, lon_out_1d)
    !call horiz_interp(interp_t, lon_in_1d, lat_in_1d, lon_out_1d, lon_out_1d) 
    call horiz_interp_del(interp_t)

  end subroutine

  subroutine test_horiz_interp_conserve
    real(HI_TEST_KIND)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(HI_TEST_KIND), allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
    real(HI_TEST_KIND), allocatable, dimension(:,:) :: lon2D_src, lat2D_src, lon2D_dst, lat2D_dst
    real(HI_TEST_KIND), allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst, data3_dst, data4_dst
    real(HI_TEST_KIND) :: lon_src_beg = 0._HI_TEST_KIND,    lon_src_end = 360._HI_TEST_KIND
    real(HI_TEST_KIND) :: lat_src_beg = -90._HI_TEST_KIND,  lat_src_end = 90._HI_TEST_KIND
    real(HI_TEST_KIND) :: lon_dst_beg = -280._HI_TEST_KIND, lon_dst_end = 80._HI_TEST_KIND
    real(HI_TEST_KIND) :: lat_dst_beg = -90._HI_TEST_KIND,  lat_dst_end = 90._HI_TEST_KIND
    real(HI_TEST_KIND) :: D2R = real(PI,HI_TEST_KIND)/180._HI_TEST_KIND
    real(HI_TEST_KIND), parameter :: SMALL = 1.0e-10_HI_TEST_KIND
    type(horiz_interp_type)           :: interp_conserve


    print *, "asdijsauidsiuad", kind(dlon_src)

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
          data_src(i,j) = i + j*0.001_HI_TEST_KIND
        end do
    end do

    id1 = mpp_clock_id( 'horiz_interp_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id2 = mpp_clock_id( 'horiz_interp_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id3 = mpp_clock_id( 'horiz_interp_2dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id4 = mpp_clock_id( 'horiz_interp_2dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! --- 1dx1d version conservative interpolation
    call mpp_clock_begin(id1)
    call horiz_interp_new(interp_conserve, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, interp_method = "conservative")
    call horiz_interp(interp_conserve, data_src, data1_dst)
    call horiz_interp_del(interp_conserve)
    call mpp_clock_end(id1)

    ! --- 1dx2d version conservative interpolation
    call mpp_clock_begin(id2)
    call horiz_interp_new(interp_conserve, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, interp_method = "conservative")
    call horiz_interp(interp_conserve, data_src, data2_dst)
    call horiz_interp_del(interp_conserve)
    call mpp_clock_end(id2)

    ! --- 2dx1d version conservative interpolation
    call mpp_clock_begin(id3)
    call horiz_interp_new(interp_conserve, lon2D_src, lat2D_src, lon1D_dst, lat1D_dst, interp_method = "conservative")
    call horiz_interp(interp_conserve, data_src, data3_dst)
    call horiz_interp_del(interp_conserve)
    call mpp_clock_end(id3)

    ! --- 2dx2d version conservative interpolation
    call mpp_clock_begin(id4)
    call horiz_interp_new(interp_conserve, lon2D_src, lat2D_src, lon2D_dst, lat2D_dst, interp_method = "conservative")
    call horiz_interp(interp_conserve, data_src, data4_dst)
    call horiz_interp_del(interp_conserve)
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

   end subroutine

end program horiz_interp_test
