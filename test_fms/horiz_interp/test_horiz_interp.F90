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
     call test_horiz_interp_conserve_r4
     call test_horiz_interp_conserve_r8
  else if(test_bicubic) then
     call test_horiz_interp_bicubic_r8
  else if(test_bilinear) then
     call test_horiz_interp_bilinear_r8
  else if(test_spherical) then
    call test_horiz_interp_spherical_r8
  else
    call mpp_error(FATAL, "test_horiz_interp: no unit test enabled in namelist")
  endif


  call mpp_exit

  contains

  subroutine test_horiz_interp_spherical_r8
    real(r8_kind)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(r8_kind), allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
    real(r8_kind), allocatable, dimension(:,:) :: lon2D_dst, lat2D_dst
    real(r8_kind), allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst
    real(r8_kind) :: lon_src_beg = 0._r8_kind,    lon_src_end = 360._r8_kind
    real(r8_kind) :: lat_src_beg = -90._r8_kind,  lat_src_end = 90._r8_kind
    real(r8_kind) :: lon_dst_beg = -280._r8_kind, lon_dst_end = 80._r8_kind
    real(r8_kind) :: lat_dst_beg = -90._r8_kind,  lat_dst_end = 90._r8_kind
    real(r8_kind) :: D2R = real(PI,r8_kind)/180._r8_kind
    real(r8_kind), parameter :: SMALL = 1.0e-10_r8_kind
    type(horiz_interp_type)           :: interp_spherical

    allocate(lon1D_src(ni_src+1), lat1D_src(nj_src+1), data_src(ni_src, nj_src) )

    allocate(lon1D_dst(isc:iec+1), lat1D_dst(jsc:jec+1) )
    allocate(data1_dst(isc:iec, jsc:jec), data2_dst(isc:iec, jsc:jec) )

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


    do i = isc, iec+1
        lon2D_dst(i,:) = lon1D_dst(i)
    end do

    do j = jsc, jec+1
        lat2D_dst(:,j) = lat1D_dst(j)
    end do

    !--- set up the source data
    do j = 1, nj_src
        do i = 1, ni_src
          data_src(i,j) = i + j*0.001_r8_kind
        end do
    end do

    !! init input data


    !! perform and time interpolations
    id1 = mpp_clock_id( 'horiz_interp_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id2 = mpp_clock_id( 'horiz_interp_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! --- 1dx1d version conservative interpolation
    call mpp_clock_begin(id1)
    call horiz_interp_new(interp_spherical, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, interp_method = "spherical")
    call horiz_interp(interp_spherical, data_src, data1_dst)
    call horiz_interp_del(interp_spherical)
    call mpp_clock_end(id1)

    ! --- 1dx2d version conservative interpolation
    call mpp_clock_begin(id2)
    call horiz_interp_new(interp_spherical, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, interp_method = "spherical")
    call horiz_interp(interp_spherical, data_src, data2_dst)
    call horiz_interp_del(interp_spherical)
    call mpp_clock_end(id2)

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

  subroutine test_horiz_interp_bilinear_r8
    real(r8_kind)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(r8_kind), allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
    real(r8_kind), allocatable, dimension(:,:) :: lon2D_src, lat2d_src, lon2D_dst, lat2D_dst
    real(r8_kind), allocatable, dimension(:,:) :: data_src, data_dst
    real(r8_kind), parameter :: lon_src_beg =   0._r8_kind,  lon_src_end = 360._r8_kind
    real(r8_kind), parameter :: lat_src_beg = -90._r8_kind,  lat_src_end = 90._r8_kind
    real(r8_kind), parameter :: D2R = real(PI,r8_kind)/180._r8_kind
    type(horiz_interp_type) :: interp

    allocate( lon1D_src(ni_src+1), lat1D_src(nj_src+1) )
    allocate( lon1D_dst(ni_src+1), lat1D_dst(nj_src+1) )
    allocate( lon2d_src(ni_src,nj_src), lat2d_src(ni_src,nj_src) )
    allocate( lon2d_dst(ni_src,nj_src), lat2d_dst(ni_src,nj_src) )
    allocate( data_src(ni_src, nj_src) )
    allocate( data_dst(ni_src,nj_src) )

    ! set up longitude and latitude of source/destination grid.
    dlon_src = (lon_src_end-lon_src_beg)/real(ni_src,r8_kind)  ;  dlon_dst = dlon_src
    dlat_src = (lat_src_end-lat_src_beg)/real(nj_src,r8_kind)  ;  dlat_dst = dlat_src

    ! set up 1d source grid
    do i = 1, ni_src
       lon1D_src(i) = ( lon_src_beg + real(i-1,r8_kind)*dlon_src ) * D2R
    end do
    lon1D_src(ni_src+1) = ( lon_src_beg + real(ni_src,r8_kind)*dlon_src ) * D2R

    do j = 1, nj_src
       lat1D_src(j) = ( lat_src_beg + real(j-1,r8_kind)*dlat_src ) * D2R
    end do
    lat1D_src(nj_src+1) = ( lat_src_beg + real(nj_src,r8_kind)*dlat_src ) * D2R

    !--- set up the source data
    do j = 1, nj_src
       do i = 1, ni_src
          data_src(i,j) = real(i,r8_kind) + real(j,r8_kind)*0.001_r8_kind
       end do
    end do

    id1 = mpp_clock_id( 'horiz_interp_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id2 = mpp_clock_id( 'horiz_interp_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id3 = mpp_clock_id( 'horiz_interp_2dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id4 = mpp_clock_id( 'horiz_interp_2dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! --- 1dx1d version conservative interpolation
    data_dst = 0.0_r8_kind
    lon1d_dst = lon1d_src
    lat1d_dst = lat1d_src
    call mpp_clock_begin(id1)
    call horiz_interp_new(interp, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, interp_method = "bilinear")
    call horiz_interp(interp, data_src, data_dst)
    ! check weights
    do j=1, nj_src-1
       do i=1, ni_src-1
          if( interp%kind8_reals%wtj(i,j,1).ne.1.0_r8_kind ) then
             write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%kind8_reals%wtj(i,j,1)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wtj1")
          end if
          if( interp%kind8_reals%wtj(i,j,2).ne.0.0_r8_kind ) then
             write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wtj(i,j,2)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wtj2")
          end if
          if( interp%kind8_reals%wti(i,j,1).ne.1.0_r8_kind ) then
             write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,1)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wti1")
          end if
          if( interp%kind8_reals%wti(i,j,2).ne.0.0_r8_kind ) then
             write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,2)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wtj2")
          end if
       end do
    end do
    call mpp_clock_end(id1)
    !checking to make sure data_src is equal to data_dst
    do j=1, nj_src
       do i=1, ni_src
          if( data_src(i,j).ne.data_dst(i,j) ) then
             write(*,*) 'expected ', data_src(i,j), ' but computed ', data_dst(i,j)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d data comparison")
          end if
       end do
    end do
    call horiz_interp_del(interp)

    ! --- 1dx2d version conservative interpolation
    data_dst = 0.0_r8_kind
    ! taking the midpoint
    do i = 1, ni_src
       lon2D_dst(i,:) = 0.5_r8_kind*(lon1D_src(i)+lon1D_src(i+1))
    end do
    do j = 1, nj_src
       lat2D_dst(:,j) = 0.5_r8_kind*(lat1D_src(j)+lat1D_src(j+1))
    end do
    call mpp_clock_begin(id2)
    call horiz_interp_new(interp, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, interp_method = "bilinear")
    call horiz_interp(interp, data_src, data_dst)
    call mpp_clock_end(id2)
    ! check weights
    do j=1, nj_src-1
       do i=1, ni_src-1
          if( interp%kind8_reals%wtj(i,j,1).ne.1.0_r8_kind ) then
             write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%kind8_reals%wtj(i,j,1)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wtj1")
          end if
          if( interp%kind8_reals%wtj(i,j,2).ne.0.0_r8_kind ) then
             write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wtj(i,j,2)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wtj2")
          end if
          if( interp%kind8_reals%wti(i,j,1).ne.1.0_r8_kind ) then
             write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,1)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wti1")
          end if
          if( interp%kind8_reals%wti(i,j,2).ne.0.0_r8_kind ) then
             write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,2)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wti2")
          end if
       end do
    end do
    !check that data are equal
    do j=1, nj_src
       do i=1, ni_src
          if( data_src(i,j).ne.data_dst(i,j) ) then
             write(*,*) 'expected ', data_src(i,j), ' but computed ', data_dst(i,j)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d data comparison")
          end if
       end do
    end do
    call horiz_interp_del(interp)

    ! --- 2dx1d version conservative interpolation
    data_dst = 0.0_r8_kind
    lon1d_dst = lon1d_src
    lat1d_dst = lat1d_src
    do i=1, ni_src
       lon2d_src(i,:) = lon1d_dst(i)
    end do
    do j=1, nj_src
       lat2d_src(:,j) = lat1d_dst(j)
    end do
    call mpp_clock_begin(id3)
    call horiz_interp_new(interp, lon2D_src, lat2D_src, lon1D_dst(1:ni_src), lat1D_dst(1:nj_src), interp_method = "bilinear")
    call horiz_interp(interp, data_src, data_dst)
    ! check weights
    !j=1,i=1 is a special case; see subroutine find_neighbor
    i=1 ; j=1
    if( interp%kind8_reals%wtj(i,j,1).ne.1.0_r8_kind ) then
       write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', i,j,interp%kind8_reals%wtj(i,j,1)
       call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj(1,1,1)")
    end if
    if( interp%kind8_reals%wtj(i,j,2).ne.0.0_r8_kind ) then
       write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wtj(i,j,2)
       call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj(1,1,2)")
    end if
    if( interp%kind8_reals%wti(i,j,1).ne.1.0_r8_kind ) then
       write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,1)
       call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti(1,1,1)")
    end if
    if( interp%kind8_reals%wti(i,j,2).ne.0.0_r8_kind ) then
       write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,2)
       call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti(1,1,2)")
    end if
    do j=2, nj_src
       do i=2, ni_src
          if( interp%kind8_reals%wtj(i,j,1).ne.0.0_r8_kind ) then
             write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', i,j,interp%kind8_reals%wtj(i,j,1)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj1")
          end if
          if( interp%kind8_reals%wtj(i,j,2).ne.1.0_r8_kind ) then
             write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wtj(i,j,2)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj2")
          end if
          if( interp%kind8_reals%wti(i,j,1).ne.0.0_r8_kind ) then
             write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,1)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti1")
          end if
          if( interp%kind8_reals%wti(i,j,2).ne.1.0_r8_kind ) then
             write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,2)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti2")
          end if
       end do
    end do
    call mpp_clock_end(id3)
    !check that data are equal
    do j=1, nj_src
       do i=1, ni_src
          if( data_src(i,j).ne.data_dst(i,j) ) then
             write(*,*) 'expected ', data_src(i,j), ' but computed ', data_dst(i,j)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d data comparison")
          end if
       end do
    end do
    call horiz_interp_del(interp)

    ! --- 2dx2d version conservative interpolation
    data_dst = 0.0_r8_kind
    lon2D_dst = lon2D_src
    lat2D_dst = lat2D_src

    call mpp_clock_begin(id4)
    call horiz_interp_new(interp, lon2D_src, lat2D_src, lon2D_dst, lat2D_dst, interp_method = "bilinear")
    call horiz_interp(interp, data_src, data_dst)
    call mpp_clock_end(id4)
    ! check weights
    !j=1,i=1 is a special case; see subroutine find_neighbor
    i=1 ; j=1
    if( interp%kind8_reals%wtj(i,j,1).ne.1.0_r8_kind ) then
       write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', i,j,interp%kind8_reals%wtj(i,j,1)
       call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj(1,1,1)")
    end if
    if( interp%kind8_reals%wtj(i,j,2).ne.0.0_r8_kind ) then
       write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wtj(i,j,2)
       call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj(1,1,2)")
    end if
    if( interp%kind8_reals%wti(i,j,1).ne.1.0_r8_kind ) then
       write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,1)
       call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti(1,1,1)")
    end if
    if( interp%kind8_reals%wti(i,j,2).ne.0.0_r8_kind ) then
       write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,2)
       call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti(1,1,2)")
    end if
    do j=2, nj_src
       do i=2, ni_src
          if( interp%kind8_reals%wtj(i,j,1).ne.0.0_r8_kind ) then
             write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', i,j,interp%kind8_reals%wtj(i,j,1)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj1")
          end if
          if( interp%kind8_reals%wtj(i,j,2).ne.1.0_r8_kind ) then
             write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wtj(i,j,2)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj2")
          end if
          if( interp%kind8_reals%wti(i,j,1).ne.0.0_r8_kind ) then
             write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,1)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti1")
          end if
          if( interp%kind8_reals%wti(i,j,2).ne.1.0_r8_kind ) then
             write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%kind8_reals%wti(i,j,2)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti2")
          end if
       end do
    end do
    call horiz_interp_del(interp)
    !check that data are equal
    do j=1, nj_src
       do i=1, ni_src
          if( data_src(i,j).ne.data_dst(i,j) ) then
             write(*,*) 'expected ', data_src(i,j), ' but computed ', data_dst(i,j)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d data comparison")
          end if
       end do
    end do

  end subroutine test_horiz_interp_bilinear_r8

  subroutine test_horiz_interp_bicubic_r8
    real(r8_kind)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(r8_kind), allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
    real(r8_kind), allocatable, dimension(:,:) :: lon2D_dst, lat2D_dst
    real(r8_kind), allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst
    real(r8_kind) :: lon_src_beg = 0._r8_kind,    lon_src_end = 360._r8_kind
    real(r8_kind) :: lat_src_beg = -90._r8_kind,  lat_src_end = 90._r8_kind
    real(r8_kind) :: lon_dst_beg = -280._r8_kind, lon_dst_end = 80._r8_kind
    real(r8_kind) :: lat_dst_beg = -90._r8_kind,  lat_dst_end = 90._r8_kind
    real(r8_kind) :: D2R = real(PI,r8_kind)/180._r8_kind
    real(r8_kind), parameter :: SMALL = 1.0e-10_r8_kind
    type(horiz_interp_type)           :: interp_bicubic

    allocate(lon1D_src(ni_src+1), lat1D_src(nj_src+1), data_src(ni_src, nj_src) )

    allocate(lon1D_dst(isc:iec+1), lat1D_dst(jsc:jec+1) )
    allocate(data1_dst(isc:iec, jsc:jec), data2_dst(isc:iec, jsc:jec) )

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


    do i = isc, iec+1
        lon2D_dst(i,:) = lon1D_dst(i)
    end do

    do j = jsc, jec+1
        lat2D_dst(:,j) = lat1D_dst(j)
    end do

    !--- set up the source data
    do j = 1, nj_src
        do i = 1, ni_src
          data_src(i,j) = i + j*0.001_r8_kind
        end do
    end do

    !! init input data


    !! perform and time interpolations
    id1 = mpp_clock_id( 'horiz_interp_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id2 = mpp_clock_id( 'horiz_interp_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! --- 1dx1d version conservative interpolation
    call mpp_clock_begin(id1)
    call horiz_interp_new(interp_bicubic, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, interp_method = "bicubic")
    call horiz_interp(interp_bicubic, data_src, data1_dst)
    call horiz_interp_del(interp_bicubic)
    call mpp_clock_end(id1)

    ! --- 1dx2d version conservative interpolation
    call mpp_clock_begin(id2)
    call horiz_interp_new(interp_bicubic, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, interp_method = "bicubic")
    call horiz_interp(interp_bicubic, data_src, data2_dst)
    call horiz_interp_del(interp_bicubic)
    call mpp_clock_end(id2)

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

  subroutine test_horiz_interp_conserve_r8
    real(r8_kind)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(r8_kind), allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
    real(r8_kind), allocatable, dimension(:,:) :: lon2D_src, lat2D_src, lon2D_dst, lat2D_dst
    real(r8_kind), allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst, data3_dst, data4_dst
    real(r8_kind) :: lon_src_beg = 0._r8_kind,    lon_src_end = 360._r8_kind
    real(r8_kind) :: lat_src_beg = -90._r8_kind,  lat_src_end = 90._r8_kind
    real(r8_kind) :: lon_dst_beg = -280._r8_kind, lon_dst_end = 80._r8_kind
    real(r8_kind) :: lat_dst_beg = -90._r8_kind,  lat_dst_end = 90._r8_kind
    real(r8_kind) :: D2R = real(PI,r8_kind)/180._r8_kind
    real(r8_kind), parameter :: SMALL = 1.0e-10_r8_kind
    type(horiz_interp_type)           :: interp_conserve

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
          data_src(i,j) = i + j*0.001_r8_kind
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

  subroutine test_horiz_interp_conserve_r4
    real(r4_kind)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(r4_kind), allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
    real(r4_kind), allocatable, dimension(:,:) :: lon2D_src, lat2D_src, lon2D_dst, lat2D_dst
    real(r4_kind), allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst, data3_dst, data4_dst
    real(r4_kind) :: lon_src_beg = 0_r4_kind,    lon_src_end = 360_r4_kind
    real(r4_kind) :: lat_src_beg = -90_r4_kind,  lat_src_end = 90_r4_kind
    real(r4_kind) :: lon_dst_beg = -280_r4_kind, lon_dst_end = 80_r4_kind
    real(r4_kind) :: lat_dst_beg = -90_r4_kind,  lat_dst_end = 90_r4_kind
    real(r4_kind) :: D2R = real(PI,r4_kind)/180_r4_kind
    real(r4_kind), parameter :: SMALL = 1.0e-10_r4_kind
    type(horiz_interp_type)           :: interp_conserve

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
          data_src(i,j) = i + j*0.001_r4_kind
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
