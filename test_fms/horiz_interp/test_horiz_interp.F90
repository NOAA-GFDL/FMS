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
!> @author Ryan Mulhall 2023
!> Original test is in test_conserve, modified to test the other 3 interp_method option and mixed precision reals
!! tests are split up by interp_method (same way the modules are broken up) and enabled via the nml flags.
!! Assignment test checks that the override is copying the data type properly
!! TODO some larger tests with different data sets

!! defaults to 8 real kind, make check will compile with both 4 and 8
#ifndef HI_TEST_KIND_
#define HI_TEST_KIND_ 8
#endif

program horiz_interp_test

use mpp_mod,          only : mpp_init, mpp_exit, mpp_error, FATAL, stdout, mpp_npes, WARNING
use mpp_mod,          only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,          only : mpp_pe, mpp_root_pe, NOTE, MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
use mpp_mod,          only : input_nml_file, mpp_sync
use mpp_domains_mod,  only : mpp_define_layout, mpp_define_domains, mpp_get_compute_domain
use mpp_domains_mod,  only : mpp_domains_init, domain2d
use fms_mod,          only : check_nml_error, fms_init
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_new, horiz_interp_del
use horiz_interp_mod, only : horiz_interp, horiz_interp_type
use horiz_interp_spherical_mod, only: horiz_interp_spherical_wght
use horiz_interp_type_mod, only: SPHERICA
use constants_mod,    only : constants_init, PI
use platform_mod

implicit none

  logical :: test_conserve = .false. , test_bicubic = .false. , test_spherical =.false. , test_bilinear =.false.
  logical :: test_assign = .false.
  logical :: test_solo = .false.!< test with the 'solo' wrappers that hide the _new and _del calls for the derived type
  integer :: ni_src = 360, nj_src = 180
  integer :: ni_dst = 144, nj_dst = 72
  integer, parameter :: max_neighbors = 400 !! took this from spherical mod
                !! max amount found neighbors to loop through in spherical search
  logical :: decreasing_lat = .False. !< latitude and longitude are decreasing instead of increasing
                                      !! This is only for the bilinear test case

  namelist /test_horiz_interp_nml/ test_conserve, test_bicubic, test_spherical, test_bilinear, test_assign, test_solo,&
                                   ni_src, nj_src, ni_dst,nj_dst, decreasing_lat


  type(domain2d)                    :: domain
  integer                           :: id1, id2, id3, id4
  integer                           :: isc, iec, jsc, jec, i, j
  integer                           :: io, ierr, layout(2)
  integer, parameter :: lkind = HI_TEST_KIND_

  call fms_init
  call mpp_init
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
  else if(test_assign) then
    call test_assignment
  else
    call mpp_error(FATAL, "test_horiz_interp: no unit test enabled in namelist")
  endif

  call mpp_exit

  contains

  !> Tests spherical module interpolation with each dimension conversion
  !! test without passing in the type when test_solo is true
  !! The spherical module has a nml option for whether using a full or radially bounded search
  !! for finding the nearest points and distances so this gets run for both
  subroutine test_horiz_interp_spherical
    !! grid data
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: lat_in_2D, lon_in_2D
    type(horiz_interp_type)                       :: interp_t
    !! input data
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: data_src, data_dst
    !! output data
    real(HI_TEST_KIND_), allocatable, dimension(:,:)   :: lat_out_2D, lon_out_2D
    real(HI_TEST_KIND_), allocatable, dimension(:,:,:) :: wghts
    !! array sizes and number of lat/lon per index
    real(HI_TEST_KIND_) :: dlon_src, dlat_src, dlon_dst, dlat_dst
    !! parameters for lon/lat setup
    real(HI_TEST_KIND_) :: lon_src_beg = 0._lkind,    lon_src_end = 360._lkind
    real(HI_TEST_KIND_) :: lat_src_beg = -90._lkind,  lat_src_end = 90._lkind
    real(HI_TEST_KIND_) :: lon_dst_beg = -280._lkind, lon_dst_end = 80._lkind
    real(HI_TEST_KIND_) :: lat_dst_beg = -90._lkind,  lat_dst_end = 90._lkind
    real(HI_TEST_KIND_) :: D2R = real(PI,HI_TEST_KIND_)/180._lkind
    real(HI_TEST_KIND_) :: R2D = 180._lkind/real(PI,HI_TEST_KIND_)
    real(HI_TEST_KIND_), parameter :: SMALL = 1.0e-10_lkind

    ! set up longitude and latitude of source/destination grid.
    dlon_src = (lon_src_end-lon_src_beg)/real(ni_src, HI_TEST_KIND_)
    dlat_src = (lat_src_end-lat_src_beg)/real(nj_src, HI_TEST_KIND_)
    dlon_dst = (lon_dst_end-lon_dst_beg)/real(ni_dst, HI_TEST_KIND_)
    dlat_dst = (lat_dst_end-lat_dst_beg)/real(nj_dst, HI_TEST_KIND_)

    ! set up 2d lon/lat
    allocate(lon_in_2D(ni_src, nj_src), lat_in_2D(ni_src, nj_src))
    do i = 1, ni_src
        lon_in_2D(i,:) = lon_src_beg + real(i-1, HI_TEST_KIND_)*dlon_src
    end do
    do j = 1, nj_src
        lat_in_2D(:,j) = lat_src_beg + real(j-1, HI_TEST_KIND_)*dlat_src
    end do
    allocate(lon_out_2D(ni_dst, nj_dst), lat_out_2D(ni_dst, nj_dst))
    do i = 1, ni_dst
        lon_out_2D(i,:) = lon_dst_beg + real(i-1, HI_TEST_KIND_)*dlon_dst
    end do
    do j = 1, nj_dst
        lat_out_2D(:,j) = lat_src_beg + real(j-1, HI_TEST_KIND_)*dlat_dst
    end do

    ! scale to radians
    lat_in_2D = lat_in_2D * D2R
    lon_in_2D = lon_in_2D * D2R
    lat_out_2D = lat_out_2D * D2R
    lon_out_2D = lon_out_2D * D2R


    allocate(data_src(ni_src, nj_src))
    allocate(data_dst(ni_dst, nj_dst))
    allocate(wghts(ni_dst, nj_dst, max_neighbors))
    data_dst = 0.0_lkind ; data_src = 1.0_lkind

    id1 = mpp_clock_id( 'horiz_interp_spherical_2dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! 2D x 2D (only one supported for spherical)
    call mpp_clock_begin(id1)
    if(.not. test_solo) then
        call horiz_interp_new(interp_t, lon_in_2d, lat_in_2d, lon_out_2d, lon_out_2d, interp_method="spherical")
        call horiz_interp(interp_t, data_src, data_dst)
        call horiz_interp_spherical_wght(interp_t, wghts, verbose=1)
    else
        call horiz_interp(data_src, lon_in_2D, lat_in_2D, lon_out_2D, lat_out_2D, data_dst, interp_method="spherical")
    endif
    call mpp_clock_end(id1)
    do i=1, ni_dst-1
        do j=1, nj_dst-1
            if(data_dst(i,j) - 1.0_lkind .gt. SMALL) then
                print *, 'data_dst(i=', i, ', j=', j, ')=', data_dst(i,j), ' Expected value: 1.0'
                call mpp_error(FATAL, "test_horiz_interp_spherical: "// &
                                                                    "invalid output data after interpolation")
            endif
        enddo
    enddo
    if(.not. test_solo) then
        call horiz_interp_del(interp_t)
        call check_dealloc(interp_t)
    endif
    deallocate(data_src, data_dst)
    deallocate(lat_in_2D, lon_in_2D)
    deallocate(lat_out_2D, lon_out_2D)

  end subroutine

  !> Tests bilinear module interpolation with each dimension conversion
  !! test without passing in the type when test_solo is true
  subroutine test_horiz_interp_bilinear
    real(HI_TEST_KIND_)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(HI_TEST_KIND_), allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: lon2D_src, lat2d_src, lon2D_dst, lat2D_dst
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: data_src, data_dst
    real(HI_TEST_KIND_) :: lon_src_beg =  0._lkind,  lon_src_end = 360.0_lkind
    real(HI_TEST_KIND_) :: lat_src_beg = -90._lkind,  lat_src_end = 90._lkind
    real(HI_TEST_KIND_), parameter :: D2R = real(PI,lkind)/180._lkind
    real(HI_TEST_KIND_), parameter :: R2D = 180._lkind/real(PI,lkind)

    type(horiz_interp_type) :: interp

    if (decreasing_lat) then
        lon_src_beg = 360.0_lkind
        lon_src_end = 0._lkind
        lat_src_beg = 90._lkind
        lat_src_end = -90._lkind
    endif

    allocate( lon1D_src(ni_src+1), lat1D_src(nj_src+1) )
    allocate( lon1D_dst(ni_src+1), lat1D_dst(nj_src+1) )
    allocate( lon2d_src(ni_src,nj_src), lat2d_src(ni_src,nj_src) )
    allocate( lon2d_dst(ni_src,nj_src), lat2d_dst(ni_src,nj_src) )
    allocate( data_src(ni_src, nj_src) )
    allocate( data_dst(ni_src,nj_src) )

    ! set up longitude and latitude of source/destination grid.
    dlon_src = (lon_src_end-lon_src_beg)/real(ni_src,HI_TEST_KIND_)  ;  dlon_dst = dlon_src
    dlat_src = (lat_src_end-lat_src_beg)/real(nj_src,HI_TEST_KIND_)  ;  dlat_dst = dlat_src

    ! set up 1d source grid
    do i = 1, ni_src
       lon1D_src(i) = ( lon_src_beg + real(i-1,HI_TEST_KIND_)*dlon_src ) * D2R
    end do
    lon1D_src(ni_src+1) = ( lon_src_beg + real(ni_src,HI_TEST_KIND_)*dlon_src ) * D2R

    do j = 1, nj_src
       lat1D_src(j) = ( lat_src_beg + real(j-1,HI_TEST_KIND_)*dlat_src ) * D2R
    end do
    lat1D_src(nj_src+1) = ( lat_src_beg + real(nj_src,HI_TEST_KIND_)*dlat_src ) * D2R

    !--- set up the source data
    do j = 1, nj_src
       do i = 1, ni_src
          data_src(i,j) = real(i,HI_TEST_KIND_) + real(j,HI_TEST_KIND_)*0.001_lkind
       end do
    end do

    id1 = mpp_clock_id( 'horiz_interp_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id2 = mpp_clock_id( 'horiz_interp_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id3 = mpp_clock_id( 'horiz_interp_2dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id4 = mpp_clock_id( 'horiz_interp_2dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! --- 1dx1d version bilinear interpolation
    data_dst = 0.0_lkind
    lon1d_dst = lon1d_src
    lat1d_dst = lat1d_src
    call mpp_clock_begin(id1)
    if (.not. test_solo) then
        call horiz_interp_new(interp, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, interp_method = "bilinear")
        call horiz_interp(interp, data_src, data_dst)
    else
        call horiz_interp(data_src, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, data_dst, interp_method = "bilinear")
    endif
    ! check weights
    if( .not. test_solo) then
        do j=1, nj_src-1
            do i=1, ni_src-1
                if(interp%horizInterpReals8_type%is_allocated) then
                    if( interp%horizInterpReals8_type%wtj(i,j,1).ne.1.0_r8_kind ) then
                        write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wtj(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wtj1")
                    end if
                    if( interp%horizInterpReals8_type%wtj(i,j,2).ne.0.0_r8_kind ) then
                        write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wtj(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wtj2")
                    end if
                        if( interp%horizInterpReals8_type%wti(i,j,1).ne.1.0_r8_kind ) then
                        write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wti1")
                    end if
                    if( interp%horizInterpReals8_type%wti(i,j,2).ne.0.0_r8_kind ) then
                        write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wtj2")
                    end if
                else
                    if( interp%horizInterpReals4_type%wtj(i,j,1).ne.1.0_r4_kind ) then
                        write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wtj(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wtj1")
                    end if
                    if( interp%horizInterpReals4_type%wtj(i,j,2).ne.0.0_r4_kind ) then
                        write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wtj(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wtj2")
                    end if
                    if( interp%horizInterpReals4_type%wti(i,j,1).ne.1.0_r4_kind ) then
                        write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wti1")
                    end if
                    if( interp%horizInterpReals4_type%wti(i,j,2).ne.0.0_r4_kind ) then
                        write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d1d with wtj2")
                    end if
                endif
            end do
        end do
    endif
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
    if(.not. test_solo) then
        call horiz_interp_del(interp)
        call check_dealloc(interp)
    endif

    ! --- 1dx2d version bilinear interpolation
    data_dst = 0.0_lkind
    ! taking the midpoint
    do i = 1, ni_src
       lon2D_dst(i,:) = (lon1D_src(i) + lon1D_src(i+1)) * 0.5_lkind
    end do
    do j = 1, nj_src
       lat2D_dst(:,j) = (lat1D_src(j) + lat1D_src(j+1)) * 0.5_lkind
    end do
    call mpp_clock_begin(id2)
    if(.not. test_solo) then
        call horiz_interp_new(interp, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, interp_method = "bilinear")
        call horiz_interp(interp, data_src, data_dst)
    else
        call horiz_interp(data_src, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, data_dst,interp_method="bilinear")
    endif
    call mpp_clock_end(id2)
    ! check weights
    if(.not. test_solo) then
        do j=1, nj_src-1
            do i=1, ni_src-1
                if(interp%horizInterpReals8_type%is_allocated) then
                    if( interp%horizInterpReals8_type%wtj(i,j,1).ne.1.0_r8_kind ) then
                        write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wtj(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wtj1")
                    end if
                    if( interp%horizInterpReals8_type%wtj(i,j,2).ne.0.0_r8_kind ) then
                        write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wtj(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wtj2")
                    end if
                    if( interp%horizInterpReals8_type%wti(i,j,1).ne.1.0_r8_kind ) then
                        write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wti1")
                    end if
                    if( interp%horizInterpReals8_type%wti(i,j,2).ne.0.0_r8_kind ) then
                        write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wti2")
                    end if
                else
                    if( interp%horizInterpReals4_type%wtj(i,j,1).ne.1.0_r4_kind ) then
                        write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wtj(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wtj1")
                    end if
                    if( interp%horizInterpReals4_type%wtj(i,j,2).ne.0.0_r4_kind ) then
                        write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wtj(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wtj2")
                    end if
                    if( interp%horizInterpReals4_type%wti(i,j,1).ne.1.0_r4_kind ) then
                        write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wti1")
                    end if
                    if( interp%horizInterpReals4_type%wti(i,j,2).ne.0.0_r4_kind ) then
                        write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d with wti2")
                    end if
                endif
            end do
        end do
    endif
    !check that data are equal
    do j=1, nj_src
        do i=1, ni_src
            if( data_src(i,j).ne.data_dst(i,j) ) then
                write(*,*) 'expected ', data_src(i,j), ' but computed ', data_dst(i,j)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_1d2d data comparison")
            end if
        end do
    end do
    if(.not. test_solo) then
        call horiz_interp_del(interp)
        call check_dealloc(interp)
    endif

    if (decreasing_lat) return
    ! --- 2dx1d version bilinear interpolation
    data_dst = 0.0_lkind
    lon1d_dst = lon1d_src
    lat1d_dst = lat1d_src
    do i=1, ni_src
       lon2d_src(i,:) = lon1d_dst(i)
    end do
    do j=1, nj_src
       lat2d_src(:,j) = lat1d_dst(j)
    end do
    call mpp_clock_begin(id3)
    if(.not. test_solo) then
        call horiz_interp_new(interp,lon2D_src,lat2D_src,lon1D_dst(1:ni_src),lat1D_dst(1:nj_src), &
                              interp_method = "bilinear")
        call horiz_interp(interp, data_src, data_dst)
    else
        call horiz_interp(data_src, lon2D_src, lat2d_src, lon1D_dst(1:ni_src),lat1D_dst(1:nj_src), data_dst, &
                          interp_method="bilinear")
    endif
    call mpp_clock_end(id3)
    ! check weights
    !j=1,i=1 is a special case; see subroutine find_neighbor
    if(.not. test_solo) then
        i=1 ; j=1
        if(interp%horizInterpReals8_type%is_allocated) then
            if( interp%horizInterpReals8_type%wtj(i,j,1).ne.1.0_r8_kind ) then
                write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', i,j,interp%horizInterpReals8_type%wtj(i,j,1)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj(1,1,1)")
            end if
            if( interp%horizInterpReals8_type%wtj(i,j,2).ne.0.0_r8_kind ) then
                write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wtj(i,j,2)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj(1,1,2)")
            end if
            if( interp%horizInterpReals8_type%wti(i,j,1).ne.1.0_r8_kind ) then
                write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,1)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti(1,1,1)")
            end if
            if( interp%horizInterpReals8_type%wti(i,j,2).ne.0.0_r8_kind ) then
                write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,2)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti(1,1,2)")
            end if
        else
            if( interp%horizInterpReals4_type%wtj(i,j,1).ne.1.0_r4_kind ) then
                write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', i,j,interp%horizInterpReals4_type%wtj(i,j,1)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj(1,1,1)")
            end if
            if( interp%horizInterpReals4_type%wtj(i,j,2).ne.0.0_r4_kind ) then
                write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wtj(i,j,2)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj(1,1,2)")
            end if
            if( interp%horizInterpReals4_type%wti(i,j,1).ne.1.0_r4_kind ) then
                write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,1)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti(1,1,1)")
            end if
            if( interp%horizInterpReals4_type%wti(i,j,2).ne.0.0_r4_kind ) then
                write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,2)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti(1,1,2)")
            end if
        endif
        do j=2, nj_src
            do i=2, ni_src
                if(interp%horizInterpReals8_type%is_allocated) then
                    if( interp%horizInterpReals8_type%wtj(i,j,1).ne.0.0_r8_kind ) then
                        write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', i,j, &
                                    interp%horizInterpReals8_type%wtj(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj1")
                    end if
                    if( interp%horizInterpReals8_type%wtj(i,j,2).ne.1.0_r8_kind ) then
                        write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wtj(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj2")
                    end if
                    if( interp%horizInterpReals8_type%wti(i,j,1).ne.0.0_r8_kind ) then
                        write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti1")
                    end if
                    if( interp%horizInterpReals8_type%wti(i,j,2).ne.1.0_r8_kind ) then
                        write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti2")
                    end if
                else
                    if( interp%horizInterpReals4_type%wtj(i,j,1).ne.0.0_r4_kind ) then
                        write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', i,j, &
                                    interp%horizInterpReals4_type%wtj(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj1")
                    end if
                    if( interp%horizInterpReals4_type%wtj(i,j,2).ne.1.0_r4_kind ) then
                        write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wtj(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wtj2")
                    end if
                    if( interp%horizInterpReals4_type%wti(i,j,1).ne.0.0_r4_kind ) then
                        write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti1")
                    end if
                    if( interp%horizInterpReals4_type%wti(i,j,2).ne.1.0_r4_kind ) then
                        write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d with wti2")
                    end if
                endif
            end do
        end do
    endif
    !check that data are equal
    do j=1, nj_src
       do i=1, ni_src
          if( data_src(i,j).ne.data_dst(i,j) ) then
             write(*,*) 'expected ', data_src(i,j), ' but computed ', data_dst(i,j)
             call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d1d data comparison")
          end if
       end do
    end do
    if(.not. test_solo) then
        call horiz_interp_del(interp)
        call check_dealloc(interp)
    endif

    ! --- 2dx2d version bilinear interpolation
    data_dst = 0.0_lkind
    lon2D_dst = lon2D_src
    lat2D_dst = lat2D_src

    call mpp_clock_begin(id4)
    if(.not. test_solo) then
        call horiz_interp_new(interp, lon2D_src, lat2D_src, lon2D_dst, lat2D_dst, interp_method = "bilinear")
        call horiz_interp(interp, data_src, data_dst)
    else
        call horiz_interp(data_src, lon2D_src, lat2d_src, lon2D_dst, lat2D_dst, data_dst, interp_method="bilinear")
    endif
    call mpp_clock_end(id4)
    ! check weights
    if(.not. test_solo) then
        !j=1,i=1 is a special case; see subroutine find_neighbor
        i=1 ; j=1
        if(interp%horizInterpReals8_type%is_allocated) then
            if( interp%horizInterpReals8_type%wtj(i,j,1).ne.1.0_r8_kind ) then
                write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', i,j,interp%horizInterpReals8_type%wtj(i,j,1)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj(1,1,1)")
            end if
            if( interp%horizInterpReals8_type%wtj(i,j,2).ne.0.0_r8_kind ) then
                write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wtj(i,j,2)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj(1,1,2)")
            end if
            if( interp%horizInterpReals8_type%wti(i,j,1).ne.1.0_r8_kind ) then
                write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,1)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti(1,1,1)")
            end if
            if( interp%horizInterpReals8_type%wti(i,j,2).ne.0.0_r8_kind ) then
                write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,2)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti(1,1,2)")
            end if
        else
            if( interp%horizInterpReals4_type%wtj(i,j,1).ne.1.0_r4_kind ) then
                write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', i,j,interp%horizInterpReals4_type%wtj(i,j,1)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj(1,1,1)")
            end if
            if( interp%horizInterpReals4_type%wtj(i,j,2).ne.0.0_r4_kind ) then
                write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wtj(i,j,2)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj(1,1,2)")
            end if
            if( interp%horizInterpReals4_type%wti(i,j,1).ne.1.0_r4_kind ) then
                write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,1)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti(1,1,1)")
            end if
            if( interp%horizInterpReals4_type%wti(i,j,2).ne.0.0_r4_kind ) then
                write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,2)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti(1,1,2)")
            end if
        endif
        do j=2, nj_src
            do i=2, ni_src
                if(interp%horizInterpReals8_type%is_allocated) then
                    if( interp%horizInterpReals8_type%wtj(i,j,1).ne.0.0_r8_kind ) then
                        write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', i,j, &
                                    interp%horizInterpReals8_type%wtj(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj1")
                    end if
                    if( interp%horizInterpReals8_type%wtj(i,j,2).ne.1.0_r8_kind ) then
                        write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wtj(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj2")
                    end if
                    if( interp%horizInterpReals8_type%wti(i,j,1).ne.0.0_r8_kind ) then
                        write(*,*) 'expected ', 1.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti1")
                    end if
                    if( interp%horizInterpReals8_type%wti(i,j,2).ne.1.0_r8_kind ) then
                        write(*,*) 'expected ', 0.0_r8_kind, ' but computed ', interp%horizInterpReals8_type%wti(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti2")
                    end if
                else
                    if( interp%horizInterpReals4_type%wtj(i,j,1).ne.0.0_r4_kind ) then
                        write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', i,j, &
                                    interp%horizInterpReals4_type%wtj(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj1")
                    end if
                    if( interp%horizInterpReals4_type%wtj(i,j,2).ne.1.0_r4_kind ) then
                        write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wtj(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wtj2")
                    end if
                    if( interp%horizInterpReals4_type%wti(i,j,1).ne.0.0_r4_kind ) then
                        write(*,*) 'expected ', 1.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,1)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti1")
                    end if
                    if( interp%horizInterpReals4_type%wti(i,j,2).ne.1.0_r4_kind ) then
                        write(*,*) 'expected ', 0.0_r4_kind, ' but computed ', interp%horizInterpReals4_type%wti(i,j,2)
                        call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d wti2")
                    end if
                endif
            end do
        end do
    endif
    if(.not. test_solo) then
        call horiz_interp_del(interp)
        call check_dealloc(interp)
    endif
    !check that data are equal
    do j=1, nj_src
        do i=1, ni_src
            if( data_src(i,j).ne.data_dst(i,j) ) then
                write(*,*) 'expected ', data_src(i,j), ' but computed ', data_dst(i,j)
                call mpp_error(FATAL, "failed at horiz_interp_bilinear_2d2d data comparison")
            end if
        end do
    end do

  end subroutine test_horiz_interp_bilinear

  !> Tests bicubic module interpolation with each dimension conversion
  !! test without passing in the type when test_solo is true
  subroutine test_horiz_interp_bicubic
    !! grid data
    real(HI_TEST_KIND_), allocatable, dimension(:) :: lat_in_1D, lon_in_1D
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: lat_in_2D, lon_in_2D
    type(horiz_interp_type)                       :: interp_t
    !! input data
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: data_src, data_dst
    !! output data
    real(HI_TEST_KIND_), allocatable, dimension(:)   :: lat_out_1D, lon_out_1D
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: lat_out_2D, lon_out_2D
    !! array sizes and number of lat/lon per index
    real(HI_TEST_KIND_) :: nlon_in, nlat_in
    real(HI_TEST_KIND_) :: nlon_out, nlat_out
    real(HI_TEST_KIND_) :: dlon_src, dlat_src, dlon_dst, dlat_dst
    !! parameters for lon/lat setup
    real(HI_TEST_KIND_) :: lon_src_beg = 0._lkind,    lon_src_end = 360._lkind
    real(HI_TEST_KIND_) :: lat_src_beg = -90._lkind,  lat_src_end = 90._lkind
    real(HI_TEST_KIND_) :: lon_dst_beg = -280._lkind, lon_dst_end = 80._lkind
    real(HI_TEST_KIND_) :: lat_dst_beg = -90._lkind,  lat_dst_end = 90._lkind
    real(HI_TEST_KIND_) :: D2R = real(PI,HI_TEST_KIND_)/180._lkind
    real(HI_TEST_KIND_) :: R2D = 180._lkind/real(PI,HI_TEST_KIND_)
    real(HI_TEST_KIND_), parameter :: SMALL = 1.0e-10_lkind

    ! set up longitude and latitude of source/destination grid.
    dlon_src = (lon_src_end-lon_src_beg)/real(ni_src, lkind)
    dlat_src = (lat_src_end-lat_src_beg)/real(nj_src, lkind)
    dlon_dst = (lon_dst_end-lon_dst_beg)/real(ni_dst, lkind)
    dlat_dst = (lat_dst_end-lat_dst_beg)/real(nj_dst, lkind)

    allocate(lon_in_1D(ni_src+1), lat_in_1D(nj_src+1))
    do i = 1, ni_src+1
        lon_in_1D(i) = lon_src_beg + real(i-1, lkind)*dlon_src
    end do
    do j = 1, nj_src+1
        lat_in_1D(j) = lat_src_beg + real(j-1, lkind)*dlat_src
    end do
    allocate(lon_out_1D(isc:iec+1), lat_out_1D(jsc:jec+1))
    do i = isc, iec+1
        lon_out_1D(i) = lon_dst_beg + real(i-1,lkind)*dlon_dst
    end do
    do j = jsc, jec+1
        lat_out_1D(j) = lat_dst_beg + real(j-1,lkind)*dlat_dst
    end do
    ! convert to rads
    lon_in_1D = lon_in_1D * D2R
    lat_in_1D = lat_in_1D * D2R
    lon_out_1D = lon_out_1D * D2R
    lat_out_1D = lat_out_1D * D2R

    ! set up 2d lon/lat
    allocate(lon_out_2D(isc:iec+1, jsc:jec+1), lat_out_2D(isc:iec+1, jsc:jec+1))
    do i = isc, iec+1
        lon_out_2D(i,:) = lon_out_1D(i)
    end do
    do j = jsc, jec+1
        lat_out_2D(:,j) = lat_out_1D(j)
    end do

    nlon_in = real(ni_src, lkind);  nlat_in = real(nj_src, lkind)
    nlon_out = real(iec - isc, lkind); nlat_out = real(jec - jsc, lkind)

    ! allocate data
    allocate(data_src(ni_src, nj_src))
    allocate(data_dst(isc:iec, jsc:jec))
    data_dst = 0.0_lkind ; data_src = 1.0_lkind

    id1 = mpp_clock_id( 'horiz_interp_bicubic_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id2 = mpp_clock_id( 'horiz_interp_bicubic_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! 1D x 1D
    call mpp_clock_begin(id1)
    if(.not. test_solo) then
        call horiz_interp_new(interp_t, lon_in_1d, lat_in_1d, lon_out_1d, lat_out_1d, interp_method="bicubic")
        call horiz_interp(interp_t, data_src, data_dst)
    else
        call horiz_interp(data_src, lon_in_1D, lat_in_1D, lon_out_1D, lat_out_1D, data_dst, interp_method="bicubic")
    endif
    call mpp_clock_end(id1)
    call mpp_sync()
    ! check weights (for last index, 1=x,2=y,3=xy derivatives)
    ! 1 radian (in degrees) at edges, 0.5 otherwise
    if( .not. test_solo) then
        do i=1, ni_src-1
            do j=1, nj_src-1
                if( interp_t%horizInterpReals4_type%is_allocated) then
                    if( interp_t%horizInterpReals4_type%wti(i,j,1) * interp_t%horizInterpReals4_type%wti(i,j,2) &
                        - interp_t%horizInterpReals4_type%wti(i,j,3) .gt. SMALL .or.    &
                        interp_t%horizInterpReals4_type%wti(i,j,3) - (57.2958_lkind * 57.2958_lkind) .gt. SMALL) then
                            print *, i, j, interp_t%horizInterpReals4_type%wti(i,j,:)
                            call mpp_error(FATAL, "test_horiz_interp: bicubic test failed 1Dx1D weight calculation")
                    endif
                else
                    if( interp_t%horizInterpReals8_type%wti(i,j,1) * interp_t%horizInterpReals8_type%wti(i,j,2) &
                        - interp_t%horizInterpReals8_type%wti(i,j,3) .gt. SMALL .and. &
                        interp_t%horizInterpReals8_type%wti(i,j,3) - (57.2958_lkind * 57.2958_lkind) .gt. SMALL) then
                            print *, i, j, interp_t%horizInterpReals8_type%wti(i,j,:)
                            call mpp_error(FATAL, "test_horiz_interp: bicubic test failed 1Dx1D weight calculation")
                    endif
                endif
            enddo
        enddo
        call horiz_interp_del(interp_t)
        call check_dealloc(interp_t)
    endif
    do i=isc, iec
        do j=jsc, jec
            if( data_dst(i,j) .ne. 1.0_lkind) call mpp_error(FATAL, "test_horiz_interp: error in 1Dx1D output data")
        enddo
    enddo

    ! 1D x 2D
    deallocate(data_src, data_dst)
    allocate(data_src(ni_src+1, nj_src+1))
    allocate(data_dst(isc:iec+1, jsc:jec+1))
    data_dst = 0.0_lkind ; data_src = 1.0_lkind

    call mpp_clock_begin(id2)
    if(.not. test_solo) then
        call horiz_interp_new(interp_t, lon_in_1d, lat_in_1d, lon_out_2d, lat_out_2d, interp_method="bicubic")
        call horiz_interp(interp_t, data_src, data_dst)
    else
        call horiz_interp(data_src, lon_in_1D, lat_in_1D, lon_out_2D, lat_out_2D, data_dst, interp_method="bicubic")
    endif
    call mpp_clock_end(id2)
    if( .not. test_solo) then
        do i=1, ni_src-1
            do j=1, nj_src-1
                if( interp_t%horizInterpReals4_type%is_allocated) then
                    if( interp_t%horizInterpReals4_type%wti(i,j,1) * interp_t%horizInterpReals4_type%wti(i,j,2) &
                        - interp_t%horizInterpReals4_type%wti(i,j,3) .gt. SMALL .or.    &
                        interp_t%horizInterpReals4_type%wti(i,j,3) - (57.2958_lkind * 57.2958_lkind) .gt. SMALL) then
                            print *, i, j, interp_t%horizInterpReals4_type%wti(i,j,:)
                            call mpp_error(FATAL, "test_horiz_interp: bicubic test failed 1Dx1D weight calculation")
                    endif
                else
                    if( interp_t%horizInterpReals8_type%wti(i,j,1) * interp_t%horizInterpReals8_type%wti(i,j,2) &
                        - interp_t%horizInterpReals8_type%wti(i,j,3) .gt. SMALL .or. &
                        interp_t%horizInterpReals8_type%wti(i,j,3) - (57.2958_lkind * 57.2958_lkind) .gt. SMALL) then
                            print *, i, j, interp_t%horizInterpReals8_type%wti(i,j,:)
                            call mpp_error(FATAL, "test_horiz_interp: bicubic test failed 1Dx1D weight calculation")
                    endif
                endif
            enddo
        enddo
        call horiz_interp_del(interp_t)
        call check_dealloc(interp_t)
    endif
    do i=isc, iec
        do j=jsc, jec
            if( data_dst(i,j) .ne. 1.0_lkind) call mpp_error(FATAL, "test_horiz_interp: error in 1Dx2D output data")
        enddo
    enddo

    deallocate(data_src, data_dst)
    deallocate(lat_in_1D, lon_in_1D)
    deallocate(lat_out_1D, lon_out_1D, lat_out_2D, lon_out_2D)

  end subroutine test_horiz_interp_bicubic

  !> Tests conservative (default) interpolation module and checks grids reproduce across 1/2d versions
  subroutine test_horiz_interp_conserve
    real(HI_TEST_KIND_)                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
    real(HI_TEST_KIND_), allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: lon2D_src, lat2D_src, lon2D_dst, lat2D_dst
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst, data3_dst, data4_dst
    real(HI_TEST_KIND_), allocatable, dimension(:,:) :: data1_solo, data2_solo, data3_solo, data4_solo
    real(HI_TEST_KIND_) :: lon_src_beg = 0._lkind,    lon_src_end = 360._lkind
    real(HI_TEST_KIND_) :: lat_src_beg = -90._lkind,  lat_src_end = 90._lkind
    real(HI_TEST_KIND_) :: lon_dst_beg = -280._lkind, lon_dst_end = 80._lkind
    real(HI_TEST_KIND_) :: lat_dst_beg = -90._lkind,  lat_dst_end = 90._lkind
    real(HI_TEST_KIND_) :: D2R = real(PI,HI_TEST_KIND_)/180._lkind
    real(HI_TEST_KIND_), parameter :: SMALL = 1.0e-10_lkind
    type(horiz_interp_type)           :: interp_conserve

    allocate(lon2D_src(ni_src+1, nj_src+1), lat2D_src(ni_src+1, nj_src+1) )
    allocate(lon1D_src(ni_src+1), lat1D_src(nj_src+1), data_src(ni_src, nj_src) )

    allocate(lon2D_dst(isc:iec+1, jsc:jec+1), lat2D_dst(isc:iec+1, jsc:jec+1) )
    allocate(lon1D_dst(isc:iec+1), lat1D_dst(jsc:jec+1) )
    allocate(data1_dst(isc:iec, jsc:jec), data2_dst(isc:iec, jsc:jec) )
    allocate(data3_dst(isc:iec, jsc:jec), data4_dst(isc:iec, jsc:jec) )

    ! set up longitude and latitude of source/destination grid.
    dlon_src = (lon_src_end-lon_src_beg)/real(ni_src, lkind)
    dlat_src = (lat_src_end-lat_src_beg)/real(nj_src, lkind)
    dlon_dst = (lon_dst_end-lon_dst_beg)/real(ni_dst, lkind)
    dlat_dst = (lat_dst_end-lat_dst_beg)/real(nj_dst, lkind)

    do i = 1, ni_src+1
        lon1D_src(i) = lon_src_beg + real(i-1, lkind)*dlon_src
    end do

    do j = 1, nj_src+1
        lat1D_src(j) = lat_src_beg + real(j-1, lkind)*dlat_src
    end do

    do i = isc, iec+1
        lon1D_dst(i) = lon_dst_beg + real(i-1, lkind)*dlon_dst
    end do

    do j = jsc, jec+1
        lat1D_dst(j) = lat_dst_beg + real(j-1, lkind)*dlat_dst
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
          data_src(i,j) = real(i,lkind) + real(j,lkind)*0.001_lkind
        end do
    end do

    id1 = mpp_clock_id( 'horiz_interp_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id2 = mpp_clock_id( 'horiz_interp_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id3 = mpp_clock_id( 'horiz_interp_2dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    id4 = mpp_clock_id( 'horiz_interp_2dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

    ! --- 1dx1d version conservative interpolation
    call mpp_clock_begin(id1)
    if(.not. test_solo) then
        call horiz_interp_new(interp_conserve, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, &
                              interp_method = "conservative")
        call horiz_interp(interp_conserve, data_src, data1_dst)
        call horiz_interp_del(interp_conserve)
        call check_dealloc(interp_conserve)
    else
        call horiz_interp(data_src, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, data1_dst, &
                          interp_method="conservative")
    endif
    call mpp_clock_end(id1)

    ! --- 1dx2d version conservative interpolation
    call mpp_clock_begin(id2)
    if(.not. test_solo) then
        call horiz_interp_new(interp_conserve, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, &
                              interp_method = "conservative")
        call horiz_interp(interp_conserve, data_src, data2_dst)
        call horiz_interp_del(interp_conserve)
        call check_dealloc(interp_conserve)
    else
        call horiz_interp(data_src, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, data2_dst, &
                          interp_method="conservative")
    endif
    call mpp_clock_end(id2)

    ! --- 2dx1d version conservative interpolation
    call mpp_clock_begin(id3)
    if(.not. test_solo) then
        call horiz_interp_new(interp_conserve, lon2D_src, lat2D_src, lon1D_dst, lat1D_dst, &
                              interp_method = "conservative")
        call horiz_interp(interp_conserve, data_src, data3_dst)
        call horiz_interp_del(interp_conserve)
        call check_dealloc(interp_conserve)
    else
        call horiz_interp(data_src, lon2D_src, lat2D_src, lon1D_dst, lat1D_dst, data3_dst, &
                          interp_method="conservative")
    endif
    call mpp_clock_end(id3)

    ! --- 2dx2d version conservative interpolation
    call mpp_clock_begin(id4)
    if(.not. test_solo) then
        call horiz_interp_new(interp_conserve, lon2D_src, lat2D_src, lon2D_dst, lat2D_dst, &
                              interp_method = "conservative")
        call horiz_interp(interp_conserve, data_src, data4_dst)
        call horiz_interp_del(interp_conserve)
        call check_dealloc(interp_conserve)
    else
        call horiz_interp(data_src, lon2D_src, lat2D_src, lon2D_dst, lat2D_dst, data4_dst, &
                          interp_method="conservative")
    endif
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

    !> Tests the assignment overload for horiz_interp_type
    !! creates some new instances of the derived type for the different methods
    !! and tests equality of fields after initial weiht calculations
    subroutine test_assignment()
        type(horiz_interp_type) :: Interp_new1, Interp_new2, Interp_cp, intp_3
        !! grid data points
        real(HI_TEST_KIND_), allocatable, dimension(:) :: lat_in_1D, lon_in_1D
        real(HI_TEST_KIND_), allocatable, dimension(:,:) :: lat_in_2D, lon_in_2D
        !! output data points
        real(HI_TEST_KIND_), allocatable, dimension(:)   :: lat_out_1D, lon_out_1D
        real(HI_TEST_KIND_), allocatable, dimension(:,:) :: lat_out_2D, lon_out_2D
        real(HI_TEST_KIND_), allocatable, dimension(:) :: lat_out_bil, lon_out_bil
        real(HI_TEST_KIND_), allocatable, dimension(:,:) :: lat_in_bil, lon_in_bil
        !! array sizes and number of lat/lon per index
        real(HI_TEST_KIND_) :: nlon_in, nlat_in
        real(HI_TEST_KIND_) :: nlon_out, nlat_out
        real(HI_TEST_KIND_) :: dlon_src, dlat_src, dlon_dst, dlat_dst
        !! parameters for lon/lat setup
        real(HI_TEST_KIND_) :: lon_src_beg = 0._lkind,    lon_src_end = 360._lkind
        real(HI_TEST_KIND_) :: lat_src_beg = -90._lkind,  lat_src_end = 90._lkind
        real(HI_TEST_KIND_) :: lon_dst_beg = 0.0_lkind, lon_dst_end = 360._lkind
        real(HI_TEST_KIND_) :: lat_dst_beg = -90._lkind,  lat_dst_end = 90._lkind
        real(HI_TEST_KIND_) :: D2R = real(PI,HI_TEST_KIND_)/180._lkind
        real(HI_TEST_KIND_) :: R2D = 180._lkind/real(PI,HI_TEST_KIND_)
        real(HI_TEST_KIND_), parameter :: SMALL = 1.0e-10_lkind

        ! set up longitude and latitude of source/destination grid.
        dlon_src = (lon_src_end-lon_src_beg)/real(ni_src, lkind)
        dlat_src = (lat_src_end-lat_src_beg)/real(nj_src, lkind)
        dlon_dst = (lon_dst_end-lon_dst_beg)/real(ni_dst, lkind)
        dlat_dst = (lat_dst_end-lat_dst_beg)/real(nj_dst, lkind)

        allocate(lon_in_1D(ni_src+1), lat_in_1D(nj_src+1))
        allocate(lon_out_1D(isc:iec+1), lat_out_1D(jsc:jec+1))
        do i = 1, ni_src+1
            lon_in_1D(i) = lon_src_beg + real(i-1,HI_TEST_KIND_)*dlon_src
        end do
        do j = 1, nj_src+1
            lat_in_1D(j) = lat_src_beg + real(j-1,HI_TEST_KIND_)*dlat_src
        end do
        do i = isc, iec+1
            lon_out_1D(i) = lon_dst_beg + real(i-1,HI_TEST_KIND_)*dlon_dst
        end do
        do j = jsc, jec+1
            lat_out_1D(j) = lat_dst_beg + real(j-1, HI_TEST_KIND_)*dlat_dst
        end do

        lon_in_1D = lon_in_1D * D2R
        lat_in_1D = lat_in_1D * D2R
        lon_out_1D = lon_out_1D * D2R
        lat_out_1D = lat_out_1D * D2R

        allocate(lon_in_2D(ni_src+1, nj_src+1), lat_in_2D(ni_src+1, nj_src+1))
        do i = 1, ni_src+1
            lon_in_2D(i,:) = lon_in_1D(i)
        end do
        do j = 1, nj_src+1
            lat_in_2D(:,j) = lat_in_1D(j)
        end do
        allocate(lon_out_2D(isc:iec+1, jsc:jec+1), lat_out_2D(isc:iec+1, jsc:jec+1))
        do i = isc, iec+1
            lon_out_2D(i,:) = lon_out_1D(i)
        end do
        do j = jsc, jec+1
            lat_out_2D(:,j) = lat_out_1D(j)
        end do

        ! conservative
        ! 1dx1d
        call horiz_interp_new(Interp_new1, lon_in_1D, lat_in_1D, lon_out_1D, lat_out_1D, interp_method="conservative")
        call horiz_interp_new(Interp_new2, lon_in_1D, lat_in_1D, lon_out_1D, lat_out_1D, interp_method="conservative")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 1x1d conservative")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)
        ! 1dx2d
        call horiz_interp_new(Interp_new1, lon_in_1D, lat_in_1D, lon_out_2D, lat_out_2D, interp_method="conservative")
        call horiz_interp_new(Interp_new2, lon_in_1D, lat_in_1D, lon_out_2D, lat_out_2D, interp_method="conservative")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 1x2d conservative")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)
        ! 2dx1d
        call horiz_interp_new(Interp_new1, lon_in_2D, lat_in_2D, lon_out_1D, lat_out_1D, interp_method="conservative")
        call horiz_interp_new(Interp_new2, lon_in_2D, lat_in_2D, lon_out_1D, lat_out_1D, interp_method="conservative")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 2x1d conservative")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)
        ! 2dx2d
        call horiz_interp_new(Interp_new1, lon_in_2D, lat_in_2D, lon_out_2D, lat_out_2D, interp_method="conservative")
        call horiz_interp_new(Interp_new2, lon_in_2D, lat_in_2D, lon_out_2D, lat_out_2D, interp_method="conservative")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 2x2d conservative")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)

        ! bicubic only works with 1d src
        ! 1dx1d
        call horiz_interp_new(Interp_new1, lon_in_1D, lat_in_1D, lon_out_1D, lat_out_1D, interp_method="bicubic")
        call horiz_interp_new(Interp_new2, lon_in_1D, lat_in_1D, lon_out_1D, lat_out_1D, interp_method="bicubic")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 1x1d bicubic")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)
        ! 1dx2d
        call horiz_interp_new(Interp_new1, lon_in_1D, lat_in_1D, lon_out_2D, lat_out_2D, interp_method="bicubic")
        call horiz_interp_new(Interp_new2, lon_in_1D, lat_in_1D, lon_out_2D, lat_out_2D, interp_method="bicubic")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 1x2d bicubic")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)

        deallocate(lon_out_2D, lat_out_2D, lon_in_2D, lat_in_2D)
        allocate(lon_out_2D(ni_dst, nj_dst), lat_out_2D(ni_dst, nj_dst))
        allocate(lon_in_2D(ni_src, nj_src), lat_in_2D(ni_src, nj_src))
        do i = 1, ni_dst
            lon_out_2D(i,:) = lon_dst_beg + real(i-1, HI_TEST_KIND_)*dlon_dst
        end do
        do j = 1, nj_dst
            lat_out_2D(:,j) = lat_dst_beg + real(j-1, HI_TEST_KIND_)*dlat_dst
        end do
        do i = 1, ni_src
            lon_in_2D(i,:) = lon_src_beg + real(i-1, HI_TEST_KIND_)*dlon_src
        end do
        do j = 1, nj_src
            lat_in_2D(:,j) = lat_src_beg + real(j-1, HI_TEST_KIND_)*dlat_src
        end do
        ! scale to radians
        lat_in_2D = lat_in_2D * D2R
        lon_in_2D = lon_in_2D * D2R
        lat_out_2D = lat_out_2D * D2R
        lon_out_2D = lon_out_2D * D2R

        ! spherical
        ! only 2dx2d
        call horiz_interp_new(Interp_new1, lon_in_2D, lat_in_2D, lon_out_2D, lat_out_2D, interp_method="spherical")
        call horiz_interp_new(Interp_new2, lon_in_2D, lat_in_2D, lon_out_2D, lat_out_2D, interp_method="spherical")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 1x2d bilinear")
        call check_type_eq(Interp_cp, Interp_new1)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)

        ! bilinear
        ! 1dx1d
        call horiz_interp_new(Interp_new1, lon_in_1D, lat_in_1D, lon_in_1D, lat_in_1D, interp_method="bilinear")
        call horiz_interp_new(Interp_new2, lon_in_1D, lat_in_1D, lon_in_1D, lat_in_1D, interp_method="bilinear")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 1x1d bilinear")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)
        ! 1dx2d
        call horiz_interp_new(Interp_new1, lon_in_1D, lat_in_1D, lon_in_2D, lat_in_2D, interp_method="bilinear")
        call horiz_interp_new(Interp_new2, lon_in_1D, lat_in_1D, lon_in_2D, lat_in_2D, interp_method="bilinear")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 1x2d bilinear")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)
        ! 2dx1d
        deallocate(lon_out_1D, lat_out_1D)
        allocate(lon_out_1D(ni_dst+1), lat_out_1D(nj_dst+1))
        do i=1, ni_dst
            lon_out_1d(i) = real(i-1, HI_TEST_KIND_) * dlon_dst + lon_dst_beg
        enddo
        do j=1, nj_dst
            lat_out_1d(j) = real(j-1, HI_TEST_KIND_) * dlat_dst + lat_dst_beg
        enddo
        lat_out_1d = lat_out_1D * D2R
        lon_out_1d = lon_out_1D * D2R
        call horiz_interp_new(Interp_new1, lon_in_2D, lat_in_2D, lon_out_1D, lat_out_1D, interp_method="bilinear")
        call horiz_interp_new(Interp_new2, lon_in_2D, lat_in_2D, lon_out_1D, lat_out_1D, interp_method="bilinear")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 1x2d bilinear")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)
        ! 2dx2d
        call horiz_interp_new(Interp_new1, lon_in_2D, lat_in_2D, lon_in_2D, lat_in_2D, interp_method="bilinear")
        call horiz_interp_new(Interp_new2, lon_in_2D, lat_in_2D, lon_in_2D, lat_in_2D, interp_method="bilinear")
        Interp_cp = Interp_new1
        call mpp_error(NOTE,"testing horiz_interp_type assignment 1x2d bilinear")
        call check_type_eq(Interp_cp, Interp_new2)
        call check_type_eq(Interp_cp, Interp_new1)
        call horiz_interp_del(Interp_new1)
        call horiz_interp_del(Interp_new2)
        call horiz_interp_del(Interp_cp)

   end subroutine
    !> helps assignment test with derived type comparisons
    subroutine check_type_eq(interp_1, interp_2)
        type(horiz_interp_type), intent(in) :: interp_1, interp_2
        integer :: k
        if(interp_1%horizInterpReals4_type%is_allocated) then
            if(allocated(interp_1%horizInterpReals4_type%faci)) then
                if( ANY(interp_2%horizInterpReals4_type%faci .ne. interp_1%horizInterpReals4_type%faci)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: faci")
            endif
            if(allocated(interp_1%horizInterpReals4_type%facj)) then
                if( ANY(interp_2%horizInterpReals4_type%facj .ne. interp_1%horizInterpReals4_type%facj)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: facj")
            endif
            if(allocated(interp_1%horizInterpReals4_type%area_src)) then
                if( ANY(interp_2%horizInterpReals4_type%area_src .ne. interp_1%horizInterpReals4_type%area_src)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: area_src")
            endif
            if(allocated(interp_1%horizInterpReals4_type%area_dst)) then
                if( ANY(interp_2%horizInterpReals4_type%area_dst .ne. interp_1%horizInterpReals4_type%area_dst)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: area_dst")
            endif
            if(allocated(interp_1%horizInterpReals4_type%wti)) then
                if( ANY(interp_2%horizInterpReals4_type%wti .ne. interp_1%horizInterpReals4_type%wti)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: wti")
            endif
            if(allocated(interp_1%horizInterpReals4_type%wtj)) then
                if( ANY(interp_2%horizInterpReals4_type%wtj .ne. interp_1%horizInterpReals4_type%wtj)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: wtj")
            endif
            if(allocated(interp_1%horizInterpReals4_type%src_dist)) then
                if( ANY(interp_2%horizInterpReals4_type%src_dist .ne. interp_1%horizInterpReals4_type%src_dist)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: src_dst")
            endif
            if(allocated(interp_1%horizInterpReals4_type%rat_x)) then
                if( ANY(interp_2%horizInterpReals4_type%rat_x .ne. interp_1%horizInterpReals4_type%rat_x)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: src_dst")
            endif
            if(allocated(interp_1%horizInterpReals4_type%rat_y)) then
                if( ANY(interp_2%horizInterpReals4_type%rat_y .ne. interp_1%horizInterpReals4_type%rat_y)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: src_dst")
            endif
            if(allocated(interp_1%horizInterpReals4_type%lon_in)) then
                if( ANY(interp_2%horizInterpReals4_type%lon_in .ne. interp_1%horizInterpReals4_type%lon_in)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: lon_in")
            endif
            if(allocated(interp_1%horizInterpReals4_type%lat_in)) then
                if( ANY(interp_2%horizInterpReals4_type%lat_in .ne. interp_1%horizInterpReals4_type%lat_in)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: lat_in")
            endif

            if(allocated(interp_1%horizInterpReals4_type%area_frac_dst)) then
                if(ANY(interp_2%horizInterpReals4_type%area_frac_dst.ne.interp_1%horizInterpReals4_type%area_frac_dst))&
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: area_frac_dst")
            endif
            if(allocated(interp_1%horizInterpReals4_type%mask_in)) then
                if(ANY(interp_2%horizInterpReals4_type%mask_in.ne.interp_1%horizInterpReals4_type%mask_in)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: mask_in")
            endif
            !! only set during spherical
            if(interp_1%interp_method .eq. SPHERICA) then
                if( interp_2%horizInterpReals4_type%max_src_dist .ne. interp_1%horizInterpReals4_type%max_src_dist) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: max_src_dist")
            endif

        else if(interp_1%horizInterpReals8_type%is_allocated) then
            !!
            if(allocated(interp_1%horizInterpReals8_type%faci)) then
                if( ANY(interp_2%horizInterpReals8_type%faci .ne. interp_1%horizInterpReals8_type%faci)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: faci")
            endif
            if(allocated(interp_1%horizInterpReals8_type%facj)) then
                if( ANY(interp_2%horizInterpReals8_type%facj .ne. interp_1%horizInterpReals8_type%facj)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: facj")
            endif
            if(allocated(interp_1%horizInterpReals8_type%area_src)) then
                if( ANY(interp_2%horizInterpReals8_type%area_src .ne. interp_1%horizInterpReals8_type%area_src)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: area_src")
            endif
            if(allocated(interp_1%horizInterpReals8_type%area_dst)) then
                if( ANY(interp_2%horizInterpReals8_type%area_dst .ne. interp_1%horizInterpReals8_type%area_dst)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: area_dst")
            endif
            if(allocated(interp_1%horizInterpReals8_type%wti)) then
                if( ANY(interp_2%horizInterpReals8_type%wti .ne. interp_1%horizInterpReals8_type%wti)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: wti")
            endif
            if(allocated(interp_1%horizInterpReals8_type%wtj)) then
                if( ANY(interp_2%horizInterpReals8_type%wtj .ne. interp_1%horizInterpReals8_type%wtj)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: wtj")
            endif
            if(allocated(interp_1%horizInterpReals8_type%src_dist)) then
                if( ANY(interp_2%horizInterpReals8_type%src_dist .ne. interp_1%horizInterpReals8_type%src_dist)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: src_dst")
            endif
            if(allocated(interp_1%horizInterpReals8_type%rat_x)) then
                if( ANY(interp_2%horizInterpReals8_type%rat_x .ne. interp_1%horizInterpReals8_type%rat_x)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: src_dst")
            endif
            if(allocated(interp_1%horizInterpReals8_type%rat_y)) then
                if( ANY(interp_2%horizInterpReals8_type%rat_y .ne. interp_1%horizInterpReals8_type%rat_y)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: src_dst")
            endif
            if(allocated(interp_1%horizInterpReals8_type%lon_in)) then
                if( ANY(interp_2%horizInterpReals8_type%lon_in .ne. interp_1%horizInterpReals8_type%lon_in)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: lon_in")
            endif
            if(allocated(interp_1%horizInterpReals8_type%lat_in)) then
                if( ANY(interp_2%horizInterpReals8_type%lat_in .ne. interp_1%horizInterpReals8_type%lat_in)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: lat_in")
            endif

            if(allocated(interp_1%horizInterpReals8_type%area_frac_dst)) then
                if(ANY(interp_2%horizInterpReals8_type%area_frac_dst.ne.interp_1%horizInterpReals8_type%area_frac_dst))&
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: area_frac_dst")
            endif
            if(allocated(interp_1%horizInterpReals8_type%mask_in)) then
                if(ANY(interp_2%horizInterpReals8_type%mask_in.ne.interp_1%horizInterpReals8_type%mask_in)) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: mask_in")
            endif

            !! only set during spherical
            if(interp_1%interp_method .eq. SPHERICA) then
                if( interp_2%horizInterpReals8_type%max_src_dist .ne. interp_1%horizInterpReals8_type%max_src_dist) &
                    call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: max_src_dist")
            endif
        else
            call mpp_error(FATAL, "check_type.ne. both real kinds unallocated")
        endif
        ! non reals
        if(allocated(interp_1%ilon)) then
            if( ANY(interp_2%ilon .ne. interp_1%ilon)) &
                call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: ilon")
        endif
        if(allocated(interp_1%jlat)) then
            if( ANY(interp_2%jlat .ne. interp_1%jlat)) &
                call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: ilon")
        endif
        if(allocated(interp_1%found_neighbors)) then
            if( ANY(interp_2%found_neighbors .neqv. interp_1%found_neighbors)) &
                call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: found_neighbors")
        endif
        if(allocated(interp_1%num_found)) then
            if( ANY(interp_2%num_found .ne. interp_1%num_found)) &
                call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: num_found")
        endif
        if(allocated(interp_1%i_src)) then
            if(ANY(interp_2%i_src .ne. interp_1%i_src)) &
                call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: i_src")
        endif
        if(allocated(interp_1%j_src)) then
            if(ANY(interp_2%j_src .ne. interp_1%j_src)) &
                call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: j_src")
        endif
        if(allocated(interp_1%i_dst)) then
            if(ANY(interp_2%i_dst .ne. interp_1%i_dst)) &
                call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: i_dst")
        endif
        if(allocated(interp_1%j_dst)) then
            if(ANY(interp_2%j_dst .ne. interp_1%j_dst)) &
                call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: j_dst")
        endif
        if(interp_2%nlon_src .ne. interp_1%nlon_src) &
            call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: nlon_src")
        if(interp_2%nlat_src .ne. interp_1%nlat_src) &
            call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: nlat_src")
        if(interp_2%nlon_dst .ne. interp_1%nlon_dst) &
            call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: nlon_dst")
        if(interp_2%nlat_dst .ne. interp_1%nlat_dst) &
                call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: nlat_dst")

        ! these checks were giving me issues with the ALL comparison (gnu), seems to work here tho
        if( allocated(interp_1%i_lon)) then
            do i=1, SIZE(interp_1%i_lon, 1)
                do j=1, SIZE(interp_1%i_lon, 2)
                    do k=1, SIZE(interp_1%i_lon, 3)
                        if(interp_1%i_lon(i,j,k) .ne. interp_2%i_lon(i,j,k)) &
                            call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: i_lon")
                    enddo
                enddo
            enddo
            do i=1, SIZE(interp_1%j_lat, 1)
                do j=1, SIZE(interp_1%j_lat, 2)
                    do k=1, SIZE(interp_1%j_lat, 3)
                        if(interp_1%j_lat(i,j,k) .ne. interp_2%j_lat(i,j,k)) &
                            call mpp_error(FATAL, "Invalid value for copied horiz_interp_type field: j_lat")
                    enddo
                enddo
            enddo
        endif
    end subroutine

    subroutine check_dealloc(hi_type)
        type(horiz_interp_type), intent(in) :: hi_type
        !! can only check the encapsulating real types, inner fields are inaccessible after deallocation
        if(hi_type%horizInterpReals4_type%is_allocated) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: horizInterpReals4_type")
        endif
        if(hi_type%horizInterpReals8_type%is_allocated) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: horizInterpReals8_type")
        endif
        !! non reals
        if(allocated(hi_type%ilon)) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: ilon")
        endif
        if(allocated(hi_type%jlat)) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: jlat")
        endif
        if(allocated(hi_type%found_neighbors)) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: found_neighbors")
        endif
        if(allocated(hi_type%num_found)) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: num_found")
        endif
        if(allocated(hi_type%i_src)) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: i_src")
        endif
        if(allocated(hi_type%j_src)) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: j_src")
        endif
        if(allocated(hi_type%i_dst)) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: i_dst")
        endif
        if(allocated(hi_type%j_dst)) then
            call mpp_error(FATAL, "horiz_interp_test: field left allocated after type deletion: j_dst")
        endif

    end subroutine


end program horiz_interp_test
