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
!> @brief  This programs tests calls to get_mosaic_ntiles, get_mosaic_ncontacts,
!! get_mosaic_grid_sizes, get_mosaic_contact.  All subroutines here are tested
!! with C1 tiles where tiles 1-6 are identical.  The tile points are made up with
!! values that result in simple answers.  See write_files module for grid details.
program test_mosaic

use mpp_mod,      only : mpp_init, mpp_error, FATAL, mpp_npes, mpp_pe, mpp_root_pe
use mpp_domains_mod, only: domain2D, domainUG, mpp_define_domains, mpp_get_compute_domain, mpp_define_unstruct_domain
use fms2_io_mod,  only : open_file, close_file, FmsNetcdfFile_t, fms2_io_init
use fms2_io_mod,  only : register_axis, register_field, write_data, read_data
use fms_mod,      only : fms_init, fms_end
use platform_mod, only : r4_kind, r8_kind
use grid2_mod
use WRITE_FILES_MOD_

implicit none

!> write out netcdf files
!! write_all sets up the grids
call fms2_io_init()
call write_all()
call fms_init()

if(mpp_pe() .eq. mpp_root_pe()) write(*,*) 'TEST GET_CELL_VERTICIES'
call test_get_cell_vertices

if(mpp_pe() .eq. mpp_root_pe()) write(*,*) 'TEST GET_CELL_CENTERS'
call test_get_cell_centers

if(mpp_pe() .eq. mpp_root_pe()) write(*,*) 'TEST GET_GRID_CELL_AREA_SG'
call test_get_grid_cell_area_sg

if(mpp_pe() .eq. mpp_root_pe()) write(*,*) 'TEST_GET_GRID_CELL_AREA_UG'
call test_get_grid_cell_area_ug

if(mpp_pe() .eq. mpp_root_pe()) write(*,*) 'TEST GET_GRID_COMP_AREA_SG'
call test_get_grid_comp_area_sg

if(mpp_pe() .eq. mpp_root_pe()) write(*,*) 'TEST GET_GRID_COMP_AREA_UG'
call test_get_grid_comp_area_ug

call fms_end()

contains
  !------------------------------------------!
  subroutine test_get_cell_vertices

    !> This subroutine tests get_cell_verticees. This
    !! subroutine only tests for vertices in tile 1.

    implicit none

    real(TEST_FMS_KIND_) :: lonb_2d(c1_nx,c1_ny)       !< returned values for lon 2d
    real(TEST_FMS_KIND_) :: latb_2d(c1_nx,c1_ny)       !< returned values for lat 2d
    real(TEST_FMS_KIND_) :: answer_lon_2d(c1_nx,c1_ny) !< answers for lon 2d
    real(TEST_FMS_KIND_) :: answer_lat_2d(c1_nx,c1_ny) !< answers for lat 2d

    integer :: i,j

    !> answers
    answer_lon_2d=x(1:c1_nxp:2, 1:c1_nxp:2)
    answer_lat_2d=y(1:c1_nxp:2, 1:c1_nxp:2)

    call get_grid_cell_vertices('ATM',1,lonb_2d,latb_2d)
    !> check
    do j=1, c1_ny
       do i=1, c1_nx
          call check_answer(answer_lon_2d(i,j), lonb_2d(i,j), 'TEST_GET_CELL_VERTICIES_2D lon')
          call check_answer(answer_lat_2d(i,j), latb_2d(i,j), 'TEST_GET_CELL_VERTICIES_2D lat')
       end do
    end do

  end subroutine test_get_cell_vertices
  !------------------------------------------!
  subroutine test_get_cell_centers

    !> This subroutine tests get_cell_centers.
    !! There is only one cell center point in a C1 tile.

    implicit none

    integer, parameter :: nx = c1_nx/2 !< number of center points
    integer, parameter :: ny = c1_ny/2 !< number of center points

    real(TEST_FMS_KIND_) :: glon_2d(nx,ny) !< results from grid_cell_centers
    real(TEST_FMS_KIND_) :: glat_2d(nx,ny) !< results from grid_cell_centers
    real(TEST_FMS_KIND_) :: answer_glon_2d(nx,ny) !< answers for glon
    real(TEST_FMS_KIND_) :: answer_glat_2d(nx,ny) !< answers for glat

    integer :: i, j

    !--- 2d ---!
    answer_glon_2d=x(2:c1_nx:2, 2:c1_nx:2)
    answer_glat_2d=y(2:c1_nx:2, 2:c1_nx:2)

    call get_grid_cell_centers('ATM', 1, glon_2d, glat_2d)
    do i=1, nx
       do j=1, ny
          call check_answer(answer_glon_2d(j,i), glon_2d(j,i), 'TEST_GET_CELL_CENTERS_2D lon')
          call check_answer(answer_glat_2d(j,i), glat_2d(j,i), 'TEST_GET_CELL_CENTERS_2D lat')
       end do
    end do

  end subroutine test_get_cell_centers
  !------------------------------------------!
  subroutine test_get_grid_cell_area_sg

    !> This subroutine tests get_grid_cell_area_SG
    !! first without the domain input argument and second
    !! with the domain input argument

    implicit none

    type(domain2D) :: SG_domain
    real(TEST_FMS_KIND_) :: area_out2(1,1)
    real(TEST_FMS_KIND_) :: answer

    answer = real(2.0_r8_kind*PI*RADIUS*RADIUS,lkind)

    !> total of 1 domain with 1 (center) point in the domain
    call mpp_define_domains((/1,1,1,1/), (/1,1/), SG_domain)

    !> The area computed by get_grid_cell_area is for the entire cell
    !! The array area, set in write_files.F90, is the area for 1/4th of the cell

    !> Test withtout SG_domain
    call get_grid_cell_area('ATM',2, area_out2)
    call check_answer(answer, area_out2(1,1), 'TEST_GRID_CELL_AREA_SG')

    !> Test with SG_domain
    call get_grid_cell_area('ATM',2, area_out2, SG_domain)
    call check_answer(answer, area_out2(1,1), 'TEST_GRID_CELL_AREA_SG with SG_domain')

  end subroutine test_get_grid_cell_area_sg
  !------------------------------------------!
  subroutine test_get_grid_cell_area_ug

    !> This subroutine tests get_grid_cell_area_ug

    implicit none
    type(domain2D) :: SG_domain
    type(domainUG) :: UG_domain !< UG_domain is the same as SG_domain
    real(TEST_FMS_KIND_) :: area_out1(1)
    real(TEST_FMS_KIND_) :: answer
    integer :: i
    integer :: npts_tile(1),grid_nlevel(1), ndivs, grid_index(1)

    npts_tile=1
    grid_nlevel=1
    ndivs=1
    grid_index=1

    answer = real( 4.0_r8_kind * area(1,1), TEST_FMS_KIND_)

    !> The unstructured grid is the same as the structured grid; there's only one center point in the tile.
    call mpp_define_domains((/1,1,1,1/), (/1,1/), SG_domain)
    call mpp_define_unstruct_domain(UG_domain, SG_domain,npts_tile,grid_nlevel,&
                                    mpp_npes(),ndivs,grid_index,name='immadeup')

    !> The area computed by get_grid_cell_area is for the entire cell
    !! The array area, set in write_files.F90, is the area for 1/4th of the cell
    call get_grid_cell_area('ATM',1, area_out1, SG_domain, UG_domain)
    call check_answer(answer, area_out1(1), 'TEST_GRID_CELL_AREA_UG')

  end subroutine test_get_grid_cell_area_ug
  !------------------------------------------!
  subroutine test_get_grid_comp_area_sg

    !> This subroutine tests get_grid_comp_area_sg
    !! first without the domain input argument and second
    !! with the domain input argument

    implicit none
    type(domain2D) :: SG_domain
    real(TEST_FMS_KIND_) :: area_out2(1,1)
    real(TEST_FMS_KIND_) :: answer

    answer = real( 4.0_r8_kind * area(1,1), TEST_FMS_KIND_)

    call mpp_define_domains((/1,1,1,1/), (/1,1/), SG_domain)

    !> The area computed by get_grid_cell_area is for the entire cell
    !! The array area, set in write_files.F90, is the area for 1/4th of the cell
    !! Test without SG_domain
    call get_grid_comp_area('ATM', 1, area_out2)
    call check_answer(answer, area_out2(1,1), 'TEST_GRID_COMP_AREA_SG')

    !> The area computed by get_grid_cell_area is for the entire cell
    !! The array area, set in write_files.F90, is the area for 1/4th of the cell
    !! Test with SG_domain
    call get_grid_comp_area('ATM', 1, area_out2, SG_domain)
    call check_answer(answer, area_out2(1,1), 'TEST_GRID_COMP_AREA_SG with SG_domain')

  end subroutine test_get_grid_comp_area_sg
  !------------------------------------------!
  subroutine test_get_grid_comp_area_ug

    !> This subroutine tests get_grid_comp_area_ug

    implicit none
    type(domain2D) :: SG_domain
    type(domainUG) :: UG_domain !< UG_domain is the same as SG_domain
    integer :: npts_tile(1), ntiles_grid(1), grid_index(1)
    real(TEST_FMS_KIND_) :: answer
    real(TEST_FMS_KIND_) :: area_out1(1)

    npts_tile=1
    ntiles_grid=1
    grid_index(1)=1
    answer = real( 4.0_r8_kind * area(1,1), TEST_FMS_KIND_)

    !> the unstructured grid is the same as the structured grid
    call mpp_define_domains((/1,1,1,1/), (/1,1/), SG_domain)
    call mpp_define_unstruct_domain(UG_domain, SG_domain, npts_tile, ntiles_grid,mpp_npes(),1,grid_index)

    !> The area computed by get_grid_cell_area is for the entire cell
    !! The array area, set in write_files.F90, is the area for 1/4th of the cell
    call get_grid_comp_area('ATM',3,area_out1,SG_domain, UG_domain)
    call check_answer(answer, area_out1(1), 'TEST_GRID_CELL_AREA_UG')

  end subroutine test_get_grid_comp_area_ug
  !------------------------------------------!
  subroutine check_answer(answer, myvalue, whoami)

    implicit none
    real(TEST_FMS_KIND_) :: answer
    real(TEST_FMS_KIND_) :: myvalue
    character(*) :: whoami

    if( answer .ne. myvalue ) then
       write(*,*) '*************************************'
       write(*,*) 'EXPECTED ', answer, 'but got ', myvalue
       write(*,*) 'difference of', abs(answer-myvalue)
       call mpp_error( FATAL,'failed '//trim(whoami) )
    end if

  end subroutine check_answer
!------------------------------------------------------!
end program test_mosaic
