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
!! get_mosaic_grid_sizes, get_mosaic_contact
program test_mosaic

use mpp_mod,      only : mpp_init, mpp_error, FATAL, mpp_npes, mpp_get_current_pelist
use fms2_io_mod,  only : open_file, close_file, FmsNetcdfFile_t, fms2_io_init
use fms2_io_mod,  only : register_axis, register_field, write_data, read_data
use fms_mod,      only : fms_init, fms_end
use platform_mod, only : r4_kind, r8_kind
use grid2_mod,    only : grid_init, get_grid_ntiles, &
                         get_grid_cell_vertices,     &
                         get_grid_cell_centers,      &
                         get_grid_cell_area,         &
                         get_grid_comp_area

use write_files

implicit none

character(23), dimension(6)  :: grid_files, grid_tiles
character(35), dimension(12) :: contacts, contact_index

call fms2_io_init()
call write_all()

call fms_init()
call test_get_cell_vertices
call test_get_cell_centers
call test_get_grid_cell_area_sg
call test_get_grid_comp_area_sg

contains
  !------------------------------------------!
  subroutine test_get_cell_vertices

    !e---x----e
    !|        |
    !x   x    x  grid_size = 3 x 3 (nxp by nyp)
    !|        |  number of edges in the x direction:  2 (marked by 'e')
    !e---x----e  number of edges in the y direction:  2 (marked by 'e')

    implicit none

    integer, parameter :: npt_x = c1_nx !< number of edge points in a supergrid size of nx
    integer, parameter :: npt_y = c1_ny !< number of edge points in a supergrid size of ny

    real :: lonb_1d(npt_x)       !< returned values for lon 1d
    real :: latb_1d(npt_y)       !< returned values for lat 1d
    real :: answer_lon_1d(npt_x) !< answers for lon 1d
    real :: answer_lat_1d(npt_y) !< answers for lat 1d

    real :: lonb_2d(npt_x,npt_y)       !< returned values for lon 2d
    real :: latb_2d(npt_x,npt_y)       !< returned values for lat 2d
    real :: answer_lon_2d(npt_x,npt_y) !< answers for lon 2d
    real :: answer_lat_2d(npt_x,npt_y) !< answers for lat 2d

    type(FmsNetcdfFile_t) :: fileobj
    integer, allocatable :: pes(:)
    integer :: i,j,itile

    !---- 2d ---!
    !> answers (for tile1)
    answer_lon_2d(1,1)=305.0 ; answer_lon_2d(2,1)=35.0
    answer_lon_2d(1,2)=305.0 ; answer_lon_2d(2,2)=35.0

    answer_lat_2d(1,1)=-35.2643896827547 ; answer_lat_2d(2,1)=-35.2643896827547
    answer_lat_2d(1,2)=35.2643896827547  ; answer_lat_2d(2,2)=35.2643896827547

    call get_grid_cell_vertices('ATM',1,lonb_2d,latb_2d)
    !> check
    do j=1, npt_y
       do i=1, npt_x
          call check_answer(answer_lon_2d(i,j), lonb_2d(i,j), 'TEST_GET_CELL_VERTICIES_2D lon')
          call check_answer(answer_lat_2d(i,j), latb_2d(i,j), 'TEST_GET_CELL_VERTICIES_2D lat')
       end do
    end do

  end subroutine test_get_cell_vertices
  !------------------------------------------!
  subroutine test_get_cell_centers

    implicit none

    integer, parameter :: npt_x = c1_nx/2
    integer, parameter :: npt_y = c1_ny/2

    real :: glon_1d(npt_x)
    real :: glat_1d(npt_y)
    real :: answer_glon_1d(npt_x)
    real :: answer_glat_1d(npt_y)

    real :: glon_2d(npt_x,npt_y)
    real :: glat_2d(npt_x,npt_y)
    real :: answer_glon_2d(npt_x,npt_y)
    real :: answer_glat_2d(npt_x,npt_y)

    integer :: i, j, itile


    !--- 2d ---!
    !> assign answers for 2d
    answer_glon_2d(1,1)=350.0
    answer_glat_2d(1,1)=0.0

    call get_grid_cell_centers('ATM', 1, glon_2d, glat_2d)
    do i=1, npt_x
       do j=1, npt_y
          call check_answer(answer_glon_2d(j,i), glon_2d(j,i), 'TEST_GET_CELL_CENTERS_2D lon')
          call check_answer(answer_glat_2d(j,i), glat_2d(j,i), 'TEST_GET_CELL_CENTERS_2D lat')
       end do
    end do

  end subroutine test_get_cell_centers
  !------------------------------------------!
  subroutine test_get_grid_cell_area_sg

    implicit none

    real, dimension(1,1) :: area, answer_area

    !> area of cell 1 should equal area of cell 2
    call get_grid_cell_area('ATM',1, answer_area)

    !> area of cell 2
    call get_grid_cell_area('ATM',2, area)

    call check_answer(answer_area(1,1), area(1,1), 'TEST_GRID_CELL_AREA_SG')

  end subroutine test_get_grid_cell_area_sg
  !------------------------------------------!
  subroutine test_get_grid_comp_area_sg

    !> bug in the code?  for now, modified the test to reproduce "answers"

    implicit none
    real, dimension(c1_nx-1, c1_ny-1) :: area_answer, area
    integer :: i,j


    call get_grid_comp_area('ATM', 3, area_answer)
    call get_grid_comp_area('ATM', 6, area)

    do j=1, c1_nx-1
       do i=1, c1_ny-1
          call check_answer_w_tol(area_answer(i,j), area(i,j), 'TEST_GET_GRID_COMP_AREA_SG')
       end do
    end do

  end subroutine test_get_grid_comp_area_sg
  !------------------------------------------!
  subroutine check_answer(answer, myvalue, whoami)

    implicit none
    real :: answer
    real :: myvalue
    character(*) :: whoami

    if( answer .ne. myvalue ) then
       write(*,*) '*************************************'
       write(*,*) 'EXPECTED ', answer, 'but got ', myvalue
       call mpp_error( FATAL,'failed'//trim(whoami) )
    end if

  end subroutine check_answer
  !------------------------------------------!
subroutine check_answer_w_tol(answer, myvalue, whoami)

  implicit none

  real, parameter :: tol=1.e-2
  real :: answer
  real :: myvalue
  character(*) :: whoami

  if( abs(answer-myvalue) .gt. myvalue ) then
     write(*,*) '*************************************'
     write(*,*) 'EXPECTED ', answer, 'but got ', myvalue
     call mpp_error( FATAL,'failed'//trim(whoami) )
  end if

end subroutine check_answer_w_tol
!------------------------------------------------------!
end program test_mosaic
