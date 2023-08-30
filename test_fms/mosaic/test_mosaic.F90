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

use mosaic2_mod
use grid2_mod
use write_files
use mpp_mod,       only : mpp_init, mpp_error, FATAL, mpp_sync, mpp_npes, mpp_get_current_pelist
use fms2_io_mod,   only : open_file, close_file, FmsNetcdfFile_t, fms2_io_init, read_data
use fms_mod,       only : fms_init, fms_end
use constants_mod, only : DEG_TO_RAD
use platform_mod,  only : r4_kind, r8_kind

implicit none

!> create mosaic and grid files
!! In orderr to create the mosaic and grid files, fms2_io needs to be initialized first
call fms2_io_init()
call write_all()

!< fms_init calls grid_init which reads in the grid_spec file
!! In this case, the grid_version is VERSION_OCN_MOSAIC_FILE.
call fms_init()

write(*,*) 'TEST GET_MOSAIC_GRID_SIZES'
call test_get_mosaic_grid_sizes()

write(*,*) 'TEST GET_MOSAIC_CONTACT'
call test_get_mosaic_contact()

write(*,*) 'TEST GET_GRID_GREAT_CIRCLE_AREA'
call test_get_grid_great_circle_area()

write(*,*) 'TEST GET_GRID_AREA'
call test_get_grid_area()

write(*,*) 'TEST GET_MOSAIC_XGRID'
call test_get_mosaic_xgrid()

write(*,*) 'TEST IS_INSIDE_POLYGON'
call test_is_inside_polygon()

call fms_end()

contains
!------------------------------------------------------!
subroutine test_get_mosaic_grid_sizes

  !> test get_mosaic_grid_sizes

  integer              :: ntiles
  integer              :: n
  integer, allocatable :: nx_out(:), ny_out(:)

  type(FmsNetcdfFile_t):: fileobj

  !-- ocean --!
  if( .not. open_file(fileobj, 'INPUT/'//trim(ocn_mosaic_file), 'read') ) &
       call mpp_error(FATAL, 'test_mosaic: error in opening file '//'INPUT/'//trim(ocn_mosaic_file))

  allocate( nx_out(ocn_ntiles), ny_out(ocn_ntiles) )
  !> get_mosaic_grid_sizes reads in the grid file
  call get_mosaic_grid_sizes(fileobj, nx_out, ny_out )
  do n=1, ocn_ntiles
     call check_answer(ocn_nx/2, nx_out(n), 'ocn TEST_GET_MOSAIC_GRID_SIZES')
     call check_answer(ocn_nY/2, ny_out(n), 'ocn TEST_GET_MOSAIC_GRID_SIZES')
  end do
  deallocate(nx_out, ny_out)
  call close_file(fileobj)

  !-- atm --!
  if( .not. open_file(fileobj, 'INPUT/'//trim(c1_mosaic_file), 'read') ) &
       call mpp_error(FATAL, 'test_mosaic: error in opening file '//'INPUT/'//trim(c1_mosaic_file))

  allocate( nx_out(c1_ntiles), ny_out(c1_ntiles) )
  call get_mosaic_grid_sizes(fileobj, nx_out, ny_out)
  do n=1, ntiles
     call check_answer(c1_nx/2, nx_out(n), 'atm TEST_GET_MOSAIC_GRID_SIZES')
     call check_answer(c1_nx/2, ny_out(n), 'atm TEST_GET_MOSAIC_GRID_SIZES')
  end do
  deallocate(nx_out, ny_out)
  call close_file(fileobj)

end subroutine test_get_mosaic_grid_sizes
!------------------------------------------------------!
subroutine test_get_mosaic_contact

  !< @uriel.ramirez

  integer              :: ntiles         !< Number of tiles
  integer              :: ncontacts      !< Number of contacts
  integer              :: n              !< For do loops
  integer, allocatable :: tile1(:)       !< tile number for first contact
  integer, allocatable :: tile2(:)       !< tile number of the second contact
  integer, allocatable :: nx(:), ny(:)   !< Number of x/y points for each tile
  integer, allocatable :: istart1(:), iend1(:), jstart1(:), jend1(:) !< Indexes of first contact point
  integer, allocatable :: istart2(:), iend2(:), jstart2(:), jend2(:) !< Indexes of second contact point

  integer              :: answers(2, 8)  !< Expected results

  type(FmsNetcdfFile_t):: ocn_fileobj

  if( .not. open_file(ocn_fileobj, 'INPUT/'//trim(ocn_mosaic_file), 'read') ) &
       call mpp_error(FATAL, 'test_mosaic: error in opening file '//'INPUT/'//trim(ocn_mosaic_file))

  answers(1,:) = (/1440, 1440, 1, 1080, 1, 1, 1, 1080 /)
  answers(2,:) = (/1, 720, 1080, 1080, 1440, 721, 1080, 1080 /)

  ntiles = get_mosaic_ntiles(ocn_fileobj)
  ncontacts = get_mosaic_ncontacts(ocn_fileobj)

  allocate(nx(ntiles), ny(ntiles))
  allocate(tile1(ncontacts), tile2(ncontacts) )
  allocate(istart1(ncontacts), iend1(ncontacts), jstart1(ncontacts), jend1(ncontacts) )
  allocate(istart2(ncontacts), iend2(ncontacts), jstart2(ncontacts), jend2(ncontacts) )

  call get_mosaic_contact(ocn_fileobj, tile1, tile2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2)

  !< Compare with expected results:
  if (ntiles .ne. 1)    call mpp_error(FATAL, "ntiles is not equal to 1")
  if (ncontacts .ne. 2) call mpp_error(FATAL, "ncontacts is not the expected result")
  do n = 1, ncontacts
     if (istart1(n) .ne. answers(n,1)) call mpp_error(FATAL, "istart1 is not the expected result")
     if (iend1(n)   .ne. answers(n,2)) call mpp_error(FATAL, "iend1 is not the expected result")

     if (jstart1(n) .ne. answers(n,3)) call mpp_error(FATAL, "jstart1 is not the expected result")
     if (jend1(n)   .ne. answers(n,4)) call mpp_error(FATAL, "jend1 is not the expected result")

     if (istart2(n) .ne. answers(n,5)) call mpp_error(FATAL, "istart2 is not the expected result")
     if (iend2(n)   .ne. answers(n,6)) call mpp_error(FATAL, "iend2 is not the expected result")

     if (jstart2(n) .ne. answers(n,7)) call mpp_error(FATAL, "jstart2 is not the expected result")
     if (jend2(n)   .ne. answers(n,8)) call mpp_error(FATAL, "jend2 is not the expected result")
  end do

  deallocate(tile1, tile2, nx, ny)
  deallocate(istart1, iend1, jstart1, jend1)
  deallocate(istart2, iend2, jstart2, jend2)

end subroutine test_get_mosaic_contact
!------------------------------------------------------!
subroutine test_get_grid_area

  !> This subroutine tests get_grid_area
  !! This subroutine does not check the correctness of the computed area.
  !! Rather, this subroutine only checks the consistency of the areas
  !! given different grid points.  Since each tile is a C1 grid,
  !! the areas of grid cells in tile 1 should equal those of tile 2, etc.
  !! This subroutine only tests calc_mosaic_grid_great_circle_area_2d
  !! For reasons that haven't been figured out, areas of tile 3 do not equal
  !! those of tile 1.  Instead, areas of tile 3 equal those of tile 6.  Tile 3
  !! and tile 6 are unique in that they contain the North and South Poles.  As
  !! to why this should give rise to different areas is to be determined
  !! by FMS in the future.  It may be due to the inaccuracy of the method.

  implicit none

  real :: x_rad(c1_nx, c1_ny), y_rad(c1_nx, c1_ny)
  real :: area_out(1,1)

  !> get answers.  Tile 1 will be the reference/benchmark data
  x_rad = x(1:2, 1:2) * DEG_TO_RAD
  y_rad = y(1:2, 1:2) * DEG_TO_RAD

  call calc_mosaic_grid_area(x_rad, y_rad, area_out)
  call check_answer(area(1,1), area_out(1,1), 'TEST_GET_GRID_AREA')

end subroutine test_get_grid_area
!------------------------------------------------------!
subroutine test_get_grid_great_circle_area

  !> This subroutine tests calc_mosaic_grid_great_circle_area
  !! This subroutine does not check the correctness of the computed area.
  !! Rather, this subroutine only checks the consistency of the areas
  !! given different grid points.  Since each tile is a C1 grid,
  !! the areas of grid cells in tile 1 should equal those of tile 2, etc.
  !! This subroutine only tests calc_mosaic_grid_great_circle_area_2d

  implicit none

  real :: x_rad(c1_nx, c1_ny), y_rad(c1_nx, c1_ny)
  real :: area_out(1,1)

  !> get answers.  Tile 1 will be the reference/benchmark data
  x_rad = x(1:2, 1:2) * DEG_TO_RAD
  y_rad = y(1:2, 1:2) * DEG_TO_RAD
  call calc_mosaic_grid_great_circle_area(x_rad, y_rad, area_out)
  call check_answer(area(1,1), area_out(1,1), 'TEST_GET_GRID_AREA')

end subroutine test_get_grid_great_circle_area
!------------------------------------------------------!
subroutine test_get_mosaic_xgrid

  implicit none

  integer, dimension(ncells) :: i1, j1, i2, j2
  real, dimension(ncells) :: area

  integer :: i
  real :: garea, get_global_area

  type(FmsNetcdfFile_t) x_fileobj

  garea = get_global_area()

  if( .not. open_file(x_fileobj, 'INPUT/'//trim(exchange_file), 'read')) &
       call mpp_error(FATAL, 'test_mosaic: error in opening file '//'INPUT/'//trim(exchange_file))

  call get_mosaic_xgrid(x_fileobj, i1, j1, i2, j2, area)

  !> check answers
  do i=1, ncells
     call check_answer( xgrid_area(i)/garea, area(i),"TEST_GET_MOSAIC_XGRID area")
     call check_answer(tile1_cell(1,i), i1(i), "TEST_GET_MOSAIC_XGRID i1")
     call check_answer(tile2_cell(1,i), i2(i), "TEST_GET_MOSAIC_XGRID i2")
     call check_answer(tile1_cell(1,i), j1(i), "TEST_GET_MOSAIC_XGRID j1")
     call check_answer(tile2_cell(1,i), j2(i), "TEST_GET_MOSAIC_XGRID j2")
  end do

  call close_file(x_fileobj)


end subroutine test_get_mosaic_xgrid
!------------------------------------------------------!
subroutine test_is_inside_polygon

  !> cheating a little.  starting with xyz coordinates (cause easier to understand)

  implicit none

  integer, parameter :: n=5
  integer :: i
  real :: lat1, lon1, x1, y1, z1, r
  real, dimension(n) :: lon2, lat2, x2, y2, z2
  logical :: answer, is_inside

  !> polygon
  x2=0.0
  y2(1)=1.0 ; y2(2)=1.0 ; y2(3)=4.0 ; y2(4)=4.0 ; y2(5)=1.0
  z2(1)=2.0 ; z2(2)=4.0 ; z2(3)=4.0 ; z2(4)=2.0 ; z2(5)=2.0
  do i=1, n
     r = sqrt( x2(i)**2 + y2(i)**2 + z2(i)**2 )
     lon2(i)=atan(y2(i)/x2(i))
     lat2(i)=asin(z2(i)/r)
  end do

  !> point outside of the polygon
  x1=2.0
  y1=5.0
  z1=4.2
  r = sqrt(x1**2+y1**2+z1**2)
  lon1=atan(y1/x1)
  lat1=asin(z1/r)

  answer=.false.
  is_inside=is_inside_polygon(lon1, lat1, lon2, lat2)
  call check_answer(answer,is_inside,' TEST_IS_INSIDE_POLYGON')

  !> point inside the polygon
  x1=0.0
  y1=3.0
  z1=2.5
  r = sqrt(x1**2+y1**2+z1**2)
  lon1=atan(y1/x1)
  lat1=asin(z1/r)

  answer=.true.
  is_inside=is_inside_polygon(lon1, lat1, lon2, lat2)
  call check_answer(answer,is_inside,'TEST_IS_INSIDE_POLYGON')

end subroutine test_is_inside_polygon
!------------------------------------------------------!
subroutine check_answer(answer, myvalue, whoami)

  implicit none
  class(*) :: answer
  class(*) :: myvalue
  character(*) :: whoami

  select type(answer)
  type is ( logical )
     select type(myvalue)
     type is( logical )
        if( answer .neqv. myvalue ) then
           write(*,*) '*************************************'
           write(*,*) 'EXPECTED ', answer, 'but got ', myvalue
           call mpp_error( FATAL,'failed'//trim(whoami) )
        end if
     end select
  type is( real(r4_kind) )
     select type( myvalue)
        type is(real(r4_kind) )
           if( answer .ne. myvalue ) then
              write(*,*) '*************************************'
              write(*,*) 'EXPECTED ', answer, 'but got ', myvalue
              call mpp_error( FATAL,'failed'//trim(whoami) )
           end if
        end select
  type is( real(r8_kind) )
     select type( myvalue)
        type is(real(r4_kind) )
           if( answer .ne. myvalue ) then
              write(*,*) '*************************************'
              write(*,*) 'EXPECTED ', answer, 'but got ', myvalue
              call mpp_error( FATAL,'failed'//trim(whoami) )
           end if
        end select
  end select

end subroutine check_answer
!------------------------------------------------------!
subroutine check_answer_w_tol(answer, myvalue, whoami)

  implicit none

  real, parameter :: tol=1.e-6
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
