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
!! Helper routines for generating input files for test_mosaic2/grid2
!! This portion of code is included in two modules, one for each real precision
  use fms2_io_mod, only: fms2_io_init, open_file, close_file, FmsNetcdfFile_t
  use fms2_io_mod, only: register_axis, register_field, write_data
  use mpp_mod,     only: mpp_init, mpp_sync, mpp_npes, mpp_get_current_pelist
  use fms_mod,     only: fms_init
  use constants_mod, only:  PI, RADIUS
  use platform_mod, only :r4_kind, r8_kind

  implicit none

  character(23), parameter :: grid_spec_file="grid_spec.nc"
  character(23), parameter :: c1_mosaic_file="C1_mosaic.nc"
  character(30), parameter :: ocn_mosaic_file="ocean_mosaic.nc"
  character(50), parameter :: exchange_file="C96_mosaic_tile1Xocean_mosaic_tile1.nc"
  character(30), parameter :: ocn_tile_file="ocean_hgrid.nc"
  character(23), parameter :: tile1_file="C1_grid.tile1.nc"
  character(23), parameter :: tile2_file="C1_grid.tile2.nc"
  character(23), parameter :: tile3_file="C1_grid.tile3.nc"
  character(23), parameter :: tile4_file="C1_grid.tile4.nc"
  character(23), parameter :: tile5_file="C1_grid.tile5.nc"
  character(23), parameter :: tile6_file="C1_grid.tile6.nc"

  ! atm and land
  integer, parameter :: c1_nx=2      !x---x----x
  integer, parameter :: c1_ny=2      !|        |
  integer, parameter :: c1_nxp=3     !x   x    x
  integer, parameter :: c1_nyp=3     !|        |
  integer, parameter :: c1_ntiles=6  !x---x----x
  integer, parameter :: c1_ncontacts=12

  !ocn
  integer, parameter :: ocn_nx=2880
  integer, parameter :: ocn_ny=2160
  integer, parameter :: ocn_ntiles=1
  integer, parameter :: ocn_ncontacts=2

  !exchange
  integer, parameter :: ncells=2

  ! variables for tile1
  character(5) :: tile
  real(TEST_FMS_KIND_), dimension(c1_nxp,c1_nyp) :: x
  real(TEST_FMS_KIND_), dimension(c1_nxp,c1_nyp) :: y
  real(TEST_FMS_KIND_), dimension(c1_nx,c1_ny)   :: area

  !variables for exchange grid cells
  real(TEST_FMS_KIND_), dimension(2,ncells) :: tile1_cell, tile2_cell
  real(TEST_FMS_KIND_), dimension(ncells) :: xgrid_area

  integer, parameter :: lkind=TEST_FMS_KIND_ !< local kind parameter

contains
  !---------------------------------!
  subroutine write_grid_spec

    implicit none
    type(FmsNetcdfFile_t) :: fileobj
    integer, allocatable :: pes(:)

    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)

    if( open_file(fileobj, 'INPUT/'//grid_spec_file, 'overwrite', pelist=pes) ) then
       call register_axis(fileobj, "string", 128)

       call register_field(fileobj, "atm_mosaic_file", "char", dimensions=(/"string"/))
       call register_field(fileobj, "lnd_mosaic_file", "char", dimensions=(/"string"/))
       call register_field(fileobj, "ocn_mosaic_file", "char", dimensions=(/"string"/))

       call write_data(fileobj, "atm_mosaic_file", "C1_mosaic.nc")
       call write_data(fileobj, "lnd_mosaic_file", "C1_mosaic.nc")
       call write_data(fileobj, "ocn_mosaic_file", "ocean_mosaic.nc")

       call close_file(fileobj)
    end if

  end subroutine write_grid_spec
  !---------------------------------!
  subroutine write_c1_mosaic

    implicit none

    type(FmsNetcdfFile_t) :: fileobj
    integer, allocatable :: pes(:)

    character(50), dimension(c1_ntiles) :: strings6
    character(50), dimension(c1_ncontacts) :: strings12

    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)


    if( open_file(fileobj, 'INPUT/'//trim(c1_mosaic_file), 'overwrite', pelist=pes) ) then

       call register_axis(fileobj, 'ntiles', c1_ntiles)
       call register_axis(fileobj, 'ncontact', c1_ncontacts)
       call register_axis(fileobj, 'string', 55)

       call register_field(fileobj, 'mosaic', 'char', dimensions=(/'string'/))
       call register_field(fileobj, 'gridfiles', 'char', dimensions=(/'string','ntiles'/))
       call register_field(fileobj, "gridtiles", "char", dimensions=(/"string","ntiles"/))
       call register_field(fileobj, "contacts", "char", dimensions=(/"string  ","ncontact"/))
       call register_field(fileobj, "contact_index", "char", dimensions=(/"string  ","ncontact"/))

       call write_data(fileobj, "mosaic", "C1_mosaic")

       strings6(1)=tile1_file
       strings6(2)=tile2_file
       strings6(3)=tile3_file
       strings6(4)=tile4_file
       strings6(5)=tile5_file
       strings6(6)=tile6_file
       call write_data(fileobj, "gridfiles", strings6)

       strings6(1)='tile1'
       strings6(2)='tile2'
       strings6(3)='tile3'
       strings6(4)='tile4'
       strings6(5)='tile5'
       strings6(6)='tile6'
       call write_data(fileobj, "gridtiles", strings6)

       strings12(1) ="C1_mosaic:tile1::C1_mosaic:tile2"
       strings12(2) ="C1_mosaic:tile1::C1_mosaic:tile3"
       strings12(3) ="C1_mosaic:tile1::C1_mosaic:tile5"
       strings12(4) ="C1_mosaic:tile1::C1_mosaic:tile6"
       strings12(5) ="C1_mosaic:tile2::C1_mosaic:tile3"
       strings12(6) ="C1_mosaic:tile2::C1_mosaic:tile4"
       strings12(7) ="C1_mosaic:tile2::C1_mosaic:tile6"
       strings12(8) ="C1_mosaic:tile3::C1_mosaic:tile4"
       strings12(9) ="C1_mosaic:tile3::C1_mosaic:tile5"
       strings12(10)="C1_mosaic:tile4::C1_mosaic:tile5"
       strings12(11)="C1_mosaic:tile4::C1_mosaic:tile6"
       strings12(12)="C1_mosaic:tile5::C1_mosaic:tile6"
       call write_data(fileobj, "contacts", strings12)

       strings12(1) ="2:2,1:2::1:1,1:2"
       strings12(2) ="1:2,2:2::1:1,2:1"
       strings12(3) ="1:1,1:2::2:1,2:2"
       strings12(4) ="1:2,1:1::1:2,2:2"
       strings12(5) ="1:2,2:2::1:2,1:1"
       strings12(6) ="2:2,1:2::2:1,1:1"
       strings12(7) ="1:2,1:1::2:2,2:1"
       strings12(8) ="2:2,1:2::1:1,1:2"
       strings12(9) ="1:2,2:2::1:1,2:1"
       strings12(10)="1:2,2:2::1:2,1:1"
       strings12(11)="2:2,1:2::2:1,1:1"
       strings12(12)="2:2,1:2::1:1,1:2"
       call write_data(fileobj, "contact_index", strings12)

       call close_file(fileobj)

    end if

  end subroutine write_c1_mosaic
  !---------------------------------!
  subroutine write_c1_tiles

    !> These are made up numbers, numbers chosen
    !! for computational convenience

    implicit none

    character(5) :: tile
    real(TEST_FMS_KIND_), parameter :: area_value = real(PI*RADIUS*RADIUS/2.0_r8_kind, TEST_FMS_KIND_)
    real(r8_kind) :: xtmp(c1_nxp, c1_nyp), ytmp(c1_nxp, c1_nyp)

    xtmp(1,1)=0.0_r8_kind ; xtmp(2,1)=90.0_r8_kind ; xtmp(3,1)=180.0_r8_kind
    xtmp(1,2)=0.0_r8_kind ; xtmp(2,2)=90.0_r8_kind ; xtmp(3,2)=180.0_r8_kind
    xtmp(1,3)=0.0_r8_kind ; xtmp(2,3)=90.0_r8_kind ; xtmp(3,3)=180.0_r8_kind

    x = real(xtmp,lkind)

    ytmp(1,1)=-90.0_r8_kind ; ytmp(2,1)=-90.0_r8_kind ; ytmp(3,1)=-90.0_r8_kind
    ytmp(1,2)=  0.0_r8_kind ; ytmp(2,2)=  0.0_r8_kind ; ytmp(3,2)=  0.0_r8_kind
    ytmp(1,3)= 90.0_r8_kind ; ytmp(2,3)= 90.0_r8_kind ; ytmp(3,3)= 90.0_r8_kind

    y = real(ytmp,lkind)

    area(1,1)=area_value ; area(2,1)=area_value
    area(1,2)=area_value ; area(2,2)=area_value

    tile='tile1' ; call call_fms2_io(tile1_file, tile, x, y, area)
    tile='tile2' ; call call_fms2_io(tile2_file, tile, x, y, area)
    tile='tile3' ; call call_fms2_io(tile3_file, tile, x, y, area)
    tile='tile4' ; call call_fms2_io(tile4_file, tile, x, y, area)
    tile='tile5' ; call call_fms2_io(tile5_file, tile, x, y, area)
    tile='tile6' ; call call_fms2_io(tile6_file, tile, x, y, area)

  end subroutine write_c1_tiles
 !-----------------------------------!
  subroutine call_fms2_io(filename, tile, x_in, y_in, area_in)

    implicit none

    character(*) :: filename
    character(*) :: tile
    real(TEST_FMS_KIND_), dimension(c1_nxp,c1_nyp), intent(in) :: x_in, y_in
    real(TEST_FMS_KIND_), dimension(c1_nx,c1_ny),   intent(in) :: area_in

    type(FmsNetcdfFile_t) :: fileobj
    integer, allocatable :: pes(:)

    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)

    if( open_file(fileobj, 'INPUT/'//trim(filename), 'overwrite', pelist=pes) ) then

       call register_axis(fileobj, "nx", c1_nx)
       call register_axis(fileobj, "ny", c1_ny)
       call register_axis(fileobj, 'nxp', c1_nxp)
       call register_axis(fileobj, 'nyp', c1_nyp)
       call register_axis(fileobj, "string", 5)

       call register_field(fileobj, 'tile', 'char', dimensions=(/'string'/))
       call register_field(fileobj, 'x', 'double', dimensions=(/'nxp', 'nyp'/))
       call register_field(fileobj, 'y', 'double', dimensions=(/'nxp', 'nyp'/))
       call register_field(fileobj, 'area', 'double', dimensions=(/'nx','ny'/))

       call write_data(fileobj, 'tile', trim(tile))
       call write_data(fileobj, 'x', x_in)
       call write_data(fileobj, 'y', y_in)
       call write_data(fileobj, 'area', area_in)

       call close_file(fileobj)

    end if

  end subroutine call_fms2_io
  !---------------------------------!
  subroutine write_ocean_mosaic()

    !> from @uriel.ramirez

    implicit none

    type(FmsNetcdfFile_t):: fileobj        !< Fileobj for the files written by the test
    integer, allocatable :: pes(:)

    character(38), dimension(ocn_ntiles) :: strings1
    character(38), dimension(ocn_ncontacts) :: strings2

    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)

    if( open_file(fileobj, 'INPUT/'//ocn_mosaic_file, 'overwrite', pelist=pes)) then
       call register_axis(fileobj, "ntiles", ocn_ntiles)
       call register_axis(fileobj, "ncontact", ocn_ncontacts)
       call register_axis(fileobj, "string", 50)

       call register_field(fileobj, "contacts", "char",  dimensions=(/"string  ","ncontact"/))
       call register_field(fileobj, "contact_index", "char",  dimensions=(/"string  ","ncontact"/))
       call register_field(fileobj, "gridfiles", "char", dimensions=(/"string", "ntiles"/))
       call register_field(fileobj, "gridtiles", "char", dimensions=(/"string", "ntiles"/))

       strings1(1)=ocn_tile_file
       call write_data(fileobj, "gridfiles",strings1)

       strings1(1)='tile1'
       call write_data(fileobj, "gridtiles",strings1)

       strings2(1)="2880:2880,1:2160::1:1,1:2160"
       strings2(2)="1:1440,2160:2160::2880:1441,2160:2160"
       call write_data(fileobj, "contact_index", strings2)

       strings2(1)="ocean_mosaic:tile1::ocean_mosaic:tile1"
       strings2(2)="ocean_mosaic:tile1::ocean_mosaic:tile1"
       call write_data(fileobj, "contacts", strings2)

       call close_file(fileobj)
    endif

  end subroutine write_ocean_mosaic
  !----------------------------------
  subroutine write_exchange

    implicit none

    type(FmsNetcdfFile_t):: fileobj        !< Fileobj for the files written by the test
    integer, allocatable :: pes(:)
    integer :: i, j, k
    real(r8_kind) :: get_global_area !< get_global_area returns a double

    !> These are made up numbers, numbers chosen
    !! for computational convenience

    do i=1,ncells
       tile1_cell(1,i) = i
       tile1_cell(2,i) = i
       tile2_cell(1,i) = i
       tile2_cell(2,i) = i
    end do

    do i=1, ncells
       xgrid_area(i) = real(get_global_area(), TEST_FMS_KIND_)
    end do

    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)
    if( open_file(fileobj, 'INPUT/'//trim(exchange_file), "overwrite", pelist=pes)) then
       call register_axis(fileobj, "ncells", ncells)
       call register_axis(fileobj, "two", 2)

       call register_field(fileobj, "tile1_cell", "double", dimensions=(/"two   ", "ncells"/))
       call register_field(fileobj, "tile2_cell", "double", dimensions=(/"two   ", "ncells"/))
       call register_field(fileobj, "xgrid_area", "double", dimensions=(/"ncells"/))

       call write_data(fileobj, "tile1_cell", tile1_cell)
       call write_data(fileobj, "tile2_cell", tile2_cell)
       call write_data(fileobj, "xgrid_area", xgrid_area)

       call close_file(fileobj)
    end if

  end subroutine write_exchange
  !----------------------------------
  subroutine write_hgrid

    !> from @uriel.ramirez

    implicit none

    type(FmsNetcdfFile_t):: fileobj        !< Fileobj for the files written by the test
    integer, allocatable :: pes(:)

    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)

    if( open_file(fileobj, 'INPUT/'//ocn_tile_file, "overwrite", pelist=pes)) then
       call register_axis(fileobj, "nx", ocn_nx)
       call register_axis(fileobj, "ny", ocn_ny)

       call close_file(fileobj)
    endif

  end subroutine write_hgrid
  !---------------------------------!
  subroutine write_all

    implicit none
    call write_grid_spec()
    call write_c1_mosaic()
    call write_ocean_mosaic()
    call write_c1_tiles()
    call write_hgrid()
    call write_exchange()

  end subroutine write_all
