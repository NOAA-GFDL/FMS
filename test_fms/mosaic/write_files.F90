module write_files

  use fms2_io_mod, only: fms2_io_init, open_file, close_file, FmsNetcdfFile_t
  use fms2_io_mod, only: register_axis, register_field, write_data
  use mpp_mod,     only: mpp_init, mpp_sync, mpp_npes, mpp_get_current_pelist
  use fms_mod,     only: fms_init
  use constants_mod, only:  PI

  implicit none

  character(23), parameter :: grid_spec_file="grid_spec.nc"
  character(23), parameter :: c1_mosaic_file="C1_mosaic.nc"
  character(30), parameter :: ocn_mosaic_file="ocean_mosaic.nc"
  character(30), parameter :: ocn_tile_file="ocean_hgrid.nc"
  character(50), parameter :: exchange_file="C96_mosaic_tile1Xocean_mosaic_tile1.nc"
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
  integer, parameter :: ncells=10

  ! variables for tile1
  character(5) :: tile
  real, dimension(c1_nxp,c1_nyp) :: x1, x2, x3, x4, x5, x6
  real, dimension(c1_nxp,c1_nyp) :: y1, y2, y3, y4, y5, y6
  real, dimension(c1_nxp,c1_nyp) :: angle_dx1, angle_dx2, angle_dx3, angle_dx4, angle_dx5, angle_dx6
  real, dimension(c1_nxp,c1_nyp) :: angle_dy1, angle_dy2, angle_dy3, angle_dy4, angle_dy5, angle_dy6
  real, dimension(c1_nx,c1_ny)   :: area1, area2, area3, area4, area5, area6
  real, dimension(c1_nx,c1_nyp)  :: dx1, dx2, dx3, dx4, dx5, dx6
  real, dimension(c1_nxp,c1_ny)  :: dy1, dy2, dy3, dy4, dy5, dy6

  !variables for exchange grid cells
  real, dimension(2,ncells) :: tile1_cell, tile2_cell
  real, dimension(ncells) :: xgrid_area

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
  subroutine write_c1_tile1

    implicit none

    character(5) :: tile

    tile='tile1'

    x1(1,1)=305.0 ; x1(2,1)=350.0 ; x1(3,1)=35.0
    x1(1,2)=305.0 ; x1(2,2)=350.0 ; x1(3,2)=35.0
    x1(1,3)=305.0 ; x1(2,3)=350.0 ; x1(3,3)=35.0

    y1(1,1)=-35.2643896827547 ; y1(2,1)=-45.0 ; y1(3,1)=-35.2643896827547
    y1(1,2)=  0.0000000000000 ; y1(2,2)=  0.0 ; y1(3,2)=  0.0000000000000
    y1(1,3)= 35.2643896827547 ; y1(2,3)= 45.0 ; y1(3,3)= 35.2643896827547

    dx1(1,1)=3921221.22393904 ; dx1(2,1)=3921221.22393904
    dx1(1,2)=5003771.69900514 ; dx1(2,2)=5003771.69900514
    dx1(1,3)=3921221.22393904 ; dx1(2,3)=3921221.22393904

    dy1(1,1)=3921221.22393904 ; dy1(2,1)=5003771.69900514 ; dy1(3,1)=3921221.22393904
    dy1(1,2)=3921221.22393904 ; dy1(2,2)=5003771.69900514 ; dy1(3,2)=3921221.22393904

    area1(1,1)=21252686329574.5 ; area1(2,1)=21252686329574.5
    area1(1,2)=21252686329574.6 ; area1(2,2)=21252686329574.5

    angle_dx1(1,1)=-12.2083047812356 ; angle_dx1(2,1)=180.0 ; angle_dx1(3,1)= 178.229638007655
    angle_dx1(1,2)= 0.0              ; angle_dx1(2,2)=180.0 ; angle_dx1(3,2)= 180.0
    angle_dx1(1,3)= 12.2083047812356 ; angle_dx1(2,3)=180.0 ; angle_dx1(3,3)=-178.229638007655

    angle_dy1 = 90.0

    call call_fms2_io(tile1_file, tile, x1, y1, dx1, dy1, area1, angle_dx1, angle_dy1)

  end subroutine write_c1_tile1
  !-----------------------------------!
  subroutine write_c1_tile2

    implicit none

    character(5) :: tile

    tile='tile2'

    x2(1,1)=35.0 ; x2(2,1)=80.0 ; x2(3,1)=125.0
    x2(1,2)=35.0 ; x2(2,2)=80.0 ; x2(3,2)=125.0
    x2(1,3)=35.0 ; x2(2,3)=80.0 ; x2(3,3)=125.0

    y2(1,1)=-35.2643896827547 ; y2(2,1)=-45.0 ;                ; y2(3,1)=-35.2643896827547
    y2(1,2)=  0.0             ; y2(2,2)=  2.75444115227293e-15 ; y2(3,2)=  3.89536803430295e-15
    y2(1,3)= 35.2643896827547 ; y2(2,3)= 45.0                  ; y2(3,3)= 35.2643896827547

    dx2(1,1)=3921221.22393904 ; dx2(2,1)=3921221.22393904
    dx2(1,2)=5003771.69900514 ; dx2(2,2)=5003771.69900514
    dx2(1,3)=3921221.22393904 ; dx2(2,3)=3921221.22393904

    dy2(1,1)=3921221.22393904 ; dy2(2,1)=5003771.69900514 ; dy2(3,1)=3921221.22393904
    dy2(1,2)=3921221.22393904 ; dy2(2,2)=5003771.69900514 ; dy2(3,2)=3921221.22393904

    area2(1,1)=21252686329574.5 ; area2(2,1)=21252686329574.5
    area2(1,2)=21252686329574.6 ; area2(2,2)= 21252686329574.5

    angle_dx2(1,1)=-12.2083047812356 ; angle_dx2(2,1)=0.000000000000000000    ; angle_dx2(3,1)=45.0016529252264
    angle_dx2(1,2)= 3.50706339871978e-15 ; angle_dx2(2,2)=2.4798683112859e-15 ; angle_dx2(3,2)=38.0841344663695
    angle_dx2(1,3)= 12.2083047812356 ; angle_dx2(2,3)=4.04998849185443e-15    ; angle_dx2(3,3)=-4.12501788256499

    angle_dy2 = 90.0

    call call_fms2_io(tile2_file, tile, x2, y2, dx2, dy2, area2, angle_dx2, angle_dy2)

  end subroutine write_c1_tile2
  !-----------------------------------!
  subroutine write_c1_tile3

    implicit none

    character(5) :: tile

    tile='tile3'

    x3(1,1)= 35.0 ; x3(2,1)=80.0  ; x3(3,1)=125.0
    x3(1,2)=350.0 ; x3(2,2)=0.0   ; x3(3,2)=170.0
    x3(1,3)=305.0 ; x3(2,3)=260.0 ; x3(3,3)=215.0

    y3(1,1)=35.2643896827547 ; y3(2,1)=45.0 ; y3(3,1)=35.2643896827547
    y3(1,2)=45.0             ; y3(2,2)=90.0 ; y3(3,2)=45.0
    y3(1,3)=35.2643896827547 ; y3(2,3)=45.0 ; y3(3,3)=35.2643896827547

    dx3(1,1)=3921221.22393904 ; dx3(2,1)=3921221.22393904
    dx3(1,2)=5003771.69900514 ; dx3(2,2)=5003771.69900514
    dx3(1,3)=3921221.22393904 ; dx3(2,3)=3921221.22393904

    dy3(1,1)=3921221.22393904 ; dy3(2,1)=5003771.69900514 ; dy3(3,1)=3921221.22393904
    dy3(1,2)=3921221.22393904 ; dy3(2,2)=5003771.69900514 ; dy3(3,2)=3921221.22393904

    area3(1,1)=21252686329574.5 ; area3(2,1)=21252686329574.5
    area3(1,2)=21252686329574.6 ; area3(2,2)=21252686329574.5

    angle_dx3(1,1)=12.2083047812356 ; angle_dx3(2,1)=4.04998849185443e-15 ; angle_dx3(3,1)=-12.2083047812356
    angle_dx3(1,2)=172.67291243803  ; angle_dx3(2,2)=180.000000000000     ; angle_dx3(3,2)=-14.827811681987
    angle_dx3(1,3)=167.791695218764 ; angle_dx3(2,3)=-180.00000000000     ; angle_dx3(3,3)=-167.791695218764

    angle_dy3(1,1)=1.77036199234526 ; angle_dy3(2,1)=150.639946160165 ; angle_dy3(3,1)=12.2083047812356
    angle_dy3(1,2)=0            ; angle_dy3(2,2)=-4.05113034111413e-15; angle_dy3(3,2)=-1.61999539674177e-14
    angle_dy3(1,3)=-167.791695218764; angle_dy3(2,3)=-9.82020528879405; angle_dy3(3,3)=-12.2083047812356

    call call_fms2_io(tile3_file, tile, x3, y3, dx3, dy3, area3, angle_dx3, angle_dy3)

  end subroutine write_c1_tile3
  !-----------------------------------!
  subroutine write_c1_tile4

    implicit none

    character(5) :: tile

    tile='tile4'

    x4(1,1)=125.0 ; x4(2,1)=125.0 ; x4(3,1)=125.0
    x4(1,2)=170.0 ; x4(2,2)=170.0 ; x4(3,2)=170.0
    x4(1,3)=215.0 ; x4(2,3)=215.0 ; x4(3,3)=215.0

    y4(1,1)=35.2643896827547 ; y4(2,1)=3.89536803430295e-15 ; y4(3,1)=-35.2643896827547
    y4(1,2)=45.0000000000000 ; y4(2,2)=8.26332345681879e-15 ; y4(3,2)=-45.0000000000000
    y4(1,3)=35.2643896827547 ; y4(2,3)=7.7907360686059e-15  ; y4(3,3)=-35.2643896827546

    dx4(1,1)=3921221.22393904 ; dx4(2,1)=3921221.22393904
    dx4(1,2)=5003771.69900515 ; dx4(2,2)=5003771.69900515
    dx4(1,3)=3921221.22393904 ; dx4(2,3)=3921221.22393904

    dy4(1,1)=3921221.22393904 ; dy4(2,1)=5003771.69900515 ; dy4(3,1)=3921221.22393904
    dy4(1,2)=3921221.22393904 ; dy4(2,2)=5003771.69900515 ; dy4(3,2)=3921221.22393904

    area4(1,1)=21252686329574.5 ; area4(2,1)=21252686329574.5
    area4(1,2)=21252686329574.6 ; area4(2,2)=21252686329574.5

    angle_dx4(1,1)=-90.0 ; angle_dx4(2,1)=-90.0 ; angle_dx4(3,1)=-45.0016529252264
    angle_dx4(1,2)=-90.0 ; angle_dx4(2,2)=-90.0 ; angle_dx4(3,2)=-38.0867480461066
    angle_dx4(1,3)=-90.0 ; angle_dx4(2,3)=-90.0 ; angle_dx4(3,3)=-21.3976877148648

    angle_dy4(1,1)=12.2083047812356 ; angle_dy4(2,1)=5.56145357358754e-15     ; angle_dy4(3,1)=-12.2083047812356
    angle_dy4(1,2)=-1.61999539674177e-14 ; angle_dy4(2,2)=2.4798683112859e-15 ; angle_dy4(3,2)=2.02499424592721e-14
    angle_dy4(1,3)=-12.2083047812356 ; angle_dy4(2,3)=-6.01716951015755e-16   ; angle_dy4(3,3)=12.2083047812356

    call call_fms2_io(tile4_file, tile, x4, y4, dx4, dy4, area4, angle_dx4, angle_dy4)

  end subroutine write_c1_tile4
  !-----------------------------------!
  subroutine write_c1_tile5

    implicit none

    character(5) :: tile

    tile='tile5'

    x5(1,1)=215.0 ; x5(2,1)=215.0 ; x5(3,1)=215.0
    x5(1,2)=260.0 ; x5(2,2)=260.0 ; x5(3,2)=260.0
    x5(1,3)=305.0 ; x5(2,3)=305.0 ; x5(3,3)=305.0

    y5(1,1)=35.2643896827547 ; y5(2,1)=7.7907360686059e-15  ; y5(3,1)=-35.2643896827546
    y5(1,2)=45.0000000000000 ; y5(2,2)=5.50888230454586e-15 ; y5(3,2)=-45.0000000000000
    y5(1,3)=35.2643896827547 ; y5(2,3)=0.00000000000000     ; y5(3,3)=-35.2643896827547

    dx5(1,1)=3921221.22393904 ; dx5(2,1)=3921221.22393904
    dx5(1,2)=5003771.69900514 ; dx5(2,2)=5003771.69900514
    dx5(1,3)=3921221.22393904 ; dx5(2,3)=3921221.22393904

    dy5(1,1)=3921221.22393904 ; dy5(2,1)=5003771.69900514 ; dy5(3,1)=3921221.22393904
    dy5(1,2)=3921221.22393904 ; dy5(2,2)=5003771.69900514 ; dy5(3,2)=3921221.22393904

    area5(1,1)=21252686329574.5 ; area5(2,1)=21252686329574.5
    area5(1,2)=21252686329574.6 ; area5(2,2)=21252686329574.5

    angle_dx5 = -90.0

    angle_dy5(1,1)=12.2083047812356  ; angle_dy5(2,1)=-2.90534644770403e-15    ; angle_dy5(3,1)=-12.2083047812356
    angle_dy5(1,2)=1.21499654755633e-14 ; angle_dy5(2,2)=-4.95973662257179e-15 ; angle_dy5(3,2)=-2.02499424592721e-14
    angle_dy5(1,3)=-12.2083047812356 ; angle_dy5(2,3)=-7.01412679743956e-15    ; angle_dy5(3,3)=12.2083047812355

    call call_fms2_io(tile5_file, tile, x5, y5, dx5, dy5, area5, angle_dx5, angle_dy5)

  end subroutine write_c1_tile5
 !-----------------------------------!
 subroutine write_c1_tile6

   implicit none

   character(5) :: tile

   tile='tile6'

   x6(1,1)=215.0 ; x6(2,1)=170.0 ; x6(3,1)=125.0
   x6(1,2)=260.0 ; x6(2,2)=0.0   ; x6(3,2)=80.0
   x6(1,3)=305.0 ; x6(2,3)=350.0 ; x6(3,3)=35.0

   y6(1,1)=-35.2643896827546 ; y6(2,1)=-45.0 ; y6(3,1)=-35.2643896827547
   y6(1,2)=-45.0000000000000 ; y6(2,2)=-90.0 ; y6(3,2)=-45.0000000000000
   y6(1,3)=-35.2643896827547 ; y6(2,3)=-45.0 ; y6(3,3)=-35.2643896827547

   dx6(1,1)=3921221.22393904 ; dx6(2,1)=3921221.22393904
   dx6(1,2)=5003771.69900514 ; dx6(2,2)=5003771.69900514
   dx6(1,3)=3921221.22393904 ; dx6(2,3)=3921221.22393904

   dy6(1,1)=3921221.22393904 ; dy6(2,1)=5003771.69900514 ; dy6(3,1)=3921221.22393904
   dy6(1,2)=3921221.22393904 ; dy6(2,2)=5003771.69900514 ; dy6(3,2)=3921221.22393904

   area6(1,1)=21252686329574.5 ; area6(2,1)=21252686329574.5
   area6(1,2)=21252686329574.6 ; area6(2,2)=21252686329574.5

   angle_dx6(1,1)=-167.791695218764 ; angle_dx6(2,1)=-180.0 ; angle_dx6(3,1)=-180.0
   angle_dx6(1,2)=-170.179794711206 ; angle_dx6(2,2)=-180.0 ; angle_dx6(3,2)=57.4060711511999
   angle_dx6(1,3)=-12.2083047812356 ; angle_dx6(2,3)=180.0  ; angle_dx6(3,3)=165.704009005321

   angle_dy6(1,1)=-12.2083047812356 ; angle_dy6(2,1)=-165.172188318013 ; angle_dy6(3,1)=-167.791695218764
   angle_dy6(1,2)=-2.02499424592721e-14 ; angle_dy6(2,2)=-4.05113034111413e-15 ; angle_dy6(3,2)=180.0
   angle_dy6(1,3)=12.2083047812355  ; angle_dy6(2,3)=7.32708756196959  ; angle_dy6(3,3)=167.791695218764

   call call_fms2_io(tile6_file, tile, x6, y6, dx6, dy6, area6, angle_dx6, angle_dy6)

 end subroutine write_c1_tile6
 !-----------------------------------!
  subroutine call_fms2_io(filename, tile, x, y, dx, dy, area, angle_dx, angle_dy)

    implicit none

    character(*) :: filename
    character(*) :: tile
    real, dimension(c1_nxp,c1_nyp), intent(in) :: x, y, angle_dx, angle_dy
    real, dimension(c1_nx,c1_ny),   intent(in) :: area
    real, dimension(c1_nx,c1_nyp),  intent(in) :: dx
    real, dimension(C1_nxp,C1_ny),  intent(in) :: dy

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
       call register_field(fileobj, 'dx', 'double', dimensions=(/'nx','nyp'/))
       call register_field(fileobj, 'dy', 'double', dimensions=(/'nxp','ny'/))
       call register_field(fileobj, 'area', 'double', dimensions=(/'nx','ny'/))
       call register_field(fileobj, 'angle_dx','double',dimensions=(/'nxp','nyp'/))
       call register_field(fileobj, 'angle_dy','double',dimensions=(/'nxp','nyp'/))

       call write_data(fileobj, 'tile', trim(tile))
       call write_data(fileobj, 'x', x)
       call write_data(fileobj, 'y', y)
       call write_data(fileobj, 'dx', dx)
       call write_data(fileobj, 'dy', dy)
       call write_data(fileobj, 'area', area)
       call write_data(fileobj, 'angle_dx', angle_dx)
       call write_data(fileobj, 'angle_dy', angle_dy)

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

    !< made up numbers
    do i=1,ncells
       tile1_cell(1,i) = i
       tile1_cell(2,i) = i
       tile2_cell(1,i) = i
       tile2_cell(2,i) = i
    end do

    !< made up numbers
    do i=1, ncells
       xgrid_area(i) = real(i) * PI
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
    call write_c1_tile1()
    call write_c1_tile2()
    call write_c1_tile3()
    call write_c1_tile4()
    call write_c1_tile5()
    call write_c1_tile6()
    call write_hgrid()
    call write_exchange()

  end subroutine write_all
  !---------------------------------!
end module write_files
