program main

  ! Test for horiz_interp_read_weights_conserve where
  ! the source and destination grids are identical

  use constants_mod, only: DEG_TO_RAD
  use fms_mod, only: fms_init, fms_end
  use horiz_interp_mod, only: horiz_interp_type, horiz_interp_read_weights_conserve, horiz_interp
  use mpp_mod, only: FATAL, mpp_error, mpp_npes, mpp_pe, mpp_root_pe, mpp_sync
  use mpp_domains_mod, only: domain2d, mpp_define_layout, mpp_define_domains
  use mpp_domains_mod, only: mpp_get_compute_domain
  use fms2_io_mod, only: FmsNetcdfFile_t, open_file, close_file
  use fms2_io_mod, only: register_axis, register_field, write_data

  use platform_mod, only: r8_kind

  implicit none

  integer, parameter :: nlon_dst = 32 !< number of lon cells in destination grid
  integer, parameter :: nlat_dst = 16 !< number of lat cells in destionation grid
  integer, parameter :: nlon_src = nlon_dst !< number of lon cells in src grid
  integer, parameter :: nlat_src = nlat_dst !< number of lat cells in src_grid
  integer, parameter :: ncells = nlon_dst*nlat_dst !< number of exchange grid cells

  type(horiz_interp_type) :: Interp !< Interp to be populated
  character(8) :: remap = "remap.nc" !< weight file name

  type(domain2d) :: domain !< domain
  integer :: global_indices(4) = [1, nlon_dst, 1, nlat_dst] !< global_indices
  integer :: layout(2) !< layout
  integer :: isc, iec, jsc, jec !< compute domain indices on dst grid
  integer :: nlon_dst_d, nlat_dst_d !< dst grid size on pe domain

  real(8), allocatable :: data_src(:,:), data_dst(:,:) !< data for remapping

  ! counter
  integer :: i,j

  ! initialize fms
  call fms_init()

  ! get domain
  call mpp_define_layout(global_indices, mpp_npes(), layout)
  call mpp_define_domains(global_indices, layout, domain, tile_id=1)
  call mpp_get_compute_domain(domain, isc, iec, jsc, jec, xsize=nlon_dst_d, ysize=nlat_dst_d)

  ! write remap file to test
  if(mpp_pe() == mpp_root_pe()) call write_remap(nlon_dst, nlat_dst, ncells)
  call mpp_sync()

  call fms_end()
  stop

  ! allocate src data
  allocate(data_src(nlon_src, nlat_src))
  do j=1, nlat_src
    do i=1, nlon_src
      data_src(i,j) = j*nlon_src + i
    end do
  end do

  ! allocate dst data on domain
  allocate(data_dst(isc:iec, jsc:jec))
  data_dst = -99.0

  ! read remap file
  call horiz_interp_read_weights_conserve(Interp, remap, "fregrid", nlon_src, nlat_src, &
    nlon_dst_d, nlat_dst_d, isc, iec, jsc, jec, 1)

  ! remap
  call horiz_interp(Interp, data_src, data_dst)

  ! check answers
  if(.not.all(data_src(isc:iec, jsc:jec) == data_dst)) then
    call mpp_error(FATAL, "error remapping with read-in exchange grid")
  end if

  call fms_end()

contains

  subroutine write_remap(nlon_dst, nlat_dst, ncells)

    implicit none
    integer, intent(in) :: nlon_dst, nlat_dst, ncells

    type(FmsNetCDFfile_t) :: fileobj
    integer, allocatable :: tile1(:), tile1_cell(:, :), tile2_cell(:, :)
    real(r8_kind), allocatable :: lon_dst(:,:), lat_dst(:,:), xgrid_area(:)

    integer :: i, j, ij

    allocate(lon_dst(nlon_dst+1, nlat_dst+1), lat_dst(nlon_dst+1, nlat_dst+1))
    ! create grid in radians
    do j=1, nlat_dst+1
      do i=1, nlon_dst+1
        lon_dst(i,j) = real(i, r8_kind) * DEG_TO_RAD
        lat_dst(i,j) = real(j, r8_kind) * DEG_TO_RAD
      end do
    end do

    ! get grid area, make xgrid_area equal to the dst grid area
    allocate(xgrid_area(ncells))
    call get_grid_area(nlon_dst, nlat_dst, lon_dst, lat_dst, xgrid_area)
    deallocate(lon_dst, lat_dst)

    ! allocate arrays to write
    allocate(tile1(ncells), tile1_cell(2, ncells), tile2_cell(2, ncells))

    ! xgrid
    ij = 1
    do j=1, nlat_src
      do i=1, nlon_src
        tile1_cell(1, ij) = i
        tile1_cell(2, ij) = j
        ij = ij + 1
      end do
    end do
    tile2_cell = tile1_cell
    tile1 = 1

   !write remap file
    if(open_file(fileobj, "remap.nc", "overwrite")) then

      ! register axis
      call register_axis(fileobj, "two", 2)
      call register_axis(fileobj, "ncells", ncells)

      !register field
      call register_field(fileobj, "tile1", "int", dimensions=["ncells"])
      call register_field(fileobj, "tile1_cell", "int", dimensions=["two   ", "ncells"])
      call register_field(fileobj, "tile2_cell", "int", dimensions=["two   ", "ncells"])
      call register_field(fileobj, "xgrid_area", "double", dimensions=["ncells"])

      call write_data(fileobj, "tile1", tile1)
      call write_data(fileobj, "tile1_cell", tile1_cell)
      call write_data(fileobj, "tile2_cell", tile2_cell)
      call write_data(fileobj, "xgrid_area", xgrid_area)

      call close_file(fileobj)
    end if

    deallocate(tile1, tile1_cell, tile2_cell, xgrid_area)

  end subroutine write_remap

end program main
