program main

  use constants_mod, only: DEG_TO_RAD
  use fms_mod, only: fms_init, fms_end
  use horiz_interp_mod, only: horiz_interp_type, horiz_interp
  use horiz_interp_type_mod, only: CONSERVE
  use mpp_mod, only: FATAL, mpp_error, mpp_npes, mpp_pe, mpp_root_pe, mpp_sync
  use mpp_domains_mod, only: domain2d, mpp_define_layout, mpp_define_domains
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_tile_id
  use fms2_io_mod, only: FmsNetcdfFile_t, open_file, get_dimension_size, read_data, close_file
  use fms2_io_mod, only: register_axis, register_field, write_data

  implicit none

  integer, parameter :: nlon_dst = 16
  integer, parameter :: nlat_dst = 16
  integer, parameter :: nlon_src = nlon_dst
  integer, parameter :: nlat_src = nlat_dst
  integer, parameter :: ncells = nlon_dst*nlat_dst

  type(horiz_interp_type) :: Interp
  character(8) :: remap = "remap.nc"

  integer :: global_indices(4) = [1, nlon_dst, 1, nlat_dst]
  integer :: layout(2)
  integer :: isc, iec, jsc, jec
  integer :: nlon_dst_d, nlat_dst_d, nlon_src_d, nlat_src_d
  type(domain2d) :: domain

  real(8), allocatable :: data_src(:,:), data_dst(:,:)

  integer :: i, j

  ! initialize fms
  call fms_init()

  call mpp_define_layout(global_indices, mpp_npes(), layout)
  call mpp_define_domains(global_indices, layout, domain, tile_id=1)
  call mpp_get_compute_domain(domain, isc, iec, jsc, jec)

  nlon_dst_d = iec - isc + 1
  nlat_dst_d = jec - jsc + 1

  if(mpp_pe() == mpp_root_pe()) call write_remap(nlon_dst, nlat_dst)
  call mpp_sync()

  call read_weights(Interp, remap, "fregrid", nlon_src, nlat_src, nlon_dst_d, nlat_dst_d, domain, src_tile=1)

  !test
  allocate(data_src(nlon_src, nlat_src))
  allocate(data_dst(isc:iec, jsc:jec))
  do j=1, nlat_src
    do i=1, nlon_src
      data_src(i,j) = j*nlon_src + i
    end do
  end do

  call horiz_interp(Interp, data_src, data_dst)

  if(all(data_src(isc:iec, jsc:jec) == data_dst)) then
    write(*,*) 'SUCCESS'
  else
    write(*,*) 'FAIL'
  end if

  call fms_end()

contains

  subroutine read_weights(Interp, weight_filename, weight_file_source, nlon_src, nlat_src, &
    nlon_dst, nlat_dst, dst_domain, src_tile)

    implicit none
    type(horiz_interp_type), intent(inout) :: Interp !< Horiz interp time to fill
    character(len=*), intent(in) :: weight_filename !< Name of the weight file
    character(len=*), intent(in) :: weight_file_source !< Source of the weight file
    integer, intent(in) :: nlon_src, nlat_src
    integer, intent(in) :: nlon_dst, nlat_dst
    type(domain2d), intent(in) :: dst_domain !< domain, only one can be present
    integer, intent(in), optional :: src_tile

    integer :: tile_id(1) !< tile id for saving xgrid for the correct tile
    integer :: ncells, mpi_ncells !< number of exchange grid cells
    integer :: isc, iec, jsc, jec !< compute indices for dst domain
    integer :: istart(1), iend(1), i_dst, j_dst
    real(8), allocatable :: dst_area(:,:), var1(:)
    integer, allocatable :: tile1(:), indices(:), var2(:,:)
    type(FmsNetcdfFile_t) :: weight_fileobj !< FMS2io fileob for the weight file

    integer :: i

    ! check if weight_file was generated from fregrid
    if(trim(weight_file_source) /= "fregrid") then
      call mpp_error(FATAL, trim(weight_file_source)//&
        &" is not a supported weight file source. fregrid is the only supported weight file source.")
    end if

    ! get compute and global domain
    call mpp_get_compute_domain(dst_domain, isc, iec, jsc, jec)

    if(open_file(weight_fileobj, trim(weight_filename), "read")) then

      ! get ncells
      call get_dimension_size(weight_fileobj, "ncells", ncells)

      !allocate temporary data
      allocate(tile1(ncells))

      call read_data(weight_fileobj, "tile1", tile1)
      if(present(src_tile)) then
        istart = FINDLOC(tile1, src_tile)
        iend = FINDLOC(tile1, src_tile, back=.true.)
        ncells = iend(1) - istart(1) + 1
      end if

      deallocate(tile1)

      ! allocate fields to read in xgrid
      allocate(var2(2, ncells))
      allocate(var1(ncells))
      allocate(indices(ncells))
      allocate(dst_area(nlon_dst, nlat_dst))

      mpi_ncells = 1
      call read_data(weight_fileobj, "tile2_cell", var2, corner=[1,istart(1)], edge_lengths=[2,iend(1)])

      do i=1, ncells
        i_dst = var2(1,i)
        j_dst = var2(2,i)
        if(j_dst < jsc .or. j_dst > jec) cycle
        if( i_dst < isc .or. i_dst > iec) cycle
        indices(mpi_ncells) = i
        mpi_ncells = mpi_ncells + 1
      end do

      mpi_ncells = mpi_ncells - 1

      allocate(Interp%i_src(mpi_ncells))
      allocate(Interp%j_src(mpi_ncells))
      allocate(Interp%i_dst(mpi_ncells))
      allocate(Interp%j_dst(mpi_ncells))
      allocate(Interp%horizInterpReals8_type%area_frac_dst(mpi_ncells))

      do i=1, mpi_ncells
        Interp%i_dst(i) = var2(1,indices(i)) - isc + 1
        Interp%j_dst(i) = var2(2,indices(i)) - jsc + 1
      end do

      call read_data(weight_fileobj, "tile1_cell", var2, corner=[1,istart(1)], edge_lengths=[2,iend(1)])
      do i=1, mpi_ncells
        Interp%i_src(i) = var2(1, indices(i))
        Interp%j_src(i) = var2(2, indices(i))
      end do


      call read_data(weight_fileobj, "xgrid_area", var1, corner=[istart(1)], edge_lengths=[iend(1)])

      !sum over exchange grid area to get destination grid area
      dst_area = 0.0
      do i = 1, mpi_ncells
        dst_area(Interp%i_dst(i), Interp%j_dst(i)) =+ var1(indices(i))
        Interp%horizInterpReals8_type%area_frac_dst(i) = var1(indices(i))
      end do

      do i=1, mpi_ncells
        Interp%horizInterpReals8_type%area_frac_dst(i) = dst_area(Interp%i_dst(i), Interp%j_dst(i))/Interp%horizInterpReals8_type%area_frac_dst(i)
      end do

      call close_file(weight_fileobj)

      deallocate(var2)
      deallocate(var1)
      deallocate(indices)
      deallocate(dst_area)

      ! if(mpp_pe()==1) then
      !   write(*,*) Interp%i_src
      !   write(*,*) Interp%j_src
      !   write(*,*) Interp%i_dst
      !   write(*,*) Interp%j_dst
      !   write(*,*) Interp%horizInterpReals8_type%area_frac_dst
      ! end if

    else
      call mpp_error(FATAL, "cannot open weight file")
    end if

    Interp%nxgrid = mpi_ncells
    Interp%nlon_src = nlon_src
    Interp%nlat_src = nlat_src
    Interp%nlon_dst = nlon_dst
    Interp%nlat_dst = nlat_dst
    Interp%horizInterpReals8_type%is_allocated = .true.
    Interp%interp_method = CONSERVE
    Interp%version = 2
    Interp%I_am_initialized = .true.

  end subroutine read_weights


  subroutine write_remap(nlon_dst, nlat_dst)

    implicit none

    integer, intent(in) :: nlon_dst, nlat_dst

    type(FmsNetCDFfile_t) :: fileobj

    integer :: tile1(ncells), tile1_cell(2, ncells), tile2_cell(2, ncells)

    real(8) :: lon_dst(nlon_dst+1, nlat_dst+1), lat_dst(nlon_dst+1, nlat_dst+1), xgrid_area(ncells)

    integer :: i, j, ij

    do j=1, nlat_dst+1
      do i=1, nlon_dst+1
        lon_dst(i,j) = i
        lat_dst(i,j) = j
      end do
    end do
    lon_dst = lon_dst * DEG_TO_RAD
    lat_dst = lat_dst * DEG_TO_RAD

    call get_grid_area(nlon_dst, nlat_dst, lon_dst, lat_dst, xgrid_area)

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
    if(mpp_pe() == mpp_root_pe()) then
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
    end if

  end subroutine write_remap

end program main
