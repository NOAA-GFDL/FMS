!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Tests offloaded IO for a cubed-sphere (6 tile) domain with one writer PE per tile.
!! Each writer should produce <name>.tileN.nc with the data of tile N.
program test_io_offloading_cubic
  use fms_mod, only: fms_init, fms_end, check_nml_error
  use mpp_mod, only: mpp_pe, mpp_set_current_pelist, mpp_sync, mpp_error, FATAL, input_nml_file
  use mpp_domains_mod, only: domain2d, mpp_define_cubic_domain, mpp_get_data_domain, mpp_get_compute_domain, &
                             mpp_get_tile_id, mpp_get_ntile_count, mpp_deallocate_domain
  use fms2_io_mod, only: FmsNetcdfFile_t, FmsNetcdfDomainFile_t, open_file, close_file, read_data
  use offloading_io_mod
  use platform_mod, only: i4_kind
  use netcdf
  implicit none

  integer, parameter :: n = 48, ntiles = 6, nsteps = 5, nextra = 4
  integer, allocatable :: all_pes(:), model_pes(:), writer_pes(:)
  type(domain2d) :: domain
  type(FmsNetcdfFile_t) :: fileobj
  type(FmsNetcdfDomainFile_t) :: offload_file
  integer :: is, ie, js, je, isc, iec, jsc, jec, tile_ids(1), tile
  integer :: i, j, k, step, global_size(2), ncid, varid, dimids(2), dimlen, i4_att
  real, allocatable :: data_domain(:,:), compute_domain(:,:), actual(:,:)
  character(len=32) :: field_name
  character(len=64) :: tile_file
  character(len=*), parameter :: filename = "offloaded_cubic.nc"
  integer :: io_status
  !> Number of writer PEs; anything but 6 must be rejected for a cubic domain.
  integer :: nwriters = ntiles
  !> Fewer pending writes than writes per step, so send buffers are reused.
  integer :: max_pending_writes = 2
  namelist /test_io_offloading_cubic_nml/ nwriters, max_pending_writes

  call fms_init()
  read(input_nml_file, nml=test_io_offloading_cubic_nml, iostat=io_status)
  io_status = check_nml_error(io_status, "test_io_offloading_cubic_nml")
  call offloading_io_init(all_pes, model_pes, writer_pes, nwriters=nwriters)

  ! model pe
  if (any(model_pes == mpp_pe())) then
    ! domain set up
    call mpp_set_current_pelist(model_pes)
    call mpp_define_cubic_domain(n, n, ntiles, (/1,1/), nhalos=1, domain_out=domain)
    if (mpp_get_ntile_count(domain) /= ntiles) call mpp_error(FATAL, "test: cubic domain is not 6 tiles")
    tile_ids = mpp_get_tile_id(domain)
    tile = tile_ids(1)
    call mpp_get_data_domain(domain, is, ie, js, je)
    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    allocate(data_domain(is:ie, js:je), compute_domain(isc:iec, jsc:jec))
    ! send domain info and open the offloaded file
    call send_domain_metadata(domain, max_pending_writes=max_pending_writes)
    call open_offloaded_file(filename, offload_file)
    ! register axes/fields and send writing data 
    call register_axes_offloaded_file((/"grid_xt", "grid_yt"/), (/"x", "y"/), (/n, n/), domain, offload_file)
    call global_metadata_offload("ntiles", int(ntiles, i4_kind), offload_file)
    call register_field_offloaded_file("temperature", (/1,2/), fileobj=offload_file)
    do k = 1, nextra
      call register_field_offloaded_file(extra_name(k), (/1,2/), fileobj=offload_file)
    enddo
    do step = 1, nsteps
      data_domain = -999.
      do j = jsc, jec
        do i = isc, iec
          data_domain(i,j) = expected(tile, step, i, j)
        enddo
      enddo
      call write_field_offloaded_file("temperature", data_domain, domain, offload_file)
    enddo
    do k = 1, nextra
      do j = jsc, jec
        do i = isc, iec
          compute_domain(i,j) = -expected(tile, k, i, j)
        enddo
      enddo
      call write_field_offloaded_file(extra_name(k), compute_domain, domain, offload_file)
    enddo
    call close_offloaded_file(offload_file)
    call offloading_io_end()
  ! writer pes
  else
    call recv_domain_metadata(domain, global_size)
    if (any(global_size /= (/n, n/))) call mpp_error(FATAL, "test: incorrect tile size metadata")
    call run_offload_writer()
  endif

  call mpp_deallocate_domain(domain)
  call mpp_set_current_pelist(all_pes, no_sync=.true.)
  call mpp_sync()

  ! check data for correctness
  if (mpp_pe() == all_pes(1)) then
    call mpp_set_current_pelist((/all_pes(1)/), no_sync=.true.)
    if (nf90_open(filename, nf90_nowrite, ncid) == nf90_noerr) &
      call mpp_error(FATAL, "test: an untiled file was written for a cubic domain")
    allocate(actual(n,n))
    do tile = 1, ntiles
      write(tile_file, '(a,i0,a)') "offloaded_cubic.tile", tile, ".nc"
      call check_nc(nf90_open(tile_file, nf90_nowrite, ncid))
      call check_nc(nf90_inq_varid(ncid, "temperature", varid))
      call check_nc(nf90_inquire_variable(ncid, varid, dimids=dimids))
      do k = 1, 2
        call check_nc(nf90_inquire_dimension(ncid, dimids(k), len=dimlen))
        if (dimlen /= n) call mpp_error(FATAL, "test: incorrect dimension size in "//trim(tile_file))
      enddo
      call check_nc(nf90_get_att(ncid, nf90_global, "ntiles", i4_att))
      if (i4_att /= ntiles) call mpp_error(FATAL, "test: incorrect global attribute in "//trim(tile_file))
      call check_nc(nf90_close(ncid))

      if (.not. open_file(fileobj, tile_file, "read")) call mpp_error(FATAL, "test: missing "//trim(tile_file))
      call read_data(fileobj, "temperature", actual)
      call check_field(actual, tile, nsteps, 1., "temperature")
      do k = 1, nextra
        field_name = extra_name(k)
        call read_data(fileobj, trim(field_name), actual)
        call check_field(actual, tile, k, -1., field_name)
      enddo
      call close_file(fileobj)
    enddo
  endif
  call mpp_set_current_pelist(all_pes, no_sync=.true.)
  call fms_end()

contains

  pure real function expected(tile, step, i, j)
    integer, intent(in) :: tile, step, i, j
    expected = real(1000000*tile + 10000*step + 100*j + i)
  end function

  character(len=32) function extra_name(k)
    integer, intent(in) :: k
    write(extra_name, '(a,i0)') "field_", k
  end function

  subroutine check_field(values, tile, step, sign, name)
    real, intent(in) :: values(:,:), sign
    integer, intent(in) :: tile, step
    character(len=*), intent(in) :: name
    integer :: i, j

    do j = 1, size(values, 2)
      do i = 1, size(values, 1)
        if (values(i,j) /= sign*expected(tile, step, i, j)) &
          call mpp_error(FATAL, "test: "//trim(name)//" data mismatch in "//trim(tile_file))
      enddo
    enddo
  end subroutine

  subroutine check_nc(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) call mpp_error(FATAL, "test: "//trim(nf90_strerror(status)))
  end subroutine
end program test_io_offloading_cubic
