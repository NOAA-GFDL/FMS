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
!> Tests offloaded IO with more than one writer PE per tile, for single-tile and cubed-sphere
!! domains. The number of writers must be a multiple of the number of tiles.
!! After the run, the writers read the files back through a domain with the same decomposition,
!! which checks both combined (writer_io_layout = 1,1) and distributed (.nc.XXXX) output.
program test_io_offloading_multi_writer
  use fms_mod, only: fms_init, fms_end, check_nml_error
  use mpp_mod, only: mpp_pe, mpp_set_current_pelist, mpp_sync, mpp_error, FATAL, input_nml_file
  use mpp_domains_mod, only: domain2d, mpp_define_cubic_domain, mpp_define_domains, mpp_define_layout, &
                             mpp_define_io_domain, mpp_define_mosaic, mpp_get_data_domain, &
                             mpp_get_compute_domain, mpp_get_tile_id, mpp_deallocate_domain
  use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file, read_data, register_axis
  use offloading_io_mod
  implicit none

  integer, parameter :: n = 48, nsteps = 3
  integer, allocatable :: all_pes(:), model_pes(:), writer_pes(:)
  type(domain2d) :: domain, check_domain
  type(FmsNetcdfDomainFile_t) :: fileobj, offload_file
  integer :: is, ie, js, je, isc, iec, jsc, jec, tile_ids(1), tile, layout(2)
  integer :: i, j, step, global_size(2), writers_per_tile, io_status
  real, allocatable :: model_data(:,:), actual(:,:)
  character(len=*), parameter :: filename = "offloaded_multi.nc"
  character(len=64) :: expected_file
  logical :: found
  !> 1 for a single tile domain, 6 for a cubed sphere
  integer :: ntiles = 6
  !> Total number of writer PEs; must be a multiple of ntiles.
  integer :: nwriters = 12
  !> Writer layout per tile; (0,0) uses the default from mpp_define_layout.
  integer :: writer_layout(2) = 0
  !> io_layout of the writers of each tile; (1,1) writes one combined file per tile.
  integer :: writer_io_layout(2) = 1
  namelist /test_io_offloading_multi_writer_nml/ ntiles, nwriters, writer_layout, writer_io_layout

  call fms_init()
  read(input_nml_file, nml=test_io_offloading_multi_writer_nml, iostat=io_status)
  io_status = check_nml_error(io_status, "test_io_offloading_multi_writer_nml")
  call offloading_io_init(all_pes, model_pes, writer_pes, nwriters=nwriters)
  writers_per_tile = nwriters / ntiles

  if (any(model_pes == mpp_pe())) then
    call mpp_set_current_pelist(model_pes)
    if (ntiles == 6) then
      call mpp_define_cubic_domain(n, n, ntiles, (/1,1/), nhalos=1, domain_out=domain)
    else
      call mpp_define_layout((/1,n,1,n/), size(model_pes), layout)
      call mpp_define_domains((/1,n,1,n/), layout, domain, xhalo=1, yhalo=1)
      call mpp_define_io_domain(domain, (/1,1/))
    endif
    tile_ids = mpp_get_tile_id(domain)
    tile = tile_ids(1)
    call mpp_get_data_domain(domain, is, ie, js, je)
    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    allocate(model_data(is:ie, js:je))
    if (all(writer_layout > 0)) then
      call send_domain_metadata(domain, writer_layout=writer_layout, writer_io_layout=writer_io_layout)
    else
      call send_domain_metadata(domain, writer_io_layout=writer_io_layout)
    endif
    call open_offloaded_file(filename, offload_file)
    call register_axes_offloaded_file((/"grid_xt", "grid_yt"/), (/"x", "y"/), (/n, n/), domain, offload_file)
    call register_field_offloaded_file("temperature", offload_file, axis_names=["grid_xt", "grid_yt"])
    do step = 1, nsteps
      model_data = -999.
      do j = jsc, jec
        do i = isc, iec
          model_data(i,j) = expected(tile, step, i, j)
        enddo
      enddo
      call write_field_offloaded_file("temperature", model_data, domain, offload_file)
    enddo
    call close_offloaded_file(offload_file)
    call offloading_io_end()
    call mpp_deallocate_domain(domain)
  else
    call recv_domain_metadata(domain, global_size)
    if (any(global_size /= (/n, n/))) call mpp_error(FATAL, "test: incorrect tile size metadata")
    call run_offload_writer()
    call mpp_deallocate_domain(domain)

    ! Read the output back with the same decomposition as the writers.
    call define_check_domain()
    tile_ids = mpp_get_tile_id(check_domain)
    tile = tile_ids(1)
    if (ntiles > 1) then
      write(expected_file, '(a,i0,a)') "offloaded_multi.tile", tile, ".nc"
    else
      expected_file = filename
    endif
    if (any(writer_io_layout > 1)) expected_file = trim(expected_file)//".0001"
    inquire(file=trim(expected_file), exist=found)
    if (.not. found) call mpp_error(FATAL, "test: missing output file "//trim(expected_file))
    call mpp_get_compute_domain(check_domain, isc, iec, jsc, jec)
    allocate(actual(isc:iec, jsc:jec))
    if (.not. open_file(fileobj, filename, "read", check_domain)) &
      call mpp_error(FATAL, "test: unable to read back the offloaded output")
    call register_axis(fileobj, "grid_xt", "x")
    call register_axis(fileobj, "grid_yt", "y")
    call read_data(fileobj, "temperature", actual)
    call close_file(fileobj)
    do j = jsc, jec
      do i = isc, iec
        if (actual(i,j) /= expected(tile, nsteps, i, j)) &
          call mpp_error(FATAL, "test: temperature data mismatch on a writer block")
      enddo
    enddo
    call mpp_deallocate_domain(check_domain)
  endif

  call mpp_set_current_pelist(all_pes, no_sync=.true.)
  call mpp_sync()
  call fms_end()

contains

  pure real function expected(tile, step, i, j)
    integer, intent(in) :: tile, step, i, j
    expected = real(1000000*tile + 10000*step + 100*j + i)
  end function

  !> Same decomposition as the offloading writers: writers_per_tile writers per tile, in order.
  subroutine define_check_domain()
    integer :: global_indices(4, ntiles), layouts(2, ntiles), pe_start(ntiles), pe_end(ntiles), t
    integer :: no_contact(0), tile_layout(2)

    tile_layout = writer_layout
    if (any(writer_layout < 1)) call mpp_define_layout((/1,n,1,n/), writers_per_tile, tile_layout)
    do t = 1, ntiles
      global_indices(:,t) = (/1, n, 1, n/)
      layouts(:,t) = tile_layout
      pe_start(t) = writer_pes((t-1)*writers_per_tile + 1)
      pe_end(t) = writer_pes(t*writers_per_tile)
    enddo
    call mpp_define_mosaic(global_indices, layouts, check_domain, ntiles, 0, no_contact, no_contact, &
                           no_contact, no_contact, no_contact, no_contact, no_contact, no_contact, no_contact, &
                           no_contact, pe_start, pe_end, name="check domain")
    call mpp_define_io_domain(check_domain, writer_io_layout)
  end subroutine
end program test_io_offloading_multi_writer
