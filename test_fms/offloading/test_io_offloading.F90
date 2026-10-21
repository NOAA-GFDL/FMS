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
program test_io_offloading
  use fms_mod, only: fms_init, fms_end, check_nml_error
  use mpp_mod, only: mpp_pe, mpp_set_current_pelist, mpp_sync, mpp_error, FATAL, input_nml_file
  use mpp_domains_mod, only: domain2d, mpp_define_layout, mpp_define_domains, mpp_define_io_domain, &
                            mpp_get_data_domain, mpp_get_compute_domain, mpp_deallocate_domain
  use fms2_io_mod, only: FmsNetcdfFile_t, open_file, close_file, read_data
  use offloading_io_mod
  use platform_mod, only: r4_kind, r8_kind, i4_kind
  use netcdf
  implicit none

  integer, parameter :: nx = 128, ny = 96
  integer, allocatable :: all_pes(:), model_pes(:), writer_pes(:)
  type(domain2d) :: domain
  type(FmsNetcdfFile_t) :: fileobj
  integer :: layout(2), is, ie, js, je, isc, iec, jsc, jec
  integer :: i, j, step, global_size(2), ncid, varid, dimids(2), dimlen, ndims
  real, allocatable :: model_data(:,:), actual(:,:)
  real(kind=r4_kind), allocatable :: model_data_r4(:,:)
  real(kind=r4_kind) :: r4_att
  real(kind=r8_kind) :: r8_att
  integer(kind=i4_kind) :: i4_att
  character(len=64) :: str_att
  character(len=9) :: axes(2) = [character(len=9) :: "latitude", "longitude"]
  character(len=*), parameter :: filename = "offloaded_output.nc"
  integer :: io_status
  ! Select normal, bad-axis, bad-dims, or bad-attribute in input.nml.
  character(len=32) :: mode = "normal"
  namelist /test_io_offloading_nml/ mode

  call fms_init()
  read(input_nml_file, nml=test_io_offloading_nml, iostat=io_status)
  io_status = check_nml_error(io_status, "test_io_offloading_nml")
  select case (trim(mode))
  case ("normal", "bad-axis", "bad-dims", "bad-attribute")
  case default
    call mpp_error(FATAL, "test_io_offloading: unknown mode "//trim(mode))
  end select
  call setup_offloading_pelists(all_pes, model_pes, writer_pes)
  ! model pes section: define domain and send metadata and i/o ops to the writer pe
  if (any(model_pes == mpp_pe())) then
    call mpp_set_current_pelist(model_pes)
    call mpp_define_layout((/1,nx,1,ny/), size(model_pes), layout)
    call mpp_define_domains((/1,nx,1,ny/), layout, domain, xhalo=1, yhalo=1)
    call mpp_define_io_domain(domain, (/1,1/))
    call mpp_get_data_domain(domain, is, ie, js, je)
    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    allocate(model_data(is:ie, js:je))
    allocate(model_data_r4(is:ie, js:je))
    call send_domain_metadata(domain)
    call open_offloaded_file(filename)
    ! Deliberately register y before x to exercise transmitted dimension indices.
    call register_axes_offloaded_file(axes, (/"y", "x"/), (/ny, nx/), domain)
    if (trim(mode) == "bad-axis") call register_axis_offloaded_file("bad", "z", nx)
    if (trim(mode) == "bad-dims") call register_field_offloaded_file("bad", (/1/))
    if (trim(mode) == "bad-attribute") call global_metadata_offload("bad", .true.)
    call global_metadata_offload("Number of times Fortran made you cry", 20_i4_kind)
    call global_metadata_offload("Number of lines of code", 19.54326541_r4_kind)
    call global_metadata_offload("pi", 3.141592653589793_r8_kind)
    call register_field_offloaded_file("temperature", (/2,1/))
    call register_field_offloaded_file("rainfall", (/2,1/))
    call register_field_offloaded_file("temperature_r4", (/2,1/), "float")
    do step = 1, 3
      model_data = -999.
      do j = jsc, jec
        do i = isc, iec
          model_data(i,j) = real(1000*step + 100*j + i)
        enddo
      enddo
      call write_field_offloaded_file("temperature", model_data, domain)
    enddo
    model_data = -999.
    do j = jsc, jec
      do i = isc, iec
        model_data(i,j) = real(-100*j - i)
      enddo
    enddo
    call write_field_offloaded_file("rainfall", model_data, domain)
    model_data_r4 = -999.0_r4_kind
    do j = jsc, jec
      do i = isc, iec
        model_data_r4(i,j) = real(2000 + 100*j + i, r4_kind)
      enddo
    enddo
    call write_field_offloaded_file("temperature_r4", model_data_r4, domain)
    ! Registered after data writes to exercise switching the file back to define mode.
    call global_metadata_offload("title", "offloading test")
    call close_offloaded_file()

  else ! writer PE

    call recv_domain_metadata(domain, global_size)
    if (any(global_size /= (/nx, ny/))) call mpp_error(FATAL, "test: incorrect global domain metadata")
    call run_offload_writer()

  endif

  call mpp_deallocate_domain(domain)
  call mpp_set_current_pelist(all_pes, no_sync=.true.)
  call mpp_sync()

  if (mpp_pe() == all_pes(1)) then
    call mpp_set_current_pelist((/all_pes(1)/), no_sync=.true.)
    call check_nc(nf90_open(filename, nf90_nowrite, ncid))
    call check_nc(nf90_inq_varid(ncid, "temperature", varid))
    call check_nc(nf90_inquire_variable(ncid, varid, ndims=ndims))
    if (ndims /= 2) call mpp_error(FATAL, "test: field rank differs from two")
    call check_nc(nf90_inquire_variable(ncid, varid, dimids=dimids))
    call check_nc(nf90_inquire_dimension(ncid, dimids(1), len=dimlen))
    if (dimlen /= nx) call mpp_error(FATAL, "test: incorrect x dimension size")
    call check_nc(nf90_inquire_dimension(ncid, dimids(2), len=dimlen))
    if (dimlen /= ny) call mpp_error(FATAL, "test: incorrect y dimension size")
    call check_nc(nf90_inq_dimid(ncid, "longitude", varid))
    if (varid /= dimids(1)) call mpp_error(FATAL, "test: incorrect x dimension name/order")
    call check_nc(nf90_inq_dimid(ncid, "latitude", varid))
    if (varid /= dimids(2)) call mpp_error(FATAL, "test: incorrect y dimension name/order")
    call check_nc(nf90_get_att(ncid, nf90_global, "Number of times Fortran made you cry", i4_att))
    if (i4_att /= 20_i4_kind) call mpp_error(FATAL, "test: incorrect i4 global attribute")
    call check_nc(nf90_get_att(ncid, nf90_global, "Number of lines of code", r4_att))
    if (r4_att /= 19.54326541_r4_kind) call mpp_error(FATAL, "test: incorrect r4 global attribute")
    call check_nc(nf90_get_att(ncid, nf90_global, "pi", r8_att))
    if (r8_att /= 3.141592653589793_r8_kind) call mpp_error(FATAL, "test: incorrect r8 global attribute")
    str_att = ""
    call check_nc(nf90_get_att(ncid, nf90_global, "title", str_att))
    if (trim(str_att) /= "offloading test") &
      call mpp_error(FATAL, "test: incorrect string global attribute")
    call check_nc(nf90_close(ncid))
    allocate(actual(nx,ny))
    if (.not. open_file(fileobj, filename, "read")) call mpp_error(FATAL, "test: output file missing")
    call read_data(fileobj, "temperature", actual)
    do j = 1, ny
      do i = 1, nx
        if (actual(i,j) /= real(3000 + 100*j + i)) call mpp_error(FATAL, "test: temperature data mismatch")
      enddo
    enddo
    call read_data(fileobj, "rainfall", actual)
    do j = 1, ny
      do i = 1, nx
        if (actual(i,j) /= real(-100*j - i)) call mpp_error(FATAL, "test: rainfall data mismatch")
      enddo
    enddo
    call read_data(fileobj, "temperature_r4", actual)
    do j = 1, ny
      do i = 1, nx
        if (actual(i,j) /= real(2000 + 100*j + i)) &
          call mpp_error(FATAL, "test: r4 temperature data mismatch")
      enddo
    enddo
    call close_file(fileobj)
  endif
  call mpp_set_current_pelist(all_pes, no_sync=.true.)
  call fms_end()
contains
  subroutine check_nc(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) call mpp_error(FATAL, "test: "//trim(nf90_strerror(status)))
  end subroutine
end program test_io_offloading
