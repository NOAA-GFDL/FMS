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
  use fms2_io_mod, only: FmsNetcdfFile_t, FmsNetcdfDomainFile_t, open_file, close_file, read_data, unlimited
  use offloading_io_mod
  use platform_mod, only: r4_kind, r8_kind, i4_kind, i8_kind
  use netcdf
  implicit none

  integer, parameter :: nx = 128, ny = 96
  integer, allocatable :: all_pes(:), model_pes(:), writer_pes(:)
  type(domain2d) :: domain
  type(FmsNetcdfFile_t) :: fileobj
  type(FmsNetcdfDomainFile_t) :: offload_file
  type(FmsNetcdfDomainFile_t) :: second_file
  integer :: layout(2), is, ie, js, je, isc, iec, jsc, jec
  integer :: i, j, step, global_size(2), ncid, varid, dimids(2), dimlen, ndims, unlimited_dimid
  integer :: record_dimids(3)
  real, allocatable :: model_data(:,:), actual(:,:)
  real(kind=r4_kind), allocatable :: model_data_r4(:,:)
  real(kind=r4_kind) :: r4_att
  real(kind=r8_kind) :: r8_att
  integer(kind=i4_kind) :: i4_att
  character(len=64) :: str_att
  character(len=9) :: axes(2) = [character(len=9) :: "latitude", "longitude"]
  character(len=*), parameter :: filename = "offloaded_output.nc"
  character(len=*), parameter :: second_filename = "offloaded_output_second.nc"
  integer :: io_status
  integer :: nwriters = 1
  ! Select normal, bad-axis, bad-dims, or bad-attribute in input.nml.
  character(len=32) :: mode = "normal"
  namelist /test_io_offloading_nml/ mode, nwriters

  call fms_init()
  read(input_nml_file, nml=test_io_offloading_nml, iostat=io_status)
  io_status = check_nml_error(io_status, "test_io_offloading_nml")
  select case (trim(mode))
  case ("normal", "bad-axis", "bad-dims", "bad-attribute", "bad-level", "bad-type")
  case default
    call mpp_error(FATAL, "test_io_offloading: unknown mode "//trim(mode))
  end select
  call offloading_io_init(all_pes, model_pes, writer_pes, nwriters=nwriters)
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
    call open_offloaded_file(filename, offload_file)
    call open_offloaded_file(second_filename, second_file)
    call register_axis_offloaded_file("level", 3, offload_file)
    ! Deliberately register y before x to exercise named-axis ordering.
    call register_axes_offloaded_file(axes, (/"y", "x"/), (/ny, nx/), domain, offload_file)
    call register_axes_offloaded_file(axes, (/"y", "x"/), (/ny, nx/), domain, second_file)
    call register_axis_offloaded_file("member", 2, offload_file, "")
    call register_axis_offloaded_file("time", axis_size=unlimited, fileobj=offload_file)
    if (trim(mode) == "bad-axis") call register_axis_offloaded_file("bad", nx, offload_file, "z")
    if (trim(mode) == "bad-dims") &
      call register_field_offloaded_file("bad", offload_file, axis_names=["missing"])
    if (trim(mode) == "bad-attribute") call global_metadata_offload("bad", .true., offload_file)
    call global_metadata_offload("Number of times Fortran made you cry", 20_i4_kind, offload_file)
    call global_metadata_offload("Number of lines of code", 19.54326541_r4_kind, offload_file)
    call global_metadata_offload("pi", 3.141592653589793_r8_kind, offload_file)
    call register_field_offloaded_file("temperature", offload_file, axis_names=axes(2:1:-1))
    call register_field_offloaded_file("independent", second_file, axis_names=axes(2:1:-1))
    call register_field_offloaded_file("rainfall", offload_file, axis_names=axes(2:1:-1))
    call register_field_offloaded_file("temperature_r4", offload_file, "float", axes(2:1:-1))
    call register_field_offloaded_file("history", offload_file, &
                       axis_names=[character(len=9) :: "longitude", "latitude", "time"])
    call register_field_offloaded_file("history_r4", offload_file, "float", &
                       [character(len=9) :: "longitude", "latitude", "time"])
    if (trim(mode) == "bad-level") &
      call write_field_offloaded_file("history", model_data, domain, offload_file, unlim_dim_level=0)
    if (trim(mode) == "bad-type") call write_field_offloaded_file("history", .true., domain, offload_file)
    do step = 1, 3
      model_data = -999.
      do j = jsc, jec
        do i = isc, iec
          model_data(i,j) = real(1000*step + 100*j + i)
          model_data_r4(i,j) = real(2000*step + 100*j + i, r4_kind)
        enddo
      enddo
      call write_field_offloaded_file("temperature", model_data, domain, offload_file)
      call write_field_offloaded_file("independent", model_data, domain, second_file)
      call write_field_offloaded_file("history", model_data, domain, offload_file, unlim_dim_level=4-step)
      call write_field_offloaded_file("history_r4", model_data_r4, domain, offload_file, unlim_dim_level=4-step)
    enddo
    model_data = -999.
    do j = jsc, jec
      do i = isc, iec
        model_data(i,j) = real(-100*j - i)
      enddo
    enddo
    call write_field_offloaded_file("rainfall", model_data, domain, offload_file)
    call write_field_offloaded_file("history", model_data, domain, offload_file)
    model_data_r4 = -999.0_r4_kind
    do j = jsc, jec
      do i = isc, iec
        model_data_r4(i,j) = real(2000 + 100*j + i, r4_kind)
      enddo
    enddo
    call write_field_offloaded_file("temperature_r4", model_data_r4, domain, offload_file)
    call generate_fields(.true.)
    call register_axis_offloaded_file("strlen", 12, offload_file)
    call register_field_offloaded_file("label", offload_file, "char", ["strlen"])
    call register_field_offloaded_file("labels", offload_file, "char", ["strlen", "member"])
    call write_field_offloaded_file("label", "offloaded", domain, offload_file)
    call write_field_offloaded_file("labels", [character(len=12) :: "first", "second"], domain, offload_file)
    ! Registered after data writes to exercise switching the file back to define mode.
    call global_metadata_offload("title", "offloading test", offload_file)
    call close_offloaded_file(second_file)
    call close_offloaded_file(offload_file)
    call offloading_io_end()

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
    call check_nc(nf90_inq_dimid(ncid, "level", varid))
    call check_nc(nf90_inquire_dimension(ncid, varid, len=dimlen))
    if (dimlen /= 3) call mpp_error(FATAL, "test: incorrect non-spatial dimension size")
    call check_nc(nf90_inq_dimid(ncid, "member", varid))
    call check_nc(nf90_inquire_dimension(ncid, varid, len=dimlen))
    if (dimlen /= 2) call mpp_error(FATAL, "test: incorrect blank-cartesian dimension size")
    call check_nc(nf90_inq_dimid(ncid, "time", varid))
    call check_nc(nf90_inquire(ncid, unlimiteddimid=unlimited_dimid))
    if (varid /= unlimited_dimid) call mpp_error(FATAL, "test: time dimension is not unlimited")
    call check_nc(nf90_inquire_dimension(ncid, varid, len=dimlen))
    if (dimlen /= 3) call mpp_error(FATAL, "test: incorrect unlimited dimension length")
    call check_nc(nf90_inq_varid(ncid, "history", varid))
    call check_nc(nf90_inquire_variable(ncid, varid, ndims=ndims))
    if (ndims /= 3) call mpp_error(FATAL, "test: incorrect record field rank")
    call check_nc(nf90_inquire_variable(ncid, varid, dimids=record_dimids))
    if (any(record_dimids(:2) /= dimids) .or. record_dimids(3) /= unlimited_dimid) &
      call mpp_error(FATAL, "test: incorrect record field dimension order")
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
    call check_nc(nf90_open(second_filename, nf90_nowrite, ncid))
    call check_nc(nf90_inq_varid(ncid, "independent", varid))
    allocate(actual(nx,ny))
    call check_nc(nf90_get_var(ncid, varid, actual))
    do j = 1, ny
      do i = 1, nx
        if (actual(i,j) /= real(3000 + 100*j + i)) &
          call mpp_error(FATAL, "test: independent file data mismatch")
      enddo
    enddo
    call check_nc(nf90_close(ncid))
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
    do step = 1, 3
      call read_data(fileobj, "history", actual, unlim_dim_level=step)
      do j = 1, ny
        do i = 1, nx
          if (step == 1) then
            if (actual(i,j) /= real(-100*j - i)) &
              call mpp_error(FATAL, "test: default record level data mismatch")
          else
            if (actual(i,j) /= real(1000*(4-step) + 100*j + i)) &
              call mpp_error(FATAL, "test: r8 record data mismatch")
          endif
        enddo
      enddo
      call read_data(fileobj, "history_r4", actual, unlim_dim_level=step)
      do j = 1, ny
        do i = 1, nx
          if (actual(i,j) /= real(2000*(4-step) + 100*j + i)) &
            call mpp_error(FATAL, "test: r4 record data mismatch")
        enddo
      enddo
    enddo
    call generate_fields(.false.)
    block
      character(len=12) :: label, labels(2)
      call read_data(fileobj, "label", label)
      if (label /= "offloaded") call mpp_error(FATAL, "test: scalar character data mismatch")
      call read_data(fileobj, "labels", labels)
      if (any(labels /= [character(len=12) :: "first", "second"])) &
        call mpp_error(FATAL, "test: vector character data mismatch")
    end block
    call close_file(fileobj)
  endif
  call mpp_set_current_pelist(all_pes, no_sync=.true.)
  call fms_end()
contains

  !> Generate the data to write, and either write it to the offloaded file or check an existing file for validity
  subroutine generate_fields(writing)
    logical, intent(in) :: writing
    character(len=6), parameter :: storage_types(4) = [character(len=6) :: "float", "double", "int", "int64"]
    integer :: kind_index, data_rank, record, dims(4), xstart, ystart
    character(len=9) :: field_axes(5)
    integer :: ix, iy, iz, member
    integer, allocatable :: level
    character(len=32) :: field_name
    class(*), allocatable :: values(:,:,:,:)
    real(kind=r8_kind), allocatable :: expected(:,:,:,:)

    do kind_index = 1, 4
      do data_rank = 0, 4
        write(field_name, '(a,i0,a,i0)') "kind", kind_index, "_rank", data_rank
        dims = 1
        xstart = 1
        ystart = 1
        if (data_rank == 1) dims(1) = 3
        if (data_rank >= 2) then
          dims(:2) = [nx,ny]
          if (writing) then
            dims(:2) = [ie-is+1,je-js+1]
            xstart = is
            ystart = js
          endif
        endif
        if (data_rank >= 3) dims(3) = 3
        if (data_rank == 4) dims(4) = 2
        select case (kind_index)
        case (1)
          allocate(real(kind=r4_kind) :: values(dims(1),dims(2),dims(3),dims(4)))
        case (2)
          allocate(real(kind=r8_kind) :: values(dims(1),dims(2),dims(3),dims(4)))
        case (3)
          allocate(integer(kind=i4_kind) :: values(dims(1),dims(2),dims(3),dims(4)))
        case (4)
          allocate(integer(kind=i8_kind) :: values(dims(1),dims(2),dims(3),dims(4)))
        end select
        allocate(expected(dims(1),dims(2),dims(3),dims(4)))
        field_axes = ""
        if (data_rank == 1) field_axes(1) = "level"
        if (data_rank >= 2) field_axes(:2) = axes(2:1:-1)
        if (data_rank >= 3) field_axes(3) = "level"
        if (data_rank == 4) field_axes(4) = "member"
        field_axes(data_rank+1) = "time"
        if (writing) &
          call register_field_offloaded_file(field_name, offload_file, &
                                             trim(storage_types(kind_index)), field_axes(:data_rank+1))
        do record = 2, 0, -1
          if (record > 0) then
            if (.not. allocated(level)) allocate(level)
            level = record
          else
            deallocate(level)
            field_name = trim(field_name)//"_static"
            if (writing) then
              if (data_rank == 0) then
                call register_field_offloaded_file(field_name, offload_file, trim(storage_types(kind_index)))
              else
                call register_field_offloaded_file(field_name, offload_file, &
                                                   trim(storage_types(kind_index)), field_axes(:data_rank))
              endif
            endif
          endif
          do member = 1, dims(4)
            do iz = 1, dims(3)
              do iy = 1, dims(2)
                do ix = 1, dims(1)
                  expected(ix,iy,iz,member) = 100000*record + 10000*member + 1000*iz + &
                                             100*(ystart+iy-1) + xstart+ix-1
                enddo
              enddo
            enddo
          enddo
          if (kind_index == 2) expected = expected + 0.1234567890123_r8_kind
          if (writing) then
            select type (values)
            type is (real(kind=r4_kind))
              values = real(expected,r4_kind)
            type is (real(kind=r8_kind))
              values = expected
            type is (integer(kind=i4_kind))
              values = int(expected,i4_kind)
            type is (integer(kind=i8_kind))
              values = int(expected,i8_kind) + 2_i8_kind**60
            end select
            select case (data_rank)
            case (0)
              call write_field_offloaded_file(field_name, values(1,1,1,1), domain, offload_file, level)
            case (1)
              call write_field_offloaded_file(field_name, values(:,1,1,1), domain, offload_file, level)
            case (2)
              call write_field_offloaded_file(field_name, values(:,:,1,1), domain, offload_file, level)
            case (3)
              call write_field_offloaded_file(field_name, values(:,:,:,1), domain, offload_file, level)
            case (4)
              call write_field_offloaded_file(field_name, values, domain, offload_file, level)
            end select
          else
            select case (data_rank)
            case (0)
              call read_data(fileobj, field_name, values(1,1,1,1), unlim_dim_level=level)
            case (1)
              call read_data(fileobj, field_name, values(:,1,1,1), unlim_dim_level=level)
            case (2)
              call read_data(fileobj, field_name, values(:,:,1,1), unlim_dim_level=level)
            case (3)
              call read_data(fileobj, field_name, values(:,:,:,1), unlim_dim_level=level)
            case (4)
              call read_data(fileobj, field_name, values, unlim_dim_level=level)
            end select
            select type (values)
            type is (real(kind=r4_kind))
              if (any(values /= real(expected,r4_kind))) call mpp_error(FATAL, "test: r4 mismatch in "//field_name)
            type is (real(kind=r8_kind))
              if (any(values /= expected)) call mpp_error(FATAL, "test: r8 mismatch in "//field_name)
            type is (integer(kind=i4_kind))
              if (any(values /= int(expected,i4_kind))) call mpp_error(FATAL, "test: i4 mismatch in "//field_name)
            type is (integer(kind=i8_kind))
              if (any(values /= int(expected,i8_kind) + 2_i8_kind**60)) &
                call mpp_error(FATAL, "test: i8 mismatch in "//field_name)
            end select
          endif
        enddo
        deallocate(values, expected)
      enddo
    enddo
  end subroutine

  subroutine check_nc(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) call mpp_error(FATAL, "test: "//trim(nf90_strerror(status)))
  end subroutine
end program test_io_offloading
