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
module offloading_io_mod
  use mpp_mod
  use mpp_domains_mod
  use fms_string_utils_mod, only: string
  use fms2_io_mod
  use metadata_transfer_mod
  use platform_mod
  use fms_mod

  implicit none

  integer :: current_files_init !< number of currently initialized offloading files
  logical :: module_is_initialized !< .true. if module has been initialized

  integer :: max_files = 10 !< amount of offloaded files to allocate space for

  namelist / offloading_io_nml / max_files

  !> Structure to hold offloading file information
  type :: offloading_obj_out
    integer :: id
    character(len=:), allocatable :: filename !< filename of the offloaded netcdf file
    class(FmsNetcdfFile_t), allocatable :: fileobj !< fms2_io file object
    type(domain2D) :: domain_out !< domain on offloading PEs
  end type

  !> Offload equivalent of register_axis in fms2_io_mod
  !! Registers an axis to a netcdf file on offloaded PEs. File must have been opened with open_file_offload.
  !! TODO: add register_unstructured_axis_offload for the unstructured grid
  interface register_axis_offload
    procedure :: register_netcdf_axis_offload
    procedure :: register_domain_axis_offload
  end interface

  !> Offload equivalent of write_data in fms2_io_mod
  !! Writes data to a netcdf file on offloaded PEs. File must have been opened with open_file_offload.
  interface write_data_offload
    procedure :: write_data_offload_2d
    procedure :: write_data_offload_3d
  end interface

  !> Array of offloading objects that have been initialized by this module.
  type(offloading_obj_out), allocatable, target :: offloading_objs(:)

  private

  public :: offloading_io_init, open_file_offload
  public :: global_metadata_offload, close_file_offload, register_axis_offload, register_field_offload
  public :: write_data_offload
  public :: create_cubic_domain, create_lat_lon_domain

  contains

  !> Initialize by allocating array used to keep track of offloading file objects.
  subroutine offloading_io_init()
    integer :: ierr, io
    if (module_is_initialized) return
    current_files_init = 0
    allocate(offloading_objs(max_files))
    module_is_initialized = .true.
    read (input_nml_file, offloading_io_nml, iostat=io)
    ierr = check_nml_error(io,'offloading_io_nml')
  end subroutine offloading_io_init

  !> Open a netcdf file and set it up for offloaded writes
  !! This routine should be called from both the model PEs and the offload PEs, with the full list
  !! of pes for each group being provided. The model PEs will broadcast the filename and domain.
  subroutine open_file_offload(fileobj, filename, domain_in, pe_in, pe_out)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj !< fms2_io file object
    character(len=*),       intent(in)    :: filename !< filename to open
    type(domain2D),         intent(inout) :: domain_in !< model domain (from model pes)
    integer,                intent(in)    :: pe_in(:) !< model pes
    integer,                intent(in)    :: pe_out(:) !< offload pes

    integer, parameter :: str_len = 255
    character(len=str_len) :: filename_out(1)
    integer :: object_id
    integer :: global_domain_size(2)
    integer :: ntile
    integer, allocatable :: all_current_pes(:)
    integer, allocatable :: broadcasting_pes(:)
    logical :: is_pe_out

    is_pe_out = ANY(pe_out .eq. mpp_pe())

    ! This should be called from the model PEs and the offload PEs
    if (.not. module_is_initialized) &
      call mpp_error(FATAL, "offloading_io_mod is not initialized")

    allocate(all_current_pes(mpp_npes()))
    call mpp_get_current_pelist(all_current_pes)

    filename_out(1) = ""
    if (mpp_pe() .eq. pe_in(1)) then
      !< The root model pe gets the domain info (ntiles and size of global domain)
      ntile = mpp_get_ntile_count(domain_in)

      !< The number of tiles must be the same as the number of offloading pes
      if ( MOD(size(pe_out), ntile) .ne. 0 ) &
        call mpp_error(FATAL, "The number of offloading PEs must be the same as the number of tiles of the domain")
      filename_out(1) = filename

      call mpp_get_global_domain(domain_in, xsize=global_domain_size(1), ysize=global_domain_size(2))
    endif

    ! A "Root" model PE broadcasts the filename, domain size, and number of tiles to the offload pes
    if (mpp_pe() .eq. pe_in(1) .or. is_pe_out) then
      allocate(broadcasting_pes(1 + size(pe_out)))
      broadcasting_pes(1) = pe_in(1) ! root pe
      broadcasting_pes(2:size(broadcasting_pes)) = pe_out ! offload pes
      call mpp_set_current_pelist( broadcasting_pes )
      ! TODO bundle these into a single derived type to reduce the number of broadcasts
      call mpp_broadcast(filename_out, str_len, pe_in(1))
      call mpp_broadcast(global_domain_size, size(global_domain_size), pe_in(1))
      call mpp_broadcast(ntile, pe_in(1))
    endif

    ! Broadcast the domain
    call mpp_set_current_pelist(all_current_pes)
    if (is_pe_out) call mpp_define_null_domain(domain_in)
    call mpp_broadcast_domain(domain_in)

    ! The offload pes inits the offloading object and return an object id
    if (is_pe_out) then
      call mpp_set_current_pelist(pe_out)
      object_id = init_offloading_object(filename_out(1), global_domain_size(1), global_domain_size(2),&
         ntile)
    endif
    call mpp_set_current_pelist(all_current_pes)
    call mpp_broadcast(object_id, pe_out(1))

    ! Init the "offloading object" in the fileobj
    call fileobj%offloading_obj_in%init(object_id, pe_out, pe_in, domain_in)
  end subroutine open_file_offload

  !> Broadcast and register a global metadata attribute on offloading PEs
  subroutine global_metadata_offload(fileobj, attribute_name, attribute_value)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj !< fms2_io file object
    character(len=*), intent(in) :: attribute_name !< name of the global attribute to register
    class(*), intent(in) :: attribute_value !< value of the global attribute to register (r4, r8, i4, i8, str)

    integer :: id
    type(offloading_obj_out), pointer :: this

    !TODO better PEs management!
    integer, allocatable :: offloading_pes(:)
    integer, allocatable :: model_pes(:)
    integer, allocatable :: all_current_pes(:)
    integer, allocatable :: broadcasting_pes(:)
    logical :: is_model_pe
    character(len=255) :: att_name(1)

    integer :: int_buf
    real(r4_kind), allocatable :: r4_tmp(:)
    real(r8_kind), allocatable :: r8_tmp(:)
    integer(i4_kind) , allocatable :: i4_tmp(:)
    integer(i8_kind) , allocatable :: i8_tmp(:)
    character(len=:), allocatable :: str_tmp
    class(metadata_class), allocatable :: transfer_obj

    select type (attribute_value)
    type is (real(kind=r8_kind))
      allocate(metadata_r8_type :: transfer_obj)
      call transfer_obj%fms_metadata_transfer_init(real8_type)
    type is (real(kind=r4_kind))
      allocate(metadata_r4_type :: transfer_obj)
      call transfer_obj%fms_metadata_transfer_init(real4_type)
    type is (integer(kind=i4_kind))
      allocate(metadata_i4_type :: transfer_obj)
      call transfer_obj%fms_metadata_transfer_init(int4_type)
    type is (integer(kind=i8_kind))
      allocate(metadata_i8_type :: transfer_obj)
      call transfer_obj%fms_metadata_transfer_init(int8_type)
    type is (character(*))
      allocate(metadata_str_type :: transfer_obj)
      call transfer_obj%fms_metadata_transfer_init(str_type)
    class default
      call mpp_error(FATAL, "Unsupported attribute type for offloading: " // string(attribute_value))
    end select

    offloading_pes = fileobj%offloading_obj_in%offloading_pes
    model_pes = fileobj%offloading_obj_in%model_pes
    is_model_pe = fileobj%offloading_obj_in%is_model_pe

    id = fileobj%offloading_obj_in%id
    this => offloading_objs(id)

    if (is_model_pe) then
      att_name(1) = attribute_name
      call transfer_obj%set_attribute_name(attribute_name)
    endif

    allocate(all_current_pes(mpp_npes()))
    call mpp_get_current_pelist(all_current_pes)

    if (mpp_pe() .eq. model_pes(1) .or. .not. is_model_pe) then

      allocate(broadcasting_pes(1 + size(offloading_pes)))
      broadcasting_pes(1) = model_pes(1)
      broadcasting_pes(2:size(broadcasting_pes)) = offloading_pes
      call mpp_set_current_pelist( broadcasting_pes )

      select type (attribute_value)
        type is (real(kind=r4_kind))
          ! TODO replace this mess with a single call if possible
          if (is_model_pe) then
            select type(transfer_obj)
              type is (metadata_r4_type)
                call transfer_obj%set_attribute_value([attribute_value])
            end select
          endif
          call transfer_obj%fms_metadata_broadcast()
          select type(transfer_obj)
            type is (metadata_r4_type)
              r4_tmp = transfer_obj%get_attribute_value()
          end select
          att_name(1) = transfer_obj%get_attribute_name()
          if (.not. is_model_pe) then
            call register_global_attribute(this%fileobj, att_name(1), r4_tmp)
          endif

        type is (real(kind=r8_kind))
          ! TODO replace this mess with a single call if possible
          if (is_model_pe) then
            select type(transfer_obj)
              type is (metadata_r8_type)
                call transfer_obj%set_attribute_value([attribute_value])
            end select
          endif
          call transfer_obj%fms_metadata_broadcast()
          select type(transfer_obj)
            type is (metadata_r8_type)
              r8_tmp = transfer_obj%get_attribute_value()
          end select
          att_name(1) = transfer_obj%get_attribute_name()
          if (.not. is_model_pe) then
            call register_global_attribute(this%fileobj, att_name(1), r8_tmp)
          endif

        type is (integer(kind=i4_kind))
          if (is_model_pe) then
            select type(transfer_obj)
              type is (metadata_i4_type)
                call transfer_obj%set_attribute_value([attribute_value])
                int_buf = attribute_value
            end select
          endif
          call transfer_obj%fms_metadata_broadcast()
          select type(transfer_obj)
            type is (metadata_i4_type)
              i4_tmp = transfer_obj%get_attribute_value()
          end select
          att_name(1) = transfer_obj%get_attribute_name()
          if (.not. is_model_pe) then
            call register_global_attribute(this%fileobj, att_name(1), i4_tmp)
          endif

        type is (integer(kind=i8_kind))
          ! TODO replace this mess with a single call if possible
          if (is_model_pe) then
            select type(transfer_obj)
              type is (metadata_i8_type)
                call transfer_obj%set_attribute_value([attribute_value])
            end select
          endif
          call transfer_obj%fms_metadata_broadcast()
          select type(transfer_obj)
            type is (metadata_i8_type)
              i8_tmp = transfer_obj%get_attribute_value()
          end select
          att_name(1) = transfer_obj%get_attribute_name()
          if (.not. is_model_pe) then
            call register_global_attribute(this%fileobj, att_name(1), i8_tmp)
          endif

        type is (character(len=*))
          ! TODO replace this mess with a single call if possible
          if (is_model_pe) then
            select type(transfer_obj)
              type is (metadata_str_type)
                call transfer_obj%set_attribute_value(attribute_value)
            end select
          endif
          call transfer_obj%fms_metadata_broadcast()
          select type(transfer_obj)
            type is (metadata_str_type)
              str_tmp = transfer_obj%get_attribute_value()
          end select
          att_name(1) = transfer_obj%get_attribute_name()
          if (.not. is_model_pe) then
            call register_global_attribute(this%fileobj, att_name(1), str_tmp)
          endif

      end select
    endif

    call mpp_set_current_pelist(all_current_pes)
  end subroutine

  !> Register a domain axis (ie. x or y) on offloading PEs
  subroutine register_domain_axis_offload(fileobj, axis_name, cart)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj !< fms2_io file object
    character(len=*), intent(in) :: axis_name !< axis name to be written to file
    character(len=1), intent(in) :: cart !< must be either 'x' or 'y' for cartesian axis

    integer :: id
    type(offloading_obj_out), pointer :: this

    integer, allocatable :: offloading_pes(:)
    integer, allocatable :: model_pes(:)
    integer, allocatable :: all_current_pes(:)
    integer, allocatable :: broadcasting_pes(:)
    logical :: is_model_pe
    character(len=ATTR_NAME_MAX_LENGTH) :: var_info(2)

    offloading_pes = fileobj%offloading_obj_in%offloading_pes
    model_pes = fileobj%offloading_obj_in%model_pes
    is_model_pe = fileobj%offloading_obj_in%is_model_pe

    id = fileobj%offloading_obj_in%id
    this => offloading_objs(id)

    if (is_model_pe) then
      var_info(1) = axis_name
      var_info(2) = cart
    endif

    allocate(all_current_pes(mpp_npes()))
    call mpp_get_current_pelist(all_current_pes)

    if (mpp_pe() .eq. model_pes(1) .or. .not. is_model_pe) then
      allocate(broadcasting_pes(1 + size(offloading_pes)))
      broadcasting_pes(1) = model_pes(1)
      broadcasting_pes(2:size(broadcasting_pes)) = offloading_pes
      call mpp_set_current_pelist( broadcasting_pes )
      call mpp_broadcast(var_info, 255, model_pes(1))
    endif

    if (.not. is_model_pe) then
      select type(file=>this%fileobj)
        type is(FmsNetcdfDomainFile_t)
          call register_axis(file, trim(var_info(1)), trim(var_info(2)))
        class default
          call mpp_error(FATAL, "offloading_io_mod::register_domain_axis_offload currently only supports FmsNetcdfDomainFile_t")
      end select
    endif

    call mpp_set_current_pelist(all_current_pes)
  end subroutine register_domain_axis_offload

  !> Register a netcdf axis on offloading PEs
  subroutine register_netcdf_axis_offload(fileobj, axis_name, length)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj !< fms2_io file object
    character(len=*), intent(in) :: axis_name !< axis name to be written to file
    integer, intent(in) :: length !< length of the axis

    integer :: id
    type(offloading_obj_out), pointer :: this

    !TODO better PEs management!
    integer, allocatable :: offloading_pes(:)
    integer, allocatable :: model_pes(:)
    integer, allocatable :: all_current_pes(:)
    integer, allocatable :: broadcasting_pes(:)
    logical :: is_model_pe
    character(len=255) :: var_axis(1)
    integer :: var_length
    integer :: axis_length

    offloading_pes = fileobj%offloading_obj_in%offloading_pes
    model_pes = fileobj%offloading_obj_in%model_pes
    is_model_pe = fileobj%offloading_obj_in%is_model_pe

    id = fileobj%offloading_obj_in%id
    this => offloading_objs(id)

    ! get var data on root for broadcasting to offload pes
    if (mpp_pe() .eq. model_pes(1)) then
      var_axis(1) = trim(axis_name)
      axis_length = len_trim(var_axis(1))
      var_length = length
    endif

    allocate(all_current_pes(mpp_npes()))
    call mpp_get_current_pelist(all_current_pes)

    ! root pe broadcasts the axis name and length to offload pes
    if (mpp_pe() .eq. model_pes(1) .or. .not. is_model_pe) then
      allocate(broadcasting_pes(1 + size(offloading_pes)))
      broadcasting_pes(1) = model_pes(1)
      broadcasting_pes(2:size(broadcasting_pes)) = offloading_pes
      call mpp_set_current_pelist( broadcasting_pes )
      ! TODO bundle these into a single derived type to reduce the number of broadcasts
      call mpp_broadcast(axis_length, model_pes(1))
      call mpp_broadcast(var_axis, axis_length, model_pes(1))
      call mpp_broadcast(var_length, model_pes(1))
    endif

    if (.not. is_model_pe) then
      select type(wut=>this%fileobj)
        type is(FmsNetcdfDomainFile_t)
          call register_axis(this%fileobj, var_axis(1)(1:axis_length), var_length)
        class default
          call mpp_error(FATAL, "offloading_io_mod::register_netcdf_axis_offload currently only supports FmsNetcdfDomainFile_t")
      end select
    endif

    call mpp_set_current_pelist(all_current_pes)
  end subroutine register_netcdf_axis_offload

  !> Register a netcdf field on offloading PEs
  subroutine register_field_offload(fileobj, varname, vartype, dimensions)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj !> fms2_io file object
    character(len=*), intent(in) :: varname !> name of the variable to be registered
    character(len=*), intent(in) :: vartype !> type of the variable to be registered
                                            !! must be one of {r4_type, r8_type, i4_type, i8_type, str_type}
    character(len=*), intent(in) :: dimensions(:) !< previously registered dimension/axis names for the variable

    integer :: id
    type(offloading_obj_out), pointer :: this
    integer :: ndim

    !TODO better PEs management!
    integer, allocatable :: offloading_pes(:)
    integer, allocatable :: model_pes(:)
    integer, allocatable :: all_current_pes(:)
    integer, allocatable :: broadcasting_pes(:)
    logical :: is_model_pe
    character(len=255), allocatable :: var_info(:)

    offloading_pes = fileobj%offloading_obj_in%offloading_pes
    model_pes = fileobj%offloading_obj_in%model_pes
    is_model_pe = fileobj%offloading_obj_in%is_model_pe

    id = fileobj%offloading_obj_in%id
    this => offloading_objs(id)

    allocate(all_current_pes(mpp_npes()))
    call mpp_get_current_pelist(all_current_pes)

    if (is_model_pe) then
      ndim = size(dimensions)

      allocate(var_info(ndim + 2))
      var_info(1) = varname
      var_info(2) = vartype
      var_info(3:) = dimensions

    endif

    if (mpp_pe() .eq. model_pes(1) .or. .not. is_model_pe) then
      allocate(broadcasting_pes(1 + size(offloading_pes)))
      broadcasting_pes(1) = model_pes(1)
      broadcasting_pes(2:size(broadcasting_pes)) = offloading_pes
      call mpp_set_current_pelist( broadcasting_pes )
      call mpp_broadcast(ndim, model_pes(1))

      if (.not. is_model_pe) allocate(var_info(ndim + 2))
      call mpp_broadcast(var_info, 255, model_pes(1))
    endif

    if (.not. is_model_pe) then
      !select type(wut=>this%fileobj)
      !  type is(FmsNetcdfDomainFile_t)
          call register_field(this%fileobj, trim(var_info(1)), trim(var_info(2)), var_info(3:))
      !end select
    endif

    call mpp_set_current_pelist(all_current_pes)
  end subroutine

  !> Write 3D data to offloaded netcdf file
  subroutine write_data_offload_3d(fileobj, varname, vardata, unlim_dim_level)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj !< fms2_io file object
    character(len=*), intent(in) :: varname !< name of the variable to be written
    real(kind=r4_kind), intent(in) :: vardata(:,:,:) !< 3D data to be written
    integer, intent(in), optional :: unlim_dim_level !< level along unlimited dimension to write to

    integer :: id
    type(offloading_obj_out), pointer :: this

    !TODO better PEs management!
    integer, allocatable :: offloading_pes(:)
    integer, allocatable :: model_pes(:)
    integer, allocatable :: all_current_pes(:)
    logical :: is_model_pe

    real(kind=r4_kind), allocatable :: var_r4_data(:,:,:)
    type(domain2D) :: domain_out
    type(domain2D) :: domain_in
    integer :: isc, iec, jsc, jec, nz, redistribute_clock
    character(len=ATTR_NAME_MAX_LENGTH) :: varname_tmp(1)

    offloading_pes = fileobj%offloading_obj_in%offloading_pes
    model_pes = fileobj%offloading_obj_in%model_pes
    is_model_pe = fileobj%offloading_obj_in%is_model_pe

    id = fileobj%offloading_obj_in%id
    this => offloading_objs(id)

    allocate(all_current_pes(mpp_npes()))
    call mpp_get_current_pelist(all_current_pes)

    redistribute_clock = mpp_clock_id( 'data_transfer' )
    call mpp_clock_begin(redistribute_clock)

    nz = size(vardata, 3)
    call mpp_broadcast(nz, model_pes(1))

    !Allocate space to store the data!
    if (.not. is_model_pe) then
      domain_out = this%domain_out
      call mpp_get_data_domain(domain_out, isc, iec, jsc, jec)
      allocate(var_r4_data(isc:iec, jsc:jec, nz))
      call mpp_define_null_domain(domain_in)
    else
      domain_in = fileobj%offloading_obj_in%domain_in
      call mpp_define_null_domain(domain_out)
    endif

    ! get domain from the other pes
    call mpp_broadcast_domain(domain_out)
    call mpp_broadcast_domain(domain_in)

    call mpp_redistribute( domain_in, vardata, domain_out, var_r4_data)
    call mpp_redistribute( domain_in, vardata, domain_out, var_r4_data, free=.true.)

    if(mpp_pe() .eq. model_pes(1)) then
      varname_tmp(1) = trim(varname)
    endif
    call mpp_broadcast(varname_tmp, 255, model_pes(1))

    call mpp_clock_end(redistribute_clock)

      if (.not. is_model_pe) then
        select type(wut=>this%fileobj)
          type is(FmsNetcdfDomainFile_t)
            if (present(unlim_dim_level)) then
              call write_data(wut, varname, var_r4_data, unlim_dim_level=unlim_dim_level)
            else
              call write_data(wut, varname, var_r4_data)
            endif
        end select
      endif
    call mpp_set_current_pelist(all_current_pes)
  end subroutine write_data_offload_3d

  !> Write 2D data to offloaded netcdf file
  subroutine write_data_offload_2d(fileobj, varname, vardata)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj !< fms2_io file object
    character(len=*), intent(in) :: varname !< name of the variable to be written
    real(kind=r4_kind), intent(in) :: vardata(:,:) !< 2D data to be written

    integer :: id
    type(offloading_obj_out), pointer :: this

    !TODO better PEs management!
    integer, allocatable :: offloading_pes(:)
    integer, allocatable :: model_pes(:)
    integer, allocatable :: all_current_pes(:)
    logical :: is_model_pe

    real(kind=r4_kind), allocatable :: var_r4_data(:,:)
    type(domain2D) :: domain_out
    type(domain2D) :: domain_in
    integer :: isc, iec, jsc, jec
    character(len=255) :: varname_tmp(1)

    offloading_pes = fileobj%offloading_obj_in%offloading_pes
    model_pes = fileobj%offloading_obj_in%model_pes
    is_model_pe = fileobj%offloading_obj_in%is_model_pe

    id = fileobj%offloading_obj_in%id
    this => offloading_objs(id)

    allocate(all_current_pes(mpp_npes()))
    call mpp_get_current_pelist(all_current_pes)

    !Allocate space to store the data!
    if (.not. is_model_pe) then
      domain_out = this%domain_out
      call mpp_get_compute_domain(domain_out, isc, iec, jsc, jec)
      allocate(var_r4_data(isc:iec, jsc:jec))
      call mpp_define_null_domain(domain_in)
    else
      domain_in = fileobj%offloading_obj_in%domain_in
      call mpp_define_null_domain(domain_out)
    endif

    call mpp_broadcast_domain(domain_out)
    call mpp_broadcast_domain(domain_in)

    ! redistribute data from model domain to offload domain and then free memory for future calls
    call mpp_redistribute( domain_in, vardata, domain_out, var_r4_data)
    call mpp_redistribute( domain_in, vardata, domain_out, var_r4_data, free=.true.)

    ! broadcast the variable name
    if(mpp_pe() .eq. model_pes(1)) then
      varname_tmp(1) = trim(varname)
    endif
    call mpp_broadcast(varname_tmp, 255, model_pes(1))

    if (.not. is_model_pe) then
      select type(wut=>this%fileobj)
        type is(FmsNetcdfDomainFile_t)
          call write_data(wut, varname_tmp(1), var_r4_data)
      end select
    endif
    call mpp_set_current_pelist(all_current_pes)
  end subroutine write_data_offload_2d

  !> Close offloaded netcdf file
  subroutine close_file_offload(fileobj)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj !> fms2_io file object to close

    integer :: id
    type(offloading_obj_out), pointer :: this
    logical :: is_model_pe

    id = fileobj%offloading_obj_in%id
    this => offloading_objs(id)
    is_model_pe = fileobj%offloading_obj_in%is_model_pe

    if (.not. is_model_pe) call close_file(this%fileobj)
  end subroutine close_file_offload

  !> Initialize an offloading object on offload PEs
  integer function init_offloading_object(filename, nx, ny, ntile)
    character(len=*), intent(in) :: filename !< filename to open
    integer, intent(in) :: nx !< x size of the global domain
    integer, intent(in) :: ny !< y size of the global domain
    integer, intent(in) :: ntile !< number of tiles, only supports 1 for lat-lon or 6 for cubed-sphere

    type(offloading_obj_out), pointer :: this
    integer, allocatable :: curr_pelist(:)

    current_files_init = current_files_init + 1
    if (current_files_init .gt. max_files) &
      call mpp_error(FATAL, "The number of files is too large")

    ! An array of offloading objects is stored in this module, the "id" is the index of the object in the array
    this => offloading_objs(current_files_init)
    this%id = current_files_init
    this%filename = trim(filename)

    ! does npes return current pelist or total??
    allocate(curr_pelist(mpp_npes()))
    call mpp_get_current_pelist(curr_pelist)

    select case (ntile)
    case (1)
      this%domain_out = create_lat_lon_domain(nx, ny)
    case (6)
      this%domain_out = create_cubic_domain(nx, ny, ntile, (/1,1/), offload_pes=curr_pelist)
    case default
      call mpp_error(FATAL, "Unsupported number of tiles for offloading: " // trim(adjustl(string(ntile))))
    end select

    allocate(FmsNetcdfDomainFile_t :: this%fileobj)
    select type (fileobj => this%fileobj)
    type is (FmsNetcdfDomainFile_t)
      if ( .not. open_file(fileobj, trim(this%filename), "overwrite", this%domain_out)) &
        call mpp_error(FATAL, "Error opening file")
    end select

    init_offloading_object = current_files_init
  end function

  !! TODO move this somewhere else
  function create_lat_lon_domain(nx_in, ny_in, halox, haloy ) &
    result(domain_out)
    integer, intent(in) :: nx_in !< number of lat-lon grid points in x direction
    integer, intent(in) :: ny_in !< number of lat-lon grid points in y direction
    integer, intent(in), optional :: halox !< number of halo points in x direction
    integer, intent(in), optional :: haloy !< number of halo points in y direction
    type(domain2d) :: domain_out
    integer :: layout(2)

    call mpp_define_layout( (/1,nx_in,1,ny_in/), mpp_npes(), layout )
    call mpp_define_domains( (/1,nx_in,1,ny_in/), layout, domain_out, xhalo=halox, yhalo=haloy)
    call mpp_define_io_domain(domain_out, (/1,1/))
  end function create_lat_lon_domain

  !! TODO move this somewhere else
  function create_cubic_domain(nx_in, ny_in, ntiles, io_layout, nhalos, offload_pes, layout) &
    result(domain_out)
    integer, intent(in) :: nx_in !< number of grid points in x direction per tile
    integer, intent(in) :: ny_in !< number of grid points in y direction per tile
    integer, intent(in) :: ntiles !< number of tiles, must be 6 for cubed-sphere
    integer, intent(in) :: io_layout(2) !< layout for I/O operations
    integer, intent(in), optional :: nhalos !< number of halo points
    integer, intent(in), optional :: offload_pes(:) !< list of PEs used for offloading write operations
    integer, optional :: layout(2) !< layout to be used for each tile)

    type(domain2d) :: domain_out

    integer :: npes
    integer :: npes_per_tile
    integer :: layout_tmp(2)
    integer, allocatable :: global_indices(:,:)
    integer, allocatable :: layout2D(:,:)
    integer, allocatable :: pe_start(:)
    integer, allocatable :: pe_end(:)
    integer :: n

    npes = mpp_npes()
    if( mod(npes, ntiles) .NE. 0 ) call mpp_error(FATAL, &
          "create_cubic_domain: npes is not divisible by ntiles")

    npes_per_tile = npes/ntiles
    if( .not. present(layout)) then
      call mpp_define_layout ((/1,nx_in,1,ny_in/), npes_per_tile, layout_tmp )
    else
      layout_tmp = layout
    endif

    allocate(global_indices(4, ntiles))
    allocate(layout2D(2, ntiles))
    allocate(pe_start(ntiles), pe_end(ntiles))

    do n = 1, ntiles
      global_indices(:,n) = (/1,nx_in,1,ny_in/)
      layout2D(:,n) = layout_tmp
      if( present(offload_pes)) then
        pe_start(n) = offload_pes((n-1)*npes_per_tile+1)
        pe_end(n) = offload_pes((n)*npes_per_tile)
      else
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n) = n*npes_per_tile-1
      endif
    end do

    print *, "pe: ", mpp_pe(), "creating mosaic with pe_start", pe_start, "pe_end", pe_end

    call define_cubic_mosaic(domain_out, (/nx_in,nx_in,nx_in,nx_in,nx_in,nx_in/), &
                            (/ny_in,ny_in,ny_in,ny_in,ny_in,ny_in/), &
                              global_indices, layout2D, pe_start, pe_end, io_layout, nhalos )
  end function create_cubic_domain

  !> @brief Initialize a cubed-sphere atomsphere domain.
  !! TODO move this somehere else
  subroutine define_cubic_mosaic(domain, ni, nj, global_indices, layout, pe_start, pe_end, &
                                      io_layout, nhalos)

    integer, dimension(:), intent(in) :: ni !< number of grid points in i direction per tile
    integer, dimension(:), intent(in) :: nj !< number of grid points in j direction per tile
    integer, dimension(:,:), intent(in) :: global_indices !< global indices for each tile
    integer, dimension(:,:), intent(in) :: layout !< array of layouts for each tile
    integer, dimension(:), intent(in) :: pe_start !< starting PE for each tile
    integer, dimension(:), intent(in) :: pe_end !< ending PE for each tile
    integer, dimension(2), intent(in) :: io_layout !< layout for I/O operations
    type(domain2d), intent(inout) :: domain !< A cubed-sphere domain.
    integer, optional, intent(in) :: nhalos !< number of halo points

    integer, dimension(12) :: tile1
    integer, dimension(12) :: tile2
    integer, dimension(12) :: istart1
    integer, dimension(12) :: iend1
    integer, dimension(12) :: jstart1
    integer, dimension(12) :: jend1
    integer, dimension(12) :: istart2
    integer, dimension(12) :: iend2
    integer, dimension(12) :: jstart2
    integer, dimension(12) :: jend2
    integer :: ntiles
    integer :: num_contact
    integer, dimension(2) :: msize
    integer :: whalo
    integer :: ehalo
    integer :: shalo
    integer :: nhalo

    ntiles = 6
    num_contact = 12

    whalo = 2
    if (present(nhalos)) whalo = nhalos
    ehalo = whalo
    shalo = whalo
    nhalo = whalo

    if (size(pe_start) .ne. 6 .or. size(pe_end) .ne. 6 ) then
      call mpp_error(FATAL, "size of pe_start and pe_end should be 6.")
    endif
    if (size(global_indices,1) .ne. 4) then
      call mpp_error(FATAL, "size of first dimension of global_indices should be 4.")
    endif
    if (size(global_indices,2) .ne. 6) then
      call mpp_error(FATAL, "size of second dimension of global_indices should be 6.")
    endif
    if (size(layout,1) .ne. 2) then
      call mpp_error(FATAL, "size of first dimension of layout should be 2.")
    endif
    if (size(layout,2) .ne. 6) then
      call mpp_error(FATAL, "size of second dimension of layout should be 6.")
    endif
    if (size(ni) .ne. 6 .or. size(nj) .ne. 6) then
      call mpp_error(FATAL, "size of ni and nj should be 6.")
    endif

    !Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1
    tile2(1) = 2
    istart1(1) = ni(1)
    iend1(1) = ni(1)
    jstart1(1) = 1
    jend1(1) = nj(1)
    istart2(1) = 1
    iend2(1) = 1
    jstart2(1) = 1
    jend2(1) = nj(2)

    !Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1
    tile2(2) = 3
    istart1(2) = 1
    iend1(2) = ni(1)
    jstart1(2) = nj(1)
    jend1(2) = nj(1)
    istart2(2) = 1
    iend2(2) = 1
    jstart2(2) = nj(3)
    jend2(2) = 1

    !Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1
    tile2(3) = 5
    istart1(3) = 1
    iend1(3) = 1
    jstart1(3) = 1
    jend1(3) = nj(1)
    istart2(3) = ni(5)
    iend2(3) = 1
    jstart2(3) = nj(5)
    jend2(3) = nj(5)

    !Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1
    tile2(4) = 6
    istart1(4) = 1
    iend1(4) = ni(1)
    jstart1(4) = 1
    jend1(4) = 1
    istart2(4) = 1
    iend2(4) = ni(6)
    jstart2(4) = nj(6)
    jend2(4) = nj(6)

    !Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2
    tile2(5) = 3
    istart1(5) = 1
    iend1(5) = ni(2)
    jstart1(5) = nj(2)
    jend1(5) = nj(2)
    istart2(5) = 1
    iend2(5) = ni(3)
    jstart2(5) = 1
    jend2(5) = 1

    !Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2
    tile2(6) = 4
    istart1(6) = ni(2)
    iend1(6) = ni(2)
    jstart1(6) = 1
    jend1(6) = nj(2)
    istart2(6) = ni(4)
    iend2(6) = 1
    jstart2(6) = 1
    jend2(6) = 1

    !Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2
    tile2(7) = 6
    istart1(7) = 1
    iend1(7) = ni(2)
    jstart1(7) = 1
    jend1(7) = 1
    istart2(7) = ni(6)
    iend2(7) = ni(6)
    jstart2(7) = nj(6)
    jend2(7) = 1

    !Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3
    tile2(8) = 4
    istart1(8) = ni(3)
    iend1(8) = ni(3)
    jstart1(8) = 1
    jend1(8) = nj(3)
    istart2(8) = 1
    iend2(8) = 1
    jstart2(8) = 1
    jend2(8) = nj(4)

    !Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3
    tile2(9) = 5
    istart1(9) = 1
    iend1(9) = ni(3)
    jstart1(9) = nj(3)
    jend1(9) = nj(3)
    istart2(9) = 1
    iend2(9) = 1
    jstart2(9) = nj(5)
    jend2(9) = 1

    !Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4
    tile2(10) = 5
    istart1(10) = 1
    iend1(10) = ni(4)
    jstart1(10) = nj(4)
    jend1(10) = nj(4)
    istart2(10) = 1
    iend2(10) = ni(5)
    jstart2(10) = 1
    jend2(10) = 1

    !Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4
    tile2(11) = 6
    istart1(11) = ni(4)
    iend1(11) = ni(4)
    jstart1(11) = 1
    jend1(11) = nj(4)
    istart2(11) = ni(6)
    iend2(11) = 1
    jstart2(11) = 1
    jend2(11) = 1

    !Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5
    tile2(12) = 6
    istart1(12) = ni(5)
    iend1(12) = ni(5)
    jstart1(12) = 1
    jend1(12) = nj(5)
    istart2(12) = 1
    iend2(12) = 1
    jstart2(12) = 1
    jend2(12) = nj(6)
    msize(1) = maxval(ni(:)/layout(1,:)) + whalo + ehalo + 1
    msize(2) = maxval(nj(:)/layout(2,:)) + shalo + nhalo + 1
    call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, &
                          tile2, istart1, iend1, jstart1, jend1, istart2, iend2, &
                          jstart2, jend2, pe_start, pe_end, symmetry = .true., &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          name=trim("Cubed-sphere"), memory_size=msize)
    call mpp_define_io_domain(domain, io_layout)
  end subroutine define_cubic_mosaic
end module offloading_io_mod
