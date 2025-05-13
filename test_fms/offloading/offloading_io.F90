module offloading_io_mod
  use mpp_mod
  use mpp_domains_mod
  use fms_string_utils_mod, only: string
  use fms2_io_mod
  use offloading_io_obj_mod
  use platform_mod

  implicit none

  integer, parameter :: domain_decomposed = 0
  integer, parameter :: non_domain_decomposed = 1
  integer, parameter :: Unstructured_grid = 2
  integer :: nfiles = 10
  integer :: current_files_init
  logical :: module_is_initialized

  type :: offloading_obj_out
    integer :: id
    character(len=:), allocatable :: filename
    class(FmsNetcdfFile_t), allocatable :: fileobj
    type(domain2D) :: domain_out
    integer :: type_of_domain
  end type

  type :: offloading_obj_in
  end type

  type(offloading_obj_out), allocatable, target :: offloading_objs(:)

  private

  !TODO create domain should be private
  public :: offloading_io_init, open_file_offload, create_lat_lon_domain, create_cubic_domain
  public :: domain_decomposed, non_domain_decomposed, Unstructured_grid
  public :: global_metadata_offload, close_file_offload, register_axis_offload, register_field_offload
  public :: write_data_offload

  contains

  subroutine offloading_io_init()
    if (module_is_initialized) return

    current_files_init = 0
    allocate(offloading_objs(nfiles))
    module_is_initialized = .true.
  end subroutine offloading_io_init

  subroutine open_file_offload(fileobj, filename, domain_in, pe_in, pe_out, is_pe_in, is_pe_out)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj
    character(len=*),       intent(in)    :: filename
    type(domain2D),         intent(inout) :: domain_in
    integer,                intent(in)    :: pe_in(:) !< model pes
    integer,                intent(in)    :: pe_out(:) !< offload pes
    logical,                intent(in)    :: is_pe_in
    logical,                intent(in)    :: is_pe_out

    integer, parameter :: str_len = 255
    character(len=str_len) :: filename_out(1)
    integer :: object_id
    integer :: type_of_domain
    integer :: global_domain_size(2)
    integer :: ntile
    integer, allocatable :: all_current_pes(:)
    integer, allocatable :: broadcasting_pes(:)
    type(domain2D) :: domain_out

    ! This should be called from the model PEs and the offload PEs
    if (.not. module_is_initialized) &
      call mpp_error(FATAL, "offloading_io_mod is not initialized")

    allocate(all_current_pes(mpp_npes()))
    call mpp_get_current_pelist(all_current_pes)

    select type(fileobj)
    type is (FmsNetcdfDomainFile_t)
      type_of_domain = non_domain_decomposed
    type is (FmsNetcdfFile_t)
      type_of_domain = domain_decomposed
    type is (FmsNetcdfUnstructuredDomainFile_t)
      type_of_domain = Unstructured_grid
    end select

    filename_out(1) = ""
    if (mpp_pe() .eq. pe_in(1)) then
      !< The root model pe gets the domain info (ntiles and size of global domain)
      ntile = mpp_get_ntile_count(domain_in)

      !< The number of tiles must be the same as the number of offloading pes
      if (size(pe_out) .ne. ntile ) &
        call mpp_error(FATAL, "The number of offloading PEs must be the same as the number of tiles of the domain")
      filename_out(1) = filename

      call mpp_get_global_domain(domain_in, xsize=global_domain_size(1), ysize=global_domain_size(2))
    endif

    ! A "Root" model PE broadcasts the filename, domain size, and number of tiles to the offload pes
    if (mpp_pe() .eq. pe_in(1) .or. is_pe_out) then
      allocate(broadcasting_pes(1 + size(pe_out)))
      broadcasting_pes(1) = pe_in(1)
      broadcasting_pes(2:size(broadcasting_pes)) = pe_out
      call mpp_set_current_pelist( broadcasting_pes )
      call mpp_broadcast(filename_out, str_len, pe_in(1))
      call mpp_broadcast(global_domain_size, size(global_domain_size), pe_in(1))
      call mpp_broadcast(ntile, pe_in(1))
    endif

    ! Broadcast the domain
    call mpp_set_current_pelist(all_current_pes)
    if (is_pe_out) call mpp_define_null_domain(domain_in)
    call mpp_broadcast_domain(domain_in)

    ! The offload pes inits the t offloading objecand return an object id
    if (is_pe_out) then
      call mpp_set_current_pelist(pe_out)
      print *, "PE", string(mpp_pe()), " knows that the filename is ", trim(filename_out(1)), &
        " and the domain size ", string(global_domain_size(1)), " by ", string(global_domain_size(2)), &
        " and the number of tiles is ", string(ntile)
      object_id = init_offloading_object(filename_out(1), global_domain_size(1), global_domain_size(2),&
         ntile)
    endif
    call mpp_set_current_pelist(all_current_pes)
    call mpp_broadcast(object_id, pe_out(1))

    ! Init the "offloading object" in the fileobj
    call fileobj%offloading_obj_in%init(object_id, pe_out, pe_in, domain_in)
  end subroutine open_file_offload

  subroutine global_metadata_offload(fileobj, attribute_name, attribute_value)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj
    character(len=*), intent(in) :: attribute_name
    class(*), intent(in) :: attribute_value

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
    real :: real_buf

    offloading_pes = fileobj%offloading_obj_in%offloading_pes
    model_pes = fileobj%offloading_obj_in%model_pes
    is_model_pe = fileobj%offloading_obj_in%is_model_pe

    id = fileobj%offloading_obj_in%id
    this => offloading_objs(id)

    if (is_model_pe) then
      att_name(1) = attribute_name
    endif

    allocate(all_current_pes(mpp_npes()))
    call mpp_get_current_pelist(all_current_pes)

    if (mpp_pe() .eq. model_pes(1) .or. .not. is_model_pe) then

      allocate(broadcasting_pes(1 + size(offloading_pes)))
      broadcasting_pes(1) = model_pes(1)
      broadcasting_pes(2:size(broadcasting_pes)) = offloading_pes
      call mpp_set_current_pelist( broadcasting_pes )
      call mpp_broadcast(att_name, 255, model_pes(1))

      select type (attribute_value)
        type is (real(kind=r4_kind))
          print *, trim(attribute_name), " is a r4"
        type is (real(kind=r8_kind))
          if (is_model_pe) real_buf = attribute_value
          call mpp_broadcast(real_buf, model_pes(1))
          if (.not. is_model_pe) call register_global_attribute(this%fileobj, att_name(1), real_buf)
        type is (integer(kind=i4_kind))
          if (is_model_pe) int_buf = attribute_value
          call mpp_broadcast(int_buf, model_pes(1))
          if (.not. is_model_pe) call register_global_attribute(this%fileobj, att_name(1), int_buf)
        type is (integer(kind=i8_kind))
          print *, trim(attribute_name), " is a i8"
        type is (character(len=*))
          print *, trim(attribute_name), " is a string"
      end select
    endif

    call mpp_set_current_pelist(all_current_pes)
  end subroutine

  !TODO Need an interface and more
  subroutine register_axis_offload(fileobj, axis_name, cart)
    class(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    character(len=*), intent(in) :: axis_name
    character(len=1), intent(in) :: cart

    integer :: id
    type(offloading_obj_out), pointer :: this

    !TODO better PEs management!
    integer, allocatable :: offloading_pes(:)
    integer, allocatable :: model_pes(:)
    integer, allocatable :: all_current_pes(:)
    integer, allocatable :: broadcasting_pes(:)
    logical :: is_model_pe
    character(len=255) :: var_info(2)

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
      select type(wut=>this%fileobj)
        type is(FmsNetcdfDomainFile_t)
          call register_axis(wut, trim(var_info(1)), trim(var_info(2)))
      end select
    endif

    call mpp_set_current_pelist(all_current_pes)
  end subroutine

  subroutine register_field_offload(fileobj, varname, vartype, dimensions)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: vartype
    character(len=*), intent(in) :: dimensions(:)

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
      select type(wut=>this%fileobj)
        type is(FmsNetcdfDomainFile_t)
          call register_field(wut, trim(var_info(1)), trim(var_info(2)), var_info(3:))
      end select
    endif

    call mpp_set_current_pelist(all_current_pes)
  end subroutine

  subroutine write_data_offload(fileobj, varname, vardata)
    class(FmsNetcdfDomainFile_t), intent(inout) :: fileobj
    character(len=*), intent(in) :: varname
    real(kind=r4_kind), intent(in) :: vardata(:,:)

    integer :: id
    type(offloading_obj_out), pointer :: this

    !TODO better PEs management!
    integer, allocatable :: offloading_pes(:)
    integer, allocatable :: model_pes(:)
    integer, allocatable :: all_current_pes(:)
    integer, allocatable :: broadcasting_pes(:)
    logical :: is_model_pe

    real(kind=r4_kind), allocatable :: var_r4_data(:,:)
    type(domain2D) :: domain_out
    type(domain2D) :: domain_in
    integer :: isc, iec, jsc, jec

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

    call mpp_redistribute( domain_in, vardata,  &
    domain_out, var_r4_data &
      )

      if (.not. is_model_pe) then
        select type(wut=>this%fileobj)
          type is(FmsNetcdfDomainFile_t)
            call write_data(wut, "mullions", var_r4_data) !TODO mullions need to be broadcasted
        end select
      endif
    call mpp_set_current_pelist(all_current_pes)
  end subroutine

  subroutine close_file_offload(fileobj)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj

    integer :: id
    type(offloading_obj_out), pointer :: this
    logical :: is_model_pe

    id = fileobj%offloading_obj_in%id
    this => offloading_objs(id)
    is_model_pe = fileobj%offloading_obj_in%is_model_pe

    if (.not. is_model_pe) call close_file(this%fileobj)
  end subroutine close_file_offload

  integer function init_offloading_object(filename, nx, ny, ntile)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: ntile

    type(offloading_obj_out), pointer :: this

    current_files_init = current_files_init + 1
    if (current_files_init .gt. nfiles) &
      call mpp_error(FATAL, "The number of files is too large")

    ! An array of offloading objects is stored in this module, the "id" is the index of the object in the array
    this => offloading_objs(current_files_init)
    this%id = current_files_init
    this%filename = trim(filename)
    this%type_of_domain = domain_decomposed

    select case (ntile)
    case (1)
      this%domain_out = create_lat_lon_domain(nx, ny)
    case (6)
      this%domain_out = create_cubic_domain(nx, ny, ntile, (/1,1/))
    end select

    allocate(FmsNetcdfDomainFile_t :: this%fileobj)
    select type (fileobj => this%fileobj)
    type is (FmsNetcdfDomainFile_t)
      if ( .not. open_file(fileobj, trim(this%filename), "overwrite", this%domain_out)) &
        call mpp_error(FATAL, "Error opening file")
    end select

    init_offloading_object = current_files_init
  end function

  function create_lat_lon_domain(nx_in, ny_in, halox, haloy ) &
    result(domain_out)
    integer, intent(in) :: nx_in
    integer, intent(in) :: ny_in
    integer, intent(in), optional :: halox
    integer, intent(in), optional :: haloy
    type(domain2d) :: domain_out
    integer :: layout(2)

    call mpp_define_layout( (/1,nx_in,1,ny_in/), mpp_npes(), layout )
    call mpp_define_domains( (/1,nx_in,1,ny_in/), layout, domain_out, xhalo=halox, yhalo=haloy)
    call mpp_define_io_domain(domain_out, (/1,1/))
  end function create_lat_lon_domain

  !Assumes all members of the domain are in the current pelist
  function create_cubic_domain(nx_in, ny_in, ntiles, io_layout, nhalos) &
    result(domain_out)
    integer, intent(in) :: nx_in
    integer, intent(in) :: ny_in
    integer, intent(in) :: ntiles
    integer, intent(in) :: io_layout(2)
    integer, intent(in), optional :: nhalos

    type(domain2d) :: domain_out

    integer :: layout(2)
    integer :: npes
    integer :: npes_per_tile
    integer, allocatable :: global_indices(:,:)
    integer, allocatable :: layout2D(:,:)
    integer, allocatable :: pe_start(:)
    integer, allocatable :: pe_end(:)
    integer :: n

    npes = mpp_npes()
    if( mod(npes, ntiles) .NE. 0 ) call mpp_error(FATAL, &
          "create_cubic_domain: npes is not divisible by ntiles")

    npes_per_tile = npes/ntiles
    call mpp_define_layout ((/1,nx_in,1,ny_in/), npes_per_tile, layout )

    allocate(global_indices(4, ntiles))
    allocate(layout2D(2, ntiles))
    allocate(pe_start(ntiles), pe_end(ntiles))

    do n = 1, ntiles
      global_indices(:,n) = (/1,nx_in,1,ny_in/)
      layout2D(:,n) = layout
      pe_start(n) = (n-1)*npes_per_tile
      pe_end(n) = n*npes_per_tile-1
    end do

    call define_cubic_mosaic(domain_out, (/nx_in,nx_in,nx_in,nx_in,nx_in,nx_in/), &
                            (/ny_in,ny_in,ny_in,ny_in,ny_in,ny_in/), &
                              global_indices, layout2D, pe_start, pe_end, io_layout, nhalos )
  end function create_cubic_domain

  !> @brief Initialize a cubed-sphere atomsphere domain.
subroutine define_cubic_mosaic(domain, ni, nj, global_indices, layout, pe_start, pe_end, &
                                    io_layout, nhalos)

  integer, dimension(:), intent(in) :: ni
  integer, dimension(:), intent(in) :: nj
  integer, dimension(:,:), intent(in) :: global_indices
  integer, dimension(:,:), intent(in) :: layout
  integer, dimension(:), intent(in) :: pe_start
  integer, dimension(:), intent(in) :: pe_end
  integer, dimension(2), intent(in) :: io_layout
  type(domain2d), intent(inout) :: domain !< A cubed-sphere domain.
  integer, optional, intent(in) :: nhalos

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