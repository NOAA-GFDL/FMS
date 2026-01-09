program test_io_offloading
  use fms_mod, only: fms_init, fms_end, string, check_nml_error
  use platform_mod
  use mpp_mod
  use mpp_domains_mod
  use offloading_io_mod
  use fms2_io_mod

  implicit none

  integer, parameter          :: lat_lon = 1     !< Using a lat lon domain
  integer, parameter          :: cube_sphere = 2 !< Using a cube sphere domain
  integer, allocatable        :: model_pes(:)    !< The model PEs
  integer, allocatable        :: offload_pes(:)  !< The PEs for offloading
  integer, allocatable        :: og_pes(:)       !< all of the PEs
  logical                     :: is_root_pe      !< The root PE from all of the PEs
  logical                     :: is_model_pe     !< .True. if current PE is a member of the model_pes
  logical                     :: is_offload_pe   !< .True. if current PE is a member of the offload_pes
  type(domain2D)              :: model_domain    !< Domain for the model PEs
  integer                     :: i               !< For do loops
  integer                     :: io              !< Error code when reading namelist
  integer                     :: ierr            !< Error code for namelist
  character(len=30)           :: filename        !< Filename for the test
  type(FmsNetcdfDomainFile_t) :: fileobj         !< Fileobj for the test
  real(kind=r4_kind), allocatable :: var_r4(:,:) !< Variable data for the "model"

  !< Namelist variables
  integer :: nx = 96                  !< Number of points in the x direction (per tile)
  integer :: ny = 96                  !< Number of points in the y direction (per tile)
  integer :: nxhalo = 2               !< Number of halo points in the x direction
  integer :: nyhalo = 2               !< Number of halo points in the x direction
  integer :: noffload_pes = 1         !< Number of PEs to use for offloading (need 1 per tile)
  integer :: nmodel_pes = 6           !< Number of PEs to use for the model
  integer :: io_layout(2)             !< Io layout to use for the model domain
  integer :: domain_type = lat_lon    !< The type of domain to use [lat_lon or cube_sphere]

  namelist /test_io_offloading_nml/ nx, ny, nxhalo, nyhalo, noffload_pes, nmodel_pes, domain_type

  call fms_init
  call offloading_io_init

  read(input_nml_file, nml=test_io_offloading_nml, iostat=io)
  ierr = check_nml_error(io, 'test_io_offloading_nml')

  if (mpp_npes() .ne. nmodel_pes + noffload_pes) &
    call mpp_error(FATAL, "The total number of PEs "//string(mpp_npes())//" is not equal to model_pes + noffload_pes")

  is_root_pe = mpp_pe() .eq. mpp_root_pe()
  allocate(og_pes(mpp_npes()))
  call mpp_get_current_pelist(og_pes)

  allocate(model_pes(nmodel_pes))
  model_pes(1) = 0
  do i = 2, nmodel_pes
    model_pes(i) = model_pes(i-1) + 1
  enddo
  if (is_root_pe) print *, "Model PEs:", model_pes
  call mpp_declare_pelist(model_pes, "model_pes")

  allocate(offload_pes(noffload_pes))
  offload_pes(1) = model_pes(nmodel_pes) + 1
  do i = 2, noffload_pes
    offload_pes(i) = offload_pes(i-1) + 1
  enddo
  if (is_root_pe) print *, "Offload PEs:", offload_pes
  call mpp_declare_pelist(offload_pes, "offload_pes")

  is_model_pe = .false.
  is_offload_pe = .false.
  if (any(model_pes .eq. mpp_pe())) is_model_pe = .true.
  if (any(offload_pes .eq. mpp_pe())) is_offload_pe = .true.

  if (is_model_pe) then
    call mpp_set_current_pelist( model_pes )
    ! Only the model pes are creating the domain and allocating the data
    select case (domain_type)
    case (lat_lon)
      model_domain = create_lat_lon_domain(nx, ny, halox=nxhalo, haloy=nyhalo)
    case (cube_sphere)
      model_domain = create_cubic_domain(nx, ny, 6, io_layout, nhalos=nxhalo)
    end select

    filename = "atmos.daily.nc"
    var_r4 = create_dummy_data(model_domain)
  endif

  call mpp_set_current_pelist(og_pes) !All of the PEs need to call the offloading stuff
  call open_file_offload(fileobj, filename, &
    model_domain, &
    model_pes, offload_pes)

  call global_metadata_offload(fileobj, "Number of times Fortran made you cry", 20)
  call global_metadata_offload(fileobj, "Number of lines of code", 19.54326541)

  call register_axis_offload(fileobj, "lon", "x")
  call register_axis_offload(fileobj, "lat", "x")

  call register_field_offload(fileobj, "mullions", "double", (/"lon", "lat"/))
  call write_data_offload(fileobj, "mullions", var_r4)
  call close_file_offload(fileobj)
  call fms_end

  contains

  function create_dummy_data(domain) &
    result(dummy_data)

    type(domain2D), intent(in) :: domain
    real(kind=r4_kind), allocatable :: dummy_data(:,:)

    integer :: is !< Starting x index
    integer :: ie !< Ending x index
    integer :: js !< Starting y index
    integer :: je !< Ending y index

    integer :: j, k

    !Allocate the data to the size of the data domain but only fill the compute domain with data
    call mpp_get_data_domain(domain, is, ie, js, je)
    allocate(dummy_data(is:ie, js:je))
    dummy_data = -999_r4_kind

    call mpp_get_compute_domain(domain, is, ie, js, je)
    do j = is, ie
      do k = js, je
        dummy_data(j, k) = real(j, kind=r4_kind)* 100_r4_kind + &
          real(k, kind=r4_kind)
      enddo
    enddo
  end function
end program test_io_offloading