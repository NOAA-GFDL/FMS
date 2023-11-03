program test
  use mpp_mod
  use fms_mod, only: fms_init, fms_end
  use mpp_domains_mod
  use offloading_io_mod
  use platform_mod

  implicit none

  integer :: noffload_pes = 1         !< Number of PEs to use for offloading (need 1 per tile)
  integer :: nmodel_pes = 6           !< Number of PEs to use for the model
  integer :: nx = 96                  !< Number of points in the x direction (per tile)
  integer :: ny = 96                  !< Number of points in the y direction (per tile)

  integer, allocatable        :: model_pes(:)    !< The model PEs
  integer, allocatable        :: offload_pes(:)  !< The PEs for offloading
  type(domain2D)              :: model_domain    !< Domain for the model PEs
  type(domain2D)              :: offload_domain    !< Domain for the offload PEs
  integer, allocatable        :: og_pes(:)       !< all of the PEs
  logical                     :: is_root_pe      !< The root PE from all of the PEs
  logical                     :: is_model_pe     !< .True. if current PE is a member of the model_pes
  logical                     :: is_offload_pe   !< .True. if current PE is a member of the offload_pes
  real(kind=r4_kind), allocatable :: var_r4_model(:,:) !< Variable data for the "model"
  real(kind=r4_kind), allocatable :: var_r4_offload(:,:) !< Variable data for the "model"

  integer :: i
  call fms_init()

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
    model_domain = create_lat_lon_domain(nx, ny)
    var_r4_model = create_dummy_data(model_domain, is_model_pe)
  endif

  if (is_offload_pe) then
    call mpp_set_current_pelist(offload_pes)
    offload_domain = create_lat_lon_domain(nx, ny)
    var_r4_offload = create_dummy_data(offload_domain, is_model_pe)
  endif

  call mpp_set_current_pelist(og_pes)
  call mpp_redistribute( model_domain, var_r4_model,  offload_domain, var_r4_offload)

  call fms_end()

  contains

  function create_dummy_data(domain, use_data_domain) &
    result(dummy_data)

    type(domain2D), intent(in) :: domain
    logical, intent(in) :: use_data_domain
    real(kind=r4_kind), allocatable :: dummy_data(:,:)

    integer :: is !< Starting x index
    integer :: ie !< Ending x index
    integer :: js !< Starting y index
    integer :: je !< Ending y index

    integer :: j, k

    !Allocate the data to the size of the data domain but only fill the compute domain with data

    if (use_data_domain) then
      call mpp_get_data_domain(domain, is, ie, js, je)
    else
      call mpp_get_compute_domain(domain, is, ie, js, je)
    endif
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
end program test