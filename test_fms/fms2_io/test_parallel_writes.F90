program test_parallel_writes
  use fms_mod,            only: fms_init, fms_end
  use platform_mod,       only: r8_kind
  use random_numbers_mod, only: randomNumberStream, initializeRandomNumberStream, getRandomNumbers
  use mpp_mod,            only: mpp_sync, mpp_broadcast, mpp_pe, mpp_root_pe, mpp_error, FATAL
  use mpp_domains_mod,    only: mpp_define_domains, mpp_define_io_domain, mpp_get_data_domain, domain2d, &
                                mpp_domains_set_stack_size
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, register_field, register_axis, unlimited, open_file, &
                                close_file, write_data, read_data
use fms_string_utils_mod, only: string
  implicit none

  integer, parameter :: nx = 96 !< Size of the "x" dimension
  integer, parameter :: ny = 96 !< Size of the "y" dimension
  integer, parameter :: nu = 5  !< Size of the "u" dimension
  integer, parameter :: nv = 4  !< Size of the "v" dimension
  integer, parameter :: nw = 3  !< Size of the "w" dimension
  integer, parameter :: nt = 2  !< Number of time levels

  integer, parameter :: layout(2) = [1,6] !< Layout of the domain
  type(domain2d) :: domain
  integer :: is, ie, js, je !< Starting and ending indices of the data domain

  character(*), parameter :: dim_names(*) = ["x", "y", "u", "v", "w", "t"]
  character(*), parameter :: filename = "test_parallel_writes.nc"

  interface init_data
    procedure :: init_data_1d
    procedure :: init_data_2d
    procedure :: init_data_3d
    procedure :: init_data_4d
    procedure :: init_data_5d
    procedure :: init_data_6d
  end interface init_data

  interface check_answers
    procedure :: check_answers_1d
    procedure :: check_answers_2d
    procedure :: check_answers_3d
    procedure :: check_answers_4d
    procedure :: check_answers_5d
    procedure :: check_answers_6d
  end interface check_answers

  real(r8_kind), allocatable :: data_0d(:,:)
  real(r8_kind), allocatable :: data_1d(:,:,:)
  real(r8_kind), allocatable :: data_2d(:,:,:,:)
  real(r8_kind), allocatable :: data_3d(:,:,:,:,:)
  real(r8_kind), allocatable :: data_4d(:,:,:,:,:,:)
  real(r8_kind), allocatable :: data_5d(:,:,:,:,:,:,:)

  type(randomNumberStream) :: random_stream

  call fms_init
  call define_domain

  random_stream = initializeRandomNumberStream(mpp_pe())

  allocate(data_0d(nt, 2))
  allocate(data_1d(nu, nt, 2))
  allocate(data_2d(is:ie, js:je, nt, 2))
  allocate(data_3d(is:ie, js:je, nu, nt, 2))
  allocate(data_4d(is:ie, js:je, nu, nv, nt, 2))
  allocate(data_5d(is:ie, js:je, nu, nv, nw, nt, 2))

  if (mpp_pe().eq.mpp_root_pe()) then
    call init_data(data_0d(:,1))
    call init_data(data_1d(:,:,1))
  endif

  call mpp_broadcast(data_0d(:,1), nt, mpp_root_pe())
  call mpp_broadcast(data_1d(:,:,1), nu*nt, mpp_root_pe())

  call init_data(data_2d(:,:,:,1))
  call init_data(data_3d(:,:,:,:,1))
  call init_data(data_4d(:,:,:,:,:,1))
  call init_data(data_5d(:,:,:,:,:,:,1))

  call write_netcdf
  call mpp_sync
  call read_netcdf

  call check_answers(data_0d)
  call check_answers(data_1d)
  call check_answers(data_2d)
  call check_answers(data_3d)
  call check_answers(data_4d)
  call check_answers(data_5d)

  call fms_end

contains

  subroutine define_domain
    call mpp_domains_set_stack_size(17280000)
    call mpp_define_domains( [1,nx,1,ny], layout, domain, xhalo=0, yhalo=0)
    call mpp_define_io_domain(domain, [1,1])
    call mpp_get_data_domain(domain, is, ie, js, je)
  end subroutine define_domain

  subroutine check_answers_1d(arr)
    real(r8_kind), intent(in) :: arr(:,:)

    if (any(arr(:, 1).ne.arr(:, 2))) then
      call mpp_error(FATAL, "The data written out do not match the data read back in")
    endif
  end subroutine check_answers_1d

  subroutine check_answers_2d(arr)
    real(r8_kind), intent(in) :: arr(:,:,:)

    if (any(arr(:, :, 1).ne.arr(:, :, 2))) then
      call mpp_error(FATAL, "The data written out do not match the data read back in")
    endif
  end subroutine check_answers_2d

  subroutine check_answers_3d(arr)
    real(r8_kind), intent(in) :: arr(:,:,:,:)

    if (any(arr(:, :, :, 1).ne.arr(:, :, :, 2))) then
      call mpp_error(FATAL, "The data written out do not match the data read back in")
    endif
  end subroutine check_answers_3d

  subroutine check_answers_4d(arr)
    real(r8_kind), intent(in) :: arr(:,:,:,:,:)

    if (any(arr(:, :, :, :, 1).ne.arr(:, :, :, :, 2))) then
      call mpp_error(FATAL, "The data written out do not match the data read back in")
    endif
  end subroutine check_answers_4d

  subroutine check_answers_5d(arr)
    real(r8_kind), intent(in) :: arr(:,:,:,:,:,:)

    if (any(arr(:, :, :, :, :, 1).ne.arr(:, :, :, :, :, 2))) then
      call mpp_error(FATAL, "The data written out do not match the data read back in")
    endif
  end subroutine check_answers_5d

  subroutine check_answers_6d(arr)
    real(r8_kind), intent(in) :: arr(:,:,:,:,:,:,:)

    if (any(arr(:, :, :, :, :, :, 1).ne.arr(:, :, :, :, :, :, 2))) then
      call mpp_error(FATAL, "The data written out do not match the data read back in")
    endif
  end subroutine check_answers_6d

  subroutine init_data_1d(arr)
    real(r8_kind), intent(out) :: arr(:)

    call getRandomNumbers(random_stream, arr)
  end subroutine init_data_1d

  subroutine init_data_2d(arr)
    real(r8_kind), intent(out) :: arr(:,:)

    call getRandomNumbers(random_stream, arr)
  end subroutine init_data_2d

  subroutine init_data_3d(arr)
    real(r8_kind), intent(out) :: arr(:,:,:)
    integer :: i, n

    n = size(arr, 3)
    do i = 1, n
      call init_data(arr(:, :, i))
    enddo
  end subroutine init_data_3d

  subroutine init_data_4d(arr)
    real(r8_kind), intent(out) :: arr(:,:,:,:)
    integer :: i, n

    n = size(arr, 4)
    do i = 1, n
      call init_data(arr(:, :, :, i))
    enddo
  end subroutine init_data_4d

  subroutine init_data_5d(arr)
    real(r8_kind), intent(out) :: arr(:,:,:,:,:)
    integer :: i, n

    n = size(arr, 5)
    do i = 1, n
      call init_data(arr(:, :, :, :, i))
    enddo
  end subroutine init_data_5d

  subroutine init_data_6d(arr)
    real(r8_kind), intent(out) :: arr(:,:,:,:,:,:)
    integer :: i, n

    n = size(arr, 6)
    do i = 1, n
      call init_data(arr(:, :, :, :, :, i))
    enddo
  end subroutine init_data_6d

  subroutine netcdf_file_register_fields(fileobj)
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj

    call register_axis(fileobj, dim_names(1), "x")
    call register_axis(fileobj, dim_names(2), "y")
    call register_axis(fileobj, dim_names(3), nu)
    call register_axis(fileobj, dim_names(4), nv)
    call register_axis(fileobj, dim_names(5), nw)
    call register_axis(fileobj, dim_names(6), unlimited)

    call register_field(fileobj, "random_0d", "double", &
                        pack(dim_names, [.false., .false., .false., .false., .false., .true.]))
    call register_field(fileobj, "random_1d", "double", &
                        pack(dim_names, [.false., .false., .true., .false., .false., .true.]))
    call register_field(fileobj, "random_2d", "double", &
                        pack(dim_names, [.true., .true., .false., .false., .false., .true.]))
    call register_field(fileobj, "random_3d", "double", &
                        pack(dim_names, [.true., .true., .true., .false., .false., .true.]))
    call register_field(fileobj, "random_4d", "double", &
                        pack(dim_names, [.true., .true., .true., .true., .false., .true.]))
    call register_field(fileobj, "random_5d", "double", &
                        pack(dim_names, [.true., .true., .true., .true., .true., .true.]))
  end subroutine netcdf_file_register_fields

  subroutine write_netcdf
    type(FmsNetcdfDomainFile_t) :: fileobj
    integer :: i

    if (.not.open_file(fileobj, filename, "overwrite", domain, use_netcdf_mpi=.true.)) then
      call mpp_error(FATAL, "Unable to open the NetCDF file for writing")
    endif

    call netcdf_file_register_fields(fileobj)

    do i = 1, nt
      call write_data(fileobj, "random_0d", data_0d(i, 1),                unlim_dim_level=i)
      call write_data(fileobj, "random_1d", data_1d(:, i, 1),             unlim_dim_level=i)
      call write_data(fileobj, "random_2d", data_2d(:, :, i, 1),          unlim_dim_level=i)
      call write_data(fileobj, "random_3d", data_3d(:, :, :, i, 1),       unlim_dim_level=i)
      call write_data(fileobj, "random_4d", data_4d(:, :, :, :, i, 1),    unlim_dim_level=i)
      call write_data(fileobj, "random_5d", data_5d(:, :, :, :, :, i, 1), unlim_dim_level=i)
    enddo

    call close_file(fileobj)
  end subroutine write_netcdf

  subroutine read_netcdf
    type(FmsNetcdfDomainFile_t) :: fileobj
    integer :: i

    if (.not.open_file(fileobj, filename, "read", domain)) then
      call mpp_error(FATAL, "Unable to open the NetCDF file for reading")
    endif

    call netcdf_file_register_fields(fileobj)

    do i = 1, nt
      call read_data(fileobj, "random_0d", data_0d(i, 2),                unlim_dim_level=i)
      call read_data(fileobj, "random_1d", data_1d(:, i, 2),             unlim_dim_level=i)
      call read_data(fileobj, "random_2d", data_2d(:, :, i, 2),          unlim_dim_level=i)
      call read_data(fileobj, "random_3d", data_3d(:, :, :, i, 2),       unlim_dim_level=i)
      call read_data(fileobj, "random_4d", data_4d(:, :, :, :, i, 2),    unlim_dim_level=i)
      call read_data(fileobj, "random_5d", data_5d(:, :, :, :, :, i, 2), unlim_dim_level=i)
    enddo

    call close_file(fileobj)
  end subroutine read_netcdf
end program test_parallel_writes
