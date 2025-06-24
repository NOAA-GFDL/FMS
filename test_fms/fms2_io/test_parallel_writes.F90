program test
  use fms_mod, only: fms_init, fms_end
  use mpp_domains_mod
  use mpp_mod
  use fms2_io_mod

  implicit none

  !< Namelist variables (configuration)
  integer, dimension(2) :: layout = (/1,6/)           !< Layout of the domain
  integer, dimension(2) :: io_layout = (/1,1/)        !< Io layout (currently overwritten to be the same as the
                                                      !! layout for reading, and 1,1 for writting)
  integer               :: nx = 96                    !< Size of the "x" dimension
  integer               :: ny = 96                    !< Size of the "y" dimension
  integer               :: nz = 65                    !< Size of the "z" dimension
  integer               :: ntimes = 2                 !< Number of time levels
  integer               :: test_case = 1              !< 1 use fms2io domain writes
                                                      !! 2 use netcdf collective writes

  character(len=10)     :: nc_format = "netcdf4"

  !< Domain stuff
  integer                               :: is, ie, js, je !< Starting and ending indices
  type(domain2d)                        :: Domain_write !< Domain of the data for when writing the file
                                                        !! (it uses an io_layout of 1,1)

  !< FMS2 io stuff
  character(len=6), dimension(4)        :: names        !< Dimension names for the dummy variables
  type(FmsNetcdfDomainFile_t)           :: fileobj, fileobj2             !< fms2io fileobj for domain decomposed

  real, allocatable, dimension(:,:)   :: sst_in         !< Data to be written
  integer                               :: i, j, k

  call fms_init()

  ! Create a domain
  call mpp_domains_set_stack_size(17280000)
  call mpp_define_domains( (/1,nx,1,ny/), layout, Domain_write, xhalo=0, yhalo=0)
  call mpp_define_io_domain(Domain_write, io_layout) !< io_layout is only relevant for fms2_io
  call mpp_get_data_domain(Domain_write, is, ie, js, je)

  ! Dummy data
  allocate(sst_in(is:ie, js:je))
  do i = is, ie
    do j = js, je
        sst_in(i,j) = i*10000. + j
    enddo
  enddo

  names(1) = "lon"
  names(2) = "lat"
  names(3) = "time"

  if (open_file(fileobj, "test_domain_io.nc", "overwrite", Domain_write, nc_format=nc_format)) then
    call register_axis(fileobj, names(1), "x")
    call register_axis(fileobj, names(2), "y")
    call register_axis(fileobj, names(3), unlimited)

    call register_field(fileobj, "sst_3d", "double", names(1:3))
    do i = 1, ntimes
      call write_data(fileobj, "sst_3d", sst_in, unlim_dim_level = i)
    enddo

    call close_file(fileobj)
  else
    call mpp_error(FATAL, "Unable to open the file for writing")
  endif

  call mpp_sync()

  if (open_file(fileobj2, "test_domain_io_efficient.nc", "overwrite", Domain_write, use_netcdf_mpi=.true.)) then
    call register_axis(fileobj2, names(1), "x")
    call register_axis(fileobj2, names(2), "y")
    call register_axis(fileobj2, names(3), unlimited)

    call register_field(fileobj2, "sst_3d", "double", names(1:3))

    do i = 1, ntimes
      call write_data(fileobj2, "sst_3d", sst_in, unlim_dim_level = i)
    enddo

    call close_file(fileobj2)
  else
    call mpp_error(FATAL, "Unable to open the file for writing")
  endif

  call fms_end()
end program