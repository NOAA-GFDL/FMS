subroutine atmosphere_restart_file(domain,nz,nt,debug_flag)
  type(domain2d),intent(in),target :: domain
  integer,intent(in) :: nz
  integer,intent(in) :: nt
  logical,intent(in),optional :: debug_flag

  type(FmsNetcdfDomainFile_t) :: fileobj
  type(domain2d),pointer :: io_domain
  integer :: nx
  integer :: ny
  real(kind=real64),dimension(:),allocatable :: lon_data
  real(kind=real64),dimension(:),allocatable :: lat_data
  real(kind=real64),dimension(:),allocatable :: lev_data
  real(kind=real64),dimension(:,:,:),allocatable :: temperature_data
  real(kind=real64),dimension(:),allocatable :: avg_t_data
  integer(kind=int32),dimension(:,:),allocatable :: surf_t_data
  integer,dimension(:),allocatable :: dim_sizes
  integer(kind=int64),dimension(3) :: out_chksum
  integer(kind=int64),dimension(3) :: in_chksum
  integer,dimension(3) :: chksum_unlim_level
  logical :: debug
  integer :: err
  integer :: i
  integer :: j
  integer :: k
  integer :: l

  !Initialize.
  call mpi_barrier(mpi_comm_world,err)
  call mpi_check(err)
  if (present(debug_flag)) then
    debug = debug_flag
  else
    debug = .false.
  endif
  call random_seed()

  if (debug) then
    !Tell which test is running.
    if (mpp_pe() .eq. 0) then
      write(error_unit,'(/a)') &
          "Running atmosphere (6-tile domain decomposed)" &
          //" restart file test ... "
    endif
  endif
  call mpi_barrier(mpi_comm_world,err)
  call mpi_check(err)

  !Get the sizes of the I/O compute domain and allocate buffer.
  io_domain => mpp_get_io_domain(domain)
  if (.not. associated(io_domain)) then
    call mpp_error(fatal, &
                   "I/O domain is not associated.")
  endif
  call mpp_get_compute_domain(io_domain, &
                              xsize=nx, &
                              ysize=ny)

  !Write a restart file.
  !Open a restart file and initialize the file object.
  call open_check(open_file(fileobj, &
                            "atmosphere.nc", &
                            "overwrite", &
                            domain, &
                            nc_format="64bit", &
                            is_restart=.true.))

  !Add some global attributes to the restart file.
  call register_global_attribute(fileobj, &
                                 "globalatt1", &
                                 real(7.,kind=real64))
  call register_global_attribute(fileobj, &
                                 "globalatt2", &
                                 real(4.,kind=real32))
  call register_global_attribute(fileobj, &
                                 "globalatt3", &
                                 int(3,kind=int32))
  call register_global_attribute(fileobj, &
                                 "globalatt4", &
                                 int(2,kind=int64))
  call register_global_attribute(fileobj, &
                                 "globalatt5", &
                                 "some text")

  !Add dimensions and corresponding variables to the file.
  !Longitude (domain "x" dimension).
  call register_axis(fileobj, &
                     "lon", &
                     "x")
  call register_field(fileobj, &
                      "lon", &
                      "double", &
                      (/"lon"/))
  allocate(lon_data(nx))
  do i = 1,nx
    lon_data(i) = 33.33 + 1.2*i
  enddo
  call write_data(fileobj, &
                  "lon", &
                  lon_data)
  call register_variable_attribute(fileobj, &
                                   "lon", &
                                   "units", &
                                   "degrees")

  !Longitude2 (domain "x" dimension).
  call register_axis(fileobj, &
                     "lon2", &
                     "x")
  call register_field(fileobj, &
                      "lon2", &
                      "double", &
                      (/"lon2"/))
  do i = 1,nx
    lon_data(i) = 2.*3.14*lon_data(i)/(360.)
  enddo
  call write_data(fileobj, &
                  "lon2", &
                  lon_data)
  call register_variable_attribute(fileobj, &
                                   "lon2", &
                                   "units", &
                                   "radians")

  !Latitude (domain "y" dimension).
  call register_axis(fileobj, &
                     "lat", &
                     "y")
  call register_field(fileobj, &
                      "lat", &
                      "double", &
                      (/"lat"/))
  allocate(lat_data(ny))
  do i = 1,ny
    lat_data(i) = 65. + 0.05*i
  enddo
  call write_data(fileobj, &
                  "lat", &
                  lat_data)
  call register_variable_attribute(fileobj, &
                                   "lat", &
                                   "units", &
                                   "degrees")

  !Latitude2 (domain "y" dimension).
  call register_axis(fileobj, &
                     "lat2", &
                     "y")
  call register_field(fileobj, &
                      "lat2", &
                      "double", &
                      (/"lat2"/))
  do i = 1,ny
    lat_data(i) = 2.*3.14*(lat_data(i))/360.
  enddo
  call write_data(fileobj, &
                  "lat2", &
                  lat_data)
  call register_variable_attribute(fileobj, &
                                   "lat", &
                                   "units2", &
                                   "radians")

  !Height.
  call register_axis(fileobj, &
                     "lev", &
                     nz)
  call register_field(fileobj, &
                      "lev", &
                      "double", &
                      (/"lev"/))
  allocate(lev_data(nz))
  do i = 1,nz
    lev_data(i) = 0.1 + 20.2*i
  enddo
  call write_data(fileobj, &
                  "lev", &
                  lev_data)
  call register_variable_attribute(fileobj, &
                                   "lev", &
                                   "units", &
                                   "mb")

  !Height2.
  call register_axis(fileobj, &
                     "lay", &
                     nz-1)
  call register_field(fileobj, &
                      "lay", &
                      "double", &
                      (/"lay"/))
  do i = 1,nz-1
    lev_data(i) = lev_data(i) - 0.5
  enddo
  call write_data(fileobj, &
                  "lay", &
                  lev_data(1:nz-1))
  call register_variable_attribute(fileobj, &
                                   "lay", &
                                   "units", &
                                   "mb")

  !Time.
  call register_axis(fileobj, &
                     "time", &
                     unlimited)
  call register_field(fileobj, &
                      "time", &
                      "float", &
                      (/"time"/))
  call register_variable_attribute(fileobj, &
                                   "time", &
                                   "units", &
                                   "years")

  !Add a restart variable.
  allocate(temperature_data(nx,ny,nz))
  call register_restart_field(fileobj, &
                              "temperature", &
                              temperature_data, &
                              (/"lon ","lat ","lev ","time"/))
  call register_variable_attribute(fileobj, &
                                   "temperature", &
                                   "units", &
                                   "K")

  !Add a non-domain decomposed variable.
  allocate(avg_t_data(nz-1))
  call register_field(fileobj, &
                      "avg layer temperature", &
                      "double", &
                      (/"lay ","time"/))
  call register_variable_attribute(fileobj, &
                                   "avg layer temperature", &
                                   "units", &
                                   "K")

  !Add a domain decomposed variable.
  allocate(surf_t_data(nx,ny))
  call register_field(fileobj, &
                      "surface temperature", &
                      "int", &
                      (/"lon2","lat2","time"/))
  call register_variable_attribute(fileobj, &
                                   "surface temperature", &
                                   "units", &
                                   "K")

  !Initialize some random temperature data.
  do j = 1,nz
    do k = 1,ny
      do l = 1,nx
        call random_number(temperature_data(l,k,j))
        temperature_data(l,k,j) = temperature_data(l,k,j) + 300.05
      enddo
    enddo
  enddo
  do i = 1,nz-1
    avg_t_data(i) = sum(temperature_data(:,:,i))/(nx*ny)
  enddo
  do i = 1,ny
    do j = 1,nx
      surf_t_data(j,i) = temperature_data(j,i,1)
    enddo
  enddo

  !Perform a "simulation".
  do i = 1,nt
    do j = 1,nz
      do k = 1,ny
        do l = 1,nx
          temperature_data(l,k,j) = temperature_data(l,k,j) + 0.5*mpp_pe()
        enddo
      enddo
    enddo
    do j = 1,nz-1
      avg_t_data(j) = sum(temperature_data(:,:,j))/(nx*ny)
    enddo
    call write_data(fileobj, &
                    "avg layer temperature", &
                    avg_t_data, &
                    unlim_dim_level=i)
    if (i .eq. 1) then
      out_chksum(2) = mpp_chksum(avg_t_data,pelist=(/mpp_pe()/))
      chksum_unlim_level(2) = i
    endif

    do j = 1,ny
      do k = 1,nx
        surf_t_data(k,j) = int(temperature_data(k,j,1))
      enddo
    enddo
    call write_data(fileobj, &
                    "surface temperature", &
                    surf_t_data, &
                    unlim_dim_level=i)
    if (i .eq. 2) then
      out_chksum(3) = mpp_chksum(surf_t_data,pelist=(/mpp_pe()/))
      chksum_unlim_level(3) = i
    endif

    !Write out the times of each "time step".
    call write_data(fileobj, &
                    "time", &
                    1950. + real(i,kind=real32), &
                    unlim_dim_level=i)
  enddo

  !Write out the restart file.
  call write_restart(fileobj, &
                     unlim_dim_level=nt)
  out_chksum(1) = mpp_chksum(temperature_data,pelist=(/mpp_pe()/))
  chksum_unlim_level(1) = nt

  !Close the file and destroy the file object.
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world,err)
  call mpi_check(err)

  !Prepare to read the data back in.
  deallocate(lon_data)
  deallocate(lat_data)
  deallocate(lev_data)
  deallocate(temperature_data)
  deallocate(avg_t_data)
  deallocate(surf_t_data)

  !Check if a non-existent file exists (just to see if this feature
  !works.
  if (open_file(fileobj, &
                "atmosphere.foobar.nc", &
                "read", &
                domain, &
                nc_format="64bit", &
                is_restart=.true.)) then
    call mpp_error(FATAL, &
                   "Found non-existent file.")
  endif

  !Re-open the restart file and re-initialize the file object.
  call open_check(open_file(fileobj, &
                            "atmosphere.nc", &
                            "read", &
                            domain, &
                            nc_format="64bit", &
                            is_restart=.true.))

  !Associate the lon and lat dimensions with the
  !"x" and "y" dimensions of the domain.
  call register_axis(fileobj, &
                     "lon", &
                     "x")
  call register_axis(fileobj, &
                     "lat", &
                     "y")
  call register_axis(fileobj, &
                     "lon2", &
                     "x")
  call register_axis(fileobj, &
                     "lat2", &
                     "y")

  !Allocate buffers associated with the temperature variables.
  call get_variable_size(fileobj, &
                         "temperature", &
                         dim_sizes, &
                         broadcast=.true.)
  allocate(temperature_data(nx,ny,dim_sizes(3)))
  temperature_data(:,:,:) = 0.
  call get_variable_size(fileobj, &
                         "avg layer temperature", &
                         dim_sizes, &
                         broadcast=.true.)
  allocate(avg_t_data(dim_sizes(1)))
  avg_t_data(:) = 0.
  deallocate(dim_sizes)
  allocate(surf_t_data(nx,ny))
  surf_t_data(:,:) = 0

  !Register the restart variable (so a pointer can point at the newly
  !allocated buffer.
  call register_restart_field(fileobj, &
                              "temperature", &
                              temperature_data)

  !Read in the restart variable and checksum.
  call read_restart(fileobj, &
                    unlim_dim_level=chksum_unlim_level(1))
  in_chksum(1) = mpp_chksum(temperature_data,pelist=(/mpp_pe()/))
  call chksum_match(out_chksum(1), &
                    in_chksum(1), &
                    "temperature", &
                    debug)

  !Read in the non-restart variables and checksum.
  call read_data(fileobj, &
                 "avg layer temperature", &
                 avg_t_data, &
                 unlim_dim_level=chksum_unlim_level(2))
! in_chksum(2) = mpp_chksum(avg_t_data,pelist=(/mpp_pe()/))
! call chksum_match(out_chksum(2), &
!                   in_chksum(2), &
!                   "avg layer temperature", &
!                   debug)
  call read_data(fileobj, &
                 "surface temperature", &
                 surf_t_data, &
                 unlim_dim_level=chksum_unlim_level(3))
  in_chksum(3) = mpp_chksum(surf_t_data,pelist=(/mpp_pe()/))
  call chksum_match(out_chksum(3), &
                    in_chksum(3), &
                    "surface temperature", &
                    debug)

  !Close the file.
  call close_file(fileobj)
  call mpi_barrier(mpi_comm_world,err)
  call mpi_check(err)

  if (mpp_pe() .eq. 0) then
    write(output_unit,'(a/)') &
      "Atmosphere (6-tile domain decomposed) file test ... " &
      //trim(green)//"passed."//trim(color_end)
  endif
  call mpi_barrier(mpi_comm_world,err)
  call mpi_check(err)
end subroutine atmosphere_restart_file
