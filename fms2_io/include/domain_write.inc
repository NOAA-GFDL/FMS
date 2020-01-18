!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_0d(fileobj, variable_name, vdata, unlim_dim_level, corner)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), intent(in) :: vdata !< Data that will
                                !! be written out
                                !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true.)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner)
    return
  else
    call error("this branch should never be reached.")
  endif
end subroutine domain_write_0d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_1d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in) :: vdata !< Data that will
                                              !! be written out
                                              !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(1), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(1), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true.)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  else
    call error("this branch should never be reached.")
  endif
end subroutine domain_write_1d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_2d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in) :: vdata !< Data that will
                                                !! be written out
                                                !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(2), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(2), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: xdim_index
  integer :: ydim_index
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: xmax
  integer :: ymax
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  logical :: buffer_includes_halos
  
  
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer :: i
  integer(kind=i4_kind), dimension(:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:), allocatable :: buf_r8_kind

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root gathers the data and writes it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_global_domain(io_domain, xend=xmax, position=xpos)
    call mpp_get_global_domain(io_domain, yend=ymax, position=ypos)
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, position=xpos)
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, ysize=pe_jcsize, position=ypos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i) - pe_isc(1) + 1
      c(ydim_index) = pe_jsc(i) - pe_jsc(1) + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=i4_kind))
          call allocate_array(buf_i4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_i4_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i4_kind, size(buf_i4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_i4_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_i4_kind)
        type is (integer(kind=i8_kind))
          call allocate_array(buf_i8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_i8_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i8_kind, size(buf_i8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_i8_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_i8_kind)
        type is (real(kind=r4_kind))
          call allocate_array(buf_r4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_r4_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r4_kind, size(buf_r4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_r4_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_r4_kind)
        type is (real(kind=r8_kind))
          call allocate_array(buf_r8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_r8_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r8_kind, size(buf_r8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_r8_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_r8_kind)
        class default
          call error("unsupported type.")
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jcsize)
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size

    select type(vdata)
      type is (integer(kind=i4_kind))
        call allocate_array(buf_i4_kind, e)
        call get_array_section(buf_i4_kind, vdata, c, e)
        call mpp_send(buf_i4_kind, size(buf_i4_kind), fileobj%io_root)
      type is (integer(kind=i8_kind))
        call allocate_array(buf_i8_kind, e)
        call get_array_section(buf_i8_kind, vdata, c, e)
        call mpp_send(buf_i8_kind, size(buf_i8_kind), fileobj%io_root)
      type is (real(kind=r4_kind))
        call allocate_array(buf_r4_kind, e)
        call get_array_section(buf_r4_kind, vdata, c, e)
        call mpp_send(buf_r4_kind, size(buf_r4_kind), fileobj%io_root)
      type is (real(kind=r8_kind))
        call allocate_array(buf_r8_kind, e)
        call get_array_section(buf_r8_kind, vdata, c, e)
        call mpp_send(buf_r8_kind, size(buf_r8_kind), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
end subroutine domain_write_2d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_3d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in) :: vdata !< Data that will
                                                  !! be written out
                                                  !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(3), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(3), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: xdim_index
  integer :: ydim_index
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: xmax
  integer :: ymax
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  logical :: buffer_includes_halos
  
  
  integer, dimension(3) :: c
  integer, dimension(3) :: e
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer :: i
  integer(kind=i4_kind), dimension(:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:), allocatable :: buf_r8_kind

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root gathers the data and writes it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_global_domain(io_domain, xend=xmax, position=xpos)
    call mpp_get_global_domain(io_domain, yend=ymax, position=ypos)
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, position=xpos)
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, ysize=pe_jcsize, position=ypos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i) - pe_isc(1) + 1
      c(ydim_index) = pe_jsc(i) - pe_jsc(1) + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=i4_kind))
          call allocate_array(buf_i4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_i4_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i4_kind, size(buf_i4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_i4_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_i4_kind)
        type is (integer(kind=i8_kind))
          call allocate_array(buf_i8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_i8_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i8_kind, size(buf_i8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_i8_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_i8_kind)
        type is (real(kind=r4_kind))
          call allocate_array(buf_r4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_r4_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r4_kind, size(buf_r4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_r4_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_r4_kind)
        type is (real(kind=r8_kind))
          call allocate_array(buf_r8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_r8_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r8_kind, size(buf_r8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_r8_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_r8_kind)
        class default
          call error("unsupported type.")
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jcsize)
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size

    select type(vdata)
      type is (integer(kind=i4_kind))
        call allocate_array(buf_i4_kind, e)
        call get_array_section(buf_i4_kind, vdata, c, e)
        call mpp_send(buf_i4_kind, size(buf_i4_kind), fileobj%io_root)
      type is (integer(kind=i8_kind))
        call allocate_array(buf_i8_kind, e)
        call get_array_section(buf_i8_kind, vdata, c, e)
        call mpp_send(buf_i8_kind, size(buf_i8_kind), fileobj%io_root)
      type is (real(kind=r4_kind))
        call allocate_array(buf_r4_kind, e)
        call get_array_section(buf_r4_kind, vdata, c, e)
        call mpp_send(buf_r4_kind, size(buf_r4_kind), fileobj%io_root)
      type is (real(kind=r8_kind))
        call allocate_array(buf_r8_kind, e)
        call get_array_section(buf_r8_kind, vdata, c, e)
        call mpp_send(buf_r8_kind, size(buf_r8_kind), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
end subroutine domain_write_3d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_4d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in) :: vdata !< Data that will
                                                    !! be written out
                                                    !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(4), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(4), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: xdim_index
  integer :: ydim_index
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: xmax
  integer :: ymax
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  logical :: buffer_includes_halos
  
  
  integer, dimension(4) :: c
  integer, dimension(4) :: e
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer :: i
  integer(kind=i4_kind), dimension(:,:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:,:), allocatable :: buf_r8_kind

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root gathers the data and writes it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_global_domain(io_domain, xend=xmax, position=xpos)
    call mpp_get_global_domain(io_domain, yend=ymax, position=ypos)
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, position=xpos)
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, ysize=pe_jcsize, position=ypos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i) - pe_isc(1) + 1
      c(ydim_index) = pe_jsc(i) - pe_jsc(1) + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=i4_kind))
          call allocate_array(buf_i4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_i4_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i4_kind, size(buf_i4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_i4_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_i4_kind)
        type is (integer(kind=i8_kind))
          call allocate_array(buf_i8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_i8_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i8_kind, size(buf_i8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_i8_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_i8_kind)
        type is (real(kind=r4_kind))
          call allocate_array(buf_r4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_r4_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r4_kind, size(buf_r4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_r4_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_r4_kind)
        type is (real(kind=r8_kind))
          call allocate_array(buf_r8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_r8_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r8_kind, size(buf_r8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_r8_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_r8_kind)
        class default
          call error("unsupported type.")
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jcsize)
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size

    select type(vdata)
      type is (integer(kind=i4_kind))
        call allocate_array(buf_i4_kind, e)
        call get_array_section(buf_i4_kind, vdata, c, e)
        call mpp_send(buf_i4_kind, size(buf_i4_kind), fileobj%io_root)
      type is (integer(kind=i8_kind))
        call allocate_array(buf_i8_kind, e)
        call get_array_section(buf_i8_kind, vdata, c, e)
        call mpp_send(buf_i8_kind, size(buf_i8_kind), fileobj%io_root)
      type is (real(kind=r4_kind))
        call allocate_array(buf_r4_kind, e)
        call get_array_section(buf_r4_kind, vdata, c, e)
        call mpp_send(buf_r4_kind, size(buf_r4_kind), fileobj%io_root)
      type is (real(kind=r8_kind))
        call allocate_array(buf_r8_kind, e)
        call get_array_section(buf_r8_kind, vdata, c, e)
        call mpp_send(buf_r8_kind, size(buf_r8_kind), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
end subroutine domain_write_4d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_5d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in) :: vdata !< Data that will
                                                      !! be written out
                                                      !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(5), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(5), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: xdim_index
  integer :: ydim_index
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: xmax
  integer :: ymax
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  logical :: buffer_includes_halos
  
  
  integer, dimension(5) :: c
  integer, dimension(5) :: e
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer :: i
  integer(kind=i4_kind), dimension(:,:,:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:,:,:), allocatable :: buf_r8_kind

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root gathers the data and writes it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_global_domain(io_domain, xend=xmax, position=xpos)
    call mpp_get_global_domain(io_domain, yend=ymax, position=ypos)
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, position=xpos)
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, ysize=pe_jcsize, position=ypos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i) - pe_isc(1) + 1
      c(ydim_index) = pe_jsc(i) - pe_jsc(1) + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=i4_kind))
          call allocate_array(buf_i4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_i4_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i4_kind, size(buf_i4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_i4_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_i4_kind)
        type is (integer(kind=i8_kind))
          call allocate_array(buf_i8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_i8_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i8_kind, size(buf_i8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_i8_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_i8_kind)
        type is (real(kind=r4_kind))
          call allocate_array(buf_r4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_r4_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r4_kind, size(buf_r4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_r4_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_r4_kind)
        type is (real(kind=r8_kind))
          call allocate_array(buf_r8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call get_array_section(buf_r8_kind, vdata, c, e)
            if (buffer_includes_halos) then
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r8_kind, size(buf_r8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Write the out the data.
          call netcdf_write_data(fileobj, variable_name, buf_r8_kind, &
                                 unlim_dim_level=unlim_dim_level, corner=c, &
                                 edge_lengths=e)
          deallocate(buf_r8_kind)
        class default
          call error("unsupported type.")
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jcsize)
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size

    select type(vdata)
      type is (integer(kind=i4_kind))
        call allocate_array(buf_i4_kind, e)
        call get_array_section(buf_i4_kind, vdata, c, e)
        call mpp_send(buf_i4_kind, size(buf_i4_kind), fileobj%io_root)
      type is (integer(kind=i8_kind))
        call allocate_array(buf_i8_kind, e)
        call get_array_section(buf_i8_kind, vdata, c, e)
        call mpp_send(buf_i8_kind, size(buf_i8_kind), fileobj%io_root)
      type is (real(kind=r4_kind))
        call allocate_array(buf_r4_kind, e)
        call get_array_section(buf_r4_kind, vdata, c, e)
        call mpp_send(buf_r4_kind, size(buf_r4_kind), fileobj%io_root)
      type is (real(kind=r8_kind))
        call allocate_array(buf_r8_kind, e)
        call get_array_section(buf_r8_kind, vdata, c, e)
        call mpp_send(buf_r8_kind, size(buf_r8_kind), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
end subroutine domain_write_5d
