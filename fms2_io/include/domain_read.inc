!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_0d(fileobj, variable_name, vdata, unlim_dim_level, corner)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), intent(inout) :: vdata !< Data that will
                                   !! be written out
                                   !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true.)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          broadcast=.true.)
    return
  else
    call error("this branch should never be reached.")
  endif
end subroutine domain_read_0d


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_1d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(inout) :: vdata !< Data that will
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
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  else
    call error("this branch should never be reached.")
  endif
end subroutine domain_read_1d


!> @brief Utility routine that retrieves domain indices.
!! @internal
subroutine domain_offsets(data_xsize, data_ysize, io_domain, dpos, &
                          isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos)

  integer, intent(in) :: data_xsize !< Size of buffer's domain "x" dimension.
  integer, intent(in) :: data_ysize !< Size of buffer's domain "y" dimension.
  type(domain2d), intent(in) :: io_domain !< I/O domain variable is decomposed over.
  integer, intent(in) :: dpos !< Variable's domain position.
  integer, intent(out) :: isd !< Starting index for x dimension of data domain.
  integer, intent(out) :: isc !< Starting index for x dimension of compute domain.
  integer, intent(out) :: xc_size !< Size of x dimension of compute domain.
  integer, intent(out) :: jsd !< Starting index for y dimension of data domain.
  integer, intent(out) :: jsc !< Starting index for y dimension of compute domain.
  integer, intent(out) :: yc_size !< Size of y dimension of compute domain.
  logical, intent(out) :: buffer_includes_halos !< Flag telling if input buffer includes space for halos.

  integer :: xd_size
  integer :: yd_size

  call mpp_get_data_domain(io_domain, xbegin=isd, xsize=xd_size, ybegin=jsd, &
                           ysize=yd_size, position=dpos)
  call mpp_get_compute_domain(io_domain, xbegin=isc, xsize=xc_size, ybegin=jsc, &
                              ysize=yc_size, position=dpos)
  buffer_includes_halos = (data_xsize .eq. xd_size) .and. (data_ysize .eq. yd_size)
  if (.not. buffer_includes_halos .and. data_xsize .ne. xc_size .and. data_ysize &
      .ne. yc_size) then
    call error("size of x dimension of input buffer does not match size" &
               //" of x dimension of data or compute domain.")
  endif
end subroutine domain_offsets


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_2d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(inout) :: vdata !< Data that will
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
  integer :: dpos
  integer :: i
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer(kind=int32), dimension(:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:), allocatable :: buf_real64
  logical :: buffer_includes_halos

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  dpos = get_domain_decomposed_variable_position(fileobj, variable_name)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), io_domain, &
                      dpos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root reads in the data and scatters it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, &
                                 ybegin=pe_jsc, ysize=pe_jcsize, position=dpos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i) - pe_isc(1) + 1
      c(ydim_index) = pe_jsc(i) - pe_jsc(1) + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int32, e)
          call netcdf_read_data(fileobj, variable_name, buf_int32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_int32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int32, size(buf_int32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int32)
        type is (integer(kind=int64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int64, e)
          call netcdf_read_data(fileobj, variable_name, buf_int64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_int64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int64, size(buf_int64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int64)
        type is (real(kind=real32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real32, e)
          call netcdf_read_data(fileobj, variable_name, buf_real32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_real32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real32, size(buf_real32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real32)
        type is (real(kind=real64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real64, e)
          call netcdf_read_data(fileobj, variable_name, buf_real64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_real64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real64, size(buf_real64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real64)
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
      type is (integer(kind=int32))
        call allocate_array(buf_int32, e)
        call mpp_recv(buf_int32, size(buf_int32), fileobj%io_root, block=.true.)
        call put_array_section(buf_int32, vdata, c, e)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, e)
        call mpp_recv(buf_int64, size(buf_int64), fileobj%io_root, block=.true.)
        call put_array_section(buf_int64, vdata, c, e)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call allocate_array(buf_real32, e)
        call mpp_recv(buf_real32, size(buf_real32), fileobj%io_root, block=.true.)
        call put_array_section(buf_real32, vdata, c, e)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call allocate_array(buf_real64, e)
        call mpp_recv(buf_real64, size(buf_real64), fileobj%io_root, block=.true.)
        call put_array_section(buf_real64, vdata, c, e)
        deallocate(buf_real64)
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine domain_read_2d


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_3d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(inout) :: vdata !< Data that will
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
  integer :: dpos
  integer :: i
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(3) :: c
  integer, dimension(3) :: e
  integer(kind=int32), dimension(:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:), allocatable :: buf_real64
  logical :: buffer_includes_halos

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  dpos = get_domain_decomposed_variable_position(fileobj, variable_name)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), io_domain, &
                      dpos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root reads in the data and scatters it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, &
                                 ybegin=pe_jsc, ysize=pe_jcsize, position=dpos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i) - pe_isc(1) + 1
      c(ydim_index) = pe_jsc(i) - pe_jsc(1) + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int32, e)
          call netcdf_read_data(fileobj, variable_name, buf_int32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_int32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int32, size(buf_int32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int32)
        type is (integer(kind=int64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int64, e)
          call netcdf_read_data(fileobj, variable_name, buf_int64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_int64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int64, size(buf_int64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int64)
        type is (real(kind=real32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real32, e)
          call netcdf_read_data(fileobj, variable_name, buf_real32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_real32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real32, size(buf_real32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real32)
        type is (real(kind=real64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real64, e)
          call netcdf_read_data(fileobj, variable_name, buf_real64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_real64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real64, size(buf_real64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real64)
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
      type is (integer(kind=int32))
        call allocate_array(buf_int32, e)
        call mpp_recv(buf_int32, size(buf_int32), fileobj%io_root, block=.true.)
        call put_array_section(buf_int32, vdata, c, e)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, e)
        call mpp_recv(buf_int64, size(buf_int64), fileobj%io_root, block=.true.)
        call put_array_section(buf_int64, vdata, c, e)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call allocate_array(buf_real32, e)
        call mpp_recv(buf_real32, size(buf_real32), fileobj%io_root, block=.true.)
        call put_array_section(buf_real32, vdata, c, e)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call allocate_array(buf_real64, e)
        call mpp_recv(buf_real64, size(buf_real64), fileobj%io_root, block=.true.)
        call put_array_section(buf_real64, vdata, c, e)
        deallocate(buf_real64)
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine domain_read_3d


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_4d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(inout) :: vdata !< Data that will
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
  integer :: dpos
  integer :: i
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(4) :: c
  integer, dimension(4) :: e
  integer(kind=int32), dimension(:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:), allocatable :: buf_real64
  logical :: buffer_includes_halos

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  dpos = get_domain_decomposed_variable_position(fileobj, variable_name)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), io_domain, &
                      dpos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root reads in the data and scatters it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, &
                                 ybegin=pe_jsc, ysize=pe_jcsize, position=dpos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i) - pe_isc(1) + 1
      c(ydim_index) = pe_jsc(i) - pe_jsc(1) + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int32, e)
          call netcdf_read_data(fileobj, variable_name, buf_int32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_int32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int32, size(buf_int32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int32)
        type is (integer(kind=int64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int64, e)
          call netcdf_read_data(fileobj, variable_name, buf_int64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_int64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int64, size(buf_int64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int64)
        type is (real(kind=real32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real32, e)
          call netcdf_read_data(fileobj, variable_name, buf_real32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_real32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real32, size(buf_real32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real32)
        type is (real(kind=real64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real64, e)
          call netcdf_read_data(fileobj, variable_name, buf_real64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_real64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real64, size(buf_real64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real64)
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
      type is (integer(kind=int32))
        call allocate_array(buf_int32, e)
        call mpp_recv(buf_int32, size(buf_int32), fileobj%io_root, block=.true.)
        call put_array_section(buf_int32, vdata, c, e)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, e)
        call mpp_recv(buf_int64, size(buf_int64), fileobj%io_root, block=.true.)
        call put_array_section(buf_int64, vdata, c, e)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call allocate_array(buf_real32, e)
        call mpp_recv(buf_real32, size(buf_real32), fileobj%io_root, block=.true.)
        call put_array_section(buf_real32, vdata, c, e)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call allocate_array(buf_real64, e)
        call mpp_recv(buf_real64, size(buf_real64), fileobj%io_root, block=.true.)
        call put_array_section(buf_real64, vdata, c, e)
        deallocate(buf_real64)
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine domain_read_4d


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_5d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(inout) :: vdata !< Data that will
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
  integer :: dpos
  integer :: i
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(5) :: c
  integer, dimension(5) :: e
  integer(kind=int32), dimension(:,:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:,:), allocatable :: buf_real64
  logical :: buffer_includes_halos

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  dpos = get_domain_decomposed_variable_position(fileobj, variable_name)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), io_domain, &
                      dpos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root reads in the data and scatters it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, &
                                 ybegin=pe_jsc, ysize=pe_jcsize, position=dpos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i) - pe_isc(1) + 1
      c(ydim_index) = pe_jsc(i) - pe_jsc(1) + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int32, e)
          call netcdf_read_data(fileobj, variable_name, buf_int32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_int32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int32, size(buf_int32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int32)
        type is (integer(kind=int64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int64, e)
          call netcdf_read_data(fileobj, variable_name, buf_int64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_int64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int64, size(buf_int64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int64)
        type is (real(kind=real32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real32, e)
          call netcdf_read_data(fileobj, variable_name, buf_real32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_real32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real32, size(buf_real32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real32)
        type is (real(kind=real64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real64, e)
          call netcdf_read_data(fileobj, variable_name, buf_real64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            endif
            call put_array_section(buf_real64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real64, size(buf_real64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real64)
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
      type is (integer(kind=int32))
        call allocate_array(buf_int32, e)
        call mpp_recv(buf_int32, size(buf_int32), fileobj%io_root, block=.true.)
        call put_array_section(buf_int32, vdata, c, e)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, e)
        call mpp_recv(buf_int64, size(buf_int64), fileobj%io_root, block=.true.)
        call put_array_section(buf_int64, vdata, c, e)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call allocate_array(buf_real32, e)
        call mpp_recv(buf_real32, size(buf_real32), fileobj%io_root, block=.true.)
        call put_array_section(buf_real32, vdata, c, e)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call allocate_array(buf_real64, e)
        call mpp_recv(buf_real64, size(buf_real64), fileobj%io_root, block=.true.)
        call put_array_section(buf_real64, vdata, c, e)
        deallocate(buf_real64)
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine domain_read_5d
